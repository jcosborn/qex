## Actions for HMC
##
## Author: Curtis T. Peterson
##
## Details:
## Implements actions/forces for use in HMC. Based on many iterations/versions of
## similar code that Curtis, James Osborn, and Xiao-Yong Jin have written over the
## years. Currently only implements gauge and staggered fermion actions without
## rooting. But the template is there for more complicated actions. Specifically,
## see https://github.com/ctpeterson/qex/tree/devel/src/mcmc for how to implement
## rooted actions.

import field
import rng
import layout

import mdevolve

import std/[os,strutils]

import algorithms/[integrator]

import base/[qexInternal]
import base/[basicOps]

import gauge/[gaugeAction]
import gauge/[gaugeUtils]

import gauge/[hypsmear]
import gauge/[stoutsmear]
import gauge/[hisqsmear]

import physics/[stagD]
import physics/[stagSolve]

when defined(HisqSmearing):
  import maths/[matproject]

import hmc/[metropolis]

export metropolis

type
  ActionField* = ref object of RootObj
  ActionRngField* = ref object of ActionField
  ActionForce* = ref object of ActionField

type ActionRoot* = ref object of RootObj
  heatbathProc*: proc(self: ActionRoot, u: ActionField, r: ActionRngField)
  actionProc*: proc(self: ActionRoot, u: ActionField): float
  forceProc*: proc(self: ActionRoot, u: ActionField, dtau: float, f: ActionForce)
  name*: string
  id*: string
  description*: string
  stats*: Table[string,Table[string,float]]

type GaugeConfiguration*[U] = ref object of ActionField
  u*: seq[U]
  when defined(HypSmearing):
    su*: seq[U]
  elif defined(StoutSmearing):
    discard
  elif defined(HisqSmearing):
    su*, sul*: seq[U]

  when defined(HypSmearing):
    sc*: HypCoefs
  elif defined(StoutSmearing):
    discard
  elif defined(HisqSmearing):
    sc*: HisqCoefs

  when defined(HypSmearing):
    deriv*: proc(f: seq[U]; chain: seq[U])
  elif defined(StoutSmearing):
    discard
  elif defined(HisqSmearing):
    deriv*: proc(dsdu: var seq[U]; dsdsu, dsdsul: seq[U])

type
  ActionPrng*[R] = ref object of ActionRngField
    r*: R

  GaugeForce*[U] = ref object of ActionForce
    f*: seq[U]

  GaugeAction*[U] = ref object of ActionRoot
    gc*: GaugeActionCoeffs
    aperf*: PerfInfo
    fperf*: PerfInfo

  StaggeredFermionAction*[U, T, S, R] = ref object of ActionRoot
    mass*: float
    stag*: Staggered[U, T]
    phi*: S
    prng*: R
    spa*: SolverParams
    spf*: SolverParams

  StaggeredPauliVillarsAction*[U, T, S] = ref object of ActionRoot
    mass*: float
    stag*: Staggered[U, T]
    phi*: S
    spa*: SolverParams
    spf*: SolverParams

  StaggeredRatioAction*[U, T, S] = ref object of ActionRoot
    massNum*: float
    massDen*: float
    stagNum*: Staggered[U, T]
    stagDen*: Staggered[U, T]
    phi*: S
    spa*: SolverParams
    spf*: SolverParams

type
  ActionLevel* = object
    actions*: seq[ActionRoot]
    multiplier*: int
    integrator*: IntegratorProc

  HmcAction*[U, R] = ref object of MetropolisRoot
    tau*: float
    uc*: GaugeConfiguration[U]
    prng*: R
    levels*: seq[ActionLevel]
    p*: seq[U]
    f*: seq[U]
    bu*: seq[U]
    heatbathProc*: proc()
    globalRandProc*: proc(): float
    hmcStats*: seq[proc():string]
    secs*: float

  #HmcEvolver*[U, R] = ref object of MetropolisRoot
  #  tau*: float
  #  uc*: GaugeConfiguration[U]
  #  prng*: R
  #  levels*: seq[ActionLevel]
  #  p*: seq[U]
  #  f*: seq[U]
  #  bu*: seq[U]
  #  heatbathProc*: proc()
  #  globalRandProc*: proc(): float
  #  hmcStats*: seq[proc():string]

# HmcAction
#   tuple(fields)
#   levels
#   stats

# HmcRunner
#   tau, globalRand

proc newHmcAction*[U,R](
  uc: GaugeConfiguration[U];
  srng: auto;
  prng: R;
  tau: float
): HmcAction[U,R] =
  new(result)
  result.uc = uc
  result.prng = prng
  let lo = uc.u[0].l
  result.p = lo.newGauge()
  result.f = lo.newGauge()
  result.bu = lo.newGauge()
  result.hmcStats.newSeq(0)
  result.tau = tau
  #template ET: untyped = evalType(lo.ColorVector()[0])
  #template FT: untyped = evalType(lo.ColorVector())

  let pHB = result.p
  result.heatbathProc = proc =
    var p = pHB
    threads:
      for mu in 0..<p.len: p[mu].randomTAH(prng)

  var sRng = srng
  result.globalRandProc = proc(): float = sRng.uniform()

proc description*(h: HmcAction): string =
  result = ""
  for i,l in h.levels.pairs:
    if i>0: result &= "\n"
    result &= "ActionLevel " & $i
    for a in l.actions:
      if a.id != "":
        result &= "\n" & a.id
        if a.description != "":
          result &= "\n" & a.description.indent(2)

#[ ActionLevel: virtual dispatch for nested integrators ]#

proc newActionLevel*(multiplier: int = 1; integrator: IntegratorProc = "2MN"): ActionLevel =
  ActionLevel(actions: @[], multiplier: multiplier, integrator: integrator)

proc add*(level: var ActionLevel; action: ActionRoot) =
  level.actions.add(action)

proc heatbath(level: ActionLevel, u: ActionField, r: ActionRngField) =
  for a in level.actions: a.heatbathProc(a, u, r)

proc action(level: ActionLevel, u: ActionField): float =
  for a in level.actions: result += a.actionProc(a, u)

proc force(level: ActionLevel, u: ActionField, dtau: float, f: ActionForce) =
  for a in level.actions: a.forceProc(a, u, dtau, f)

proc add*(hmc: var HmcAction; level: ActionLevel) =
  hmc.levels.add(level)

proc heatbath(hmc: var HmcAction; rng: auto) =
  var p = hmc.p
  threads:
    for mu in 0..<p.len: p[mu].randomTAH(rng)
  for level in hmc.levels: level.heatbath()

proc kineticAction*(hmc: HmcAction): float =
  var p2: float
  threads:
    var p2t = 0.0
    for mu in 0..<hmc.p.len: p2t += hmc.p[mu].norm2
    threadBarrier()
    threadMaster: p2 = p2t
  return 0.5*p2 - 16.0*float(hmc.p[0].l.physVol)

proc action*(hmc: var HmcAction): float =
  result = 0.0
  for level in hmc.levels: result += level.action(hmc.uc)

proc hamiltonian*(hmc: var HmcAction): float = hmc.kineticAction() + hmc.action()

proc integrator*(hmc: var HmcAction): Integrator =
  let uc = hmc.uc
  let pp = addr hmc.p
  let fp = addr hmc.f
  let gf = GaugeForce[hmc.U](f: fp[])
  let levels = addr hmc.levels
  let nlevels = hmc.levels.len

  proc mdt(dtau: float) =
    threads:
      for mu in 0..<uc.u.len:
        for s in uc.u[mu]:
          uc.u[mu][s] := exp(dtau * pp[][mu][s]) * uc.u[mu][s]

  proc mdv(dtau: openArray[float]) =
    threads:
      for mu in 0..<fp[].len: fp[][mu] := 0
    for i in 0..<nlevels:
      if dtau[i] != 0.0: levels[][i].force(uc, dtau[i], gf)
    threads:
      for mu in 0..<pp[].len: pp[][mu] -= fp[][mu]

  let (V, T) = newIntegratorPair(mdv, mdt)
  var integrator = T
  for i in 0..<nlevels:
    integrator = levels[][i].integrator(
      T = integrator,
      V = V[i],
      steps = levels[][i].multiplier
    )

  return integrator

proc setGauge[U](g: seq[U]; u: seq[U]) =
  threads:
    for mu in 0..<g.len: g[mu] := u[mu]

proc reunit(g: auto) =
  threads:
    let d = g.checkSU
    threadBarrier()
    echo "unitary deviation avg: ",d.avg," max: ",d.max
    g.projectSU
    threadBarrier()
    let dd = g.checkSU
    echo "new unitary deviation avg: ",dd.avg," max: ",dd.max

proc reunit*(hmc: var HmcAction) = reunit(hmc.uc.u)

#[ "virtual" MetropolisRootObj procedures ]#

proc getH*(hmc: var HmcAction): float = hmc.hamiltonian()

proc start*(hmc: var HmcAction) =
  hmc.heatbathProc()
  let r = ActionPrng[hmc.R](r: hmc.prng)
  for level in hmc.levels: level.heatbath(hmc.uc, r)
  setGauge(hmc.bu, hmc.uc.u)

proc generate*(hmc: var HmcAction) =
  var integ = hmc.integrator()
  integ.evolve(hmc.tau)
  integ.finish()

proc globalRand*(hmc: var HmcAction): float = hmc.globalRandProc()

proc accept*(hmc: var HmcAction) = hmc.reunit()

proc reject*(hmc: var HmcAction) = setGauge(hmc.uc.u, hmc.bu)

proc run*(hmc: var HmcAction) =
  tic("HmcAction:run")
  let nup = hmc.nUpdates + 1
  echo &"== Begin HMC update {nup} =========="
  metropolis.update(hmc)
  let dt = getElapsedTime()
  hmc.secs += dt
  echo &"End HMC update {nup}: {dt:.2f} seconds ({hmc.secs:.2f} total)"
  #for s in hmc.hmcStats:
  #  echo s()
  #for level in hmc.levels:
  #  for a in level.actions:
      #echo &"{a.name}"
      #let s = a.stats(a)
      #if s != "": echo s
  #echo &"  {}"
  toc("end")

#[ gauge configuration implementation ]#

when defined(HypSmearing):
  proc newGaugeConfiguration*[U](
    u: seq[U];
    alpha1: float = 0.4;
    alpha2: float = 0.5;
    alpha3: float = 0.5
  ): GaugeConfiguration[U] =
    let lo = u[0].l
    result = GaugeConfiguration[U](u: u)
    result.sc = HypCoefs(alpha1: alpha1, alpha2: alpha2, alpha3: alpha3)
    result.su = lo.newGauge()
elif defined(StoutSmearing):
  proc newGaugeConfiguration*[U](u: seq[U]): GaugeConfiguration[U] =
    qexError "Stout smearing for HMC not yet implemented"
elif defined(HisqSmearing):
  proc newGaugeConfiguration*[U](
    u: seq[U];
    naik: float = 1.0;
    lepage: float = 0.0;
    unitaryItr: int = 10;
    unitaryEps: float = 1e-20;
    unitaryPrj: ProjectionMethod = CayleyHamilton
  ): GaugeConfiguration[U] =
    let lo = u[0].l
    result = GaugeConfiguration[U](u: u)

    result.sc = HisqCoefs()
    result.sc.fat7first.setHisqFat7(lepage, 0.0)
    result.sc.fat7second.setHisqFat7(2.0 - lepage, naik)
    result.sc.naik = -naik/24.0
    result.sc.projection = newUnitaryProjection(unitaryPrj, unitaryEps, unitaryItr)

    result.su = lo.newGauge()
    result.sul = lo.newGauge()
else:
  proc newGaugeConfiguration*[U](u: seq[U]): GaugeConfiguration[U] =
    result = GaugeConfiguration[U](u: u)

proc cold*(c: var GaugeConfiguration) = c.u.unit()

proc read*(c: var GaugeConfiguration; filename: string) =
  if fileExists(filename):
    if 0 != c.u.loadGauge(filename):
      qexError "error loading gauge configuration from file: ", filename
    qexLog "gauge configuration loaded from file: ", filename
  else: qexError "gauge configuration file does not exist: ", filename

proc smear(self: GaugeConfiguration) =
  # HISQ is a staggered smearing, so I don't mind having to explicitly insert the
  # the staggered rephasing here. For other smearings, rephasing must be decoupled
  # from the smearing, which is the case for nHYP
  when defined(HypSmearing):
    var info: PerfInfo
    discard self.sc.smearGetForce(self.u, self.su, info)
  elif defined(StoutSmearing):
    qexError "Stout smearing for HMC not yet implemented"
  elif defined(HisqSmearing):
    threads:
      self.u.setBC()
      threadBarrier()
      self.u.stagPhase()
    discard self.sc.smearGetForce(self.u, self.su, self.sul)
    threads:
      self.u.setBC()
      threadBarrier()
      self.u.stagPhase()

proc smearGetForce[U](self: GaugeConfiguration[U]) =
  when defined(HypSmearing):
    var info: PerfInfo
    self.deriv = self.sc.smearGetForce(self.u, self.su, info)
  elif defined(StoutSmearing):
    qexError "Stout smearing for HMC not yet implemented"
  elif defined(HisqSmearing):
    threads:
      self.u.setBC()
      threadBarrier()
      self.u.stagPhase()
    self.deriv = self.sc.smearGetForce(self.u, self.su, self.sul)
    threads:
      self.u.setBC()
      threadBarrier()
      self.u.stagPhase()

proc setBC[U](self: GaugeConfiguration[U]) =
  threads:
    when defined(HypSmearing): self.su.setBC()
    elif defined(StoutSmearing): qexError "Stout smearing for HMC not yet implemented"
    elif defined(HisqSmearing): discard
    else: self.u.setBC()

proc stagPhase[U](self: GaugeConfiguration[U]) =
  threads:
    when defined(HypSmearing): self.su.stagPhase()
    elif defined(StoutSmearing): qexError "Stout smearing for HMC not yet implemented"
    elif defined(HisqSmearing): discard
    else: self.u.stagPhase()

#[ gauge action implementation ]#

proc action*(self: GaugeAction, u: GaugeConfiguration): float =
  tic("GaugeAction:action")
  result = self.gc.gaugeAction1(u.u)
  self.aperf.secs += getElapsedTime()
  toc("end")

proc force*(self: GaugeAction, u: GaugeConfiguration): auto =
  tic("GaugeAction:force")
  let lo = u.u[0].l
  var f = lo.newGauge()
  self.gc.gaugeForce(u.u, f)
  self.fperf.secs += getElapsedTime()
  toc("end")
  return f

proc gaugeActionHeatbathProc(self: ActionRoot, u: ActionField, r: ActionRngField) = discard

proc gaugeActionActionProc[U](self: ActionRoot, u: ActionField): float =
  let ga = GaugeAction[U](self)
  let gc = GaugeConfiguration[U](u)
  ga.action(gc)

proc gaugeActionForceProc[U](self: ActionRoot, u: ActionField, dtau: float, f: ActionForce) =
  let ga = GaugeAction[U](self)
  let gc = GaugeConfiguration[U](u)
  let gf = GaugeForce[U](f)
  let fi = ga.force(gc)
  threads:
    for mu in 0..<gf.f.len: gf.f[mu] += dtau * fi[mu]

proc gaugeActionStatsA(a: GaugeAction): string =
  result = &"GaugeAction: {a.aperf.secs:.2f}"
  a.aperf.clear

proc gaugeActionStatsF(a: GaugeAction): string =
  result = &"GaugeForce:  {a.fperf.secs:.2f}"
  a.fperf.clear

var GaugeActionCount = 0
proc newGaugeAction*[U](h: HmcAction, gc: GaugeActionCoeffs, u: GaugeConfiguration[U]): GaugeAction[U] =
  result = GaugeAction[U](gc: gc)
  result.name = "GaugeAction"
  result.id = "GA" & $GaugeActionCount; inc GaugeActionCount
  result.description = result.name & $gc
  result.heatbathProc = gaugeActionHeatbathProc
  result.actionProc = gaugeActionActionProc[U]
  result.forceProc = gaugeActionForceProc[U]
  let a = result
  h.hmcStats.add proc():string = gaugeActionStatsA(a)
  h.hmcStats.add proc():string = gaugeActionStatsF(a)

#[ fermion force helpers ]#

proc stagForceSolve[T](
  stag: auto;
  psi: auto;
  phi: T;
  mass: float;
  sp0: var SolverParams
) =
  tic()
  var
    varphi = newOneOf(psi)
    r = newOneOf(phi)
    r2, b2: float
    sp = sp0

  # Prepare solver parameters
  sp.resetStats()
  dec sp.verbosity
  sp.usePrevSoln = false
  threads:
    psi := 0
    varphi := 0
    r := 0
    let b2t = phi.norm2
    threadBarrier()
    threadMaster: b2 = b2t

  # Get solution
  stag.solveEE(varphi, phi, mass, sp)
  threads:
    # Get residual
    var r2t: float
    stagD2ee(stag.se, stag.so, r, stag.g, varphi, mass*mass)
    threadBarrier()
    r := r - phi
    threadBarrier()
    r2t = r.norm2
    threadBarrier()
    threadMaster: r2 = r2t

    # Get solution for force
    varphi.even := 4.0*varphi
    threadBarrier()
    stagD2(stag.so, psi, stag.g, varphi, 0, 0)
    threadBarrier()
    psi.even := varphi

  # Get information about solve
  sp.r2.init r2/b2
  sp.calls = 1
  sp.seconds = getElapsedTime()
  sp.flops += float((stag.g.len*4*72+24)*psi.l.nEven) # ???
  if sp0.verbosity>0: echo "stagSolve: ", sp.getStats
  sp0.addStats(sp)

proc fermForce[U, S](f: seq[U]; psi: S; u: GaugeConfiguration[U]) =
  let nc = psi[0].len
  let nd = psi.l.nDim
  var
    ff = f.newOneOf()
    f1 = f.newOneOf()
    t1 = newSeq[Shifter[typeOf(psi), typeOf(psi[0])]](nd)
  when defined(HisqSmearing):
    var
      f3 = f.newOneOf()
      t3 = newSeq[Shifter[typeOf(psi), typeOf(psi[0])]](nd)

  for mu in 0..<nd:
    t1[mu] = newShifter(psi, mu, 1)
    discard t1[mu] ^* psi
    when defined(HisqSmearing):
      t3[mu] = newShifter(psi, mu, 3)
      discard t3[mu] ^* psi

  # dslash
  threads:
    for mu in 0..<nd:
      for n in f[mu]:
        forO a, 0, nc-1:
          forO b, 0, nc-1:
            f1[mu][n][a, b] := psi[n][a] * t1[mu].field[n][b].adj
            when defined(HisqSmearing):
              f3[mu][n][a, b] := psi[n][a] * t3[mu].field[n][b].adj
        f[mu][n] := 0
        ff[mu][n] := 0
    threadBarrier()

    # rephase
    when defined(HypSmearing):
      f1.setBC()
      threadBarrier()
      f1.stagPhase()
    elif defined(StoutSmearing):
      qexError "Stout smearing for HMC not yet implemented"
    elif defined(HisqSmearing):
      u.u.setBC()
      threadBarrier()
      u.u.stagPhase()
    # else: u rephasing already takes care of phases
    threadBarrier()

    # correct odd sites
    for mu in 0..<nd:
      for n in f[mu].odd:
        f1[mu][n] *= -1
        when defined(HisqSmearing):
          f3[mu][n] *= -1

  # smear
  when defined(HypSmearing): u.deriv(ff, f1)
  elif defined(StoutSmearing): qexError "Stout smearing for HMC not yet implemented"
  elif defined(HisqSmearing): u.deriv(ff, f1, f3)
  else: # no smearing - identity chain rule
    threads:
      for mu in 0..<ff.len: ff[mu] := f1[mu]

  # traceless/anti-Hermitian projection
  threads:
    for mu in 0..<f.len:
      for n in f[mu]: f1[mu][n] := ff[mu][n] * u.u[mu][n].adj
    threadBarrier()
    when defined(HisqSmearing):
      u.u.setBC()
      threadBarrier()
      u.u.stagPhase()
    for mu in 0..<f.len:
      for n in f[mu]: f[mu][n].projectTAH(f1[mu][n])

#[ staggered fermion action ]#

proc heatbath*(self: StaggeredFermionAction; u: GaugeConfiguration) =
  let lo = self.phi.l
  var psi = lo.ColorVector()
  u.smear()
  u.setBC()
  u.stagPhase()
  threads:
    psi.gaussian(self.prng)
    threadBarrier()
    self.stag.D(self.phi, psi, -self.mass)
    threadBarrier()
    self.phi.odd := 0
  u.setBC()
  u.stagPhase()

proc action*(self: StaggeredFermionAction; u: GaugeConfiguration): float =
  let lo = self.phi.l
  var psi = lo.ColorVector()
  var psi2: float
  threads: psi := 0
  u.smear()
  u.setBC()
  u.stagPhase()
  self.stag.solve(psi, self.phi, -self.mass, self.spa)
  threads:
    var psi2t = psi.norm2()
    threadBarrier()
    threadMaster: psi2 = psi2t
  u.setBC()
  u.stagPhase()
  return 0.5*psi2

proc force*(self: StaggeredFermionAction; u: GaugeConfiguration): auto =
  let lo = self.phi.l
  var f = lo.newGauge()
  var psi = lo.ColorVector()
  u.smearGetForce()
  u.setBC()
  u.stagPhase()
  self.stag.stagForceSolve(psi, self.phi, self.mass, self.spf)
  f.fermForce(psi, u)
  u.setBC()
  u.stagPhase()
  threads:
    for mu in 0..<f.len:
      for n in f[mu]: f[mu][n] *= 0.25
  return f

proc staggeredFermionActionHeatbathProc[U,ET,FT,R](self: ActionRoot, u: ActionField,
                                                   r: ActionRngField) =
  let fa = StaggeredFermionAction[U,ET,FT,R](self)
  let gc = GaugeConfiguration[U](u)
  fa.heatbath(gc)

proc staggeredFermionActionActionProc[U,ET,FT,R](self: ActionRoot, u: ActionField): float =
  let fa = StaggeredFermionAction[U, ET, FT, R](self)
  let gc = GaugeConfiguration[U](u)
  fa.action(gc)

proc staggeredFermionActionForceProc[U,ET,FT,R](self: ActionRoot, u: ActionField, dtau: float,
                                              f: ActionForce) =
  let fa = StaggeredFermionAction[U, ET, FT, R](self)
  let gc = GaugeConfiguration[U](u)
  let gf = GaugeForce[U](f)
  let fi = fa.force(gc)
  threads:
    for mu in 0..<gf.f.len: gf.f[mu] += dtau * fi[mu]

var StaggeredFermionActionCount = 0
proc newStaggeredFermionAction*[U, T, R](
  stag: Staggered[U, T];
  mass: float;
  spa: SolverParams;
  spf: SolverParams;
  r: R
): auto =
  let lo = stag.g[0].l
  var phi = lo.newField(T)
  type FT = evalType(phi)
  result = StaggeredFermionAction[U, T, FT, R](
    mass: mass,
    stag: stag,
    phi: phi,
    prng: r,
    spa: spa,
    spf: spf
  )
  result.name = "StaggeredFermionAction"
  result.id = "SFA" & $StaggeredFermionActionCount; inc StaggeredFermionActionCount
  result.description = result.name & &"(mass: {mass})"
  result.heatbathProc = staggeredFermionActionHeatbathProc[U,T,FT,R]
  result.actionProc = staggeredFermionActionActionProc[U,T,FT,R]
  result.forceProc = staggeredFermionActionForceProc[U,T,FT,R]

#[ staggered Pauli-Villars action ]#

proc heatbath*(self: StaggeredPauliVillarsAction; u: GaugeConfiguration; rng: auto) =
  let lo = self.phi.l
  var psi = lo.ColorVector()
  u.smear()
  u.setBC()
  u.stagPhase()
  threads:
    psi.gaussian(rng)
    self.phi := 0
  self.stag.solve(self.phi, psi, self.mass, self.spa)
  threads: self.phi.odd := 0
  u.setBC()
  u.stagPhase()

proc action*(self: StaggeredPauliVillarsAction; u: GaugeConfiguration): float =
  let lo = self.phi.l
  var psi = lo.ColorVector()
  var psi2: float
  u.smear()
  u.setBC()
  u.stagPhase()
  threads:
    psi := 0
    threadBarrier()
    self.stag.D(psi, self.phi, self.mass)
    threadBarrier()
    var psi2t = psi.norm2()
    threadBarrier()
    threadMaster: psi2 = psi2t
  u.setBC()
  u.stagPhase()
  return 0.5*psi2

proc force*[U](self: StaggeredPauliVillarsAction; u: GaugeConfiguration[U]): auto =
  let lo = self.phi.l
  var f = lo.newGauge()
  var psi = lo.ColorVector()
  u.smearGetForce()
  u.setBC()
  u.stagPhase()
  threads:
    psi := 0
    threadBarrier()
    stagD2(self.stag.so, psi, self.stag.g, self.phi, 0, 0)
    threadBarrier()
    psi.even := self.phi
  f.fermForce(psi, u)
  u.setBC()
  u.stagPhase()
  threads:
    for mu in 0..<f.len:
      for n in f[mu]: f[mu][n] *= -0.25
  return f

proc staggeredPauliVillarsActionHeatbathProc[U,ET,FT,R](self: ActionRoot, u: ActionField,
                                                        r: ActionRngField) =
  let pva = StaggeredPauliVillarsAction[U, ET, FT](self)
  let gc = GaugeConfiguration[U](u)
  let ap = ActionPrng[R](r)
  pva.heatbath(gc, ap.r)

proc staggeredPauliVillarsActionActionProc[U,ET,FT](self: ActionRoot, u: ActionField): float =
  let pva = StaggeredPauliVillarsAction[U, ET, FT](self)
  let gc = GaugeConfiguration[U](u)
  pva.action(gc)

proc staggeredPauliVillarsActionForceProc[U,ET,FT](self: ActionRoot, u: ActionField, dtau: float,
                                              f: ActionForce) =
  let pva = StaggeredPauliVillarsAction[U, ET, FT](self)
  let gc = GaugeConfiguration[U](u)
  let gf = GaugeForce[U](f)
  let fi = pva.force(gc)
  threads:
    for mu in 0..<gf.f.len: gf.f[mu] += dtau * fi[mu]

var StaggeredPauliVillarsActionCount = 0
proc newStaggeredPauliVillarsAction*[U, T, R](
  stag: Staggered[U, T];
  mass: float;
  spa: SolverParams;
  spf: SolverParams;
  r: R
): auto =
  let lo = stag.g[0].l
  var phi = lo.newField(T)
  type FT = evalType(phi)
  result = StaggeredPauliVillarsAction[U, T, FT](
    mass: mass,
    stag: stag,
    phi: phi,
    spa: spa,
    spf: spf
  )
  result.name = "StaggeredPauliVillarsAction"
  result.id = "SPVA" & $StaggeredPauliVillarsActionCount; inc StaggeredPauliVillarsActionCount
  result.description = result.name & &"(mass: {mass})"
  result.heatbathProc = staggeredPauliVillarsActionHeatbathProc[U,T,FT,R]
  result.actionProc = staggeredPauliVillarsActionActionProc[U,T,FT]
  result.forceProc = staggeredPauliVillarsActionForceProc[U,T,FT]

#[ staggered "ratio" (Pauli-Villars/fermion) action ]#

proc heatbath*(self: StaggeredRatioAction; u: GaugeConfiguration; rng: auto) =
  let lo = self.phi.l
  var psi = lo.ColorVector()
  var phi = lo.ColorVector()
  u.smear()
  u.setBC()
  u.stagPhase()
  threads:
    self.phi := 0
    psi.gaussian(rng)
    threadBarrier()
    self.stagNum.D(phi, psi, -self.massNum)
  self.stagDen.solve(self.phi, phi, -self.massDen, self.spa)
  threads: self.phi.odd := 0
  u.setBC()
  u.stagPhase()

proc action*(self: StaggeredRatioAction; u: GaugeConfiguration): float =
  let lo = self.phi.l
  var psi = lo.ColorVector()
  var phi = lo.ColorVector()
  var psi2: float
  u.smear()
  u.setBC()
  u.stagPhase()
  threads:
    self.stagDen.D(phi, self.phi, -self.massDen)
    psi := 0
  self.stagNum.solve(psi, phi, -self.massNum, self.spa)
  threads:
    var psi2t = psi.norm2()
    threadBarrier()
    threadMaster: psi2 = psi2t
  u.setBC()
  u.stagPhase()
  return 0.5*psi2

proc force*[U](self: StaggeredRatioAction; u: GaugeConfiguration[U]): auto =
  let lo = self.phi.l
  var f = lo.newGauge()
  var psi = lo.ColorVector()
  u.smearGetForce()
  u.setBC()
  u.stagPhase()
  self.stagNum.stagForceSolve(psi, self.phi, self.massNum, self.spf)
  f.fermForce(psi, u)
  u.setBC()
  u.stagPhase()
  let ffac = 0.25 * (self.massDen*self.massDen - self.massNum*self.massNum)
  threads:
    for mu in 0..<f.len:
      for n in f[mu]: f[mu][n] *= ffac
  return f

proc staggeredRatioActionHeatbathProc[U,ET,FT,R](self: ActionRoot, u: ActionField,
                                                 r: ActionRngField) =
  let ra = StaggeredRatioAction[U, ET, FT](self)
  let gc = GaugeConfiguration[U](u)
  let ap = ActionPrng[R](r)
  ra.heatbath(gc, ap.r)

proc staggeredRatioActionActionProc[U,ET,FT](self: ActionRoot, u: ActionField): float =
  let ra = StaggeredRatioAction[U, ET, FT](self)
  let gc = GaugeConfiguration[U](u)
  ra.action(gc)

proc staggeredRatioActionForceProc[U,ET,FT](self: ActionRoot, u: ActionField, dtau: float,
                                            f: ActionForce) =
  let ra = StaggeredRatioAction[U, ET, FT](self)
  let gc = GaugeConfiguration[U](u)
  let gf = GaugeForce[U](f)
  let fi = ra.force(gc)
  threads:
    for mu in 0..<gf.f.len: gf.f[mu] += dtau * fi[mu]

var StaggeredRatioActionCount = 0
proc newStaggeredRatioAction*[U, T, R](
  stagNum: Staggered[U, T];
  stagDen: Staggered[U, T];
  massNum: float;
  massDen: float;
  spa: SolverParams;
  spf: SolverParams;
  r: R
): auto =
  let lo = stagNum.g[0].l
  var phi = lo.newField(T)
  type FT = evalType(phi)
  result = StaggeredRatioAction[U, T, FT](
    massNum: massNum,
    massDen: massDen,
    stagNum: stagNum,
    stagDen: stagDen,
    phi: phi,
    spa: spa,
    spf: spf
  )
  result.name = "StaggeredRatioAction"
  result.id = "SRA" & $StaggeredRatioActionCount; inc StaggeredRatioActionCount
  result.description = result.name & &"(massNum: {massNum}, massDen: {massDen})"
  result.heatbathProc = staggeredRatioActionHeatbathProc[U,T,FT,R]
  result.actionProc = staggeredRatioActionActionProc[U,T,FT]
  result.forceProc = staggeredRatioActionForceProc[U,T,FT]
