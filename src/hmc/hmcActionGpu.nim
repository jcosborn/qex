## Actions for HMC on the GPU, as hmc/hmcAction with StaggeredSmearing HYP.
##
## The links, momenta, pseudofermions and forces stay on the device, and each
## step of the molecular dynamics runs in kernels: gauge/gaugeGpu for the
## gauge action, its force and the link update, gauge/hypGpu for the
## smearing and its force, physics/stagGpu for the staggered operators and
## solves.  A host HmcAction without actions keeps the random number
## generators, the files and the measurements: the momenta and the
## pseudofermion sources are drawn there in the order of hmcAction and
## uploaded, and the links come back once per trajectory.  The operators and
## solves follow stag.D, stag.solve and stagForceSolve of hmcAction with their
## stopping rules, so trajectories agree with hmcAction to rounding.  The
## sums run in a fixed order (gpuSum, and the CG of stagGpu with fixed), so
## a trajectory repeats bit for bit, as the checkpoint test needs.  One boundary condition holds for all
## fermion actions.
##
## The smearing force is linear in the one link force, so the fermion
## actions of a level add their one link forces and pull the sum back once.
## The force stats of each action need its own pullback: with forceStats = 1
## they come from the first fermion force of each trajectory, whose actions
## are pulled back one by one, with forceStats = 2 from every force, as
## hmcAction does.

import qex
import gauge
import physics/[qcdTypes, stagD, stagSolve, stagGpu]
import gauge/[gaugeGpu, hypGpu, hypsmear]
import hmc/[hmcAction, metropolis]
import algorithms/[integrator]
import backend/accel
import base/metaUtils
import std/[sequtils, strformat, strutils, tables]

export hmcAction

type
  GpuActionKind* = enum gkGauge, gkFermion, gkRatio, gkPV
  GpuAction* = ref object
    kind*: GpuActionKind
    name*, id*, description*: string
    gc*: GaugeActionCoeffs
    mass*: float  ## fermion and PV mass, the numerator mass of a ratio
    massDen*: float  ## the denominator mass of a ratio
    phi*, x*, src*, rhs*: ptr UncheckedArray[float]  ## 6n reals: pseudofermion, solution, source, right-hand side
    spa*, spf*: SolverParams
    stats*: Table[string, ActionStats]
  GpuLevel* = object
    actions*: seq[GpuAction]
    multiplier*: int
    integrator*: IntegratorProc
  HmcGpu*[U,R,S] = ref object of MetropolisRoot
    cpu*: HmcAction[U,R,S]  ## host links, random numbers, files, measurements
    tau*: float
    levels*: seq[GpuLevel]
    coef*: HypCoefs
    gg*, sgo*: GpuGauge[VLEN]  ## thin and smeared links
    hg*: HypGpu[VLEN]
    s*: StagGpu[VLEN,float]
    ss*: StagGpu[VLEN,float32]  ## single precision solver of the force solves, with mixed
    sg*: ptr UncheckedArray[float]  ## signs of the phases and boundary conditions, stagSigns
    p*, f*, bu*, ft*, fl*, g1*, p1*: ptr UncheckedArray[float]  ## momenta, force, saved links, force of an action, one link force, reverse check copies
    w*: ptr UncheckedArray[float]  ## work vector, 6n reals
    red*, hred*: ptr UncheckedArray[float]  ## force stats sums and maxima, device and pinned host
    smeared*: bool  ## sgo and s hold the smeared links of gg.u
    forceStats*: int  ## 1: the per action force stats of the first fermion force of a trajectory, 2: of all
    mixed*: bool  ## force solves in mixed precision
    sample: bool  ## pull the actions of the next fermion force back one by one
    hostNew*: bool  ## the host links are newer than gg.u
    hmcStats*: Table[string, ActionStats]
    secs*: float
    revCheckFreq*: int
    forceAccept*: bool
    atEnd: bool

const nSlot = 128  # slots of the second stage of the force stats sums
const V = VLEN

proc newVec(h: HmcGpu): ptr UncheckedArray[float] =
  cast[ptr UncheckedArray[float]](gpuMalloc(6*h.s.n*sizeof(float)))

proc newHmcGpu*[U,R,S](uc: GaugeConfiguration[U]; srng: S; prng: R; tau: float; revCheckFreq = 0;
                       bc: openArray[bool] = [true, true, true, false]; reals = 18;
                       forceStats = 1; mixed = false; fixed = true): HmcGpu[U,R,S] =
  ## the HMC of hmcAction on the GPU, the fermion boundary conditions bc
  ## (periodic when true) for all fermion actions; reals per link of the
  ## solvers, 18 or 14 (rows 0 and 1 and the determinant); with mixed the
  ## force solves restart single precision CGs in double (solveM), which
  ## changes the forces within f_tol; with fixed the CG dot products add in a
  ## fixed order, so a trajectory repeats bit for bit
  new(result)
  result.forceStats = forceStats
  result.mixed = mixed
  result.cpu = uc.newHmcAction(srng, prng, tau, revCheckFreq)
  result.tau = tau
  result.revCheckFreq = revCheckFreq
  result.coef = uc.sc
  let lo = uc.u[0].l
  var gs = lo.newGauge()
  gs.unit
  gs.setBC(bc)
  gs.stagPhase
  result.sg = stagSigns(gs)
  result.s = newStagGpu(gs, float64, reals, batch = true)
  result.s.fixed = fixed
  if mixed:
    result.ss = newStagGpu(gs, float32, reals, batch = true)
    result.ss.fixed = fixed
  result.gg = newGpuGauge(lo)
  result.sgo = newGpuGauge(lo)
  result.hg = newHypGpu(lo)
  for d in [addr result.p, addr result.f, addr result.bu, addr result.ft, addr result.fl,
            addr result.g1, addr result.p1]:
    d[] = result.gg.newLinks
  result.w = result.newVec
  let nt = (4*result.gg.n + 15) div 16
  result.red = cast[ptr UncheckedArray[float]](gpuMalloc((3*nSlot + 3*nt)*sizeof(float)))
  result.hred = cast[ptr UncheckedArray[float]](gpuMallocHost(3*nSlot*sizeof(float)))
  result.hmcStats = initTable[string, ActionStats]()
  result.hostNew = true

proc description*(h: HmcGpu): string =
  var sp = ""
  let n = h.levels.len
  for i in 0..<n:
    let l = h.levels[n-1-i]
    if i > 0: result &= "\n"
    result &= sp & "ActionLevel " & $i & ": "
    sp &= "  "
    for a in l.actions:
      if a.id != "" and a.description != "":
        result &= "\n" & sp & a.id & ": " & a.description.indent(2)

proc newGpuLevel*(multiplier = 1; integrator: IntegratorProc = "2MN"): GpuLevel =
  GpuLevel(multiplier: multiplier, integrator: integrator)

proc add*(level: var GpuLevel; a: GpuAction) = level.actions.add a

proc add*(h: HmcGpu; level: GpuLevel) =
  ## the levels from the innermost, as hmcAction
  h.levels.insert level

var gaugeCount, fermCount, ratioCount, pvCount = 0

proc newGaugeAction*(h: HmcGpu; gc: GaugeActionCoeffs): GpuAction =
  result = GpuAction(kind: gkGauge, gc: gc, name: "GaugeAction", id: "GA" & $gaugeCount)
  inc gaugeCount
  result.description = result.name & $gc

proc newFermVecs(h: HmcGpu; a: GpuAction) =
  a.phi = h.newVec
  a.x = h.newVec
  a.src = h.newVec
  a.rhs = h.newVec

proc newStaggeredFermionAction*(h: HmcGpu; mass: float; spa, spf: SolverParams): GpuAction =
  ## |M(-m)^-1 phi|^2/2 for phi on the even sites, as newStaggeredFermionAction
  result = GpuAction(kind: gkFermion, mass: mass, spa: spa, spf: spf,
                     name: "StaggeredFermionAction", id: "SFA" & $fermCount)
  inc fermCount
  result.description = result.name & &"(mass: {mass})"
  h.newFermVecs result

proc newStaggeredRatioAction*(h: HmcGpu; massNum, massDen: float; spa, spf: SolverParams): GpuAction =
  ## |M(-m_num)^-1 M(-m_den) phi|^2/2, as newStaggeredRatioAction
  result = GpuAction(kind: gkRatio, mass: massNum, massDen: massDen, spa: spa, spf: spf,
                     name: "StaggeredRatioAction", id: "SRA" & $ratioCount)
  inc ratioCount
  result.description = result.name & &"(massNum: {massNum}, massDen: {massDen})"
  h.newFermVecs result

proc newStaggeredPauliVillarsAction*(h: HmcGpu; mass: float; spa, spf: SolverParams): GpuAction =
  ## |M(m) phi|^2/2, as newStaggeredPauliVillarsAction
  result = GpuAction(kind: gkPV, mass: mass, spa: spa, spf: spf,
                     name: "StaggeredPauliVillarsAction", id: "SPVA" & $pvCount)
  inc pvCount
  result.description = result.name & &"(mass: {mass})"
  h.newFermVecs result

#[ device helpers ]#

proc zeroOdd(h: HmcGpu; v: ptr UncheckedArray[float]) =
  let o = 6*h.s.ne
  gpuFor(i, 6*(h.s.n - h.s.ne)): v[o + i] = 0.0

proc zeroLinks(h: HmcGpu; v: ptr UncheckedArray[float]) =
  gpuFor(i, 4*18*h.gg.n): v[i] = 0.0

proc normEO(h: HmcGpu; a: ptr UncheckedArray[float]): array[2, float] =
  ## global |a_e|^2, |a_o|^2
  let ne6 = 6*h.s.ne
  result = gpuSum(i, 6*h.s.n, 2):
    let q = a[i]*a[i]
    if i < ne6: [q, 0.0] else: [0.0, q]
  getDefaultComm().allReduce(addr result[0], 2)

proc resid(h: HmcGpu; x, b: ptr UncheckedArray[float]; m: float): float =
  ## |b - M(m) x|^2/|b|^2 over all sites
  let w = h.w
  h.s.applyMfull(w, x, m)
  var r = gpuSum(i, 6*h.s.n, 2):
    let d = b[i] - w[i]
    [d*d, b[i]*b[i]]
  getDefaultComm().allReduce(addr r[0], 2)
  r[0]/r[1]

proc upload(h: HmcGpu; d: ptr UncheckedArray[float]; v: Field) =
  ## all sites of a color vector field, the layout of stagGpu
  gpuMemCpyToGpu(d, addr v[0], 6*h.s.n*sizeof(float))

proc smear(h: HmcGpu): float =
  ## the smeared links of gg.u in sgo and s; the seconds it took
  if h.smeared: return 0.0
  tic()
  h.hg.smear(h.coef, h.gg, h.sgo.u)
  h.sgo.fresh = false
  h.s.setLinks(h.sgo, h.sg)
  if h.mixed: h.ss.setLinks(h.sgo, h.sg)
  h.smeared = true
  result = getElapsedTime()
  toc()

proc smear(h: HmcGpu; a: GpuAction; key: string) =
  ## smear, timed as key of a
  let t = h.smear
  a.stats[key]["n"] += 1
  a.stats[key]["secs"] += t

proc solveStats(a: GpuAction; key: string; sp: SolverParams; secs: float) =
  a.stats[a.id & key]["n"] += 1
  a.stats[a.id & key]["secs"] += secs
  a.stats[a.id & key]["flops"] += sp.flops
  a.stats[a.id & key]["its"] += float sp.iterations
  a.stats[a.id & key]["r2"] += sp.r2.mean
  a.stats[a.id & key]["r2max"].maxeq sp.r2.max

type
  SysKind = enum skEE, skR, skL
  Sys = object
    ## a solve of solveAll: x from b with mass m, rhs a work vector
    x, b, rhs: ptr UncheckedArray[float]
    m: float
    kind: SysKind
    sp: SolverParams

proc solveAll(h: HmcGpu; sys: var seq[Sys]; mixed = false) =
  ## The solves of sys as those of hmcAction, together by solveM of several
  ## systems, M(m) = m + D with D = stag.D:
  ##   skEE  stagForceSolve: x = M(m)^-1 b for b_o = 0; its psi_e = x_e/m,
  ##         psi_o = -2 x_o, so the one link force of psi is -2/m that of x
  ##   skR   stag.solve for b_o = 0 (reconR): x = M(m)^-1 b
  ##   skL   stag.solve (reconL): the even sites of solveM for rhs_e =
  ##         d_e/m, d = M(m)^+ b, to the tolerance
  ##         0.99 r2req (|b_e|^2 + |b_o|^2) m^2/|d_e|^2, then x_o += b_o/m
  if sys.len == 0: return
  tic()
  let ne6 = 6*h.s.ne
  let no6 = 6*(h.s.n - h.s.ne)
  var xs, bs = newSeq[ptr UncheckedArray[float]](sys.len)
  var ms = newSeq[float](sys.len)
  var sps = newSeq[SolverParams](sys.len)
  for j in 0..<sys.len:
    let b = sys[j].b
    let m = sys[j].m
    xs[j] = sys[j].x
    ms[j] = m
    sps[j] = sys[j].sp
    sps[j].resetStats
    bs[j] = b
    if sys[j].kind == skL:
      let b2 = h.normEO(b)
      let d = sys[j].rhs
      h.s.applyM(d, b, -m)  # -m b_e + D_eo b_o
      let sc = -1.0/m
      gpuFor(i, ne6): d[i] *= sc
      sps[j].r2req = 0.99*sps[j].r2req*(b2[0] + b2[1])/h.normEO(d)[0]
      bs[j] = d
  if mixed: h.s.solveM(xs, bs, ms, sps, addr h.ss)
  else: h.s.solveM(xs, bs, ms, sps)
  let secs = getElapsedTime()/float(sys.len)
  for j in 0..<sys.len:
    let x = sys[j].x
    let b = sys[j].b
    let m = sys[j].m
    if sys[j].kind == skL:
      let c = 1.0/m
      gpuFor(i, no6): x[ne6 + i] += c*b[ne6 + i]
    var sp = sps[j]
    sp.r2req = sys[j].sp.r2req
    sp.r2.init h.resid(x, b, m)
    sp.flops = float((4*4*72 + 60)*h.s.ne*sp.iterations)
    sp.seconds = secs
    sys[j].sp = sp
  toc()

proc addForce(h: HmcGpu; a: GpuAction; secs: float) =
  ## f += ft and the stats of ft as addForce of hmcAction: the sums of |ft|^2
  ## and |ft|^4 over the links over 4V, and the largest |ft|^2 (summed over
  ## the ranks as rankSum does).  Thread t takes links t, t+T, ..., slot k
  ## of nSlot threads k, k+nSlot, ..., the host the slots.
  let n = h.gg.n
  let nl = 4*n
  let nt = (nl + 15) div 16
  let f = h.f
  let ft = h.ft
  let rd = h.red  # [3 nSlot] slots, then [3][nt] thread sums
  gpuFor(t, nt):
    var s2, s4, mx = 0.0
    for j in 0..<16:
      let i = t + j*nt
      if i < nl:
        let mu = i div n
        let k = i - mu*n
        let o = 18*mu*n + lo18(V, k)
        var q = 0.0
        forStatic e, 0, 17:
          let v = ft[o + e*V]
          f[o + e*V] += v
          q += v*v
        s2 += q
        s4 += q*q
        mx = max(mx, q)
    rd[3*nSlot + t] = s2
    rd[3*nSlot + nt + t] = s4
    rd[3*nSlot + 2*nt + t] = mx
  gpuFor(k, nSlot):
    var s2, s4, mx = 0.0
    var t = k
    while t < nt:
      s2 += rd[3*nSlot + t]
      s4 += rd[3*nSlot + nt + t]
      mx = max(mx, rd[3*nSlot + 2*nt + t])
      t += nSlot
    rd[k] = s2
    rd[nSlot + k] = s4
    rd[2*nSlot + k] = mx
  let hr = h.hred
  gpuMemCpyToCpu(hr, rd, 3*nSlot*sizeof(float))
  var fs: array[3, float]
  for k in 0..<nSlot:
    fs[0] += hr[k]
    fs[1] += hr[nSlot + k]
    fs[2] = max(fs[2], hr[2*nSlot + k])
  getDefaultComm().allReduce(addr fs[0], 3)
  let v = 1.0/float(4*h.gg.lo.physVol)
  a.stats[a.id & "F"]["n"] += 1
  a.stats[a.id & "F"]["secs"] += secs
  a.stats[a.id & "F"]["f2"] += v*fs[0]
  a.stats[a.id & "F"]["f4"] += v*fs[1]
  a.stats[a.id & "F"]["finf"].maxeq fs[2]

proc pullback(h: HmcGpu) =
  ## ft = TAH(F U^+) for the chain F of the one link force fl through the
  ## smearing, with the phases, boundary signs and -1 on the odd sites, as
  ## fermForce of hmcAction
  h.zeroLinks h.ft
  h.hg.force(h.coef, h.gg, h.sgo.u, h.fl, h.sg, h.s.ne, h.ft)

#[ heatbath, action, force ]#

proc heatbath(h: HmcGpu; level: GpuLevel) =
  ## the sources in the order of hmcAction, then the solves of the level
  var sys: seq[Sys]
  var owner: seq[GpuAction]
  for a in level.actions:
    case a.kind
    of gkGauge:
      a.stats[a.id & "A"] = baseStats0.newTable
      a.stats[a.id & "F"] = forceStats0.newTable
    of gkFermion, gkRatio, gkPV:
      a.stats[a.id & "F"] = forceStats0.newTable
      if a.kind == gkFermion:
        a.stats[a.id & "AS"] = solveStats0.newTable
        a.stats[a.id & "FS"] = solveStats0.newTable
      if a.kind == gkPV: a.stats[a.id & "RS"] = solveStats0.newTable
      if a.kind == gkRatio: discard h.smear
      else:
        a.stats["SS"] = baseStats0.newTable
        a.stats["SF"] = baseStats0.newTable
        h.smear(a, "SS")
      var psi = h.cpu.uc.u[0].l.ColorVector()
      let r = h.cpu.prng
      threads: psi.gaussian(r)
      h.upload(a.src, psi)
      case a.kind
      of gkFermion:  # phi = M(-m) psi on the even sites
        h.s.applyMfull(a.phi, a.src, -a.mass)
        h.zeroOdd a.phi
      of gkRatio:  # phi = M(-m_den)^-1 M(-m_num) psi on the even sites
        h.s.applyMfull(a.x, a.src, -a.mass)
        sys.add Sys(x: a.phi, b: a.x, rhs: a.rhs, m: -a.massDen, kind: skL, sp: a.spa)
        owner.add a
      else:  # phi = M(m)^-1 psi on the even sites
        sys.add Sys(x: a.phi, b: a.src, rhs: a.rhs, m: a.mass, kind: skL, sp: a.spa)
        owner.add a
  h.solveAll(sys)
  for j, a in owner:
    h.zeroOdd a.phi
    if a.kind == gkPV: a.solveStats("RS", sys[j].sp, sys[j].sp.seconds)

proc action(h: HmcGpu; level: GpuLevel): seq[float] =
  ## the action of each action of the level
  result.setLen level.actions.len
  var sys: seq[Sys]
  var owner: seq[int]
  for i, a in level.actions:
    case a.kind
    of gkGauge:
      tic()
      result[i] = h.gg.actionA(a.gc)
      a.stats[a.id & "A"]["n"] += 1
      a.stats[a.id & "A"]["secs"] += getElapsedTime()
      toc()
    of gkFermion:  # |M(-m)^-1 phi|^2/2
      discard h.smear
      sys.add Sys(x: a.x, b: a.phi, rhs: a.rhs, m: -a.mass, kind: skR, sp: a.spa)
      owner.add i
    of gkRatio:  # |M(-m_num)^-1 M(-m_den) phi|^2/2
      discard h.smear
      h.s.applyMfull(a.src, a.phi, -a.massDen)
      sys.add Sys(x: a.x, b: a.src, rhs: a.rhs, m: -a.mass, kind: skL, sp: a.spa)
      owner.add i
    of gkPV:  # |M(m) phi|^2/2
      discard h.smear
      h.s.applyMfull(a.x, a.phi, a.mass)
      result[i] = 0.5*h.s.norm2(a.x)
  h.solveAll(sys)
  for j, i in owner:
    let a = level.actions[i]
    result[i] = 0.5*h.s.norm2(a.x)
    if a.kind == gkFermion: a.solveStats("AS", sys[j].sp, sys[j].sp.seconds)

proc force(h: HmcGpu; level: GpuLevel; dtau: float) =
  ## f += the forces of the level with dtau, as the force procs of hmcAction:
  ##   fermion  0.25 dtau F(psi), psi_e = 4 A^-1 phi_e, psi_o = D_oe psi_e
  ##   ratio    0.25 dtau (m_den^2 - m_num^2) F(psi), psi as fermion with m_num
  ##   PV       -0.25 dtau F(psi), psi_e = phi_e, psi_o = D_oe phi_e
  ## for F the force of the one link force psi(x) psi(x+mu)^+.
  var sys: seq[Sys]
  var owner: seq[int]
  for i, a in level.actions:
    if a.kind == gkFermion: h.smear(a, "SF")
    elif a.kind in {gkRatio, gkPV}: discard h.smear
    if a.kind in {gkFermion, gkRatio}:
      sys.add Sys(x: a.x, b: a.phi, rhs: a.rhs, m: a.mass, kind: skEE, sp: a.spf)
      owner.add i
  h.solveAll(sys, h.mixed)
  var fsecs = newSeq[float](level.actions.len)
  for j, i in owner:
    let a = level.actions[i]
    fsecs[i] = sys[j].sp.seconds
    if a.kind == gkFermion: a.solveStats("FS", sys[j].sp, fsecs[i])
  let each = h.forceStats == 2 or h.sample
  # the PV vectors psi_e = phi_e, psi_o = D_oe phi_e for the hop D of dslash (twice
  # applyMfull(psi, phi, 0) on the odd sites), together
  var pv: seq[int]
  for i, a in level.actions:
    if a.kind == gkPV: pv.add i
  if pv.len > 0:
    tic()
    let ne6 = 6*h.s.ne
    for i in pv:
      let y = level.actions[i].x
      let ph = level.actions[i].phi
      gpuFor(k, ne6): y[k] = ph[k]
    h.s.hopOE(pv.mapIt(level.actions[it].x), pv.mapIt(level.actions[it].phi), 1.0)
    let t = getElapsedTime()/float(pv.len)
    for i in pv: fsecs[i] += t
    toc("PV vectors")
  var fx: seq[ptr UncheckedArray[float]]  # the one link forces summed for one pullback
  var ft: seq[float]
  var fa: seq[int]
  for i, a in level.actions:
    tic()
    case a.kind
    of gkGauge:
      h.zeroLinks h.ft
      h.gg.forceA(a.gc, h.ft, -dtau)
      h.addForce(a, getElapsedTime())
    of gkFermion, gkRatio, gkPV:
      var t: float
      if a.kind == gkPV: t = -0.25*dtau
      else:
        let sc = if a.kind == gkFermion: 0.25*dtau
                 else: 0.25*dtau*(a.massDen*a.massDen - a.mass*a.mass)
        t = -2.0*sc/a.mass
      if each:
        h.zeroLinks h.fl
        h.s.outerM(h.fl, a.x, t)
        h.pullback
        h.addForce(a, getElapsedTime() + fsecs[i])
      else:
        fx.add a.x
        ft.add t
        fa.add i
    toc()
  if fx.len > 0:  # f += the force of the sum of the one link forces, the outer products together
    tic()
    h.zeroLinks h.fl
    h.s.outerM(h.fl, fx, ft)
    let t = getElapsedTime()/float(fx.len)
    for i in fa:  # the time of the action, its stats from the sampled forces
      let a = level.actions[i]
      a.stats[a.id & "F"]["secs"] += t + fsecs[i]
    toc("outer")
    h.hg.force(h.coef, h.gg, h.sgo.u, h.fl, h.sg, h.s.ne, h.f)
    h.hmcStats["SF"]["n"] += 1
    h.hmcStats["SF"]["secs"] += getElapsedTime()
    toc()
  if level.actions.anyIt(it.kind != gkGauge): h.sample = false

#[ the molecular dynamics ]#

proc kinetic(h: HmcGpu): float =
  0.5*h.gg.norm2(h.p) - 16.0*float(h.gg.lo.physVol)

proc action*(h: HmcGpu): float =
  for level in h.levels:
    for v in h.action(level): result += v

proc hamiltonian*(h: HmcGpu): float =
  let t = h.kinetic
  let v = h.action
  result = t + v
  if h.atEnd: echo fmt"Ending H: {result} T: {t} V: {v}"
  else: echo fmt"Beginning H: {result} T: {t} V: {v}"

proc integrator*(h: HmcGpu): Integrator =
  let nlevels = h.levels.len
  proc mdt(dtau: float) =
    tic("mdt")
    h.gg.expUpdate(h.p, dtau)
    h.smeared = false
    h.hmcStats["GU"]["n"] += 1
    h.hmcStats["GU"]["secs"] += getElapsedTime()
    toc("end")
  proc mdv(dtau: openArray[float]) =
    h.zeroLinks h.f
    for i in 0..<nlevels:
      if dtau[i] != 0.0: h.force(h.levels[i], dtau[i])
    let p = h.p
    let f = h.f
    gpuFor(i, 4*18*h.gg.n): p[i] -= f[i]
  let (V, T) = newIntegratorPair(mdv, mdt)
  result = T
  for i in 0..<nlevels:
    result = h.levels[i].integrator(T = result, V = V[i], steps = h.levels[i].multiplier)

proc toHost(h: HmcGpu) = h.gg.download(h.cpu.uc.u, h.gg.u)

#[ "virtual" MetropolisRoot procedures ]#

proc getH*(h: HmcGpu): float = h.hamiltonian

proc start*(h: HmcGpu) =
  h.atEnd = false
  h.hmcStats["GU"] = baseStats0.newTable
  h.hmcStats["SF"] = baseStats0.newTable
  h.sample = h.forceStats == 1
  if h.hostNew:
    h.gg.upload(h.gg.u, h.cpu.uc.u)
    h.smeared = false
    h.hostNew = false
  h.cpu.heatbathProc()
  h.gg.upload(h.p, h.cpu.p)
  for level in h.levels: h.heatbath(level)
  h.gg.copy(h.bu, h.gg.u)

proc generate*(h: HmcGpu) =
  var integ = h.integrator
  integ.evolve h.tau
  integ.finish
  h.atEnd = true

proc checkReverse*(h: HmcGpu): bool =
  h.revCheckFreq > 0 and (h.nUpdates mod h.revCheckFreq == 0)

proc generateReverse*(h: HmcGpu) =
  h.gg.copy(h.g1, h.gg.u)
  h.gg.copy(h.p1, h.p)
  let p = h.p
  gpuFor(i, 4*18*h.gg.n): p[i] = -p[i]
  var integ = h.integrator
  integ.evolve h.tau
  integ.finish
  let t = h.kinetic
  let v = h.action
  h.hReverse = t + v
  echo fmt"Reverse H: {h.hReverse} T: {t} V: {v}"
  h.gg.copy(h.gg.u, h.g1)
  h.gg.copy(h.p, h.p1)
  h.gg.fresh = false
  h.smeared = false

proc finishReverse*(h: HmcGpu) =
  echo "REVERSE ",
    "  dH (from hNew): ", h.hReverse - h.hNew,
    "  dH0 (from hOld): ", h.hReverse - h.hOld

proc globalRand*(h: HmcGpu): float =
  if h.forceAccept: 0.0 else: h.cpu.globalRandProc()

proc accept*(h: HmcGpu) =
  ## reunit on the host as hmcAction, back to the device
  h.toHost
  h.cpu.reunit
  h.gg.upload(h.gg.u, h.cpu.uc.u)
  h.smeared = false

proc reject*(h: HmcGpu) =
  h.gg.copy(h.gg.u, h.bu)
  h.gg.fresh = false
  h.smeared = false
  h.toHost

proc run*(h: HmcGpu; forceAccept = false) =
  tic("HmcGpu:run")
  let nup = h.nUpdates + 1
  echo &"== Begin HMC update {nup} =========="
  h.forceAccept = forceAccept
  metropolis.update(h)
  let dt = getElapsedTime()
  h.secs += dt
  var parts = @[h.hmcStats]
  for level in h.levels:
    for a in level.actions: parts.add a.stats
  showStats(parts, dt, h.secs, nup, h.verbosity)
  echo "===================================="
  toc("end")

#[ host links: start, files, measurements ]#

proc cold*(h: HmcGpu) =
  h.cpu.cold
  h.hostNew = true

proc hot*(h: HmcGpu) =
  h.cpu.hot
  h.hostNew = true

proc read*(h: HmcGpu; readParallelRNG = false; readSerialRNG = false;
           gaugeFilename, parallelRNGFilename, serialRNGFilename: string) =
  h.cpu.read(readParallelRNG = readParallelRNG, readSerialRNG = readSerialRNG,
             gaugeFilename = gaugeFilename, parallelRNGFilename = parallelRNGFilename,
             serialRNGFilename = serialRNGFilename)
  h.hostNew = true

proc write*(h: HmcGpu; writeParallelRNG = false; writeSerialRNG = false;
            gaugeFilename, parallelRNGFilename, serialRNGFilename: string) =
  h.cpu.write(writeParallelRNG = writeParallelRNG, writeSerialRNG = writeSerialRNG,
              gaugeFilename = gaugeFilename, parallelRNGFilename = parallelRNGFilename,
              serialRNGFilename = serialRNGFilename)

proc measurePlaquette*(h: HmcGpu) = h.cpu.measurePlaquette
proc measurePolyakovLoop*(h: HmcGpu) = h.cpu.measurePolyakovLoop
proc measurePlaquetteS4*(h: HmcGpu) = h.cpu.measurePlaquetteS4
