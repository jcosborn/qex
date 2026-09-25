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
## sums run in a fixed order (gpuSum and the CG here), so a trajectory
## repeats bit for bit.  One boundary condition holds for all fermion
## actions.

import qex
import gauge
import physics/[qcdTypes, stagD, stagSolve, stagGpu]
import gauge/[gaugeGpu, hypGpu, hypsmear]
import hmc/[hmcAction, metropolis]
import algorithms/[integrator]
import backend/accel
import base/metaUtils
import std/[algorithm, strformat, strutils, tables]

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
    sg*: ptr UncheckedArray[float]  ## signs of the phases and boundary conditions, stagSigns
    p*, f*, bu*, ft*, fl*, g1*, p1*: ptr UncheckedArray[float]  ## momenta, force, saved links, force of an action, one link force, reverse check copies
    w*, cr*, cp*, cap*, ct*: ptr UncheckedArray[float]  ## work vectors, 6n reals
    red*, hred*: ptr UncheckedArray[float]  ## force stats partial sums and maxima, device and pinned host
    smeared*: bool  ## sgo and s hold the smeared links of gg.u
    hostNew*: bool  ## the host links are newer than gg.u
    hmcStats*: Table[string, ActionStats]
    secs*: float
    revCheckFreq*: int
    forceAccept*: bool
    atEnd: bool

const nSlot = 128  # atomic slots per sum of the force stats
const V = VLEN

proc newVec(h: HmcGpu): ptr UncheckedArray[float] =
  cast[ptr UncheckedArray[float]](gpuMalloc(6*h.s.n*sizeof(float)))

proc newHmcGpu*[U,R,S](uc: GaugeConfiguration[U]; srng: S; prng: R; tau: float; revCheckFreq = 0;
                       bc: openArray[bool] = [true, true, true, false]; reals = 18): HmcGpu[U,R,S] =
  ## the HMC of hmcAction on the GPU, the fermion boundary conditions bc
  ## (periodic when true) for all fermion actions; reals per link of the
  ## solvers, 18 or 14 (rows 0 and 1 and the determinant)
  new(result)
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
  result.gg = newGpuGauge(lo)
  result.sgo = newGpuGauge(lo)
  result.hg = newHypGpu(lo)
  for d in [addr result.p, addr result.f, addr result.bu, addr result.ft, addr result.fl,
            addr result.g1, addr result.p1]:
    d[] = result.gg.newLinks
  for d in [addr result.w, addr result.cr, addr result.cp, addr result.cap, addr result.ct]:
    d[] = result.newVec
  let nt = (4*result.gg.n + 15) div 16
  result.red = cast[ptr UncheckedArray[float]](gpuMalloc((2*nSlot + nt)*sizeof(float)))
  result.hred = cast[ptr UncheckedArray[float]](gpuMallocHost((2*nSlot + nt)*sizeof(float)))
  let rd = result.red
  gpuFor(i, 2*nSlot): rd[i] = 0.0
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

proc dot(h: HmcGpu; a, b: ptr UncheckedArray[float]): float =
  ## global a_e.b_e
  result = gpuSum(i, 6*h.s.ne, 1, [a[i]*b[i]])[0]
  getDefaultComm().allReduce(result)

proc cg(h: HmcGpu; y, d: ptr UncheckedArray[float]; m: float; sp: var SolverParams) =
  ## y_e = A^-1 d_e, A = 4 (m^2 - D_eo D_oe), by the iteration of cgSolve
  ## from y = 0 until |r|^2 <= r2req |d_e|^2
  let n6 = 6*h.s.ne
  let r = h.cr
  let p = h.cp
  let ap = h.cap
  gpuFor(i, n6):
    y[i] = 0.0
    r[i] = d[i]
  var r2 = h.dot(r, r)
  let r2stop = sp.r2req*r2
  var r2o = 0.0
  var itn = 0
  while itn < sp.maxits and r2 > r2stop:
    if itn == 0:
      gpuFor(i, n6): p[i] = r[i]
    else:
      let beta = r2/r2o
      gpuFor(i, n6): p[i] = r[i] + beta*p[i]
    inc itn
    h.s.applyD2ee(ap, p, h.ct, m*m)
    let alpha = r2/h.dot(p, ap)
    var q2 = gpuSum(i, n6, 1):
      y[i] += alpha*p[i]
      let q = r[i] - alpha*ap[i]
      r[i] = q
      [q*q]
    getDefaultComm().allReduce(q2[0])
    r2o = r2
    r2 = q2[0]
  sp.iterations += itn

proc solveAll(h: HmcGpu; sys: var seq[Sys]) =
  ## The solves of sys as those of hmcAction, with M(m) = m + D, D = stag.D,
  ## and y_e = A^-1 d_e by cg:
  ##   skEE  stagForceSolve: y for d = b, x_e = 4 y_e, x_o = 2 D_oe x_e
  ##   skR   stag.solve for b_o = 0 (reconR): y for d = b, x_e = 4m y_e,
  ##         x_o = -D_oe x_e/m
  ##   skL   stag.solve (reconL): y for d = M(m)^+ b to the tolerance
  ##         0.99 r2req (|b_e|^2 + |b_o|^2) m^2/|d_e|^2, x_e = 4 y_e,
  ##         x_o = (b_o - D_oe x_e)/m
  let ne6 = 6*h.s.ne
  let no6 = 6*(h.s.n - h.s.ne)
  for j in 0..<sys.len:
    tic()
    let x = sys[j].x
    let b = sys[j].b
    let m = sys[j].m
    let kind = sys[j].kind
    var sp = sys[j].sp
    sp.resetStats
    var d = b
    if kind == skL:
      let b2 = h.normEO(b)
      d = sys[j].rhs
      let dd = d
      h.s.applyM(dd, b, -m)  # -m b_e + D_eo b_o
      gpuFor(i, ne6): dd[i] = -dd[i]
      sp.r2req = 0.99*sp.r2req*(b2[0] + b2[1])*m*m/h.normEO(dd)[0]
    h.cg(x, d, m, sp)
    var rr = 0.0
    if kind == skEE:  # |A y - d|^2/|d|^2 as stagForceSolve
      let ap = h.cap
      h.s.applyD2ee(ap, x, h.ct, m*m)
      var q = gpuSum(i, ne6, 2):
        let e = ap[i] - d[i]
        [e*e, d[i]*d[i]]
      getDefaultComm().allReduce(addr q[0], 2)
      rr = q[0]/q[1]
    let sc = if kind == skR: 4.0*m else: 4.0
    gpuFor(i, ne6): x[i] *= sc
    h.zeroOdd x
    let w = h.w
    h.s.applyMfull(w, x, 0.0)  # w_o = D_oe x_e
    case kind
    of skEE:
      gpuFor(i, no6): x[ne6 + i] = 2.0*w[ne6 + i]
    of skR:
      let c = -1.0/m
      gpuFor(i, no6): x[ne6 + i] = c*w[ne6 + i]
    of skL:
      let c = 1.0/m
      gpuFor(i, no6): x[ne6 + i] = c*(b[ne6 + i] - w[ne6 + i])
    if kind != skEE: rr = h.resid(x, b, m)
    sp.r2.init rr
    sp.flops = float((4*4*72 + 60)*h.s.ne*sp.iterations)
    sp.seconds = getElapsedTime()
    sys[j].sp = sp
    toc()

proc addForce(h: HmcGpu; a: GpuAction; secs: float) =
  ## f += ft and the stats of ft as addForce of hmcAction: the sums of |ft|^2
  ## and |ft|^4 over the links over 4V, and the largest |ft|^2 (summed over
  ## the ranks as rankSum does)
  let n = h.gg.n
  let nl = 4*n
  let nt = (nl + 15) div 16
  let f = h.f
  let ft = h.ft
  let rd = h.red
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
    gpuAtomicAdd(rd, t mod nSlot, s2)
    gpuAtomicAdd(rd, nSlot + t mod nSlot, s4)
    rd[2*nSlot + t] = mx
  let hr = h.hred
  gpuMemCpyToCpu(hr, rd, (2*nSlot + nt)*sizeof(float))
  gpuFor(i, 2*nSlot): rd[i] = 0.0
  var fs: array[3, float]
  for k in 0..<nSlot:
    fs[0] += hr[k]
    fs[1] += hr[nSlot + k]
  for t in 0..<nt: fs[2] = max(fs[2], hr[2*nSlot + t])
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
  h.solveAll(sys)
  var fsecs = newSeq[float](level.actions.len)
  for j, i in owner:
    let a = level.actions[i]
    fsecs[i] = sys[j].sp.seconds
    if a.kind == gkFermion: a.solveStats("FS", sys[j].sp, fsecs[i])
  for i, a in level.actions:
    tic()
    case a.kind
    of gkGauge:
      h.zeroLinks h.ft
      h.gg.forceA(a.gc, h.ft, -dtau)
    of gkFermion, gkRatio:
      let sc = if a.kind == gkFermion: 0.25*dtau
               else: 0.25*dtau*(a.massDen*a.massDen - a.mass*a.mass)
      h.zeroLinks h.fl
      h.s.outerM(h.fl, a.x, sc)
      h.pullback
    of gkPV:  # psi_e = phi_e, psi_o = 2 D_oe phi_e
      let y = a.x
      let ph = a.phi
      let ne6 = 6*h.s.ne
      h.s.applyMfull(y, ph, 0.0)
      gpuFor(k, 6*h.s.n):
        y[k] = if k < ne6: ph[k] else: 2.0*y[k]
      h.zeroLinks h.fl
      h.s.outerM(h.fl, y, -0.25*dtau)
      h.pullback
    h.addForce(a, getElapsedTime() + fsecs[i])
    toc()

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
