## Staggered HMC with Hasenbusch mass ratios as examples/staghmc_sh, with the
## links, momenta and fermion fields on the GPU and every step in kernels:
## the solves and one link forces (physics/stagGpu), the gauge action and
## force (gauge/gaugeGpu), the HYP smearing and its force (gauge/hypGpu).
## The random numbers come from the generators of staghmc_sh (rng/rngGpu)
## drawn in its order, so
##   staghmc_sh -gintalg:2MN -fintalg:2MN -pbpreps:0 -revCheckFreq:0 -savefreq:0
## gives the same trajectories as this program with the same lattice, seed,
## masses, steps, tau and HYP coefficients.  The smeared links are U(3), so
## with smearing the solvers keep the determinant of each link with rows 0
## and 1, 14 reals, and rebuild row 2 = det U conj(row 0 x row 1).
##   -rg: ranks per dimension, lanes as in bestagcg; -hot:1 random links
##   -mass, -hmasses: the mass and the Hasenbusch masses, heavier
##   -gsteps, -fsteps, -hfsteps: steps of the gauge, mass and Hasenbusch terms
##   -arsq, -frsq, -hfrsq: CG tolerances of the action and force solves
##   -alpha1, -alpha2, -alpha3: HYP coefficients, 0.4, 0.5, 0.5 as staghmc_sh
##   -mixed:1 mixed precision solves; -reals: reals per link of the solvers, 12 for
##   SU(3) links, 14 by default with smearing (rows 0 and 1 and the determinant), 18
##   -check:1 first compares the smeared links and the smearing force with
##   gauge/hypsmear2 on random links
import qex, gauge, gauge/[hypsmear, hypsmear2, gaugeGpu, hypGpu], physics/[qcdTypes, stagSolve, stagGpu]
import backend/accel, rng/rngGpu
import mdevolve
import times, macros, sequtils

qexinit()

let
  lat = intSeqParam("lat", @[8,8,8,8])
  rg = intSeqParam("rg")
  beta = floatParam("beta", 6.0)
  adjFac = floatParam("adjFac", -0.25)
  tau = floatParam("tau", 2.0)
  gsteps = intParam("gsteps", 4)
  fsteps = intParam("fsteps", 4)
  trajs = intParam("trajs", 10)
  seed = intParam("seed", int(1000*epochTime())).uint64
  mass = floatParam("mass", 0.1)
  hmasses = floatSeqParam("hmasses", @[0.2, 0.4])
  hfsteps = intSeqParam("hfsteps", repeat(fsteps, hmasses.len))
  arsq = floatParam("arsq", 1e-20)
  frsq = floatParam("frsq", 1e-12)
  hfrsq = floatSeqParam("hfrsq", repeat(frsq, hmasses.len))
  nt = hmasses.len + 1  # fermion terms: the ratios and the heaviest determinant
  alpha1 = floatParam("alpha1", 0.4)  # HYP smearing coefficients, as staghmc_sh
  alpha2 = floatParam("alpha2", 0.5)
  alpha3 = floatParam("alpha3", 0.5)
  smear = alpha1 != 0 or alpha2 != 0 or alpha3 != 0
  coef = HypCoefs(alpha1: alpha1, alpha2: alpha2, alpha3: alpha3)

macro echoparam(x: typed): untyped =
  let n = x.repr
  result = quote do:
    echo `n`, ": ", `x`

echoparam(beta)
echoparam(adjFac)
echoparam(tau)
echoparam(gsteps)
echoparam(fsteps)
echoparam(hfsteps)
echoparam(trajs)
echoparam(seed)
echoparam(mass)
echoparam(hmasses)
echoparam(arsq)
echoparam(frsq)
echoparam(hfrsq)
echo "smear = ", coef
if hfsteps.len != hmasses.len or hfrsq.len != hmasses.len:
  qexError "hfsteps and hfrsq need one entry per Hasenbusch mass"

let
  gc = GaugeActionCoeffs(plaq: beta, adjplaq: beta*adjFac)
  lo = if rg.len > 0: newLayout(lat, VLEN, rg, innerGeom(lat, rg, static(VLEN))) else: lat.newLayout
  vol = lo.physVol

var r = lo.newRNGFieldV(RngMilc6, seed)
var R: RngMilc6  # global RNG
R.seed(seed, 987654321)

var g = lo.newgauge
if intParam("hot", 0) != 0: g.random r
else: g.unit
var gs = lo.newgauge  # unit links with the phases of the fermions
gs.unit
threads:
  gs.setBC
  threadBarrier()
  gs.stagPhase

var gg = newGpuGauge(lo)
# the smeared links are U(3), projectU without a determinant condition, so
# the solvers keep their determinant; unsmeared links are SU(3) and 12 do
let reals = intParam("reals", if smear: 14 else: 12)
if smear and reals < 14: qexError "the smeared links are U(3): -reals:14 or 18"
var s = newStagGpu(gs, float64, reals)
var ss: StagGpu[VLEN,float32]
let ssp = if intParam("mixed", 0) != 0: (ss = newStagGpu(gs, float32, reals); addr ss) else: nil
echo "GPU links: ", s.nl, " reals", if s.lb == nil: ", forward only" else: ""
let sg = stagSigns(gs)
let mom = gg.newLinks
let g0 = gg.newLinks
# with smearing the solvers take the smeared links sgo
var sgo: GpuGauge[VLEN]
var hg: HypGpu[VLEN]
if smear:
  sgo = newGpuGauge(lo)
  hg = newHypGpu(lo)
let sgg = if smear: addr sgo else: addr gg  # the links of the solver
let fF = if smear: gg.newLinks else: nil  # sum of the one link forces of a step
let n6 = 6*s.n
proc newVec(): ptr UncheckedArray[float] = cast[ptr UncheckedArray[float]](gpuMalloc(n6*sizeof(float)))
var phi = newSeq[ptr UncheckedArray[float]](nt)
for i in 0..<nt: phi[i] = newVec()
let psi = newVec()
let ftmp = newVec()
let x = newVec()

proc mt(i: int): float =
  ## the lighter mass of term i
  if i == 0: mass else: hmasses[i-1]
proc mh(i: int): float =
  ## the heavier mass of term i
  if i < hmasses.len: hmasses[i] else: mt(i)
proc fscale(i: int, t: float): float =
  ## as staghmc_sh: the force of term i is d/dU of |M(m_i)^-1 M(h_i) phi_i|^2/2
  if hmasses.len == 0: 0.5*t/mass
  elif i == 0: 0.5*t*(hmasses[0]*hmasses[0] - mass*mass)/mass
  elif i < hmasses.len: 0.5*t*(hmasses[i]*hmasses[i] - hmasses[i-1]*hmasses[i-1])/hmasses[i-1]
  else: 0.5*t/hmasses[i-1]

var spa = initSolverParams()
spa.r2req = arsq
spa.maxits = 1000000
spa.verbosity = 0
var spf = newSeq[SolverParams](nt)
for i in 0..<nt:
  spf[i] = initSolverParams()
  spf[i].r2req = if i == 0: frsq else: hfrsq[i-1]
  spf[i].maxits = 1000000
  spf[i].verbosity = 0

proc zeroOdd(v: ptr UncheckedArray[float]) =
  let o = 6*s.ne
  gpuFor(i, n6 - o): v[o + i] = 0.0

var linksSet = false
proc setLinks() =
  if linksSet: return
  if smear:
    hg.smear(coef, gg, sgo.u)
    sgo.fresh = false
  s.setLinks(sgg[], sg)
  if ssp != nil: ssp[].setLinks(sgg[], sg)
  linksSet = true

proc faction(): seq[float] =
  ## |M(m_i)^-1 M(h_i) phi_i|^2/2 per term, the last |M(h)^-1 phi|^2/2
  for i in 0..<nt:
    tic("faction")
    if i != nt-1:
      s.applyMfull(ftmp, phi[i], mh(i))
      s.solveM(psi, ftmp, mt(i), spa, ssp, full = true)
    else:
      s.solveM(psi, phi[i], mh(i), spa, ssp)
    toc("solve")
    result.add 0.5*s.norm2(psi)
    toc("norm")

proc mdt(t: float) =
  tic()
  gg.expUpdate(mom, t)
  linksSet = false
  toc("mdt")

proc mdv(t: float) =
  tic()
  gg.forceA(gc, mom, t)
  toc("mdv")

proc mdvf(i: int, t: float) =
  tic()
  setLinks()
  toc("fforce links")
  s.solveM(x, phi[i], mt(i), spf[i], ssp)
  toc("fforce solve")
  if smear: s.outerM(fF, x, fscale(i, t))
  else: s.forceM(sgg[], sg, mom, x, -fscale(i, t))
  toc("mdvf")

proc pullback() =
  ## p += TAH(smearing pullback of the phased force sum), as
  ## smearedOneLinkForce of staghmc_sh
  tic("pullback")
  hg.force(coef, gg, sgo.u, fF, sg, s.ne, mom)
  toc("force")

proc mdvAll(t: openarray[float]) =
  if t[0] != 0: mdv t[0]
  var any = false
  for i in 0..<nt:
    if t[i+1] != 0:
      if smear and not any:
        let f = fF
        gpuFor(i, 4*18*s.n): f[i] = 0.0
      any = true
      mdvf(i, t[i+1])
  if smear and any: pullback()

let (VAll, T) = newIntegratorPair(mdvAll, mdt)
let H = newParallelEvolution(mkOmelyan2MN(steps = gsteps, V = VAll[0], T = T))
H.add mkOmelyan2MN(steps = fsteps, V = VAll[1], T = T)
for i in 0..<hmasses.len:
  H.add mkOmelyan2MN(steps = hfsteps[i], V = VAll[2+i], T = T)

proc rel(a, b: auto): float =
  ## |a-b|^2/|b|^2 summed over the fields
  var d, n = 0.0
  threads:
    var dt, nt = 0.0
    for mu in 0..<a.len:
      dt += norm2(a[mu] - b[mu])
      nt += b[mu].norm2
    threadMaster:
      d = dt
      n = nt
  d/n

if intParam("check", 0) != 0 and smear:
  # random links and chain: the smeared links and p = TAH(F U^+) for the
  # force F of the phased chain against gauge/hypsmear2 and contractProjectTAH
  var rc = lo.newRNGField(RngMilc6, seed + 1)
  var gr = lo.newgauge
  var ch = lo.newgauge
  var cp = lo.newgauge
  var fl = lo.newgauge
  var f = lo.newgauge
  var fg = lo.newgauge
  var z = lo.newgauge
  gr.random rc
  ch.gaussian rc
  threads:
    for mu in 0..<cp.len:
      cp[mu] := ch[mu]
      z[mu] := 0
    cp.setBC
    threadBarrier()
    cp.stagPhase
    threadBarrier()
    for mu in 0..<cp.len:
      for i in cp[mu].odd:
        cp[mu][i] *= -1
  let htc = newHypTemps(gr)
  htc.smear(coef, fl)
  htc.force(coef, f, cp)
  threads: contractProjectTAH(f, f, gr)
  gg.upload(gg.u, gr)
  hg.smear(coef, gg, sgo.u)
  gg.download(fg, sgo.u)
  echo "check smear  |GPU-CPU|^2/|CPU|^2: ", rel(fg, fl)
  gg.upload(fF, ch)
  gg.upload(mom, z)
  hg.force(coef, gg, sgo.u, fF, sg, s.ne, mom)
  gg.download(fg, mom)
  echo "check force  |GPU-CPU|^2/|CPU|^2: ", rel(fg, f)

gg.upload(gg.u, g)
echo "plaq: ", gg.plaq
echo "actionA: ", gg.actionA(gc)

let rgpu = newRngGpu(r)
for n in 1..trajs:
  tic()
  let t0 = epochTime()
  rgpu.randomTAH mom
  gg.copy(g0, gg.u)
  setLinks()
  toc("random, links")
  # phi_i = M(-h_i)^-1 M(-m_i) psi_i on the even sites, the last M(-h) psi,
  # the minus signs as staghmc_sh (bsm.lua convention)
  for i in 0..<nt:
    rgpu.gaussian(psi, 6)
    if i != nt-1:
      s.applyMfull(ftmp, psi, -mt(i))
      s.solveM(phi[i], ftmp, -mh(i), spa, ssp, full = true)
    else:
      s.applyMfull(phi[i], psi, -mh(i))
    zeroOdd phi[i]
  toc("init")
  let fa0 = faction()
  toc("fa solve 1")
  let
    p2 = gg.norm2(mom)
    ga0 = gg.actionA(gc)
    k0 = 0.5*p2 - (16*vol).float
    h0 = ga0 + fa0.sum + k0
  toc("init gauge action")
  echo "Begin H: ",h0,"  Sg: ",ga0,"  Sf: ",fa0,"  T: ",k0

  H.evolve tau
  H.finish
  toc("evolve")

  let p2e = gg.norm2(mom)
  setLinks()
  let fa1 = faction()
  toc("fa solve 2")
  let
    ga1 = gg.actionA(gc)
    k1 = 0.5*p2e - (16*vol).float
    h1 = ga1 + fa1.sum + k1
  toc("final gauge action")
  echo "End H: ",h1,"  Sg: ",ga1,"  Sf: ",fa1,"  T: ",k1

  let
    dH = h1 - h0
    acc = exp(-dH)
    accr = R.uniform
  if accr <= acc:
    echo "ACCEPT:  dH: ",dH,"  exp(-dH): ",acc,"  r: ",accr
  else:
    echo "REJECT:  dH: ",dH,"  exp(-dH): ",acc,"  r: ",accr
    gg.copy(gg.u, g0)
    gg.fresh = false
    linksSet = false
  echo "plaq: ", gg.plaq, "  traj secs: ", epochTime() - t0

echoTimers()
qexfinalize()
