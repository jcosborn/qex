## Staggered HMC as examples/staghmc, with the links, momenta and fermion
## fields resident on the GPU and every step of the trajectories in GPU
## kernels.  The momenta and the fermion source come from the same host
## random number fields as in staghmc and go to the device once per
## trajectory, so the two programs give the same trajectories.
##   -rg: ranks per dimension, with the lanes as in bestagcg
##   -hot:1 starts from random links instead of unit ones, as staghmc -hot:1
##   -mixed:1 solves in mixed precision, restarting single precision CGs
##   -recon:0 keeps 18 reals per link in the solvers, which otherwise take 12
##   and so need SU(3) links (the random ones of -hot:1 are so to 1e-10)
##   -check:1 first compares the gauge action, gauge force, link update and
##   fermion force with the CPU ones on a random gauge field
import qex, gauge, physics/[qcdTypes, stagSolve, stagGpu], gauge/gaugeGpu
import backend/accel
import mdevolve
import times, macros

qexinit()

let
  lat = intSeqParam("lat", @[8,8,8,8])
  rg = intSeqParam("rg")
  beta = floatParam("beta", 6.0)
  adjFac = floatParam("adjFac", -0.25)
  tau = floatParam("tau", 2.0)
  gsteps = intParam("gsteps", 32)
  fsteps = intParam("fsteps", 32)
  trajs = intParam("trajs", 10)
  seed = intParam("seed", int(1000*epochTime())).uint64
  mass = floatParam("mass", 0.1)
  arsq = floatParam("arsq", 1e-20)
  frsq = floatParam("frsq", 1e-12)

macro echoparam(x: typed): untyped =
  let n = x.repr
  result = quote do:
    echo `n`, ": ", `x`

echoparam(beta)
echoparam(adjFac)
echoparam(tau)
echoparam(gsteps)
echoparam(fsteps)
echoparam(trajs)
echoparam(seed)
echoparam(mass)
echoparam(arsq)
echoparam(frsq)

let
  gc = GaugeActionCoeffs(plaq: beta, adjplaq: beta*adjFac)
  lo = if rg.len > 0: newLayout(lat, VLEN, rg, innerGeom(lat, rg, static(VLEN))) else: lat.newLayout
  vol = lo.physVol

var r = lo.newRNGField(RngMilc6, seed)
var R: RngMilc6  # global RNG
R.seed(seed, 987654321)

var g = lo.newgauge
if intParam("hot", 0) != 0: g.random r  # random links instead of unit ones
else: g.unit
var p = lo.newgauge
var psi = lo.ColorVector()
var gs = lo.newgauge  # unit links with the phases of the fermions
gs.unit
threads:
  gs.setBC
  threadBarrier()
  gs.stagPhase

var gg = newGpuGauge(lo)
let recon = intParam("recon", 1) != 0
var s = newStagGpu(gs, float64, recon)
var ss: StagGpu[VLEN,float32]  # single precision solver with -mixed:1
let ssp = if intParam("mixed", 0) != 0: (ss = newStagGpu(gs, float32, recon); addr ss) else: nil
echo "GPU links: ", s.nl, " reals", if s.lb == nil: ", forward only" else: ""
let sg = stagSigns(gs)
let mom = gg.newLinks  # momenta
let g0 = gg.newLinks  # links at the start of the trajectory
let phi = cast[ptr UncheckedArray[float]](gpuMalloc(6*s.n*sizeof(float)))
let x = cast[ptr UncheckedArray[float]](gpuMalloc(6*s.n*sizeof(float)))

var spa = initSolverParams()
spa.r2req = arsq
spa.maxits = 10000
spa.verbosity = 0
var spf = initSolverParams()
spf.r2req = frsq
spf.maxits = 10000
spf.verbosity = 0

proc olf(f: var auto, v1: auto, v2: auto) =
  var t {.noInit.}: type(f)
  for i in 0..<v1.len:
    for j in 0..<v2.len:
      t[i,j] := v1[i] * v2[j].adj
  projectTAH(f, t)

proc oneLinkForce(f: auto, p: auto, g: auto) =
  let t = newTransporters(g, p, 1)
  for mu in 0..<g.len:
    discard t[mu] ^* p
  for mu in 0..<g.len:
    for i in f[mu]:
      olf(f[mu][i], p[i], t[mu].field[i])
    for i in f[mu].odd:
      f[mu][i] *= -1

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

if intParam("check", 0) != 0:
  var rc = lo.newRNGField(RngMilc6, seed + 1)
  var gr = lo.newgauge
  var f = lo.newgauge
  var fg = lo.newgauge
  var z = lo.newgauge
  gr.random rc
  threads:
    gr.projectSU  # SU(3) to rounding, as the solvers with links of 12 reals need
    for mu in 0..<z.len: z[mu] := 0
  gg.upload(gg.u, gr)
  echo "check actionA  CPU: ", gc.actionA(gr), "  GPU: ", gg.actionA(gc)
  gc.forceA(gr, f)
  gg.upload(mom, z)
  gg.forceA(gc, mom, -1.0)
  gg.download(fg, mom)
  echo "check forceA  |GPU-CPU|^2/|CPU|^2: ", rel(fg, f)
  var pr = lo.newgauge
  threads: pr.randomTAH rc
  threads: axexpmuly(f, 0.3, pr, gr)
  gg.upload(mom, pr)
  gg.expUpdate(mom, 0.3)
  gg.download(fg, gg.u)
  echo "check exp  |GPU-CPU|^2/|CPU|^2: ", rel(fg, f)
  gg.upload(gg.u, gr)
  threads:
    gr.setBC
    threadBarrier()
    gr.stagPhase
    psi.gaussian rc
  let stag = newStag(gr)
  var ph = lo.ColorVector()
  var ps = lo.ColorVector()
  threads:
    stag.D(ph, psi, mass)
    threadBarrier()
    ph.odd := 0
  var sp = initSolverParams()
  sp.r2req = 1e-24
  sp.maxits = 10000
  sp.verbosity = 0
  stag.solve(ps, ph, mass, sp)
  f.oneLinkForce(ps, gr)
  s.setLinks(gg, sg)
  gpuMemCpyToGpu(x, addr psi[0], 6*s.n*sizeof(float))
  s.applyM(phi, x, mass)
  s.solveM(x, phi, mass, sp)
  gg.upload(mom, z)
  s.forceM(gg, sg, mom, x, -1.0)
  gg.download(fg, mom)
  var xg = lo.ColorVector()
  gpuMemCpyToCpu(addr xg[0], x, 6*s.n*sizeof(float))
  var dx, nx = 0.0
  threads:
    let a = norm2(xg - ps)
    let b = ps.norm2
    threadMaster:
      dx = a
      nx = b
  echo "check solve  |GPU-CPU|^2/|CPU|^2: ", dx/nx
  echo "check fermion force  |GPU-CPU|^2/|CPU|^2: ", rel(fg, f)

gg.upload(gg.u, g)
echo "plaq: ", gg.plaq
echo "actionA: ", gg.actionA(gc)

proc mdt(t: float) =
  tic()
  gg.expUpdate(mom, t)
  toc("mdt")

proc mdv(t: float) =
  tic()
  gg.forceA(gc, mom, t)
  toc("mdv")

proc setLinks() =
  s.setLinks(gg, sg)
  if ssp != nil: ssp[].setLinks(gg, sg)

proc mdvf(t: float) =
  tic()
  setLinks()
  toc("fforce links")
  s.solveM(x, phi, mass, spf, ssp)
  toc("fforce solve")
  s.forceM(gg, sg, mom, x, -0.5*t/mass)
  toc("mdvf")

proc mdvAll(t: openarray[float]) =
  if t[0] != 0: mdv t[0]
  if t[1] != 0: mdvf t[1]

let
  (VAll,T) = newIntegratorPair(mdvAll, mdt)
  H = newParallelEvolution(
    mkOmelyan2MN(steps = gsteps, V = VAll[0], T = T),
    mkOmelyan2MN(steps = fsteps, V = VAll[1], T = T))

for n in 1..trajs:
  tic()
  let t0 = epochTime()
  threads:
    p.randomTAH r
    psi.gaussian r
  gg.upload(mom, p)
  gg.copy(g0, gg.u)
  gpuMemCpyToGpu(x, addr psi[0], 6*s.n*sizeof(float))
  setLinks()
  s.applyM(phi, x, mass)
  toc("init traj")
  s.solveM(x, phi, mass, spa, ssp)
  toc("fa solve 1")
  let
    p2 = gg.norm2(mom)
    ga0 = gg.actionA(gc)
    fa0 = 0.5*s.norm2(x)
    k0 = 0.5*p2 - (16*vol).float
    h0 = ga0 + fa0 + k0
  toc("init gauge action")
  echo "Begin H: ",h0,"  Sg: ",ga0,"  Sf: ",fa0,"  T: ",k0

  H.evolve tau
  H.finish
  toc("evolve")

  let p2e = gg.norm2(mom)
  setLinks()
  s.solveM(x, phi, mass, spa, ssp)
  toc("fa solve 2")
  let
    ga1 = gg.actionA(gc)
    fa1 = 0.5*s.norm2(x)
    k1 = 0.5*p2e - (16*vol).float
    h1 = ga1 + fa1 + k1
  toc("final gauge action")
  echo "End H: ",h1,"  Sg: ",ga1,"  Sf: ",fa1,"  T: ",k1

  let
    dH = h1 - h0
    acc = exp(-dH)
    accr = R.uniform
  if accr <= acc:  # accept
    echo "ACCEPT:  dH: ",dH,"  exp(-dH): ",acc,"  r: ",accr
  else:  # reject
    echo "REJECT:  dH: ",dH,"  exp(-dH): ",acc,"  r: ",accr
    gg.copy(gg.u, g0)
    gg.fresh = false
  echo "plaq: ", gg.plaq, "  traj secs: ", epochTime() - t0

echoTimers()
qexfinalize()
