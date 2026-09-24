## Staggered even-odd CG on the GPU, physics/stagGpu, against the CPU solver.
##   -lat:x,y,z,t -rg: ranks per dimension -ig: lanes per dimension (default innerGeom)
##   -mass -r2req -maxits; -ncpu, -ngpu: number of CPU and GPU solves
##   -mixed:1 also the mixed precision solve, -r2in its restart tolerance
##   -recon:0 links of 18 reals, -reals:14 rows 0, 1 and the determinant; -fwd: 1 forward links only, 0 both, -1 by size
##   -ipc:0 MPI for peers on the node too; -split: second hop split, 1, 0, -1 off-node
##   -nd:n times n applications of A on the CPU and the GPU; -prof:1 profile
## Prints the true residual of the GPU solution and its difference from the CPU one.
import qex
import physics/[qcdTypes, stagSolve, stagGpu]
import comms/halogpu
import backend/accel
import times

qexInit()
let lat = intSeqParam("lat", @[8,8,8,8])
let rg = intSeqParam("rg", @[1,1,1,1])
let ig = intSeqParam("ig", innerGeom(lat, rg, static(VLEN)))
let lo = newLayout(lat, VLEN, rg, ig)
var g = newSeq[type(lo.ColorMatrix())](lat.len)
for mu in 0..<lat.len: g[mu] = lo.ColorMatrix()
echo "lattice ", lat, "  ranks ", rg, "  lanes ", ig
let mass = floatParam("mass", 0.1)
let r2req = floatParam("r2req", 1e-16)
let maxits = intParam("maxits", 10000)
let ncpu = intParam("ncpu", 1)
let ngpu = intParam("ngpu", 2)
var rng = newRngField(lo, MRG32k3a, 987654321'u)
g.random rng
if intParam("reunit", 1) != 0:  # SU(3) to rounding, for links of 12 reals
  threads:
    g.projectSU
threads:
  g.setBC
  threadBarrier()
  g.stagPhase
var s = newStag(g)
var src = lo.ColorVector()
var xc = lo.ColorVector()
var xg = lo.ColorVector()
var d = lo.ColorVector()
src.gaussian rng
threads:
  src.odd := 0
  xg := 0

var sp = initSolverParams()
sp.r2req = r2req
sp.maxits = maxits
sp.verbosity = 1
for k in 0..<ncpu:
  resetTimers()
  s.solveEE(xc, src, mass, sp)
  echoProf()

haloIpc = intParam("ipc", 1) != 0
hopSplit = intParam("split", -1)
let reals = intParam("reals", if intParam("recon", 1) != 0: 12 else: 18)
let fwd = intParam("fwd", -1)
var sg = newStagGpu(g, float64, reals, fwd)
echo "GPU links: ", sg.nl, " reals", if sg.lb == nil: ", forward only" else: ""
for k in 0..<ngpu:
  var spg = initSolverParams()
  spg.r2req = r2req
  spg.maxits = maxits
  spg.verbosity = intParam("verb", 1)
  resetTimers()
  sg.solveEE(xg, src, mass, spg)
  echoProf()

let nd = intParam("nd", 0)
if nd > 0:  # A = 4m^2 - D_eo D_oe applications, flops as stagD2xx
  let gf = 1e-9*float(1158*lo.nEven*nRanks*nd)
  var t0 = epochTime()
  threads:
    for k in 0..<nd:
      stagD2ee(s.se, s.so, d, s.g, src, mass*mass)
      threadBarrier()
  var dt = epochTime() - t0
  echo "CPU A: ", 1e6*dt/nd.float, " us  ", gf/dt, " Gflops"
  let x = cast[ptr UncheckedArray[float64]](gpuMalloc(6*sg.n*sizeof(float64)))
  let r = cast[ptr UncheckedArray[float64]](gpuMalloc(6*sg.n*sizeof(float64)))
  let t = cast[ptr UncheckedArray[float64]](gpuMalloc(6*sg.n*sizeof(float64)))
  sg.upload(x, src)
  sg.applyD2ee(r, x, t, mass*mass)
  getDefaultComm().barrier
  t0 = epochTime()
  for k in 0..<nd:
    sg.applyD2ee(r, x, t, mass*mass)
  dt = epochTime() - t0
  echo "GPU A: ", 1e6*dt/nd.float, " us  ", gf/dt, " Gflops"
  for q in [x, r, t]: gpuFree(q)

if intParam("mixed", 0) != 0:
  var sgs = newStagGpu(g, float32, reals, fwd)
  for k in 0..<ngpu:
    var spg = initSolverParams()
    spg.r2req = r2req
    spg.maxits = maxits
    spg.verbosity = intParam("verb", 1)
    resetTimers()
    solveEE(sg, sgs, xg, src, mass, spg, floatParam("r2in", 1e-6))
    echoProf()
  sgs.free

threads:
  stagD2ee(s.se, s.so, d, s.g, xg, mass*mass)
  threadBarrier()
  d.even -= src
  let rg = d.even.norm2
  let sn = src.even.norm2
  d.even := xg - xc
  let dx = d.even.norm2
  let xn = xc.even.norm2
  echo "GPU true r2/b2: ", rg/sn, "  |xg-xc|^2/|xc|^2: ", dx/xn
sg.free
qexFinalize()
