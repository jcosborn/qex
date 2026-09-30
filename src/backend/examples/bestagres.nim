## Checks the stopping rules of the solvers of physics/stagGpu against the
## CPU operators, stag.D and stagD2ee, on random SU(3) links, for masses of
## both signs and sources on all sites (full) or on the even ones, whose
## odd sites hold 1e30 to show they are not read:
##   solveM   |b - M x|^2 <= r2req |b|^2, M = m + D/2, of one system and of
##            several (nBatch per hop), in double and mixed precision
##   solveEE  |b_e - A x_e|^2 <= r2req |b_e|^2
## with fixed and atomic dot products.  solveM records the residual it
## reaches, which the CPU one matches, leaves the source unchanged, stops
## within maxits at a residual above the request when maxits is too small,
## and gives x = 0 for b = 0.  Exits with 1 if a check fails.
##   -lat: lattice (default 4^4 per rank), -reals: links of 18 or 12 reals,
##   -ipc:0 MPI for peers on the node too, -split: second hop split around
##   the exchange, 1, 0, -1 with off-node neighbors, -verb: verbosity
import qex
import physics/[qcdTypes, stagSolve, stagGpu]
import comms/halogpu
import backend/accel
import strformat

qexInit()
let lat = intSeqParam("lat", latticeFromLocalLattice(@[4,4,4,4], nRanks))
let lo = lat.newLayout
echo "lattice ", lo.physGeom, "  ranks ", lo.rankGeom, "  lanes ", lo.innerGeom, "  nBatch ", nBatch
var g = lo.newGauge
var rng = newRngField(lo, MRG32k3a, 987654321'u)
g.random rng
threads:
  g.projectSU
  threadBarrier()
  g.setBC
  threadBarrier()
  g.stagPhase
var s = newStag(g)
haloIpc = intParam("ipc", 1) != 0
hopSplit = intParam("split", -1)
var sd = newStagGpu(g, float64, intParam("reals", 18), batch = true)
var sf = newStagGpu(g, float32, intParam("reals", 18), batch = true)
echo "GPU links: ", sd.nl, " reals", "  split hop: ", sd.split, "  peer IPC: ", haloIpc
let n6 = 6*sd.n
let ne6 = 6*sd.ne
const big = 100000  # maxits of the solves that converge
let verb = intParam("verb", 0)
var fails = 0

proc dvec(): ptr UncheckedArray[float] = cast[ptr UncheckedArray[float]](gpuMalloc(n6*sizeof(float)))

type Src = object
  ## a source on the device and its host copy
  d: ptr UncheckedArray[float]
  h: seq[float]
  full: bool

proc newSrc(scale: float; full: bool): Src =
  ## scale times a gaussian source, b_o = 1e30 unless full
  var f = lo.ColorVector()
  threads: f.gaussian rng
  result.h = newSeq[float](n6)
  copyMem(addr result.h[0], addr f[0], n6*sizeof(float))
  for i in 0..<n6:
    result.h[i] = if full or i < ne6: scale*result.h[i] else: 1e30
  result.full = full
  result.d = dvec()
  gpuMemCpyToGpu(result.d, addr result.h[0], n6*sizeof(float))

proc resid(x: ptr UncheckedArray[float]; b: Src; m: float): float =
  ## |b - M x|^2/|b|^2 by stag.D, b_o = 0 unless full; |M x|^2 for b = 0
  var hx, hb, r = lo.ColorVector()
  gpuMemCpyToCpu(addr hx[0], x, n6*sizeof(float))
  copyMem(addr hb[0], unsafeAddr b.h[0], n6*sizeof(float))
  var rr, bb = 0.0
  threads:
    if not b.full: hb.odd := 0
    threadBarrier()
    s.D(r, hx, m)
    threadBarrier()
    r := hb - r
    let t = r.norm2
    let u = hb.norm2
    threadMaster:
      rr = t
      bb = u
  if bb > 0: rr/bb else: rr

proc unchanged(b: Src): bool =
  ## the device source still holds its host copy
  var h = newSeq[float](n6)
  gpuMemCpyToCpu(addr h[0], b.d, n6*sizeof(float))
  var d = 0.0
  for i in 0..<n6:
    if h[i] != b.h[i]: d += 1
  getDefaultComm().allReduce(d)
  d == 0

proc check(label: string; x: ptr UncheckedArray[float]; b: Src; m: float; sp: SolverParams;
           r2req: float; maxits: int; zero = false) =
  ## the solve of sp, its only one: with maxits = big it meets r2req, with
  ## fewer it misses it within maxits; its recorded residual matches the CPU one
  let r2 = resid(x, b, m)
  let rec = sp.r2.max
  var ok = unchanged(b) and sp.calls == 1 and sp.iterations <= maxits
  if zero: ok = ok and r2 == 0 and rec == 0 and sp.iterations == 0
  else: ok = ok and abs(rec - r2) <= 1e-3*max(rec, r2) + 1e-6*r2req and (rec <= r2req) == (maxits == big)
  echo &"{label:<28} m {m:6} its {sp.iterations:5}  |b-Mx|^2/|b|^2 recorded {rec:9.3e} CPU {r2:9.3e}  requested {r2req:7.1e}",
       if ok: "  ok" else: "  FAILED"
  if not ok: inc fails

proc solveOne(label: string; m, r2req: float; maxits: int; full, mixed: bool; scale = 1.0) =
  let b = newSrc(scale, full)
  let x = dvec()
  var sp = initSolverParams()
  sp.r2req = r2req
  sp.maxits = maxits
  sp.verbosity = verb
  sd.solveM(x, b.d, m, sp, if mixed: addr sf else: nil, full = full)
  check(label, x, b, m, sp, r2req, maxits, scale == 0.0)
  gpuFree(x)
  gpuFree(b.d)

proc solveMany(label: string; ns: int; full, mixed: bool) =
  ## ns systems of different masses and requests, system 3 with 5
  ## iterations only, system 5 with b = 0
  const ms = [0.1, -0.02, 0.05, 0.01, -0.2, 0.03, 0.005, 0.07, -0.04]
  const rq = [1e-12, 1e-10, 1e-14, 1e-12, 1e-16, 1e-12, 1e-20, 1e-12, 1e-14]
  var bs: seq[Src]
  var xs, bd: seq[ptr UncheckedArray[float]]
  var sps = newSeq[SolverParams](ns)
  for j in 0..<ns:
    bs.add newSrc(if j == 5: 0.0 else: 1.0, full)
    bd.add bs[j].d
    xs.add dvec()
    sps[j] = initSolverParams()
    sps[j].r2req = rq[j]
    sps[j].maxits = if j == 3: 5 else: big
    sps[j].verbosity = verb
  sd.solveM(xs, bd, ms[0..<ns], sps, if mixed: addr sf else: nil, full = full)
  for j in 0..<ns:
    check(&"{label} {j}/{ns}", xs[j], bs[j], ms[j], sps[j], rq[j], sps[j].maxits, j == 5)
    gpuFree(xs[j])
    gpuFree(bs[j].d)

proc solveEEcheck(label: string; m, r2req: float; mixed: bool) =
  ## the true residual of solveEE by stagD2ee
  var src = lo.ColorVector()
  var xg = lo.ColorVector()
  var d = lo.ColorVector()
  threads:
    src.gaussian rng
    src.odd := 0
    xg := 0
  var sp = initSolverParams()
  sp.r2req = r2req
  sp.maxits = big
  sp.verbosity = verb
  if mixed: solveEE(sd, sf, xg, src, m, sp)
  else: sd.solveEE(xg, src, m, sp)
  var r2 = 0.0
  threads:
    stagD2ee(s.se, s.so, d, s.g, xg, m*m)
    threadBarrier()
    d.even -= src
    let t = d.even.norm2
    let u = src.even.norm2
    threadMaster: r2 = t/u
  let ok = r2 <= 1.001*r2req and sp.iterations <= big
  echo &"{label:<28} m {m:6} its {sp.iterations:5}  |b-Ax|^2/|b|^2 CPU {r2:9.3e}  requested {r2req:7.1e}",
       if ok: "  ok" else: "  FAILED"
  if not ok: inc fails

for fixed in [true, false]:
  sd.fixed = fixed
  sf.fixed = fixed
  let f = if fixed: "fixed" else: "atomic"
  for mixed in [false, true]:
    let p = if mixed: "mixed" else: "double"
    for full in [true, false]:
      let e = if full: "full" else: "even"
      for m in [0.1, 0.01, 0.001, -0.01]:
        solveOne(&"{f} {p} {e}", m, 1e-12, big, full, mixed)
      solveOne(&"{f} {p} {e} tight", 0.001, 1e-24, big, full, mixed)
      solveOne(&"{f} {p} {e} maxits 5", 0.01, 1e-12, 5, full, mixed)
      solveOne(&"{f} {p} {e} zero", 0.01, 1e-12, big, full, mixed, 0.0)
      for ns in [2, nBatch, nBatch + 1, 2*nBatch + 1]:
        solveMany(&"{f} {p} {e} batch", ns, full, mixed)
    solveEEcheck(&"{f} {p} solveEE", 0.01, 1e-14, mixed)
    solveEEcheck(&"{f} {p} solveEE", 0.001, 1e-20, mixed)

echo if fails == 0: "all solves ok" else: $fails & " solves FAILED"
qexExit(if fails == 0: 0 else: 1)
