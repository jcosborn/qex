## Checks the generators of RNG fields in the SIMD layout, RNGFieldV, and
## on the GPU, rng/rngGpu, against the host RNG fields: for each RNG,
## gaussian into a color vector field and randomTAH into 4 matrix fields
## from each, then the generators.
##   -lat: lattice, -rg: ranks per dimension
import qex, gauge, physics/qcdTypes
import backend/accel, rng/rngGpu

qexInit()
let lat = intSeqParam("lat", @[8,8,8,8])
let rg = intSeqParam("rg")
let lo = if rg.len > 0: newLayout(lat, VLEN, rg) else: lat.newLayout
let n = lo.nSites
echo "lattice ", lat, "  lanes ", lo.innerGeom
let dv = cast[ptr UncheckedArray[float]](gpuMalloc(6*n*sizeof(float)))
let dp = cast[ptr UncheckedArray[float]](gpuMalloc(4*18*n*sizeof(float)))
var fails = 0

proc rel(a, b: auto): float =
  ## |a-b|^2/|b|^2, the gauge fields as one
  var d, m = 0.0
  threads:
    var dt, mt = 0.0
    when a is seq:
      for mu in 0..<a.len:
        dt += norm2(a[mu] - b[mu])
        mt += b[mu].norm2
    else:
      dt = norm2(a - b)
      mt = b.norm2
    threadMaster:
      d = dt
      m = mt
  d/m

proc gpuDraws(g: auto; v: auto; p: auto) =
  ## v and p from the generators g on the GPU
  g.gaussian(dv, 6)
  g.randomTAH dp
  gpuMemCpyToCpu(addr v[0], dv, 6*n*sizeof(float))
  for mu in 0..3:
    gpuMemCpyToCpu(addr p[mu][0], addr dp[18*mu*n], 18*n*sizeof(float))

proc check(R: typedesc; name: string) =
  var r = lo.newRNGField(R, 987654321)
  var rv = lo.newRNGFieldV(R, 987654321)
  var g = newRngGpu(lo, r)
  var gv = newRngGpu(rv)
  var v, vv, vg, vgv = lo.ColorVector()
  var p = lo.newgauge
  var pv = lo.newgauge
  var pg = lo.newgauge
  var pgv = lo.newgauge
  threads:
    v.gaussian r
    p.randomTAH r
    vv.gaussian rv
    pv.randomTAH rv
  g.gpuDraws(vg, pg)
  gv.gpuDraws(vgv, pgv)
  var rd = lo.newRNGField(R, 1)
  g.download rd
  var rvd = lo.newRNGFieldV(R, 1)
  gv.download rvd
  var nd = 0.0  # generators differing from the host ones
  for j in rd.l.sites:
    if rd[j] != r[j]: nd += 1
  if rvd.s != rv.s: nd += 1
  var ld = 0.0  # RNGFieldV draws differing from the RNGField ones
  for j in rd.l.sites:
    if r[j] != rv[lo.rankIndex(rd.l.coords, j).index]: ld += 1
  let d = [rel(vg, v), rel(pg, p), rel(vv, v), rel(pv, p), rel(vgv, vg), rel(pgv, pg)]
  getDefaultComm().allReduce(nd)
  getDefaultComm().allReduce(ld)
  let ok = nd == 0 and ld == 0 and d[0] < 1e-28 and d[1] < 1e-28 and d[2] == 0 and d[3] == 0 and d[4] == 0 and d[5] == 0
  if not ok: inc fails
  echo name, "  GPU-host |d|^2/|x|^2 gaussian: ", d[0], "  randomTAH: ", d[1],
    "  RNGFieldV-RNGField host: ", d[2], " ", d[3], "  GPU: ", d[4], " ", d[5],
    "  generators differing: ", nd, " ", ld, if ok: "  ok" else: "  FAILED"
  g.free
  gv.free

check(RngMilc6, "RngMilc6")
check(MRG32k3a, "MRG32k3a")
check(Philox4x64, "Philox4x64")
check(Threefry4x64, "Threefry4x64")
echo if fails == 0: "all generators ok" else: $fails & " generators FAILED"
qexFinalize()
