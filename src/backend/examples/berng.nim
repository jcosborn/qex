## Checks that the generators of RNG fields on the GPU, rng/rngGpu, draw
## the numbers of the host ones: for each RNG, gaussian into a color vector
## field and randomTAH into 4 matrix fields on both, then the generators.
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

proc check(R: typedesc; name: string) =
  var r = lo.newRNGField(R, 987654321)
  var g = newRngGpu(lo, r)
  var v = lo.ColorVector()
  var p = lo.newgauge
  threads:
    v.gaussian r
    p.randomTAH r
  g.gaussian(dv, 6)
  g.randomTAH dp
  var vg = lo.ColorVector()
  var pg = lo.newgauge
  gpuMemCpyToCpu(addr vg[0], dv, 6*n*sizeof(float))
  for mu in 0..3:
    gpuMemCpyToCpu(addr pg[mu][0], addr dp[18*mu*n], 18*n*sizeof(float))
  var rd = lo.newRNGField(R, 1)
  g.download rd
  var nd = 0.0  # generators differing from the host ones
  for j in rd.l.sites:
    if rd[j] != r[j]: nd += 1
  var dvv, nvv, dpp, npp = 0.0
  threads:
    let a = norm2(vg - v)
    let b = v.norm2
    var c, d = 0.0
    for mu in 0..3:
      c += norm2(pg[mu] - p[mu])
      d += p[mu].norm2
    threadMaster:
      dvv = a
      nvv = b
      dpp = c
      npp = d
  getDefaultComm().allReduce(nd)
  let ok = nd == 0 and dvv/nvv < 1e-28 and dpp/npp < 1e-28
  if not ok: inc fails
  echo name, "  |GPU-CPU|^2/|CPU|^2 gaussian: ", dvv/nvv, "  randomTAH: ", dpp/npp,
    "  generators differing: ", nd, if ok: "  ok" else: "  FAILED"
  g.free

check(RngMilc6, "RngMilc6")
check(MRG32k3a, "MRG32k3a")
check(Philox4x64, "Philox4x64")
check(Threefry4x64, "Threefry4x64")
echo if fails == 0: "all generators ok" else: $fails & " generators FAILED"
qexFinalize()
