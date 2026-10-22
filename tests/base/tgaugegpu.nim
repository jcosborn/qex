import testutils
import qex, gauge, physics/qcdTypes
import gauge/gaugeGpu
import backend/accel

# gauge/gaugeGpu: the action and force of plaq + adjplaq agree with actionA
# and forceA, and the coefficients they leave out, rect and pgm, raise
# ValueError as there, before the force changes its momenta.

qexInit()
let lo = latticeFromLocalLattice([4,4,4,8], nRanks).newLayout
var g = lo.newGauge
var r = lo.newRNGField(RngMilc6, 987654321)
g.random r
threads: g.projectSU
var gg = newGpuGauge(lo)
gg.upload(gg.u, g)
let p = gg.newLinks

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

suite "GPU gauge action":
  let gc = GaugeActionCoeffs(plaq: 5.4, adjplaq: -1.35)

  test "actionA":
    let a = gc.actionA(g)
    let b = gg.actionA(gc)
    check abs(b - a) <= 1e-12*abs(a)

  test "forceA":
    var f = lo.newGauge
    var z = lo.newGauge
    var fg = lo.newGauge
    threads:
      for mu in 0..<z.len: z[mu] := 0
    gc.forceA(g, f)
    gg.upload(p, z)
    gg.forceA(gc, p, -1.0)  # p = F
    gg.download(fg, p)
    check rel(fg, f) <= 1e-26

  test "rect and pgm raise ValueError":
    var q = lo.newGauge
    var qg = lo.newGauge
    q.gaussian r
    for c in [GaugeActionCoeffs(plaq: 5.4, rect: -0.3), GaugeActionCoeffs(plaq: 5.4, pgm: 0.1),
              GaugeActionCoeffs(rect: -0.3), GaugeActionCoeffs(plaq: 5.4, rect: -0.3, pgm: 0.1, adjplaq: -1.35)]:
      expect ValueError: discard c.actionA(g)
      expect ValueError: discard gg.actionA(c)
      gg.upload(p, q)
      expect ValueError: gg.forceA(c, p, 1.0)
      gg.download(qg, p)
      check rel(qg, q) == 0.0

qexFinalize()
