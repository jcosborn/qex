import testutils
import qex, os

qexInit()

threads: echo "thread ",threadNum," / ",numThreads
const fn = "tmplat.lime"
var
  lat = latticeFromLocalLattice([8,8,8,8], nRanks)
  (l,g,_) = setupLattice(lat)
  p = g.plaq
  gt = g.newGaugeS
  gs = gt.newGauge
  ps = gs.plaq
echo "plaq64: ", p
echo "plaq32: ", ps

suite "Test gauge IO":

  test "save double precision (default)":
    let err = g.saveGauge(fn)
    check(err == 0)

  test "load double precision":
    var gg = l.newGauge
    let err = gg.loadGauge(fn)
    check(err == 0)
    var pp = gg.plaq
    check(p == pp)
    for i in 0..<g.len:
      gg[i] -= g[i]
      let n2 = gg[i].norm2
      check(n2 == 0)

  test "save single precision":
    let err = g.saveGauge(fn,"F")
    check(err == 0)

  test "load single precision":
    var gg = l.newGauge
    let err = gg.loadGauge(fn)
    check(err == 0)
    var pp = gg.plaq
    check(ps == pp)
    for i in 0..<g.len:
      gg[i] -= gs[i]
      let n2 = gg[i].norm2
      check(n2 == 0)

removeFile fn
qexFinalize()
