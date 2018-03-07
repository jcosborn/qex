import qex
import testutils

suite "Test gauge IO":
  qexInit()
  threads: echo "thread ",threadNum," / ",numThreads
  const fn = "tmplat.lime"
  var
    (l,g,r) = setupLattice([8,8,8,8])
    p = g.plaq

  test "save double precision (default)":
    if 0 != g.saveGauge(fn):
      echo "Error: failed to save gauge to ",fn
      qexExit 1
    var gg = l.newGauge
    if 0 != gg.loadGauge(fn):
      echo "Error: failed to load gauge from ",fn
      qexExit 1
    #threads: gg.projectSU
    #CT = 1E-10    # Larger error would result from projectSU due to random gauge.
    var pp = gg.plaq
    check(p ~ pp)

  test "save single precision":
    if 0 != g.saveGauge(fn,"F"):
      echo "Error: failed to save gauge to ",fn
      qexExit 1
    var gg = l.newGauge
    if 0 != gg.loadGauge(fn):
      echo "Error: failed to load gauge from ",fn
      qexExit 1
    #threads: gg.projectSU
    var pp = gg.plaq
    CT = 1E-5
    check(p ~ pp)

  qexFinalize()
