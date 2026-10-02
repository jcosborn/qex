import testutils
import qex, physics/[qcdTypes, stagSolve]

qexInit()
let lo = latticeFromLocalLattice([4,4,4,4], nRanks).newLayout
var g = lo.newGauge()
g.unit()
var s = newStag(g)
var b, x = lo.ColorVector()
threads: b := 0

suite "zero staggered source":
  for mixed in [false,true]:
    test "single mass, mixed " & $mixed:
      for which in 0..2:
        var sp = initSolverParams()
        sp.backend = sbQex
        sp.verbosity = 0
        sp.sloppySolve = if mixed: SloppySingle else: SloppyNone
        sp.usePrevSoln = true
        for call in 1..2:
          threads: x := 1
          case which
          of 0: s.solveEE(x,b,0.1,sp)
          of 1: s.solveOO(x,b,0.1,sp)
          else: s.solve(x,b,0.1,sp)
          var x2 = 0.0
          threads:
            let t = x.norm2
            threadMaster: x2 = t
          check x2 == 0
          check sp.calls == call
          check sp.iterations == 0
          check sp.iterationsMax == 0
          check sp.r2.n == call
          check sp.r2.mean == 0
          check sp.r2.max == 0
  test "multiple masses":
    var xs = @[lo.ColorVector(),lo.ColorVector()]
    var sp = initSolverParams()
    sp.backend = sbQex
    sp.sloppySolve = SloppyNone
    sp.verbosity = 0
    threads:
      for y in xs: y := 1
    s.solve(xs,b,@[0.1,0.2],sp)
    for y in xs:
      var y2 = 0.0
      threads:
        let t = y.norm2
        threadMaster: y2 = t
      check y2 == 0
    check sp.calls == 1
    check sp.iterations == 0
    check sp.r2.mean == 0
    check sp.r2.max == 0
qexFinalize()
