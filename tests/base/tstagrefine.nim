import testutils
import qex, physics/[qcdTypes,stagSolve]

qexInit()
disableParamFiltering()
let lo = @[8,8,8,8].newLayout
var rng = lo.newRngField(MRG32k3a,13579'u)
var g = lo.newGauge()
g.random rng
threads:
  g.projectSU
  threadBarrier()
  g.setBC
  threadBarrier()
  g.stagPhase
var s = newStag(g)
var b,x,r = lo.ColorVectorD()
threads: b.gaussian rng

proc residual(mass: float): float =
  var rr = 0.0
  threads:
    s.D(r,x,mass)
    threadBarrier()
    r := b-r
    let rn = r.norm2
    let bn = b.norm2
    threadMaster: rr = rn/bn
  rr

suite "mixed full equation refinement":
  for acc64 in [false,true]:
    test "inner precision floor, wide accumulator " & $acc64:
      var sp = initSolverParams()
      sp.backend = sbQex
      sp.sloppySolve = SloppySingle
      sp.cg.acc64 = acc64
      sp.r2req = 1e-24
      sp.maxits = 30000
      sp.verbosity = 0
      s.solve(x,b,0.001,sp)
      let rr = residual(0.001)
      check rr <= 1.01*sp.r2req
      check sp.iterations < sp.maxits
      check abs(rr-sp.r2.mean) <= 1e-8*sp.r2req
  test "failed mixed reconstruction retries the accepted residual":
    var sp = initSolverParams()
    sp.backend = sbQex
    sp.sloppySolve = SloppySingle
    sp.cg.kind = 0
    sp.r2req = 1e-16
    sp.maxits = 20000
    sp.verbosity = 0
    s.solve(x,b,0.000001,sp)
    let rr = residual(0.000001)
    check rr <= 1.01*sp.r2req
    check sp.iterations < sp.maxits
    check abs(rr-sp.r2.mean) <= 1e-8*sp.r2req
    check sp.sloppySolve == SloppySingle
  test "limited mixed reconstruction retains the initial residual":
    var sp = initSolverParams()
    sp.backend = sbQex
    sp.sloppySolve = SloppySingle
    sp.cg.kind = 0
    sp.r2req = 1e-16
    sp.maxits = 200
    sp.verbosity = 0
    s.solve(x,b,0.000001,sp)
    let rr = residual(0.000001)
    check rr <= 1.0
    check rr > sp.r2req
    check sp.iterations<=sp.maxits
    check abs(rr-sp.r2.mean) <= 1e-12
qexFinalize()
