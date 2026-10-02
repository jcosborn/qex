## Mixed staggered solves retain the full equation's accuracy across CG controls.
import testutils
import qex, physics/[qcdTypes,stagSolve]

qexInit()
disableParamFiltering()
let lo = latticeFromLocalLattice(@[4,4,4,8],nRanks).newLayout
var rng = lo.newRngField(MRG32k3a,918273645'u)
var g = lo.newGauge()
g.random rng
threads:
  g.projectSU
  threadBarrier()
  g.setBC
  threadBarrier()
  g.stagPhase
var s = newStag(g)
var b,x,r,save = lo.ColorVector()

proc residual(m: float): float =
  var value = 0.0
  threads:
    s.D(r,x,m)
    threadBarrier()
    r := b-r
    let rr = r.norm2
    let bb = b.norm2
    threadMaster: value = if bb > 0.0: rr/bb else: rr
  value

suite "mixed staggered CG":
  for mode in 0..6:
    for par in ["all","even","odd"]:
      test "mode " & $mode & " source " & par:
        threads:
          b.gaussian rng
          threadBarrier()
          if par=="even": b.odd := 0
          if par=="odd": b.even := 0
          threadBarrier()
          save := b
        var sp = initSolverParams()
        sp.backend = sbQex
        sp.sloppySolve = SloppySingle
        sp.r2req = 1e-16
        sp.maxits = 20000
        sp.verbosity = 0
        sp.cg.kind = if mode==0: 0 elif mode==6: -1 else: 1
        sp.cg.acc64 = mode!=2
        sp.cg.beta = mode!=3
        sp.cg.keep = mode!=4
        sp.cg.delta = if mode==5: 0.0 else: 0.1
        if mode==5: sp.cg.period = 0
        if mode==6:
          sp.cg.acc64 = false
          sp.cg.beta = false
          sp.cg.delta = 0.1
        let m = if par=="odd": -0.05 else: 0.05
        s.solve(x,b,m,sp)
        let rr = residual(m)
        check rr <= 1.01*sp.r2req
        check abs(rr-sp.r2.mean) <= 1e-8*sp.r2req
        check sp.calls==1 and sp.r2.n==1 and sp.iterations<=sp.maxits
        if mode==6:
          check sp.cg.kind == -1
          if par=="all": check sp.reliable>0
          else: check sp.reliable==0
        var dx = 0.0
        threads:
          let d = norm2diff(b,save)
          threadMaster: dx = d
        check dx==0.0
        sp.resetStats()
        sp.usePrevSoln = true
        s.solve(x,b,m,sp)
        check sp.iterations==0
        check residual(m) <= 1.01*sp.r2req
  test "reliable CG iteration budget":
    var sp = initSolverParams()
    sp.backend = sbQex
    sp.sloppySolve = SloppySingle
    sp.cg.kind = 1
    sp.r2req = 1e-20
    sp.maxits = 5
    sp.verbosity = 0
    s.solve(x,b,0.01,sp)
    check sp.iterations==5 and sp.iterationsMax<=5
    check residual(0.01)>sp.r2req
  test "precise residual floor leaves room for full M refinement":
    let ll = lo
    var rg = ll.newRngField(MRG32k3a,987654321'u)
    var gg = ll.newGauge()
    gg.random rg
    threads:
      gg.projectSU
      threadBarrier()
      gg.setBC
      threadBarrier()
      gg.stagPhase
    var ss = newStag(gg)
    var bb,xx,rr = ll.ColorVector()
    # Source 191 of bestagres stalled above the FP64 A-residual floor.
    for i in 1..191:
      threads: bb.gaussian rg
    var sp = initSolverParams()
    sp.backend = sbQex
    sp.sloppySolve = SloppySingle
    sp.cg.kind = 1
    sp.cg.acc64 = false
    sp.cg.beta = false
    sp.cg.delta = 0.1
    sp.cg.period = 1000
    sp.r2in = 0.0
    sp.r2req = 1e-24
    sp.maxits = 50000
    sp.verbosity = 0
    ss.solve(xx,bb,0.001,sp)
    var value = 0.0
    threads:
      ss.D(rr,xx,0.001)
      threadBarrier()
      rr := bb-rr
      let r2 = rr.norm2
      let b2 = bb.norm2
      threadMaster: value = r2/b2
    check value <= 1.001*sp.r2req
    check sp.iterations < sp.maxits and sp.reliable>0
qexFinalize()
