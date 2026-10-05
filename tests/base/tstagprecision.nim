import testutils
import qex, physics/[qcdTypes,stagSolve]
import std/strformat

qexInit()
disableParamFiltering()
# Staggered operator scratch is shared per field type: one layout per process.
let lo = @[8,8,8,8].newLayout
let mass = 0.0001
let seed = 987654321'u
let source = 3
var rng = lo.newRngField(MRG32k3a,seed)
var g = lo.newGauge()
g.random rng
threads:
  g.projectSU
  threadBarrier()
  g.setBC
  threadBarrier()
  g.stagPhase
var s = newStag(g)
var ri = lo.newRngField(MRG32k3a,seed xor 192837465'u)
var bs = lo.ColorVectorS()
var b,x,r,save = lo.ColorVectorD()
for i in 1..source:
  threads: bs.gaussian ri

proc residual(even: bool): float =
  var rr = 0.0
  let par = if even: "even" else: "odd"
  threads:
    if even: stagD2ee(s.se,s.so,r,s.g,x,mass*mass)
    else: stagD2oo(s.se,s.so,r,s.g,x,mass*mass)
    threadBarrier()
    r[par] := b-r
    let rn = r[par].norm2
    let bn = b[par].norm2
    threadMaster: rr = rn/bn
  rr

echo &"Residual fixture: lattice {lo.physGeom}, mass {mass}, seed {seed}, source {source}"
suite "CPU CG true residual":
  for sloppy in [SloppyNone,SloppySingle]:
    let name = if sloppy==SloppyNone: "double precision CG" else: "CG with a single precision inner solve"
    test name:
      threads:
        b := bs
        threadBarrier()
        b.odd := 0
        threadBarrier()
        save := b
      var sp = initSolverParams()
      sp.backend = sbQex
      sp.sloppySolve = sloppy
      sp.r2req = 1e-12
      sp.maxits = 100000
      sp.verbosity = 0
      s.solveEE(x,b,mass,sp)
      let rr = residual(true)
      echo &"  {name}: iterations {sp.iterations}, true r2/b2 {rr:.6e}, requested {sp.r2req:.1e}"
      check rr <= 1.01*sp.r2req
      check sp.iterations>0 and sp.iterations<=sp.maxits
      var dx = 0.0
      threads:
        let d = norm2diff(b,save)
        threadMaster: dx = d
      check dx==0.0

suite "mixed checkerboard residuals":
  for even in [true,false]:
    for mode in 0..3:
      test "parity " & $even & " mode " & $mode:
        threads:
          b := bs
          threadBarrier()
          if even: b.odd := 0
          else: b.even := 0
          threadBarrier()
          save := b
        var sp = initSolverParams()
        sp.backend = sbQex
        sp.sloppySolve = SloppySingle
        sp.cg.kind = if mode==0: 0 elif mode==1: -1 else: 1
        sp.cg.acc64 = mode==3
        sp.r2req = 1e-12
        sp.maxits = 100000
        sp.verbosity = 0
        let kind = sp.cg.kind
        if even: s.solveEE(x,b,mass,sp)
        else: s.solveOO(x,b,mass,sp)
        let rr = residual(even)
        check rr <= 1.01*sp.r2req
        check abs(rr-sp.r2.mean) <= 1e-8*sp.r2req
        check sp.calls==1 and sp.r2.n==1
        check sp.iterations>0 and sp.iterations<=sp.maxits
        check sp.cg.kind==kind and sp.sloppySolve==SloppySingle
        var dx = 0.0
        threads:
          let d = norm2diff(b,save)
          threadMaster: dx = d
        check dx==0.0
  test "periodic checks allow a nondecreasing CG residual":
    var ug = lo.newGauge()
    ug.unit()
    threads:
      ug.setBC
      threadBarrier()
      ug.stagPhase
      b := 1
      threadBarrier()
      b.odd := 0
    var us = newStag(ug)
    var sp = initSolverParams()
    sp.backend = sbQex
    sp.sloppySolve = SloppySingle
    sp.cg.kind = 1
    sp.cg.delta = 0.0
    sp.cg.period = 1
    sp.cg.maxInc = 0
    sp.cg.maxTotal = 0
    sp.r2req = 1e-24
    sp.maxits = 2
    sp.verbosity = 0
    # This source has two eigenvalues: the first CG residual has the same
    # norm as b; the second vanishes. A restart cannot fit in this budget.
    us.solveEE(x,b,0.0,sp)
    var rr = 0.0
    threads:
      stagD2ee(us.se,us.so,r,us.g,x,0.0)
      threadBarrier()
      r.even := b-r
      let rn = r.even.norm2
      let bn = b.even.norm2
      threadMaster: rr = rn/bn
    check rr <= sp.r2req
    check sp.iterations==2
  for kind in [0,1]:
    test "iteration limit retains the best residual, kind " & $kind:
      threads:
        b := bs
        threadBarrier()
        b.odd := 0
      var sp = initSolverParams()
      sp.backend = sbQex
      sp.sloppySolve = SloppySingle
      sp.cg.kind = kind
      sp.r2req = 1e-12
      sp.maxits = 5000
      sp.verbosity = 0
      s.solveEE(x,b,mass,sp)
      let rr = residual(true)
      check rr <= 1.0
      check rr > sp.r2req
      check sp.iterations>0 and sp.iterations<=sp.maxits
      check sp.iterationsMax<=sp.maxits
      check abs(rr-sp.r2.mean) <= 1e-12
qexFinalize()
