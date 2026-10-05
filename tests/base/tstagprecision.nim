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

qexFinalize()
