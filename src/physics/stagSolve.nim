import times, strformat
export strformat  # workaround #7632
import base
import layout
import field
import stagD
export stagD
import solvers/bicgstab
import solvers/gcr
import maths
import quda/qudaWrapper
import grid/Grid

#var precon = false

proc solveEO*(s: Staggered; r,x: Field; m: SomeNumber; sp0: var SolverParams) =
  var sp = sp0
  sp.subset.layoutSubset(r.l, sp.subsetName)
  var t = newOneOf(r)
  var top = 0.0
  proc op2(a,b: Field) =
    threadBarrier()
    if threadNum==0: top -= epochTime()
    stagD2ee(s.se, s.so, a, s.g, b, m*m)
    if threadNum==0: top += epochTime()
    #threadBarrier()
  var oa = (apply: op2)
  let t0 = epochTime()
  cgSolve(r, x, oa, sp)
  let t1 = epochTime()
  let secs = t1-t0
  let flops = (s.g.len*4*72+60)*r.l.nEven*sp.finalIterations
  sp0.finalIterations = sp.finalIterations
  sp0.seconds = secs
  if sp0.verbosity>0:
    echo "op time: ", top
    echo "solve time: ", secs, "  Gflops: ", 1e-9*flops.float/secs

# multimass (trivial version with multiple single mass calls for now)
proc solveEO*(s: Staggered; r: seq[Field]; x: Field; m: seq[float];
              sp: seq[SolverParams]) =
  let n = m.len
  doAssert(r.len == n)
  doAssert(sp.len == n)
  for i in 0..<n:
    solveEO(s, r[i], x, m[i], sp[i])

proc solveEO*(s: Staggered; r: seq[Field]; x: Field; m: seq[float];
              sp: var SolverParams) =
  let n = m.len
  doAssert(r.len == n)
  doAssert(sp.len == n)
  for i in 0..<n:
    solveEO(s, r[i], x, m[i], sp)

proc solveXX*(s: Staggered; r,x: Field; m: SomeNumber; sp0: var SolverParams;
              parEven = true) =
  tic("solveXX")
  var sp = sp0
  sp.resetStats()
  dec sp.verbosity
  threads:
    r := 0
  case sp.backend
  of sbQex:
    tic("sbQex")
    proc opSolveXX(a,b: Field) =
      tic("solveXX>sbQex>op")
      threadBarrier()
      if parEven:
        stagD2ee(s.se, s.so, a, s.g, b, m*m)
        toc("stagD2ee")
      else:
        stagD2oo(s.se, s.so, a, s.g, b, m*m)
        toc("stagD2oo")
      #threadBarrier()
    var cg = newCgState(r, x)
    if parEven:
      sp.subset.layoutSubset(r.l, "even")
    else:
      sp.subset.layoutSubset(r.l, "odd")
    #if precon:
      #var oap = (apply: op, applyPrecon: oppre)
      #cg.solve(oap, sp)
    #else:
    var oa = (apply: opSolveXX, precon: cpNone)
    cg.solve(oa, sp)
    toc("cg.solve")
    sp.calls = 1
    sp.seconds = getElapsedTime()
    let flops = (s.g.len*4*72+60)*r.l.nEven*sp.iterations
    sp.flops = flops.float
    if sp0.verbosity>0:
      if parEven:
        echo "solveEE(QEX): ", sp.getStats
      else:
        echo "solveOO(QEX): ", sp.getStats
    toc("end sbQex")
  of sbQuda:
    tic("sbQuda")
    if parEven:
      #echo x.even.norm2, " ", sp.r2req
      s.qudaSolveEE(r,x,m,sp)
      toc("qudaSolveEE")
    else:
      s.qudaSolveOO(r,x,m,sp)
      toc("qudaSolveOO")
    sp.calls = 1
    sp.seconds = getElapsedTime()
    let flopsquda = (s.g.len*4*72+60)*r.l.nEven*sp.iterations
    sp.flops = flopsquda.float
    if sp0.verbosity>0:
      echo "solveXX(QUDA): ", sp.getStats
  of sbGrid:
    tic("sbQuda")
    if parEven:
      s.gridSolveEE(r,x,m,sp)
      toc("gridSolveEE")
    else:
      s.gridSolveOO(r,x,m,sp)
      toc("gridSolveOO")
    sp.calls = 1
    sp.seconds = getElapsedTime()
    let flops = (s.g.len*4*72+60)*r.l.nEven*sp.iterations
    sp.flops = flops.float
    if sp0.verbosity>0:
      echo "solveXX(Grid): ", sp.getStats
  sp.iterationsMax = sp.iterations
  sp.r2.push 0.0
  sp0.addStats(sp)
  toc("end solveXX")

proc solveEE*(s: Staggered; r,x: Field; m: SomeNumber; sp0: var SolverParams) =
  solveXX(s, r, x, m, sp0, parEven=true)

proc solveOO*(s: Staggered; r,x: Field; m: SomeNumber; sp0: var SolverParams) =
  solveXX(s, r, x, m, sp0, parEven=false)

# right-preconditioned
proc solveReconR(s:Staggered; x,b:Field; m:SomeNumber; sp: var SolverParams;
                 b2e,b2o: float) =
  tic("solveReconR")
  let b2 = b2e + b2o
  let r2stop = sp.r2req * b2
  let r2stop2 = 0.5 * r2stop
  var r2stope = (if b2o <= r2stop2: r2stop-b2o else: r2stop2)
  if b2e > r2stope:
    var y = newOneOf(x)
    threads:
      y := 0
    sp.r2req = r2stope / b2e
    toc("setup")
    s.solveEE(y, b, m, sp)
    toc("solveEE")
    threads:
      y.even *= 4
      threadBarrier()
      s.Ddag(x, y, m)
    toc("reconstruct")
    return
  var r2stopo = (if b2e <= r2stop2: r2stop-b2e else: r2stop2)
  if b2o > r2stopo:
    var y = newOneOf(x)
    threads:
      y := 0
    sp.r2req = r2stopo / b2o
    toc("setup")
    s.solveOO(y, b, m, sp)
    toc("solveOO")
    threads:
      y.odd *= 4
      threadBarrier()
      s.Ddag(x, y, m)
    toc("reconstruct")
    return

# left-preconditioned with odd reconstruction
proc solveReconL(s:Staggered; x,b:Field; m:SomeNumber; sp: var SolverParams;
                 b2e,b2o: float) =
  tic("solveReconL")
  #if b2e == 0.0 or b2o == 0.0:
  #solveR(s, y, r, m, sp, r2e, r2o)
  var d = newOneOf(b)
  var d2e = 0.0
  threads:
    #s.eoReduce(t, x, m)
    s.Ddag(d, b, m)
    x := 0
    threadBarrier()
    let
      d2et = d.even.norm2
      #d2ot = d.odd.norm2
    threadMaster:
      d2e = d2et
      #d2o = d2ot
  #echo "d2e: ", d2e, "  d2o: ", d2o
  sp.r2req = 0.99 * sp.r2req * (b2e+b2o) * m*m / d2e
  toc("setup")
  s.solveEE(x, d, m, sp)
  toc("solveEE")
  #echo "solveReconL ", d2e, "  ", sp.r2req
  threads:
    x.even *= 4
    threadBarrier()
    s.eoReconstruct(x, b, m)
    threadBarrier()
  toc("reconstruct")

# solveInner {init, run, free}
proc solveInner(s:Staggered; x,b:Field; m:SomeNumber; sp: var SolverParams;
                b2e,b2o: float) =
  let b2 = b2e + b2o
  let r2stop = sp.r2req * b2
  let r2stop2 = 0.5 * r2stop
  var r2stope = (if b2o <= r2stop2: r2stop-b2o else: r2stop2)
  var r2stopo = (if b2e <= r2stop2: r2stop-b2e else: r2stop2)
  if b2e <= r2stope or b2o <= r2stopo or m == 0.0:
    solveReconR(s, x, b, m, sp, b2e, b2o)
  else:
    solveReconL(s, x, b, m, sp, b2e, b2o)
    #solveReconR(s, x, b, m, sp, b2e, b2o)

proc solve*(s:Staggered; x,b:Field; m:SomeNumber; sp0: var SolverParams) =
  tic()
  var b2,r2e,r2o,r2 = 0.0
  threads:
    let
      b2t = b.norm2
    threadMaster:
      b2 = b2t
  let r2stop = sp0.r2req * b2
  var r = newOneOf(b)
  if sp0.usePrevSoln:
    threads:
      s.D(r, x, m)
      threadBarrier()
      r := b - r
  else:
    threads:
      x := 0
      r := b
  threads:
    let
      r2et = r.even.norm2
      r2ot = r.odd.norm2
    threadMaster:
      r2e = r2et
      r2o = r2ot
  r2 = r2e + r2o
  if sp0.verbosity>1:
    echo &"stagSolve b2: {b2:.6g}  r2: {r2/b2:.6g}  r2stop: {r2stop:.6g}"

  var y = newOneOf(x)
  var sp = sp0
  sp.resetStats()
  dec sp.verbosity
  sp.usePrevSoln = false
  while r2 > r2stop:
    sp.maxits = sp0.maxits - sp.iterations
    if sp.maxits <= 0: break
    sp.r2req = r2stop / r2;

    solveInner(s, y, r, m, sp, r2e, r2o)

    threads:
      x += y
      threadBarrier()
      s.D(r, x, m)
      threadBarrier()
      r := b - r
      threadBarrier()
      let
        r2et = r.even.norm2
        r2ot = r.odd.norm2
      threadMaster:
        r2e = r2et
        r2o = r2ot
    r2 = r2e + r2o
    if sp.verbosity>0:
      echo "stagSolve r2/b2: ", r2/b2

  sp.r2.init r2/b2
  sp.calls = 1
  sp.seconds = getElapsedTime()
  sp.flops += float((s.g.len*4*72+24)*x.l.nEven) # ???
  if sp0.verbosity>0:
    #let its = sp.iterations
    #let s = sp.seconds
    #let f = sp.flops
    #let gf = 1e-9*f/s
    #echo "op time: ", top
    echo "stagSolve: ", sp.getStats
  sp0.addStats(sp)

# multimass (trivial version with multiple single mass calls for now)
proc solve*(s: Staggered; r: seq[Field]; x: Field; m: seq[float];
            sp: seq[SolverParams]) =
  let n = m.len
  doAssert(r.len == n)
  doAssert(sp.len == n)
  for i in 0..<n:
    solve(s, r[i], x, m[i], sp[i])

proc solve*(s: Staggered; r: seq[Field]; x: Field; m: seq[float];
            sp: var SolverParams) =
  let n = m.len
  doAssert(r.len == n)
  doAssert(sp.len == n)
  for i in 0..<n:
    solve(s, r[i], x, m[i], sp)

proc solve*(s:Staggered; r,x:Field; m:SomeNumber; res:float;
            cpuonly = false; sloppySolve = SloppyNone) =
  var sp = newSolverParams()
  sp.r2req = res
  #sp.maxits = 1000
  sp.verbosity = 1
  sp.subsetName = "even"
  sp.sloppySolve = sloppySolve
  if cpuonly:
    sp.backend = sbQex
  solve(s, r, x, m, sp)


type S2oa*[T] = object
  s: T
  s2: T
  m: float
proc apply*(oa: S2oa; Dt,t: Field) =
  threadBarrier()
  oa.s.D(Dt, t, oa.m)
  threadBarrier()
proc preconditioner*(oa: S2oa; z: Field; gs: GcrState) =
  threadBarrier()
  z := gs.r
  #oa.s2.Ddag(z, gs.r, oa.m)
  threadBarrier()
  #for e in z:
  #  #let r2 = gs.r[e].norm2
  #  let x2 = gs.x[e].norm2
  #  let s = 1.0/(1e10+x2)
  #  z[e] := asReal(s)*gs.r[e]
proc solve2*(s:Staggered; x,b:Field; m:SomeNumber;
             sp:var SolverParams; s2: Staggered) =
  #var t = newOneOf(r)
  var gcr = newGcrState(x=x, b=b)
  #proc op(Dt,t: Field) =
  #  threadBarrier()
  #  s.D(Dt, t, m)
  #  threadBarrier()
  #proc id(z: Field, c: type(gcr)) =
  #  assign(z, c.r)
  #var oa = (apply: op, preconditioner: id)
  let s2oa = S2oa[type(s)](s: s, s2: s2, m: m)
  gcr.solve(s2oa, sp)

proc solve2*(s:Staggered; r,x:Field; m:SomeNumber; res:float; s2: Staggered;
             cpuonly = false) =
  var sp = newSolverParams()
  sp.r2req = res
  #sp.maxits = 1000
  sp.verbosity = 1
  solve2(s, r, x, m, sp, s2, cpuonly)


when isMainModule:
  import qex
  qexInit()
  #var defaultLat = [4,4,4,4]
  #var defaultLat = [8,8,8,8]
  var defaultLat = @[8,8,8,8]
  #var defaultLat = @[12,12,12,12]
  defaultSetup()
  var v1 = lo.ColorVector()
  var v2 = lo.ColorVector()
  var v3 = lo.ColorVector()
  var r = lo.ColorVector()
  var rs = newRNGField(RngMilc6, lo, intParam("seed", 987654321).uint64)

  if fn == "":
    var warm = floatParam("warm", 0.15)
    threads:
      #g.random rs
      g.warm warm, rs
  let plaq = g.plaq
  echo "plaq: ", plaq.sum, " ", plaq

  threads:
    g.setBC
    g.stagPhase
    #v1 := 0
    v1.gaussian rs
  #if myRank==0:
  #  v1{0}[0] := 1
  #  #v1{2*1024}[0] := 1
  echo v1.norm2

  var s = newStag(g)
  var m = floatParam("m", 0.01)
  var sp = newSolverParams()
  sp.verbosity = intParam("verb", 2)
  sp.subset.layoutSubset(lo, "all")
  sp.maxits = int(1e9/lo.physVol.float)
  sp.r2req = floatParam("rsq", 1e-12)

  proc test =
    v2 := 0
    s.solve(v2, v1, m, sp)
    threads:
      s.D(v3, v2, m)
      v1 := 0
    resetTimers()
    s.solve(v1, v3, m, sp)
    threads:
      r := v1 - v2
      echo "err2: ", r.norm2
    threads:
      v1 := 0
    resetTimers()
    #precon = true
    s.solve(v1, v3, m, sp)
    threads:
      r := v1 - v2
      echo "err2: ", r.norm2

  block:
    v1 := 0
    let p = lo.rankIndex([0,0,0,0])
    if myRank==p.rank:
      v1{p.index}[0] := 1
    echo "even point"
    test()
    echo sp.getStats()
#[
  block:
    v1 := 0
    let p = lo.rankIndex([0,0,0,1])
    if myRank==p.rank:
      v1{p.index}[0] := 1
    echo "odd point"
    test()
    echo sp.getStats()

  block:
    v1.gaussian rs
    echo "random"
    test()
    echo sp.getStats()
]#
  if intParam("timers", 0)!=0:
    echoTimers()

  qexFinalize()
