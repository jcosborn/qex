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
when defined(gridDir):
  import grid/GridDefs

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

proc solveEE*(s: Staggered; r,x: Field; m: SomeNumber; sp0: var SolverParams) =
  tic()
  var sp = sp0
  sp.resetStats()
  dec sp.verbosity
  threads:
    r := 0
  case sp.backend
  of sbQex:
    tic()
    proc op(a,b: Field) =
      tic()
      threadBarrier()
      stagD2ee(s.se, s.so, a, s.g, b, m*m)
      #threadBarrier()
      toc("stagD2ee")
    var oa = (apply: op)
    var cg = newCgState(r, x)
    sp.subset.layoutSubset(r.l, sp.subsetName)
    if sp.subsetName=="all": sp.subset.layoutSubset(r.l, "even")
    cg.solve(oa, sp)
    toc("cg.solve")
    sp.calls = 1
    sp.seconds = getElapsedTime()
    let flops = (s.g.len*4*72+60)*r.l.nEven*sp.iterations
    sp.flops = flops.float
    if sp0.verbosity>0:
      echo "solveEE: ", sp.getStats
  of sbQuda:
    tic()
    s.qudaSolveEE(r,x,m,sp)
    toc("qudaSolveEE")
    sp.calls = 1
    sp.seconds = getElapsedTime()
    let flopsquda = (s.g.len*4*72+60)*r.l.nEven*sp.iterations
    sp.flops = flopsquda.float
    if sp0.verbosity>0:
      echo "solveEE(QUDA): ", sp.getStats
  of sbGrid:
    tic()
    s.qudaSolveEE(r,x,m,sp)
    toc("qudaSolveEE")
    sp.calls = 1
    sp.seconds = getElapsedTime()
    let flopsquda = (s.g.len*4*72+60)*r.l.nEven*sp.iterations
    sp.flops = flopsquda.float
    if sp0.verbosity>0:
      echo "solveEE(QUDA): ", sp.getStats
  #[
  else:  # remove?
    tic()
    proc op(a,b: Field) =
      tic()
      threadBarrier()
      s.D(a, b, m*m)
      threadBarrier()
      a.odd += (1-m*m)*b
      toc("stagDee2")
    #var oa = (apply: op)
    #var cg = newCgState(r, x)
    sp.subsetName = "all"
    sp.subset.layoutSubset(r.l, sp.subsetName)
    x.odd := 0
    #cg.solve(oa, sp)
    bicgstabSolve(r, x, op, sp)
    toc("bicg.solve")
    r.even := 0.25*r
    r.odd := 0
    sp.calls = 1
    sp.seconds = getElapsedTime()
    let flops = (s.g.len*4*72+60)*r.l.nEven*sp.iterations
    sp.flops = flops.float
    if sp0.verbosity>0:
      echo "solveEE: ", sp.getStats
  ]#
  sp.iterationsMax = sp.iterations
  sp.r2.push 0.0
  sp0.addStats(sp)

# solveInner {init, run, free}

proc solve*(s:Staggered; x,b:Field; m:SomeNumber; sp0: var SolverParams) =
  tic()
  var c = newOneOf(b)
  var b2,r2 = 0.0
  if sp0.usePrevSoln:
    threads:
      s.D(c, x, m)
      threadBarrier()
      c := b - c
      let
        b2t = b.norm2
        c2t = c.norm2
      threadMaster:
        b2 = b2t
        r2 = c2t
  else:
    threads:
      x := 0
      c := b
      let
        b2t = b.norm2
      threadMaster:
        b2 = b2t
        r2 = b2t
  var r2stop = sp0.r2req * b2
  if sp0.verbosity>1:
    echo &"stagSolve b2: {b2:.6g}  r2: {r2/b2:.6g}  r2stop: {r2stop:.6g}"
  var d = newOneOf(b)
  var y = newOneOf(x)
  var sp = sp0
  sp.resetStats()
  dec sp.verbosity
  sp.usePrevSoln = false
  while r2 > r2stop:
    sp.maxits = sp0.maxits - sp.iterations
    if sp.maxits <= 0: break
    var d2e = 0.0
    threads:
      #s.eoReduce(t, x, m)
      s.Ddag(d, c, m)
      y := 0
      threadBarrier()
      let
        d2et = d.even.norm2
        #d2ot = d.odd.norm2
      threadMaster:
        d2e = d2et
        #d2o = d2ot
    #echo "d2e: ", d2e, "  d2o: ", d2o
    sp.r2req = 0.99 * sp0.r2req * b2 * m*m / d2e
    toc("setup")
    s.solveEE(y, d, m, sp)
    toc("solveEE")
    threads:
      y.even := 4*y
      threadBarrier()
      s.eoReconstruct(y, c, m)
      threadBarrier()
      x += y
      threadBarrier()
      s.D(c, x, m)
      threadBarrier()
      c := b - c
      let
        c2t = c.norm2
      threadMaster:
        r2 = c2t
    toc("reconstruct")
    if sp.verbosity>0:
      echo "stagSolve r2: ", r2/b2

  #echo "r2: ", r2
  #sp.r2sum = r2/b2
  #sp.r2max = r2/b2
  sp.r2.init r2/b2
  sp.calls = 1
  sp.seconds = getElapsedTime()
  sp.flops += float((s.g.len*4*72+24)*x.l.nEven)
  if sp0.verbosity>0:
    #let its = sp.iterations
    #let s = sp.seconds
    #let f = sp.flops
    #let gf = 1e-9*f/s
    #echo "op time: ", top
    echo "stagSolve: ", sp.getStats
  sp0.addStats(sp)

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

  threads:
    g.random rs
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

  block:
    v1 := 0
    let p = lo.rankIndex([0,0,0,0])
    if myRank==p.rank:
      v1{p.index}[0] := 1
    echo "even point"
    test()
    echo sp.getStats()

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

  if intParam("timers", 0)!=0:
    echoTimers()
