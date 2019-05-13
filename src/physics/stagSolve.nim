import times
import base
import layout
import field
import stagD
export stagD
import solvers/gcr
import maths

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

# Late import to avoid problem with circular dependence.
when defined(qudaDir):
  import quda/qudaWrapper

proc solve*(s:Staggered; r,x:Field; m:SomeNumber; sp0: var SolverParams;
            cpuonly = false) =
  ## When QUDA is available, we use QUDA unless `cpuonly` is true.
  var sp = sp0
  var its = 0
  var t = newOneOf(r)
  var top = 0.0
  proc op(a,b: Field) =
    threadBarrier()
    if threadNum==0: top -= epochTime()
    stagD2ee(s.se, s.so, a, s.g, b, m*m)
    if threadNum==0: top += epochTime()
    #threadBarrier()
  var oa = (apply: op)
  threads:
    #echo "x2: ", x.norm2
    s.eoReduce(t, x, m)
    #echo "te2: ", t.even.norm2
    r := 0
  let t0 = epochTime()
  when defined(qudaDir):
    if not cpuonly:
      let tquda0 = epochTime()
      s.qudaSolveEE(r,t,m,sp)
      let tquda1 = epochTime()
      let squda = tquda1-tquda0
      its = sp.finalIterations
      sp.finalIterations = 0
      let flopsquda = (s.g.len*4*72+60)*r.l.nEven*its
      if sp0.verbosity>0:
        echo "quda time: ",squda,"  Gflops: ", 1e-9*flopsquda.float/squda
    # After QUDA, we still run through our solver.
  #cgSolve(r, t, oa, sp)
  var cg = newCgState(r, t)
  sp.subset.layoutSubset(r.l, sp.subsetName)
  if sp.subsetName=="all": sp.subset.layoutSubset(r.l, "even")
  cg.solve(oa, sp)
  its += sp.finalIterations
  let t1 = epochTime()
  #var u = newOneOf(r)
  #stagD2ee(s.se, s.so, u, s.g, r, m*m)
  #u := t - u
  #echo "u2: ", u.even.norm2
  threads:
    r[s.se.sub] := 4*r
    threadBarrier()
    s.eoReconstruct(r, x, m)
  sp0.finalIterations += its
  let secs = t1-t0
  let flops = (s.g.len*4*72+60)*r.l.nEven*sp.finalIterations
  if sp0.verbosity>0:
    echo "op time: ", top
    echo "solve time: ", secs, "  Gflops: ", 1e-9*flops.float/secs
proc solve*(s:Staggered; r,x:Field; m:SomeNumber; res:float;
            cpuonly = false; sloppySolve = SloppyNone) =
  var sp = initSolverParams()
  sp.r2req = res
  #sp.maxits = 1000
  sp.verbosity = 1
  sp.subsetName = "even"
  sp.sloppySolve = sloppySolve
  solve(s, r, x, m, sp, cpuonly)


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
  var sp = initSolverParams()
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
  var r = lo.ColorVector()
  var rs = newRNGField(RngMilc6, lo, intParam("seed", 987654321).uint64)
  threads:
    #g.random rs
    g.setBC
    g.stagPhase
    v1 := 0
    #for e in v1:
    #  template x(d:int):untyped = lo.vcoords(d,e)
    #  v1[e][0].re := foldl(x, 4, a*10+b)
    #  #echo v1[e][0]
  #echo v1.norm2
  var g2 = lo.newGauge
  for mu in 0..<g2.len:
    for s in 0..<lo.nSites:
      if lo.coords[mu][s] mod 2 == 1:
        g2[mu]{s} := 0
      else:
        g2[mu]{s} := g[mu]{s}
  var s2 = newStag(g2)
  if myRank==0:
    v1{0}[0] := 1
    #v1{2*1024}[0] := 1
  echo v1.norm2
  #var gs = lo.newGaugeS
  #for i in 0..<gs.len: gs[i] := g[i]
  var s = newStag(g)
  var m = 0.1
  threads:
    v2 := 0
    echo v2.norm2
    threadBarrier()
    s.D(v2, v1, m)
    threadBarrier()
    #echoAll v2
    echo v2.norm2
  #echo v2
  var sp = initSolverParams()
  sp.subset.layoutSubset(lo, "all")
  sp.maxits = int(1e9/lo.physVol.float)
  s.solve(v2, v1, m, sp)
  #resetTimers()
  s.solve2(v2, v1, m, sp, s2)
  threads:
    echo "v2: ", v2.norm2
    echo "v2.even: ", v2.even.norm2
    echo "v2.odd: ", v2.odd.norm2
    s.D(r, v2, m)
    threadBarrier()
    r := v1 - r
    threadBarrier()
    echo r.norm2
  #echo v2
  #echoTimers()
