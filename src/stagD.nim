import os
import strUtils
import qex
import qcdTypes
import stdUtils
import times
import cg
import types
import matrixConcept
import profile
import metaUtils

#{.emit:"#define memset(a,b,c)".}

type StaggeredD*[T] = object
  sf*:array[4,ShiftB[T]]
  sb*:array[4,ShiftB[T]]
  sub*:string
  subset*:Subset
type Staggered*[G,T] = object
  se*,so*:StaggeredD[T]
  g*:seq[G]

proc initSB(s:var ShiftB; x:Field; dir,ln:int; sub="all") =
  initShiftB(s, x.l, type(x[0]), dir, ln, sub)

template initStagDT*(l:var Layout; T:typedesc; ss:string):expr =
  var sd:StaggeredD[T]
  for mu in 0..<4:
    #echoRank "init: ", mu
    initShiftB(sd.sf[mu], l, T, mu, 1, ss)
    initShiftB(sd.sb[mu], l, T, mu,-1, ss)
  sd.sub = ss
  sd.subset.layoutSubset(l, ss)
  sd

#proc initStagD*[T](l:var Layout
proc initStagD*(x:Field; sub:string):auto =
  result = initStagDT(x.l, type(x[0]), sub)
  #type t = type(x[0])
  #var sd:StaggeredD[t]
  #initShiftB(sd.sf[3], x.l, type(x[0]), 3, 1, sub)
  #initSB(sd.sf[3], x, 3, 1, sub)
  #sd

proc stagD*(sd:StaggeredD; r:Field; g:openArray[Field2];
            x:Field; m:SomeNumber; sc:SomeNumber=1.0) =
  #{.emit:"#define memset(a,b,c)".}
  template sf0:expr = sd.sf
  template sb0:expr = sd.sb
  #mixin imsub
  tic()
  let sch = 0.5*sc
  for mu in 0..<4:
    #startSB(sf0[mu], x[ix])
    startSB(sf0[mu], sch*x[ix])
  toc("startShiftF")
  for mu in 0..<4:
    startSB(sb0[mu], g[mu][ix].adj*x[ix])
  toc("startShiftB")
  for ir in r[sd.subset]:
    var rir{.noInit.}:type(r[ir])
    #r[ir] := m * x[ir]
    mul(rir, m, x[ir])
    for mu in 0..<4:
      #localSB(sf[mu], ir, r[ir] += g[mu][ir]*it, x[ix])
      #localSB(sf[mu], ir, imadd(r[ir], g[mu][ir], it), x[ix])
      #localSB(sb[mu], ir, isub(r[ir], it), g[mu][ix].adj*x[ix])
      localSB(sf0[mu], ir, imadd(rir, g[mu][ir], it), sch*x[ix])
      localSB(sb0[mu], ir, imsub(rir, sch, it), g[mu][ix].adj*x[ix])
    assign(r[ir], rir)
  toc("local", flops=(6+(4*(6+72+66+12)))*sd.subset.len)
  for mu in 0..<4:
    boundarySB(sf0[mu], imadd(r[ir], g[mu][ir], it))
  toc("boundaryF")
  for mu in 0..<4:
    boundarySB(sb0[mu], imsub(r[ir], sch, it))
  #threadBarrier()
  toc("boundaryB")
  #{.emit:"#undef memset".}

proc stagD2ee*(sde,sdo:StaggeredD; r:Field; g:openArray[Field2];
               x:Field; m2:SomeNumber) =
  var t{.global.}:type(x)
  if t==nil:
    threadBarrier()
    if threadNum==0:
      t = newOneOf(x)
    threadBarrier()
  #threadBarrier()
  stagD(sdo, t, g, x, 0.0)
  threadBarrier()
  stagD(sde, r, g, t, 0.0)
  #threadBarrier()
  #r[sde.sub] := m2*x - r
  for ir in r[sde.subset]:
    msubVSVV(r[ir], m2, x[ir], r[ir])

proc setBC*(g:openArray[Field]) =
  let gt = g[3]
  tfor i, 0..<gt.l.nSites:
    #let e = i div gt.l.nSitesInner
    if gt.l.coords[3][i] == gt.l.physGeom[3]-1:
      gt{i} *= -1
      #echoAll isMatrix(gt{i})
      #echoAll i, " ", gt[e][0,0]
proc stagPhase*(g:openArray[Field]) =
  const phases = [8,9,11,0]
  let l = g[0].l
  for mu in 0..<4:
    tfor i, 0..<l.nSites:
      var s = 0
      for k in 0..<4:
        s += (phases[mu] shr k) and l.coords[k][i].int
      if (s and 1)==1:
        g[mu]{i} *= -1
        #echoAll i, " ", gt[e][0,0]

proc newStag*[G](g:openArray[G]):auto =
  var l = g[0].l
  template t:expr =
    type(l.ColorVector()[0])
  var r:Staggered[G,t]
  r.se = initStagDT(l, t, "even")
  r.so = initStagDT(l, t, "odd")
  r.g = @g
  r

proc D*(s:Staggered; r,x:Field; m:SomeNumber) =
  stagD(s.se, r, s.g, x, m)
  stagD(s.so, r, s.g, x, m)
proc Ddag*(s:Staggered; r,x:Field; m:SomeNumber) =
  stagD(s.se, r, s.g, x, m, -1)
  stagD(s.so, r, s.g, x, m, -1)
proc eoReduce*(s:Staggered; r,b:Field; m:SomeNumber) =
  # r.even = (D^+ b).even
  #dump: "b.even.norm2"
  #dump: "b.odd.norm2"
  stagD(s.se, r, s.g, b, m, -1)
  #dump: r.even.norm2
  #dump: r.odd.norm2
proc eoReconstruct*(s:Staggered; r,b:Field; m:SomeNumber) =
  # r.odd = (b.odd - Doe r.even)/m
  stagD(s.so, r, s.g, r, 0.0, -1.0/m)
  r.odd += b/m

proc solve*(s:Staggered; r,x:Field; m:SomeNumber; res:float) =
  var sp:SolverParams
  sp.r2req = res
  sp.maxits = 1000
  sp.verbosity = 1
  sp.subset.layoutSubset(r.l, "even")
  var t = newOneOf(r)
  var top = 0.0
  proc op(a,b:Field) =
    threadBarrier()
    if threadNum==0: top -= epochTime()
    stagD2ee(s.se, s.so, a, s.g, b, m*m)
    if threadNum==0: top += epochTime()
    #threadBarrier()
  threads:
    #echo "x2: ", x.norm2
    s.eoReduce(t, x, m)
    #echo "te2: ", t.even.norm2
  let t0 = epochTime()
  cgSolve(r, t, op, sp)
  let t1 = epochTime()
  threads:
    s.eoReconstruct(r, x, m)
  let secs = t1-t0
  let flops = (1152+60)*r.l.nEven*sp.finalIterations
  echo top
  echo "solve time: ", secs, "  Gflops: ", 1e-9*flops.float/secs
proc solve2*(s:Staggered; r,x:Field; m:SomeNumber; res:float) =
  var sp:SolverParams
  sp.r2req = res
  sp.maxits = 100
  sp.verbosity = 1
  sp.subset.layoutSubset(r.l, "all")
  var t = newOneOf(r)
  proc op(a,b:Field) =
    #stagD2ee(s.se, s.so, r, s.g, x, m*m)
    threadBarrier()
    s.Ddag(t, b, m)
    #threadBarrier()
    #echo t.norm2
    s.D(a, t, m)
    #a := b
    #threadBarrier()
  #echo r.norm2
  cgSolve(r, x, op, sp)
  threads:
    s.Ddag(t, r, m)
    r := t

template foldl*(f,n,op:untyped):expr =
  var r:type(f(0))
  r = f(0)
  for i in 1..<n:
    let
      a {.inject.} = r
      b {.inject.} = f(i)
    r = op
  r

when isMainModule:
  qexInit()
  echo "rank ", myRank, "/", nRanks
  let cp = commandLineParams()
  #var lat = [4,4,4,4]
  var lat = [8,8,8,8]
  #var lat = [8,8,8,16]
  #var lat = [8,8,16,16]
  #var lat = [16,16,16,8]
  #var lat = [16,16,16,16]
  #var lat = [16,16,16,32]
  if cp.len>0:
    for i in 0..<lat.len:
      lat[i] = (if i<cp.len: parseInt(cp[i]) else: lat[i-1])
  var lo = newLayout(lat)
  var g:array[4,type(lo.ColorMatrix())]
  for i in 0..<4:
    g[i] = lo.ColorMatrix()
    threads: g[i] := 1
  var v1 = lo.ColorVector()
  var v2 = lo.ColorVector()
  #g.loadGauge("l88.scidac")
  var sdAll = initStagD(v1, "all")
  var sdEven = initStagD(v1, "even")
  var sdOdd = initStagD(v1, "odd")
  var s = newStag(@g)
  var m = 0.1
  threads:
    g.setBC
    threadBarrier()
    for i in 0..<4:
      echo g[i].norm2
    threadBarrier()
    g.stagPhase
    threadBarrier()
    for i in 0..<4:
      echo g[i].norm2
    v1 := 0
    #v2 := 1
    if myRank==0 and threadNum==0:
      v1{0}[0] := 1
    threadBarrier()
    echo v1.norm2

    stagD(sdAll, v2, g, v1, m)
    threadBarrier()
    echo v2.norm2
    #echo v2
    s.D(v2, v1, m)
    threadBarrier()
    echo v2.norm2

    for e in v1:
      template x(d:int):expr = lo.vcoords(d,e)
      v1[e][0].re := foldl(x, 4, a*10+b)
      #echo v1[e][0]
    threadBarrier()
    stagD(sdAll, v2, g, v1, 0.5)
    echo v1[0][0]
    echo v2[0][0]

  let nrep = int(1e7/lo.physVol.float)
  #let nrep = int(1e9/lo.physVol.float)
  #let nrep = 1
  template makeBench(name:untyped; bar:bool):untyped =
    proc name(sd:var any, ss="all") =
      resetTimers()
      var t0 = epochTime()
      threads(sd):
        for rep in 1..nrep:
          stagD(sd, v2, g, v1, 0.5)
          when bar: threadBarrier()
      var t1 = epochTime()
      let dt = t1-t0
      #var vol = lo.physVol.float
      var vol = lo.nSites.float
      if sd.sub != "all": vol *= 0.5
      let flops = (6.0+8.0*72.0) * vol
      echo ss & "secs: ", dt, "  mf: ", (nrep.float*flops)/(1e6*dt)
      echoTimers()
  makeBench(bench, false)
  makeBench(benchB, true)
  bench(sdAll, "all  ")
  benchB(sdAll, "all  ")
  bench(sdEven, "even ")
  benchB(sdEven, "even ")
  bench(sdOdd, "odd  ")
  benchB(sdOdd, "odd  ")
  proc benchEO() =
    var t0 = epochTime()
    threads:
      for rep in 1..nrep:
        stagD2ee(sdEven, sdOdd, v2, g, v1, 0.1)
    var t1 = epochTime()
    let dt = t1-t0
    #var vol = 0.5 * lo.physVol.float
    var vol = 0.5 * lo.nSites.float
    let flops = (6.0+2.0*8.0*72.0) * vol
    echo "EO   secs: ", dt, "  mf: ", (nrep.float*flops)/(1e6*dt)
  benchEO()
  qexFinalize()
