static:
  echo "Starting wilsonD imports: ", staticExec("date")

import ../base/globals
setForceInline(false)
#setForceInline(true)
#setStaticUnroll(false)
#setStaticUnroll(true)
#setNoAlias(false)
#setNoAlias(true)

import os
import times
import strUtils

import base
import layout
import field
import qcdTypes
#import stdUtils
import solvers/cg
export cg
import solvers/bicgstab
#import types
#import profile
#import metaUtils
import gauge/gaugeUtils

static:
  echo "Starting wilsonD: ", staticExec("date")

type WilsonD*[T] = object
  sf*: seq[ShiftB[T]]
  sb*: seq[ShiftB[T]]
  sub*: string
  subset*: Subset
type Wilson*[G,T] = object
  se*,so*: WilsonD[T]
  g*: seq[G]

template initWilsonDT*(l: var Layout; T: typedesc; ss: string): untyped =
  var sd: WilsonD[T]
  sd.sf.newSeq(4)
  sd.sb.newSeq(4)
  for mu in 0..<4:
    initShiftB(sd.sf[mu], l, T, mu, 1, ss)
    initShiftB(sd.sb[mu], l, T, mu,-1, ss)
  sd.sub = ss
  sd.subset.layoutSubset(l, ss)
  sd

proc initWilsonD*(x: Field; sub: string): auto =
  result = initWilsonDT(x.l, type(spproj1p(x[0])), sub)

template optimizeAstX(x: untyped): untyped = x #optimizeAst(x)
# normalized to 2*D_w
template wilsonDP*(sd: WilsonD; r: Field; g: openArray[Field2];
                   x: Field3; expFlops: int; exp: untyped) =
  static:
    echo "Starting wilsonDP: ", staticExec("date")
  tic()
  optimizeAstX:
    startSB(sd.sf[0], spproj1p(x[ix]))
    startSB(sd.sf[1], spproj2p(x[ix]))
    startSB(sd.sf[2], spproj3p(x[ix]))
    startSB(sd.sf[3], spproj4p(x[ix]))
  toc("startShiftF")
  optimizeAstX:
    startSB(sd.sb[0], g[0][ix].adj*spproj1m(x[ix]))
    startSB(sd.sb[1], g[1][ix].adj*spproj2m(x[ix]))
    startSB(sd.sb[2], g[2][ix].adj*spproj3m(x[ix]))
    startSB(sd.sb[3], g[3][ix].adj*spproj4m(x[ix]))
  toc("startShiftB")
  optimizeAstX:
    for ir{.inject.} in r[sd.subset]:
      var rir{.inject,noInit.}: type(load1(r[ir]))
      exp
      localSB(sd.sf[0], ir, rir-=sprecon1p(g[0][ir]*it), spproj1p(x[ix]))
      localSB(sd.sf[1], ir, rir-=sprecon2p(g[1][ir]*it), spproj2p(x[ix]))
      localSB(sd.sf[2], ir, rir-=sprecon3p(g[2][ir]*it), spproj3p(x[ix]))
      localSB(sd.sf[3], ir, rir-=sprecon4p(g[3][ir]*it), spproj4p(x[ix]))
      localSB(sd.sb[0], ir, rir-=sprecon1m(it), g[0][ix].adj*spproj1m(x[ix]))
      localSB(sd.sb[1], ir, rir-=sprecon2m(it), g[1][ix].adj*spproj2m(x[ix]))
      localSB(sd.sb[2], ir, rir-=sprecon3m(it), g[2][ix].adj*spproj3m(x[ix]))
      localSB(sd.sb[3], ir, rir-=sprecon4m(it), g[3][ix].adj*spproj4m(x[ix]))
      assign(r[ir], rir)
  toc("local", flops=(expFlops+2*g.len*(12+2*66+24))*sd.subset.len)
  optimizeAstX:
    #template bsb0(ir0,it: typed): untyped =
    #  r[ir0] -= sprecon1p(g[0][ir0]*it)
    #boundarySB2(sd.sf[0], bsb0)
    boundarySB(sd.sf[0], r[ir]-=sprecon1p(g[0][ir]*it))
    boundarySB(sd.sf[1], r[ir]-=sprecon2p(g[1][ir]*it))
    boundarySB(sd.sf[2], r[ir]-=sprecon3p(g[2][ir]*it))
    boundarySB(sd.sf[3], r[ir]-=sprecon4p(g[3][ir]*it))
  toc("boundaryF")
  optimizeAstX:
    boundarySB(sd.sb[0], r[ir]-=sprecon1m(it))
    boundarySB(sd.sb[1], r[ir]-=sprecon2m(it))
    boundarySB(sd.sb[2], r[ir]-=sprecon3m(it))
    boundarySB(sd.sb[3], r[ir]-=sprecon4m(it))
  #threadBarrier()
  toc("boundaryB")
  static:
    echo "Done wilsonDP: ", staticExec("date")

# normalized to 2*D_w
template wilsonDM*(sd: WilsonD; r: Field; g: openArray[Field2];
                   x: Field3; expFlops: int; exp: untyped) =
  static:
    echo "Starting wilsonDM: ", staticExec("date")
  tic()
  optimizeAstX:
    startSB(sd.sf[0], spproj1m(x[ix]))
    startSB(sd.sf[1], spproj2m(x[ix]))
    startSB(sd.sf[2], spproj3m(x[ix]))
    startSB(sd.sf[3], spproj4m(x[ix]))
  static:
    echo " startShiftF: ", staticExec("date")
  toc("startShiftF")
  optimizeAstX:
    startSB(sd.sb[0], g[0][ix].adj*spproj1p(x[ix]))
    startSB(sd.sb[1], g[1][ix].adj*spproj2p(x[ix]))
    startSB(sd.sb[2], g[2][ix].adj*spproj3p(x[ix]))
    startSB(sd.sb[3], g[3][ix].adj*spproj4p(x[ix]))
  static:
    echo " startShiftB: ", staticExec("date")
  toc("startShiftB")
  for ir{.inject.} in r[sd.subset]:
    var rir{.inject,noInit.}: evalType(load1(r[ir]))
    exp
    localSB(sd.sf[0], ir, rir-=sprecon1m(g[0][ir]*it), spproj1m(x[ix]))
    localSB(sd.sf[1], ir, rir-=sprecon2m(g[1][ir]*it), spproj2m(x[ix]))
    localSB(sd.sf[2], ir, rir-=sprecon3m(g[2][ir]*it), spproj3m(x[ix]))
    localSB(sd.sf[3], ir, rir-=sprecon4m(g[3][ir]*it), spproj4m(x[ix]))
    localSB(sd.sb[0], ir, rir-=sprecon1p(it), g[0][ix].adj*spproj1p(x[ix]))
    localSB(sd.sb[1], ir, rir-=sprecon2p(it), g[1][ix].adj*spproj2p(x[ix]))
    localSB(sd.sb[2], ir, rir-=sprecon3p(it), g[2][ix].adj*spproj3p(x[ix]))
    localSB(sd.sb[3], ir, rir-=sprecon4p(it), g[3][ix].adj*spproj4p(x[ix]))
    assign(r[ir], rir)
  static:
    echo " local: ", staticExec("date")
  toc("local", flops=(expFlops+2*g.len*(12+2*66+24))*sd.subset.len)
  optimizeAstX:
    boundarySB(sd.sf[0], r[ir]-=sprecon1m(g[0][ir]*it))
    boundarySB(sd.sf[1], r[ir]-=sprecon2m(g[1][ir]*it))
    boundarySB(sd.sf[2], r[ir]-=sprecon3m(g[2][ir]*it))
    boundarySB(sd.sf[3], r[ir]-=sprecon4m(g[3][ir]*it))
  static:
    echo " boundaryF: ", staticExec("date")
  toc("boundaryF")
  optimizeAstX:
    boundarySB(sd.sb[0], r[ir]-=sprecon1p(it))
    boundarySB(sd.sb[1], r[ir]-=sprecon2p(it))
    boundarySB(sd.sb[2], r[ir]-=sprecon3p(it))
    boundarySB(sd.sb[3], r[ir]-=sprecon4p(it))
  #threadBarrier()
  toc("boundaryB")
  static:
    echo "Done wilsonDM: ", staticExec("date")

# r = m*x + sc*D*x
proc wilsonD*(sd:WilsonD; r:Field; g:openArray[Field2];
              x:Field; m:SomeNumber; sc:SomeNumber=1.0) =
  #wilsonD2(sd, r, g, x, 0, m/(0.5*sc))
  #r[sd.subset] := (0.5*sc)*r
  #wilsonDP2(sd, r, g, x, 6):
  #  #for i in 0..<n:
  #  rir := m*getVec(x[ir], ic)
  wilsonDM(sd, r, g, x, 6):
    rir := (2.0*(4.0 + m)/sc) * x[ir]
    #rir := 0
  r[sd.subset] := (0.5*sc)*r

# r = m*x + sc*D*x
proc wilsonDx*(sd:WilsonD; r:Field; g:openArray[Field2];
              x:Field; m:SomeNumber; sc:SomeNumber=1.0) =
  #wilsonD2(sd, r, g, x, 0, m/(0.5*sc))
  #r[sd.subset] := (0.5*sc)*r
  #wilsonDP2(sd, r, g, x, 6):
  #  #for i in 0..<n:
  #  rir := m*getVec(x[ir], ic)
  wilsonDP(sd, r, g, x, 6):
    rir := (2.0*(4.0 + m)) * x[ir]
    #rir := 0
  r[sd.subset] := 0.5*r

#[
proc wilsonD1*(sd:WilsonD; r:Field; g:openArray[Field2];
             x:Field; m:SomeNumber) =
  wilsonDP(sd, r, g, x, 6):
    rir := 0

proc wilsonD1x*(sd:WilsonD; r:Field; g:openArray[Field2];
              x:Field; m:SomeNumber) =
  wilsonDM(sd, r, g, x, 6):
    rir := 0

# r = m*x + sc*D*x
proc wilsonDb*(sd:WilsonD; r:Field; g:openArray[Field2];
               x:Field; m:SomeNumber; sc:SomeNumber=1.0) =
  #wilsonD2(sd, r, g, x, 0, m/(0.5*sc))
  #r[sd.subset] := (0.5*sc)*r
  #wilsonDP2(sd, r, g, x, 6):
  #  #for i in 0..<n:
  #  rir := m*getVec(x[ir], ic)
  wilsonDP(sd, r, g, x, 6):
    rir := (4.0 + m) * x[ir]
]#

# r = m2 - Deo * Doe
proc wilsonD2ee*(sde,sdo:WilsonD; r:Field; g:openArray[Field2];
               x:Field; m2:SomeNumber) =
  tic()
  var t{.global.}:type(x)
  if t==nil:
    threadBarrier()
    if threadNum==0:
      t = newOneOf(x)
    threadBarrier()
  threadBarrier()
  #wilsonD(sdo, t, g, x, 0.0)
  toc("wilsonD2ee init")
  block:
    wilsonDM(sdo, t, g, x, 0):
      rir := 0
  toc("wilsonD2ee DP")
  threadBarrier()
  toc("wilsonD2ee barrier")
  #wilsonD(sde, r, g, t, 0.0)
  block:
    wilsonDM(sde, r, g, t, 6):
      rir := (-4.0*m2)*x[ir]
  toc("wilsonD2ee DM")
  threadBarrier()
  #r[sde.sub] := m2*x - r
  #for ir in r[sde.subset]:
  #  msubVSVV(r[ir], m2, x[ir], r[ir])
  r[sde.sub] := -0.25*r

# r = m2 - [Deo * Doe]'
proc wilsonD2eex*(sde,sdo:WilsonD; r:Field; g:openArray[Field2];
               x:Field; m2:SomeNumber) =
  tic()
  var t{.global.}:type(x)
  if t==nil:
    threadBarrier()
    if threadNum==0:
      t = newOneOf(x)
    threadBarrier()
  threadBarrier()
  #wilsonD(sdo, t, g, x, 0.0)
  toc("wilsonD2ee init")
  block:
    wilsonDP(sdo, t, g, x, 0):
      rir := 0
  toc("wilsonD2ee DP")
  threadBarrier()
  toc("wilsonD2ee barrier")
  #wilsonD(sde, r, g, t, 0.0)
  block:
    wilsonDP(sde, r, g, t, 6):
      rir := (-4.0*m2)*x[ir]
  toc("wilsonD2ee DM")
  threadBarrier()
  #r[sde.sub] := m2*x - r
  #for ir in r[sde.subset]:
  #  msubVSVV(r[ir], m2, x[ir], r[ir])
  r[sde.sub] := -0.25*r

proc wilsonD2oo*(sde,sdo:WilsonD; r:Field; g:openArray[Field2];
                 x:Field; m2:SomeNumber) =
  wilsonD2ee(sdo, sde, r, g, x, m2)

#[
# r = m2 - Deo * Doe
proc wilsonD2eeN*(sde,sdo:WilsonD; r:Field; g:openArray[Field2];
                x:Field; m2:SomeNumber) =
  block:
    wilsonDPN(sdo, t, g, x, 0):
      rir := 0
  threadBarrier()
  block:
    wilsonDMN(sde, r, g, t, 6):
      rir := (4.0*m2)*x[ir]
]#

proc newWilson*[G,T](g:openArray[G]; v:T):auto =
  var l = g[0].l
  template t: untyped =
    type(spproj1p(v[0]))
  var r:Wilson[G,t]
  r.se = initWilsonDT(l, t, "even")
  r.so = initWilsonDT(l, t, "odd")
  r.g = @g
  r

proc newWilson*[G](g:openArray[G]):auto =
  var l = g[0].l
  template t:untyped =
    type(spproj1p(l.DiracFermion()[0]))
    #SColorVectorV
  var r:Wilson[G,t]
  r.se = initWilsonDT(l, t, "even")
  r.so = initWilsonDT(l, t, "odd")
  r.g = @g
  r

proc newWilsonS*[G](g:openArray[G]):auto =
  var l = g[0].l
  template t:untyped =
    type(spproj1p(l.DiracFermionS()[0]))
    #SColorVectorV
  var r:Wilson[G,t]
  r.se = initWilsonDT(l, t, "even")
  r.so = initWilsonDT(l, t, "odd")
  r.g = @g
  r

proc D*(s:Wilson; r,x:Field; m:SomeNumber) =
  wilsonD(s.se, r, s.g, x, m)
  wilsonD(s.so, r, s.g, x, m)
proc Ddag*(s:Wilson; r,x:Field; m:SomeNumber) =
  wilsonDx(s.se, r, s.g, x, m)
  wilsonDx(s.so, r, s.g, x, m)
proc eoReduce*(s:Wilson; r,b:Field; m:SomeNumber) =
  # r.even = (D^+ b).even
  # re = Dee be - Deo bo
  #dump: "b.even.norm2"
  #dump: "b.odd.norm2"
  wilsonD(s.se, r, s.g, b, m, -1)
  #dump: r.even.norm2
  #dump: r.odd.norm2
proc eoReconstruct*(s:Wilson; r,b:Field; m:SomeNumber) =
  # r.odd = (b.odd - Doe r.even)/m
  let mi = 1.0/(4.0+m)
  wilsonD(s.so, r, s.g, r, 0.0, -mi)
  r.odd += mi*b

proc eoReduceO*(s:Wilson; r,b:Field; m:SomeNumber) =
  # r.odd = (D^+ b).odd
  # ro = Doo bo - Doe be
  wilsonD(s.so, r, s.g, b, m, -1)
proc eoReconstructO*(s:Wilson; r,b:Field; m:SomeNumber) =
  # r.even = (b.even - Deo r.odd)/m
  wilsonD(s.se, r, s.g, r, 0.0, -1.0/m)
  let mi = 1.0/m
  r.even += mi*b

proc solveEO*(s: Wilson; r,x: Field; m: SomeNumber; sp0: var SolverParams) =
  var sp = sp0
  sp.subset.layoutSubset(r.l, sp.subsetName)
  if sp.subsetName=="all": sp.subset.layoutSubset(r.l, "even")
  var t = newOneOf(r)
  var x2 = newOneOf(x)
  var top = 0.0
  let m4 = 4.0 + m
  let m2 = m4*m4
  proc op(a,b:Field) =
    threadBarrier()
    if threadNum==0: top -= epochTime()
    wilsonD2eex(s.se, s.so, t, s.g, b, m2)
    threadBarrier()
    wilsonD2ee(s.se, s.so, a, s.g, t, m2)
    if threadNum==0: top += epochTime()
    threadBarrier()
  threads:
    #echo "x2: ", x.norm2
    s.eoReduce(x2, x, m)
  let t0 = epochTime()
  var oa = (apply: op)
  cgSolve(r, x2, oa, sp)
  threads:
    wilsonD2eex(s.se, s.so, t, s.g, r, m2)
    threadBarrier()
    r[s.se.sub] := t
    threadBarrier()
    s.eoReconstruct(r, x, m)
  let t1 = epochTime()
  let secs = t1-t0
  let flops = (2*2*2*s.g.len*(12+2*66+24)+2*60)*r.l.nEven*sp.finalIterations
  sp0.finalIterations = sp.finalIterations
  sp0.seconds = secs
  echo "op time: ", top
  echo "solve time: ", secs, "  Gflops: ", 1e-9*flops.float/secs
proc solveEO*(s:Wilson; r,x:Field; m:SomeNumber; res:float) =
  var sp = initSolverParams()
  sp.r2req = res
  #sp.maxits = 1000
  sp.verbosity = 1
  solveEO(s, r, x, m, sp)
proc solve*(s:Wilson; r,x:Field; m:SomeNumber; sp0:SolverParams) =
  var sp = sp0
  sp.subset.layoutSubset(r.l, sp.subsetName)
  if sp.subsetName=="all": sp.subset.layoutSubset(r.l, "even")
  var t = newOneOf(r)
  var top = 0.0
  let m4 = 4.0 + m
  let m2 = m4*m4
  proc op(a,b:Field) =
    threadBarrier()
    if threadNum==0: top -= epochTime()
    wilsonD2ee(s.se, s.so, a, s.g, b, m2)
    if threadNum==0: top += epochTime()
    #threadBarrier()
  threads:
    #echo "x2: ", x.norm2
    s.eoReduce(t, x, m)
    #echo "te2: ", t.even.norm2
  let t0 = epochTime()
  bicgstabSolve(r, t, op, sp)
  let t1 = epochTime()
  threads:
    s.eoReconstruct(r, x, m)
  let secs = t1-t0
  #let flops = (s.g.len*4*72+60)*r.l.nEven*sp.finalIterations
  let flops = (2*2*s.g.len*(12+2*66+24)+2*60)*r.l.nEven*sp.finalIterations
  echo "op time: ", top
  echo "solve time: ", secs, "  Gflops: ", 1e-9*flops.float/secs
proc solve*(s:Wilson; r,x:Field; m:SomeNumber; res:float) =
  var sp = initSolverParams()
  sp.r2req = res
  #sp.maxits = 1000
  sp.verbosity = 1
  solve(s, r, x, m, sp)

proc solve2*(s:Wilson; r,x:Field; m:SomeNumber; res2:float) =
  var sp:SolverParams
  sp.r2req = res2
  sp.maxits = 100
  sp.verbosity = 1
  sp.subset.layoutSubset(r.l, "all")
  var t = newOneOf(r)
  proc op(a,b:Field) =
    #wilsonD2ee(s.se, s.so, r, s.g, x, m*m)
    threadBarrier()
    s.Ddag(t, b, m)
    threadBarrier()
    #echo t.norm2
    s.D(a, t, m)
    #a := b
    threadBarrier()
  #echo r.norm2
  var oa = (apply: op)
  cgSolve(r, x, oa, sp)
  threads:
    s.Ddag(t, r, m)
    threadBarrier()
    r := t

static:
  echo "Finished wilsonD: ", staticExec("date")

when isMainModule:
  import rng
  template foldl(f,n,op:untyped):untyped =
    var r:type(f(0))
    r = f(0)
    for i in 1..<n:
      let
        a {.inject.} = r
        b {.inject.} = f(i)
      r = op
    r

  proc runtest(v1,v2,sdAll,sdEven,sdOdd,s,m:any) =
    let g = s.g
    let lo = g[0].l
    #const nv = nVecs(v1[0])
    const nv = 1
    threads:
      v1 := 0
      v2 := 1
      threadBarrier()
      if myRank==0 and threadNum==0:
        when compiles(v1[0].len):
          v1{0}[0][0] := 1
        else:
          v1{0} := 1
      threadBarrier()
      echo v1.norm2

      wilsonD(sdAll, v2, g, v1, m)
      threadBarrier()
      echo v2.norm2
      #echo v2
      s.D(v2, v1, m)
      threadBarrier()
      echo v2.norm2
      echo "expect 20.81"

      for e in v1:
        template x(d:int):untyped = lo.vcoords(d,e)
        when compiles(v1[e].len):
          v1[e][0][0].re := foldl(x, 4, a*10.int16+b)
        else:
          for i in 0..<v1[e].ncols:
            v1[e][0,i].re := foldl(x, 4, a*10.int16+b)
        #echo v1[e][0]
      threadBarrier()
      wilsonD(sdAll, v2, g, v1, 0.5)
      echo v1[0][0]
      echo v2[0][0]

    template makeBench(name:untyped; bar:untyped):untyped {.dirty.} =
      proc `name T`(sd,v1,v2:any, ss="all") =
        var nrep = 1
        var dt = 0.0
        while true:
          resetTimers()
          threads:
            threadBarrier()
            let t0 = getTics()
            for rep in 1..nrep:
              wilsonD(sd, v2, g, v1, 0.5)
              when bar: threadBarrier()
            let t1 = getTics()
            var dtt = ticDiffSecs(t1,t0)
            threadSum(dtt)
            threadMaster: dt = dtt/numThreads.float
          if dt>1: break
          let nnrep = 1 + int(1.1*nrep.float/(dt+1e-9))
          nrep = min(10*nrep, nnrep)

        #var vol = lo.physVol.float
        var vol = lo.nSites.float
        if sd.sub != "all": vol *= 0.5
        let flops = nv * (24.0+g.len*2.0*(12+2*66+24)) * vol
        echo ss & "secs: ", dt, "  mf: ", (nrep.float*flops)/(1e6*dt)
        echoTimers()
      template name(sd:any, ss="all") = `name T`(sd, v1, v2, ss)
    subst(bench,_,benchB,_):
      makeBench(bench, false)
      makeBench(benchB, true)
      bench(sdAll, "all  ")
      benchB(sdAll, "all  ")
      bench(sdEven, "even ")
      benchB(sdEven, "even ")
      bench(sdOdd, "odd  ")
      benchB(sdOdd, "odd  ")
    proc benchEO() =
      var nrep = 1
      var dt = 0.0
      while true:
        resetTimers()
        threads:
          threadBarrier()
          let t0 = getTics()
          for rep in 1..nrep:
            wilsonD2ee(sdEven, sdOdd, v2, g, v1, 0.1)
            #s.D(v2, v1, 0.1)
            #when bar: threadBarrier()
          let t1 = getTics()
          var dtt = ticDiffSecs(t1,t0)
          threadSum(dtt)
          threadMaster: dt = dtt/numThreads.float
        if dt>1: break
        let nnrep = 1 + int(1.1*nrep.float/(dt+1e-9))
        nrep = min(10*nrep, nnrep)

      var vol = 0.5 * lo.nSites.float
      let flops = nv * (24.0+g.len*4.0*(12+2*66+24)) * vol
      echo "EO   secs: ", dt, "  mf: ", (nrep.float*flops)/(1e6*dt)
      echoTimers()
    benchEO()
    resetTimers()
    s.solve(v2, v1, 0.1, 1e-12)
    echoTimers()

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
    var i0 = 0
    if cp[0][0] notin {'0'..'9'}: inc i0
    for i in 0..<lat.len:
      lat[i] = (if (i0+i)<cp.len: parseInt(cp[i0+i]) else: lat[i-1])
  var lo = newLayout(lat)
  var v1 = lo.DiracFermion()
  var v2 = lo.DiracFermion()
  var rs = newRNGField(RngMilc6, lo, intParam("seed", 987654321).uint64)
  var g: array[4,type(lo.ColorMatrix())]
  for i in 0..<4:
    g[i] = lo.ColorMatrix()
    threads:
      g[i] := 1
  threads:
    g.setBC
    threadBarrier()
    for i in 0..<4:
      echo g[i].norm2

  #g.loadGauge("l88.scidac")
  var sdAll = initWilsonD(v1, "all")
  var sdEven = initWilsonD(v1, "even")
  var sdOdd = initWilsonD(v1, "odd")
  var s = newWilson(@g, v1)
  var m = 0.1
  echo "done newWilson"

  wilsonD(sdAll, v2, g, v1, m)

  runtest(v1, v2, sdAll, sdEven, sdOdd, s, m)

  #[
  var sdAll3 = initWilsonD3(v1, "all")
  var sdEven3 = initWilsonD3(v1, "even")
  var sdOdd3 = initWilsonD3(v1, "odd")
  var g3:array[8,type(lo.ColorMatrix())]
  for i in 0..3:
    g3[2*i  ] = g[i]
    #g3[2*i+1] = g[i]
    g3[2*i+1] = lo.ColorMatrix()
    g3[2*i+1].randomU rs
  var s3 = newWilson3(@g3)

  runtest(v1, v2, sdAll3, sdEven3, sdOdd3, s3, m)

  const nc = v1[0].len
  const nr = 8
  type MX* = Field[VLEN,MatrixArray[nr,nc,SComplexV]]
  #type MX* = Field[VLEN,MatrixArray[nr,nc,DComplexV]]
  var m1,m2: MX
  m1.new(lo)
  m2.new(lo)
  var sdAllM = initWilsonD(m1, "all")
  var sdEvenM = initWilsonD(m1, "even")
  var sdOddM = initWilsonD(m1, "odd")
  var sM = newWilson(@g,m1)
  echo "testing multi matrix: ", nr
  wilsonD(sdAllM, m2, g, m1, m)
  runtest(m1, m2, sdAllM, sdEvenM, sdOddM, sM, m)
  echoTimers()
  ]#

  #[
  #const n = 4
  var n = 4
  if cp.len>0:
    for i in 0..<cp.len:
      if cp[i][0]=='n':
        n = parseInt(cp[i][1..^1])
        break
  echo "n: ", n
  var v1a = newSeq[type(v1)](n)
  var v2a = newSeq[type(v2)](n)
  var sda = newSeq[type(sdAll)](n)
  var sda3 = newSeq[type(sdAll3)](n)
  #var sa = array[n,type(s)]
  v1a[0] = v1
  v2a[0] = v2
  sda[0] = sdAll
  sda3[0] = sdAll3
  #sda[0] = sdEven
  #sa[0] = s
  for i in 1..<n:
    v1a[i] = lo.ColorVector()
    v1a[i] := 1
    v2a[i] = lo.ColorVector()
    sda[i] = initWilsonD(v1, "all")
    sda3[i] = initWilsonD3(v1, "all")
    #sa[i] = newWilson(@g)

  let nrep = int(2e7/lo.physVol.float)
  template makeBenchN(name:untyped; bar:bool):untyped =
    proc name(sd,g:any, ss="all") =
      resetTimers()
      var t0 = epochTime()
      threads:
        for rep in 1..nrep:
          wilsonDN(sd, v2a, g, v1a, 0.5)
          when bar: threadBarrier()
      var t1 = epochTime()
      let dt = t1-t0
      #var vol = lo.physVol.float
      var vol = lo.nSites.float
      if sd[0].sub != "all": vol *= 0.5
      let flops = n*(6.0+g.len*2.0*72.0) * vol
      echo ss & "secs: ", dt, "  mf: ", (nrep.float*flops)/(1e6*dt)
      #echoTimers()

  makeBenchN(benchN, false)

  wilsonDN(sda, v2a, g, v1a, m)
  benchN(sda, g)
  #echoTimers()

  #wilsonDN(sda3, v2a, g3, v1a, m)
  #benchN(sda3, g3)
  ]#

  qexFinalize()
