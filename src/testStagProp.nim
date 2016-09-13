import qex
import stdUtils
import field
import qcdTypes
import gaugeUtils
import stagD
import profile
import os

template corner(l, i):untyped =
  (l.coords[0][i].int and 1) + ((l.coords[1][i].int and 1) shl 1) +
   ((l.coords[2][i].int and 1) shl 2)

template addCorner(l, s, i):untyped =
  ((s + l.coords[0][i].int) and 1) +
   ((((s shr 1) + l.coords[1][i].int) and 1) shl 1) +
   ((((s shr 2) + l.coords[2][i].int) and 1) shl 2)

proc mesons(v:any) =
  let l = v.l
  let nt = l.physGeom[3]
  var c = newSeq[array[8,float]](nt)
  when true:
  #when false:
    #var x:VectorArray[3,DComplex]
    for i in l.sites:
      let t = l.coords[3][i]
      let s = l.corner(i)
      c[t][s] += v{i}.norm2()
      #assign(x, v{i})
      #c[t][s] += x.norm2
  else:
    threads:
      #var x:VectorArray[3,DComplex]
      for i in 0..<l.nSites:
        let t = l.coords[3][i]
        let s = l.corner(i)
        let tpar = (8*t+s) mod numThreads
        if tpar==threadNum:
          c[t][s] += v{i}.norm2()
          #assign(x, v{i})
          #c[t][s] += x.norm2
  rankSum(c)
  for s in 0..<8:
    echo "corner: ", s
    for t in 0..<nt:
      let r = c[t][s]
      echo t, " ", r
      #echo t, " ", c[t][s]
  echo "sum:"
  for t in 0..<nt:
    var r = c[t][0]
    for s in 1..<8:
      r += c[t][s]
    echo t, " ", r

proc mesonsV(v:any) =
  let l = v.l
  let nt = l.physGeom[3]
  var cv = newAlignedMem[array[8,type(v[0].norm2())]](nt)
  var c = newSeq[array[8,float64]](nt)
  threads:
    for e in 0..<l.nSitesOuter:
      let i = l.nSitesInner * e
      let t = l.coords[3][i]
      let s = l.corner(i)
      let tpar = (8*t+s) mod numThreads
      if tpar==threadNum:
        cv[t][s] += v[e].norm2
  for tt in 0..<nt:
    let tt0 = tt - l.coords[3][0] + nt
    for ss in 0..<8:
      for i in 0..<l.nSitesInner:
        let t = (tt0 + l.coords[3][i]) mod nt
        let s = l.addCorner(ss, i)
        c[t][s] += cv[tt][ss][i]
  rankSum(c)
  for s in 0..<8:
    echo "corner: ", s
    for t in 0..<nt:
      let r = c[t][s]
      echo t, " ", r
      #echo t, " ", c[t][s]
  echo "sum:"
  for t in 0..<nt:
    var r = c[t][0]
    for s in 1..<8:
      r += c[t][s]
    echo t, " ", r

when isMainModule:
  qexInit()
  var defaultGaugeFile = "l88.scidac"
  #var defaultLat = [4,4,4,4]
  var defaultLat = [8,8,8,8]
  #var defaultLat = @[8,8,8,16]
  #var defaultLat = @[12,12,12,12]
  defaultSetup()
  var v1 = lo.ColorVector()
  var v2 = lo.ColorVector()
  var r = lo.ColorVector()
  threads:
    g.setBC
    g.stagPhase
    v1 := 0
    #for e in v1:
    #  template x(d:int):untyped = lo.vcoords(d,e)
    #  v1[e][0].re := foldl(x, 4, a*10+b)
    #  #echo v1[e][0]
  #echo v1.norm2
  if myRank==0:
    v1{0}[0] := 1
    #v1{2*1024}[0] := 1
  echo v1.norm2
  var s = newStag(g)
  var m = 0.001
  threads:
    v2 := 0
    echo v2.norm2
    threadBarrier()
    s.D(v2, v1, m)
    threadBarrier()
    #echoAll v2
    echo v2.norm2
  #echo v2
  s.solve(v2, v1, m, 1e-8)
  resetTimers()
  s.solve(v2, v1, m, 1e-8)
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
  mesons(v2)
  mesonsV(v2)
  qexFinalize()
