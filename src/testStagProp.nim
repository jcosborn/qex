import qex
import stdUtils
import field
import qcdTypes
import gaugeUtils
import stagD
import profile
import os

template corner(l, i):expr =
  (l.coords[0][i].int and 1) + ((l.coords[1][i].int and 1) shl 1) +
   ((l.coords[2][i].int and 1) shl 2)

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

when isMainModule:
  qexInit()
  var defaultGaugeFile = "l88.scidac"
  #var defaultLat = [4,4,4,4]
  #var defaultLat = [8,8,8,8]
  var defaultLat = @[8,8,8,16]
  defaultSetup()
  var v1 = lo.ColorVector()
  var v2 = lo.ColorVector()
  var r = lo.ColorVector()
  threads:
    g.setBC
    g.stagPhase
    v1 := 0
    #for e in v1:
    #  template x(d:int):expr = lo.vcoords(d,e)
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
  qexFinalize()
