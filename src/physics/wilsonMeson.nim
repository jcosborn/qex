import ../base/globals
setForceInline(false)

import base
import gauge
import os

proc wilsonMesons*(v: any) =
  let l = v.l
  let nt = l.physGeom[3]
  var c = newSeq[array[16,float]](nt)
  when true:
  #when false:
    #var x:VectorArray[3,DComplex]
    for i in l.sites:
      let t = l.coords[3][i]
      c[t][0] += v{i}.norm2()
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
  for s in 0..0:
    echo "gamma: ", 15-s
    for t in 0..<nt:
      let r = c[t][s]
      echo t, " ", r
      #echo t, " ", c[t][s]
  echo "sum:"
  for t in 0..<nt:
    var r = c[t][0]
    #for s in 1..<8:
    #  r += c[t][s]
    echo t, " ", r

when isMainModule:
  import qex
  import physics/qcdTypes
  import physics/wilsonD
  qexInit()
  var defaultGaugeFile = "l88.scidac"
  #var defaultLat = [4,4,4,4]
  var defaultLat = [8,8,8,8]
  #var defaultLat = @[8,8,8,16]
  #var defaultLat = @[12,12,12,12]
  defaultSetup()
  var v1 = lo.DiracFermionS()
  var v2 = lo.DiracFermionS()
  var v3 = lo.DiracFermionS()
  var r = lo.DiracFermionS()
  threads:
    g.setBC
    v1 := 0
    #for e in v1:
    #  template x(d:int):untyped = lo.vcoords(d,e)
    #  v1[e][0].re := foldl(x, 4, a*10+b)
    #  #echo v1[e][0]
  #echo v1.norm2
  if myRank==0:
    v1{0}[0][0] := 1
    v1{0}[1][1] := 1
    #v1{2*1024}[0] := 1
  echo "v1sq: ", v1.norm2
  var gs = newGaugeS(g)
  var s = newWilsonS(gs)
  #var s = newWilson(g)
  var m = floatParam("mass", 0.001)
  echo "mass: ", m
  threads:
    #v2 := 0
    #v3 := 0
    #echo v2.norm2
    #threadBarrier()
    s.D(v2, v1, m)
    threadBarrier()
    echo "v1'Dv1: ", v1.dot(v2)
    s.D(v3, v2, m)
    threadBarrier()
    echo "v1'DDv1: ", v1.dot(v3)
    s.D(v2, v3, m)
    threadBarrier()
    echo "v1'DDDv1: ", v1.dot(v2)
    s.D(v3, v2, m)
    threadBarrier()
    echo "v1'DDDDv1: ", v1.dot(v3)
  let rsq = 1e-20
  #echo v2
  v2 := 0
  #s.solve2(v2, v1, m, 1e-12)
  s.solveEO(v2, v1, m, rsq)
  #s.solve(v2, v1, m, 1e-12)
  resetTimers()
  v2 := 0
  #s.solve2(v2, v1, m, 1e-12)
  s.solveEO(v2, v1, m, rsq)
  #s.solve(v2, v1, m, 1e-12)
  threads:
    echo "v2: ", v2.norm2
    echo "v2.even: ", v2.even.norm2
    echo "v2.odd: ", v2.odd.norm2
    s.D(r, v2, m)
    threadBarrier()
    r := v1 - r
    threadBarrier()
    echo "r2: ", r.norm2
  #echo v2
  #s.Ddag(v3, v2, m)
  wilsonMesons(v2)
  for i in 0..3:
    echo i, ": ", v2[0][i][0].re[0], " ", v2[0][i][0].im[0]

  echoTimers()
  qexFinalize()
