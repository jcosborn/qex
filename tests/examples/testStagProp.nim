import qex
import physics/qcdTypes
import gauge
import physics/stagSolve
import os
import physics/stagMesonLocal
import physics/hisqLinks

qexInit()
#var defaultGaugeFile = "l88.scidac"
#var defaultLat = [4,4,4,4]
#var defaultLat = [8,8,8,8]
#var defaultLat = [8,8,8,16]
#var defaultLat = @[12,12,12,12]
var defaultLat = latticeFromLocalLattice([8,8,8,8], nRanks)
defaultSetup()
var v1 = lo.ColorVector()
var v2 = lo.ColorVector()
var r = lo.ColorVector()
var hc: HisqCoefs
hc.init()
var fl = lo.newGauge()
var ll = lo.newGauge()
threads:
  g.setBC
  g.stagPhase
  v1 := 0
  #for e in v1:
  #  template x(d:int):untyped = lo.vcoords(d,e)
  #  v1[e][0].re := foldl(x, 4, a*10+b)
  #  #echo v1[e][0]
#echo v1.norm2
hc.smear(g, fl, ll)
if myRank==0:
  discard v1{0}
  discard v1{0}[0]
  v1{0}[0] := 1
  #v1{2*1024}[0] := 1
echo v1.norm2
var s = newStag3(fl, ll)
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
stagMesons(v2)
stagMesonsV(v2)
qexFinalize()
