import qex
#import stdUtils
import field
import physics/qcdTypes
import gauge/gaugeUtils
import physics/stagSolve
#import profile
import os
import rng

qexInit()
installHelpParam()
processHelpParam()
#var defaultLat = [4,4,4,4]
#var defaultLat = [8,8,8,8]
var defaultLat = @[8,8,8,8]
#var defaultLat = @[12,12,12,12]
defaultSetup()
var v1 = lo.ColorVector()
var v2 = lo.ColorVector()
var r = lo.ColorVector()
var rs = newRNGField(RngMilc6, lo, intParam("seed", 987654321).uint64)
var mass = floatParam("mass", 0.000001)
var warm = floatParam("warm", 0.5)
var maxtime = floatParam("maxtime", 2.0)
var comm = getDefaultComm()
#sp.maxits = intParam("maxits",int(1e9/lo.nSitesOuter.float))
echo "mass: ", mass
echo "warm: ", warm
echo "maxtime(seconds): ", maxtime
threads:
  #g.random rs
  g.warm warm, rs
  threadBarrier()
  g.setBC
  threadBarrier()
  g.stagPhase
  v1 := 0
  #for e in v1:
  #  template x(d:int):untyped = lo.vcoords(d,e)
  #  v1[e][0].re := foldl(x, 4, a*10+b)
  #  #echo v1[e][0]
echo "g.norm2: ", g.norm2
if myRank==0:
  v1{0}[0] := 1
  #v1{2*1024}[0] := 1
echo v1.norm2
#var gs = lo.newGaugeS
#for i in 0..<gs.len: gs[i] := g[i]
var s = newStag(g)
threads:
  v2 := 0
  echo v2.norm2
  threadBarrier()
  s.D(v2, v1, mass)
  threadBarrier()
  #echoAll v2
  echo v2.norm2
#echo v2
var sp = initSolverParams()
sp.maxits = 10
s.solve(v2, v1, mass, sp)
sp.maxits = comm.broadcast int(1.0+(0.1*sp.iterations*maxtime)/sp.seconds)
#echo "maxits: ", sp.maxits
sp.resetStats
s.solve(v2, v1, mass, sp)
sp.maxits = comm.broadcast int(1.0+(sp.iterations*maxtime)/sp.seconds)
#echo "maxits: ", sp.maxits
sp.resetStats
resetTimers()
s.solve(v2, v1, mass, sp)
threads:
  echo "v2: ", v2.norm2
  echo "v2.even: ", v2.even.norm2
  echo "v2.odd: ", v2.odd.norm2
  s.D(r, v2, mass)
  threadBarrier()
  r := v1 - r
  threadBarrier()
  echo r.norm2
#echo v2
echoProf()

var g3:array[8,type(g[0])]
for i in 0..3:
  g3[2*i] = g[i]
  g3[2*i+1] = lo.ColorMatrix()
  g3[2*i+1].warmSU warm, rs
  g3[2*i+1] *= 0.1
for i in 0..<g3.len:
  echo &"g3[{i}]2: {g3[i].norm2}"
var s3 = newStag3(g3)
#s3.D(v2, v1, m)
sp.resetStats
sp.maxits = 10
s3.solve(v2, v1, mass, sp)
sp.maxits = comm.broadcast int(1.0+(0.1*sp.iterations*maxtime)/sp.seconds)
#echo "maxits: ", sp.maxits
sp.resetStats
s3.solve(v2, v1, mass, sp)
sp.maxits = comm.broadcast int(1.0+(sp.iterations*maxtime)/sp.seconds)
#echo "maxits: ", sp.maxits
sp.resetStats
resetTimers()
s3.solve(v2, v1, mass, sp)
echoProf()

qexFinalize()
finalizeParams()
