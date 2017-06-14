import os, strutils
import qex, physics/stagD

template ff(x:untyped):auto = formatFloat(x,ffScientific,17)
template apply(op:Staggered, x,y:Field) =
  threadBarrier()
  op.so.stagD(x, op.g, y, 0.0)
template applyAdj(op:Staggered, x,y:Field) =
  threadBarrier()
  op.se.stagD(x, op.g, y, 0.0, -1.0)
template newVector(op:Staggered): untyped =
  op.g[0].l.ColorVector()

qexInit()
threads: echo "thread $# / $#"%[$threadNum, $numThreads]
var (lo, g, r) = setupLattice([8,8,8,8])

threads:
  g.setBC
  g.stagPhase
var hc: HisqCoefs
hc.init
echo hc
var
  fl = lo.newGauge
  ll = lo.newGauge
hc.smear(g, fl, ll)
let
  ae = floatParam("abserr", 1e-6)
  re = floatParam("relerr", 1e-4)
var
  s = newStag3(fl, ll)
  pp = lo.primmeSVDInitialize(
    s, abserr=ae, relerr=re, nVals = intParam("nv", 16),
    presetStage1 = intParam("methodStage1", 2).primme_preset_method)
pp.prepare
pp.run
for i in 0..<pp.p.initSize:
  echo "$#  $#  $#"%[$i, pp.vals[i].ff, pp.rnorms[i].ff]
qexFinalize()
