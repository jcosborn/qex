import os, strutils
import qex
template ff(x:untyped):auto = formatFloat(x,ffScientific,17)
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
var
  s = newStag3(fl, ll)
  opInfo = newOpInfo(s.addr)
  pp = lo.primmeSVDInitialize opInfo
pp.numSvals = intParam("nsv", 16).cint
pp.eps = floatParam("eps", 1e-10)
let pevs = pp.run
for i in 0..<pp.initSize:
  echo "$#  $#  $#"%[$i, pevs.vals[i].ff, pevs.rnorms[i].ff]
qexFinalize()
