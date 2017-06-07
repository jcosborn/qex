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
  pp = lo.primmeInitialize opInfo
let Nev = intParam("nev", 16)
pp.numEvals = Nev.cint
pp.eps = floatParam("eps", 1e-12)
let pevs = pp.run
for i in 0..<pp.initSize:
  echo "$#  $#  $#"%[$i, pevs.evals[i].ff, pevs.rnorms[i].ff]
pp.free

type MyOp = object
  s: type(s)
  r: type(r)
  lo: type(lo)
template rand(op: var MyOp, v: any) =
  gaussian(v, op.r)
template newVector(op: MyOp): untyped =
  op.lo.ColorVector()
template apply(op: MyOp, r,v: typed) =
  threadBarrier()
  stagD(op.s.so, r.field, op.s.g, v.field, 0.0)
template applyAdj(op: MyOp, r,v: typed) =
  threadBarrier()
  stagD(op.s.se, r.field, op.s.g, v.field, 0.0, -1)
template newRightVec(op: MyOp): untyped = newVector(op).even
template newLeftVec(op: MyOp): untyped = newVector(op).odd
if 0 < intParam("hisqev", 0):
  var op = MyOp(r:r,s:s,lo:lo)
  var opts: EigOpts
  opts.initOpts
  opts.nev = Nev
  opts.nvecs = intParam("nvecs", (opts.nev*11) div 10)
  opts.rrbs = intParam("rrbs", opts.nvecs)
  opts.relerr = 1e-4
  opts.abserr = 1e-6
  #opts.relerr = 1e-6
  #opts.abserr = 1e-8
  opts.svdits = intParam("svdits", 5000)
  opts.maxup = 10
  var evs = hisqev(op, opts)
  echo "HISQEV evs:"
  for i in 0..<evs.len:
    echo "$#  $#  $#"%[$i, evs[i].sv.ff, evs[i].err.ff]
qexFinalize()
