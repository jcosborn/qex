import strutils
import qex
import qex/eigens/qexPrimmeInternal

template ff(x:untyped):auto = formatFloat(x,ffScientific,17)

template apply(op:Staggered, x,y:Field) =
  threadBarrier()
  op.so.stagD(x, op.g, y, 0.0)
template applyAdj(op:Staggered, x,y:Field) =
  threadBarrier()
  op.se.stagD(x, op.g, y, 0.0, -1.0)
template newVector(op:Staggered): untyped =
  op.g[0].l.ColorVector()
template set(p:typed, n:untyped, runifset:untyped) =
  let
    o = p.n
    s = astToStr n
  when compiles(strParam(s, p.n)):
    p.n = type(p.n)strParam(s, p.n)
  elif compiles(intParam(s, p.n)):
    p.n = type(p.n)intParam(s, p.n)
  elif compiles(floatParam(s, p.n)):
    p.n = type(p.n)floatParam(s, p.n)
  else:
    {.fatal:"Cannot set argument "&s&" of "&astToStr(p)&" for command line.".}
  if o != p.n:
    runifset
    echo "Customize $# : $# -> $#"%[s, $o, $p.n]
template set(p:typed, n:untyped) =
  p.set n:
    discard

qexInit()
threads: echo "thread $# / $#"%[$threadNum, $numThreads]
var (lo, g, _) = setupLattice([8,8,8,8])

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
  pp = lo.primmeInitialize(
    s, abserr=ae, relerr=re, nVals = intParam("nv", 16),
    preset = intParam("method", 2).primme_preset_method)
pp.p.set locking
pp.p.set maxBlockSize
pp.p.set maxBasisSize
pp.p.set minRestartSize
pp.p.restartingParams.set maxPrevRetain
pp.p.set eps:
  echo "Ignoring abserr and relerr."
  pp.p.convTestFun = nil
pp.prepare
pp.run
for i in 0..<pp.p.initSize:
  echo "$#  $#  $#"%[$i, pp.vals[i].ff, pp.rnorms[i].ff]
var evs = newseq[type(s.newVector)](pp.p.initSize)
for i in 0..<evs.len: evs[i] = s.newVector
threads:
  for i in 0..<evs.len:
    evs[i] := 0
    threadBarrier()
    evs[i].fromPrimmeArray(pp.vecs[i*pp.p.nLocal].addr)

if intParam("dosolve", 1) > 0:
  # Following copied from hisqev with modifications.
  var m = floatParam("mass", 0.01)
  var m2 = m*m
  var src = s.newVector
  var src2 = s.newVector
  var d1 = s.newVector
  var d2 = s.newVector
  var d3 = s.newVector
  var r = s.newVector
  var t = s.newVector
  var t2 = s.newVector
  src := 0
  if myRank==0:
    src{0}[0] := 1

  proc getResid(rr: any, dd,ss: any) =
    apply(s, rr.odd.field, dd.even.field)
    applyAdj(s, t.even.field, rr.odd.field)
    rr.even := ss - 4.0*(t + m2*dd)
  proc rsolve(dt: any, sc: any, m: float, sp: var SolverParams) =
    getResid(t2, dt, sc)
    let s1 = sc.even.norm2
    let s2 = t2.even.norm2
    let r2req = sp.r2req
    sp.r2req = r2req * (s1/s2)
    s.solveEO(t, t2, m, sp)
    dt += t
    sp.r2req = r2req

  var sp = initSolverParams()
  sp.maxits = 20_000
  sp.r2req = 1e-16

  d1 := 0
  s.solveEO(d1, src, m, sp)
  getResid(r, d1, src)
  echo "r1: ", r.even.norm2

  let ng = pp.p.numEvals.int
  let nv = pp.p.initSize.int
  #var nps = @[0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 20, 30]
  var nps = @[0, nv]
  var step = ng
  var stepMin = intParam("stepMin", 1)
  while step>=stepMin:
    for i in countup(step,nv,step):
      if not nps.contains(i): nps.add i
    if step<=stepMin: break
    let st = step div stepMin
    let f = factor(st)[^1]
    step = max(stepMin, (st*stepMin) div f)
  var np0 = nv+1
  for ip in 0..<nps.len:
    var np = nps[ip]
    if np<np0:
      src2 := src
      d2 := 0
      np0 = 0
    for i in np0..<np:
      let c = evs[i].even.dot(src2)
      src2 -= c*evs[i]
      let s = 0.25/(pp.vals[i] + m2)
      d2 += (s*c)*evs[i]
    #d2 := 0.9999*d1
    d3 := d2
    rsolve(d3, src, m, sp)
    getResid(r, d3, src)
    echo "np: ", np, " its: ", sp.finalIterations,
         " time: ", sp.seconds, " r2: ", r.even.norm2
    np0 = np
qexFinalize()
