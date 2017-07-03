import strutils
import qex, physics/stagD
import eigens/qexPrimmeInternal
import chebyshev

qexInit()
threads: echo "thread $# / $#"%[$threadNum, $numThreads]

template ff(x:untyped):auto = formatFloat(x,ffScientific,17)
template apply(op:Staggered, x,y:Field) =
  threads: op.so.stagD(x, op.g, y, 0.0)
template applyAdj(op:Staggered, x,y:Field) =
  threads: op.se.stagD(x, op.g, y, 0.0, -1.0)
template newVector(op:Staggered): untyped =
  op.g[0].l.ColorVector()

type
  Cheby = object
    a,b,s:float
    n:int
    count:int
var cheby = Cheby(a:1e-4, b:8.0, s:2.0, n:128, count:0)
cheby.CLIset a, "Cheby"
cheby.CLIset b, "Cheby"
cheby.CLIset s, "Cheby"
cheby.CLIset n, "Cheby"
type
  ChebyOp[O,T] = object
    o:O
    t:T
template apply(op:ChebyOp, y,x:Field) =
  inc cheby.count
  let
    i = 2.0/(cheby.b-cheby.a)
    c = (cheby.b+cheby.a)/2.0
  op.o.apply(op.t, x)
  op.o.applyAdj(y, op.t)
  threads: y := i*(y - c*x)
template assign(y:Field, x:FieldAddSub) =
  threads: assignM(y,x,x.type)
template assign(y:Field, x:Field) =
  threads: assignM(y,x,x.type)

# Augment Primme Operators
type
  MyOp[S,T] = object
    o:S
    c:ChebyOp[S,T]
# Requires knowledge of the internal applyD in qexPrimme.nim
# Note that there are pp.x and pp.y used in applyD,
# but only pp.x is used both when pass in and pass out.
# Also note that the apply(Adj) is called inside threads.
template apply(op:MyOp, x,y:Field) =
  # The above comment means we only need to modify the y for applyD.
  x.chebyshevT(cheby.n, op.c, y)
  threads: x += cheby.s*y
template applyAdj(op:MyOp, x,y:Field) =
  # and do a simple copy here.
  threads: x := y
template newVector(op:MyOp): untyped = op.o.newVector

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
  ae = floatParam("abserr", 1e-10)
  re = floatParam("relerr", 1e-10)
var
  s = newStag3(fl, ll)
  c = ChebyOp[type(s),type(s.newVector)](o:s,t:s.newVector)
  myO = MyOp[c.O,c.T](o:s,c:c)
  pp = lo.primmeInitialize(
    myO, abserr=ae, relerr=re, nVals = intParam("nv", 16),
    preset = intParam("method", 2).primme_preset_method)
pp.p.target = primme_largest
pp.p.CLIset maxBlockSize
pp.p.CLIset maxBasisSize
pp.p.CLIset minRestartSize
pp.p.CLIset printLevel
pp.p.restartingParams.CLIset maxPrevRetain
pp.p.CLIset eps, "":
  echo "Ignoring abserr and relerr."
  pp.p.convTestFun = nil
pp.prepare
pp.run
echo "ChebyMatvec : ",cheby.count
for i in 0..<pp.p.initSize:
  echo "$#  $#  $#"%[$i, pp.vals[i].ff, pp.rnorms[i].ff]
var
  evs = newseq[type(s.newVector)](pp.p.initSize)
  val = newseq[float](pp.p.initSize)
  err = newseq[float](pp.p.initSize)
  u = s.newVector
  v = s.newVector
for i in 0..<evs.len: evs[i] = s.newVector
for i in 0..<evs.len:
  var vn,r2:float
  threads:
    evs[i] := 0
    threadBarrier()
    evs[i].fromPrimmeArray(pp.vecs[i*pp.p.nLocal].addr)
    vn = evs[i].even.norm2
  s.apply(u, evs[i])
  # echo vn," ",ff(u.odd.norm2/vn)
  s.applyAdj(v, u)
  val[i] = sqrt(v.even.norm2/vn)
  threads:
    u.even := v - val[i]*evs[i]
    r2 = u.even.norm2/vn
  err[i] = abs(val[i] - sqrt(abs(val[i]*val[i]-sqrt(r2))))
  echo "$#  eval: $#  r: $#  err: $#"%[$i, val[i].ff, r2.sqrt.ff, err[i].ff]
qexFinalize()
