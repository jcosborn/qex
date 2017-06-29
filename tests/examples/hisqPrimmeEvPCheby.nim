import strutils
import qex, physics/stagD
import eigens/qexPrimmeInternal
import chebyshev

qexInit()
threads: echo "thread $# / $#"%[$threadNum, $numThreads]
let
  chebyA = floatParam("chebyA", 1e-4)
  chebyB = floatParam("chebyB", 10.0)
  chebyN = intParam("chebyN", 128)
var cheby = newChebyshev(chebyA..chebyB, chebyN, 1.0/x)

template ff(x:untyped):auto = formatFloat(x,ffScientific,17)
template apply(op:Staggered, x,y:Field) =
  threadBarrier()
  op.so.stagD(x, op.g, y, 0.0)
template applyAdj(op:Staggered, x,y:Field) =
  threadBarrier()
  op.se.stagD(x, op.g, y, 0.0, -1.0)
template newVector(op:Staggered): untyped =
  op.g[0].l.ColorVector()
template set(p:typed, n:untyped) =
  let
    o = p.n
    s = astToStr n
  p.n = intParam(s, p.n).cint
  if o != p.n:
    echo "Customize $# : $# -> $#"%[s, $o, $p.n]

type
  Op[O,T] = object
    o:O
    t:T
template apply(op:Op, y,x:Field) =
  threads:
    op.o.apply(op.t, x)
    op.o.applyAdj(y, op.t)
template assign(y:Field, x:FieldMul) = assignM(y,x,x.type)
template assign(y:Field, x:FieldAddSub) = assignM(y,x,x.type)
template assign(y:Field, x:Field) = assignM(y,x,x.type)
proc applyApproxDInv(pp:ptr Primme; xi,yo:ptr complex[float]) =
  # x in, y out
  let o = Op[type(pp.o),type(pp.x)](o:pp.o,t:newOneOf(pp.x))
  threads:
    pp.x.fromPrimmeArray xi
  cheby.apply(pp.y, o, pp.x)
  threads:
    pp.y.toPrimmeArray yo

proc precond[O](x:pointer, ldx:ptr PRIMME_INT,
                y:pointer, ldy:ptr PRIMME_INT,
                blocksize:ptr cint,
                primme:ptr primme_params, err:ptr cint) {.noconv.} =
  var
    x = asarray[complex[cdouble]] x
    dx = ldx[]
    y = asarray[complex[cdouble]] y
    dy = ldy[]
  for i in 0..<blocksize[]:
    let
      xp = x[i*dx].addr         # Input vector
      yp = y[i*dy].addr         # Output
    applyApproxDInv(cast[ptr O](primme.matrix), xp, yp)
  err[] = 0

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
  pp = lo.primmeInitialize(
    s, abserr=ae, relerr=re, nVals = intParam("nv", 16),
    preset = intParam("method", 2).primme_preset_method)
pp.p.set maxBlockSize
pp.p.set maxBasisSize
pp.p.set minRestartSize
pp.p.set printLevel
pp.p.restartingParams.set maxPrevRetain
pp.p.applyPreconditioner = precond[type(pp)]
pp.p.correctionParams.precondition = 1
pp.prepare
pp.run
for i in 0..<pp.p.initSize:
  echo "$#  $#  $#"%[$i, pp.vals[i].ff, pp.rnorms[i].ff]
qexFinalize()
