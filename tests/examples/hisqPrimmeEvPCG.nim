import strutils
import qex, qex/physics/stagD
import qex/eigens/qexPrimmeInternal

qexInit()
threads: echo "thread $# / $#"%[$threadNum, $numThreads]

template ff(x:untyped):auto = formatFloat(x,ffScientific,17)
template apply(op:Staggered, x,y:Field) =
  threadBarrier()
  op.so.stagD(x, op.g, y, 0.0)
template applyAdj(op:Staggered, x,y:Field) =
  threadBarrier()
  op.se.stagD(x, op.g, y, 0.0, -1.0)
template newVector(op:Staggered): untyped =
  op.g[0].l.ColorVector()
template set(p:typed, n:untyped, prefix:string, runifset:untyped) =
  let
    o = p.n
    s = prefix & astToStr(n)
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
template set(p:typed, n:untyped, prefix = "") =
  p.set n, prefix:
    discard

proc setPrecondSolver:SolverParams =
  echo "Set SolverParams for the preconditioner"
  result.r2req = 1e-6
  result.maxits = 5000
  result.verbosity = 1
  result.set r2req, "PCG"
  result.set maxits, "PCG"
  result.set verbosity, "PCG"
var precondSP = setPrecondSolver()
proc applyApproxDInv(pp:ptr Primme; xi,yo:ptr complex[float]) =
  # x in, y out
  threads:
    pp.x.fromPrimmeArray xi
  precondSP.subset.layoutSubset(pp.y.l, "even")
  let s {.global.} = floatParam("PCGs", 1e-4)
  var t = newOneOf(pp.y)
  proc op(a,b:pp.F) =
    pp.o.apply(t, b)
    pp.o.applyAdj(a, t)
    threadBarrier()
    a.even += s*b.even
  cgSolve(pp.y, pp.x, op, precondSP)
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
pp.p.restartingParams.set maxPrevRetain
pp.p.applyPreconditioner = precond[type(pp)]
pp.p.correctionParams.precondition = 1
pp.prepare
pp.run
for i in 0..<pp.p.initSize:
  echo "$#  $#  $#"%[$i, pp.vals[i].ff, pp.rnorms[i].ff]
qexFinalize()
