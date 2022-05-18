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

proc convTest[O](val:ptr cdouble; vec:pointer; rNorm:ptr cdouble; isconv:ptr cint;
                 primme:ptr primme_params; ierr:ptr cint) {.noconv.} =
  let
    pp = cast[ptr O](primme.matrix)
    re = pp.relerr
    ae = pp.abserr
  var
    v = val[].float
    r = rNorm[].float
  if vec != nil:
    echo "convTest: filtered v: ",v,"    r: ",r
    var u = pp.o.newVector
    threads: pp.x.fromPrimmeArray cast[ptr complex[float]](vec)
    pp.o.o.apply(u, pp.x)
    pp.o.o.applyAdj(pp.y, u)
    threads:
      v = sqrt(pp.y.even.norm2/pp.x.even.norm2)
      pp.y.even -= v*pp.x
      r = sqrt(pp.y.even.norm2)
    echo "convTest: original v: ",v,"    r: ",r
  if r < 2*ae*sqrt(v) or r < 2*re*v:
    isconv[] = 1
  else: isconv[] = 0
  ierr[] = 0

type
  PChebyOp[O,T] = object
    o:O
    t:T
    s:float
template apply(op:PChebyOp, y,x:Field) =
  op.o.apply(op.t, x)
  op.o.applyAdj(y, op.t)
  threads: y += op.s*x
template assign(y:Field, x:FieldMul) = assignM(y,x,x.type)
template assign(y:Field, x:FieldAddSub) = assignM(y,x,x.type)
template assign(y:Field, x:Field) = assignM(y,x,x.type)

proc cg(x:Field; b:Field2; A:proc; sp:var SolverParams) =
  # A copy of cgSolve with A running outside of the thread block.
  tic()
  let vrb = sp.verbosity
  template verb(n:int; body:untyped):untyped =
    if vrb>=n: body
  let sub = sp.subset
  template subset(body:untyped):untyped =
    onNoSync(sub):
      body
  template mythreads(body:untyped):untyped =
    threads:
      onNoSync(sub):
        body

  var b2: float
  mythreads:
    x := 0
    b2 = b.norm2
  verb(1):
    echo("input norm2: ", b2)
  if b2 == 0.0:
    sp.finalIterations = 0
    return

  var r = newOneOf(x)
  var p = newOneOf(x)
  var Ap = newOneOf(x)
  let r2stop = sp.r2req * b2;
  let maxits = sp.maxits
  var finalIterations = 0

  threads:
    subset:
      p := 0
      r := b
      verb(3):
        echo("p2: ", p.norm2)
        echo("r2: ", r.norm2)

  var itn = 0
  var r2 = b2
  var r2o = r2
  verb(1):
    #echo(-1, " ", r2)
    echo(itn, " ", r2/b2)
  toc("cg setup")

  while itn<maxits and r2>=r2stop:
    tic()
    inc itn
    let beta = r2/r2o;
    r2o = r2
    threads:
      subset:
        p := r + beta*p
    toc("p update", flops=2*numNumbers(r[0])*sub.lenOuter)
    A(Ap, p)
    toc("Ap")
    threads:
      subset:
        let pAp = p.redot(Ap)
        toc("pAp", flops=2*numNumbers(p[0])*sub.lenOuter)
        let alpha = r2/pAp
        x += alpha*p
        toc("x", flops=2*numNumbers(p[0])*sub.lenOuter)
        r -= alpha*Ap
        toc("r", flops=2*numNumbers(r[0])*sub.lenOuter)
        r2 = r.norm2
        toc("r2", flops=2*numNumbers(r[0])*sub.lenOuter)
    verb(2):
      #echo(itn, " ", r2)
      echo(itn, " ", r2/b2)
    verb(3):
      threads:
        subset:
          let pAp = p.redot(Ap)
          echo "beta: ", beta
          echo "p2: ", p.norm2
          echo "Ap2: ", Ap.norm2
          echo "pAp: ", pAp
          echo "alpha: ", r2o/pAp
          echo "x2: ", x.norm2
          echo "r2: ", r2
      A(Ap, x)
      var fr2: float
      threads:
        subset:
          fr2 = (b - Ap).norm2
      echo "   ", fr2/b2
  toc("cg iterations")
  if threadNum==0: finalIterations = itn

  var fr2: float
  A(Ap, x)
  threads:
    subset:
      r := b - Ap
      fr2 = r.norm2
  verb(1):
    echo finalIterations, " acc r2:", r2/b2
    echo finalIterations, " tru r2:", fr2/b2

  sp.finalIterations = finalIterations
  toc("cg final")

type
  PKind = enum
    pNone,
    pCheby,
    pCG
  Precond = object
    case kind:PKind
    of pNone: discard
    of pCheby:
      cheby:Chebyshev
      chebys:float
    of pCG:
      cg:SolverParams
      cgs:float
let precondKind = intParam("PKind", 0)
var precond:Precond
if 0 == precondKind:
  echo "Preconditioner: None"
  precond.kind = pNone
elif 1 == precondKind:
  echo "Preconditioner: Chebyshev"
  precond.kind = pCheby
  precond.CLIset chebys, "Pshift"
  let
    chebyA = floatParam("PChebya", precond.chebys+cheby.s-1)
    chebyB = floatParam("PChebyb", precond.chebys+cheby.s+2)
    chebyN = intParam("PChebyn", 6)
  echo "PChebya: ",chebyA
  echo "PChebyb: ",chebyB
  echo "PChebyn: ",chebyN
  precond.cheby = newChebyshev(chebyA..chebyB, chebyN, 1.0/x)
  echo "Pshift: ",precond.chebys
elif 2 == precondKind:
  echo "Preconditioner: CG"
  precond.kind = pCG
  precond.CLIset cgs, "Pshift"
  precond.cg.r2req = 1e-8
  precond.cg.maxits = 1000
  precond.cg.verbosity = 1
  precond.cg.CLIset r2req, "PCG"
  precond.cg.CLIset maxits, "PCG"
  precond.cg.CLIset verbosity, "PCG"
  echo "Pshift: ",precond.cgs
else:
  echo "ERROR: PKind can only be 0..2"
  qexAbort()

proc applyApproxDInv(pp:ptr Primme; xi,yo:ptr complex[float]) =
  # x in, y out
  threads:
    pp.x.fromPrimmeArray xi
  case precond.kind:
  of pNone:
    echo "Internal ERROR: applyApproxDInv shouldn't be called with PKind of pNone."
    qexAbort()
  of pCheby:
    let o = PChebyOp[type(pp.o),type(pp.x)](o:pp.o,t:newOneOf(pp.x),s:precond.chebys)
    precond.cheby.apply(pp.y, o, pp.x)
  of pCG:
    precond.cg.subset.layoutSubset(pp.y.l, "even")
    var t = newOneOf(pp.y)
    proc op(a,b:pp.F) =
      pp.o.apply(t, b)
      pp.o.applyAdj(a, t)
      threads: a.even += precond.cgs*b.even
    cg(pp.y, pp.x, op, precond.cg)
  threads:
    pp.y.toPrimmeArray yo

proc precondFun[O](x:pointer, ldx:ptr PRIMME_INT,
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
  c = ChebyOp[type(s),type(s.newVector)](o:s,t:s.newVector)
  myO = MyOp[c.O,c.T](o:s,c:c)
  pp = lo.primmeInitialize(
    myO, abserr=ae, relerr=re, nVals = intParam("nv", 16),
    preset = intParam("method", 2).primme_preset_method)
pp.p.convTestFun = convTest[type(pp)]
pp.p.target = primme_largest
pp.p.CLIset maxBlockSize
pp.p.CLIset maxBasisSize
pp.p.CLIset minRestartSize
pp.p.CLIset printLevel
pp.p.restartingParams.CLIset maxPrevRetain
pp.p.CLIset eps, "":
  echo "Ignoring abserr and relerr."
  pp.p.convTestFun = nil
if precondKind > 0:
  pp.p.applyPreconditioner = precondFun[type(pp)]
  pp.p.correctionParams.precondition = 1
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
