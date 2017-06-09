import strutils
import primme
import base, layout, comms/qmp, physics/qcdTypes, physics/stagD

type
  OpInfo*[S,F] = object
    s:ptr S
    m:float
    x,y:F
proc newOpInfo*[S](s:ptr S, m:float = 0):auto =
  type F = type(s.g[0].l.ColorVectorD())
  var r:OpInfo[S,F]
  r.s = s
  r.m = m
  r.x.new(s.g[0].l)
  r.y.new(s.g[0].l)
  return r

# WARNING: low level implementation details follow.
proc toPrimmeArray(f:Field, a:ptr float, l:int) =
  const
    vl = f.V
    cl = 2*vl
  var
    a = asarray[float]a
    f = asarray[float]f.s.data
  let n = l div cl
  tfor i, 0..<n:
    forO j, 0, <vl:
      a[cl*i+2*j] = f[cl*i+j]
      a[cl*i+2*j+1] = f[cl*i+j+vl]
template toPrimmeArray(f:Field, a:ptr float) =
  toPrimmeArray(f, a, 6*f.l.nEven)
proc fromPrimmeArray(f:Field, a:ptr float, l:int) =
  const
    vl = f.V
    cl = 2*vl
  var
    a = asarray[float]a
    f = asarray[float]f.s.data
  let n = l div cl
  tfor i, 0..<n:
    forO j, 0, <vl:
      f[cl*i+j] = a[cl*i+2*j]
      f[cl*i+j+vl] = a[cl*i+2*j+1]
template fromPrimmeArray(f:Field, a:ptr float) =
  fromPrimmeArray(f, a, 6*f.l.nEven)

proc applyD(op:ptr OpInfo; x,y:ptr float) =
  # x in, y out
  threads:
    op.x.fromPrimmeArray x
    threadBarrier()
    stagD(op.s.so, op.y, op.s.g, op.x, 0.0)
    threadBarrier()
    stagD(op.s.se, op.x, op.s.g, op.y, 0.0, -1.0)
    threadBarrier()
    op.x.toPrimmeArray y

proc matvec[O](x:pointer, ldx:ptr PRIMME_INT,
               y:pointer, ldy:ptr PRIMME_INT,
               blocksize:ptr cint,
               primme:ptr primme_params, err:ptr cint) {.noconv.} =
  var
    x = asarray[cdouble] x
    dx = ldx[]
    y = asarray[cdouble] y
    dy = ldy[]
  for i in 0..<blocksize[]:
    let
      xp = x[i*dx].addr         # Input vector
      yp = y[i*dy].addr         # Output
    applyD(cast[ptr O](primme.matrix), xp, yp)
  err[] = 0

proc sumReal(sendBuf: pointer; recvBuf: pointer; count: ptr cint;
               primme: ptr primme_params; ierr: ptr cint) {.noconv.} =
  for i in 0..<count[]:
    asarray[float](recvBuf)[i] = asarray[float](sendBuf)[i]
  QMP_sum_double_array(cast[ptr cdouble](recvBuf), count[])
  ierr[] = 0

proc primmeInitialize*(lo: Layout, op: var OpInfo): primme_params =
  result = primme_initialize()
  result.n = 3*lo.physVol div 2
  result.matrixMatvec = matvec[type(op)]
  result.numEvals = 16
  result.target = primme_smallest
  result.eps = 1e-9
  result.numProcs = nRanks.cint
  result.procId = myRank.cint
  result.nLocal = 3*lo.nEven
  result.globalSumReal = sumReal
  result.matrix = op.addr
  result.printLevel = 3
  let ret = result.set_method PRIMME_DYNAMIC
  if 0 != ret:
    echo "ERROR: set_method returned with nonzero exit status: ", ret
    quit QuitFailure

template ff(x:untyped):auto = formatFloat(x,ffScientific,17)
type
  PrimmeEvs* = object
    evecs: alignedMem[complex[float]] # Wrap it in Field later.
    intWork: alignedMem[char]
    realWork: alignedMem[char]
    evals*: seq[float]
    rnorms*: seq[float]
proc run*(param: var primme_params): PrimmeEvs =
  result.evecs = newAlignedMemU[complex[float]]int(param.numEvals*param.nLocal)
  result.evals = newseq[float]param.numEvals
  result.rnorms = newseq[float]param.numEvals
  block primmeSetSize:
   let ret = zprimme(nil,nil,nil,param.addr)
   if 1 != ret:
     echo "Error: zprimme(nil) returned with exit status: ", ret
     quit QuitFailure
  result.intWork = newAlignedMemU[char]param.intWorkSize
  result.realWork = newAlignedMemU[char]param.realWorkSize
  param.intWork = cast[ptr cint](result.intWork.data)
  param.realWork = result.realWork.data
  if myRank == 0: param.display_params
  let ret = param.run(result.evals,
                      asarray[complex[float]](result.evecs.data)[],
                      result.rnorms)
  if ret != 0:
    echo "Error: primme returned with nonzero exit status: ", ret
    quit QuitFailure
  echo "Neigens    : ",param.initSize
  echo "Tolerance  : ",ff param.aNorm*param.eps
  echo "Iterations : ",param.stats.numOuterIterations
  echo "Restarts   : ",param.stats.numRestarts
  echo "Matvecs    : ",param.stats.numMatvecs
  echo "Preconds   : ",param.stats.numPreconds
  echo "GlobalSums : ",param.stats.numGlobalSum
  echo "VGlobalSum : ",param.stats.volumeGlobalSum
  echo "OrthoIProd : ",param.stats.numOrthoInnerProds
  echo "ElapsedT   : ",param.stats.elapsedTime
  echo "MatvecT    : ",param.stats.timeMatvec
  echo "PrecondT   : ",param.stats.timePrecond
  echo "OrthoT     : ",param.stats.timeOrtho
  echo "GlobalSumT : ",param.stats.timeGlobalSum
  echo "EstMinEv   : ",param.stats.estimateMinEVal
  echo "EstMaxEv   : ",param.stats.estimateMaxEVal
  echo "EstMaxSv   : ",param.stats.estimateLargestSVal
  echo "EstResid   : ",param.stats.estimateResidualError
  echo "MaxConvTol : ",param.stats.maxConvTol
  if param.locking != 0 and param.intWork != nil and param.intWork[] == 1:
    echo "\nA locking problem has occurred."
    echo "Some eigenpairs do not have a residual norm less than the tolerance."
    echo "However, the subspace of evecs is accurate to the required tolerance."
  case param.dynamicMethodSwitch:
  of -1: echo "Recommended method for next run: DEFAULT_MIN_MATVECS"
  of -2: echo "Recommended method for next run: DEFAULT_MIN_TIME"
  of -3: echo "Recommended method for next run: DYNAMIC (close call)"
  else: discard

export primme

when isMainModule:
  import qex, gauge, rng
  qexInit()
  var lat = [4,4,4,4]
  threads:
    echo "thread ", threadNum, "/", numThreads
  var
    lo = lat.newLayout
    g = newSeq[type(lo.ColorMatrix())](lat.len)
    r = RngMilc6.newRNGField(lo, 987654321)
  for mu in 0..<lat.len: g[mu] = lo.Colormatrix
  threads:
    g.random r
    g.setBC
    g.stagPhase
  var s = g.newStag
  var m = 0.1
  var opInfo = newOpInfo(s.addr, m)
  var pp = lo.primme_initialize(opInfo)
  var pevs = pp.run
  for i in 0..<pp.initSize:
    echo "Eval[",i,"]: ",pevs.evals[i].ff," rnorm: ",pevs.rnorms[i].ff
  #pp.free
  import hisqev
  type MyOp = object
    s: type(s)
    r: type(r)
    lo: type(lo)
  var op = MyOp(r:r,s:s,lo:lo)
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
  var opts: EigOpts
  opts.initOpts
  opts.nev = intParam("nev", 16)
  opts.nvecs = intParam("nvecs", (opts.nev*11) div 10)
  opts.rrbs = intParam("rrbs", opts.nvecs)
  opts.relerr = 1e-4
  opts.abserr = 1e-6
  #opts.relerr = 1e-6
  #opts.abserr = 1e-8
  opts.svdits = intParam("svdits", 500)
  opts.maxup = 10
  var evals0 = hisqev(op, opts)
  import unittest
  var CT = 1e-10                  # comparison tolerance
  const eigenResults = [
    0.0055851198854,
    0.0113277827524,
    0.0190303521552,
    0.0246828813692,
    0.0278549950745,
    0.0336481575597,
    0.0395180223843,
    0.0448678361203,
    0.0535298681860,
    0.0599862410362,
    0.0651594072052,
    0.0697624694074,
    0.0726489452823,
    0.0772123224085,
    0.0920417369866,
    0.0951845085945]
  proc `~=`(x,y:float):bool = abs(x-y)/max(abs(x),abs(y)) < CT
  if myrank == 0:
    suite "primme vs. hisqev":
      test "First 16 evs":
        forStatic i, 0, 15:
          check eigenResults[i] ~= sqrt(pevs.evals[i])
          check eigenResults[i] ~= evals0[i].sv
          if not(eigenResults[i] ~= sqrt(pevs.evals[i])) or
             not(eigenResults[i] ~= evals0[i].sv):
            echo "expect: ", $eigenResults[i]
            echo "primme: ", sqrt(pevs.evals[i]).ff
            echo "hisqev: ", evals0[i].sv.ff
  qexFinalize()
