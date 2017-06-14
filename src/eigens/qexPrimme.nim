import primme
import qexPrimmeInternal
import base, layout

type
  Operator* = concept o
    o.newVector
    #o.apply(Field, Field)
    #o.applyAdj(Field, Field)
  Primme*[Op,F,P] = object
    o*:Op
    x*,y*:F
    abserr*,relerr*:float
    m*:float
    p*:P
    vecs*: alignedMem[complex[float]] # Wrap it in Field later.
    realWork*: alignedMem[char]
    intWork*: alignedMem[char]
    vals*: seq[float]
    rnorms*: seq[float]

proc applyD(pp:ptr Primme; xi,yo:ptr complex[float]) =
  # x in, y out
  threads:
    pp.x.fromPrimmeArray xi
    pp.o.apply(pp.y, pp.x)
    pp.o.applyAdj(pp.x, pp.y)
    threadBarrier()
    pp.x.toPrimmeArray yo

proc matvec[O](x:pointer, ldx:ptr PRIMME_INT,
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
    applyD(cast[ptr O](primme.matrix), xp, yp)
  err[] = 0

proc convTest[O](val:ptr cdouble; vec:pointer; rNorm:ptr cdouble; isconv:ptr cint;
                 primme:ptr primme_params; ierr:ptr cint) {.noconv.} =
  let
    r = rNorm[].float
    pp = cast[ptr O](primme.matrix)
    re = pp.relerr
    ae = pp.abserr
  if r < 2*ae*sqrt(val[]) or r < 2*re*val[]:
    isconv[] = 1
  else: isconv[] = 0
  ierr[] = 0

proc primmeInitialize*(lo: Layout, op: Operator,
                       relerr:float = 1e-4, abserr:float = 1e-6, m:float = 0.0,
                       nVals:int = 16,
                       printLevel:int = 3,
                       preset:primme_preset_method = PRIMME_DEFAULT_METHOD): auto =
  var pp:Primme[type(op), type(op.newVector), primme_params]
  pp.o = op
  pp.x = op.newVector
  pp.y = op.newVector
  pp.abserr = abserr
  pp.relerr = relerr
  pp.m = m
  echo "RelErr target: ",relerr
  echo "AbsErr target: ",abserr
  const nc = pp.x[0].len
  pp.p = primme_initialize()
  pp.p.n = nc*lo.physVol div 2
  pp.p.matrixMatvec = matvec[type(pp)]
  pp.p.numEvals = nVals.cint
  pp.p.target = primme_smallest
  pp.p.numProcs = nRanks.cint
  pp.p.procId = myRank.cint
  pp.p.nLocal = nc*lo.nEven
  pp.p.globalSumReal = sumReal[primme_params]
  pp.p.printLevel = printLevel.cint
  pp.p.convTestFun = convTest[type(pp)]
  block primmeSetMethod:
    let ret = pp.p.set_method preset
    if 0 != ret:
      echo "ERROR: set_method returned with nonzero exit status: ", ret
      qexAbort()
  pp
  
proc prepare*[Op,F](pp:var Primme[Op,F,primme_params]) =
  let ret = zprimme(nil,nil,nil,pp.p.addr)
  if 1 != ret:
    echo "Error: zprimme(nil) returned with exit status: ", ret
    qexAbort()
  pp.vecs = newAlignedMemU[complex[float]]int(pp.p.numEvals*pp.p.nLocal)
  pp.vals = newseq[float]pp.p.numEvals
  pp.rnorms = newseq[float]pp.p.numEvals
  pp.intWork = newAlignedMemU[char]pp.p.intWorkSize
  pp.realWork = newAlignedMemU[char]pp.p.realWorkSize

proc run*[Op,F](pp:var Primme[Op,F,primme_params]) =
  pp.p.intWork = cast[ptr cint](pp.intWork.data)
  pp.p.realWork = pp.realWork.data
  pp.p.matrix = pp.addr
  if myRank == 0: pp.p.display_params
  let ret = pp.p.run(pp.vals,
                     asarray[complex[float]](pp.vecs.data)[],
                     pp.rnorms)
  if ret != 0:
    echo "Error: primme returned with nonzero exit status: ", ret
    qexAbort()
  echo "Neigens    : ",pp.p.initSize
  echo "Iterations : ",pp.p.stats.numOuterIterations
  echo "Restarts   : ",pp.p.stats.numRestarts
  echo "Matvecs    : ",pp.p.stats.numMatvecs
  echo "Preconds   : ",pp.p.stats.numPreconds
  echo "GlobalSums : ",pp.p.stats.numGlobalSum
  echo "VGlobalSum : ",pp.p.stats.volumeGlobalSum
  echo "OrthoIProd : ",pp.p.stats.numOrthoInnerProds
  echo "ElapsedT   : ",pp.p.stats.elapsedTime
  echo "MatvecT    : ",pp.p.stats.timeMatvec
  echo "PrecondT   : ",pp.p.stats.timePrecond
  echo "OrthoT     : ",pp.p.stats.timeOrtho
  echo "GlobalSumT : ",pp.p.stats.timeGlobalSum
  echo "EstMinEv   : ",pp.p.stats.estimateMinEVal
  echo "EstMaxEv   : ",pp.p.stats.estimateMaxEVal
  echo "EstMaxSv   : ",pp.p.stats.estimateLargestSVal
  echo "EstResid   : ",pp.p.stats.estimateResidualError
  echo "MaxConvTol : ",pp.p.stats.maxConvTol
  if pp.p.locking != 0 and pp.p.intWork != nil and pp.p.intWork[] == 1:
    echo "\nA locking problem has occurred."
    echo "Some eigenpairs do not have a residual norm less than the tolerance."
    echo "However, the subspace of evecs is accurate to the required tolerance."
  case pp.p.dynamicMethodSwitch:
  of -1: echo "Recommended method for next run: DEFAULT_MIN_MATVECS"
  of -2: echo "Recommended method for next run: DEFAULT_MIN_TIME"
  of -3: echo "Recommended method for next run: DYNAMIC (close call)"
  else: discard

export primme

when isMainModule:
  import qex, gauge, rng, physics/stagD
  template apply(op:Staggered, x,y:Field) =
    threadBarrier()
    op.so.stagD(x, op.g, y, 0.0)
  template applyAdj(op:Staggered, x,y:Field) =
    threadBarrier()
    op.se.stagD(x, op.g, y, 0.0, -1.0)
  template newVector(op:Staggered): untyped =
    op.g[0].l.ColorVector()

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
  var pp = lo.primmeInitialize(s, relerr=1e-6, abserr=1e-8,
                               preset = intParam("method", 2).primme_preset_method)
  pp.p.maxBasisSize = 32
  pp.prepare
  pp.run
  for i in 0..<pp.p.initSize:
    echo "Eval[",i,"]: ",pp.vals[i].ff," rnorm: ",pp.rnorms[i].ff
  # Must avoid calling free, because we allocate memory ourselves.
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
  opts.nvecs = intParam("nvecs", 32)
  opts.rrbs = intParam("rrbs", opts.nvecs)
  opts.relerr = 1e-4
  opts.abserr = 1e-6
  #opts.relerr = 1e-6
  #opts.abserr = 1e-8
  opts.svdits = intParam("svdits", 500)
  opts.maxup = 10
  var ev = hisqev(op, opts)
  import unittest
  var CT = 1e-8
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
          check eigenResults[i] ~= sqrt(pp.vals[i])
          check eigenResults[i] ~= ev[i].sv
          if not(eigenResults[i] ~= sqrt(pp.vals[i])) or
             not(eigenResults[i] ~= ev[i].sv):
            echo "expect: ", $eigenResults[i]
            echo "primme: ", sqrt(pp.vals[i]).ff
            echo "hisqev: ", ev[i].sv.ff
  qexFinalize()
