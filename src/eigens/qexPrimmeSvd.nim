import primme
import base, layout, physics/stagD
import qexPrimmeInternal
import qexPrimme

proc applyD(op:ptr OpInfo; x,y:ptr float) =
  # x in, y out
  threads:
    op.x.fromPrimmeArray x      # even
    threadBarrier()
    stagD(op.s.so, op.y, op.s.g, op.x, 0.0)
    threadBarrier()
    op.y.toPrimmeArray(y, op.y.l.nEven) # odd
proc applyDdag(op:ptr OpInfo; x,y:ptr float) =
  # x in, y out
  threads:
    op.x.fromPrimmeArray(x, op.x.l.nEven) # odd
    threadBarrier()
    stagD(op.s.se, op.y, op.s.g, op.x, 0.0, -1.0)
    threadBarrier()
    op.y.toPrimmeArray y        # even
proc matvec[O](x:pointer, ldx:ptr PRIMME_INT,
               y:pointer, ldy:ptr PRIMME_INT,
               blocksize:ptr cint, transpose:ptr cint,
               primme:ptr primme_svds_params, err:ptr cint) {.noconv.} =
  var
    x = asarray[cdouble] x
    dx = ldx[]
    y = asarray[cdouble] y
    dy = ldy[]
    op = cast[ptr O](primme.matrix)
  if 0 == transpose[]:          # Do y <- A x
    for i in 0..<blocksize[]:
      let
        xp = x[i*dx].addr       # Input vector
        yp = y[i*dy].addr       # Output
      op.applyD(xp, yp)
  else:                         # Do y <- A^dag x
    for i in 0..<blocksize[]:
      let
        xp = x[i*dx].addr       # Input vector
        yp = y[i*dy].addr       # Output
      op.applyDdag(xp, yp)
  err[] = 0

proc primmeSVDInitialize*(lo: Layout, op: var OpInfo): primme_svds_params =
  const nc = op.nc
  result = primme_svds_initialize()
  result.n = nc*lo.physVol div 2
  result.m = result.n
  result.matrixMatvec = matvec[type(op)]
  result.numSvals = 16
  result.target = primme_svds_smallest
  result.eps = 1e-9
  result.numProcs = nRanks.cint
  result.procId = myRank.cint
  result.nLocal = nc*lo.nEven
  result.mLocal = result.nLocal
  result.globalSumReal = sumReal[primme_svds_params]
  result.matrix = op.addr
  result.printLevel = 3
  let ret = result.set_method(primme_svds_default,
                              PRIMME_DEFAULT_METHOD, PRIMME_DEFAULT_METHOD)
  if 0 != ret:
    echo "ERROR: set_method returned with nonzero exit status: ", ret
    quit QuitFailure

proc run*(param: var primme_svds_params): PrimmeResults =
  result.vecs = newAlignedMemU[complex[float]]int(param.numSvals*(param.nLocal+param.mLocal))
  result.vals = newseq[float]param.numSvals
  result.rnorms = newseq[float]param.numSvals
  block primmeSetSize:
   let ret = zprimme_svds(nil,nil,nil,param.addr)
   if 1 != ret:
     echo "Error: zprimme_svds(nil) returned with exit status: ", ret
     quit QuitFailure
  result.intWork = newAlignedMemU[char]param.intWorkSize
  result.realWork = newAlignedMemU[char]param.realWorkSize
  param.intWork = cast[ptr cint](result.intWork.data)
  param.realWork = result.realWork.data
  if myRank == 0: param.display_params
  let ret = param.run(result.vals,
                      asarray[complex[float]](result.vecs.data)[],
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
  if param.locking != 0 and param.intWork != nil and param.intWork[] == 1:
    echo "\nA locking problem has occurred."
    echo "Some eigenpairs do not have a residual norm less than the tolerance."
    echo "However, the subspace of evecs is accurate to the required tolerance."
  case param.primme.dynamicMethodSwitch:
  of -1: echo "Recommended method for next run 1st stage: DEFAULT_MIN_MATVECS"
  of -2: echo "Recommended method for next run 1st stage: DEFAULT_MIN_TIME"
  of -3: echo "Recommended method for next run 1st stage: DYNAMIC (close call)"
  else: discard
  case param.primmeStage2.dynamicMethodSwitch:
  of -1: echo "Recommended method for next run 2nd stage: DEFAULT_MIN_MATVECS"
  of -2: echo "Recommended method for next run 2nd stage: DEFAULT_MIN_TIME"
  of -3: echo "Recommended method for next run 2nd stage: DYNAMIC (close call)"
  else: discard

export primme
export qexPrimme

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
  var opInfo = newOpInfo(s.addr)
  var pp = lo.primmeSVDInitialize(opInfo)
  var pevs = pp.run
  for i in 0..<pp.initSize:
    echo "Sval[",i,"]: ",pevs.vals[i].ff," rnorm: ",pevs.rnorms[i].ff
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
          check eigenResults[i] ~= pevs.vals[i]
          check eigenResults[i] ~= evals0[i].sv
          if not(eigenResults[i] ~= pevs.vals[i]) or
             not(eigenResults[i] ~= evals0[i].sv):
            echo "expect: ", $eigenResults[i]
            echo "primme: ", pevs.vals[i].ff
            echo "hisqev: ", evals0[i].sv.ff
  qexFinalize()
