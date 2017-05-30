import primme
import base, comms/qmp, physics/qcdTypes, physics/stagD

type
  OpInfo[S,F] = object
    s:ptr S
    m:float
    x,y,tmp:F
proc newOpInfo[S](s:ptr S, m:float):auto =
  type F = type(s.g[0].l.ColorVector())
  var r:OpInfo[S,F]
  r.s = s
  r.m = m
  r.tmp.new(s.g[0].l)
  r.x.new
  r.x[] = r.tmp[]
  r.y.new
  r.y[] = r.tmp[]
  return r
proc makeX(op:OpInfo, x:ptr float) = op.x.s.data = cast[type(op.x.s.data)](x)
proc makeY(op:OpInfo, y:ptr float) = op.y.s.data = cast[type(op.y.s.data)](y)
proc applyD(op:ptr OpInfo; x,y:ptr float) =
  # x in, y out
  let op = op[]
  op.makeX x
  op.makeY y
  op.s[].D(op.tmp, op.x, op.m)
  op.s[].Ddag(op.y, op.tmp, op.m)
  
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

when isMainModule:
  import qex, gauge
  qexInit()
  var lat = [4,4,4,4]
  threads:
    echo "thread ", threadNum, "/", numThreads
  var
    lo = lat.newLayout
    g = newSeq[type(lo.ColorMatrix())](lat.len)
  for i in 0..<lat.len:
    g[i] = lo.Colormatrix()
    threads: g[i] := 1
  gauge.random(g)
  for mu in 0..<lat.len:
    #var t, s: DColorMatrixV   # FIXME: get vectorized code to work with projectU
    var t, s: DColorMatrix
    tfor i, 0..<lo.nSites:
      for a in 0..2:
        for b in 0..2:
          s[a,b].re := g[mu]{i}[a,b].re
          s[a,b].im := g[mu]{i}[a,b].im
      t.projectU s
      for a in 0..2:
        for b in 0..2:
          g[mu]{i}[a,b].re := t[a,b].re
          g[mu]{i}[a,b].im := t[a,b].im
  threads:
    g.setBC
    g.stagPhase
  var s = g.newStag
  var m = 0.1
  var opInfo = newOpInfo(s.addr, m)
  var pp = primme_initialize()
  block primmeSetup:
    pp.n = 3*lo.physVol
    pp.nLocal = 3*lo.nSites
    pp.globalSumReal = sumReal
    pp.matrixMatvec = matvec[type(opInfo)]
    pp.matrix = opInfo.addr
    pp.numEvals = 16
    pp.eps = 1e-8
    pp.target = primme_smallest
    let ret = pp.set_method PRIMME_DYNAMIC
    if 0 != ret:
      echo "Error: set_method returned with nonzero exit status: ", ret
      quit QuitFailure
    pp.display_params
  var
    evecs = newAlignedMemU[complex[float]](pp.n * pp.numEvals)
    evals = newseq[float](pp.numEvals)
    rnorms = newseq[float](pp.numEvals)
  block primmeRun:
    let ret = pp.run(evals, evecs.data, rnorms)
    if ret != 0:
      echo "Error: primme returned with nonzero exit status: ", ret
      quit QuitFailure
  block primmeReport:
    for i in 0..<pp.initSize:
      echo "Eval[",i,"]: ",evals[i].ff," rnorm: ",rnorms[i].ff
    echo " ",pp.initSize," eigenpairs converged"
    echo "Tolerance  : ",ff pp.aNorm*pp.eps
    echo "Iterations : ",pp.stats.numOuterIterations
    echo "Restarts   : ",pp.stats.numRestarts
    echo "Matvecs    : ",pp.stats.numMatvecs
    echo "Preconds   : ",pp.stats.numPreconds
    if pp.locking != 0 and pp.intWork != nil and pp.intWork[] == 1:
      echo "\nA locking problem has occurred."
      echo "Some eigenpairs do not have a residual norm less than the tolerance."
      echo "However, the subspace of evecs is accurate to the required tolerance."
    case pp.dynamicMethodSwitch:
    of -1: echo "Recommended method for next run: DEFAULT_MIN_MATVECS"
    of -2: echo "Recommended method for next run: DEFAULT_MIN_TIME"
    of -3: echo "Recommended method for next run: DYNAMIC (close call)"
    else: discard
  pp.free
