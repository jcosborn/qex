import quda_milc_interface
import quda
import enum_quda

#import base
import layout
import physics/qcdTypes
import physics/stagD
import gauge/gaugeAction
#import solvers/cg

import os, times

when not defined(qudaDir):
  {.fatal:"Must define qudaDir to use QUDA.".}
const qudaDir {.strdefine.} = ""
{.passC: "-I" & qudaDir & "/include".}

const cudaLibDir {.strdefine.} = ""
const cudaMathLibDir {.strdefine.} = ""
when cudaLibDir != "":
  const cudaLib0 = "-L" & cudaLibDir & " -Wl,-rpath," & cudaLibDir
  const cudaLib1 =
    when cudaMathLibDir == "": cudaLib0
    else: cudaLib0 & " -L" & cudaMathLibDir & " -Wl,-rpath," & cudaMathLibDir
  const cudaLib = cudaLib1 & " -lcudart -lcublas -lcufft -L" & cudaLibDir & "/stubs -lcuda"
  {.passL: "-L" & qudaDir & "/lib -lquda -lstdc++ " & cudaLib .}
  {.passL: "-Wl,-rpath," & qudaDir & "/lib".}

const nvhpcDir {.strdefine.} = ""
when nvhpcDir != "":
  const cudaLib = "-L" & nvhpcDir & "/cuda/lib64 -L" & nvhpcDir & "/math_libs/lib64 -lcudart -lcublas -lcufft -Wl,-rpath," & cudaLibDir & " -L" & cudaLibDir & "/stubs -lcuda"
  {.passL: "-L" & qudaDir & "/lib -lquda -lstdc++ " & cudaLib .}
  {.passL: "-Wl,-rpath," & qudaDir & "/lib".}

when cudaLibDir=="" and nvhpcDir=="":
  {.passL: "-L" & qudaDir & "/lib -lquda -lstdc++ ".}
  {.passL: "-Wl,-rpath," & qudaDir & "/lib".}

const qmpDir {.strdefine.} = getEnv("QMPDIR")
const qioDir {.strdefine.} = getEnv("QIODIR")

when qioDir.len > 0:
  {.passL: "-L" & qioDir & "/lib -lqio -llime" .}

when qmpDir.len > 0:
  {.passL: "-L" & qmpDir & "/lib -lqmp".}

#proc cudaGetDeviceCount(n:ptr cint):cint {.importc,nodecl.}
#proc cudaGetDeviceCount*: int =
#  var c: cint
#  discard cudaGetDeviceCount(c.addr)
#  c.int

type
  D4ColorMatrix = array[4, DColorMatrix]
  D4LatticeColorMatrix = Field[1, D4ColorMatrix]
  D4Mom = array[4, array[10, float]]
  D4LatticeMom = Field[1, D4Mom]
  QudaParam = object
    ## For the global quda parameter.
    initialized: bool
    initArg: QudaInitArgs_t
    physGeom,rankGeom: array[4,cint]    ## 4D only, used in QudaInitArgs_t
    layout: Layout[1]
    #longlinkG: D4LatticeColorMatrix  ## all 0 for staggered

var qudaParam: QudaParam    ## Global quda parameter.

proc qudaInit* =
  ## Just to initialize the global parameter.
  ## Assumes single GPU per rank.
  #let n = max(1,cudaGetDeviceCount())
  #echo "cudaGetDeviceCount: ",n
  #let n = 1
  #qudaParam.initArg.layout.device = cint(myrank mod n)
  #qudaParam.initArg.layout.device = -1.cint
  qudaParam.initArg.layout.device = 0.cint
  qudaParam.initArg.layout.latsize = qudaParam.physGeom[0].addr
  qudaParam.initArg.layout.machsize = qudaParam.rankGeom[0].addr
  #qudaParam.initArg.verbosity = QUDA_SUMMARIZE
  qudaParam.initArg.verbosity = QUDA_SILENT
  qudaParam.initialized = false

qexGlobalInitializers.add qudaInit
qexGlobalFinalizers.add qudaFinalize

var qudaLayout: Layout[1]

proc setQudaLayout*(l: Layout) =
  qudaLayout = l

proc getQudaLayout*(): Layout[1] =
  qudaLayout

proc qudaSetup*(l: Layout, verbosity = QUDA_SILENT): Layout[1] =
  ## Actually initialize QUDA given the specific layout.
  var updated = false
  template update(a,b:typed) =
    if a != b:
      updated = true
      a = b
  qudaParam.initArg.verbosity.update verbosity
  for i in 0..3:
    qudaParam.physGeom[i].update l.physGeom[i].cint
    qudaParam.rankGeom[i].update l.rankGeom[i].cint
  if updated or (not qudaParam.initialized):
    if qudaParam.initialized: qudaFinalize()
    proc qudaCommsMap(coords0: ptr ConstInt; fdata: pointer): cint {.cdecl.} =
      let pl = cast[ptr type(l)](fdata)
      let coords = cast[ptr UncheckedArray[cint]](coords0)
      let r = pl[].rankFromRankCoords(coords)
      r.cint
    initCommsGridQuda(qudaParam.rankGeom.len.cint,
                      cast[ptr ConstInt](qudaParam.rankGeom[0].addr),
                      qudaCommsMap, unsafeAddr(l))
    qudaInit(qudaParam.initArg)
    qudaParam.layout = l.physGeom.newLayout 1
    #qudaParam.longlinkG.new qudaParam.layout
    #threads:
    #  for i in qudaParam.longlinkG:
    #    forO mu, 0, 3:
    #      qudaParam.longlinkG[i][mu] := 0   # zero out for hacking asqtad to do naive
    qudaParam.initialized = true
  setQudaLayout(qudaParam.layout)
  qudaParam.layout

proc setGaugeParam(gauge_param: var QudaGaugeParam) =
  gauge_param.location = QUDA_CPU_FIELD_LOCATION
  #gauge_param.type = QUDA_GENERAL_LINKS
  gauge_param.type = QUDA_SU3_LINKS
  gauge_param.gauge_order = QUDA_MILC_GAUGE_ORDER
  gauge_param.t_boundary = QUDA_PERIODIC_T
  for i in 0..3: gauge_param.X[i] = cint qudaLayout.localGeom[i]
  var cpu_prec = QUDA_DOUBLE_PRECISION
  gauge_param.cpu_prec = cpu_prec
  #var cuda_prec = QUDA_DOUBLE_PRECISION
  var cuda_prec = QUDA_SINGLE_PRECISION
  gauge_param.cuda_prec = cuda_prec
  gauge_param.cuda_prec_sloppy = cuda_prec
  gauge_param.cuda_prec_precondition = cuda_prec
  gauge_param.cuda_prec_eigensolver = cuda_prec
  #var link_recon = QUDA_RECONSTRUCT_NO
  var link_recon = QUDA_RECONSTRUCT_12
  gauge_param.reconstruct = link_recon
  gauge_param.reconstruct_sloppy = link_recon
  gauge_param.reconstruct_precondition = link_recon
  gauge_param.reconstruct_eigensolver = link_recon
  gauge_param.reconstruct_refinement_sloppy = link_recon
  #gauge_param.scale = 1.0
  gauge_param.scale = 0.0
  gauge_param.anisotropy = 1.0
  gauge_param.tadpole_coeff = 1.0
  gauge_param.ga_pad = 0
  gauge_param.site_ga_pad = 0
  gauge_param.mom_ga_pad = 0
  gauge_param.llfat_ga_pad = 0
  gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO
  gauge_param.use_resident_mom = 0
  gauge_param.make_resident_mom = 0
  gauge_param.return_result_mom = 1
  gauge_param.overwrite_mom = 0
  #gauge_param.struct_size = sizeof(gauge_param)

proc qudaSolveXX*(s:Staggered; r,t:Field; m:SomeNumber; sp: var SolverParams;
                  parEven = true) =
  tic()
  let lo1 = r.l.qudaSetup
  toc("QUDA one time setup")
  var
    # t2: float
    t1, r1: DLatticeColorVector
    g1: D4LatticeColorMatrix
    g3: D4LatticeColorMatrix
  t1.new lo1
  r1.new lo1
  g1.new lo1
  g3.new lo1
  toc("QUDA alloc")
  var
    invargs: QudaInvertArgs_t
    precision = 2   # 2 - double, 1 - single
    res = sqrt sp.r2req
    relRes = 0
    rres:cdouble = 0.0
    rrelRes:cdouble = 0.0
    iters:cint = 1
    fatlink: pointer = g1.dataPtr
    longlink: pointer = g3.dataPtr
    srcGpu: pointer = t1.dataPtr
    destGpu: pointer = r1.dataPtr
  invargs.maxIter = sp.maxits.cint
  invargs.evenodd = (if parEven: QUDA_EVEN_PARITY else: QUDA_ODD_PARITY)
  invargs.mixedPrecision = case sp.sloppySolve:
    of SloppyNone: 0
    of SloppySingle: 1
    of SloppyHalf: 2
  threads:
    # t2 = t.norm2
    for i in r.sites:
      var cv: array[4,cint]
      r.l.coord(cv,(r.l.myRank,i))
      let ri1 = lo1.rankIndex(cv)
      # assert(ri1.rank == r.l.myRank)
      forO a, 0, 2:
        t1[ri1.index][a].re := t{i}[a].re
        t1[ri1.index][a].im := t{i}[a].im
    for i in r.sites:
      var cv: array[4,cint]
      r.l.coord(cv,(r.l.myRank,i))
      let ri1 = lo1.rankIndex(cv)
      # assert(ri1.rank == r.l.myRank)
      forO a, 0, 2:
        r1[ri1.index][a].re := r{i}[a].re
        r1[ri1.index][a].im := r{i}[a].im
    if s.g.len == 4: # plain staggered
      longlink = nil
      for i in r.sites:
        var cv: array[4,cint]
        r.l.coord(cv,(r.l.myRank,i))
        let ri1 = lo1.rankIndex(cv)
        # assert(ri1.rank == r.l.myRank)
        forO mu, 0, 3:
          forO a, 0, 2:
            forO b, 0, 2:
              g1[ri1.index][mu][a,b].re := s.g[mu]{i}[a,b].re
              g1[ri1.index][mu][a,b].im := s.g[mu]{i}[a,b].im
    elif s.g.len == 8: # Naik staggered
      for i in r.sites:
        var cv: array[4,cint]
        r.l.coord(cv,(r.l.myRank,i))
        let ri1 = lo1.rankIndex(cv).index
        # assert(ri1.rank == r.l.myRank)
        forO mu, 0, 3:
          forO a, 0, 2:
            forO b, 0, 2:
              g1[ri1][mu][a,b].re := s.g[2*mu]{i}[a,b].re
              g1[ri1][mu][a,b].im := s.g[2*mu]{i}[a,b].im
              g3[ri1][mu][a,b].re := s.g[2*mu+1]{i}[a,b].re
              g3[ri1][mu][a,b].im := s.g[2*mu+1]{i}[a,b].im
    else:
      echo "unknown s.g.len: ", s.g.len
      quit(-1)
  # echo "input norm2: ",t2
  toc("QUDA setup")
  qudaInvert(precision.cint, precision.cint,   # host, QUDA
    m.cdouble, invargs, res.cdouble, relRes.cdouble,
    fatlink, longlink, srcGpu, destGpu,
    rres.addr, rrelRes.addr, iters.addr)
  toc("QUDA invert")
  sp.iterations = iters.int
  threads:
    for i in r.sites:
      var cv: array[4,cint]
      r.l.coord(cv,(r.l.myRank,i))
      let ri1 = lo1.rankIndex(cv)
      # assert(ri1.rank == r.l.myRank)
      forO a, 0, 2:
        r{i}[a].re = r1[ri1.index][a].re
        r{i}[a].im = r1[ri1.index][a].im
  toc("QUDA teardown")

proc qudaSolveEE*(s:Staggered; r,t:Field; m:SomeNumber; sp: var SolverParams) =
  qudaSolveXX(s, r, t, m, sp, parEven=true)

proc qudaSolveOO*(s:Staggered; r,t:Field; m:SomeNumber; sp: var SolverParams) =
  qudaSolveXX(s, r, t, m, sp, parEven=false)

#  Compute the gauge force and update the mometum field
#
#  @param mom The momentum field to be updated
#  @param sitelink The gauge field from which we compute the force
#  @param input_path_buf[dim][num_paths][path_length]
#  @param path_length One less that the number of links in a loop (e.g., 3 for a staple)
#  @param loop_coeff Coefficients of the different loops in the Symanzik action
#  @param num_paths How many contributions from path_length different "staples"
#  @param max_length The maximum number of non-zero of links in any path in the action
#  @param dt The integration step size (for MILC this is dt*beta/3)
#  @param param The parameters of the external fields and the computation settings
#proc computeGaugeForceQuda*(mom: pointer; sitelink: pointer;
#                           input_path_buf: ptr ptr ptr cint; path_length: ptr cint;
#                           loop_coeff: ptr cdouble; num_paths: cint;
#                           max_length: cint; dt: cdouble;
#                           qudaGaugeParam: ptr QudaGaugeParam): cint {.
proc qudaGaugeForce*[G,F](c: GaugeActionCoeffs, g0: openArray[G], f0: openArray[F]) =
  tic()
  let g = cast[ptr UncheckedArray[G]](unsafeAddr(g0[0]))
  let f = cast[ptr UncheckedArray[F]](unsafeAddr(f0[0]))
  let lo = g[0].l
  let lo1 = lo.qudaSetup
  toc("QUDA one time setup")
  var
    g1: D4LatticeColorMatrix
    #f1: D4LatticeColorMatrix
    f1: D4LatticeMom
  g1.new lo1
  f1.new lo1
  toc("QUDA alloc")
  var
    mom = pointer f1.dataPtr
    sitelink = pointer g1.dataPtr
    input_path_buf: array[4, ptr ptr cint]
    input_path_buf2: array[4, array[6, ptr cint]]
    path_length = [3.cint, 3, 3, 3, 3, 3]
    cp = c.plaq / 3.0
    loop_coeff = [cp, cp, cp, cp, cp, cp]
    num_paths = cint 6
    max_length = cint 3
    dt = cdouble 1
    qudaGaugeParam = newQudaGaugeParam()
    plaq_path: array[4,array[6,array[3,cint]]]
    plaq_path2 = [[[1, 7, 6], [6, 7, 1], [2, 7, 5], [5, 7, 2], [3, 7, 4], [4, 7, 3]],
                  [[2, 6, 5], [5, 6, 2], [3, 6, 4], [4, 6, 3], [0, 6, 7], [7, 6, 0]],
                  [[3, 5, 4], [4, 5, 3], [0, 5, 7], [7, 5, 0], [1, 5, 6], [6, 5, 1]],
                  [[0, 4, 7], [7, 4, 0], [1, 4, 6], [6, 4, 1], [2, 4, 5], [5, 4, 2]]]
  for mu in 0..3:
    input_path_buf[mu] = addr input_path_buf2[mu][0]
    for i in 0..<num_paths:
      input_path_buf2[mu][i] = addr plaq_path[mu][i][0]
      for j in 0..2:
        plaq_path[mu][i][j] = cint plaq_path2[mu][i][j]
  setGaugeParam qudaGaugeParam
  threads:
    for i in lo.sites:
      var cv: array[4,cint]
      lo.coord(cv,(lo.myRank,i))
      let ri1 = lo1.rankIndex(cv)
      # assert(ri1.rank == r.l.myRank)
      for mu in 0..3:
        g1[ri1.index][mu] := g[mu]{i}
        #f1[ri1.index][mu] := 0
  discard computeGaugeForceQuda(mom, sitelink, addr input_path_buf[0],
                                addr path_length[0], addr loop_coeff[0], num_paths,
                                max_length, dt, addr qudaGaugeParam)
  #let s2 = 0.70710678118654752440;  # sqrt(1/2)
  #let s3 = 0.57735026918962576450;  # sqrt(1/3)
  threads:
    for i in lo.sites:
      var cv: array[4,cint]
      lo.coord(cv,(lo.myRank,i))
      let ri1 = lo1.rankIndex(cv)
      # assert(ri1.rank == r.l.myRank)
      for mu in 0..3:
        #f[mu]{i} := f1[ri1.index][mu]
        let t = getPtr f1[ri1.index][mu]
        f[mu]{i}[0,1].set -t[0], -t[1]
        f[mu]{i}[1,0].set  t[0], -t[1]
        f[mu]{i}[0,2].set -t[2], -t[3]
        f[mu]{i}[2,0].set  t[2], -t[3]
        f[mu]{i}[1,2].set -t[4], -t[5]
        f[mu]{i}[2,1].set  t[4], -t[5]
        f[mu]{i}[0,0].set     0, -t[6]
        f[mu]{i}[1,1].set     0, -t[7]
        f[mu]{i}[2,2].set     0, -t[8]

when isMainModule:
  import qex
  import gauge
  import physics/stagSolve
  qexInit()
  var defaultLat = @[8,8,8,8]
  defaultSetup()
  echo "rank ", myRank, "/", nRanks
  threads:
    echo "thread ", threadNum, "/", numThreads
  var
    #lo = lat.newLayout
    src = lo.ColorVector()
    dest = lo.ColorVector()
    destG = lo.ColorVector()
    r = lo.ColorVector()
    #g = lo.newGauge
    rng = RngMilc6.newRNGField lo
  if fn == "":
    g.random rng
  threads:
    g.setBC
    g.stagPhase
    dest := 0
    destG := 0
    src.z4 rng
  var s = g.newStag
  let m = floatParam("mass", 0.02)
  let res = floatParam("cg_prec", 1e-9)
  let sloppy = intParam("sloppy", 2).SloppyType
  s.solve(dest, src, m, res, cpuonly = true)
  var n2: tuple[a,e,o:float]
  threads:
    echo "src.norm2: ", src.norm2
    echo "src.even: ", src.even.norm2
    echo "src.odd: ", src.odd.norm2
    n2.a = dest.norm2
    n2.e = dest.even.norm2
    n2.o = dest.odd.norm2
    echo "dest.norm2: ", n2.a
    echo "dest.even: ", n2.e
    echo "dest.odd: ", n2.o
    s.D(r, dest, m)
    threadBarrier()
    r := src - r
    threadBarrier()
    echo "r.norm2: ", r.norm2
    echo "r.even: ", r.even.norm2
    echo "r.odd: ", r.odd.norm2

  s.solve(destG, src, m, res, sloppySolve = sloppy)
  s.solve(destG, src, m, res, sloppySolve = sloppy)
  s.solve(destG, src, m, res, sloppySolve = sloppy)

  threads:
    echo "destG.norm2: ", destG.norm2
    echo "destG.even: ", destG.even.norm2
    echo "destG.odd: ", destG.odd.norm2
    s.D(r, destG, m)
    threadBarrier()
    r := src - r
    threadBarrier()
    echo "r.norm2: ", r.norm2
    echo "r.even: ", r.even.norm2
    echo "r.odd: ", r.odd.norm2
    r := destG - dest
    threadBarrier()
    echo "gpu-cpu: ", r.norm2 / n2.a
    echo "gpu-cpu:even ", r.even.norm2 / n2.e
    echo "gpu-cpu:odd ", r.odd.norm2 / n2.o

  qexFinalize()
  echoTimers()
