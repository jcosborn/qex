import quda_milc_interface
import quda
import enum_quda

import base
import physics/qcdTypes
import physics/stagD
import solvers/cg
import layout

import times

when not defined(qudaDir):
  {.fatal:"Must define qudaDir to use QUDA.".}
when not defined(cudaLibDir):
  {.fatal:"Must define cudaLibDir to use QUDA.".}
const qudaDir {.strdefine.} = ""
const cudaLibDir {.strdefine.} = ""
const cudaLib = "-L" & cudaLibDir & " -lcudart -lcufft -Wl,-rpath," & cudaLibDir
{.passC: "-I" & qudaDir & "/include".}
{.passL: qudaDir & "/lib/libquda.a -lstdc++ " & cudaLib.}

type
  D4ColorMatrix = array[4, DColorMatrix]
  D4LatticeColorMatrix = Field[1, D4ColorMatrix]
  QudaParam = object
    ## For the global quda parameter.
    initialized: bool
    initArg: QudaInitArgs_t
    physGeom,rankGeom: array[4,cint]    ## 4D only, used in QudaInitArgs_t
    layout: Layout[1]
    longlinkG: D4LatticeColorMatrix  ## all 0 for staggered

var qudaParam: QudaParam    ## Global quda parameter.

proc qudaInit* =
  ## Just to initialize the global parameter.
  qudaParam.initArg.layout.device = 0    # Single GPU per rank.
  qudaParam.initArg.layout.latsize = qudaParam.physGeom[0].addr
  qudaParam.initArg.layout.machsize = qudaParam.rankGeom[0].addr
  qudaParam.initArg.verbosity = QUDA_SUMMARIZE
  qudaParam.initialized = false

qexGlobalInitializers.add qudaInit
qexGlobalFinalizers.add qudaFinalize

proc qudaSetup*(l:Layout, verbosity = QUDA_SUMMARIZE):Layout[1] =
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
    qudaInit(qudaParam.initArg)
    qudaParam.layout = l.physGeom.newLayout 1
    qudaParam.longlinkG.new qudaParam.layout
    threads:
      for i in qudaParam.longlinkG:
        forO mu, 0, 3:
          qudaParam.longlinkG[i][mu] := 0   # zero out for hacking asqtad to do naive
    qudaParam.initialized = true
  qudaParam.layout

proc qudaSolveEE*(s:Staggered; r,t:Field; m:SomeNumber; sp:SolverParams) =
  tic()
  let lo1 = r.l.qudaSetup
  toc("QUDA one time setup")
  var
    # t2: float
    t1, r1: DLatticeColorVector
    g1: D4LatticeColorMatrix
  t1.new lo1
  r1.new lo1
  g1.new lo1
  toc("QUDA alloc")
  var
    invargs: QudaInvertArgs_t
    precision = 2   # 2 - double, 1 - single
    res = sqrt sp.r2req
    relRes = 0
    u0 = 1.0
    rres:cdouble = 0.0
    rrelRes:cdouble = 0.0
    iters:cint = 1
    fatlink: pointer = g1.s.data
    longlink: pointer = qudaParam.longlinkG.s.data
    srcGpu: pointer = t1.s.data
    destGpu: pointer = r1.s.data
  invargs.maxIter = sp.maxits.cint
  invargs.evenodd = QUDA_EVEN_PARITY
  invargs.mixedPrecision = 1    # 0: NO, 1: YES
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
  # echo "input norm2: ",t2
  toc("QUDA setup")
  # FIX ME and FIX QUDA interface: this is for asqtad, we use zero longlink
  qudaInvert(precision.cint, precision.cint,   # host, QUDA
    m.cdouble, invargs, res.cdouble, relRes.cdouble,
    fatlink, longlink, u0.cdouble, srcGpu, destGpu,
    rres.addr, rrelRes.addr, iters.addr)
  toc("QUDA invert")
  threads:
    for i in r.sites:
      var cv: array[4,cint]
      r.l.coord(cv,(r.l.myRank,i))
      let ri1 = lo1.rankIndex(cv)
      # assert(ri1.rank == r.l.myRank)
      forO a, 0, 2:
        r{i}[a].re := r1[ri1.index][a].re
        r{i}[a].im := r1[ri1.index][a].im
  toc("QUDA teardown")

when isMainModule:
  import qex
  import gauge
  qexInit()
  #var lat = [4,4,4,4]
  #var lat = [8,8,8,8]
  var lat = [16,8,4,32]
  threads:
    echo "thread ", threadNum, "/", numThreads
  var
    lo = lat.newLayout
    src = lo.ColorVector()
    dest = lo.ColorVector()
    destG = lo.ColorVector()
    r = lo.ColorVector()
    g = lo.newGauge
    rng = RngMilc6.newRNGField lo
  g.random rng
  threads:
    g.setBC
    g.stagPhase
    dest := 0
    destG := 0
    src.z4 rng
  var s = g.newStag
  var m = 0.0123
  var res = 1e-12
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

  s.solve(destG, src, m, res)

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
