import quda_milc_interface
import quda
import enum_quda

#import base
import layout
import physics/qcdTypes
import physics/stagD
#import solvers/cg

import os, times

when not defined(qudaDir):
  {.fatal:"Must define qudaDir to use QUDA.".}
when not defined(cudaLibDir):
  {.fatal:"Must define cudaLibDir to use QUDA.".}
const qudaDir {.strdefine.} = ""
const cudaLibDir {.strdefine.} = ""
const cudaLib = "-L" & cudaLibDir & " -lcudart -lcublas -lcufft -Wl,-rpath," & cudaLibDir & " -L" & cudaLibDir & "/stubs -lcuda"
{.passC: "-I" & qudaDir & "/include".}

const qmpDir {.strdefine.} = getEnv("QMPDIR")
const qioDir {.strdefine.} = getEnv("QIODIR")

when qioDir.len > 0:
  when qmpDir.len > 0:
    # Assume quda is built with QIO and QMP.
    #{.passL: qudaDir & "/lib/libquda.a -lstdc++ " & cudaLib & " -L" & qioDir & "/lib -lqio -llime -L" & qmpDir & "/lib -lqmp".}
    {.passL: "-L" & qudaDir & "/lib -lquda -lstdc++ " & cudaLib & " -L" & qioDir & "/lib -lqio -llime -L" & qmpDir & "/lib -lqmp".}
  else:
    # Assume QUDA is built with QIO.
    {.passL: qudaDir & "/lib/libquda.a -lstdc++ " & cudaLib & " -L" & qioDir & "/lib -lqio -llime".}
else:
  {.passL: qudaDir & "/lib/libquda.a -lstdc++ " & cudaLib.}

{.passL: "-Wl,-rpath," & qudaDir & "/lib".}

#proc cudaGetDeviceCount(n:ptr cint):cint {.importc,nodecl.}
#proc cudaGetDeviceCount*: int =
#  var c: cint
#  discard cudaGetDeviceCount(c.addr)
#  c.int

type
  D4ColorMatrix = array[4, DColorMatrix]
  D4LatticeColorMatrix = Field[1, D4ColorMatrix]
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
    initCommsGridQuda(qudaParam.rankGeom.len.cint, qudaParam.rankGeom[0].addr,
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

proc qudaSolveEE*(s:Staggered; r,t:Field; m:SomeNumber; sp: var SolverParams) =
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
    fatlink: pointer = g1.s.data
    longlink: pointer = g3.s.data
    srcGpu: pointer = t1.s.data
    destGpu: pointer = r1.s.data
  invargs.maxIter = sp.maxits.cint
  invargs.evenodd = QUDA_EVEN_PARITY
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
        r{i}[a].re := r1[ri1.index][a].re
        r{i}[a].im := r1[ri1.index][a].im
  toc("QUDA teardown")

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
