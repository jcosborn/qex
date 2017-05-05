import quda_milc_interface
import quda
import enum_quda

import base
import physics/qcdTypes
import physics/stagD
import solvers/cg

import times
import random

import os
when existsEnv("QUDADIR"):
  const qudaDir = getEnv("QUDADIR")
else:
  const homeDir = getHomeDir()
  const qudaDir = homeDir & "/lqcdM/build/quda"
when existsEnv("CUDADIR"):
  const cudaDir = getEnv("CUDADIR")
else:
  const cudaDir = "/usr/local/cuda/lib64"
const cudaLib = "-L" & cudaDir & " -lcudart -lcufft -Wl,-rpath," & cudaDir
{.passC: "-I" & qudaDir & "/include".}
{.passL: qudaDir & "/lib/libquda.a -lstdc++ " & cudaLib.}

type
  D4ColorMatrix = array[4, DColorMatrix]
  D4LatticeColorMatrix = Field[1, D4ColorMatrix]
proc qudaSolve*(s:Staggered; r,x:Field; m:SomeNumber; sp0:SolverParams) =
  var sp = sp0
  sp.subset.layoutSubset(r.l, sp.subsetName)
  var t = newOneOf r
  threads:
    s.eoReduce(t, x, m)
  let tpresetup = epochTime()
  let lo1 = r.l.physGeom.newLayout 1
  var
    t1, r1: DLatticeColorVector
    g0, g1: D4LatticeColorMatrix
  t1.new lo1
  r1.new lo1
  g0.new lo1
  g1.new lo1
  tfor i, 0..<g0.l.nSites:
    for mu in 0..3:
      g0[i][mu] := 0   # zero out for hacking asqtad to do naive
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
    longlink: pointer = g0.s.data
    srcGpu: pointer = t1.s.data
    destGpu: pointer = r1.s.data
  invargs.maxIter = sp.maxits.cint
  invargs.evenodd = QUDA_EVEN_PARITY   # QUDA_ODD_PARITY
  invargs.mixedPrecision = 0
  tfor i, 0..<r.l.nSites:
    var cv: array[4,cint]
    lo.coord(cv,(r.l.myRank,i))
    let ri1 = lo1.rankIndex(cv)
    assert(ri1.rank == r.l.myRank)
    for a in 0..2:
      t1[ri1.index][a].re := t{i}[a].re
      t1[ri1.index][a].im := t{i}[a].im
      r1[ri1.index][a].re := r{i}[a].re
      r1[ri1.index][a].im := r{i}[a].im
    for mu in 0..3:
      for a in 0..2:
        for b in 0..2:
          g1[ri1.index][mu][a,b].re := s.g[mu]{i}[a,b].re
          g1[ri1.index][mu][a,b].im := s.g[mu]{i}[a,b].im
  let secpresetup = epochTime() - tpresetup
  let tsolve = epochTime()
  # FIX ME and FIX QUDA interface: this is for asqtad, we use zero longlink
  qudaInvert(precision.cint, precision.cint,   # host, QUDA
    m.cdouble, invargs, res.cdouble, relRes.cdouble,
    fatlink, longlink, u0.cdouble, srcGpu, destGpu,
    rres.addr, rrelRes.addr, iters.addr)
  let secsolve = epochTime() - tsolve
  let tpostproc = epochTime()
  tfor i, 0..<r.l.nSites:
    var cv: array[4,cint]
    lo.coord(cv,(r.l.myRank,i))
    let ri1 = lo1.rankIndex(cv)
    assert(ri1.rank == r.l.myRank)
    for a in 0..2:
      r{i}[a].re := r1[ri1.index][a].re
      r{i}[a].im := r1[ri1.index][a].im
  let secpostproc = epochTime() - tpostproc
  echo "----- QUDA results -----"
  echo "rres: ", rres
  echo "rrelRes: ", rrelRes
  echo "iters: ", iters
  threads:
    echo "r1.norm2: ", r1.norm2
    echo "r1.even: ", r1.even.norm2
    echo "r1.odd: ", r1.odd.norm2
    echo "r.norm2: ", r.norm2
    echo "r.even: ", r.even.norm2
    echo "r.odd: ", r.odd.norm2
  echo "GPU presetup time: ", secpresetup
  echo "GPU solve time: ", secsolve
  echo "GPU postproc time: ", secpostproc
  threads:
    r[s.se.sub] := 4*r
    threadBarrier()
    s.eoReconstruct(r, x, m)
proc qudaSolve*(s:Staggered; r,x:Field; m:SomeNumber; res:float) =
  var sp = initSolverParams()
  sp.r2req = res
  sp.verbosity = 1
  qudaSolve(s, r, x, m, sp)

when isMainModule:
  import gauge
  import layout
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
    src := 0
    dest := 0
    destG := 0
  tfor i, 0..<lo.nSites:
    for a in 0..2:
      src{i}[a].re := (if random(1.0)>0.5: 1.0 else: -1.0)
      src{i}[a].im := (if random(1.0)>0.5: 1.0 else: -1.0)
  var s = g.newStag
  var m = 0.0123
  var res = 1e-12
  s.solve(dest, src, m, res)
  threads:
    echo "src.norm2: ", src.norm2
    echo "src.even: ", src.even.norm2
    echo "src.odd: ", src.odd.norm2
    echo "dest.norm2: ", dest.norm2
    echo "dest.even: ", dest.even.norm2
    echo "dest.odd: ", dest.odd.norm2
    s.D(r, dest, m)
    threadBarrier()
    r := src - r
    threadBarrier()
    echo "r.norm2: ", r.norm2
    echo "r.even: ", r.even.norm2
    echo "r.odd: ", r.odd.norm2

  var
    initArg: QudaInitArgs_t
    latGpu: array[4,cint]
    msGpu: array[4,cint]
  for i in 0..3:
    latGpu[i] = lat[i].cint
    msGpu[i] = 1   # SINGLE GPU machine for now
  #initArg.verbosity = QUDA_DEBUG_VERBOSE
  initArg.verbosity = QUDA_SUMMARIZE
  initArg.layout.device = 0   # single gpu
  initArg.layout.latsize = latGpu[0].addr
  initArg.layout.machsize = msGpu[0].addr
  qudaInit(initArg)
  s.qudaSolve(destG, src, m, res)
  qudaFinalize()

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
    echo "gpu-cpu: ", r.norm2
    echo "gpu-cpu:even ", r.even.norm2
    echo "gpu-cpu:odd ", r.odd.norm2

  qexFinalize()
