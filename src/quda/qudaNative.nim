import quda

import os

# TODO select cpp when import this file?

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
    {.passL: " -L" & qudaDir & "/lib -lquda -Wl,-rpath," & qudaDir & "/lib " & cudaLib & " -L" & qioDir & "/lib -lqio -llime -L" & qmpDir & "/lib -lqmp".}
  else:
    # Assume QUDA is built with QIO.
    {.passL: " -L" & qudaDir & "/lib -lquda -Wl,-rpath," & qudaDir & "/lib " & cudaLib & " -L" & qioDir & "/lib -lqio -llime".}
else:
  {.passL: " -L" & qudaDir & "/lib -lquda -Wl,-rpath," & qudaDir & "/lib " & cudaLib.}


when isMainModule:
  import qex
  letParam:
    gaugefile = ""
    lat =
      if existsFile(gaugefile):
        getFileLattice gaugefile
      else:
        if gaugefile.len > 0:
          qexWarn "Nonexistent gauge file: ", gaugefile
        @[8,8,8,8]
  qexInit()
  echoParams()
  echo "rank ", myRank, "/", nRanks
  threads: echo "thread ", threadNum, "/", numThreads

  let
    lo = lat.newLayout
  var
    rgC: array[4,cint]

  proc qudaCommsMap(coords0: ptr cint; fdata: pointer): cint {.cdecl.} =
    let pl = cast[ptr type(lo)](fdata)
    let coords = cast[ptr UncheckedArray[cint]](coords0)
    let r = pl[].rankFromRankCoords(coords)
    r.cint
  for i in 0..<4: rgC[i] = lo.rankGeom[i].cint
  initCommsGridQuda(rgC.len.cint, rgC[0].addr,
                    qudaCommsMap, unsafeAddr(lo))
  initQuda(0)
  endQuda()
  qexFinalize()
  echoTimers()
