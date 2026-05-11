import qex
import physics/qcdTypes
import backend/[accel,cpugpu,cgfield]
import bench/commonBench
#import strformat
import macros
import base/metaUtils
import parseUtils
import sequtils, strutils
import comms/[halo, halogpu]
import std/decls # for byaddr

type
  StagHalo[L,F,T] = object
    hl*: HaloLayout[L]
    hmg*: seq[HaloMap[L]]
    hmf*: HaloMap[L]
    hg*: seq[Halo[L,F,T]]

proc newStagHalo*(g:auto, L,F,T: typedesc): StagHalo[L,F,T] =
  tic "newStagHalo"
  let lo = g[0].l
  let nd = lo.nDim
  let hl = lo.makeHaloLayout([1,1,1,1],[1,1,1,1])
  result.hl = hl
  toc "makeHaloLayout"
  type HM = HaloMap[type lo]
  let comm = getDefaultComm()
  var hmg = newSeq[HM](nd)
  for mu in 0..<nd:
    var offsets = newSeq[seq[int32]](0)
    var t = newSeq[int32](nd)
    t[mu] = -1
    offsets.add t
    hmg[mu] = hl.makeHaloMap(comm, offsets)
  result.hmg = hmg
  var offsets = newSeq[seq[int32]](0)
  for mu in 0..<nd:
    var t = newSeq[int32](nd)
    t[mu] = 1
    offsets.add t
    t[mu] = -1
    offsets.add t
  result.hmf = hl.makeHaloMap(comm, offsets)
  toc "makeHaloMap"
  type GH = type makeHalo(hl, g[0])
  var hg = newSeq[GH](nd)
  for d in 0..<nd:
    hg[d] = makeHalo(hl, g[d])
    hg[d].gpuFlagsExcl {gmGpuWrite}
    hg[d].gpuFlagsIncl {gmCpuWriteOnce}
    g[d].gpuFlagsExcl {gmGpuWrite}
    g[d].gpuFlagsIncl {gmCpuWriteOnce}
    #f[d].gpuFlagsExcl {gmGpuRead}  # force doesn't need copyIn
  #f.gpuFlagsIncl {gmCpuWriteOnce}
  #f.gpuFlagsExcl {gmGpuWrite}
  result.hg = hg
  toc "makeHalo"

proc newStagHalo*(g:auto): auto =
  newStagHalo(g, typeof(g[0].l), typeof(g[0]), typeof(g[0][0]))

proc newFermion*(sh: StagHalo, f: auto): auto =
  result = makeHalo(sh.hl, f)

proc updateFermion*(sh: StagHalo, hf: auto) =
  let comm = getDefaultComm()
  hf.update sh.hmf, comm

proc updateLinks*(sh: StagHalo) =
  let comm = getDefaultComm()
  for d in 0..<sh.hg.len:
    sh.hg[d].update sh.hmg[d], comm

proc stagD(sh: StagHalo, dst: auto, src: auto) =
  tic("stagD")
  dst.gpuFlagsExcl {gmGpuRead}
  src.gpuFlagsExcl {gmGpuWrite}
  dst.field.gpuFlagsExcl {gmGpuRead}
  src.field.gpuFlagsExcl {gmGpuWrite}
  let lo = src.field.l
  let hl = sh.hl
  var hg {.byaddr.} = sh.hg
  #var hg: openArray[Halo]
  #echo "sh.hg: ", cast[int](addr sh.hg[0])
  #echo "hg:    ", cast[int](addr hg[0])
  var dst = dst
  var src = src
  #threads:
  #  for s in lo:
  let gpuWaitFlops = lo.nSites*8*(66+6+6)
  onGpu(lo):
    for s in gpuSites(dst.field):
      dst[s] := 0
      for mu in 0..<4:
        #let fmu = hl.neighborFwd[mu][s]
        let fmu = hl.nbrFwd(mu, s)
        dst[s] += 0.5 * (hg[mu].field[s] * src[fmu])
      #for mu in 0..<4:
        #let bmu = hl.neighborBck[mu][s]
        let bmu = hl.nbrBck(mu, s)
        dst[s] -= 0.5 * (hg[mu][bmu].adj * src[bmu])
  toc("end",flops=gpuWaitFlops)

proc testStagHalo(g:auto, sd:auto, src,v1,v2: auto) =
  let sh = newStagHalo(g)
  let hsrc = sh.newFermion(src)
  let hv2 = sh.newFermion(v2)
  sh.updateLinks
  sh.updateFermion hsrc
  sh.stagD(hv2, hsrc)
  #toc "update"
  #var gadtime = getElapsedTime()
  #stagD(sdAll, v1, g, src, 1.0)
  #gadtime = getElapsedTime() - gadtime
  resetTimers()
  tic "testStagHalo"
  for i in 0..<10:
    threads:
      stagD(sd, v1, g, src, 0.0)
  for i in 0..<10:
    sh.updateLinks
    sh.updateFermion hsrc
    sh.stagD(hv2, hsrc)
  echo "v1: ", v1.norm2
  echo "v2: ", v2.norm2
  v2 -= v1
  echo "d2: ", v2.norm2
  #echo "nf: ", nf, "  nf2: ", nf2, "  diff2: ", f2[mu].norm2
  #echo "haloTime: ", halotime, "  gadtime: ", gadtime, "  ratio: ", halotime/gadtime
  echo dumpGpuMem()
  freeGpuMem(sh.hl)

when isMainModule:
  qexInit()
  tic("main")
  var defaultLat = @[4,4,4,4]
  defaultSetup()
  let nd = lo.nDim
  var seed = 987654321'u
  #var rng = newRngField(lo, RngMilc6, seed)
  var rng = newRngField(lo, MRG32k3a, seed)
  g.gaussian rng
  #g.unit
  toc "gaussian"
  echo 6.0 * g.plaq
  toc "plaq"
  var src = lo.ColorVector()
  var v1 = lo.ColorVector()
  var v2 = lo.ColorVector()
  var sdAll = initStagD(src, "all")
  threads:
    g.setBC
    threadBarrier()
    g.stagPhase
    threadBarrier()
    if myRank==0 and threadNum==0:
      when compiles(src[0].len):
        src{0}[0] := 1
      else:
        src{0} := 1
    threadBarrier()
    stagD(sdAll, v1, g, src, 0.0)
    threadBarrier()
    echo "src: ", src.norm2
    echo "v1:  ", v1.norm2
  resetTimers()
  testStagHalo(g, sdAll, src, v1, v2)
  echoProf()
  freeAllGpuMem()
  var gs = eval(toSingle(g))
  var sdAlls = toSingle(sdAll)
  var srcs = eval(toSingle(src))
  var v1s = eval(toSingle(v1))
  var v2s = eval(toSingle(v2))
  testStagHalo(gs, sdAlls, srcs, v1s, v2s)
  echoProf()
  qexFinalize()
