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

proc testPlaq(g:auto) =
  tic "testPlaq"
  let lo = g[0].l
  let nd = lo.nDim
  let hl = lo.makeHaloLayout([1,1,1,1],[0,0,0,0])
  toc "makeHaloLayout"
  type HM = HaloMap[type lo]
  let comm = getDefaultComm()
  var hm = newSeq[HM](nd)
  for mu in 0..<nd:
    var offsets = newSeq[seq[int32]](0)
    for nu in 0..<nd:
      if nu == mu: continue
      var t = newSeq[int32](nd)
      t[nu] = 1
      offsets.add t
    hm[mu] = hl.makeHaloMap(comm, offsets)
  toc "makeHaloMap"
  type H = type makeHalo(hl, g[0])
  var h = newSeq[H](nd)
  for d in 0..<nd:
    h[d] = makeHalo(hl, g[d])
    h[d].gpuFlagsExcl {gmGpuWrite}
    h[d].gpuFlagsIncl {gmCpuWriteOnce}
    g[d].gpuFlagsExcl {gmGpuWrite}
    g[d].gpuFlagsIncl {gmCpuWriteOnce}
  toc "makeHalo"
  #var p = newSeq[typeof g[0]](6)
  #for i in 0..<6:
  #  p[i] = g[0].newOneOf
  #threads:
  #  for i in 0..<6:
  #    p[i] := 0
  var pl = newSeq[float](6)
  var gs = newGpuSum[array[6,float]](lo.nSites)
  pushGpuMemTag("testPlaq")
  toc "create fields"
  proc gpuSite(x: GpuField): auto =
    for i in gpuSites(x): return x[i]
  for nreps in [2,10]:
    resetTimers()
    for rep in 0..<nreps:
      tic "rep"
      for d in 0..<nd:
        h[d].update hm[d], comm
      toc "update"
      when false:
        threads:
          for i in g[0]:
            var k = 0
            for mu in 1..<4:
              let n0 = hl.neighborFwd[mu][i]
              for nu in 0..<mu:
                let n1 = hl.neighborFwd[nu][i]
                let a = g[mu][i] * h[nu][n0]
                let b = g[nu][i] * h[mu][n1]
                p[k][i] += a.adj * b
                inc k
      else:
        onGpu(lo):
          template g(i:int):auto = h[i].field
          var tpl: array[6,float]
          #var tpl: array[6,typeof redot(g(0).gpuSite,g(0).gpuSite)]
          for s in gpuSites(g(0)):
            var k = 0
            for mu in 1..<4:
              let smu = hl.nbrFwd(mu, s)
              for nu in 0..<mu:
                let snu = hl.nbrFwd(nu, s)
                let a = g(mu)[s] * h[nu][smu]
                let b = g(nu)[s] * h[mu][snu]
                #let a = h[1][s]
                #let b = h[1][s]
                tpl[k] += redot(a, b).simdSum
                #tpl[k] += redot(a, b)
                inc k
          gs.reduce tpl
          #var tplf: array[6,float]
          #for k in 0..<6: tplf[k] = tpl[k].simdSum
          #gs.reduce tplf
      let nc = g[0][0].getNc
      let mm = nc*nc*(8*nc-2)
      let rd = 8*nc*nc
      toc("plaq",flops=lo.nSites*6*(2*mm+rd))
  #toc "p"
  #threads:
  for k in 0..<6:
    #pl[k] = p[k].trace.re
    pl[k] = gs.value[k]
  #toc "pl"
  rankSum pl
  let vf = 1.0/(g[0][0].nRows*lo.physVol)
  let ph = pl * vf
  let pp = 6.0 * g.plaq
  let d = ph - pp
  echo ph
  echo d
  echo "norm2 diff: ", sum(d*d)
  #echo pl * vf
  #echo 6.0 * g.plaq
  echo dumpGpuMem()

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
  #var r0 = lo.Real()
  #var cv1 = lo.ColorVector()
  #var cv2 = lo.ColorVector()
  echo 6.0 * g.plaq
  toc "plaq"
  #cv0.gaussian rng
  resetTimers()
  testPlaq(g)
  echoProf()
  qexFinalize()
