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

proc testPlaqForce(g:auto, force:static bool = true) =
  tic "testPlaqForce"
  let lo = g[0].l
  let nd = lo.nDim
  let hl = lo.makeHaloLayout([1,1,1,1],[1,1,1,1])
  toc "makeHaloLayout"
  type HM = HaloMap[type lo]
  let comm = getDefaultComm()
  var hm = newSeq[HM](nd)
  for mu in 0..<nd:
    # g[mu] needs offsets: +nu, -nu, -mu, -mu+nu
    var offsets = newSeq[seq[int32]](0)
    for nu in 0..<nd:
      if nu == mu: continue
      template addOffset(omu, onu: int) =
        var t = newSeq[int32](nd)
        t[mu] = omu
        t[nu] = onu
        offsets.add t
      addOffset( 0, 1)
      addOffset( 0,-1)
      addOffset(-1, 0)
      addOffset(-1, 1)
    hm[mu] = hl.makeHaloMap(comm, offsets)
  toc "makeHaloMap"
  var f = g.newOneOf
  var f2 = f.newOneOf
  type H = type makeHalo(hl, g[0])
  var h = newSeq[H](nd)
  for d in 0..<nd:
    h[d] = makeHalo(hl, g[d])
    h[d].gpuFlagsExcl {gmGpuWrite}
    #h[d].gpuFlagsIncl {gmCpuWriteOnce}
    g[d].gpuFlagsExcl {gmGpuWrite}
    #g[d].gpuFlagsIncl {gmCpuWriteOnce}
    f[d].gpuFlagsExcl {gmGpuRead}  # force doesn't need copyIn
  f.gpuFlagsIncl {gmCpuWriteOnce}
  f.gpuFlagsExcl {gmGpuWrite}
  toc "makeHalo"
  let gac = GaugeActionCoeffs(plaq:3)  # coeff is devided by Nc
  let nc = g[0][0].getNc
  let mm = nc*nc*(8*nc-2)
  let pm = 2*nc*nc
  when force:
    let flops = lo.nSites*(6*(7*mm+4*pm)+4*(mm+2*nc*nc+4*nc))
    gac.gaugeForce(g, f2)  # warmup
  else:
    let flops = lo.nSites*(6*(7*mm+4*pm))
    gac.gaugeActionDeriv(g, f2)  # warmup
  toc "create fields"
  pushGpuMemTag("testPlaqForce")
  var halotime = 0.0
  let gpuWaitFlops = flops
  for nreps in [2,10]:
    resetTimers()
    halotime = 0
    for rep in 0..<nreps:
      tic "rep"
      for d in 0..<nd:
        h[d].update hm[d], comm
      toc "update"
      onGpu(lo):
        template g(i:int):auto = h[i].field
        for s in gpuSites(g(0)):
          for mu in 0..<4:
            f[mu][s] := 0
          for mu in 1..<4:
            let fmu = hl.nbrFwd(mu, s)
            let bmu = hl.nbrBck(mu, s)
            for nu in 0..<mu:
              let fnu = hl.nbrFwd(nu, s)
              let bnu = hl.nbrBck(nu, s)
              let bmufnu = hl.nbrFwd(nu, bmu)
              let bnufmu = hl.nbrFwd(mu, bnu)
              let a = h[mu][fnu] * h[nu][fmu].adj
              f[mu][s] += g(nu)[s] * a
              f[nu][s] += g(mu)[s] * a.adj
              f[mu][s] += h[nu][bnu].adj * h[mu][bnu] * h[nu][bnufmu]
              f[nu][s] += h[mu][bmu].adj * h[nu][bmu] * h[mu][bmufnu]
          when force:
            for mu in 0..<4:
              let t = g(mu)[s] * f[mu][s].adj
              f[mu][s].projectTAH t
      halotime += getElapsedTime() / nreps
      toc("gpu",flops=flops)
  toc("gpu")
  var gadtime = getElapsedTime()
  when force:
    gac.gaugeForce(g, f2)  # warmup
  else:
    gac.gaugeActionDeriv(g, f2)  # warmup
  gadtime = getElapsedTime() - gadtime
  for mu in 0..<4:
    let nf = f[mu].norm2
    let nf2 = f2[mu].norm2
    f2[mu] -= f[mu]
    echo "nf: ", nf, "  nf2: ", nf2, "  diff2: ", f2[mu].norm2
  echo "haloTime: ", halotime, "  gadtime: ", gadtime, "  ratio: ", halotime/gadtime
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
  echo 6.0 * g.plaq
  toc "plaq"
  resetTimers()
  #testPlaqForce(g, false)
  testPlaqForce(g)
  echoProf()
  qexFinalize()
