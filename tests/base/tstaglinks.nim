## Adjoint identity of the link representation, including quantized caches.
## Promote the FP32 cache to FP64 before applying D: FP32 stencil roundoff
## would otherwise mask an inconsistent forward/backward representation.
import testutils
import qex, physics/qcdTypes, physics/stagGpu
import gauge/gaugeGpu
import backend/accel
import std/[math, strformat]

qexInit()
const V = static(VLEN)
let lo = intSeqParam("lat", latticeFromLocalLattice(@[4,4,4,8], nRanks)).newLayout
let n = lo.nSites
var rng = lo.newRngField(MRG32k3a, 987654321'u)
var g, gs, signs = lo.newGauge()
signs.unit()
threads:
  signs.setBC
  threadBarrier()
  signs.stagPhase
let sg = stagSigns(signs)
var gg = newGpuGauge(lo)
var u, v = lo.ColorVectorD()
threads:
  u.gaussian rng
  v.gaussian rng

proc fresh(unitary: bool) =
  g.random rng
  threads:
    g.projectSU
    if unitary:
      for mu in 0..<g.len:
        let z = newComplex(cos(0.17*(mu+1).float), sin(0.17*(mu+1).float))
        for i in g[mu]: g[mu][i] := z * g[mu][i]
    for mu in 0..<g.len: gs[mu] := g[mu]
    threadBarrier()
    gs.setBC
    threadBarrier()
    gs.stagPhase
  gg.upload(gg.u, g)

proc promote(d: StagGpu[V,float64]; s: StagGpu[V,float32]) =
  template copy(a,b: untyped; count: int) =
    let dst = a
    let src = b
    gpuFor(i, count): dst[i] = float(src[i])
  copy(d.lf, s.lf, 4*s.nl*s.n)
  if s.lb != nil: copy(d.lb, s.lb, 4*s.nl*s.n)
  for q in 0..1:
    if s.lh[q] != nil: copy(d.lh[q], s.lh[q], s.nl*s.ex[1-q].nrecv)

proc defect(s: StagGpu[V,float64]): float =
  let du = s.vec[0]
  let dv = s.vec[1]
  let su = s.vec[5]
  let sv = s.vec[6]
  gpuMemCpyToGpu(su, addr u[0], 6*n*sizeof(float))
  gpuMemCpyToGpu(sv, addr v[0], 6*n*sizeof(float))
  s.applyMfull(du, su, 0.0)
  s.applyMfull(dv, sv, 0.0)
  var hu, hv, hdu, hdv = newSeq[float](6*n)
  copyMem(addr hu[0], addr u[0], 6*n*sizeof(float))
  copyMem(addr hv[0], addr v[0], 6*n*sizeof(float))
  gpuMemCpyToCpu(addr hdu[0], du, 6*n*sizeof(float))
  gpuMemCpyToCpu(addr hdv[0], dv, 6*n*sizeof(float))
  var a: array[6,float]
  for k in 0..<n:
    for c in 0..2:
      let r = (k div V)*6*V + 2*c*V + k mod V
      let i = r + V
      a[0] += hu[r]*hdv[r] + hu[i]*hdv[i] + hdu[r]*hv[r] + hdu[i]*hv[i]
      a[1] += hu[r]*hdv[i] - hu[i]*hdv[r] + hdu[r]*hv[i] - hdu[i]*hv[r]
      a[2] += hu[r]*hu[r] + hu[i]*hu[i]
      a[3] += hv[r]*hv[r] + hv[i]*hv[i]
      a[4] += hdu[r]*hdu[r] + hdu[i]*hdu[i]
      a[5] += hdv[r]*hdv[r] + hdv[i]*hdv[i]
  getDefaultComm().allReduce(addr a[0], a.len)
  sqrt(a[0]*a[0]+a[1]*a[1]) / (sqrt(a[2]*a[5])+sqrt(a[3]*a[4]))

suite "staggered cached links":
  for nl in [12,14,18]:
    for fw in [0,1]:
      test &"{nl} reals, forward {fw}":
        fresh(nl == 14)
        var sd = newStagGpu(gs, float64, nl, fw)
        var ss = newStagGpu(gs, float32, nl, fw)
        check sd.nl == nl
        check ss.nl == nl
        for update in 0..2:
          if update > 0:
            fresh(nl == 14)
            sd.setLinks(gg, sg)
            ss.setLinks(gg, sg)
          let dd = sd.defect()
          sd.promote(ss)
          let ds = sd.defect()
          echo &"links {nl} fwd {fw} update {update} FP64 {dd:.3e} promoted FP32 {ds:.3e}"
          check dd < 1e-13
          check ds < 1e-13
        sd.free()
        ss.free()
gpuFree(sg)
gg.free()
qexFinalize()
