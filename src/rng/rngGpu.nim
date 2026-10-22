## The generators of a host RNG field on the GPU, stored like a field of a
## SIMD layout lo of V lanes as RNGFieldV: 32-bit word w of the generator
## of site k at (k div V)*W*V + w*V + k mod V, W words per generator.  Site k draws the
## numbers of the host generator at its coordinates, so the kernels below
## fill device fields as their host versions fill fields of lo.  Works for
## generators of plain data whose draws are inline procs, which then
## compile for the device: RngMilc6, MRG32k3a, Philox4x64, Threefry4x64.
import qex
import backend/accel
import std/math

type RngGpu*[V: static int; R] = object
  lo*: Layout[V]
  s*: ptr UncheckedArray[uint32]  # [n div V][W][V]

proc upload*[V: static int; R](g: RngGpu[V,R]; r: Field[1,R]) =
  ## the generators of r to g, site k of g.lo from the generator of its
  ## coordinates
  doAssert r.l.rankGeom == g.lo.rankGeom, "the RNG field needs the rank grid of the layout"
  let n = g.lo.nSites
  var h = newSeq[uint32](rngWords(R)*n)
  let hp = cast[ptr UncheckedArray[uint32]](addr h[0])
  for k in 0..<n: rngStore(hp, V, k, r[r.l.rankIndex(g.lo.coords, k).index])
  gpuMemCpyToGpu(g.s, hp, h.len*sizeof(uint32))

proc download*[V: static int; R](g: RngGpu[V,R]; r: Field[1,R]) =
  ## the generators of g to r
  doAssert r.l.rankGeom == g.lo.rankGeom, "the RNG field needs the rank grid of the layout"
  let n = g.lo.nSites
  var h = newSeq[uint32](rngWords(R)*n)
  let hp = cast[ptr UncheckedArray[uint32]](addr h[0])
  gpuMemCpyToCpu(hp, g.s, h.len*sizeof(uint32))
  for k in 0..<n:
    let a = rngLoad(hp, V, k, R)
    copyMem(addr r[r.l.rankIndex(g.lo.coords, k).index], unsafeAddr a, sizeof(R))  # RNGs other than RngMilc6 lack :=

proc newRngGpu*[V: static int; R](lo: Layout[V]; r: Field[1,R]): RngGpu[V,R] =
  ## the generators of r on the device, for fields of lo
  result.lo = lo
  result.s = cast[ptr UncheckedArray[uint32]](gpuMalloc(rngWords(R)*lo.nSites*sizeof(uint32)))
  result.upload r

proc upload*[V: static int; R](g: RngGpu[V,R]; r: RNGFieldV[V,R]) =
  ## the generators of r to g, a copy
  gpuMemCpyToGpu(g.s, addr r.s[0], r.s.len*sizeof(uint32))

proc download*[V: static int; R](g: RngGpu[V,R]; r: RNGFieldV[V,R]) =
  ## the generators of g to r, a copy
  gpuMemCpyToCpu(addr r.s[0], g.s, r.s.len*sizeof(uint32))

proc newRngGpu*[V: static int; R](r: RNGFieldV[V,R]): RngGpu[V,R] =
  ## the generators of r on the device
  result.lo = r.l
  result.s = cast[ptr UncheckedArray[uint32]](gpuMalloc(r.s.len*sizeof(uint32)))
  result.upload r

proc free*[V: static int; R](g: var RngGpu[V,R]) =
  gpuFree(g.s)
  g.s = nil

proc randomTAH*[V: static int; R](g: RngGpu[V,R]; p: ptr UncheckedArray[float]) =
  ## p as randomTAH of 4 fields of 3x3 matrices, p in the layout of
  ## GpuGauge.u: real e = 6a+2b of M_ab (e+1 its imaginary part) of
  ## direction mu at 18*mu*n + (k div V)*18V + e*V + k mod V
  let n = g.lo.nSites
  let s = g.s
  gpuFor(k, n):
    const s2 = 0.70710678118654752440  # sqrt(1/2)
    const s3 = 0.57735026918962576450  # sqrt(1/3)
    var a = rngLoad(s, V, k, R)
    for mu in 0..3:
      let o = 18*mu*n + (k div V)*(18*V) + k mod V
      let r3 = s2 * gaussian(a)
      let r8 = s2 * s3 * gaussian(a)
      let r01 = s2 * gaussian(a)
      let r02 = s2 * gaussian(a)
      let r12 = s2 * gaussian(a)
      let i01 = s2 * gaussian(a)
      let i02 = s2 * gaussian(a)
      let i12 = s2 * gaussian(a)
      p[o] = 0.0
      p[o + V] = r8 + r3
      p[o + 8*V] = 0.0
      p[o + 9*V] = r8 - r3
      p[o + 16*V] = 0.0
      p[o + 17*V] = -2*r8
      p[o + 2*V] = r01
      p[o + 3*V] = i01
      p[o + 6*V] = -r01
      p[o + 7*V] = i01
      p[o + 4*V] = r02
      p[o + 5*V] = i02
      p[o + 12*V] = -r02
      p[o + 13*V] = i02
      p[o + 10*V] = r12
      p[o + 11*V] = i12
      p[o + 14*V] = -r12
      p[o + 15*V] = i12
    rngStore(s, V, k, a)

proc gaussian*[V: static int; R](g: RngGpu[V,R]; x: ptr UncheckedArray[float]; ne: int) =
  ## x as gaussian of a field of ne reals per site in the SIMD layout, real
  ## c of site k at (k div V)*ne*V + c*V + k mod V
  let n = g.lo.nSites
  let s = g.s
  gpuFor(k, n):
    var a = rngLoad(s, V, k, R)
    let o = (k div V)*(ne*V) + k mod V
    for c in 0..<ne: x[o + c*V] = gaussian(a)
    rngStore(s, V, k, a)

proc u1*[V: static int; R](g: RngGpu[V,R]; x: ptr UncheckedArray[float]; ne: int) =
  ## x as u1 of a field of ne reals per site, in the layout of gaussian:
  ## complex c = exp(2 pi i u) for each pair of reals, u uniform, as the u1
  ## of the host generators other than RngFuel
  let n = g.lo.nSites
  let s = g.s
  gpuFor(k, n):
    var a = rngLoad(s, V, k, R)
    let o = (k div V)*(ne*V) + k mod V
    for c in 0..<ne div 2:
      let t = 2.0*PI*float(uniform(a))
      x[o + 2*c*V] = cos(t)
      x[o + (2*c+1)*V] = sin(t)
    rngStore(s, V, k, a)
