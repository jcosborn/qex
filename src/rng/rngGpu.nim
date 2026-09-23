## The generators of a host RNG field on the GPU, stored like a field of a
## SIMD layout lo of V lanes: 32-bit word w of the generator of site k at
## (k div V)*W*V + w*V + k mod V, W words per generator.  Site k draws the
## numbers of the host generator at its coordinates, so the kernels below
## fill device fields as their host versions fill fields of lo.  Works for
## generators of plain data whose draws are inline procs, which then
## compile for the device: RngMilc6, MRG32k3a, Philox4x64, Threefry4x64.
import qex
import backend/accel

type RngGpu*[V: static int; R] = object
  lo*: Layout[V]
  s*: ptr UncheckedArray[uint32]  # [n div V][W][V]

template words(R: typedesc): int =
  static: doAssert sizeof(R) mod sizeof(uint32) == 0
  sizeof(R) div sizeof(uint32)

template load(s: ptr UncheckedArray[uint32]; V: static int; k: int; R: typedesc): untyped =
  block:
    const W = words(R)
    var a {.noInit.}: array[W, uint32]
    let o = (k div V)*(W*V) + k mod V
    for w in 0..<W: a[w] = s[o + w*V]
    cast[R](a)

template store(s: ptr UncheckedArray[uint32]; V: static int; k: int; x: typed) =
  block:
    const W = words(typeof(x))
    let a = cast[array[W, uint32]](x)
    let o = (k div V)*(W*V) + k mod V
    for w in 0..<W: s[o + w*V] = a[w]

proc hostIndex[V: static int; R](lo: Layout[V]; r: Field[1,R]): seq[int32] =
  ## the index in r of each site of lo
  result = newSeq[int32](lo.nSites)
  var c = newSeq[int32](lo.nDim)
  for k in 0..<lo.nSites:
    lo.coord(c, k)
    result[k] = int32 r.l.rankIndex(c).index

proc upload*[V: static int; R](g: RngGpu[V,R]; r: Field[1,R]) =
  ## the generators of r to g
  let n = g.lo.nSites
  let ix = g.lo.hostIndex(r)
  var h = newSeq[uint32](words(R)*n)
  let hp = cast[ptr UncheckedArray[uint32]](addr h[0])
  for k in 0..<n: store(hp, V, k, r[ix[k]])
  gpuMemCpyToGpu(g.s, hp, h.len*sizeof(uint32))

proc download*[V: static int; R](g: RngGpu[V,R]; r: Field[1,R]) =
  ## the generators of g to r
  let n = g.lo.nSites
  let ix = g.lo.hostIndex(r)
  var h = newSeq[uint32](words(R)*n)
  let hp = cast[ptr UncheckedArray[uint32]](addr h[0])
  gpuMemCpyToCpu(hp, g.s, h.len*sizeof(uint32))
  for k in 0..<n:
    let a = load(hp, V, k, R)
    copyMem(addr r[ix[k]], unsafeAddr a, sizeof(R))  # RNGs other than RngMilc6 lack :=

proc newRngGpu*[V: static int; R](lo: Layout[V]; r: Field[1,R]): RngGpu[V,R] =
  ## the generators of r on the device, for fields of lo
  result.lo = lo
  result.s = cast[ptr UncheckedArray[uint32]](gpuMalloc(words(R)*lo.nSites*sizeof(uint32)))
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
    var a = load(s, V, k, R)
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
    store(s, V, k, a)

proc gaussian*[V: static int; R](g: RngGpu[V,R]; x: ptr UncheckedArray[float]; ne: int) =
  ## x as gaussian of a field of ne reals per site in the SIMD layout, real
  ## c of site k at (k div V)*ne*V + c*V + k mod V
  let n = g.lo.nSites
  let s = g.s
  gpuFor(k, n):
    var a = load(s, V, k, R)
    let o = (k div V)*(ne*V) + k mod V
    for c in 0..<ne: x[o + c*V] = gaussian(a)
    store(s, V, k, a)
