## The RngMilc6 generators of a host RNG field on the GPU, stored like a
## field of a SIMD layout lo of V lanes: word w of the generator of site k
## at (k div V)*9V + w*V + k mod V, words r0..r6, icState, multiplier.
## Site k draws the numbers of the host generator at its coordinates, so
## the kernels below fill device fields as their host versions fill fields
## of lo; a kernel loads a generator, draws with the RngMilc6 procs, which
## are inline so that they compile for the device, and stores it back.
import qex
import backend/accel

type RngGpu*[V: static int] = object
  lo*: Layout[V]
  s*: ptr UncheckedArray[uint32]  # [n div V][9][V]

template wo(V, k, w: untyped): untyped = (k div V)*(9*V) + w*V + k mod V

template load(s: ptr UncheckedArray[uint32]; V: static int; k: int): RngMilc6 =
  RngMilc6(r0: s[wo(V,k,0)], r1: s[wo(V,k,1)], r2: s[wo(V,k,2)], r3: s[wo(V,k,3)],
           r4: s[wo(V,k,4)], r5: s[wo(V,k,5)], r6: s[wo(V,k,6)],
           icState: s[wo(V,k,7)], multiplier: s[wo(V,k,8)])

template store(s: ptr UncheckedArray[uint32]; V: static int; k: int; a: RngMilc6) =
  s[wo(V,k,0)] = a.r0
  s[wo(V,k,1)] = a.r1
  s[wo(V,k,2)] = a.r2
  s[wo(V,k,3)] = a.r3
  s[wo(V,k,4)] = a.r4
  s[wo(V,k,5)] = a.r5
  s[wo(V,k,6)] = a.r6
  s[wo(V,k,7)] = a.icState
  s[wo(V,k,8)] = a.multiplier

proc hostIndex[V: static int](lo: Layout[V]; r: Field[1,RngMilc6]): seq[int32] =
  ## the index in r of each site of lo
  result = newSeq[int32](lo.nSites)
  var c = newSeq[int32](lo.nDim)
  for k in 0..<lo.nSites:
    lo.coord(c, k)
    result[k] = int32 r.l.rankIndex(c).index

proc upload*[V: static int](g: RngGpu[V]; r: Field[1,RngMilc6]) =
  ## the generators of r to g
  let n = g.lo.nSites
  let ix = g.lo.hostIndex(r)
  var h = newSeq[uint32](9*n)
  let hp = cast[ptr UncheckedArray[uint32]](addr h[0])
  for k in 0..<n: store(hp, V, k, r[ix[k]])
  gpuMemCpyToGpu(g.s, hp, 9*n*sizeof(uint32))

proc download*[V: static int](g: RngGpu[V]; r: Field[1,RngMilc6]) =
  ## the generators of g to r
  let n = g.lo.nSites
  let ix = g.lo.hostIndex(r)
  var h = newSeq[uint32](9*n)
  let hp = cast[ptr UncheckedArray[uint32]](addr h[0])
  gpuMemCpyToCpu(hp, g.s, 9*n*sizeof(uint32))
  for k in 0..<n: r[ix[k]] = load(hp, V, k)

proc newRngGpu*[V: static int](lo: Layout[V]; r: Field[1,RngMilc6]): RngGpu[V] =
  ## the generators of r on the device, for fields of lo
  result.lo = lo
  result.s = cast[ptr UncheckedArray[uint32]](gpuMalloc(9*lo.nSites*sizeof(uint32)))
  result.upload r

proc free*[V: static int](g: var RngGpu[V]) =
  gpuFree(g.s)
  g.s = nil

proc randomTAH*[V: static int](g: RngGpu[V]; p: ptr UncheckedArray[float]) =
  ## p as randomTAH of 4 fields of 3x3 matrices, p in the layout of
  ## GpuGauge.u: real e = 6a+2b of M_ab (e+1 its imaginary part) of
  ## direction mu at 18*mu*n + (k div V)*18V + e*V + k mod V
  let n = g.lo.nSites
  let s = g.s
  gpuFor(k, n):
    const s2 = 0.70710678118654752440  # sqrt(1/2)
    const s3 = 0.57735026918962576450  # sqrt(1/3)
    var a = load(s, V, k)
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

proc gaussian*[V: static int](g: RngGpu[V]; x: ptr UncheckedArray[float]; ne: int) =
  ## x as gaussian of a field of ne reals per site in the SIMD layout, real
  ## c of site k at (k div V)*ne*V + c*V + k mod V
  let n = g.lo.nSites
  let s = g.s
  gpuFor(k, n):
    var a = load(s, V, k)
    let o = (k div V)*(ne*V) + k mod V
    for c in 0..<ne: x[o + c*V] = gaussian(a)
    store(s, V, k, a)
