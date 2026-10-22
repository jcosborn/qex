## Times the draw kernels of rng/rngGpu with each thread's generator loaded
## and stored in chunks of c 32-bit words, objects aligned to 4c bytes: word
## w of site k at (k div V)*M*c*V + (w div c)*c*V + (k mod V)*c + w mod c,
## M = ceil(W/c) chunks of the W words (padding when c does not divide W).
## c = 1 is the layout of rngGpu.  Kernels over the n sites of the layout:
##   copy: the loads and stores of gaussian, no draws, the generators unchanged
##   uniform: two uniforms per real, the generator work of gaussian
##   gaussian: ne reals per site, rngGpu.gaussian
##   randomTAH: 72 reals per site, rngGpu.randomTAH
## For each: the best time per call over reps batches of 4 calls, and the
## bytes moved per second, also as a fraction of -peak (GB/s).  The gaussian
## draws of each c are checked against those of rngGpu.
##   -lat  -reps  -ne  -peak
import qex
import backend/accel, rng/rngGpu
import base/metaUtils
import times, strformat

type
  K1 = object
    x: array[1, uint32]
  K2 = object
    x {.align(8).}: array[2, uint32]
  K4 = object
    x {.align(16).}: array[4, uint32]

template kT(c: static int): typedesc =
  when c == 1: K1
  elif c == 2: K2
  else: K4

template chunks(R: typedesc; c: static int): int = (rngWords(R) + c - 1) div c

template ldc(s: ptr UncheckedArray[uint32]; V, c: static int; k: int; R: typedesc): untyped =
  ## the generator of site k in s, in chunks of c words
  block:
    const W = rngWords(R)
    const M = chunks(R, c)
    var a {.noInit.}: array[W, uint32]
    let p = cast[ptr UncheckedArray[kT(c)]](s)
    let o = (k div V)*(M*V) + k mod V
    forStatic m, 0, M-1:
      let y = p[o + m*V]
      forStatic j, 0, c-1:
        when m*c + j < W: a[m*c + j] = y.x[j]
    cast[R](a)

template stc(s: ptr UncheckedArray[uint32]; V, c: static int; k: int; g: typed) =
  ## g to the generator of site k in s, in chunks of c words
  block:
    const W = rngWords(typeof(g))
    const M = chunks(typeof(g), c)
    let a = cast[array[W, uint32]](g)
    let p = cast[ptr UncheckedArray[kT(c)]](s)
    let o = (k div V)*(M*V) + k mod V
    forStatic m, 0, M-1:
      var y {.noInit.}: kT(c)
      forStatic j, 0, c-1:
        when m*c + j < W: y.x[j] = a[m*c + j]
        else: y.x[j] = 0
      p[o + m*V] = y

proc copyK[V, c: static int; R](s: ptr UncheckedArray[uint32]; x: ptr UncheckedArray[float]; n, ne: int; z: uint32) =
  ## z = 0 at run time: the generators are stored unchanged
  gpuForAsync(k, n):
    const W = rngWords(R)
    let g = ldc(s, V, c, k, R)
    var b = cast[array[W, uint32]](g)
    var t = 0'u32
    forStatic w, 0, W-1:
      t = t xor b[w]
      b[w] = b[w] xor z
    let o = (k div V)*(ne*V) + k mod V
    for i in 0..<ne: x[o + i*V] = float(t + uint32(i))
    stc(s, V, c, k, cast[R](b))

proc unifK[V, c: static int; R](s: ptr UncheckedArray[uint32]; x: ptr UncheckedArray[float]; n, ne: int) =
  gpuForAsync(k, n):
    var a = ldc(s, V, c, k, R)
    let o = (k div V)*(ne*V) + k mod V
    for i in 0..<ne:
      let u = float(uniform(a))
      x[o + i*V] = u + float(uniform(a))
    stc(s, V, c, k, a)

proc gaussK[V, c: static int; R](s: ptr UncheckedArray[uint32]; x: ptr UncheckedArray[float]; n, ne: int) =
  gpuForAsync(k, n):
    var a = ldc(s, V, c, k, R)
    let o = (k div V)*(ne*V) + k mod V
    for i in 0..<ne: x[o + i*V] = gaussian(a)
    stc(s, V, c, k, a)

proc tahK[V, c: static int; R](s: ptr UncheckedArray[uint32]; p: ptr UncheckedArray[float]; n: int) =
  gpuForAsync(k, n):
    const s2 = 0.70710678118654752440  # sqrt(1/2)
    const s3 = 0.57735026918962576450  # sqrt(1/3)
    var a = ldc(s, V, c, k, R)
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
    stc(s, V, c, k, a)

proc report(name, label: string; t: float; n: int; bytes, peak: float) =
  let gbs = 1e-9*bytes*float(n)/t
  echo &"{name:<13} {label:<14} {1e6*t:9.1f} us {1e9*t/float(n):7.3f} ns/site {gbs:7.0f} GB/s {100*gbs/peak:5.1f}%"

proc run[V: static int](lo: Layout[V]; R: typedesc; name: string; reps, ne: int; peak: float) =
  let n = lo.nSites
  let rv = lo.newRNGFieldV(R, 987654321)
  let x = cast[ptr UncheckedArray[float]](gpuMalloc(max(ne, 72)*n*sizeof(float)))
  template bench(label: string; bytes: float; call: untyped) =
    var best = 0.0
    for r in 0..<reps:
      let t0 = epochTime()
      for b in 0..<4: call
      gpuWaitAsync()
      let dt = (epochTime() - t0)/4
      if r == 0 or dt < best: best = dt
    report(name, label, best, n, bytes, peak)
  var y0, y = newSeq[float](ne*n)
  block:  # the draws of rngGpu
    var g = newRngGpu(rv)
    g.gaussian(x, ne)
    gpuMemCpyToCpu(addr y0[0], x, y0.len*sizeof(float))
    let sb = float(2*4*rngWords(R))
    bench("rngGpu gauss", sb + float(8*ne)): g.gaussian(x, ne)
    bench("rngGpu TAH", sb + float(8*72)): g.randomTAH x
    g.free
  template variant(c: static int) =
    const M = chunks(R, c)
    var h = newSeq[uint32](M*c*n)
    let hp = cast[ptr UncheckedArray[uint32]](addr h[0])
    for k in 0..<n: stc(hp, V, c, k, rv[k])
    let s = cast[ptr UncheckedArray[uint32]](gpuMalloc(h.len*sizeof(uint32)))
    gpuMemCpyToGpu(s, hp, h.len*sizeof(uint32))
    gaussK[V, c, R](s, x, n, ne)
    gpuWaitAsync()
    gpuMemCpyToCpu(addr y[0], x, y.len*sizeof(float))
    echo name, " c", c, ": ", 4*M*c, " bytes per generator, gaussian draws ", (if y == y0: "same as rngGpu" else: "DIFFER from rngGpu")
    let sb = float(2*4*M*c)
    bench("c" & $c & " copy", sb + float(8*ne)): copyK[V, c, R](s, x, n, ne, zero)
    bench("c" & $c & " uniform", sb + float(8*ne)): unifK[V, c, R](s, x, n, ne)
    bench("c" & $c & " gaussian", sb + float(8*ne)): gaussK[V, c, R](s, x, n, ne)
    bench("c" & $c & " randomTAH", sb + float(8*72)): tahK[V, c, R](s, x, n)
    gpuFree(s)
  variant(1)
  variant(2)
  variant(4)
  gpuFree(x)

qexInit()
let lat = intSeqParam("lat", @[48,48,48,48])
let reps = intParam("reps", 10)
let ne = intParam("ne", 6)
let peak = floatParam("peak", 3350)
let zero = uint32(intParam("zero", 0))  # unknown to the compiler
let lo = lat.newLayout
echo "lattice ", lat, "  lanes ", lo.innerGeom, "  sites ", lo.nSites, "  ne ", ne, "  peak ", peak, " GB/s"
run(lo, RngMilc6, "RngMilc6", reps, ne, peak)
run(lo, MRG32k3a, "MRG32k3a", reps, ne, peak)
run(lo, Philox4x64, "Philox4x64", reps, ne, peak)
run(lo, Threefry4x64, "Threefry4x64", reps, ne, peak)
qexFinalize()
