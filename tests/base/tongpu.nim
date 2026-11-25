import testutils
import qex
import backend/[accel, cgfield]
import base/metaUtils

# What a kernel body may name and how onGpu captures it (backend/expr).  The
# kernels launch inside procs: for code at the top level of a module Nim
# emits the pragmas of a device launch at file scope.

getOptimPragmas()
qexInit()

type
  Pt = object
    u, v: int
  Sc = object  # plain data without a toGpu: the launch copies it
    a, b: float
  Reals = ptr UncheckedArray[float]

const K = 3
var glob = 5.0  # copied into the launch

proc axpy(a, x, y: float; scale = 1.0): float = scale*(a*x + y)  # no mark: copied into device kernels
proc twice(x: float): float {.gpuInline.} = 2.0*x  # called

proc dev(n: int): Reals = cast[Reals](gpuMalloc(n*sizeof(float)))

proc host(p: Reals; n: int): seq[float] =
  ## the n reals at device address p
  result = newSeq[float](n)
  gpuMemCpyToCpu(addr result[0], p, n*sizeof(float))

proc names(n: int): int =
  ## b[i] = 2 a[i] + glob through a global, a constant, a named argument, a
  ## constructor and a field, with variables of the names of the argument
  ## and of the field: the elements of b that differ
  var a = newSeq[float](n)
  var b = newSeq[float](n)
  let scale = 3.0
  let u = 7
  let sc = Sc(a: 1.5, b: 0.5)
  for i in 0..<n: a[i] = float(i)
  onGpu(n):
    for i in gpuRange(n):
      let q = Pt(u: i, v: K)
      b[i] = axpy(1.0, twice(a[i]), glob, scale = 1.0) + float(q.v - K + q.u - i + u - 7) + scale - 3.0 + sc.a - 3.0*sc.b
  for i in 0..<n:
    if b[i] != 2.0*float(i) + 5.0: inc result

proc generic[V: static int; T](x: seq[T]; y: var seq[T]; s: T): int =
  ## y[i] = s x[i] + 0 + 1 + ... + (V-1) in a generic proc, the sum over the
  ## index of a macro: the elements of y that differ
  let n = x.len
  onGpu(n):
    for i in gpuRange(n):
      var t = T(0)
      forStatic c, 0, V-1:
        t += T(c)
      y[i] = s*x[i] + t
  for i in 0..<n:
    if y[i] != s*x[i] + T((V*(V-1)) div 2): inc result

proc nested[V: static int; T](y: var seq[T]): int =
  ## A kernel in a template of a generic proc, with a template of its own
  ## that loops with forStatic, as the hops of stagGpu: the elements of y
  ## that are not 0 + 1 + ... + (V-1)
  let n = y.len
  template kern(K: static int) =
    onGpu(n):
      template tot(r: untyped) =
        forStatic c, 0, K-1: r += T(c)
      for i in gpuRange(n):
        var t = T(0)
        tot(t)
        y[i] = t
  kern(V)
  for i in 0..<n:
    if y[i] != T((V*(V-1)) div 2): inc result

proc wide(n: int): int =
  ## a loop over twice the indices of the launch count: the elements of p
  ## that are not their index
  let p = dev(2*n)
  onGpu(n):
    for i in gpuRange(2*n): p[i] = float(i)
  let h = host(p, 2*n)
  for i in 0..<2*n:
    if h[i] != float(i): inc result
  gpuFree(p)

proc cases(n: int): int =
  ## A kernel that names a constant W and a seq w, whose names Nim tells
  ## apart by the case of the first letter only: the elements of w that are
  ## not 3 i
  const W = 3
  var w = newSeq[float](n)
  onGpu(n):
    for i in gpuRange(n): w[i] = float(W*i)
  for i in 0..<n:
    if w[i] != float(3*i): inc result

template fill(f: untyped; n: int; v: float) =
  ## a kernel in a template, on the caller's variable
  onGpu(n):
    for i in gpuRange(n): f[i] = v

proc inTemplate(n: int): int =
  var z = newSeq[float](n)
  fill(z, n, 4.0)
  for i in 0..<n:
    if z[i] != 4.0: inc result

proc outside(n: int): int =
  ## Device pointers named by a template outside the body, one of them by
  ## the body as well: the elements of p that are not 3.
  let p = dev(n)
  let q = dev(n)
  template pq(i: untyped): untyped = p[i] + q[i]
  onGpu(n):
    for i in gpuRange(n):
      p[i] = 1.0
      q[i] = 2.0
  onGpu(n):
    for i in gpuRange(n): p[i] = pq(i)
  let h = host(p, n)
  for i in 0..<n:
    if h[i] != 3.0: inc result
  gpuFree(p)
  gpuFree(q)

proc owned(n: int): int =
  ## A name the body declares is its own: the host seq x is no part of the
  ## kernel y[i] = x with x = i.
  var x = newSeq[float](n)
  var y = newSeq[float](n)
  for i in 0..<n: x[i] = -1.0
  onGpu(n):
    for i in gpuRange(n):
      let x = float(i)
      y[i] = x
  for i in 0..<n:
    if y[i] != float(i) or x[i] != -1.0: inc result

proc inThreads(n: int): int =
  ## every thread of a threads block calls onGpu on a device array
  let p = dev(n)
  threads:
    onGpu(n):
      for i in gpuRange(n): p[i] = float(i)
  let h = host(p, n)
  for i in 0..<n:
    if h[i] != float(i): inc result
  gpuFree(p)

proc fields(lo: Layout): float =
  ## |b|^2, over all ranks, after b[s] = 2 a[s] on the sites of a field of
  ## ones
  var a = lo.ColorVector()
  var b = lo.ColorVector()
  threads:
    a := 1.0
    b := 0.0
  onGpu(lo):
    for s in gpuSites(b):
      b[s] := 2.0*a[s]
  var r: float
  threads:
    let t = b.norm2
    threadSingle: r = t
  r

proc fieldCode(lo: Layout): float =
  ## |r - s|^2/|s|^2 of the same field code in a kernel (r) and in a threads
  ## block (s): statements on fields and subsets, loops over a field and over
  ## a subset
  var rng = lo.newRNGField(RngMilc6, 4321)
  var a = lo.ColorVector()
  var b = lo.ColorVector()
  var m = lo.ColorMatrix()
  var r = lo.ColorVector()
  var s = lo.ColorVector()
  a.gaussian rng
  b.gaussian rng
  m.gaussian rng
  onGpu(lo):
    r := 2.0*a + b
    r += m*a
    r -= m.adj*b
    r.odd := 0.5*r
    for e in r.even:
      r[e] := r[e] - a[e]
    for e in r:
      r[e] *= 3.0
  threads:
    s := 2.0*a + b
    s += m*a
    s -= m.adj*b
    threadBarrier()
    s.odd := 0.5*s
    threadBarrier()
    for e in s.even:
      s[e] := s[e] - a[e]
    threadBarrier()
    for e in s:
      s[e] *= 3.0
  var d, n: float
  threads:
    let x = norm2(r - s)
    let y = s.norm2
    threadSingle:
      d = x
      n = y
  d/n

proc fieldView(lo: Layout): float =
  ## |r - 2 a|^2/|a|^2 for a field view of device memory: a copied in, the
  ## view doubled, r copied out
  var rng = lo.newRNGField(RngMilc6, 1234)
  var a = lo.ColorVector()
  var r = lo.ColorVector()
  a.gaussian rng
  let bytes = lo.nSitesOuter*sizeof(a[0])
  let p = gpuMalloc(bytes)
  gpuMemCpyToGpu(p, addr a[0], bytes)
  let v = gpuFieldView(lo, typeof(a[0]), p)
  onGpu(lo):
    v := 2.0*v
  gpuMemCpyToCpu(addr r[0], p, bytes)
  gpuFree(p)
  var d, n: float
  threads:
    let x = norm2(r - 2.0*a)
    let y = a.norm2
    threadSingle:
      d = x
      n = y
  d/n

proc scaleBy(x: auto; a: float) =
  ## a share of a loop in a proc that a kernel calls; inlineProcs copies it
  ## into device kernels
  for e in x: x[e] *= a

proc fieldParts(lo: Layout): float =
  ## |r - s|^2/|s|^2 and the same for the fields of two seqs, kernel (r, fr)
  ## against threads block (s, fs): zero, negation, a proc that takes a
  ## share, loops over the fields of a seq
  var rng = lo.newRNGField(RngMilc6, 5678)
  var a = lo.ColorVector()
  var b = lo.ColorVector()
  a.gaussian rng
  b.gaussian rng
  var r = lo.ColorVector()
  var s = lo.ColorVector()
  var fr = @[lo.ColorVector(), lo.ColorVector()]
  var fs = @[lo.ColorVector(), lo.ColorVector()]
  onGpu(lo):
    r := 0
    r -= a
    r += -b
    scaleBy(r, 3.0)
    for k in 0..1:
      for e in fr[k]: fr[k][e] := float(k+1)*a[e]
  threads:
    s := 0
    s -= a
    s += -b
    for e in s: s[e] *= 3.0
    for k in 0..1:
      for e in fs[k]: fs[k][e] := float(k+1)*a[e]
  var d, n: float
  threads:
    var x = norm2(r - s)
    var y = s.norm2
    for k in 0..1:
      x += norm2(fr[k] - fs[k])
      y += fs[k].norm2
    threadSingle:
      d = x
      n = y
  d/n

proc total(n: int): float =
  ## the sum of 0..n-1 through a GpuSum
  var a = newSeq[float](n)
  for i in 0..<n: a[i] = float(i)
  var gs = newGpuSum[float](n)
  onGpu(n):
    var r = 0.0
    for i in gpuRange(n): r += a[i]
    gs.reduce r
  gs.value

proc queue(n: int): int =
  ## three queued kernels on a device array and one wait, the second
  ## reading an array (with OpenMP offload a kernel that completes at its
  ## launch, after the queue): the elements that are not 2 i + 1
  let p = dev(n)
  let c = [2.0, 1.0]
  onGpuQueue(n):
    for i in gpuRange(n): p[i] = float(i)
  onGpuQueue(n):
    for i in gpuRange(n): p[i] = c[0]*p[i]
  onGpuQueue(n):
    for i in gpuRange(n): p[i] = p[i] + c[1]
  gpuWaitAsync()
  let h = host(p, n)
  for i in 0..<n:
    if h[i] != 2.0*float(i) + 1.0: inc result
  gpuFree(p)

proc queueInThreads(n: int): int =
  ## every thread of a threads block queues two kernels, the second reading
  ## what the other threads wrote in the first: the elements that are not
  ## 2 (n - 1 - i)
  let p = dev(n)
  let q = dev(n)
  threads:
    onGpuQueue(n):
      for i in gpuRange(n): p[i] = float(i)
    onGpuQueue(n):
      for i in gpuRange(n): q[i] = 2.0*p[n-1-i]
    gpuWaitAsync()
  let h = host(q, n)
  for i in 0..<n:
    if h[i] != 2.0*float(n-1-i): inc result
  gpuFree(p)
  gpuFree(q)

proc queueSeq(): bool =
  ## whether a queued kernel on a seq, which is copied back, compiles
  var a = newSeq[float](4)
  compiles((block:
    onGpuQueue(4):
      for i in gpuRange(4): a[i] = 1.0))

proc noDevice(s: string): bool =
  ## whether a kernel on a string, which has no device form, compiles
  compiles((block:
    onGpu(4):
      for i in gpuRange(4): discard s.len))

let lo = latticeFromLocalLattice([4,4,4,8], nRanks).newLayout

suite "onGpu captures":
  test "globals, constants, named arguments, constructors, plain data":
    check names(1000) == 0

  test "a kernel in a generic proc with a static parameter":
    var x = newSeq[float32](100)
    var y = newSeq[float32](100)
    for i in 0..<x.len: x[i] = float32(i)
    check generic[4, float32](x, y, 3'f32) == 0
    var xd = newSeq[float](100)
    var yd = newSeq[float](100)
    for i in 0..<xd.len: xd[i] = float(i)
    check generic[2, float](xd, yd, 0.5) == 0

  test "a body template with forStatic in a template of a generic proc":
    var y = newSeq[float](100)
    check nested[5, float](y) == 0

  test "a loop over more indices than threads":
    check wide(1000) == 0

  test "names that differ in the case of the first letter":
    check cases(1000) == 0

  test "a kernel in a template":
    check inTemplate(1000) == 0

  test "device pointers of a template outside the body":
    check outside(1000) == 0

  test "a name the body declares":
    check owned(1000) == 0

  test "a kernel in a threads block":
    check inThreads(1000) == 0

  test "fields":
    check fields(lo) == 4.0*3.0*float(lo.physVol)

  test "field code as in a threads block":
    let d = fieldCode(lo)
    when backendIsGpu: check d <= 1e-28
    else: check d == 0.0

  test "field code in a proc and over the fields of a seq":
    let d = fieldParts(lo)
    when backendIsGpu: check d <= 1e-28
    else: check d == 0.0

  test "a field view of device memory":
    check fieldView(lo) == 0.0

  test "a sum":
    check total(1000) == 499500.0

  test "a host variable without a device form stops the compile":
    check not noDevice("abc")

  test "queued kernels":
    check queue(1000) == 0

  test "queued kernels in a threads block":
    var bad = 0
    for r in 0..<100: bad += queueInThreads(1000)
    check bad == 0

  test "a queued kernel on a variable that is copied back stops the compile":
    check not queueSeq()

qexFinalize()
