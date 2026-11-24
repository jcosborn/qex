import macros, strutils, bitops
import accelbase
export accelbase
import gpumem
export gpumem
import expr
import qex

proc gpuFlagsIncl*(x: seq, f: set[gmFlags]) =
  ## The flags belong to the entry of x's data address: while they skip
  ## copies, keep x alive and in place, or free the entry (freeGpuMem) first.
  gpuMemFlagsIncl(addr x[0], f)
proc gpuFlagsExcl*(x: seq, f: set[gmFlags]) =
  ## as gpuFlagsIncl
  gpuMemFlagsExcl(addr x[0], f)

type
  GpuSeq*[T] = object
    n*: int
    p*: ptr UncheckedArray[T]
proc bytes*[T:GpuSeq](x: T): int = x.n * sizeof(T.T)
template gpuType*[T](x: typedesc[seq[T]]):typedesc = GpuSeq[gpuType(T)]
proc newGpuSeq*[T](g: var GpuSeq[T], n: int) =
  g.n = n
  g.p = cast[type g.p](gpuMalloc(g.bytes))
proc newGpuSeq*[T](n: int): GpuSeq[T] =
  result.newGpuSeq(n)
template `[]`*(x: GpuSeq, i: auto): auto =
  #doAssert(i>=0)
  #doAssert(i<x.n)
  x.p[i]

proc displayName*(x: typedesc[SomeNumber]): string =
  $x
proc displayName*[T](x: typedesc[ptr UncheckedArray[T]]): string =
  result = "Ptr" & capitalizeAscii($T)



template sharesMem(G, C: typedesc): bool =
  ## Values of C in host memory serve as their device copies of type G: no
  ## table entry, no copy.
  beSharedMem and G is C

proc toGpu*[T](x: seq[T]): auto =
  mixin gpuType, toGpu, displayName
  var g: GpuSeq[gpuType(T)]
  g.n = x.len
  when sharesMem(gpuType(T), T):
    g.p = cast[typeof g.p](addr x[0])
  else:
    let t = displayName(typeof(g[0]))
    pushGpuMemTag("Seq" & t)
    let pgm = getGpuMem(addr x[0], g.bytes)
    popGpuMemTag()
    g.p = cast[typeof g.p](pgm.p)
    g.toGpu(x, pgm)
  g

template getGpu*(x: seq, g: GpuSeq): auto = g

proc fromGpu*[C,G](x: var seq[C], g: GpuSeq[G]) =
  #when backendIsGpu:
  #  for i in 0..<x.len:
  #    x[i].fromGpu(g[i])
  #x.copyFromGpu(g)
  when not sharesMem(G, C):
    let pgm = getGpuMem(addr x[0], g.bytes)
    fromGpu(x, g, pgm)

proc toGpu*[G,C](g: var GpuSeq[G], x: seq[C]) =
  mixin toGpu, displayName
  g.n = x.len
  when sharesMem(G, C):
    g.p = cast[typeof g.p](addr x[0])
  else:
    let t = displayName(typeof(g[0]))
    pushGpuMemTag("Seq" & t)
    let pgm = getGpuMem(addr x[0], g.bytes)
    popGpuMemTag()
    g.p = cast[typeof g.p](pgm.p)
    g.toGpu(x, pgm)

proc toGpu*[T:SomeNumber](g: var GpuSeq[T], c: seq[T], pgm: ptr GpuMem) =
  #if pgm.needsCopyIn:
  #  gpuMemCpyToGpu(g.p, addr c[0], c.bytes)
  pgm.copyIn(addr c[0])

proc fromGpu*[T:SomeNumber](c: var seq[T], g: GpuSeq[T], pgm: ptr GpuMem) =
  pgm.copyOut(addr c[0])

proc toGpu*[G,C](g: var GpuSeq[ptr UncheckedArray[G]], c: seq[seq[C]], pgm: ptr GpuMem) =
  if pgm.needsCopyIn:
    var t = newSeq[ptr UncheckedArray[G]](g.n)
    for i in 0..<g.n:
      t[i].toGpu(c[i])
    pgm.copyIn(addr t[0])
  else:
    for i in 0..<g.n:
      discard toGpu(c[i])

proc fromGpu*[C,G](c: var seq[seq[C]], g: GpuSeq[ptr UncheckedArray[G]], pgm: ptr GpuMem) =
  if pgm.needsCopyOut:
    var t = newSeq[ptr UncheckedArray[G]](g.n)
    pgm.copyOut(addr t[0])
    for i in 0..<g.n:
      c[i].fromGpu(t[i])
  else:
    for i in 0..<g.n:
      fromGpu(c[i])

proc toGpu*[T:SomeNumber](g: var ptr UncheckedArray[T], x: seq[T]) =
  mixin toGpu
  when beSharedMem:
    g = cast[typeof g](addr x[0])
  else:
    pushGpuMemTag("Ptr"&capitalizeAscii(displayName(T)))
    let pgm = getGpuMem(addr x[0], x.bytes)
    popGpuMemTag()
    g = cast[typeof g](pgm.p)
    g.toGpu(x, pgm)

proc toGpu*[T:SomeNumber](g: var ptr UncheckedArray[T], c: seq[T], pgm: ptr GpuMem) =
  #if pgm.needsCopyIn:
  #  gpuMemCpyToGpu(g, addr c[0], c.bytes)
  pgm.copyIn(addr c[0])

proc fromGpu*[T:SomeNumber](c: var seq[T], g: ptr UncheckedArray[T]) =
  mixin fromGpu
  when not beSharedMem:
    let pgm = getGpuMem(addr c[0])
    pgm.copyOut(addr c[0])

proc fromGpu*[C](x: var seq[C]) =
  mixin gpuType
  when not sharesMem(gpuType(C), C):
    let pgm = getGpuMem(addr x[0])
    pgm.copyOut(addr x[0])


proc toGpu*[G,C](g: var GpuSeq[G], c: alignedMem[C]) =
  g.n = c.len
  when sharesMem(G, C):
    g.p = cast[typeof g.p](addr c[0])
  else:
    let pgm = getGpuMem(addr c[0], g.bytes)
    g.p = cast[typeof g.p](pgm.p)
    pgm.copyIn(addr c[0])

proc fromGpu*[C](c: var alignedMem[C]) =
  mixin gpuType
  when not sharesMem(gpuType(C), C):
    let pgm = getGpuMem(addr c[0])
    pgm.copyOut(addr c[0])

template fromGpu*(c: var alignedMem, g: GpuSeq) = fromGpu(c)


# A kernel body runs once in each of its threads.  The launch declares the
# number of the thread and their count for the text of the body, so procs
# cannot name them: code that takes its share of a loop is a template.

template gpuThreadNum*: untyped =
  ## in a kernel: the number of this thread, below gpuNumThreads
  gpuTn
template gpuNumThreads*: untyped =
  ## in a kernel: its threads, the launch count on a device, the host threads
  ## with the CPU backend
  gpuNt

iterator gpuRangeOf*(n, t, nt: int): int =
  ## the indices below n of thread t of nt: every nt-th on a device, which
  ## keeps neighboring threads on neighboring memory, a contiguous share on
  ## the host
  when backendIsGpu:
    var i = t
    while i < n:
      yield i
      i += nt
  else:
    for i in (n*t) div nt ..< (n*(t+1)) div nt:
      yield i

template gpuRange*(n: int): untyped =
  ## in a kernel: the share of this thread of 0..<n
  gpuRangeOf(n, int(gpuThreadNum), int(gpuNumThreads))

type
  SiteV*[V:static int] = distinct int
template `[]`*(x: SiteV): int = int(x)

when backendIsGpu:
  iterator gpuSitesOf*(n: int, V: static int, t, nt: int): SiteV[1] =
    for s in gpuRangeOf(n*V, t, nt):
      yield SiteV[1](s)
else:
  iterator gpuSitesOf*(n: int, V: static int, t, nt: int): SiteV[V] =
    for s in gpuRangeOf(n, t, nt):
      yield SiteV[V](s)

template gpuSites*(n: int, V: static int): untyped =
  ## in a kernel: the share of this thread of n outer sites of V lanes, one
  ## lane at a time on a device, whole sites on the host
  gpuSitesOf(n, V, int(gpuThreadNum), int(gpuNumThreads))

template gpuType*[T](t: typedesc[Simd[T]]): typedesc =
  ## the lanes as an array on a device; the host type with the CPU backend,
  ## whose kernels take whole SIMD sites
  when backendIsCpu: Simd[T]
  else: Simd[array[T.numNumbers,T.numberType]]
template gpuType*[T](t: typedesc[ComplexType[T]]): typedesc =
  ComplexType[gpuType(T)]
template gpuType*[N:static int; T](t: typedesc[VectorArray[N,T]]): typedesc =
  VectorArray[N,gpuType(T)]
template gpuType*[N,M:static int; T](t: typedesc[MatrixArray[N,M,T]]): typedesc =
  MatrixArray[N,M,gpuType(T)]
template gpuType*[T](t: typedesc[Color[T]]): typedesc =
  Color[gpuType(T)]
template gpuType*[V:static int, T](t: typedesc[Field[V,T]]): typedesc =
  GpuField[V,gpuType(T)]

template gpuSites*(lo: Layout): int = lo.nSites
template gpuSites*(n: SomeInteger): int =
  ## the thread count of a launch given as a number
  int(n)

proc fmadd*(a, b, c: float32): float32 {.importc: "fmaf", header: "<math.h>", noSideEffect.}
  ## a b + c rounded once
proc fmadd*(a, b, c: float64): float64 {.importc: "fma", header: "<math.h>", noSideEffect.}
  ## a b + c rounded once

template cmac*(re, im, ar, ai, br, bi: untyped) =
  ## re + i im += (ar + i ai)(br + i bi), one fmadd per real product; re and
  ## im must not alias ar, ai, br, bi (im reads them after re is updated)
  re = fmadd(ar, br, re)
  re = fmadd(-ai, bi, re)
  im = fmadd(ar, bi, im)
  im = fmadd(ai, br, im)

#import gpumem
#export gpumem

const sumTerms = 16  # values per thread of a level of sumFixed
const sumHost* = 512  # values per sum sumFixed leaves to the host

proc sumFixed*(r: var openArray[float]; a, w, h: ptr UncheckedArray[float]; m, n: int) =
  ## r[c] = the sum of a[c*n + i] over i < n for c < m, in an order fixed by
  ## n: each level adds the values t, t+T, ..., T = ceil(n/16), in thread t,
  ## until at most sumHost values remain, which the host adds.  w holds
  ## m*(n div 15 + 16) reals, the pinned h m*sumHost.  The levels and the
  ## copy to h queue after the kernels already submitted; one wait covers all.
  var src = a
  var cnt = n
  var off = 0
  while cnt > sumHost:
    let t1 = (cnt + sumTerms - 1) div sumTerms
    let s = src
    let d = cast[ptr UncheckedArray[float]](addr w[off])
    let c0 = cnt
    gpuForAsync(k, m*t1):
      let c = k div t1
      let t = k - c*t1
      var acc = 0.0
      for j in 0..<sumTerms:
        let i = t + j*t1
        if i < c0: acc += s[c*c0 + i]
      d[c*t1 + t] = acc
    src = d
    off += m*t1
    cnt = t1
  gpuMemCpyToCpuAsync(h, src, m*cnt*sizeof(float))
  gpuWaitAsync()
  for c in 0..<m:
    r[c] = 0.0
    for i in 0..<cnt: r[c] += h[c*cnt + i]

# Exact sums.  A kernel adds each term to an integer accumulator with atomic
# adds: integer adds commute, so the sum does not depend on the order of the
# device threads, and a run repeats bit for bit.
#   t = +-m 2^(e-1075), m < 2^53:  m 2^(e-1) = lo 2^(52 row) + hi 2^(52 (row+1)),
#   row = (e-1) div 52, lo and hi below 2^52
# A row is one uint64 per slot of weight 2^(52 row - 1074); term i adds lo and
# hi to rows row and row+1 of slot i mod s.  A word holds 2048 pieces between
# two reads, and newSumX gives a slot up to 1024 terms.  The words are never
# reset: a read takes their change since the last read, modulo 2^64, and
# rounds the exact sum once.

const
  sumRows = 41  # rows that cover the doubles
  sumWin* = 3  # rows a read copies: terms within 2^52 of each other fit

template sumXAdd*(a: ptr UncheckedArray[uint64]; o, st, w: int; t: float) =
  ## In kernels: adds t to rows row and row+1 of a SumX, at a[o + row*st],
  ## with the window w of the SumX.  A term outside the rows w, w+1 marks
  ## word w*st, so that the read knows to copy every row.  Zeros and
  ## subnormal terms add nothing, as where the processor flushes them; a non
  ## finite t adds a magnitude beyond the doubles.
  let tt = t
  let b = cast[uint64](tt)
  let e = int((b shr 52) and 0x7ff'u64)
  if e != 0:
    let m = (b and 0xF_FFFF_FFFF_FFFF'u64) or (1'u64 shl 52)
    let row = (e - 1) div 52
    let sh = e - 1 - 52*row
    var lo = (m shl sh) and 0xF_FFFF_FFFF_FFFF'u64
    var hi = m shr (52 - sh)
    if (b shr 63) != 0'u64 and e != 0x7ff:
      lo = 0'u64 - lo
      hi = 0'u64 - hi
    if cast[uint](row - w) > uint(sumWin - 2): a[w*st] = 1'u64
    let k = o + row*st
    gpuAtomicAdd(a, k, lo)
    gpuAtomicAdd(a, k + st, hi)

const sumDigits = (52*sumRows + 64 + 31) div 32  # base 2^32 digits of the row sums

proc sumXValue(d: openArray[int64]): float =
  ## The double nearest to sum_k d[k] 2^(32 k - 1074), ties to even;
  ## d.len <= sumDigits, |d[k]| < 2^62
  var g: array[sumDigits + 2, uint64]  # base 2^32 digits of the magnitude
  var neg = false
  var c = 0'i64
  for pass in 0..1:
    c = 0
    for k in 0..<d.len:
      let v = (if neg: -d[k] else: d[k]) + c
      g[k] = uint64(v and 0xFFFF_FFFF'i64)
      c = ashr(v, 32)
    if c >= 0: break
    neg = true  # a negative sum: the digits of its negative
  g[d.len] = uint64(c) and 0xFFFF_FFFF'u64
  g[d.len + 1] = uint64(c) shr 32
  var h = d.len + 1
  while h >= 0 and g[h] == 0: dec h
  if h < 0: return 0.0
  let bl = 64 - countLeadingZeroBits(g[h])
  var top = 32*h + bl - 1  # the leading bit, of weight 2^(top - 1074)
  var bits: uint64
  if top <= 52:  # below 2^53 units: exact, the bits of the double
    bits = g[0] or (g[1] shl 32)
  else:
    template get(i: int): uint64 = (if i >= 0: g[i] else: 0'u64)
    let m = (g[h] shl (64 - bl)) or (get(h-1) shl (32 - bl)) or (get(h-2) shr bl)  # the leading 64 bits
    var sticky = (get(h-2) and ((1'u64 shl bl) - 1)) != 0
    for i in 0..h-3:
      if g[i] != 0: sticky = true
    var sig = m shr 11
    let rem = m and 0x7FF'u64
    if rem > 0x400'u64 or (rem == 0x400'u64 and (sticky or (sig and 1) == 1)):
      inc sig
      if sig == 1'u64 shl 53:
        sig = 1'u64 shl 52
        inc top
    let e = top - 51
    bits = if e >= 2047: 0x7FF'u64 shl 52 else: (uint64(e) shl 52) or (sig and 0xF_FFFF_FFFF_FFFF'u64)
  result = cast[float](bits)
  if neg: result = -result

type SumX* = ref object
  ## m exact sums on the device: kernels add terms with sumXAdd, read gives
  ## the sums of the terms since the last read.
  m*, s*: int  # sums, slots per sum (a power of 2)
  st*: int  # words per row: the mark, s slots of each sum, padding to 64 bytes
  w*: int  # first row of the window
  n: int  # terms a slot got since the last read, at most
  a*: ptr UncheckedArray[uint64]  # [sumRows][st] on the device
  h: ptr UncheckedArray[uint64]  # pinned host copy
  p: seq[uint64]  # the words at the last read

proc newSumX*(m, n: int; slots = 0): SumX =
  ## m sums of up to n terms between two reads.  Slots: `slots`, or by
  ## default n/64 within 64 to 512, enough that the device threads of a hop
  ## kernel rarely wait for a word, and at most 256 where the copy of a read
  ## with more would exceed gpuCopyFast.  Always enough for 1024 terms per
  ## slot.  Rows start on 64 bytes: a copy from the device that starts
  ## elsewhere takes microseconds longer.
  ## Word of term i of sum c in row r: r*st + 1 + c*s + (i and (s-1))
  template rowWords(s: int): int = (1 + m*s + 7) and not 7
  result = SumX(m: m, s: 1)
  var want = if slots > 0: slots else: clamp(n div 64, 64, 512)
  if slots == 0:
    while want > 256 and sumWin*rowWords(want)*sizeof(uint64) > gpuCopyFast: want = want div 2
  while result.s < want or 1024*result.s < n: result.s *= 2
  result.st = rowWords(result.s)
  let nw = sumRows*result.st
  result.a = cast[ptr UncheckedArray[uint64]](gpuMalloc(nw*sizeof(uint64)))
  result.h = cast[ptr UncheckedArray[uint64]](gpuMallocHost(nw*sizeof(uint64)))
  result.p = newSeq[uint64](nw)
  let a = result.a
  gpuFor(i, nw): a[i] = 0'u64

proc free*(x: SumX) =
  ## the device and host words; a second call does nothing
  if x.a != nil:
    gpuFree(x.a)
    gpuFreeHost(x.h)
    x.a = nil

proc terms*(x: SumX; n: int) =
  ## A kernel is about to add n terms to each sum, those of n indices in a
  ## row: up to ceil(n/s) to a slot.  A word holds 2048 pieces below 2^52
  ## between two reads; with more it could wrap unseen.
  x.n += (n + x.s - 1) div x.s
  if x.n > 2048: qexError "SumX: up to ", x.n, " terms per slot between reads, a word holds 2048"

proc read*(x: SumX; r: var openArray[float]) =
  ## r[c] = sum c of the terms added since the last read, the exact sum
  ## rounded once.  Waits for the kernels.  Moves the window to the rows of
  ## these terms.
  x.n = 0
  let st = x.st
  let s = x.s
  let h = x.h
  let p = cast[ptr UncheckedArray[uint64]](addr x.p[0])
  var r0 = x.w
  var nr = sumWin
  gpuMemCpyToCpuAsync(addr h[r0*st], addr x.a[r0*st], nr*st*sizeof(uint64))  # after the kernels
  gpuWaitAsync()
  if h[r0*st] != 0:  # terms outside the window: every row counts
    let a = x.a
    let k = r0*st
    gpuFor(i, 1): a[k] = 0'u64
    r0 = 0
    nr = sumRows
    gpuMemCpyToCpu(h, x.a, nr*st*sizeof(uint64))
  var top = -1  # the highest row that changed
  for c in 0..<x.m:
    var d: array[sumDigits, int64]
    for row in r0..<r0+nr:
      let hh = cast[ptr UncheckedArray[uint64]](addr h[row*st + 1 + c*s])
      let pp = cast[ptr UncheckedArray[uint64]](addr p[row*st + 1 + c*s])
      # the change of the row, sum of v as int64 = hi 2^32 + lo - ng 2^64;
      # shifts and masks only, so that the loop vectorizes
      var lo, hi, ng, ch = 0'u64
      for k in 0..<s:
        let v = hh[k] - pp[k]
        ch = ch or v
        lo += v and 0xFFFF_FFFF'u64
        hi += v shr 32
        ng += v shr 63
      if ch != 0:
        top = max(top, row)
        # d += (hi 2^32 + lo - ng 2^64) 2^(52 row), each piece below 2^32 2^(bit mod 32)
        for (v, bit) in [(lo, 52*row), (hi, 52*row + 32)]:
          let q = bit div 32
          let u = (v and 0xFFFF_FFFF'u64) shl (bit mod 32)
          d[q] += int64(u and 0xFFFF_FFFF'u64)
          d[q+1] += int64(u shr 32) + int64((v shr 32) shl (bit mod 32))
        d[(52*row + 64) div 32] -= int64(ng shl ((52*row + 64) mod 32))
    r[c] = sumXValue(d)
  copyMem(addr p[r0*st], addr h[r0*st], nr*st*sizeof(uint64))
  if top >= 0: x.w = clamp(top - sumWin + 1, 0, sumRows - sumWin)

template gpuSumX*(x: SumX; i: untyped; n: SomeInteger; body: untyped) =
  ## Adds, for i in 0..<n, the values body gives as an array[x.m, float] to
  ## the sums of x.  Returns before the kernel completes; read waits for it.
  block:
    let sxN = int(n)
    x.terms(sxN)
    let sxA = x.a
    let sxS = x.s
    let sxSt = x.st
    let sxW = x.w
    gpuForAsync(i, sxN):
      let v = body
      for c in 0..<v.len: sumXAdd(sxA, 1 + c*sxS + (i and (sxS-1)), sxSt, sxW, v[c])

var gpuSumHost: ptr UncheckedArray[float]  # [8 sumHost] pinned host buffer of gpuSum
var gpuSumPart: ptr UncheckedArray[float]  # [m][T] sums of the threads, then the workspace of sumFixed
var gpuSumPartLen = 0

template gpuSum*(i: untyped; n: SomeInteger; m: static int; body: untyped): array[m, float] =
  ## The m <= 8 sums over i in 0..<n of the values body gives as an
  ## array[m, float], on this rank, in an order fixed by n: thread t sums
  ## i = t, t+T, ... of 16 terms, sumFixed the threads, so equal inputs give
  ## equal sums.
  static: doAssert m <= 8, "gpuSumHost holds 8 sums"
  block:
    if gpuSumHost == nil:
      gpuSumHost = cast[ptr UncheckedArray[float]](gpuMallocHost(8*sumHost*sizeof(float)))
    let gpuN = int(n)
    let gpuT = (gpuN + sumTerms - 1) div sumTerms
    let gpuL = m*gpuT + m*(gpuT div 15 + 16)
    if gpuSumPartLen < gpuL:
      if gpuSumPart != nil: gpuFree(gpuSumPart)
      gpuSumPartLen = max(gpuL, 2*gpuSumPartLen)
      gpuSumPart = cast[ptr UncheckedArray[float]](gpuMalloc(gpuSumPartLen*sizeof(float)))
    let sp = gpuSumPart
    gpuForAsync(t, gpuT):
      var a {.noInit.}: array[m, float]
      for c in 0..<m: a[c] = 0.0
      for j in 0..<sumTerms:
        let i = t + j*gpuT
        if i < gpuN:
          let v: array[m, float] = body
          for c in 0..<m: a[c] += v[c]
      for c in 0..<m: sp[c*gpuT + t] = a[c]
    var r: array[m, float]
    sumFixed(r, sp, cast[ptr UncheckedArray[float]](addr sp[m*gpuT]), gpuSumHost, m, gpuT)
    r

when backendIsGpu:
  type
    GpuSum*[T] = ref object
      ## The sum of the values the threads of a kernel hand to reduce, each
      ## once; T is float or an array of floats.  The device keeps the value
      ## of every thread and value adds them in an order fixed by the thread
      ## count (sumFixed).  A kernel that names a GpuSum launches outside
      ## threads blocks.
      n: int  # threads it has room for
      nt: int  # threads of the last kernel
      p: ptr UncheckedArray[float]  # device: [m][nt] values, then the workspace of sumFixed
      h: ptr UncheckedArray[float]  # pinned host block of sumFixed
      val: T
      valid: bool
    GpuSumView*[T] = object
      ## a GpuSum in a kernel
      p*: ptr UncheckedArray[float]

  template sumLen(T: typedesc): int =
    when T is array: len(T) else: 1

  proc newGpuSum*[T](n: int): GpuSum[T] =
    ## for kernels of up to n threads
    const m = sumLen(T)
    result = GpuSum[T](n: n)
    result.p = cast[ptr UncheckedArray[float]](gpuMalloc((m*n + m*(n div 15 + 16))*sizeof(float)))
    result.h = cast[ptr UncheckedArray[float]](gpuMallocHost(m*sumHost*sizeof(float)))

  proc free*(x: GpuSum) =
    ## the device and host blocks; a second call does nothing
    if x.p != nil:
      gpuFree(x.p)
      gpuFreeHost(x.h)
      x.p = nil

  proc toGpu*[T](x: GpuSum[T]): GpuSumView[T] =
    if gpuLaunchThreads > x.n:
      qexError "GpuSum: a kernel of ", gpuLaunchThreads, " threads, room for ", x.n
    x.nt = gpuLaunchThreads
    x.valid = false
    GpuSumView[T](p: x.p)
  template getGpu*(x: GpuSum, g: GpuSumView): auto = g
  template fromGpu*(x: GpuSum, g: GpuSumView) = discard
  template gpuResident*[T](x: typedesc[GpuSum[T]]): bool = true

  template reduce*[T](gs: GpuSumView[T], x: T) =
    ## in a kernel, once in every thread: the value of the thread
    let gpuSumV = x
    when T is array:
      for c in 0..<gpuSumV.len:
        gs.p[c*int(gpuNumThreads) + int(gpuThreadNum)] = float(gpuSumV[c])
    else:
      gs.p[int(gpuThreadNum)] = float(gpuSumV)

  proc value*[T](x: GpuSum[T]): T =
    ## the sum of the last kernel; waits for it
    if not x.valid:
      const m = sumLen(T)
      var r: array[m, float]
      sumFixed(r, x.p, cast[ptr UncheckedArray[float]](addr x.p[m*x.nt]), x.h, m, x.nt)
      when T is array:
        for c in 0..<m: x.val[c] = r[c]
      else:
        x.val = r[0]
      x.valid = true
    x.val

when isMainModule:
  #import qex
  #qexInit()
  proc test1 =
    var x = 1.0'f32
    #var yp = cast[ptr float32](gpuMalloc(sizeof(float32)))
    echo "x: ", x
    #threads:
    onGpu:
      x = 2.0
      #if getThreadNum()==0:
      #  printf("test\n")
    echo "x: ", x

  test1()
