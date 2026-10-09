import macros, strutils
import accelbase
export accelbase
import gpumem
export gpumem
import qex

proc gpuFlagsIncl*(x: seq, f: set[gmFlags]) =
  gpuMemFlagsIncl(addr x[0], f)
proc gpuFlagsExcl*(x: seq, f: set[gmFlags]) =
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



proc toGpu*[T](x: seq[T]): auto =
  mixin gpuType, toGpu, displayName
  var g: GpuSeq[gpuType(T)]
  g.n = x.len
  let t = displayName(typeof(g[0]))
  pushGpuMemTag("Seq" & t)
  let pgm = getGpuMem(addr x[0], g.bytes)
  popGpuMemTag()
  g.p = cast[typeof g.p](pgm.p)
  g.toGpu(x, pgm)
  g

template getGpu*(x: seq, g: GpuSeq): auto = g

proc fromGpu*(x: var seq, g: GpuSeq) =
  #when backendIsGpu:
  #  for i in 0..<x.len:
  #    x[i].fromGpu(g[i])
  #x.copyFromGpu(g)
  let pgm = getGpuMem(addr x[0], g.bytes)
  fromGpu(x, g, pgm)

proc toGpu*[G,C](g: var GpuSeq[G], x: seq[C]) =
  mixin toGpu, displayName
  g.n = x.len
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

proc toGpu*[G,C](g: var GpuSeq[GpuSeq[G]], c: seq[seq[C]], pgm: ptr GpuMem) =
  if pgm.isNew:  # newly created
    var t = newSeq[GpuSeq[G]](g.n)
    for i in 0..<g.n:
      t[i] = toGpu(c[i])
    gpuMemCpyToGpu(g.p, addr t[0], t.bytes)
  elif pgm.needsCopyIn:
    for i in 0..<g.n:
      #toGpu(addr g.p[i], c[i])
      let t = getGpuMem(addr c[i][0], c[i].bytes)
      if t.useCount == 1:
        var x = toGpu(c[i])
        gpuMemCpyToGpu(addr g.p[i], addr x, sizeof(x))
      else:
        if t.needsCopyIn:
          gpuMemCpyToGpu(t.p, addr c[i][0], c[i].bytes)

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
  let pgm = getGpuMem(addr c[0])
  pgm.copyOut(addr c[0])

proc fromGpu*(x: var seq) =
  let pgm = getGpuMem(addr x[0])
  pgm.copyOut(addr x[0])


proc toGpu*(g: var GpuSeq, c: alignedMem) =
  g.n = c.len
  let pgm = getGpuMem(addr c[0], g.bytes)
  g.p = cast[typeof g.p](pgm.p)
  pgm.copyIn(addr c[0])

proc fromGpu*(c: var alignedMem) =
  let pgm = getGpuMem(addr c[0])
  pgm.copyOut(addr c[0])

template fromGpu*(c: var alignedMem, g: GpuSeq) = fromGpu(c)


iterator gpuRange*(n: int): int =
  when backendIsGpu:
    let s = int gpuNumThreads()
    var i = int gpuThreadNum()
    while i < n:
      yield i
      i += s
  else:
    let s = gpuNumThreads()
    let id = gpuThreadNum()
    let i0 = (n*id) div s
    let i1 = (n*(id+1)) div s
    for i in i0 ..< i1:
      yield i

type
  SiteV*[V:static int] = distinct int
template `[]`*(x: SiteV): int = int(x)

when backendIsGpu:
  iterator gpuSites*(n:int, V:static int): SiteV[1] =
    for s in gpuRange(n*V):
      yield SiteV[1](s)
else:
  iterator gpuSites*(n:int, V:static int): SiteV[V] =
    for s in gpuRange(n):
      yield SiteV[V](s)

template gpuType*[T](t: typedesc[Simd[T]]): typedesc =
  Simd[array[T.numNumbers,T.numberType]]
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

#import gpumem
#export gpumem

const sumTerms = 16  # values per thread of a level of sumFixed
const sumHost* = 512  # values per sum sumFixed leaves to the host

proc sumFixed*(r: var openArray[float]; a, w, h: ptr UncheckedArray[float]; m, n: int) =
  ## r[c] = the sum of a[c*n + i] over i < n for c < m, in an order fixed by
  ## n: each level adds the values t, t+T, ..., T = ceil(n/16), in thread t,
  ## until at most sumHost values remain, which the host adds.  w holds
  ## m*(n div 15 + 16) reals, the pinned h m*sumHost.  The levels queue
  ## after the kernels already submitted; the copy to h waits for all.
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
  gpuWaitAsync()
  gpuMemCpyToCpu(h, src, m*cnt*sizeof(float))
  for c in 0..<m:
    r[c] = 0.0
    for i in 0..<cnt: r[c] += h[c*cnt + i]

var gpuSumHost: ptr UncheckedArray[float]  # [8 sumHost] pinned host buffer of gpuSum
var gpuSumPart: ptr UncheckedArray[float]  # [m][T] sums of the threads, then the workspace of sumFixed
var gpuSumPartLen = 0

template gpuSum*(i: untyped; n: SomeInteger; m: static int; body: untyped): array[m, float] =
  ## The m <= 8 sums over i in 0..<n of the values body gives as an
  ## array[m, float], on this rank, in an order fixed by n: thread t sums
  ## i = t, t+T, ... of 16 terms, sumFixed the threads, so equal inputs give
  ## equal sums.
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
