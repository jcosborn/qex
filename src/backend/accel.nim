import macros, strutils
import qex
import base/qexInternal

#const Backend {.strdefine.} = "OpenMP"
const Backend {.strdefine.} = "CUDA"
#const Backend {.strdefine.} = "SYCL"
#const Backend {.strdefine.} = "CPU"

when Backend == "OpenMP":
  const backendIsGpu* = true
  import openmp
  export openmp
elif Backend == "CUDA":
  const backendIsGpu* = true
  import cuda
  export cuda
  proc init = gpuInit(myRank)
  qexGlobalInitializers.add init
elif Backend == "SYCL":
  const backendIsGpu* = true
  #const backendIsGpu* = false
  import syclbe
  export syclbe
  proc init = gpuInit(myRank)
  qexGlobalInitializers.add init
else:
  when Backend != "CPU":
    static: echo "Backend: ", Backend
    {.warning: "Backend unknown, using CPU only.".}
  const backendIsGpu* = false
  import cpu
  export cpu
const backendIsCpu* = not backendIsGpu

var kernelCallCount* = 0  # count of kernel calls

proc gpuMalloc*[T](x: var ptr T) =
  let n = sizeof(T)
  x = cast[ptr T](gpuMalloc(n))

template toGpu*(x: typedesc): auto = false
template getGpu*(x: typedesc, g: auto): typedesc = x
template fromGpu*(x: typedesc, g: auto) = discard

template gpuType*[T:SomeNumber](x: typedesc[T]): typedesc = T
template toGpu*(x:SomeNumber):auto = x
#template toGpu*(x:var SomeNumber):auto = addr x
#template getGpu*(x:SomeNumber, g:SomeNumber):auto = g
macro getGpu*(x:SomeNumber, g:SomeNumber):auto =
  #echo x.treerepr
  #if x.kind == nnkSym and x.symKind == nskConst:
  if x.kind in nnkLiterals:
    result = x
  else:
    result = g
template getGpu*(x:SomeNumber, g:ptr SomeNumber):auto = g[]
template fromGpu*(x:SomeNumber, g:SomeNumber) = discard
template fromGpu*(x:var SomeNumber, g:ptr SomeNumber) = (x = g[])

import gpumem
export gpumem

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

template fromGpu*(x: seq, g: GpuSeq) =
  #when backendIsGpu:
  #  for i in 0..<x.len:
  #    x[i].fromGpu(g[i])
  x.copyFromGpu(g)

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
  if pgm.needsCopyIn:
    gpuMemCpyToGpu(g.p, addr c[0], c.bytes)

proc toGpu*[G,C](g: var GpuSeq[GpuSeq[G]], c: seq[seq[C]], pgm: ptr GpuMem) =
  if pgm.useCount == 1:  # newly created
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
  if pgm.useCount == 1:  # newly created
    var t = newSeq[ptr UncheckedArray[G]](g.n)
    for i in 0..<g.n:
      t[i].toGpu(c[i])
    gpuMemCpyToGpu(g.p, addr t[0], t.bytes)
  elif pgm.needsCopyIn:
    for i in 0..<g.n:
      #toGpu(addr g.p[i], c[i])
      pushGpuMemTag("SeqPtr" & displayName(G))
      let t = getGpuMem(addr c[i][0], c[i].bytes)
      popGpuMemTag()
      if t.useCount == 1:
        var x = toGpu(c[i])
        gpuMemCpyToGpu(addr g.p[i], x.p, sizeof(x.p))
      else:
        if t.needsCopyIn:
          gpuMemCpyToGpu(t.p, addr c[i][0], c[i].bytes)

proc toGpu*[T:SomeNumber](g: var ptr UncheckedArray[T], x: seq[T]) =
  mixin toGpu
  pushGpuMemTag("Ptr"&capitalizeAscii(displayName(T)))
  let pgm = getGpuMem(addr x[0], x.bytes)
  popGpuMemTag()
  g = cast[typeof g](pgm.p)
  g.toGpu(x, pgm)

proc toGpu*[T:SomeNumber](g: var ptr UncheckedArray[T], c: seq[T], pgm: ptr GpuMem) =
  if pgm.needsCopyIn:
    gpuMemCpyToGpu(g, addr c[0], c.bytes)


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
