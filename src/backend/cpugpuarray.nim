const Backend {.strdefine.} = "OpenMP"
when Backend == "CUDA":
  const useGPU = true
  import cuda
  export cuda
elif Backend == "OpenCL":
  const useGPU = true
  import opencl
  export opencl
elif Backend == "OpenMP":
  const useGPU = true
  import openmp
  export openmp
else:
  {.warning: "Backend unknown, use CPU only.".}
  const useGPU = false
when useGPU:
  import expr
  import gpuarray
  export gpuarray
else:
  template onGpu*(x:untyped) = threads: x
  template onGpu*(n,x:untyped) = threads: x
  template onGpu*(n,t,x:untyped) = threads: x
  template packVarsStmt(x,y:untyped) = discard

import coalesced
export coalesced
import macros
import linalg
export linalg
include system/ansi_c
import base/threading
export threading

#template onGpu*(x: untyped): untyped = x
#template onGpu*(a,b,x: untyped): untyped = x

import vectorized
export vectorized
# XXX quick and dirty hack to simd
# ONLY works with V=4, M=1, elementType=float32
# const VLEN {.intdefine.} = 4
# import qexLite/simd
# type
#   SVec = object
#     p*:ptr array[0,SimdS4] # 4=VLEN
# template `[]`*(x:SVec, i:int):auto = x.p[i]

when useGPU:
  type
    ArrayObj*[V,M:static[int],T] = object
      p*: Coalesced[V,M,T]
      # p*: SVec
      n*: int
      g*: GpuArrayObj[V,M,T]
      lastOnGpu*: bool
      unifiedMem*: bool
      mem:pointer ## Pointer to the allocated memory.
else:
  type
    ArrayObj*[V,M:static[int],T] = object
      p*: Coalesced[V,M,T]
      # p*: SVec
      n*: int
      mem:pointer ## Pointer to the allocated memory.
  # Array*[V,M:static[int],T] = ref ArrayObj[V,M,T]
  # Array*[V,M:static[int],T] = ArrayRef[V,M,T]
  # Arrays* = ArrayObj | ArrayRef
  # Arrays2* = ArrayObj | ArrayRef
  # Arrays3* = ArrayObj | ArrayRef

proc initArrayObj*(r: var ArrayObj, n: int) =
  const align = 64
  type T = r.T
  var p: ptr T
  when useGPU:
    when Backend == "CUDA":
      r.unifiedMem = true
      if r.unifiedMem:
        let err = cudaMallocManaged(cast[ptr pointer](addr p), n*sizeof(T))
        # Somehow == and != doesn't work as expected here??!
        if err:
          if cast[cint](err) == cast[cint](cudaErrorNotSupported):
            echo "WARNING: cudaMallocManaged not supported.  Fall back to non-unified memory."
            r.unifiedMem = false
          else:
            echo "ERROR: cudaMallocManaged ", n*sizeof(T)
            quit cast[cint](err)
    else:  # wait until OpenMP 5.0
      r.unifiedMem = false
    if not r.unifiedMem:
      # p = createSharedU(T, n)
      p = cast[ptr T](allocShared(n*sizeof(T)+align))
    r.g.n = 0          # Meaning uninitialized
    r.lastOnGpu = false
  else:
    p = cast[ptr T](allocShared(n*sizeof(T)+align))
  r.mem = p
  let x = cast[ByteAddress](p)
  const a1 = align - 1
  p = cast[ptr T](x + (a1-((x+a1) mod align)))
  # echo "ArrayObj.mem: ",cast[ByteAddress](r.mem)
  # echo "ArrayObj.p: ",cast[ByteAddress](p)
  r.n = n
  r.p.initCoalesced(p, n)
  # r.p.p = cast[type(r.p.p)](p)
# proc init[T](r: var ArrayRef[T], n: int) =
#   r.new
#   r[].init(n)

proc free*(r: var ArrayObj) =
  when useGPU:
    if r.unifiedMem:
      r.p.p.gpuFree
    else:
      # r.p.p.freeShared
      r.mem.deallocShared
      r.g.free # Same as `toGpu`, r.g is not passed to init with unifiedMem.
  else:
    r.mem.deallocShared
# proc free*[T](r: ArrayRef[T]) =
#   if r.unifiedMem:
#     discard r.p.p.cudaFree
#   else:
#     r.p.p.freeShared
#     r.g.free

proc newArrayObj*(V,M:static[int], n:int, T:typedesc): auto {.noinit.} =
  var z {.noinit.}: ArrayObj[V,M,T]
  z.initArrayObj(n)
  z
# proc newArrayObj*[T](n: int): ArrayObj[T] =
#   result.init(n)

# proc newArrayRef*[T](r: var ArrayRef[T], n: int) =
#   r.init(n)
# proc newArrayRef*[T](n: int): ArrayRef[T] =
#   result.init(n)

template getThreadNum*: untyped = threadNum
template getNumThreads*: untyped = numThreads
template offloadUseVar*(x:ArrayObj):bool = true
template offloadUsePtr*(x:ArrayObj):bool = true
template rungpuPrepareOffload*(x:ArrayObj):bool = true
template runcpuFinalizeOffload*(x:ArrayObj):bool = false
template offloadPtr*(x:var ArrayObj):untyped =
  x.toGpu
  x.g.p.p
template offloadVar*(x:ArrayObj,p:untyped):untyped = x.g

proc toGpu*(x: var ArrayObj) =
  when useGPU:
    # echo ">>> toGpu"
    if x.unifiedMem:
      if x.g.n==0:
        x.g.n = x.n
        x.g.p = x.p
        # x.g.p.initCoalesced(cast[ptr x.T](x.p.p),x.n)
    else:
      if not x.lastOnGpu:
        if x.g.n==0: x.g.initGpuArrayObj(x.n)
        let err = gpuMemCpyToGpu(x.g.p.p, x.p.p, x.n*sizeof(x.T))
        if err != 0:
          echo "gpuMemCpyToGpu: ", err
          quit cast[cint](err)
        x.lastOnGpu = true

proc toCpu*(x: var ArrayObj) =
  when useGPU:
    if (not x.unifiedMem) and x.lastOnGpu:
      threadSingle:
        let err = gpuMemCpyToCpu(x.p.p, x.g.p.p, x.n*sizeof(x.T))
        if err != 0:
          echo "gpuMemCpyToGpu: ", err
          quit cast[cint](err)
      threadSingle:
        x.lastOnGpu = false

template getGpuPtr*(x: var ArrayObj): untyped =
  when useGPU:
    x.toGpu
    x.g

type ArrayIndex* = SomeInteger or ShortVectorIndex
template indexArray*(x: ArrayObj, i: ArrayIndex): untyped =
  x.p[i]
#template `[]=`(x: ArrayObj, i: SomeInteger, y: untyped): untyped =
#  x.p[][i] = y

macro indexArray*(x: ArrayObj{call}, y: ArrayIndex): untyped =
  # proc cleanUp(n:NimNode):NimNode =
  #   if n.kind in {nnkStmtListExpr,nnkStmtList} and n.len == 1:
  #     result = n[0]
  #   else:
  #     result = n
  # echo ">>>>>> indexArray"
  # echo "call[", y.repr, "]"
  # echo x.treerepr
  #if siteLocalsField.contains($x[0]):
  result = newCall(ident($x[0]))
  for i in 1..<x.len:
    let xi = x[i]
    #result.add cleanUp( quote do:
    result.add( quote do:
      indexArray(`xi`,`y`) )
  #else:
  #  result = quote do:
  #    let tt = `x`
  #    tt.d[`y`]
  # echo result.treerepr
  # echo "<<<<<< indexArray"
  #echo result.repr

template `[]`*(x: ArrayObj, i: ArrayIndex): untyped = indexArray(x, i)
#template `[]=`(x: ArrayObj, i: SomeInteger, y: untyped): untyped =
#  x.p[][i] = y

# template `[]`*(x: ArrayRef, i: SomeInteger): untyped = indexArray(x, i)
#template `[]=`(x: ArrayRef, i: SomeInteger, y: untyped): untyped =
#  x.p[][i] = y

template veclen(x:ArrayObj):untyped = x.p.veclen
iterator vIndicesT*(x:ArrayObj):ShortVectorIndex =
  mixin getThreadNum, getNumThreads
  let
    tid = getThreadNum()
    nid = getNumThreads()
    n0 = 0
    n1 = x.veclen
    n = n1 - n0
    ti0 = n0 + ((tid*n) div nid)
    ti1 = n0 + (((tid+1)*n) div nid)
  # echo "ti0: ", ti0, "  ti1: ", ti1
  var i = ti0
  while i < ti1:
    yield i.ShortVectorIndex
    inc(i)
  # var i = tid
  # while i<x.veclen:
  #   yield i.ShortVectorIndex
  #   i += nid

# var threadNum* = 0
# var numThreads* = 1
template `:=`*(x: ArrayObj, y: ArrayObj) =
  #cprintf("t %i/%i  b %i/%i\n", getThreadIdx(), getThreadDim(), getBlockIdx(), getBlockDim())
  #let i = getBlockDim().x * getBlockIdx().x + getThreadIdx().x
  # let tid = getThreadNum()
  # let nid = getNumThreads()
  packVarsStmt((x,y), toCpu)
  # var i = tid
  # # while i<x.n:
  # # while i<x.n div VLEN:
  # while i<x.veclen:
  #   x[i.ShortVectorIndex] := y[i.ShortVectorIndex]
  #   i += nid
  for i in vIndicesT(x):
    x[i] := y[i]

template `:=`*(x: ArrayObj, y: SomeNumber) =
  #cprintf("t %i/%i  b %i/%i\n", getThreadIdx(), getThreadDim(), getBlockIdx(), getBlockDim())
  #let i = getBlockDim().x * getBlockIdx().x + getThreadIdx().x
  # let tid = getThreadNum()
  # let nid = getNumThreads()
  packVarsStmt(x, toCpu)
  # var i = tid
  # # while i<x.n:
  # # while i<x.n div VLEN:
  # while i<x.veclen:
  #   x[i.ShortVectorIndex] := y
  #   i += nid
  for i in vIndicesT(x):
    x[i] := y

template `+=`*(x: ArrayObj, y: SomeNumber) =
  #cprintf("t %i/%i  b %i/%i\n", getThreadIdx(), getThreadDim(), getBlockIdx(), getBlockDim())
  #let i = getBlockDim().x * getBlockIdx().x + getThreadIdx().x
  # let tid = getThreadNum()
  # let nid = getNumThreads()
  packVarsStmt((x,y), toCpu)
  # var i = tid
  # # while i<x.n:
  # # while i<x.n div VLEN:
  # while i<x.veclen:
  #   x[i.ShortVectorIndex] += y
  #   i += nid
  for i in vIndicesT(x):
    x[i] += y

#template `+=`*(x: ArrayObj, y: ArrayObj) =
template `+=`*[VX,VY,MX,MY:static[int],TX,TY](x: ArrayObj[VX,MX,TX], y: ArrayObj[VY,MY,TY]) =
  #cprintf("t %i/%i  b %i/%i\n", getThreadIdx(), getThreadDim(), getBlockIdx(), getBlockDim())
  #let i = getBlockDim().x * getBlockIdx().x + getThreadIdx().x
  # let tid = getThreadNum()
  # let nid = getNumThreads()
  packVarsStmt((x,y), toCpu)
  # var i = tid
  # # while i<x.n:
  # # while i<x.n div VLEN:
  # while i<x.veclen:
  #   x[i.ShortVectorIndex] += y[i.ShortVectorIndex]
  #   i += nid
  for i in vIndicesT(x):
    x[i] += y[i]

template `*=`*(x: ArrayObj, y: SomeNumber) =
  packVarsStmt(x, toCpu)
  for i in vIndicesT(x):
    x[i] *= y

template norm2*(x:ArrayObj): untyped =
  packVarsStmt(x, toCpu)
  var r:type(toDouble(norm2(x[0.ShortVectorIndex])))
  for i in vIndicesT(x):
    r += x[i].norm2
  var z = simdSum(r)
  threadSum(z)
  z

proc `+`*[VX,VY,MX,MY:static[int],TX,TY](x: ArrayObj[VX,MX,TX], y: ArrayObj[VY,MY,TY]): auto =
  var r: ArrayObj[x.V,x.M,type(x[0]+y[0])]
  # when x is ArrayObj:
  #   var r: ArrayObj[x.V,x.M,type(x[0]+y[0])]
  # else:
  #   var r: ArrayRef[x.V,x.M,type(x[0]+y[0])]
  echo "+\n"
  r
proc `*`*[VX,VY,MX,MY:static[int],TX,TY](x: ArrayObj[VX,MX,TX], y: ArrayObj[VY,MY,TY]): auto =
  var r: ArrayObj[x.V,x.M,type(x[0]*y[0])]
  # when x is ArrayObj:
  #   var r: ArrayObj[x.V,x.M,type(x[0]*y[0])]
  # else:
  #   var r: ArrayRef[x.V,x.M,type(x[0]*y[0])]
  echo "*\n"
  r

template newColorMatrixArray*(V,M:static[int], n:int): untyped =
  newArrayObj(V,M,n,Colmat[3,float32])
template newComplexArray*(V,M:static[int], n:int): untyped =
  newArrayObj(V,M,n,Complex[float32])
template newFloatArray*(V,M:static[int], n:int): untyped =
  newArrayObj(V,M,n,float32)

template vectorType*(t:typedesc[float32],v:typedesc):untyped = v
template vectorType*(t:typedesc[Complex[float32]],v:typedesc):untyped = Complex[v]
template vectorType*[C:static[int]](t:typedesc[Colmat[C,float32]],v:typedesc):untyped = Colmat[C,v]

template elementType*(t:typedesc[float32]):untyped = float32
template elementType*(t:typedesc[Complex[float32]]):untyped = float32
template elementType*[C:static[int]](t:typedesc[Colmat[C,float32]]):untyped = float32

proc printf*(frmt: cstring): cint {.
  importc: "printf", header: "<stdio.h>", varargs, discardable.}

when isMainModule:
  threads:
    if getThreadNum() == 0: echo "Threads ",getThreadNum(),"/",getNumThreads()
  var N = 128
  type T = float32
  const V = structsize(SVec) div sizeof(T)

  macro dump(n:typed):typed =
    echo n.repr
    n
  proc testfloat =
    echo "### float"
    var x = newArrayObj(V,1,N,T)
    var y = newArrayObj(V,1,N,T)
    var z = newArrayObj(V,1,N,T)
    threads:
      x := 1
      y := 2
      z := 3
      x += y * z
      #if (x.n-1) mod getNumThreads() == getThreadNum():
      if (x.veclen-1) mod getNumThreads() == getThreadNum():
        cprintf("thread %i/%i\n", getThreadNum(), getNumThreads())
        cprintf("x[%i]: %g\n", x.n-1, x[x.n-1][])
    dump:
      onGpu(1,32):
        x += y * z
        if (x.n-1) mod getNumThreads() == getThreadNum():
          cprintf("thread %lld/%lld\n", getThreadNum(), getNumThreads())
          cprintf("x[%i]: %g\n", x.n-1, x[x.n-1][])
    x.toCpu
    if (x.veclen-1) mod getNumThreads() == getThreadNum():
      cprintf("thread %i/%i\n", getThreadNum(), getNumThreads())
      cprintf("x[%i]: %g\n", x.n-1, x[x.n-1][])
    x.free
    y.free
    z.free
  testfloat()

  proc testcomplex =
    echo "### complex"
    var x = newArrayObj(V,2,N,Complex[T])
    var y = newArrayObj(V,2,N,Complex[T])
    var z = newArrayObj(V,2,N,Complex[T])
    threads:
      x := 1
      y := 2
      z := 3
      x += y * z
      if (x.n-1) mod getNumThreads() == getThreadNum():
        cprintf("thread %i/%i\n", getThreadNum(), getNumThreads())
        cprintf("x[%i]: %g\n", x.n-1, x[x.n-1][].re)

    onGpu:
      x += y * z
      x += 1
      if (x.n-1) mod getNumThreads() == getThreadNum():
        cprintf("thread %lld/%lld\n", getThreadNum(), getNumThreads())
        cprintf("x[%i]: %g\n", x.n-1, x[x.n-1][].re)

    threads:
      x += y * z
      if (x.n-1) mod getNumThreads() == getThreadNum():
        cprintf("thread %i/%i\n", getThreadNum(), getNumThreads())
        cprintf("x[%i]: %g\n", x.n-1, x[x.n-1][].re)
    x.free
    y.free
    z.free
  testcomplex()

  proc testcolmat =
    echo "### colmat"
    var x = newArrayObj(V,2,N,Colmat[3,T])
    var y = newArrayObj(V,2,N,Colmat[3,T])
    var z = newArrayObj(V,2,N,Colmat[3,T])
    threads:
      x := 1
      y := 2
      z := 3
      x += y * z
      if (x.n-1) mod getNumThreads() == getThreadNum():
        cprintf("thread %i/%i\n", getThreadNum(), getNumThreads())
        cprintf("x[%i][0,0]: %g\n", x.n-1, x[x.n-1][].d[0][0].re)

    onGpu(N):
      x += y * z
      if (x.n-1) mod getNumThreads() == getThreadNum():
        cprintf("thread %lld/%lld\n", getThreadNum(), getNumThreads())
        cprintf("x[%i][0,0]: %g\n", x.n-1, x[x.n-1][].d[0][0].re)

    threads:
      x += y * z
      if (x.n-1) mod getNumThreads() == getThreadNum():
        cprintf("thread %i/%i\n", getThreadNum(), getNumThreads())
        cprintf("x[%i][0,0]: %g\n", x.n-1, x[x.n-1][].d[0][0].re)
      var n = x.norm2
      threadSingle:
        cprintf("x.norm2: %g\n", n)
    x.free
    y.free
    z.free
  testcolmat()
