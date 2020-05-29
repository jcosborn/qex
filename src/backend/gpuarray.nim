import coalesced

const Backend {.strdefine.} = "OpenMP"

when Backend == "CUDA":
  import cuda
elif Backend == "OpenCL":
  import opencl
elif Backend == "OpenMP":
  import openmp
else:
  {.warning: "Backend unknown, use OpenMP by default.".}
  import openmp

import macros
include system/ansi_c
import linalg

type
  GpuArrayObj*[V,M:static[int],T] = object
    p*: Coalesced[V,M,T]
    n*: int
  # GpuArrayRef*[V,M:static[int],T] = ref GpuArrayObj[V,M,T]
  # GpuArray*[V,M:static[int],T] = GpuArrayRef[V,M,T]
  # GpuArrays* = GpuArrayObj | GpuArrayRef
  # GpuArrays2* = GpuArrayObj | GpuArrayRef
  # GpuArrays3* = GpuArrayObj | GpuArrayRef

# Nim Bug, cannot overload this function with generic static parameters.
# proc init*(r: var GpuArrayObj, n: int) =
#   type T = r.T
#   var p: ptr T
#   when haveCuda:
#     let err = cudaMalloc(cast[ptr pointer](addr p), n*sizeof(T))
#     if err:
#       echo "alloc err: ", err
#       quit(-1)
#   else:
#     p = createSharedU(T, n)
#   r.n = n
#   r.p.newCoalesced(r.V, r.M, p, n)
# proc init[V,M:static[int],T](r: var GpuArrayRef[V,M,T], n: int) =
#   r.new
#   r[].init(n)

proc free*(r: var GpuArrayObj) =
  if r.n > 0:
    r.p.p.gpuFree
    r.n = 0
# proc free*[V,M:static[int],T](r: GpuArrayRef[V,M,T]) =
#   when haveCuda: discard r.p.p.cudaFree

proc initGpuArrayObj*(r: var GpuArrayObj, n: int) =
  type T = r.T
  var p: ptr T = cast[ptr T](gpuMalloc(csize_t(n*sizeof(T))))
  r.n = n
  r.p.initCoalesced(p, n)
  # echo "GpuArray init done."
proc newGpuArrayObj*(V,M:static[int], n:int, T:typedesc): auto {.noinit.} =
  var z {.noinit.}: GpuArrayObj[V,M,T]
  z.initGpuArrayObj(n)
  z

# proc newGpuArrayRef*[V,M:static[int],T](r: var GpuArrayRef[V,M,T], n: int) =
#   r.init(n)
# proc newGpuArrayRef*[T](V,M:static[int], n: int): auto {.noinit.} =
#   var z {.noinit.}: GpuArrayRef[V,M,T]
#   z.init(n)
#   z

template getGpuPtr*(x: SomeNumber): untyped = x
template getGpuPtr*(x: GpuArrayObj): untyped = x
# template getGpuPtr*(x: GpuArrayRef): untyped = x[]
#template getGpuPtr(x: GpuArrayRef): untyped = x.p
#template getGpuPtr(x: GpuArrayRef): untyped = (p:x.p,n:x.n)

template offloadUseVar*(x:GpuArrayObj):bool = true
template offloadUsePtr*(x:GpuArrayObj):bool = true
template rungpuPrepareOffload*(x:GpuArrayObj):bool = true
template runcpuFinalizeOffload*(x:GpuArrayObj):bool = false
template gpuVarPtr*(v:GpuArrayObj,p:untyped):untyped = v
template offloadPtr(x:GpuArrayObj):untyped = x.p.p
template offloadVar*(x:GpuArrayObj,p:untyped):untyped = x
template gpuPrepareOffload*(v:GpuArrayObj,pp:untyped):untyped = v.p.p = pp

template indexGpuArray*(x: GpuArrayObj, i: SomeInteger): untyped =
  x.p[i]

macro indexGpuArray*(x: GpuArrayObj{call}, y: SomeInteger): untyped =
  #echo "call[", y.repr, "]"
  #echo x.treerepr
  #if siteLocalsField.contains($x[0]):
  result = newCall(ident($x[0]))
  for i in 1..<x.len:
    let xi = x[i]
    result.add( quote do:
      indexGpuArray(`xi`,`y`) )
  #else:
  #  result = quote do:
  #    let tt = `x`
  #    tt.d[`y`]
  #echo result.treerepr
  #echo result.repr

template `[]`*(x: GpuArrayObj, i: SomeInteger): untyped = indexGpuArray(x, i)
# template `[]=`*(x: GpuArrayObj, i: SomeInteger, y: untyped): untyped =
#   x.p[][i] = y

# template `[]`*(x: GpuArrayRef, i: SomeInteger): untyped =
#   echo "GAR[]"
#   x.p[][i]
# template `[]=`*(x: GpuArrayRef, i: SomeInteger, y: untyped): untyped =
#   x.p[][i] = y

# var threadNum = 0
# var numThreads = 1
# template getThreadNum: untyped = threadNum
# template getNumThreads: untyped = numThreads
template `:=`*(x: GpuArrayObj, y: GpuArrayObj) =
  #cprintf("t %i/%i  b %i/%i\n", getThreadIdx(), getThreadDim(), getBlockIdx(), getBlockDim())
  #let i = getBlockDim().x * getBlockIdx().x + getThreadIdx().x
  mixin getThreadNum, getNumThreads
  let tid = getThreadNum()
  let nid = getNumThreads()
  var i = tid
  while i<x.n:
    x[i] := y[i]
    i += nid

template `:=`*(x: GpuArrayObj, y: SomeNumber) =
  #cprintf("t %i/%i  b %i/%i\n", getThreadIdx(), getThreadDim(), getBlockIdx(), getBlockDim())
  #let i = getBlockDim().x * getBlockIdx().x + getThreadIdx().x
  mixin getThreadNum, getNumThreads
  let tid = getThreadNum()
  let nid = getNumThreads()
  var i = tid
  while i<x.n:
    x[i] := y
    #echo i, "/", x.n
    i += nid

template `+=`*(x: GpuArrayObj, y: SomeNumber) =
  #cprintf("t %i/%i  b %i/%i\n", getThreadIdx(), getThreadDim(), getBlockIdx(), getBlockDim())
  #let i = getBlockDim().x * getBlockIdx().x + getThreadIdx().x
  mixin getThreadNum, getNumThreads
  let tid = getThreadNum()
  let nid = getNumThreads()
  var i = tid
  #cprintf("%i/%i\n", i, x.n)
  while i<x.n:
    x[i] += y
    #cprintf("%i/%i\n", i, x.n)
    i += nid

template `+=`*(x: GpuArrayObj, y: GpuArrayObj) =
  #cprintf("t %i/%i  b %i/%i\n", getThreadIdx(), getThreadDim(), getBlockIdx(), getBlockDim())
  #let i = getBlockDim().x * getBlockIdx().x + getThreadIdx().x
  mixin getThreadNum, getNumThreads
  let tid = getThreadNum()
  let nid = getNumThreads()
  var i = tid
  #cprintf("%i/%i\n", i, x.n)
  while i<x.n:
    x[i] += y[i]
    #cprintf("%i/%i\n", i, x.n)
    i += nid

proc `+`*[VX,VY,MX,MY:static[int],TX,TY](x: GpuArrayObj[VX,MX,TX], y: GpuArrayObj[VY,MY,TY]): auto =
  var r: GpuArrayObj[x.V,x.M,type(x[0]+y[0])]
  # when x is GpuArrayObj:
  #   var r: GpuArrayObj[x.V,x.M,type(x[0]+y[0])]
  # else:
  #   var r: GpuArrayRef[x.V,x.M,type(x[0]+y[0])]
  cprintf("+\n")
  r
proc `*`*[VX,VY,MX,MY:static[int],TX,TY](x: GpuArrayObj[VX,MX,TX], y: GpuArrayObj[VY,MY,TY]): auto =
  var r: GpuArrayObj[x.V,x.M,type(x[0]*y[0])]
  # when x is GpuArrayObj:
  #   var r: GpuArrayObj[x.V,x.M,type(x[0]*y[0])]
  # else:
  #   var r: GpuArrayRef[x.V,x.M,type(x[0]*y[0])]
  cprintf("*\n")
  r

when isMainModule:
  var N = 1000

  proc testfloat =
    # var x,y,z:GpuArrayObj[4,1,float32]
    # x.initGpuArrayObj(N)
    # y.initGpuArrayObj(N)
    # z.initGpuArrayObj(N)
    var x = newGpuArrayObj(4,1,N,float32)
    var y = newGpuArrayObj(4,1,N,float32)
    var z = newGpuArrayObj(4,1,N,float32)
    #cprintf("x.n: %i\n", x.n)
    onGpu(1,32):
      x += y * z
  testfloat()

when false:
  proc testcomplex =
    var x = newGpuArrayRef[Complex[float32]](N)
    var y = newGpuArrayRef[Complex[float32]](N)
    var z = newGpuArrayRef[Complex[float32]](N)
    onGpu(N):
      x += y * z
  testcomplex()

  proc testcolmat =
    var x = newGpuArrayRef[Colmat[3,float32]](N)
    var y = newGpuArrayRef[Colmat[3,float32]](N)
    var z = newGpuArrayRef[Colmat[3,float32]](N)
    #y := 1
    #z := 2
    onGpu(N):
      x += y * z
  testcolmat()
