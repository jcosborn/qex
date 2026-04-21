import macros, strutils
import base/qexInternal

#const Backend {.strdefine.} = "OpenMP"
#const Backend {.strdefine.} = "CUDA"
#const Backend {.strdefine.} = "SYCL"
const Backend {.strdefine.} = "CPU"

when Backend == "OpenMP":
  const backendIsGpu* = true
  const beSharedMem* {.booldefine.} = false
  import openmp
  export openmp
elif Backend == "CUDA":
  const backendIsGpu* = true
  const beSharedMem* {.booldefine.} = false
  import cudabe
  export cudabe
  import cudasum
  export cudasum
  proc init = gpuInit(myRank)
  qexGlobalInitializers.add init
elif Backend == "SYCL":
  const backendIsGpu* = true
  const beSharedMem* {.booldefine.} = false
  #const backendIsGpu* = false
  import syclbe
  export syclbe
  proc init =
    gpuInit(myRank)
    echo "SYCL device: ", gpuDeviceName()
  qexGlobalInitializers.add init
else:
  when Backend != "CPU":
    static: echo "Backend: ", Backend
    {.warning: "Backend unknown, using CPU only.".}
  const backendIsGpu* = false
  const beSharedMem* {.booldefine.} = true
  import cpu
  export cpu
  proc init =
    echo "Using CPU backend"
  qexGlobalInitializers.add init
const backendIsCpu* = not backendIsGpu

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
