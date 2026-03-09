const Backend {.strdefine.} = "OpenMP"
#const Backend {.strdefine.} = "CUDA"
#const Backend {.strdefine.} = "SYCL"
#const Backend {.strdefine.} = "CPU"

when Backend == "OpenMP":
  const useGPU* = true
  import openmp
  export openmp
elif Backend == "CUDA":
  const useGPU* = true
  import cuda
  export cuda
elif Backend == "SYCL":
  const useGPU* = true
  #const useGPU* = false
  import syclbe
  export syclbe
else:
  when Backend != "CPU":
    static: echo "Backend: ", Backend
    {.warning: "Backend unknown, using CPU only.".}
  const useGPU* = false
  import cpu
  export cpu

proc printf*(frmt: cstring): cint {.
  importc: "printf", header: "<stdio.h>", varargs, discardable.}

proc gpuMalloc*[T](x: var ptr T) =
  let n = sizeof(T)
  x = cast[ptr T](gpuMalloc(n))

template toGpu*(x:SomeNumber):auto = x
template toGpu*(x:var SomeNumber):auto = unsafeAddr x
template getGpu*(x:SomeNumber, g:SomeNumber):auto = g
template getGpu*(x:SomeNumber, g:ptr SomeNumber):auto = g[]
template fromGpu*(x:SomeNumber, g:SomeNumber) = discard
template fromGpu*(x:var SomeNumber, g:ptr SomeNumber) = (x = g[])

iterator gpuRange*(n: int): int =
  when useGPU:
    let s = getNumThreads()
    var i = getThreadNum()
    while i < n:
      yield i
      i += s
  else:
    let s = getNumThreads()
    let id = getThreadNum()
    let i0 = (n*id) div s
    let i1 = (n*(id+1)) div s
    for i in i0 .. i1:
      yield i

#when useGPU:
#  import expr
#  import gpuarray
#  export gpuarray
#else:
#  template onGpu*(x:untyped) = threads: x
#  template onGpu*(n,x:untyped) = threads: x
#  template onGpu*(n,t,x:untyped) = threads: x
#  template packVarsStmt(x,y:untyped) = discard

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
