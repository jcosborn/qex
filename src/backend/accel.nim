#const Backend {.strdefine.} = "OpenMP"
#const Backend {.strdefine.} = "CUDA"
#const Backend {.strdefine.} = "SYCL"
const Backend {.strdefine.} = "CPU"

when Backend == "OpenMP":
  #const useGPU = true
  import openmp
  export openmp
elif Backend == "CUDA":
  #const useGPU = true
  import cuda
  export cuda
elif Backend == "SYCL":
  #const useGPU = true
  import syclbe
  export syclbe
else:
  {.warning: "Backend unknown, use CPU only.".}
  #const useGPU = false
  import cpu
  export cpu

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
    echo "x: ", x

  test1()
