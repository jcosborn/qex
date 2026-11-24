{.pragma: cudah, header:"cuda_runtime.h".}

type
  CudaDim3* {.importc:"dim3",header:"cuda_runtime.h".} = object
    x*, y*, z*: cuint
  cudaError_t* {.importc,header:"cuda_runtime.h".} = object
  cudaMemcpyKind* {.importc,header:"cuda_runtime.h".} = object
var
  cudaSuccess*{.importC,header:"cuda_runtime.h".}: cudaError_t
  cudaErrorNotSupported*{.importC,header:"cuda_runtime.h".}: cudaError_t
  cudaMemcpyHostToDevice*{.importC,header:"cuda_runtime.h".}: cudaMemcpyKind
  cudaMemcpyDeviceToHost*{.importC,header:"cuda_runtime.h".}: cudaMemcpyKind

#template toPointer*(x: pointer): pointer = x
#template toPointer*[T](x: ptr T): pointer = pointer(x)
#template toPointer*(x: seq): pointer = toPointer(x[0])
#template toPointer*(x: not (pointer|seq)): pointer = pointer(unsafeAddr(x))
template toPointer*(x: typed): pointer =
  #dumpType: x
  when x is pointer: x
  elif x is ptr: x
  elif x is seq: toPointer(x[0])
  else: pointer(unsafeAddr(x))

proc cudaSetDevice*(device: cint): cudaError_t
  {.importC,header:"cuda_runtime.h".}
proc cudaGetLastError*(): cudaError_t
  {.importC,header:"cuda_runtime.h".}
proc cudaGetErrorStringX*(error: cudaError_t): ptr char
  {.importC:"cudaGetErrorString",header:"cuda_runtime.h".}
proc cudaGetErrorString*(error: cudaError_t): cstring =
  var s {.codegendecl:"const $# $#".} = cudaGetErrorStringX(error)
  result = cast[cstring](s)
proc `$`*(error: cudaError_t): string =
  let s = cudaGetErrorString(error)
  result = $s
converter toBool*(e: cudaError_t): bool =
  cast[cint](e) != cast[cint](cudaSuccess)

proc cudaMalloc*(p:ptr pointer, size: csize_t): cudaError_t
  {.importC,header:"cuda_runtime.h".}
template cudaMalloc*(p:pointer, size: csize_t): cudaError_t =
  cudaMalloc(p.addr, size)
proc cudaFree*(p: pointer): cudaError_t
  {.importC,header:"cuda_runtime.h".}
proc cudaMallocManaged*(p: ptr pointer, size: csize_t): cudaError_t
  {.importC,header:"cuda_runtime.h".}
proc cudaMallocHost*(p: ptr pointer, size: csize_t): cudaError_t
  {.importC,header:"cuda_runtime.h".}
proc cudaFreeHost*(p: pointer): cudaError_t
  {.importC,header:"cuda_runtime.h".}

proc cudaMemset*(devPtr: pointer, value: cint, count: csize_t):
  cudaError_t {.importC:"cudaMemset",header:"cuda_runtime.h".}
#template cudaMemsetX*(devPtr: pointer, value: SomeInteger, count: SomeInteger)

proc cudaMemcpyAsync*(dst, src: pointer; count: csize_t; kind: cudaMemcpyKind; stream: pointer): cudaError_t
  {.importcpp: "cudaMemcpyAsync(#, #, #, #, (cudaStream_t)#)", header: "cuda_runtime.h".}
proc cudaMemcpyX*(dst,src: pointer, count: csize_t, kind: cudaMemcpyKind):
  cudaError_t {.importC:"cudaMemcpy",header:"cuda_runtime.h".}
template cudaMemcpy*(dst,src: typed, count: csize_t,
                     kind: cudaMemcpyKind): cudaError_t =
  let pdst = toPointer(dst)
  let psrc = toPointer(src)
  cudaMemcpyX(pdst, psrc, count, kind)

proc cudaGetDeviceCount*(deviceCount: ptr cint):
  cudaError_t {.importC,header:"cuda_runtime.h".}

proc syncThreads*() {.importc:"__syncthreads",header:"cuda_runtime.h".}
proc threadFence*() {.importc:"__threadfence",header:"cuda_runtime.h".}
proc atomicInc*(address: ptr cuint, val: cuint): cuint {.importc,header:"cuda_runtime.h".}
template atomicInc*(address: ptr cuint, val: SomeInteger): auto =
  atomicInc(address, cuint val)
proc cudaDeviceReset*(): cudaError_t
  {.importC,header:"cuda_runtime.h".}
proc cudaDeviceSynchronize*(): cudaError_t
  {.importC,header:"cuda_runtime.h".}

var gridDim*{.importC,header:"cuda_runtime.h".}: CudaDim3
var blockDim*{.importC,header:"cuda_runtime.h".}: CudaDim3
var blockIdx*{.importC,header:"cuda_runtime.h".}: CudaDim3
var threadIdx*{.importC,header:"cuda_runtime.h".}: CudaDim3
template getGridDim*: auto = gridDim
template getBlockIdx*: auto = blockIdx
template getBlockDim*: auto = blockDim
template getThreadIdx*: auto = threadIdx
