{.pragma: hiph, header:"hip/hip_runtime.h".}

type
  HipDim3* {.importc:"dim3",hiph.} = object
    x*, y*, z*: cuint
  hipError_t* {.importc,hiph.} = object
  hipStream_t* {.importc,hiph.} = cint
  hipMemcpyKind* {.importc,hiph.} = object
var
  hipSuccess*{.importC,hiph.}: hipError_t
  hipErrorNotSupported*{.importC,hiph.}: hipError_t
  hipMemcpyHostToDevice*{.importC,hiph.}: hipMemcpyKind
  hipMemcpyDeviceToHost*{.importC,hiph.}: hipMemcpyKind

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

proc hipSetDevice*(device: cint): hipError_t
  {.importC,hiph.}
proc hipGetLastError*(): hipError_t
  {.importC,hiph.}
proc hipGetErrorStringX*(error: hipError_t): ptr char
  {.importC:"hipGetErrorString",hiph.}
proc hipGetErrorString*(error: hipError_t): cstring =
  var s {.codegendecl:"const $# $#".} = hipGetErrorStringX(error)
  result = cast[cstring](s)
proc `$`*(error: hipError_t): string =
  let s = hipGetErrorString(error)
  result = $s
converter toBool*(e: hipError_t): bool =
  cast[cint](e) != cast[cint](hipSuccess)

proc hipMalloc*(p:ptr pointer, size: csize_t): hipError_t
  {.importC,hiph.}
template hipMalloc*(p:pointer, size: csize_t): hipError_t =
  hipMalloc(p.addr, size)
proc hipFree*(p: pointer): hipError_t
  {.importC,hiph.}
proc hipMallocManaged*(p: ptr pointer, size: csize_t): hipError_t
  {.importC,hiph.}
proc hipMallocHost*(p: ptr pointer, size: csize_t): hipError_t
  {.importC,hiph.}
proc hipHostMalloc*(p: ptr pointer, size: csize_t): hipError_t
  {.importC,hiph.}
proc hipHostFree*(p: pointer): hipError_t
  {.importC,hiph.}

proc hipMemset*(devPtr: pointer, value: cint, count: csize_t):
  hipError_t {.importC:"hipMemset",hiph.}
#template hipMemsetX*(devPtr: pointer, value: SomeInteger, count: SomeInteger)

proc hipMemcpyAsync*(dst, src: pointer; count: csize_t; kind: hipMemcpyKind; stream: pointer): hipError_t
  {.importcpp: "hipMemcpyAsync(#, #, #, #, (hipStream_t)#)", hiph.}
proc hipMemcpyX*(dst,src: pointer, count: csize_t, kind: hipMemcpyKind):
  hipError_t {.importC:"hipMemcpy",hiph.}
template hipMemcpy*(dst,src: typed, count: csize_t,
                     kind: hipMemcpyKind): hipError_t =
  let pdst = toPointer(dst)
  let psrc = toPointer(src)
  hipMemcpyX(pdst, psrc, count, kind)

proc hipGetDeviceCount*(deviceCount: ptr cint):
  hipError_t {.importC,hiph.}

proc syncThreads*() {.importc:"__syncthreads",hiph.}
proc threadFence*() {.importc:"__threadfence",hiph.}
proc atomicInc*(address: ptr cuint, val: cuint): cuint {.importc,hiph.}
template atomicInc*(address: ptr cuint, val: SomeInteger): auto =
  atomicInc(address, cuint val)
proc hipDeviceReset*(): hipError_t
  {.importC,hiph.}
proc hipDeviceSynchronize*(): hipError_t
  {.importC,hiph.}

var gridDim*{.importC,hiph.}: HipDim3
var blockDim*{.importC,hiph.}: HipDim3
var blockIdx*{.importC,hiph.}: HipDim3
var threadIdx*{.importC,hiph.}: HipDim3
template getGridDim*: auto = gridDim
template getBlockIdx*: auto = blockIdx
template getBlockDim*: auto = blockDim
template getThreadIdx*: auto = threadIdx
