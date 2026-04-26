import macros, strutils
import base/[metaUtils,profile]
import expr

{.pragma: hiph, header:"hip/hip_runtime.h".}

# vector types
type
  float2* {.importc,hiph,completeStruct.} = object
    x,y: float32
  double2* {.importc,hiph,completeStruct.} = object
    x,y: float64
  float3* {.importc,hiph,completeStruct.} = object
    x,y,z: float32
  double3* {.importc,hiph,completeStruct.} = object
    x,y,z: float64
  float4* {.importc,hiph,completeStruct.} = object
    x,y,z,w: float32
  double4* {.importc,hiph,completeStruct.} = object
    x,y,z,w: float64
  HipVecTypes* = float2 | double2 | float3 | double3 | float4 | double4
  HipTypes* = int32 | int64 | float32 | float64 | HipVecTypes
template vecType*[N:static int,T](t: typedesc[array[N,T]]): typedesc =
  when T is float32:
    when N == 1: float32
    elif N == 2: float2
    elif N == 3: float3
    elif N == 4: float4
    else: t
  elif T is float64:
    when N == 1: float64
    elif N == 2: double2
    elif N == 3: double3
    elif N == 4: double4
    else: t
  else:
    t
template `+=`*(a,b: float2 | double2) =
  a.x += b.x
  a.y += b.y
template `+=`*(a,b: float3 | double3) =
  a.x += b.x
  a.y += b.y
  a.z += b.z
template `+=`*(a,b: float4 | double4) =
  a.x += b.x
  a.y += b.y
  a.z += b.z
  a.w += b.w

proc addChildrenFrom*(dst,src: NimNode): NimNode =
  for c in src: dst.add(c)
  result = dst
macro procInst*(p: typed): auto =
  #echo "begin procInst:"
  #echo p.treerepr
  result = p[0]

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
#template dataAddr*(x: typed): pointer =
#  #dumpType: x
#  when x is seq: dataAddr(x[0])
#  elif x is array: dataAddr(x[0])
#  #elif x is ptr: x
#  else: pointer(unsafeAddr(x))
#  #else: x
template dataAddr*(x: typed): pointer =
  when x is seq:
    var a = addr x[0]
    #dataAddr(a)
    a
  elif x is array:
    vara a = addr x[0]
    #dataAddr(a)
    a
  else: pointer(addr x)

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

proc hipMemset*(devPtr: pointer, value: cint, count: csize_t):
  hipError_t {.importC:"hipMemset",hiph.}
#template hipMemsetX*(devPtr: pointer, value: SomeInteger, count: SomeInteger)

proc hipMemcpyX*(dst,src: pointer, count: csize_t, kind: hipMemcpyKind):
  hipError_t {.importC:"hipMemcpy",hiph.}
template hipMemcpy*(dst,src: typed, count: csize_t,
                     kind: hipMemcpyKind): hipError_t =
  let pdst = toPointer(dst)
  let psrc = toPointer(src)
  hipMemcpyX(pdst, psrc, count, kind)

proc hipGetDeviceCount*(deviceCount: ptr cint):
  hipError_t {.importC,hiph.}

proc hipLaunchKernel(p:pointer, gd,bd: HipDim3, args: ptr pointer, sharedMemBytes: csize_t, stream: hipStream_t):
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

template hipDefs(body: untyped): untyped {.dirty.} =
  bind inlineProcs
  {.emit:["#define nimZeroMem(b,len) memset((b),0,(len))"].}
  {.emit:["#define nimCopyMem(a,b,len) memcpy((a),(b),(len))"].}
  {.pragma: shared, noInit, codegendecl:"__shared__ $# $#".}
  inlineProcs:
    body
  {.emit:["#undef nimZeroMem"].}
  {.emit:["#undef nimCopyMem"].}

template hipLaunch*(p: proc {.cdecl.}; blocksPerGrid,threadsPerBlock: SomeInteger;
                     arg: varargs[pointer,dataAddr]) =
  var pp = pointer p
  var gridDim, blockDim: HipDim3
  gridDim.x = cuint blocksPerGrid
  gridDim.y = 1
  gridDim.z = 1
  blockDim.x = cuint threadsPerBlock
  blockDim.y = 1
  blockDim.z = 1
  var args: array[arg.len, pointer]
  for i in 0..<arg.len: args[i] = arg[i]
  #echo "really launching kernel"
  let err = hipLaunchKernel(pp, gridDim, blockDim, addr args[0], 0, 0)
  if err:
    echo "hipLaunch: ", err
    quit cast[cint](err)

template `<<`*(p: proc, x: tuple): untyped = (p,x)
template getInst*(p: untyped): untyped =
  #when compiles((var t=p; t)): p
  #else:
  procInst(p)
    #var t =
    #t
macro `>>`*(px: tuple, y: auto): auto =
  #echo "begin >>:"
  #echo px.treerepr
  #echo "kernel type:"
  #echo px[0].getTypeImpl.treerepr
  #echo "kernel args:"
  #echo y.treerepr
  #var a = y
  #if y.kind != nnkPar: a = newNimNode(nnkPar).addChildrenFrom(y)
  result = newCall(ident("hipLaunch"))
  let krnl = newCall(px[0]).addChildrenFrom(y)
  #echo "kernel inst call:"
  #echo krnl.treerepr
  result.add getAst(getInst(krnl))[0]
  result.add px[1][0]
  result.add px[1][1]
  for c in y:
    result.add c
  #echo "kernel launch body:"
  #echo result.treerepr

proc hipproc(s:string, p:NimNode):NimNode =
  #echo "begin hip:"
  #echo s
  #let ss = s.strVal
  #echo "proc:"
  #echo p.treerepr
  p.expectKind nnkProcDef
  result = p
  # if p.kind == nnkProcDef:
  #   result = p
  # else:
  #   result = p[0]
  result.addPragma parseExpr("{.codegenDecl:\""&s&" $# $#$#\".}")[0]
  result.body = getAst(hipDefs(result.body))
  var sl = newStmtList()
  #sl.add( quote do:
  #  {.push checks: off.}
  #  {.push stacktrace: off.} )
  sl.add result
  result = sl
  #echo "end hip:"
  #echo result.treerepr
macro hipGlobal*(p: untyped): untyped = hipproc("__global__",p)

proc genCpuPrepare(n:seq[NimNode]):NimNode =
  result = newNimNode(nnkTupleConstr)
  for c in n:
    result.add newCall(ident"toGpu", c[0])

proc genCpuFinalize(n:seq[NimNode], a: NimNode):NimNode =
  template r(a,x,i:untyped):untyped =
    fromGpu(x,a[i])
  result = newstmtlist()
  for c in n:
    result.add getast r(a,c[0],c[2])


when isMainModule:
  type FltArr = ptr UncheckedArray[float32]

  proc vectorAdd(A: FltArr; B: FltArr; C: FltArr; n: uint32) {.cdecl,hipGlobal.} =
    var i = blockDim.x * blockIdx.x + threadIdx.x
    if i < n:
      C[i] = A[i] + B[i]

  proc test =
    var n = 50000
    var
      a = newSeq[float32](n)
      b = newSeq[float32](n)
      c = newSeq[float32](n)
    for i in 0..<n:
      a[i] = 1
      b[i] = 2

    var threadsPerBlock = 256
    var blocksPerGrid = (n + threadsPerBlock - 1) div threadsPerBlock

    let pa = addr a[0]
    let pb = addr b[0]
    let pc = addr c[0]
    let un = uint32 n
    hipLaunch(vectorAdd, blocksPerGrid, threadsPerBlock, pa, pb, pc, un)
    discard hipDeviceSynchronize()

    var errcnt = 0
    for i in 0..<n:
      let d = a[i] + b[i]
      if errcnt < 10 and c[i] != d:
        echo "error: ", i, "  ", c[i], "  ", d
        inc errcnt
    if errcnt == 0:
      echo "Test passed"

  test()
