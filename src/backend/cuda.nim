import macros, strutils
import base/[metaUtils,profile]
import expr

{.pragma: cudah, header:"cuda_runtime.h".}

# vector types
type
  float2* {.importc,cudah,completeStruct.} = object
    x,y: float32
  double2* {.importc,cudah,completeStruct.} = object
    x,y: float64
  float3* {.importc,cudah,completeStruct.} = object
    x,y,z: float32
  double3* {.importc,cudah,completeStruct.} = object
    x,y,z: float64
  float4* {.importc,cudah,completeStruct.} = object
    x,y,z,w: float32
  double4* {.importc,cudah,completeStruct.} = object
    x,y,z,w: float64
  CudaVecTypes* = float2 | double2 | float3 | double3 | float4 | double4
  CudaTypes* = int32 | int64 | float32 | float64 | CudaVecTypes
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

proc cudaMemset*(devPtr: pointer, value: cint, count: csize_t):
  cudaError_t {.importC:"cudaMemset",header:"cuda_runtime.h".}
#template cudaMemsetX*(devPtr: pointer, value: SomeInteger, count: SomeInteger)

proc cudaMemcpyX*(dst,src: pointer, count: csize_t, kind: cudaMemcpyKind):
  cudaError_t {.importC:"cudaMemcpy",header:"cuda_runtime.h".}
template cudaMemcpy*(dst,src: typed, count: csize_t,
                     kind: cudaMemcpyKind): cudaError_t =
  let pdst = toPointer(dst)
  let psrc = toPointer(src)
  cudaMemcpyX(pdst, psrc, count, kind)

proc cudaGetDeviceCount*(deviceCount: ptr cint):
  cudaError_t {.importC,header:"cuda_runtime.h".}

proc cudaLaunchKernel(p:pointer, gd,bd: CudaDim3, args: ptr pointer):
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

template cudaDefs(body: untyped): untyped {.dirty.} =
  bind inlineProcs
  {.emit:["#define nimZeroMem(b,len) memset((b),0,(len))"].}
  {.emit:["#define nimCopyMem(a,b,len) memcpy((a),(b),(len))"].}
  {.pragma: shared, noInit, codegendecl:"__shared__ $# $#".}
  inlineProcs:
    body
  {.emit:["#undef nimZeroMem"].}
  {.emit:["#undef nimCopyMem"].}

template cudaLaunch*(p: proc {.cdecl.}; blocksPerGrid,threadsPerBlock: SomeInteger;
                     arg: varargs[pointer,dataAddr]) =
  var pp = pointer p
  var gridDim, blockDim: CudaDim3
  gridDim.x = cuint blocksPerGrid
  gridDim.y = 1
  gridDim.z = 1
  blockDim.x = cuint threadsPerBlock
  blockDim.y = 1
  blockDim.z = 1
  var args: array[arg.len, pointer]
  for i in 0..<arg.len: args[i] = arg[i]
  #echo "really launching kernel"
  let err = cudaLaunchKernel(pp, gridDim, blockDim, addr args[0])
  if err:
    echo "cudaLaunch: ", err
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
  result = newCall(ident("cudaLaunch"))
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

proc cudaproc(s:string, p:NimNode):NimNode =
  #echo "begin cuda:"
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
  result.body = getAst(cudaDefs(result.body))
  var sl = newStmtList()
  #sl.add( quote do:
  #  {.push checks: off.}
  #  {.push stacktrace: off.} )
  sl.add result
  result = sl
  #echo "end cuda:"
  #echo result.treerepr
macro cudaGlobal*(p: untyped): untyped = cudaproc("__global__",p)

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

  proc vectorAdd(A: FltArr; B: FltArr; C: FltArr; n: uint32) {.cdecl,cudaGlobal.} =
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
    cudaLaunch(vectorAdd, blocksPerGrid, threadsPerBlock, pa, pb, pc, un)
    discard cudaDeviceSynchronize()

    var errcnt = 0
    for i in 0..<n:
      let d = a[i] + b[i]
      if errcnt < 10 and c[i] != d:
        echo "error: ", i, "  ", c[i], "  ", d
        inc errcnt

  test()
