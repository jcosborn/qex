import macros, strutils
import base/[metaUtils,profile]
import expr

const dumpKernels {.intdefine.} = 0

#{.passC:"--with-arch=sm_86".}
#{.passL:"--with-arch=sm_86".}
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
template vecType[N:static int,T](t: typedesc[array[N,T]]): typedesc =
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
    dataAddr(a)
  elif x is array:
    vara a = addr x[0]
    dataAddr(a)
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

proc gpuInit*(device: int) =
  let err = cudaSetDevice(cint device)
  if err:
    echo "CUDA error: ", cast[cint](err)
  doAssert(not err)
template gpuMalloc*(size: SomeInteger):pointer =
  var p:pointer
  let err = cudaMalloc(p, csize_t size)
  if err:
    echo "gpuMalloc: ", err
    p = cast[pointer](0)
  p
template gpuMalloc[T](x: var ptr UncheckedArray[T], n: int) =
  x = cast[typeof x](gpuMalloc(n*sizeof(T)))
template gpuMalloc[T](x: ptr T) =
  x = cast[typeof x](gpuMalloc(sizeof(T)))
template gpuFree*(p:pointer) =
  let err = cudaFree(p)
  if err:
    echo err
    quit cast[cint](err)
proc gpuMemset*(devPtr: pointer, value: SomeInteger, count: SomeInteger) =
  let err = cudaMemset(devPtr, cint value, csize_t count)
  if err:
    echo "gpuMemset: ", err
#proc gpuMemCpyToGpu*(dst,src: pointer, count: SomeInteger):cint {.discardable.} =
proc gpuMemCpyToGpu*(dst,src: pointer, count: SomeInteger) =
  let err = cudaMemcpy(dst,src,csize_t count,cudaMemcpyHostToDevice)
  if err:
    echo "gpuMemCpyToGpu: ", err
#proc gpuMemCpyToCpu*(dst,src: pointer, count: SomeInteger):cint {.discardable.} =
template gpuMemCpyToCpu*(dst,src: pointer, count: SomeInteger) =
  let err = cudaMemcpy(dst,src,csize_t count,cudaMemcpyDeviceToHost)
  if err:
    echo instantiationInfo()
    echo "gpuMemCpyToCpu: ", err

proc cudaLaunchKernel(p:pointer, gd,bd: CudaDim3, args: ptr pointer):
  cudaError_t {.importC,header:"cuda_runtime.h".}

proc syncThreads() {.importc:"__syncthreads",header:"cuda_runtime.h".}
proc threadFence() {.importc:"__threadfence",header:"cuda_runtime.h".}
proc atomicInc(address: ptr cuint, val: cuint): cuint {.importc,header:"cuda_runtime.h".}
template atomicInc(address: ptr cuint, val: SomeInteger): auto =
  atomicInc(address, cuint val)
proc cudaDeviceReset*(): cudaError_t
  {.importC,header:"cuda_runtime.h".}
proc cudaDeviceSynchronize*(): cudaError_t
  {.importC,header:"cuda_runtime.h".}

#proc printf*(fmt:cstring):cint {.importc,varargs,header:"<stdio.h>",discardable.}
#proc fprintf*(stream:ptr FILE,fmt:cstring):cint {.importc,varargs,header:"<stdio.h>".}
#proc malloc*(size: csize_t):pointer {.importc,header:"<stdlib.h>".}

var gridDim*{.importC,header:"cuda_runtime.h".}: CudaDim3
var blockDim*{.importC,header:"cuda_runtime.h".}: CudaDim3
var blockIdx*{.importC,header:"cuda_runtime.h".}: CudaDim3
var threadIdx*{.importC,header:"cuda_runtime.h".}: CudaDim3
template getGridDim*: auto = gridDim
template getBlockIdx*: auto = blockIdx
template getBlockDim*: auto = blockDim
template getThreadIdx*: auto = threadIdx
template gpuThreadNum*: auto =
  blockDim.x * blockIdx.x + threadIdx.x
template gpuNumThreads*: auto =
  gridDim.x * blockDim.x

template cudaDefs(body: untyped): untyped {.dirty.} =
  #var gridDim{.global,importC,noDecl.}: CudaDim3
  #var blockIdx{.global,importC,noDecl.}: CudaDim3
  #var blockDim{.global,importC,noDecl.}: CudaDim3
  #var threadIdx{.global,importC,noDecl.}: CudaDim3
  #template getGridDim: untyped {.used.} = gridDim
  #template getBlockIdx: untyped {.used.} = blockIdx
  #template getBlockDim: untyped {.used.} = blockDim
  #template getThreadIdx: untyped {.used.} = threadIdx
  #template getThreadNum: untyped {.used.} = blockDim.x * blockIdx.x + threadIdx.x
  #template getNumThreads: untyped {.used.} = gridDim.x * blockDim.x
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

#macro onGpuX(nn0,tpb0: untyped, body: untyped): auto =  # return callback to wait and finalize

macro onGpuNowait(nn0,tpb0: untyped, body: untyped): auto =
  let li = body.lineinfo
  let lis = li.split({'/','.','(',','})
  let fl = lis[^4] & "(" & lis[^2] & ")"
  template target(fl, nn, tpb, v, arg, cpuPrepare, cpuFinalize, body: untyped): untyped =
    mixin toGpu, getGpu, fromGpu
    block:
      tic(fl)
      var v = cpuPrepare  # kernel argument tuple
      toc("cpuPrepare")
      type ByCopy[T] {.bycopy.} = object
        d: T
      proc kern(arg: ByCopy[type(v)]) {.cdecl,cudaGlobal.} =
        const inOnGpu {.inject,used.} = true
        {.push checks: off.}
        {.push stacktrace: off.}
        body
      let ni = nn.int32
      let threadsPerBlock = tpb.int32
      let blocksPerGrid = (ni+threadsPerBlock-1) div threadsPerBlock
      #echo "launching kernel"
      threadSingle:
        cudaLaunch(kern, blocksPerGrid, threadsPerBlock, v)
      toc("launch")
      proc finalize {.gensym.} =
        tic(fl)
        threadSingle:
          discard cudaDeviceSynchronize()
        toc("wait")
        cpuFinalize
        #threadBarrier()
        toc("cpuFinalize")
      finalize
  let
    varg = gensym(nskVar, "varg")
    arg = gensym(nskParam, "arg")
  proc deref(x,g,i:NimNode):auto =
    newCall("getGpu",x,
            newNimNode(nnkBracketExpr).add(
              newNimNode(nnkDotExpr).add(arg,ident"d"), i))
  let
    v = prepareVars(body,deref)  # gather gpu pointers in symbols, body is changed accordingly
    cpuPrepare = genCpuPrepare v
    cpuFinalize = genCpuFinalize(v, varg)
  result = getast(target(fl, nn0, tpb0, varg, arg, cpuPrepare, cpuFinalize, body))
  case dumpKernels
  of 1:
    echo li
    echo result.repr
  of 2:
    echo li
    echo result.treerepr
  else:
    if dumpKernels > 2:
      echo li
      var sl = newNimNode(nnkStmtListExpr)
      sl.add newCall(bindsym"echoTyped", result)
      sl.add result
      result = sl

var gpuBlockSizeRequest* = 64
var gpuNumThreadsRequest* = 32*1024
template gpuSites(n: int): int = n
template onGpuNowait*(body: untyped): auto =
  onGpuNoWait(gpuNumThreadsRequest, gpuBlockSizeRequest, body)
template onGpuNowait*(n0,body: untyped): auto =
  mixin gpuSites
  let n = gpuSites(n0)
  var b = gpuBlockSizeRequest
  while b > n: b = b div 2
  onGpuNoWait(n, b, body)
#template onGpuNowait*(n,b,body: untyped): auto =
#  onGpuNoWait(n, b, body)

template onGpu*(body: untyped) =
  let finalize = onGpuNoWait(body)
  finalize()
template onGpu*(n,body: untyped) =
  mixin gpuSites
  let finalize = onGpuNoWait(gpuSites(n), body)
  finalize()
template onGpu*(n,b,body: untyped) =
  mixin gpuSites
  let finalize = onGpuNoWait(gpuSites(n), b, body)
  finalize()

#{.emit:"/*INCLUDESECTION*/#include <cooperative_groups.h>".}
proc importCG() {.header:"<cooperative_groups.h>",importcpp:"@".}
proc warpSumSmall*[T:CudaTypes](x: T): T =  # for up to 32 bytes
  importCG()
  result = x
  var t {.noInit.}: T
  {.emit:["""
  namespace cg = cooperative_groups;
  const int warp_size = 32;
  cg::thread_block cta = cg::this_thread_block();
  cg::thread_block_tile<warp_size> tile = cg::tiled_partition<warp_size>(cta);
  for (int offset = warp_size / 2; offset >= 1; offset /= 2) {""",
  t,""" = tile.shfl_down(""",result,""", offset);"""].}
  result += t
  {.emit:["}"].}

template atomic_type*(t: CudaTypes): typedesc = t
template atomic_type*[N:static int,T](t: typedesc[array[N,T]]): typedesc =
  mixin atomic_type
  type vt = vecType(t)
  when vt is CudaTypes:
    vt
  else:
    when N mod 2 == 0:
      atomic_type(array[N div 2, T])
    else:
      T
proc cmemcpy(dest,src: pointer, count: csize_t): pointer {.importc:"memcpy",header:"string.h".}
proc warpSumLarge*[T](x: T): T =  # for over 32 bytes
  importCG()
  type atomic_t = atomic_type(T)
  static: echo "type: ", $T, "  atomic_t: ", $atomic_t
  const n = sizeof(T) div sizeof(atomic_t)
  #doAssert(sizeof(T) == n * sizeof(atomic_t))
  var sum_tmp {.noInit.}: array[n, atomic_t]
  discard cmemcpy(addr sum_tmp, addr x, csize_t sizeof(T))
  for i in 0..<n:
    sum_tmp[i] = warpSumSmall(sum_tmp[i])
  discard cmemcpy(addr result, addr sum_tmp, csize_t sizeof(T))

template warpSum*[T](x: T): T =
  when sizeof(T) <= 32 and T is CudaTypes:
    warpSumSmall(x)
  else:
    warpSumLarge(x)

#proc warpBroadcast*[T](x: T): T =
#  cg::thread_block cta = cg::this_thread_block();
#  cg::thread_block_tile<warp_size> tile = cg::tiled_partition<warp_size>(cta);
#  result = tile.shfl(x, 0);

proc blockSumSmall*[T](x: T): T = # only thread 0 gets result
  const max_block_size = 1024
  const warp_size = 32
  const max_items = max_block_size div warp_size
  let thread_idx = threadIdx.x
  let block_size = blockDim.x
  let warp_idx = thread_idx div warp_size
  let warp_items = (block_size + warp_size - 1) div warp_size
  # first do warp reduce
  result = warpSum(x)
  if warp_items == 1: return
  # now do reduction between warps
  syncThreads()
  var storage {.shared.}: array[max_items, T]
  # if first thread in warp, write result to shared memory
  if thread_idx mod warp_size == 0: storage[warp_idx] = result
  syncThreads()
  if warp_idx == 0:
    result = if thread_idx < warp_items: storage[thread_idx] else: default(T)
    result = warpSum(result)

proc blockSum*[T](x: T): T = # only thread 0 gets result
  const min_shared_mem = 48*1024 - sizeof(bool)  # bool used in GpuSum reduce
  const max_items = 32
  const max_size = min_shared_mem div max_items
  when sizeof(T) <= max_size:
    blockSumSmall(x)
  else:
    static: echo $x.type, "  ", sizeof(T)
    {.error:"blockSum: type size too large".}  # FIXME later

type GpuSum*[T] = object
    partial: ptr UncheckedArray[T]
    npartial: cuint
    #maxblock: int
    val: ptr T
    count: ptr cuint
proc newGpuSum*[T](ns: int): GpuSum[T] =
  let n = (ns + 31) div 32  # divide by warp size
  result.partial.gpuMalloc(n)
  result.partial.gpuMemset(0, n*sizeof(T))
  result.npartial = cuint n
  let err = cudaMallocHost((ptr pointer)(addr result.val), csize_t sizeof(T))
  result.count.gpuMalloc()
  result.count.gpuMemset(0, sizeof(result.count[]))
template value*(x: GpuSum): auto = x.val[]
template toGpu*(x: GpuSum): auto = x
template getGpu*(x,g: GpuSum): auto = g
template fromGpu*(x,g: GpuSum): auto = discard

proc reduce*[T](gs: GpuSum[T], x: T) =
  var isLastBlockDone {.shared.}: bool
  var aggregate = blockSum(x)
  if threadIdx.x == 0:
    if blockIdx.x < gs.npartial:
      gs.partial[blockIdx.x] = aggregate;
    threadFence() # flush result
    # increment global block counter
    let value = atomicInc(gs.count, gridDim.x)
    # determine if last block
    isLastBlockDone = (value == (gridDim.x - 1))
  syncThreads()
  # finish the reduction if last block
  if isLastBlockDone:
    var i = threadIdx.x
    var sum = default(T)
    let n = min(gs.npartial, gridDim.x)
    while i < n:
      sum += gs.partial[i]
      i += blockDim.x;
    sum = blockSum(sum)
    # write out the final reduced value
    if threadIdx.x == 0:
      gs.val[] = sum
      gs.count[] = 0  # set to zero for next time

when isMainModule:
  type FltArr = UncheckedArray[float32]

  #proc vectorAdd(A: FltArr; B: FltArr; C: var FltArr; n: int32) {.cdecl,cudaGlobal.} =
  #  var i = blockDim.x * blockIdx.x + threadIdx.x
  #  if i < n:
  #    C[i] = A[i] + B[i]

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

    #cudaLaunch(vectorAdd, blocksPerGrid, threadsPerBlock, a, b, c, n)
    #discard cudaDeviceSynchronize()

    template toGpu(x: SomeNumber): auto = x
    template getGpu(x,g: SomeNumber): auto = g
    template fromGpu(x,g: SomeNumber) = discard

    proc toGpu[T](x: seq[T]): ptr UncheckedArray[T] =
      let n = x.len
      result = cast[typeof result](gpuMalloc(n*sizeof(T)))
      #echo "p: ", cast[int](result)
      gpuMemCpyToGpu(result, addr x[0], n*sizeof(T))
    template getGpu(x: seq, g: ptr UncheckedArray): auto = g
    template fromGpu[T](x: seq[T], g: ptr UncheckedArray[T]) =
      let n = x.len
      gpuMemCpyToCpu(addr x[0], g, n*sizeof(T))
      # should free g

    var gs = newGpuSum[float](n)
    onGpu(n):
      let i = int getBlockDim().x * getBlockIdx().x + getThreadIdx().x
      var r: typeof c[0]
      if i < n:
        c[i] = a[i] + b[i]
        r = c[i]
      gs.reduce r
      #let s = warpSumSmall(c[i])
    echo gs.val[]
    #echo 3*n

  test()
