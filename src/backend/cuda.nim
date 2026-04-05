import macros, strutils
import base/metaUtils
import expr

const dumpKernels {.intdefine.} = 0

#{.passC:"--with-arch=sm_86".}
#{.passL:"--with-arch=sm_86".}

proc addChildrenFrom*(dst,src: NimNode): NimNode =
  for c in src: dst.add(c)
  result = dst
macro procInst*(p: typed): auto =
  #echo "begin procInst:"
  #echo p.treerepr
  result = p[0]

type
  CudaDim3* {.importc:"dim3",header:"cuda_runtime.h".} = object
    x*, y*, z*: cint
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
template dataAddr*(x: typed): pointer =
  #dumpType: x
  when x is seq: dataAddr(x[0])
  elif x is array: dataAddr(x[0])
  #elif x is ptr: x
  else: pointer(unsafeAddr(x))
  #else: x

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
    echo err
    p = cast[pointer](0)
  p
template gpuFree*(p:pointer) =
  let err = cudaFree(p)
  if err:
    echo err
    quit cast[cint](err)
#proc gpuMemCpyToGpu*(dst,src: pointer, count: SomeInteger):cint {.discardable.} =
proc gpuMemCpyToGpu*(dst,src: pointer, count: SomeInteger) =
  let err = cudaMemcpy(dst,src,csize_t count,cudaMemcpyHostToDevice)
  if err:
    echo err
#proc gpuMemCpyToCpu*(dst,src: pointer, count: SomeInteger):cint {.discardable.} =
proc gpuMemCpyToCpu*(dst,src: pointer, count: SomeInteger) =
  let err = cudaMemcpy(dst,src,csize_t count,cudaMemcpyDeviceToHost)
  if err:
    echo err

proc cudaLaunchKernel(p:pointer, gd,bd: CudaDim3, args: ptr pointer):
  cudaError_t {.importC,header:"cuda_runtime.h".}

proc cudaDeviceReset*(): cudaError_t
  {.importC,header:"cuda_runtime.h".}
proc cudaDeviceSynchronize*(): cudaError_t
  {.importC,header:"cuda_runtime.h".}

#proc printf*(fmt:cstring):cint {.importc,varargs,header:"<stdio.h>",discardable.}
#proc fprintf*(stream:ptr FILE,fmt:cstring):cint {.importc,varargs,header:"<stdio.h>".}
#proc malloc*(size: csize_t):pointer {.importc,header:"<stdlib.h>".}

#template getGridDim*: untyped {.used.} =
#  var gridDim{.global,importC,noDecl.}: CudaDim3
#  gridDim
#template getBlockIdx*: untyped {.used.} =
#  var blockIdx{.global,importC,noDecl.}: CudaDim3
#  blockIdx
#template getBlockDim*: untyped {.used.} =
#  var blockDim{.global,importC,noDecl.}: CudaDim3
#  blockDim
#template getThreadIdx*: untyped {.used.} =
#  var threadIdx{.global,importC,noDecl.}: CudaDim3
#  threadIdx
var gridDim*{.importC,header:"cuda_runtime.h".}: CudaDim3
var blockDim*{.importC,header:"cuda_runtime.h".}: CudaDim3
var blockIdx*{.importC,header:"cuda_runtime.h".}: CudaDim3
var threadIdx*{.importC,header:"cuda_runtime.h".}: CudaDim3
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
  inlineProcs:
    body
  {.emit:["#undef nimZeroMem"].}
  {.emit:["#undef nimCopyMem"].}

template cudaLaunch*(p: proc {.cdecl.}; blocksPerGrid,threadsPerBlock: SomeInteger;
                     arg: varargs[pointer,dataAddr]) =
  var pp = pointer p
  var gridDim, blockDim: CudaDim3
  gridDim.x = blocksPerGrid
  gridDim.y = 1
  gridDim.z = 1
  blockDim.x = threadsPerBlock
  blockDim.y = 1
  blockDim.z = 1
  var args: array[arg.len, pointer]
  for i in 0..<arg.len: args[i] = arg[i]
  #echo "really launching kernel"
  let err = cudaLaunchKernel(pp, gridDim, blockDim, addr args[0])
  if err:
    echo err
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
template onGpuNowait*(body: untyped): auto =
  onGpuNoWait(gpuNumThreadsRequest, gpuBlockSizeRequest, body)
template onGpuNowait*(n,body: untyped): auto =
  var b = gpuBlockSizeRequest
  while b > n: b = b div 2
  onGpuNoWait(n, b, body)
#template onGpuNowait*(n,b,body: untyped): auto =
#  onGpuNoWait(n, b, body)

template onGpu*(body: untyped) =
  let finalize = onGpuNoWait(body)
  finalize()
template onGpu*(n,body: untyped) =
  let finalize = onGpuNoWait(n, body)
  finalize()
template onGpu*(n,b,body: untyped) =
  let finalize = onGpuNoWait(n, b, body)
  finalize()

when isMainModule:
  type FltArr = UncheckedArray[float32]

  proc vectorAdd(A: FltArr; B: FltArr; C: var FltArr; n: int32) {.cudaGlobal.} =
    var i = blockDim.x * blockIdx.x + threadIdx.x
    if i < n:
      C[i] = A[i] + B[i]

  proc test =
    var n = 50000.cint
    var
      a = newSeq[float32](n)
      b = newSeq[float32](n)
      c = newSeq[float32](n)
    var threadsPerBlock: cint = 256
    var blocksPerGrid: cint = (n + threadsPerBlock - 1) div threadsPerBlock

    cudaLaunch(vectorAdd, blocksPerGrid, threadsPerBlock, a, b, c, n)

    template getGpuPtr(x: int): untyped = x
    template getGpuPtr[T](x: seq[T]): untyped = addr(x[0])
    template `[]`(x: ptr SomeNumber, i: SomeInteger): untyped {.used.} =
      cast[ptr UncheckedArray[type(x[])]](x)[][i]
    template `[]=`(x: ptr SomeNumber, i: SomeInteger, y:untyped): untyped {.used.} =
      cast[ptr UncheckedArray[type(x[])]](x)[][i] = y

    onGpu(n):
      let i = getBlockDim().x * getBlockIdx().x + getThreadIdx().x
      if i < n:
        c[i] = a[i] + b[i]

  test()
