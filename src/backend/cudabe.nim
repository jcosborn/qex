import macros, strutils
import base/[metaUtils,profile]
import expr
import cuda
export cuda.dataAddr

var kernelCallCount* = 0  # _twice_ kernel calls, incremented before setup and after finalize
#var kernelSetupCount* = 0  # number of kernel calls, incremented before setup
#var kernelFinalizeCount* = 0  # number of kernel calls, incremented after finalize
const dumpKernels {.intdefine.} = 0

proc gpuNumDevices*: int =
  var deviceCount = cint 0
  discard cudaGetDeviceCount(addr deviceCount)
  result = deviceCount

proc gpuInit*(device: int) =
  let dev = device mod gpuNumDevices()
  let err = cudaSetDevice(cint dev)
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
    echo instantiationInfo()
    echo "  gpuMemCpyToGpu: ", err
#proc gpuMemCpyToCpu*(dst,src: pointer, count: SomeInteger):cint {.discardable.} =
template gpuMemCpyToCpu*(dst,src: pointer, count: SomeInteger) =
  let err = cudaMemcpy(dst,src,csize_t count,cudaMemcpyDeviceToHost)
  if err:
    echo instantiationInfo()
    echo "  gpuMemCpyToCpu: ", err

template gpuThreadNum*: auto =
  blockDim.x * blockIdx.x + threadIdx.x
template gpuNumThreads*: auto =
  gridDim.x * blockDim.x

template cudaDefs(body: untyped): untyped {.dirty.} =
  bind inlineProcs
  {.emit:["#define nimZeroMem(b,len) memset((b),0,(len))"].}
  {.emit:["#define nimCopyMem(a,b,len) memcpy((a),(b),(len))"].}
  {.pragma: shared, noInit, codegendecl:"__shared__ $# $#".}
  inlineProcs:
    body
  {.emit:["#undef nimZeroMem"].}
  {.emit:["#undef nimCopyMem"].}

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

macro onGpuNowait(nn0,tpb0: untyped, body: untyped): auto =
  let li = body.lineinfo
  let lis = li.split({'/','.','(',','})
  let fl = lis[^4] & "(" & lis[^2] & ")"
  template target(fl, nn, tpb, v, arg, cpuPrepare, cpuFinalize, body: untyped): untyped =
    mixin toGpu, getGpu, fromGpu
    inc kernelCallCount  # increment before setup
    block:
      let thisKernelCallCount = kernelCallCount  # save current value for finalizer
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
        var countSave = 0
        tic(fl)
        threadSingle:
          countSave = kernelCallCount
          kernelCallCount = thisKernelCallCount
          discard cudaDeviceSynchronize()
        toc("wait")
        cpuFinalize
        #threadBarrier()
        threadSingle:
          kernelCallCount = countSave
        toc("cpuFinalize")
      inc kernelCallCount  # increment after launch
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

    onGpu(n):
      let i = int getBlockDim().x * getBlockIdx().x + getThreadIdx().x
      if i < n:
        c[i] = a[i] + b[i]

    var errcnt = 0
    for i in 0..<n:
      let d = a[i] + b[i]
      if errcnt < 10 and c[i] != d:
        echo "error: ", i, "  ", c[i], "  ", d
        inc errcnt

  test()
