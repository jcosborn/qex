import base/threading
import backend/expr
import base/metaUtils
import macros, strutils

var kernelCallCount* = 0  # _twice_ kernel calls, incremented before setup and after finalize
const dumpKernels {.intdefine.} = 0

template gpuMalloc*(size:SomeInteger):pointer = alloc(size)
template gpuFree*(device_ptr:pointer) = dealloc(device_ptr)
proc gpuMemCpyToCPU*(dst: pointer, src: pointer; length: SomeInteger) =
  copymem(dst, src, length)
proc gpuMemCpyToGPU*(dst: pointer, src: pointer; length: SomeInteger) =
  copymem(dst, src, length)

template gpuNumThreads*:auto = numThreads
template gpuThreadNum*:auto = threadNum

proc genCpuPrepare(n:seq[NimNode]):NimNode =
  mixin toGpu
  template r(x,v:untyped):untyped =
    var v = toGpu(x)
  result = newstmtlist()
  for c in n:
    result.add getast r(c[0],c[1])

proc genCpuFinalize(n:seq[NimNode]):NimNode =
  template r(x,v:untyped):untyped =
    fromGpu(x,v)
  result = newstmtlist()
  for c in n:
    result.add getast r(c[0],c[1])

#macro echoTyped(body: auto): auto =
#  echo body.repr
#  #echo body.treerepr
#  result = body

proc gpuDefaultNumThreads*(): int =
  var nt = 0
  threads:
    threadSingle:
      nt = numThreads
  result = nt

macro onGpuNowait*(n,b,body: untyped): auto =
  let li = body.lineinfo
  let lis = li.split({'/','.','(',','})
  let fl = lis[^4] & "(" & lis[^2] & ")"
  proc deref(x,g,i:NimNode):auto = newCall("getGpu",x,g)
  template target(fl, cpuPrepare, cpuFinalize, body: untyped) =
    mixin toGpu, getGpu, fromGpu
    inc kernelCallCount  # increment before setup
    block:
      tic(fl)
      let thisKernelCallCount = kernelCallCount  # save current value for finalizer
      cpuPrepare  # a let section declare and save device pointers
      toc("cpuPrepare")
      #proc gpuProc {.gensym.} =
      block:
        const inOnGpu {.inject,used.} = true
        if numThreads == 1:
          threads:
            body
        else:
          body
      #gpuProc()
      when declared gpuWaitFlops:
        toc("wait",flops=gpuWaitFlops)
      else:
        toc("wait")
      proc finalize {.gensym.} =
        tic(fl)
        var countSave = 0
        threadBarrier()
        threadSingle:
          countSave = kernelCallCount
          kernelCallCount = thisKernelCallCount
        cpuFinalize
        threadSingle:
          kernelCallCount = countSave
        #threadBarrier()
        toc("cpuFinalize")
      inc kernelCallCount  # increment after launch
      finalize
  let
    v = prepareVars(body, deref)  # gather gpu pointers in symbols, body is changed accordingly
    cpuPrepare = genCpuPrepare v
    cpuFinalize = genCpuFinalize v
  result = getast(target(fl,cpuPrepare, cpuFinalize, body))
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

var gpuNumThreadsRequest* = 0
var gpuBlockSizeRequest* = 0
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

#template onGpu*(n,x:untyped) = onGpu(x)
#template onGpu*(n,t,x:untyped) = onGpu(x)

type GpuSum*[T] = object
    val: T
proc newGpuSum*[T](n: int): GpuSum[T] =
  discard
template value*(x: GpuSum): auto = x.val
template toGpu*(x: GpuSum): auto = addr x
template getGpu*(x: GpuSum, g: ptr GpuSum): auto = g[]
template fromGpu*(x: GpuSum, g: ptr GpuSum): auto = discard

proc reduce*[T](gs: var GpuSum[T], x: T) =
  var t = x
  t.threadSum
  gs.val = t
