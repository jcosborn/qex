import base/[threading, profile]
import backend/expr
import base/metaUtils
import macros, strutils

var kernelCallCount* = 0  # _twice_ kernel calls, incremented before setup and after finalize
const dumpKernels {.intdefine.} = 0

template gpuMalloc*(size:SomeInteger):pointer = alloc(size)
template gpuFree*(device_ptr:pointer) =
  ## nil is a no-op, as with the device frees; dealloc(nil) crashes
  let p = device_ptr
  if p != nil: dealloc(p)
proc gpuMemCpyToCPU*(dst: pointer, src: pointer; length: SomeInteger) =
  copymem(dst, src, length)
proc gpuMemCpyToGPU*(dst: pointer, src: pointer; length: SomeInteger) =
  copymem(dst, src, length)
template gpuMemCpyToCpuAsync*(dst, src: pointer; length: SomeInteger) =
  ## the copy to the host after the kernels already submitted, complete here
  copyMem(dst, src, length)
const gpuCopyFast* = high(int)  ## bytes up to which a queued copy has its low cost

template gpuNumThreads*:auto = numThreads
template gpuThreadNum*:auto = threadNum

# The kernels of the GPU backends run here as loops on the calling host
# thread, complete on return, the reference for the device code

template gpuFor*(i: untyped; n: SomeInteger; sub: untyped; body: untyped) =
  ## body for each i in 0..<n in order, n evaluated once; sub, the
  ## sub-group size of the GPU backends, does not apply
  let gpuN = int(n)
  for i in 0..<gpuN:
    body
template gpuFor*(i: untyped; n: SomeInteger; body: untyped) = gpuFor(i, n, 0, body)
template gpuForAsync*(i: untyped; n: SomeInteger; sub: untyped; body: untyped) =
  ## gpuFor, complete on return as well
  gpuFor(i, n, sub, body)
template gpuForAsync*(i: untyped; n: SomeInteger; body: untyped) = gpuFor(i, n, 0, body)
template gpuWaitAsync* = discard
template gpuAtomicAdd*(r: ptr UncheckedArray[float]; k: int; v: float) =
  ## r[k] += v, in kernels, which one thread runs
  r[k] += v
template gpuAtomicAdd*(r: ptr UncheckedArray[uint64]; k: int; v: uint64) =
  ## r[k] += v modulo 2^64, in kernels
  r[k] += v
template gpuMallocHost*(size: SomeInteger): pointer = alloc(size)
template gpuFreeHost*(p: pointer) = dealloc(p)

proc genCpuPrepare(n:seq[NimNode]):NimNode =
  mixin toGpu
  template r(x,v:untyped):untyped =
    var v = toGpu(x)
  result = newstmtlist()
  for c in n:
    result.add getast r(c[0],c[1])

proc gpuFinalizer*[P](fl: static string; count: int; ps: P): proc() =
  ## The finalizer of a kernel launch: copies what the kernel wrote back to the
  ## host variables of ps (genHostPairs).  The closure has an environment of its
  ## own, so every thread of a threads block gets one and may keep it.
  mixin fromGpu
  result = proc() =
    tic(fl)
    var countSave = 0
    threadBarrier()
    threadSingle:
      countSave = kernelCallCount
      kernelCallCount = count
    for p in fields(ps): hostFromGpu(p[0], p[1])
    threadSingle:
      kernelCallCount = countSave
    #threadBarrier()
    toc("cpuFinalize")

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
  template target(fl, cpuPrepare, hostPairs, body: untyped) =
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
      inc kernelCallCount  # increment after launch
      gpuFinalizer(fl, thisKernelCallCount, hostPairs)
  let
    v = prepareVars(body, deref)  # gather gpu pointers in symbols, body is changed accordingly
    cpuPrepare = genCpuPrepare v
    hostPairs = genHostPairs(v, proc(c: NimNode): NimNode = c[1])
  result = getast(target(fl,cpuPrepare, hostPairs, body))
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
