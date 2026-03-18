import base/threading
import backend/expr
import macros

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

macro echoTyped(body: auto): auto =
  echo body.repr
  #echo body.treerepr
  result = body

proc gpuDefaultNumThreads*(): int =
  var nt = 0
  threads:
    threadSingle:
      nt = numThreads
  result = nt

macro onGpuNowait*(n,b,body: untyped): auto =
  proc deref(x,g,i:NimNode):auto = newCall("getGpu",x,g)
  template target(cpuPrepare, cpuFinalize, body: untyped) =
    mixin toGpu, getGpu, fromGpu
    cpuPrepare  # a let section declare and save device pointers
    proc gpuProc {.gensym.} =
      if numThreads == 1:
        threads:
          body
      else:
        body
    gpuProc()
    proc finalize {.gensym.} =
      threadBarrier()
      cpuFinalize
      #threadBarrier()
    finalize
  let
    v = prepareVars(body, deref)  # gather gpu pointers in symbols, body is changed accordingly
    cpuPrepare = genCpuPrepare v
    cpuFinalize = genCpuFinalize v
  result = getast(target(cpuPrepare, cpuFinalize, body))
  case dumpKernels
  of 1:
    echo result.repr
  of 2:
    echo result.treerepr
  else:
    if dumpKernels > 2:
      result = newCall(bindsym"echoTyped", result)

var gpuNumThreadsRequest* = 0
var gpuBlockSizeRequest* = 0
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

#template onGpu*(n,x:untyped) = onGpu(x)
#template onGpu*(n,t,x:untyped) = onGpu(x)
