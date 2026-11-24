import macros, strutils
import base/[metaUtils,profile,threading]
import expr
import hip
export hip.dataAddr

var kernelCallCount* = 0  # _twice_ kernel calls, incremented before setup and after finalize
#var kernelSetupCount* = 0  # number of kernel calls, incremented before setup
#var kernelFinalizeCount* = 0  # number of kernel calls, incremented after finalize
const dumpKernels {.intdefine.} = 0

proc gpuNumDevices*: int =
  var deviceCount = cint 0
  discard hipGetDeviceCount(addr deviceCount)
  result = deviceCount

proc gpuInit*(device: int) =
  ## Selects the device of this rank in every host thread.  HIP keeps the
  ## current device per thread, device 0 until set, and the thread that takes
  ## the single of a launch in a threads block launches the kernel.
  let dev = cint(device mod gpuNumDevices())
  var bad = cint 0
  threads:
    let err = hipSetDevice(dev)
    if err:
      threadCritical: bad = cast[cint](err)
  if bad != 0:
    echo "HIP error: ", bad
  doAssert(bad == 0)

template gpuMalloc*(size: SomeInteger):pointer =
  var p:pointer
  let err = hipMalloc(p, csize_t size)
  if err:
    echo "gpuMalloc: ", $err  # $ of this module: callers need not import hip
    p = cast[pointer](0)
  p
template gpuMalloc[T](x: var ptr UncheckedArray[T], n: int) =
  x = cast[typeof x](gpuMalloc(n*sizeof(T)))
template gpuMalloc[T](x: ptr T) =
  x = cast[typeof x](gpuMalloc(sizeof(T)))

template gpuFree*(p:pointer) =
  let err = hipFree(p)
  if err:
    echo $err
    quit cast[cint](err)

proc gpuMemset*(devPtr: pointer, value: SomeInteger, count: SomeInteger) =
  let err = hipMemset(devPtr, cint value, csize_t count)
  if err:
    echo "gpuMemset: ", $err
    quit cast[cint](err)
#proc gpuMemCpyToGpu*(dst,src: pointer, count: SomeInteger):cint {.discardable.} =
proc gpuMemCpyToGpu*(dst,src: pointer, count: SomeInteger) =
  let err = hipMemcpy(dst,src,csize_t count,hipMemcpyHostToDevice)
  if err:
    echo instantiationInfo()
    echo "  gpuMemCpyToGpu: ", $err
    quit cast[cint](err)
#proc gpuMemCpyToCpu*(dst,src: pointer, count: SomeInteger):cint {.discardable.} =
template gpuMemCpyToCpu*(dst,src: pointer, count: SomeInteger) =
  let err = hipMemcpy(dst,src,csize_t count,hipMemcpyDeviceToHost)
  if err:
    echo instantiationInfo()
    echo "  gpuMemCpyToCpu: ", $err
    quit cast[cint](err)

template gpuMemCpyToCpuAsync*(dst, src: pointer; count: SomeInteger) =
  ## Queues the copy to the pinned host memory dst after the kernels already
  ## submitted; gpuWaitAsync waits for it too.
  let err = hipMemcpyAsync(dst, src, csize_t count, hipMemcpyDeviceToHost, nil)
  if err:
    echo "gpuMemCpyToCpuAsync: ", $err
    quit cast[cint](err)

const gpuCopyFast* = 16384  ## bytes up to which a queued copy has its low cost: 2.9 us on MI300X, 8.2 us from 16400 bytes on

const gpuThreads {.intdefine.} = 128  ## threads per block of gpuFor kernels
var gpuI {.importc, nodecl.}: clong  # the index in the gpuFor kernels

template gpuForAsync*(i: untyped; n: SomeInteger; sub: untyped; body: untyped) =
  ## One kernel over i in 0..<n, returning before it completes; gpuWaitAsync
  ## waits for it.  sub, the sub-group size on Intel GPUs, does not apply.
  ## Captured pointers must be device pointers.
  let gpuN = int(n)
  if gpuN > 0:
    let gpuB = cuint((gpuN + gpuThreads - 1) div gpuThreads)
    {.emit: ["qexFor<", gpuThreads, "><<<", gpuB, ", ", gpuThreads, ">>>(", gpuN, ", [=] __device__ (long gpuI) {"].}
    # Nim zeroes and copies arrays with its host procs nimZeroMem and nimCopyMem
    {.emit: "\n#define nimZeroMem(b,len) memset((b),0,(len))\n#define nimCopyMem(a,b,len) memcpy((a),(b),(len))".}
    block:
      let i = int(gpuI)
      body
    {.emit: "\n#undef nimZeroMem\n#undef nimCopyMem\n});".}
    let err = hipGetLastError()
    if err:
      echo "gpuForAsync: ", $err
      quit 1
template gpuForAsync*(i: untyped; n: SomeInteger; body: untyped) = gpuForAsync(i, n, 0, body)

template gpuWaitAsync* =
  let err = hipDeviceSynchronize()
  if err:
    echo "gpuWaitAsync: ", $err
    quit 1

template gpuFor*(i: untyped; n: SomeInteger; sub: untyped; body: untyped) =
  ## gpuForAsync waiting for the kernel
  gpuForAsync(i, n, sub, body)
  gpuWaitAsync()
template gpuFor*(i: untyped; n: SomeInteger; body: untyped) = gpuFor(i, n, 0, body)

template gpuAtomicAdd*(r: ptr UncheckedArray[float]; k: int; v: float) =
  ## r[k] += v atomically, in kernels
  let kk = k
  let vv = v
  {.emit: ["atomicAdd(&", r, "[", kk, "], ", vv, ");"].}

template gpuAtomicAdd*(r: ptr UncheckedArray[uint64]; k: int; v: uint64) =
  ## r[k] += v modulo 2^64 atomically, in kernels
  let kk = k
  let vv = v
  {.emit: ["atomicAdd((unsigned long long*)&", r, "[", kk, "], (unsigned long long)", vv, ");"].}

proc gpuMallocHost*(size: SomeInteger): pointer =
  ## pinned host memory, for fast copies from the device
  let err = hipHostMalloc(addr result, csize_t size)
  if err:
    echo "gpuMallocHost: ", $err
    quit 1
template gpuFreeHost*(p: pointer) = discard hipHostFree(p)

proc gpuFinish*[P](fl: static string; count: int; flops: float; ps: P; back: proc(ps: P) {.nimcall.}) =
  ## The end of a kernel launch: waits for the device, then back copies what
  ## the kernel wrote to the host variables of ps (hostBack).
  var countSave = -1  # the count before, in the thread that lowers it
  tic(fl)
  threadSingle:
    countSave = kernelCallCount
    kernelCallCount = count
    let err = hipDeviceSynchronize()
    if err:
      echo "onGpu kernel: ", $err
      quit cast[cint](err)
  toc("wait",flops=flops)
  back(ps)
  if countSave >= 0: kernelCallCount = countSave  # the same thread restores it
  threadBarrier()
  toc("cpuFinalize")

proc gpuFinalizer*[P](fl: static string; count: int; flops: float; ps: P; back: proc(ps: P) {.nimcall.}): proc() =
  ## gpuFinish as a closure, for onGpuNowait.  The closure has an
  ## environment of its own, so every thread of a threads block gets one
  ## and may keep it.
  result = proc() = gpuFinish(fl, count, flops, ps, back)

macro gpuKernel(fin: static int; n: untyped; sub: static int; body: untyped): untyped =
  ## Queues body as a kernel of n threads, n a count or what gpuSites turns
  ## into one.  The kernel is that of gpuForAsync, its captures the device
  ## values of the variables that body names (backend/expr), locals of the
  ## launch that the kernel reads in place.  fin: 0 queues the kernel, 1
  ## gives the finalizer, which waits for the device and copies back, 2 does
  ## that at once (a wait only where nothing is copied back).
  let li = body.lineinfo
  let lis = li.split({'/','.','(',','})
  let fl = lis[^4] & "(" & lis[^2] & ")"
  template target(fl, n, sub, pre, shadow, hostPairs, queued, copies, body, fin: untyped): untyped =
    mixin toGpu, getGpu, fromGpu, gpuSites
    inc kernelCallCount  # increment before setup
    block:
      let thisKernelCallCount = kernelCallCount  # save current value for finalizer
      tic(fl)
      let gpuNt {.inject.} = int(gpuSites(n))
      gpuLaunchThreads = gpuNt
      pre  # the device values
      when fin == 1 or (fin == 2 and copies):
        let hp = hostPairs  # before the kernel values take the names of the host variables
      shadow  # the kernel values: locals of the launch that the kernel captures and reads in place
      threadSingle:
        gpuForAsync(gpuTn, gpuNt, sub):
          const inOnGpu {.inject,used.} = true
          inlineProcs:
            body
      inc kernelCallCount  # increment after launch
      # one profile timer per launch: "launch", or "wait" where onGpu only waits
      when fin == 2 and not copies:
        if numThreads > 1:
          threadSingle: gpuWaitAsync()
        else: gpuWaitAsync()
        toc("wait", flops=(when declared gpuWaitFlops: float(gpuWaitFlops) else: 0.0))
      else:
        toc("launch")
        when fin == 1: gpuFinalizer(fl, thisKernelCallCount, (when declared gpuWaitFlops: float(gpuWaitFlops) else: 0.0), hp, hostBack(hp))
        elif fin == 2: gpuFinish(fl, thisKernelCallCount, (when declared gpuWaitFlops: float(gpuWaitFlops) else: 0.0), hp, hostBack(hp))
        else: queued
  let c = kernelCap(body, n)
  result = getast(target(fl, n, newLit(sub), c.pre, c.shadow, c.pairs, c.queued, c.copies, c.body, newLit(fin)))
  case dumpKernels
  of 1:
    echo li
    echo result.repr
  of 2:
    echo li
    echo result.treerepr
  else: discard

var gpuNumThreadsRequest* = 32*1024  ## threads of a kernel launched without a count
template onGpuNowait*(n, body: untyped): auto = gpuKernel(1, n, 0, body)
template onGpuNowait*(n, sub, body: untyped): auto =
  ## sub: the size of the sub-groups on Intel GPUs, as the compiler chooses for 0
  gpuKernel(1, n, sub, body)
template onGpuNowait*(body: untyped): auto = gpuKernel(1, gpuNumThreadsRequest, 0, body)
template onGpuQueue*(n, body: untyped) =
  ## onGpuNowait for a kernel whose variables stay on the device
  ## (gpuResident): nothing is copied back, and gpuWaitAsync waits for it.
  gpuKernel(0, n, 0, body)
template onGpuQueue*(n, sub, body: untyped) = gpuKernel(0, n, sub, body)

template onGpu*(body: untyped) = gpuKernel(2, gpuNumThreadsRequest, 0, body)
template onGpu*(n, body: untyped) = gpuKernel(2, n, 0, body)
template onGpu*(n, sub, body: untyped) = gpuKernel(2, n, sub, body)
