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
  let dev = device mod gpuNumDevices()
  let err = hipSetDevice(cint dev)
  if err:
    echo "HIP error: ", cast[cint](err)
  doAssert(not err)

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

proc gpuFinalizer*[P](fl: static string; count: int; flops: float; ps: P; back: proc(ps: P) {.nimcall.}): proc() =
  ## The finalizer of a kernel launch: waits for the device, then back copies
  ## what the kernel wrote to the host variables of ps (genHostPairs,
  ## hostBack).  The closure has an environment of its own, so every thread of
  ## a threads block gets one and may keep it.
  result = proc() =
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

macro gpuKernel(fin: static bool; n, body: untyped): untyped =
  ## Queues body as a kernel of n threads, n a count or what gpuSites turns
  ## into one.  The kernel is that of gpuForAsync, its captures the device
  ## values of the variables that body names (backend/expr).  With fin the
  ## result is the finalizer, which waits for the device and copies back.
  let li = body.lineinfo
  let lis = li.split({'/','.','(',','})
  let fl = lis[^4] & "(" & lis[^2] & ")"
  template target(fl, n, pre, shadow, hostPairs, queued, body, fin: untyped): untyped =
    mixin toGpu, getGpu, fromGpu, gpuSites
    inc kernelCallCount  # increment before setup
    block:
      let thisKernelCallCount = kernelCallCount  # save current value for finalizer
      tic(fl)
      let gpuNt {.inject.} = int(gpuSites(n))
      gpuLaunchThreads = gpuNt
      pre  # the device values
      toc("cpuPrepare")
      threadSingle:
        gpuForAsync(gpuTn, gpuNt):
          const inOnGpu {.inject,used.} = true
          shadow
          inlineProcs:
            body
      toc("launch")
      inc kernelCallCount  # increment after launch
      when fin:
        let hp = hostPairs
        gpuFinalizer(fl, thisKernelCallCount, (when declared gpuWaitFlops: float(gpuWaitFlops) else: 0.0), hp, hostBack(hp))
      else:
        queued
  let c = kernelCap(body)
  result = getast(target(fl, n, c.pre, c.shadow, c.pairs, c.queued, c.body, newLit(fin)))
  case dumpKernels
  of 1:
    echo li
    echo result.repr
  of 2:
    echo li
    echo result.treerepr
  else: discard

var gpuNumThreadsRequest* = 32*1024  ## threads of a kernel launched without a count
template onGpuNowait*(n, body: untyped): auto = gpuKernel(true, n, body)
template onGpuNowait*(body: untyped): auto = gpuKernel(true, gpuNumThreadsRequest, body)
template onGpuQueue*(n, body: untyped) =
  ## onGpuNowait for a kernel whose variables stay on the device
  ## (gpuResident): nothing is copied back, and gpuWaitAsync waits for it.
  gpuKernel(false, n, body)

template onGpu*(body: untyped) =
  let finalize = onGpuNowait(body)
  finalize()
template onGpu*(n,body: untyped) =
  let finalize = onGpuNowait(n, body)
  finalize()
