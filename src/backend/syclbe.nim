import macros, strutils
#import base/metaUtils
import backend/expr
import base/[metaUtils, threading, profile]
import sycl

var kernelCallCount* = 0  # count of kernel calls
const dumpKernels {.intdefine.} = 0

#let dev = defaultDevice()
#let q = dev.queue
var platform: Platform  # default platform
var dev: Device
var q: Queue

proc gpuNumDevices*: int =
  result = int platform.getDevices.size()

proc gpuDeviceName*: string =
  $dev.name

proc gpuInit*(device: int) =
  let devs = platform.getDevices()
  let n = int devs.size()
  let d = device mod n
  dev = devs[d]
  q = dev.inOrderQueue

template gpuMalloc*(size: SomeInteger): pointer = mallocDevice(size, q)
template gpuMalloc[T](x: var ptr UncheckedArray[T], n: int) =
  x = cast[typeof x](gpuMalloc(n*sizeof(T)))
template gpuMalloc[T](x: ptr T) =
  x = cast[typeof x](gpuMalloc(sizeof(T)))
template gpuFree*(device_ptr: pointer) = freeDevice(device_ptr, q)
template gpuMemCpyToCPU*(dst: pointer, src: pointer; length: SomeInteger) =
  memcpy(q, dst, src, length)
template gpuMemCpyToGPU*(dst: pointer, src: pointer; length: SomeInteger) =
  memcpy(q, dst, src, length)
template gpuMemCpyToCpuAsync*(dst, src: pointer; length: SomeInteger) =
  ## Queues the copy to the host after the kernels already submitted (the
  ## queue is in order); gpuWaitAsync waits for it too.
  memcpyAsync(q, dst, src, length)
const gpuCopyFast* = high(int)  ## bytes up to which a queued copy has its low cost

{.emit: """/*INCLUDESECTION*/
#include <level_zero/ze_api.h>
#include <sycl/ext/oneapi/backend/level_zero.hpp>
""".}

const gpuThreads {.intdefine.} = 128  ## work group size of gpuFor kernels
var gpuIt {.importc, nodecl.}: Nd1  # the work item of the gpuFor kernels

template gpuForAsync*(i: untyped; n: SomeInteger; sub: untyped; body: untyped) =
  ## One kernel over i in 0..<n, returning before it completes; gpuWaitAsync
  ## waits for it.  Sub-groups of sub threads, or as the compiler chooses for
  ## sub = 0.  Captured pointers must be device pointers.
  const gpuA = when sub > 0: " [[sycl::reqd_sub_group_size(" & $sub & ")]]" else: ""
  let gpuN = int(n)
  if gpuN > 0:
    let gpuM = csize_t((gpuN + gpuThreads - 1) div gpuThreads * gpuThreads)
    {.emit: [q, ".parallel_for(sycl::nd_range<1>(sycl::range<1>(", gpuM, "), sycl::range<1>(", gpuThreads, ")), [=](sycl::nd_item<1> gpuIt)", gpuA, " {"].}
    block:
      let i = int(gpuIt[])
      if i < gpuN:
        body
    {.emit: "});".}

template gpuForAsync*(i: untyped; n: SomeInteger; body: untyped) = gpuForAsync(i, n, 0, body)

template gpuFor*(i: untyped; n: SomeInteger; sub: untyped; body: untyped) =
  ## gpuForAsync waiting for the kernel
  gpuForAsync(i, n, sub, body)
  q.wait
template gpuFor*(i: untyped; n: SomeInteger; body: untyped) = gpuFor(i, n, 0, body)

template gpuWaitAsync* = q.wait

template gpuMallocHost*(size: SomeInteger): pointer =
  ## pinned host memory, for fast copies from the device
  mallocHost(size, q)
template gpuFreeHost*(p: pointer) = freeHost(p, q)

template gpuAtomicAdd*(r: ptr UncheckedArray[float]; k: int; v: float) =
  ## r[k] += v atomically, in kernels
  let kk = k
  let vv = v
  {.emit: ["sycl::atomic_ref<double, sycl::memory_order::relaxed, sycl::memory_scope::device, sycl::access::address_space::global_space>(", r, "[", kk, "]).fetch_add(", vv, ");"].}

template gpuAtomicAdd*(r: ptr UncheckedArray[uint64]; k: int; v: uint64) =
  ## r[k] += v modulo 2^64 atomically, in kernels
  let kk = k
  let vv = v
  {.emit: ["sycl::atomic_ref<unsigned long long, sycl::memory_order::relaxed, sycl::memory_scope::device, sycl::access::address_space::global_space>(*(unsigned long long*)&", r, "[", kk, "]).fetch_add((unsigned long long)", vv, ");"].}

proc gpuZeContext*(): tuple[ctx, dev: pointer] =
  ## Level Zero context and device of the queue
  var ctx, dev: pointer
  {.emit: [ctx, " = sycl::get_native<sycl::backend::ext_oneapi_level_zero>(", q, ".get_context());\n",
           dev, " = sycl::get_native<sycl::backend::ext_oneapi_level_zero>(", q, ".get_device());"].}
  (ctx, dev)

proc gpuFinish*[P](fl: static string; q: ptr Queue; flops: float; ps: P; back: proc(ps: P) {.nimcall.}) =
  ## The end of a kernel launch: waits for the queue, then back copies what
  ## the kernel wrote to the host variables of ps (hostBack).
  tic(fl)
  threadSingle:
    q[].wait
  toc("wait",flops=flops)
  back(ps)
  toc("cpuFinalize")

proc gpuFinalizer*[P](fl: static string; q: ptr Queue; flops: float; ps: P; back: proc(ps: P) {.nimcall.}): proc() =
  ## gpuFinish as a closure, for onGpuNowait.  The closure has an
  ## environment of its own, so every thread of a threads block gets one
  ## and may keep it.
  result = proc() = gpuFinish(fl, q, flops, ps, back)

proc gpuDefaultNumThreads*(): int =
  32 * q.device.maxComputeUnits.int * q.device.subGroupSizes.max_element.int
  #64 * q.device.maxComputeUnits.int * q.device.subGroupSizes.max_element.int

macro gpuKernel(fin: static int; n: untyped; sub: static int; body: untyped): untyped =
  ## Queues body as a kernel of n threads, n a count or what gpuSites turns
  ## into one.  The kernel is that of gpuForAsync, its captures the device
  ## values of the variables that body names (backend/expr), locals of the
  ## launch that the kernel reads in place.  fin: 0 queues the kernel, 1
  ## gives the finalizer, which waits for the queue and copies back, 2 does
  ## that at once (a wait only where nothing is copied back).
  let li = body.lineinfo
  let lis = li.split({'/','.','(',','})
  let fl = lis[^4] & "(" & lis[^2] & ")"
  template target(fl, n, sub, pre, shadow, hostPairs, queued, copies, body, fin: untyped) =
    mixin toGpu, getGpu, fromGpu, gpuSites
    block:
      tic(fl)
      inc kernelCallCount
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
      # one profile timer per launch: "launch", or "wait" where onGpu only waits
      when fin == 2 and not copies:
        if numThreads > 1:
          threadSingle: gpuWaitAsync()
        else: gpuWaitAsync()
        toc("wait", flops=(when declared gpuWaitFlops: float(gpuWaitFlops) else: 0.0))
      else:
        toc("launch")
        when fin == 1: gpuFinalizer(fl, addr q, (when declared gpuWaitFlops: float(gpuWaitFlops) else: 0.0), hp, hostBack(hp))
        elif fin == 2: gpuFinish(fl, addr q, (when declared gpuWaitFlops: float(gpuWaitFlops) else: 0.0), hp, hostBack(hp))
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

template onGpuNowait*(n, body: untyped): auto = gpuKernel(1, n, 0, body)
template onGpuNowait*(n, sub, body: untyped): auto =
  ## sub: the size of the sub-groups on Intel GPUs, as the compiler chooses for 0
  gpuKernel(1, n, sub, body)
template onGpuNowait*(body: untyped): auto = gpuKernel(1, gpuDefaultNumThreads(), 0, body)
template onGpuQueue*(n, body: untyped) =
  ## onGpuNowait for a kernel whose variables stay on the device
  ## (gpuResident): nothing is copied back, and gpuWaitAsync waits for it.
  gpuKernel(0, n, 0, body)
template onGpuQueue*(n, sub, body: untyped) = gpuKernel(0, n, sub, body)

template onGpu*(body: untyped) = gpuKernel(2, gpuDefaultNumThreads(), 0, body)
template onGpu*(n, body: untyped) = gpuKernel(2, n, 0, body)
template onGpu*(n, sub, body: untyped) = gpuKernel(2, n, sub, body)
