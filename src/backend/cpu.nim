import base/[threading, profile]
import backend/expr
import base/metaUtils
import macros, strutils

var kernelCallCount* = 0  # _twice_ kernel calls, incremented before setup and after finalize
const dumpKernels {.intdefine.} = 0

proc c_aligned_alloc(a, n: csize_t): pointer {.importc: "aligned_alloc", header: "<stdlib.h>".}
proc c_free(p: pointer) {.importc: "free", header: "<stdlib.h>".}

proc gpuMalloc*(size: SomeInteger): pointer =
  ## host memory on 64 bytes, as kernels load whole SIMD sites
  c_aligned_alloc(64, csize_t((int(size) + 63) and not 63))
template gpuFree*(device_ptr: pointer) = c_free(device_ptr)  ## nil is a no-op
proc gpuMemCpyToCPU*(dst: pointer, src: pointer; length: SomeInteger) =
  copymem(dst, src, length)
proc gpuMemCpyToGPU*(dst: pointer, src: pointer; length: SomeInteger) =
  copymem(dst, src, length)
template gpuMemCpyToCpuAsync*(dst, src: pointer; length: SomeInteger) =
  ## the copy to the host after the kernels already submitted, complete here
  copyMem(dst, src, length)
const gpuCopyFast* = high(int)  ## bytes up to which a queued copy has its low cost

# The kernels of the GPU backends run here as loops over the host threads,
# complete on return

template gpuFor*(i: untyped; n: SomeInteger; sub: untyped; body: untyped) =
  ## body for each i in 0..<n, n evaluated once.  One thread calls it: the
  ## host threads share the loop, or the calling thread runs it where the
  ## call is inside a threads block.  sub, the sub-group size of the GPU
  ## backends, does not apply.
  let gpuN = int(n)
  when declared(qexThreadsBlock):
    # no nested proc here: it would move the locals of every thread that
    # body names into one shared environment
    for i in 0..<gpuN:
      body
  else:
    proc gpuLoop(gpuLo, gpuHi: int) {.gensym.} =
      for i in gpuLo..<gpuHi:
        body
    if numThreads == 1:
      threads:
        gpuLoop((gpuN*threadNum) div numThreads, (gpuN*(threadNum+1)) div numThreads)
    else:
      gpuLoop(0, gpuN)
template gpuFor*(i: untyped; n: SomeInteger; body: untyped) = gpuFor(i, n, 0, body)
template gpuForAsync*(i: untyped; n: SomeInteger; sub: untyped; body: untyped) =
  ## gpuFor, complete on return as well
  gpuFor(i, n, sub, body)
template gpuForAsync*(i: untyped; n: SomeInteger; body: untyped) = gpuFor(i, n, 0, body)
template gpuWaitAsync* = discard
template gpuAtomicAdd*(r: ptr UncheckedArray[float]; k: int; v: float) =
  ## r[k] += v atomically, in kernels
  let kk = k
  let vv = v
  {.emit: ["#pragma omp atomic update\n", r, "[", kk, "] += ", vv, ";"].}
template gpuAtomicAdd*(r: ptr UncheckedArray[uint64]; k: int; v: uint64) =
  ## r[k] += v modulo 2^64 atomically, in kernels
  let kk = k
  let vv = v
  {.emit: ["#pragma omp atomic update\n", r, "[", kk, "] += ", vv, ";"].}
template gpuMallocHost*(size: SomeInteger): pointer = gpuMalloc(size)
template gpuFreeHost*(p: pointer) = c_free(p)

proc gpuFinish*[P](fl: static string; count: int; ps: P; back: proc(ps: P) {.nimcall.}) =
  ## The end of a kernel launch: meets the team, then back copies what
  ## the kernel wrote to the host variables of ps (hostBack).
  tic(fl)
  var countSave = -1  # the count before, in the thread that lowers it
  threadBarrier()
  threadSingle:
    countSave = kernelCallCount
    kernelCallCount = count
  back(ps)
  if countSave >= 0: kernelCallCount = countSave  # the same thread restores it
  threadBarrier()
  toc("cpuFinalize")

proc gpuFinalizer*[P](fl: static string; count: int; ps: P; back: proc(ps: P) {.nimcall.}): proc() =
  ## gpuFinish as a closure, for onGpuNowait.  The closure has an
  ## environment of its own, so every thread of a threads block gets one
  ## and may keep it.
  result = proc() = gpuFinish(fl, count, ps, back)

proc gpuDefaultNumThreads*(): int =
  var nt = 0
  threads:
    threadSingle:
      nt = numThreads
  result = nt

macro gpuKernel(fin: static int; n: untyped; sub: static int; body: untyped): untyped =
  ## Runs body as a kernel in the host threads: outside a threads block the
  ## kernel opens one, inside one every thread of the team calls it and
  ## the team meets at a barrier after the kernel, so that each thread finds
  ## it complete on return.  n, the thread count of the devices, does not
  ## apply.  fin: 0 copies nothing back, 1 gives the finalizer, which copies
  ## back, 2 copies back at once where there is something to copy.
  let li = body.lineinfo
  let lis = li.split({'/','.','(',','})
  let fl = lis[^4] & "(" & lis[^2] & ")"
  template target(fl, pre, shadow, hostPairs, queued, copies, body, fin: untyped) =
    mixin toGpu, getGpu, fromGpu
    inc kernelCallCount  # increment before setup
    block:
      tic(fl)
      let thisKernelCallCount = kernelCallCount  # save current value for finalizer
      pre  # the device values
      block:
        const inOnGpu {.inject,used.} = true
        when declared(qexThreadsBlock):
          # no nested proc here: it would move the device values of every
          # thread into one shared environment
          let gpuTn {.inject, used.} = threadNum
          let gpuNt {.inject, used.} = numThreads
          shadow
          body
        else:
          proc gpuKern() {.gensym.} =
            let gpuTn {.inject, used.} = threadNum
            let gpuNt {.inject, used.} = numThreads
            shadow
            body
          if numThreads == 1:
            threads:
              gpuKern()
          else:
            gpuKern()
      when declared gpuWaitFlops:
        toc("wait",flops=gpuWaitFlops)
      else:
        toc("wait")
      inc kernelCallCount  # increment after launch
      when fin == 1 or (fin == 2 and copies):
        let hp = hostPairs
      when fin == 1: gpuFinalizer(fl, thisKernelCallCount, hp, hostBack(hp))
      elif fin == 2 and copies: gpuFinish(fl, thisKernelCallCount, hp, hostBack(hp))
      else:
        when fin == 0: queued
        if numThreads > 1: threadBarrier()
  let c = kernelCap(body)
  result = getast(target(fl, c.pre, c.shadow, c.pairs, c.queued, c.copies, c.body, newLit(fin)))
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
template onGpuNowait*(body: untyped): auto = gpuKernel(1, 0, 0, body)
template onGpuQueue*(n, body: untyped) =
  ## onGpu for a kernel whose variables stay on the device (gpuResident):
  ## nothing is copied back.  Complete on return here; gpuWaitAsync waits
  ## for it with a GPU backend.
  gpuKernel(0, n, 0, body)
template onGpuQueue*(n, sub, body: untyped) = gpuKernel(0, n, sub, body)

template onGpu*(body: untyped) = gpuKernel(2, 0, 0, body)
template onGpu*(n, body: untyped) = gpuKernel(2, n, 0, body)
template onGpu*(n, sub, body: untyped) = gpuKernel(2, n, sub, body)

type GpuSum*[T] = object
    ## The sum of the values the threads of a kernel hand to reduce, each
    ## once: here the sum over the host threads.
    val: T
proc newGpuSum*[T](n: int): GpuSum[T] =
  discard
template free*(x: GpuSum) = discard
template value*(x: GpuSum): auto = x.val
template toGpu*(x: GpuSum): auto = addr x
template getGpu*(x: GpuSum, g: ptr GpuSum): auto = g
template fromGpu*(x: GpuSum, g: ptr GpuSum): auto = discard
template gpuResident*[T](x: typedesc[GpuSum[T]]): bool = true

proc reduce*[T](gs: ptr GpuSum[T], x: T) =
  var t = x
  t.threadSum
  gs.val = t
