import macros, strutils
import base/[metaUtils,profile]
import hip, hipbe
export hipbe
import gpumem

#{.emit:"/*INCLUDESECTION*/#include <cooperative_groups.h>".}
proc importCG() {.header:"<cooperative_groups.h>",importcpp:"@".}
proc warpSumSmall*[T:HipTypes](x: T): T =  # for up to 32 bytes
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

template atomic_type*(t: HipTypes): typedesc = t
template atomic_type*[N:static int,T](t: typedesc[array[N,T]]): typedesc =
  mixin atomic_type
  type vt = vecType(t)
  when vt is HipTypes:
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
  #static: echo "type: ", $T, "  atomic_t: ", $atomic_t
  const n = sizeof(T) div sizeof(atomic_t)
  #doAssert(sizeof(T) == n * sizeof(atomic_t))
  var sum_tmp {.noInit.}: array[n, atomic_t]
  discard cmemcpy(addr sum_tmp, addr x, csize_t sizeof(T))
  for i in 0..<n:
    sum_tmp[i] = warpSumSmall(sum_tmp[i])
  discard cmemcpy(addr result, addr sum_tmp, csize_t sizeof(T))

template warpSum*[T](x: T): T =
  when sizeof(T) <= 32 and T is HipTypes:
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

type
  GpuSumDev*[T] = object
    count: cuint
    partial: UncheckedArray[T]
  GpuSumObj*[T] = object
    dev: ptr GpuSumDev[T]
    npartial: cuint
    #maxblock: int
    val: ptr T
  GpuSum*[T] = ref GpuSumObj[T]
# TODO: finalizer
proc newGpuSum*[T](ns: int): GpuSum[T] =
  result.new
  let n = (ns + 31) div 32  # divide by warp size
  #result.partial.gpuMalloc(n)
  pushGpuMemTag("GpuSumPartial")
  let bytes = sizeof(result.dev.count) + n*sizeof(T)
  let pgmdev = getGpuMem(addr result[], bytes)
  popGpuMemTag()
  result.dev = cast[typeof result.dev](pgmdev.p)
  result.dev.gpuMemset(0, bytes)
  result.npartial = cuint n
  #let err = hipMallocHost((ptr pointer)(addr result.val), csize_t sizeof(T))
  let err = hipHostMalloc((ptr pointer)(addr result.val), csize_t sizeof(T))
template value*(x: GpuSum | GpuSumObj): auto = x.val[]
template toGpu*(x: GpuSum): auto =
  let bytes = sizeof(x.dev.count) + x.npartial.int*sizeof(x.T)
  let pgmdev = getGpuMem(addr x[], bytes)
  x[]
template getGpu*(x: GpuSum, g: GpuSumObj): auto = g
template fromGpu*(x: GpuSum, g: GpuSumObj): auto = discard

proc reduce*[T](gs: GpuSumObj[T], x: T) =
  var isLastBlockDone {.shared.}: bool
  var aggregate = blockSum(x)
  if threadIdx.x == 0:
    if blockIdx.x < gs.npartial:
      gs.dev.partial[blockIdx.x] = aggregate;
    threadFence() # flush result
    # increment global block counter
    let value = atomicInc(addr gs.dev.count, gridDim.x)
    # determine if last block
    isLastBlockDone = (value == (gridDim.x - 1))
  syncThreads()
  # finish the reduction if last block
  if isLastBlockDone:
    var i = threadIdx.x
    var sum = default(T)
    let n = min(gs.npartial, gridDim.x)
    while i < n:
      sum += gs.dev.partial[i]
      i += blockDim.x;
    sum = blockSum(sum)
    # write out the final reduced value
    if threadIdx.x == 0:
      gs.value = sum
      gs.dev.count = 0  # set to zero for next time

when isMainModule:
  type FltArr = UncheckedArray[float32]

  #proc vectorAdd(A: FltArr; B: FltArr; C: var FltArr; n: int32) {.cdecl,hipGlobal.} =
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

    #hipLaunch(vectorAdd, blocksPerGrid, threadsPerBlock, a, b, c, n)
    #discard hipDeviceSynchronize()

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
    echo 3*n

  test()
