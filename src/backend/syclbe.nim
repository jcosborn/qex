import macros, strutils
#import base/metaUtils
import backend/expr
import base/metaUtils
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
  q = dev.queue

template gpuMalloc*(size: SomeInteger): pointer = mallocDevice(size, q)
template gpuMalloc[T](x: var ptr UncheckedArray[T], n: int) =
  x = cast[typeof x](gpuMalloc(n*sizeof(T)))
template gpuMalloc[T](x: ptr T) =
  x = cast[typeof x](gpuMalloc(sizeof(T)))
template gpuFree*(device_ptr: pointer) = freeDevice(device_ptr, q)
template gpuMemCpyToCPU*(dst: pointer, src: pointer; length: SomeInteger) =
  memcpy(q, dst, src, length)
  q.wait
template gpuMemCpyToGPU*(dst: pointer, src: pointer; length: SomeInteger) =
  memcpy(q, dst, src, length)
  q.wait

template gpuThreadNum*: auto =
  #let item = getNdItem1()
  #item[]
  #getNdItem1()[]
  int getGlobalId1()

template gpuNumThreads*: auto =
  #let item = getNdItem1()
  #item.getRange
  int getGlobalRange1()

template syclDefs(body: untyped) =
  setupSycl()
  #var item {.item1.}: Item1
  #template getThreadNum: auto {.used.} = item[]
  #template getNumThreads: auto {.used.} = item.getRange
  #{.emit:["#define nimZeroMem(b,len) memset((b),0,(len))"].}
  inlineProcs:
    body
  #{.emit:["#undef nimZeroMem"].}

proc genCpuPrepare(n:seq[NimNode]):NimNode =
  template r(x,v:untyped):untyped =
    var v = toGpu(x)
    var `v xx` = v
  result = newstmtlist()
  for c in n:
    result.add getast r(c[0],c[1])

proc genCpuFinalize(n:seq[NimNode]):NimNode =
  template r(x,v:untyped):untyped =
    fromGpu(x,v)
    #fromGpu(x)
    #discard
  result = newstmtlist()
  for c in n:
    result.add getast r(c[0],c[1])

proc gpuDefaultNumThreads*(): int =
  32 * q.device.maxComputeUnits.int * q.device.subGroupSizes.max_element.int
  #64 * q.device.maxComputeUnits.int * q.device.subGroupSizes.max_element.int

macro onGpuQ*(q: Queue, n,b,body: untyped): auto =
  let li = body.lineinfo
  let lis = li.split({'/','.','(',','})
  let fl = lis[^4] & "(" & lis[^2] & ")"
  #proc deref(x,g,i:NimNode):auto = newCall("getGpu",x,g)
  proc deref(x,g,i:NimNode):auto = newCall("getGpu",x,newTree(nnkAccQuoted,g,ident"xx"))
  template target(q, fl, n, cpuPrepare, cpuFinalize, body: untyped) =
    mixin toGpu, getGpu, fromGpu
    {.push checks: off.}
    {.push stacktrace: off.}
    block:
      tic(fl)
      inc kernelCallCount
      cpuPrepare  # a let section declare and save device pointers
      toc("cpuPrepare")
      #proc gpuProc {.gensym.} =
      threadSingle:
        let nth = n
        #echo "Launching threads:", nth
        q.submit:
          parallelFor(nth):
            const inOnGpu {.inject,used.} = true
            syclDefs:
              body
      #gpuProc()
      toc("launch")
      proc finalize {.gensym.} =
        tic(fl)
        threadSingle:
          q.wait
        when declared gpuWaitFlops:
          toc("wait",flops=gpuWaitFlops)
        else:
          toc("wait")
        cpuFinalize
        toc("cpuFinalize")
      finalize
  let
    v = prepareVars(body, deref)  # gather gpu pointers in symbols, body is changed accordingly
    cpuPrepare = genCpuPrepare v
    cpuFinalize = genCpuFinalize v
  result = getast(target(q, fl, n, cpuPrepare, cpuFinalize, body))
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

#var gpuNumThreadsRequest* = 0
#var gpuBlockSizeRequest* = 0
template gpuSites(n: int): int = n
template onGpuNowait*(n,b,body: untyped): auto =
  mixin gpuSites
  onGpuQ(q, gpuSites(n), b, body)
template onGpuNowait*(body: untyped): auto =
  let n = gpuDefaultNumThreads()
  onGpuNoWait(n, 0, body)
template onGpuNowait*(n,body: untyped): auto =
  mixin gpuSites
  onGpuNoWait(gpuSites(n), 0, body)

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

proc subgroupSum*[T](x: T): T =  # only thread 0 gets result
  result = x
  let sg = getSubgroup()
  let n = sg.size()
  var offset = n div 2
  while offset >= 1:
    result += sg.shift_left(result, offset)
    offset = offset div 2

template `:=`*(r: LocalPtr[SomeNumber], x: SomeNumber) =
  {.emit:["*",r," = ", x, ";"].}
template `+=`*(r: LocalPtr[SomeNumber], x: SomeNumber) =
  {.emit:["*",r," += ", x, ";"].}
template `+=`*(r: SomeNumber, x: LocalPtr[SomeNumber]) =
  {.emit:[r," += *", x, ";"].}
template `:=`*[N:static int,T](r: LocalPtr[array[N,T]], x: array[N,T]) =
  for i in 0..<N:
    r[i] = x[i]
template `+=`*[N:static int,T](r: LocalPtr[array[N,T]], x: array[N,T]) =
  for i in 0..<N:
    r[i] += x[i]
template `:=`*[N:static int,T](r: array[N,T], x: LocalPtr[array[N,T]]) =
  for i in 0..<N:
    r[i] = x[i]
template `+=`*[N:static int,T](r: array[N,T], x: LocalPtr[array[N,T]]) =
  for i in 0..<N:
    r[i] += x[i]
proc `[]`*[N:static int,T](x: LocalPtr[array[N,T]], i: int): LocalPtr[T] =
  {.emit:[result," = &",x,"[0][",i,"];"].}
proc `[]`*[N,M:static int,T,I,J](x: LocalPtr[array[N,array[M,T]]], i: I, j: J): T =
  {.emit:[result," = ",x,"[0][",i,"][",j,"];"].}
#template `[]`*[N:static int,T](x: LocalPtr[array[N,T]], i: int): LocalPtr[T] =
#  {.emit:["&",x,"[0][",i,"];"].}
#proc `[]`*[N:static int,T](x: LocalPtr[array[N,T]], i: int): LocalPtr[T]
#  {.importcpp:"(&((#)[0][#]))".}
template `[]=`*[N:static int,T](x: LocalPtr[array[N,T]], i: int, y: auto) =
  var t = x[i]
  t := y
template`[]=`*[N:static int,T](x: array[N,T], i: int, y: LocalPtr[T]) =
  {.emit:[x,"[",i,"] = *",y,";"].}

template`setIndexed`*[N,M:static int,T](x: array[N,T], y: LocalPtr[array[M,array[N,T]]], i: SomeInteger) =
  for j in 0..<N:
    x[j] = y[i,j]

proc blockSumSmall*[T](x: T): T = # only thread 0 gets result
  const max_block_size = 1024
  const min_warp_size = 16
  const max_items = max_block_size div min_warp_size
  let g = getGroup()
  let sg = getSubgroup()
  let thread_idx = g.localId #threadIdx.x;
  let block_size = g.size #blockDim.x;
  let warp_size = sg.size
  let warp_idx = thread_idx div warp_size
  let warp_items = (block_size + warp_size - 1) div warp_size
  let warp_thread = thread_idx mod warp_size
  # first do warp reduce
  result = subgroupSum(x)
  #if warp_items == 1: return
  # now do reduction between warps
  g.barrier()
  #var storage {.shared.}: array[max_items, T]
  #var storage = localMem[array[max_items, T]](g)
  #var storage = localMem(array[max_items, T], g)
  var lp = localMem(array[max_items, T], g)
  var storage = lp.get
  # if first thread in warp, write result to shared memory
  var store_idx = warp_idx
  var store_items = warp_items
  while store_items >= 0:
    if store_idx >= 0 and store_idx < max_items and warp_thread == 0:
      if store_idx != warp_idx:
        var t{.noInit.}: T
        #t := storage[store_idx]
        setIndexed(t, storage, store_idx)
        result += t
      storage[store_idx] = result   # apparently can only store to local in one code location
    store_items -= max_items
    store_idx -= max_items
    g.barrier()
  if warp_idx == 0:
    if thread_idx < warp_items:
      result := storage[thread_idx]
      var i = thread_idx + warp_size
      while i < min(warp_items, max_items):
        result += storage[i]
        i += warp_size
    else:
      result = default(T)
    result = subgroupSum(result)

proc blockSum*[T](x: T): T = # only thread 0 gets result
  const min_shared_mem = 48*1024 - sizeof(bool)  # bool used in GpuSum reduce
  const max_items = 64
  const max_size = min_shared_mem div max_items
  when sizeof(T) <= max_size:
    blockSumSmall(x)
  else:
    static: echo $x.type, "  ", sizeof(T)
    {.error:"blockSum: type size too large".}  # FIXME later

type GpuSum*[T] = object
    partial: ptr UncheckedArray[T]
    npartial: cint
    #maxblock: int
    val: ptr T
    count: ptr cint
    #stats: ptr array[3,int]
proc newGpuSum*[T](ns: int): GpuSum[T] =
  let n = (ns + 15) div 16  # divide by warp size
  result.partial.gpuMalloc(n)
  #result.partial.gpuMemset(0, n*sizeof(T))
  q.memset(result.partial, 0, n*sizeof(T))
  result.npartial = cint n
  result.val = cast[ptr T](mallocHost(sizeof(T), q))
  result.count.gpuMalloc()
  #result.count.gpuMemset(0, sizeof(result.count[]))
  q.memset(result.count, 0, sizeof(result.count[]))
  #result.stats.gpuMalloc()  # nthreads, ngroups
template value*(x: GpuSum): auto =
  #var s: array[3,int]
  #gpuMemCpyToCpu(addr s, x.stats, sizeof(x.stats[]))
  #echo "subgroup: ", s[2], "  blockDim: ", s[0], "  gridDim: ", s[1]
  x.val[]
template toGpu*(x: GpuSum): auto = x
template getGpu*(x,g: GpuSum): auto = g
template fromGpu*(x,g: GpuSum): auto = discard

proc atomicInc[T](x: ptr T): T =
  var a = makeAtomicRef(x)
  a.fetchAdd(1)

proc reduce*[T](gs: GpuSum[T], x: T) =
  let g = getGroup()
  let threadIdx = g.localId #threadIdx.x;
  let blockIdx = g.groupId #
  let blockDim = g.size #blockDim.x;
  let gridDim = g.range
  #if threadIdx == 0 and blockIdx == 0:
  #  gs.stats[0] = int blockDim
  #  gs.stats[1] = int gridDim
  #  gs.stats[2] = int getSubgroup().size()
  #if blockDim > 8: gs.stats[2] = blockDim
  var isLastBlockDone = false
  var aggregate = blockSum(x)
  if threadIdx == 0:
    if blockIdx < gs.npartial:
      gs.partial[blockIdx] = aggregate;
    #threadFence() # flush result
    {.emit:["sycl::atomic_fence(sycl::memory_order::release, sycl::memory_scope::device);"].}
    # increment global block counter
    let value = atomicInc(gs.count)
    # determine if last block
    isLastBlockDone = (value == (gridDim - 1))
  isLastBlockDone = g.anyOf(isLastBlockDone)
  #g.barrier()
  # finish the reduction if last block
  if isLastBlockDone:
    var i = threadIdx
    var sum = default(T)
    let n = min(gs.npartial, gridDim)
    {.emit:["sycl::atomic_fence(sycl::memory_order::acquire, sycl::memory_scope::device);"].}
    while i < n:
      sum += gs.partial[i]
      i += blockDim;
    sum = blockSum(sum)
    # write out the final reduced value
    if threadIdx == 0:
      gs.val[] = sum
      gs.count[] = 0  # set to zero for next time


when isMainModule:
  type FltArr = object
    a:ptr UncheckedArray[float32]

  proc test =
    var n = 50000.cint
    var
      a = newSeq[float32](n)
      b = newSeq[float32](n)
      c = newSeq[float32](n)

    template `[]`(x: FltArr, i: SomeInteger): untyped = x.a[][i]
    template `[]=`(x: FltArr, i: SomeInteger, y:untyped):untyped = x.a[][i] = y

    template offloadUseVar(x:seq):bool = true
    template offloadUsePtr(x:seq):bool = true
    template rungpuPrepareOffload(x:seq):bool = true
    template runcpuFinalizeOffload(x:seq):bool = true
    template gpuVarPtr(v:FltArr,p:untyped):untyped = v
    template offloadPtr(x:seq):untyped =
      #let size = x.len * sizeof(x[0])
      #let xp = omp_target_alloc(size)
      #discard omp_target_memcpy_togpu(xp, x[0].addr, size)
      cast[ptr UncheckedArray[type(x[0])]](addr x[0])
    template offloadVar(x:seq,p:untyped):untyped = FltArr(a:p)
    template gpuPrepareOffload(v:FltArr,p:untyped):untyped = v.a=p
    template cpuFinalizeOffload(x:seq,v,p:untyped):untyped =
      discard #omp_target_free(p)

    let sel = HostSelector()
    let dev = sel.selectDevice()
    let q = dev.queue()

    macro dump(n:typed):typed =
      echo n.repr
      n
    dump:
      onGpu(q):
        let tid = getThreadNum()
        let nid = getNumThreads()
        let i0 = (n*tid) div nid
        let i1 = (n*(tid+1)) div nid
        for i in i0..<i1:
          c[i] = a[i] + b[i]
        #for i in countup(tid, n-1, nid):
        #  c[i] = a[i] + b[i]

  test()
