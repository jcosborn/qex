import macros
#import base/metaUtils
import backend/expr
import sycl

let dev = defaultDevice()
let q = dev.queue

template gpuMalloc*(size: SomeInteger): pointer = mallocDevice(size, q)
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
  {.emit:["#define nimZeroMem(b,len) memset((b),0,(len))"].}
  #inlineProcs:
  body
  {.emit:["#undef nimZeroMem"].}

proc genCpuPrepare(n:seq[NimNode]):NimNode =
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

proc gpuDefaultNumThreads*(): int =
  64 * q.device.maxComputeUnits.int * q.device.subGroupSizes.max_element.int

macro onGpuQ*(q: Queue, body: untyped): auto =
  proc deref(x,g,i:NimNode):auto = newCall("getGpu",x,g)
  template target(q, cpuPrepare, cpuFinalize, body: untyped) =
    mixin toGpu, getGpu, fromGpu
    {.push checks: off.}
    {.push stacktrace: off.}
    proc gpuProc {.gensym.} =
      cpuPrepare  # a let section declare and save device pointers
      #let nth = 32 * q.device.maxComputeUnits.int * q.device.preferredVectorWidthFloat.int
      #let nth = 32 * q.device.maxComputeUnits.int
      let nth = gpuDefaultNumThreads()
      #echo "Launching threads:", nth
      q.submit:
        parallelFor(nth):
          syclDefs:
            body
      q.wait
      cpuFinalize
    gpuProc()
  let
    v = prepareVars(body, deref)  # gather gpu pointers in symbols, body is changed accordingly
    cpuPrepare = genCpuPrepare v
    cpuFinalize = genCpuFinalize v
  result = getast(target(q, cpuPrepare, cpuFinalize, body))
  echo result.repr


# XXX fix the following
template onGpu*(body: untyped) = onGpuQ(q, body)
#template onGpu*(totalNumThreads, body: untyped): untyped = onGpu(body)
#template onGpu*(totalNumThreads, numThreadsPerTeam, body: untyped): untyped = onGpu(body)


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
