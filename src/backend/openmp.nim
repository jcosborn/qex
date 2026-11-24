import macros, strutils
import base/[metaUtils, threading, profile]
import base/omp
import backend/expr

var kernelCallCount* = 0  # count of kernel calls
const dumpKernels {.intdefine.} = 0

{.pragma: omp, header:"omp.h".}
{.passC:"-fcf-protection=none -no-pie -fno-stack-protector" .}
{.passL:"-fcf-protection=none -no-pie -fno-stack-protector" .}

template mkMemoryPragma*:untyped =
  {.pragma: restrict, codegenDecl: "$# __restrict__ $#".}
  {.pragma: aligned, codegenDecl: "$# $# __attribute__((aligned))".}
  {.pragma: aligned1, codegenDecl: "$# $# __attribute__((aligned(1)))".}
  {.pragma: aligned2, codegenDecl: "$# $# __attribute__((aligned(2)))".}
  {.pragma: aligned4, codegenDecl: "$# $# __attribute__((aligned(4)))".}
  {.pragma: aligned8, codegenDecl: "$# $# __attribute__((aligned(8)))".}
  {.pragma: aligned16, codegenDecl: "$# $# __attribute__((aligned(16)))".}
  {.pragma: aligned32, codegenDecl: "$# $# __attribute__((aligned(32)))".}
  {.pragma: aligned64, codegenDecl: "$# $# __attribute__((aligned(64)))".}
  {.pragma: aligned128, codegenDecl: "$# $# __attribute__((aligned(128)))".}
  {.pragma: aligned256, codegenDecl: "$# $# __attribute__((aligned(256)))".}
  {.pragma: alignedType, codegenDecl: "$# $# __attribute__((aligned(sizeof($1))))".}

proc alignatImpl(n:NimNode, byte:int): NimNode =
  result = n.copyNimNode
  if n.kind == nnkIdentDefs:
    let a = ident("aligned" & $byte)
    for i in 0..<n.len-2:
      if n[i].kind == nnkPragmaExpr:
        result.add n[i]
        result[i][1].expectKind nnkPragma
        result[i][1].add a
      else:
        result.add newNimNode(nnkPragmaExpr).add(n[i], newNimNode(nnkPragma).add a)
    for i in n.len-2..<n.len:
      result.add n[i]
  else:
    for c in n:
      result.add c.alignatImpl byte
macro alignat*(byte:static[int], n:untyped): untyped =
  if byte notin [1,2,4,8,16,32,64,128,256]:
    error("macro alignat: unsupported alignment: " & $byte, n)
  #echo "alignatImpl ", byte
  #echo n.treerepr
  result = n.alignatImpl byte
  #error result.treerepr

proc omp_target_alloc*(size: csize_t, device_num: cint): pointer {.omp.}
proc omp_target_free*(device_ptr: pointer, device_num: cint) {.omp.}
proc omp_target_alloc_host*(size: csize_t, device_num: cint): pointer {.omp.}
proc omp_target_free_host*(device_ptr: pointer, device_num: cint) {.omp.}
proc omp_target_alloc_shared*(size: csize_t, device_num: cint): pointer {.omp.}
proc omp_target_free_shared*(device_ptr: pointer, device_num: cint) {.omp.}
proc omp_target_memcpy*(dst: pointer, src: pointer;
    length, dst_offset, src_offset: csize_t;
    dst_device_num, src_device_num: cint): cint {.omp.}
proc omp_get_default_device*: cint {.omp.}
proc omp_get_initial_device*: cint {.omp.}
proc omp_get_num_teams*: cint {.omp.}
proc omp_get_team_num*: cint {.omp.}
proc omp_alloc*(size: SomeInteger): pointer =
    {.emit:["omp_alloc(",size,", omp_default_mem_alloc);"].}

template omp_target_alloc*(size: SomeNumber): pointer =
  omp_target_alloc(csize_t size, omp_get_default_device())
template omp_target_alloc_host*(size: SomeNumber): pointer =
  omp_target_alloc_host(csize_t size, omp_get_default_device())
template omp_target_alloc_shared*(size: SomeNumber): pointer =
  omp_target_alloc_shared(csize_t size, omp_get_default_device())
template omp_target_memcpy_tocpu*(dst: pointer, src: pointer; length: csize_t): cint =
  omp_target_memcpy(dst, src, length, 0, 0, omp_get_initial_device(), omp_get_default_device())
template omp_target_memcpy_togpu*(dst: pointer, src: pointer; length: csize_t): cint =
  omp_target_memcpy(dst, src, csize_t length, 0, 0, omp_get_default_device(), omp_get_initial_device())
template omp_target_free*(device_ptr: pointer) =
  omp_target_free(device_ptr, omp_get_default_device())

template gpuMalloc*(size:SomeInteger):pointer = omp_target_alloc(size)
template gpuFree*(device_ptr:pointer) = omp_target_free(device_ptr)
proc gpuMemCpyToCPU*(dst: pointer, src: pointer; length: SomeInteger): cint {.discardable.} =
  omp_target_memcpy_tocpu(dst, src, csize_t length)
proc gpuMemCpyToGPU*(dst: pointer, src: pointer; length: SomeInteger): cint {.discardable.} =
  omp_target_memcpy_togpu(dst, src, csize_t length)
template gpuMalloc[T](x: var ptr UncheckedArray[T], n: int) =
  x = cast[typeof x](gpuMalloc(n*sizeof(T)))
template gpuMalloc[T](x: ptr T) =
  x = cast[typeof x](gpuMalloc(sizeof(T)))

template toPointer*(x: typed): pointer =
  #dumpType: x
  when x is pointer: x
  elif x is ptr: x
  elif x is seq: toPointer(x[0])
  else: pointer(unsafeAddr(x))

const gpuThreads {.intdefine.} = 128  ## thread_limit of the kernels, 0 for the runtime default; 128 was best for the staggered kernels on PVC at 24^4 sites per tile
const gpuForMap {.strdefine.} = "defaultmap(firstprivate:pointer)"  ## pointers by value, no lookup in the mapping table at each launch
const gpuForClause* = gpuForMap & (if gpuThreads > 0: " thread_limit(" & $gpuThreads & ")" else: "")
template gpuSubClause(sub: static int): string =
  when sub > 0: " ompx_sub_group_size(" & $sub & ")" else: ""

var gpuAsyncDep {.exportc.}: int  # orders the gpuForAsync kernels before gpuWaitAsync
var gpuQueued = false  # a kernel queued outside parallel regions may wait in the queue

template gpuForAsync(i: untyped; n: SomeInteger; sub: untyped; body: untyped) =
  ## The queued kernel of an onGpu launch: one SPMD kernel over i in 0..<n,
  ## in sub-groups of sub threads, or as the compiler chooses for sub = 0,
  ## returning before it completes; gpuWaitAsync waits for it.  Captured
  ## pointers must be device pointers.
  discard addr(gpuAsyncDep)  # declares it in the C file of the caller
  const gpuC = "target teams distribute parallel for " & gpuForClause & gpuSubClause(sub) & " nowait depend(inout:gpuAsyncDep)"
  let gpuN = int(n)
  for i in `||`(0, gpuN-1, gpuC):
    body
  gpuQueued = true

template gpuWaitAsync* =
  ## a task wait for the queued kernels, skipped outside parallel regions
  ## where none was queued since the last wait.  Inside a region (a threads
  ## block, also of one thread) the wait covers the kernels queued in it, and
  ## the flag of the kernels queued before it stays.
  discard addr(gpuAsyncDep)
  let gpuL = omp_get_level()
  if gpuL > 0 or gpuQueued:
    {.emit: "#pragma omp taskwait depend(inout:gpuAsyncDep)".}
    if gpuL == 0: gpuQueued = false

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

template gpuMemCpyToCpuAsync*(dst, src: pointer; length: SomeInteger) =
  ## the copy to the host after the kernels already submitted, complete here
  gpuWaitAsync()
  gpuMemCpyToCPU(dst, src, length)
const gpuCopyFast* = high(int)  ## bytes up to which a queued copy has its low cost

template gpuMallocHost*(size: SomeInteger): pointer =
  ## pinned host memory, for fast copies from the device
  omp_target_alloc_host(size)
template gpuFreeHost*(p: pointer) = omp_target_free(p)

proc ompPlain(t: NimNode; depth = 0): bool =
  ## whether a value of type t is a number, bool, enum, character or
  ## pointer, which a deferred target region gets intact
  case t.typeKind
  of ntyBool, ntyChar, ntyEnum, ntyRange, ntyPtr, ntyPointer, ntyCString,
     ntyInt, ntyInt8, ntyInt16, ntyInt32, ntyInt64, ntyUInt, ntyUInt8, ntyUInt16, ntyUInt32, ntyUInt64,
     ntyFloat, ntyFloat32, ntyFloat64, ntyFloat128: true
  of ntyDistinct, ntyGenericInst, ntyAlias, ntyVar:
    if depth > 8: return false
    let i = t.getTypeImpl
    case i.kind
    of nnkDistinctTy, nnkVarTy: ompPlain(i[0].getType, depth+1)
    of nnkPtrTy: true
    of nnkSym: ompPlain(i.getType, depth+1)
    else: false  # objects, tuples, arrays, seqs, refs
  else: false

template gpuForSync(i: untyped; n: SomeInteger; sub: static int; body: untyped) =
  ## A kernel that completes at its launch, after the queued kernels: a
  ## wait where the queue holds work, then a target without nowait and
  ## without depend, which costs the host less than a target that depends
  ## on the queue (on PVC 0.5 to 1.2 us less per application of the
  ## staggered operator at 16^4).  The kernels of onGpu, and the queued ones
  ## that read objects; objects go to the kernel by value, which icx does
  ## for a kernel that completes at the launch (a mapped object costs ten
  ## times the launch: 270 against 25 us on PVC).
  gpuWaitAsync()
  const gpuC = "target teams distribute parallel for " & gpuForClause & gpuSubClause(sub) & " defaultmap(firstprivate:aggregate)"
  let gpuN = int(n)
  for i in `||`(0, gpuN-1, gpuC):
    body

macro ompKernel(i, n, sub, kernel: untyped; loop: typed): untyped =
  ## loop: the queued kernel of an onGpu launch (gpuForAsync over i of n),
  ## kernel: its body untyped.  When the typed body reads an object, array or
  ## tuple from outside the kernel, the result is gpuForSync over kernel: icx
  ## 2026.1 reads a captured object from the stack of the launch when a
  ## deferred kernel runs on PVC, and the captures of a launch in top level
  ## code are globals that the next launch rewrites.  Numbers and pointers
  ## stay queued.  The typed body decides.
  proc firstFor(n: NimNode): NimNode =
    ## the loop of the kernel, the one over ||
    if n.kind == nnkParForStmt: return n
    for c in n:
      let f = firstFor(c)
      if f != nil: return f
  let fs = firstFor(loop)
  doAssert fs != nil, "ompKernel: no loop"
  var decl: seq[NimNode]
  proc defs(n: NimNode) =
    case n.kind
    of nnkIdentDefs, nnkConstDef:
      for c in n[0..^3]: decl.add(if c.kind == nnkPragmaExpr: c[0] else: c)
    of nnkVarTuple:  # let (a, b) = v
      for c in n[0..^2]: decl.add c
    of nnkForStmt, nnkParForStmt:
      for c in n[0..^3]:
        if c.kind == nnkVarTuple:  # for (a, b) in ...
          for d in c:
            if d.kind == nnkSym: decl.add d
        else: decl.add c
      defs(n[^1])
      return
    else: discard
    for c in n: defs(c)
  defs(fs)
  var obj = false
  proc reads(n: NimNode) =
    case n.kind
    of RoutineNodes: discard
    of nnkSym:
      if n.symKind in {nskLet, nskVar, nskParam, nskForVar, nskResult, nskTemp} and n notin decl and
         not ompPlain(n.getType):
        obj = true
    else:
      for c in n: reads(c)
  reads(fs[^1])
  if not obj: return loop
  result = newCall(bindSym"gpuForSync", i, n, sub, kernel)

proc gpuFinish*[P](fl: static string; count: int; flops: float; ps: P; back: proc(ps: P) {.nimcall.}; wait: static bool) =
  ## The end of a kernel launch: waits for the queue where the kernel was
  ## queued (wait), then back copies what the kernel wrote to the host
  ## variables of ps (hostBack).
  tic(fl)
  var countSave = -1  # the count before, in the thread that lowers it
  threadSingle:
    countSave = kernelCallCount
    kernelCallCount = count
    when wait: gpuWaitAsync()
  toc("wait",flops=flops)
  back(ps)
  if countSave >= 0: kernelCallCount = countSave  # the same thread restores it
  threadBarrier()
  toc("cpuFinalize")

proc gpuFinalizer*[P](fl: static string; count: int; flops: float; ps: P; back: proc(ps: P) {.nimcall.}): proc() =
  ## gpuFinish as a closure, for onGpuNowait.  The closure has an
  ## environment of its own, so every thread of a threads block gets one
  ## and may keep it.
  result = proc() = gpuFinish(fl, count, flops, ps, back, true)

macro gpuKernel(fin: static int; n: untyped; sub: static int; body: untyped): untyped =
  ## Launches body as a kernel of n threads, n a count or what gpuSites
  ## turns into one.  Its captures are the device values of the variables
  ## that body names (backend/expr), locals of the launch that the kernel
  ## reads in place.  fin: 0 queues the kernel, 1 gives the finalizer, which
  ## waits for the queue and copies back, 2 runs the kernel at the launch,
  ## after the queued kernels (gpuForSync), and copies back.
  ## A kernel of fin 0 or 1 is that of gpuForAsync when it reads numbers and
  ## pointers only, else it completes at the launch (ompKernel).
  ## The launch stands in a single without its barrier, at which a queued
  ## kernel would complete.  In a threads block the launching thread waits
  ## for the kernel, as only it can, and the team meets at a barrier.
  let li = body.lineinfo
  let lis = li.split({'/','.','(',','})
  let fl = lis[^4] & "(" & lis[^2] & ")"
  template target(fl, n, sub, pre, shadow, hostPairs, queued, copies, body, fin: untyped) =
    mixin toGpu, getGpu, fromGpu, gpuSites
    inc kernelCallCount
    block:
      let thisKernelCallCount = kernelCallCount  # save current value for finalizer
      tic(fl)
      let gpuNt {.inject.} = int(gpuSites(n))
      gpuLaunchThreads = gpuNt
      pre  # the device values
      when fin == 1 or (fin == 2 and copies):
        let hp = hostPairs  # before the kernel values take the names of the host variables
      shadow  # the kernel values: locals of the launch that the kernel captures and reads in place
      ompBlock("single nowait"):
        when fin == 2:
          gpuForSync(gpuTn, gpuNt, sub):
            const inOnGpu {.inject,used.} = true
            inlineProcs:
              body
        else:
          ompKernel(gpuTn, gpuNt, sub, (const inOnGpu {.inject,used.} = true; inlineProcs: body)):
            gpuForAsync(gpuTn, gpuNt, sub):
              const inOnGpu {.inject,used.} = true
              inlineProcs:
                body
          if numThreads > 1: gpuWaitAsync()
      if numThreads > 1: threadBarrier()
      inc kernelCallCount  # increment after launch
      # one profile timer per launch: "wait" for onGpu, whose kernel completes here, else "launch"
      when fin == 2:
        toc("wait", flops=(when declared gpuWaitFlops: float(gpuWaitFlops) else: 0.0))
        when copies: gpuFinish(fl, thisKernelCallCount, (when declared gpuWaitFlops: float(gpuWaitFlops) else: 0.0), hp, hostBack(hp), false)
      else:
        toc("launch")
        when fin == 1: gpuFinalizer(fl, thisKernelCallCount, (when declared gpuWaitFlops: float(gpuWaitFlops) else: 0.0), hp, hostBack(hp))
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

{.emit: "/*INCLUDESECTION*/\n#include <omp.h>".}
proc gpuZeContext*(): tuple[ctx, dev: pointer] =
  ## Level Zero context and device of the default OpenMP device
  var ctx, dev: pointer
  {.emit: """
  omp_interop_t obj = omp_interop_none;
  int dn = omp_get_default_device(), err;
  #pragma omp interop init(targetsync: obj) device(dn)
  `ctx` = omp_get_interop_ptr(obj, omp_ipr_device_context, &err);
  `dev` = omp_get_interop_ptr(obj, omp_ipr_device, &err);
  #pragma omp interop destroy(obj)
  """.}
  (ctx, dev)
