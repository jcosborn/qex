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

proc addChildrenFrom*(dst,src: NimNode): NimNode =
  for c in src: dst.add(c)
  result = dst
macro procInst*(p: typed): auto =
  #echo "begin procInst:"
  #echo p.treerepr
  result = p[0]
macro makeCall*(p: proc, x: tuple): NimNode =
  result = newCall(p).addChildrenFrom(x)

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

proc gpuMemset*[T](p: ptr UncheckedArray[T], val: T, count: int) =
  var v = val  # an array parameter is a host pointer in C, unmapped on the device
  {.emit:["#pragma omp target teams distribute parallel for map(to:",v,")"].}
  {.emit:["for (int i = 0; i < ",count,"; i++)"].}
  block:
    var i {.importc,codegendecl:"".}: cint
    p[i] = v
proc gpuMemset*[T](p: ptr T, val: T) =
  var v = val
  {.emit:["#pragma omp target teams map(to:",v,")"].}
  p[] = v

template toPointer*(x: typed): pointer =
  #dumpType: x
  when x is pointer: x
  elif x is ptr: x
  elif x is seq: toPointer(x[0])
  else: pointer(unsafeAddr(x))
template dataAddr*(x: typed): pointer =
  #dumpType: x
  when x is seq: dataAddr(x[0])
  elif x is array: dataAddr(x[0])
  #elif x is ptr: x
  else: pointer(unsafeAddr(x))
  #else: x

const gpuThreads {.intdefine.} = 128  ## thread_limit of gpuFor kernels, 0 for the runtime default; 128 was best for the staggered kernels on PVC at 24^4 sites per tile
const gpuForMap {.strdefine.} = "defaultmap(firstprivate:pointer)"  ## pointers by value, no lookup in the mapping table at each launch
const gpuForClause* = gpuForMap & (if gpuThreads > 0: " thread_limit(" & $gpuThreads & ")" else: "")
template gpuSubClause(sub: static int): string =
  when sub > 0: " ompx_sub_group_size(" & $sub & ")" else: ""

template gpuFor*(i: untyped; n: SomeInteger; sub: untyped; body: untyped) =
  ## One SPMD kernel over i in 0..<n, in sub-groups of sub threads, or as
  ## the compiler chooses for sub = 0.  Captured pointers must be device
  ## pointers.  n is evaluated on the host, as with the other backends.
  const gpuC = "target teams distribute parallel for " & gpuForClause & gpuSubClause(sub)
  let gpuN = int(n)
  for i in `||`(0, gpuN-1, gpuC):
    body
template gpuFor*(i: untyped; n: SomeInteger; body: untyped) = gpuFor(i, n, 0, body)

var gpuAsyncDep {.exportc.}: int  # orders the gpuForAsync kernels before gpuWaitAsync

template gpuForAsync*(i: untyped; n: SomeInteger; sub: untyped; body: untyped) =
  ## gpuFor returning before the kernel completes; gpuWaitAsync waits for it.
  discard addr(gpuAsyncDep)  # declares it in the C file of the caller
  const gpuC = "target teams distribute parallel for " & gpuForClause & gpuSubClause(sub) & " nowait depend(inout:gpuAsyncDep)"
  let gpuN = int(n)
  for i in `||`(0, gpuN-1, gpuC):
    body
template gpuForAsync*(i: untyped; n: SomeInteger; body: untyped) = gpuForAsync(i, n, 0, body)

template gpuWaitAsync* =
  discard addr(gpuAsyncDep)
  {.emit: "#pragma omp taskwait depend(inout:gpuAsyncDep)".}

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

template gpuForObj(i: untyped; n: SomeInteger; body: untyped) =
  ## gpuFor for a kernel of onGpu that captures objects: they go to the
  ## kernel by value, which icx does for a kernel that completes at the
  ## launch (a mapped object costs ten times the launch: 270 against 25 us
  ## on PVC)
  const gpuC = "target teams distribute parallel for " & gpuForClause & " defaultmap(firstprivate:aggregate)"
  let gpuN = int(n)
  for i in `||`(0, gpuN-1, gpuC):
    body

proc gpuFinalizer*[P](fl: static string; count: int; flops: float; ps: P; back: proc(ps: P) {.nimcall.}): proc() =
  ## The finalizer of a kernel launch: waits for the queue, then back copies
  ## what the kernel wrote to the host variables of ps (hostBack).  The
  ## closure has an environment of its own, so every thread of a threads
  ## block gets one and may keep it.
  result = proc() =
    tic(fl)
    var countSave = -1  # the count before, in the thread that lowers it
    threadSingle:
      countSave = kernelCallCount
      kernelCallCount = count
      gpuWaitAsync()
    toc("wait",flops=flops)
    back(ps)
    if countSave >= 0: kernelCallCount = countSave  # the same thread restores it
    threadBarrier()
    toc("cpuFinalize")

macro gpuKernel(fin: static bool; n, body: untyped): untyped =
  ## Launches body as a kernel of n threads, n a count or what gpuSites
  ## turns into one.  Its captures are the device values of the variables
  ## that body names (backend/expr).  With fin the result is the finalizer,
  ## which waits for the queue and copies back.
  ## The kernel is that of gpuForAsync when every captured value is a number
  ## or a pointer.  With an object among them it completes at the launch,
  ## after the queued kernels: icx reads a captured object from the stack of
  ## the launch when the kernel runs, with every clause (icx 2026.1 on PVC), so a
  ## queued kernel would read it after the launch has returned.
  ## The launch stands in a single without its barrier, at which a queued
  ## kernel would complete.  In a threads block the launching thread waits
  ## for the kernel, as only it can, and the team meets at a barrier.
  let li = body.lineinfo
  let lis = li.split({'/','.','(',','})
  let fl = lis[^4] & "(" & lis[^2] & ")"
  template target(fl, n, pre, shadow, hostPairs, queued, body, fin: untyped) =
    mixin toGpu, getGpu, fromGpu, gpuSites
    inc kernelCallCount
    block:
      let thisKernelCallCount = kernelCallCount  # save current value for finalizer
      tic(fl)
      let gpuNt {.inject.} = int(gpuSites(n))
      gpuLaunchThreads = gpuNt
      pre  # the device values
      when fin:
        let hp = hostPairs  # before the kernel values take the names of the host variables
      shadow  # the kernel values: locals of the launch that the kernel captures and reads in place
      toc("cpuPrepare")
      ompBlock("single nowait"):
        when declared(gpuObjects):
          gpuWaitAsync()
          gpuForObj(gpuTn, gpuNt):
            const inOnGpu {.inject,used.} = true
            inlineProcs:
              body
        else:
          gpuForAsync(gpuTn, gpuNt):
            const inOnGpu {.inject,used.} = true
            inlineProcs:
              body
          if numThreads > 1: gpuWaitAsync()
      if numThreads > 1: threadBarrier()
      toc("launch")
      inc kernelCallCount  # increment after launch
      when fin:
        gpuFinalizer(fl, thisKernelCallCount, (when declared gpuWaitFlops: float(gpuWaitFlops) else: 0.0), hp, hostBack(hp))
      else:
        queued
  let c = kernelCap(body, n)
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

template toUArray(a:untyped):untyped = cast[ptr UncheckedArray[typeof(a[0])]](a[0].unsafeaddr)
proc cleanAst(n:NimNode):NimNode =
  if n.kind in {nnkHiddenDeref,nnkHiddenCallConv,nnkHiddenStdConv}:
    result = n[0].cleanAst
  else:
    result = n.copyNimNode
    for c in n:
      result.add c.cleanAst
proc identStr(n:NimNode):string =
  result = n.repr
  for i in 0..<result.len:
    if result[i] in {'.','[',']',':'}: result[i] = '_'
#proc isIndex(n,i:NimNode):bool =
#  result = n.eqident i
#  if n.kind == nnkHiddenStdConv:
#    result = n[1].eqident i
macro simdForImpl(n:typed):untyped =
  proc getIndexedPtrs(n,i:NimNode):(NimNode,seq[NimNode]) =
    #echo "### getIndexedPtrs: ", i.repr
    #echo n.treerepr
    var ptrs = newseq[NimNode]()
    proc get(n:NimNode):NimNode =
          var m = -1
          for j in 0..<ptrs.len:
            if ptrs[j][1] == n:
              m = j
              break
          if m < 0:
            let v = gensym(nskVar, n.cleanAst.identStr)
            #ptrs.add newPar(v, n)
            ptrs.add newNimNode(nnkTupleConstr).add(v, n)
            return v
          else:
            return ptrs[m][0]
    proc go(n:NimNode):NimNode =
      result = n.copyNimNode
      if n.kind in CallNodes and ($n[0] == "[]" or $n[0] == "[]="):
        if n.len > 2: # and n[2].isIndex i:
          result.add n[0].go
          result.add n[1].get
          for i in 2..<n.len: result.add n[i].go
        else:
          for c in n: result.add c.go
      elif n.kind == nnkBracketExpr:
        result.add n[0].get
        for i in 1..<n.len: result.add n[i].go
      else:
        for c in n: result.add c.go
    var nn = n.go
    (nn, ptrs)
  template res(setup, i, lo, hi, body: untyped): untyped =
    block:
      var i {.codegendecl:"/* $# $# */",noinit.}: cint
      setup
      {.emit:
        ["\n#pragma omp simd aligned(","\n",
          "for(int ",
          i,"=",lo,";",
          i,"<=",hi,";",
          i,"++){\n"
        ].}
      body
      {.emit:["\n}\n"].}

  #echo n.treerepr
  n.expectkind nnkForStmt
  #echo n[1][0].getimpl.treerepr
  let (nn,ptrs) = n[2].getIndexedPtrs(n[0])
  if ptrs.len == 0:
    echo "simdForImpl finds no pointers: ",n.treerepr
    quit 1
  let setup = newNimNode nnkVarSection
  for p in ptrs:
    setup.add newIdentDefs(
      #p[0],
      newNimNode(nnkPragmaExpr).add(
        p[0],
        newNimNode(nnkPragma).add(
          newNimNode(nnkExprColonExpr).add(ident"codegenDecl", newLit"$# __restrict__ $#"))),
      newEmptyNode(), newcall(bindsym"toUArray", p[1]))
  result = getast res(setup, n[0], n[1][1], n[1][2], nn)
  #echo result.treerepr
  let e = result[1][2]
  e.expectkind nnkPragma
  for i in 0..<ptrs.len:
    if i == 0: e[0][1].insert(1,newLit")")
    else: e[0][1].insert(1,newLit",")
    e[0][1].insert(1,ptrs[i][0])
  let i = gensym(nskvar, $n[0])
  result = result.replace(n[0], i).rebuild
  #echo result.repr
  #echo result.treerepr
  #quit 1

macro simdFor*(n:untyped):untyped =
  proc p(n:NimNode):NimNode =
    if n.kind == nnkForStmt:
      n[2] = newCall(bindsym"inlineProcs", n[2])
      #echo n.treerepr
      return newCall(bindsym"simdForImpl", n)
    elif n.kind == nnkStmtList:
      for i in 0..<n.len:
        n[i] = p n[i]
      return n
    else:
      echo "simdFor cannot handle:"
      echo n.treerepr
      quit 1
  p n

