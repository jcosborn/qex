import macros
import base/metaUtils
import base/omp
import backend/expr

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
  {.emit:["#pragma omp target teams distribute parallel for"].}
  {.emit:["for (int i = 0; i < ",count,"; i++)"].}
  block:
    var i {.importc,codegendecl:"".}: cint
    p[i] = val
proc gpuMemset*[T](p: ptr T, val: T) =
  {.emit:["#pragma omp target teams"].}
  p[] = val

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

template gpuThreadNum*: untyped =
  let teamNum = omp_get_team_num()
  let numThreads = omp_get_num_threads()
  let threadNum = omp_get_thread_num()
  teamNum.int * numThreads.int + threadNum.int
template gpuNumThreads*: untyped =
  let numTeams = omp_get_num_teams()
  let numThreads = omp_get_num_threads()
  numTeams.int * numThreads.int

template openmpDefs(n,body: untyped) =
  # XXX check if GC matters
  #let
  #  numTeams = omp_get_num_teams()
  #  teamNum = omp_get_team_num()
  #ompBlock("parallel num_threads(512)"):
  #ompBlock("parallel"):
  ompBlock2("parallel num_threads(", n, ")"):
    #let
    #  numThreads = omp_get_num_threads()
    #  threadNum = omp_get_thread_num()
    #template getThreadNum: untyped {.used.} = teamNum.int * numThreads.int + threadNum.int
    #template getNumThreads: untyped {.used.} = numTeams.int * numThreads.int
    #{.emit:["#define nimZeroMem(b,len) memset((b),0,(len))"].}
    inlineProcs:
      body
    #{.emit:["#undef nimZeroMem"].}


proc genCpuPrepare(n:seq[NimNode]):NimNode =
  mixin toGpu
  template r(x,v:untyped):untyped =
    var v = toGpu(x)
    var `v xx` = v
  result = newstmtlist()
  for c in n:
    result.add getast r(c[0],c[1])

proc genCpuFinalize(n:seq[NimNode]):NimNode =
  template r(x,v:untyped):untyped =
    fromGpu(x,v)
  result = newstmtlist()
  for c in n:
    result.add getast r(c[0],c[1])

proc declarePtrTuple(n:seq[NimNode]):NimNode =
  mixin offloadPtr
  var ps = newNimNode(nnkTupleConstr)
  ps.add newLit"is_device_ptr("
  for c in n:
    when compiles(offloadPtr(c[0])):
      ps.add c[2]
  if ps.len == 1:
    result = newNimNode(nnkTupleConstr)
  else:
    ps.add newLit")"
    result = ps
  #echo result.treerepr

#macro isDevicePtr(x: typed): untyped =
#  let n = $x
#  result = newLit(" is_device_ptr("&n&")")

#macro useDevicePtr(x: auto): auto =
  #echo x.treerepr
  #let n = x.strVal
  #echo "useDevicePtr: ", n
  #let p = newLit("#pragma omp target data use_device_ptr("&n&")")
#  result = quote do:
#    {.emit: ["#pragma omp target data use_device_ptr(",`x`,")"].}

#macro getrepr(x: untyped): auto =
#  echo x.treerepr
#  result = x

template useDevicePtr*(x: auto) =
  #getrepr:
  {.emit: ["#pragma omp target data use_device_ptr(",x,")"].}

#macro mapto(x: typed): untyped =
#  let n = $x
#  result = newLit(" map(to:"&n&")")
#macro mapto(x: typed): untyped =

macro onGpuNowait*(n,b,body: untyped): auto =
  let li = body.lineinfo
  #proc deref(x,g,i:NimNode):auto = newCall("getGpu",x,g)
  proc deref(x,g,i:NimNode):auto = newCall("getGpu",x,newTree(nnkAccQuoted,g,ident"xx"))
  template target(n,b,cpuPrepare, cpuFinalize, devicePtrDeclare, body: untyped) =
    mixin toGpu, getGpu, fromGpu
    {.push checks: off.}
    {.push stacktrace: off.}
    block:
      cpuPrepare  # a let section declare and save device pointers
      #proc gpuProc {.gensym.} =
      threadSingle:
        let nthreads = n
        let nteams = n div b
        #ompBlock2("target teams num_teams(1024)", devicePtrDeclare):
        #ompBlock2("target teams", devicePtrDeclare):
        ompBlock2("target teams num_teams(", nteams, ")", devicePtrDeclare):
          openmpDefs(nthreads):
            const inOnGpu {.inject,used.} = true
            body
      #gpuProc()
      proc finalize {.gensym.} =
        cpuFinalize
        #threadBarrier()
      finalize
  let
    v = prepareVars(body, deref)  # gather gpu pointers in symbols, body is changed accordingly
    cpuPrepare = genCpuPrepare v
    cpuFinalize = genCpuFinalize v
    isDevicePtrs = declarePtrTuple v
  result = getast(target(n,b,cpuPrepare, cpuFinalize, isDevicePtrs, body))
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

var gpuNumThreadsRequest* = 32*1024
var gpuBlockSizeRequest* = 64
template gpuSites(n: int): int = n
template onGpuNowait*(body: untyped): auto =
  onGpuNoWait(gpuNumThreadsRequest, gpuBlockSizeRequest, body)
template onGpuNowait*(n0,body: untyped): auto =
  mixin gpuSites
  let n = gpuSites(n0)
  var b = gpuBlockSizeRequest
  while b > n: b = b div 2
  onGpuNoWait(n, b, body)
#template onGpuNowait*(n,b,body: untyped): auto =
#  onGpuNoWait(n, b, body)

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

proc blockSumSmall*[T](x: T): T = # only thread 0 gets result
  const max_block_size = 1024
  const max_items = max_block_size
  let thread_idx = omp_get_thread_num()
  let block_size = omp_get_num_threads()
  var storage {.noInit,codegendecl:"static $# $#".}: array[max_items, T]
  {.emit:["#pragma omp groupprivate(",storage,")"].}
  storage[thread_idx] = x
  {.emit:["#pragma omp barrier"].}
  if thread_idx == 0:
    result = x
    for i in 1..<block_size:
      result += storage[i]

proc blockSum*[T](x: T): T = # only thread 0 gets result
  const min_shared_mem = 48*1024 - sizeof(bool)  # bool used in GpuSum reduce
  const max_items = 512 #1024
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
    valh: T
    valvalid: bool
    count: ptr cint
proc newGpuSum*[T](ns: int): GpuSum[T] =
  let n = (ns + 15) div 16  # divide by warp size
  result.partial.gpuMalloc(n)
  result.partial.gpuMemset(default(T), n)
  result.npartial = cint n
  #result.val = cast[ptr T](omp_target_alloc_host(csize_t sizeof(T)))
  #result.val = cast[ptr T](omp_target_alloc_shared(csize_t sizeof(T)))
  result.val = cast[ptr T](omp_target_alloc(csize_t sizeof(T)))
  result.count.gpuMalloc()
  #result.count.gpuMemset(0, sizeof(result.count[]))
  #q.memset(result.count, 0, sizeof(result.count[]))
  result.count.gpuMemset(0)
#template value*(x: GpuSum): auto = x.val[]
template value*(x: GpuSum): auto =
  if not x.valvalid:
    x.valvalid = true
    gpuMemCpyToCPU(addr x.valh, x.val, sizeof(x.valh))
  x.valh
template toGpu*(x: GpuSum): auto =
  x.valvalid = false
  x
template getGpu*(x,g: GpuSum): auto = g
template fromGpu*(x,g: GpuSum): auto = discard

proc reduce*[T](gs: GpuSum[T], x: T) =
  var aggregate = blockSum(x)
  let threadIdx = omp_get_thread_num()
  let blockIdx = omp_get_team_num()
  let blockDim = omp_get_num_threads()
  let gridDim = omp_get_num_teams()
  block:
    var isLastBlockDone{.noInit,codegendecl:"static $# $#".}: bool
    {.emit:["#pragma omp groupprivate(",isLastBlockDone,")"].}
    #isLastBlockDone = false
    if threadIdx == 0:
      if blockIdx < gs.npartial:
        gs.partial[blockIdx] = aggregate;
      #threadFence() # flush result
      {.emit:["#pragma omp flush release"].}
      # increment global block counter
      #let value = atomicInc(gs.count)
      var value: typeof gs.count[]
      {.emit:["#pragma omp atomic capture"].}
      block:
        value = gs.count[]
        gs.count[] += 1
      # determine if last block
      isLastBlockDone = (value == (gridDim - 1))
    {.emit:["#pragma omp barrier"].}
    # finish the reduction if last block
    if isLastBlockDone:
      var i = threadIdx
      var sum = default(T)
      let n = min(gs.npartial, gridDim)
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
  {.emit:"#pragma omp requires unified_shared_memory".}

  proc test =
    var n = 50000.cint
    var
      a = newSeq[float32](n)
      b = newSeq[float32](n)
      c = newSeq[float32](n)

    template `[]`(x: FltArr, i: SomeInteger): untyped = x.a[][i]
    template `[]=`(x: FltArr, i:SomeInteger, y:untyped):untyped = x.a[][i] = y

    #template offloadUseVar(x:seq):bool = true
    #template offloadUsePtr(x:seq):bool = true
    #template rungpuPrepareOffload(x:seq):bool = true
    #template runcpuFinalizeOffload(x:seq):bool = true
    #template gpuVarPtr(v:FltArr,p:untyped):untyped = v
    #template offloadPtr(x:seq):untyped =
    #  let size = x.len * sizeof(x[0])
    #  let xp = omp_target_alloc(size)
    #  discard omp_target_memcpy_togpu(xp, x[0].addr, size)
    #  cast[ptr UncheckedArray[type(x[0])]](xp)
    #template offloadVar(x:seq,p:untyped):untyped = FltArr(a:p)
    #template gpuPrepareOffload(v:FltArr,p:untyped):untyped = v.a=p
    #template cpuFinalizeOffload(x:seq,v,p:untyped):untyped = omp_target_free(p)

    template toGpu(x: cint): auto = x
    template getGpu(x: cint, g: cint): auto = x
    template fromGpu(x: cint, g: cint) = discard

    template toGpu(x: seq): auto =
      let size = x.len * sizeof(x[0])
      #let xp = omp_target_alloc(size)
      let xp = gpuMalloc(size)
      #discard omp_target_memcpy_togpu(xp, x[0].addr, size)
      gpuMemCpyToGpu(xp, x[0].addr, size)
      FltArr(a:cast[ptr UncheckedArray[type(x[0])]](xp))

    template getGpu(x: seq, g: FltArr): auto = g

    template fromGpu(x:seq, g:FltArr) =
      let size = x.len * sizeof(x[0])
      gpuMemCpyToCpu(x[0].addr, g.a, size)
      gpuFree(g.a)

    var x = 1.0'f32
    var y = cast[ptr UncheckedArray[float32]](omp_target_alloc(sizeof(float32)))
    useDevicePtr(y)
    #discard omp_target_memcpy_togpu(y, addr x, sizeof(float32))
    gpuMemCpyToGPU(y, addr x, sizeof(float32))
    #ompBlock("target teams"&isDevicePtr(x)):
    #ompBlock("target teams"&mapto(x)):
    ompBlock2("target teams", " map(to:", x, ")"):
      {.emit:"#pragma omp parallel".}
      {.emit:"for(int ii=0; ii<1; ii++)".}
      block:
        x = 1.0

    #macro dump(n:auto):auto =
    #  echo n.repr
    #  n
    onGpu:
      let i = getThreadNum()
      if i < n:
        c[i] = a[i] + b[i]

  test()
