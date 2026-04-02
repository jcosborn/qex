import macros
import base/metaUtils
import base/omp
import backend/expr

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
proc omp_target_memcpy*(dst: pointer, src: pointer;
    length, dst_offset, src_offset: csize_t;
    dst_device_num, src_device_num: cint): cint {.omp.}
proc omp_get_default_device*: cint {.omp.}
proc omp_get_initial_device*: cint {.omp.}
proc omp_get_num_teams*: cint {.omp.}
proc omp_get_team_num*: cint {.omp.}

template omp_target_alloc*(size: SomeNumber): pointer =
  omp_target_alloc(csize_t size, omp_get_default_device())
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

template openmpDefs(body: untyped) =
  # XXX check if GC matters
  #let
  #  numTeams = omp_get_num_teams()
  #  teamNum = omp_get_team_num()
  #ompBlock("parallel num_threads(512)"):
  ompBlock("parallel"):
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

template useDevicePtr(x: auto) =
  #getrepr:
  {.emit: ["#pragma omp target data use_device_ptr(",x,")"].}

#macro mapto(x: typed): untyped =
#  let n = $x
#  result = newLit(" map(to:"&n&")")
#macro mapto(x: typed): untyped =

macro onGpuNowait*(n,b,body: untyped): auto =
  #proc deref(x,g,i:NimNode):auto = newCall("getGpu",x,g)
  proc deref(x,g,i:NimNode):auto = newCall("getGpu",x,newTree(nnkAccQuoted,g,ident"xx"))
  template target(cpuPrepare, cpuFinalize, devicePtrDeclare, body: untyped) =
    mixin toGpu, getGpu, fromGpu
    {.push checks: off.}
    {.push stacktrace: off.}
    block:
      cpuPrepare  # a let section declare and save device pointers
      #proc gpuProc {.gensym.} =
      threadSingle:
        #ompBlock2("target teams num_teams(1024)", devicePtrDeclare):
        ompBlock2("target teams", devicePtrDeclare):
          openmpDefs:
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
  result = getast(target(cpuPrepare, cpuFinalize, isDevicePtrs, body))
  echo result.repr

var gpuNumThreadsRequest* = 0
var gpuBlockSizeRequest* = 0
template onGpuNowait*(body: untyped): auto =
  onGpuNoWait(gpuNumThreadsRequest, gpuBlockSizeRequest, body)
template onGpuNowait*(n,body: untyped): auto =
  var b = gpuBlockSizeRequest
  while b > n: b = b div 2
  onGpuNoWait(n, b, body)
#template onGpuNowait*(n,b,body: untyped): auto =
#  onGpuNoWait(n, b, body)

template onGpu*(body: untyped) =
  let finalize = onGpuNoWait(body)
  finalize()
template onGpu*(n,body: untyped) =
  let finalize = onGpuNoWait(n, body)
  finalize()
template onGpu*(n,b,body: untyped) =
  let finalize = onGpuNoWait(n, b, body)
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
