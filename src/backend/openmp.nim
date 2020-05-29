import macros
import base/metaUtils
import base/omp

{.pragma: omp, header:"omp.h".}
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
  if byte notin {1,2,4,8,16,32,64,128,256}:
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

template omp_target_alloc*(size: csize_t): pointer =
  omp_target_alloc(size, omp_get_default_device())
template omp_target_memcpy_tocpu*(dst: pointer, src: pointer; length: csize_t): cint =
  omp_target_memcpy(dst, src, length, 0, 0, omp_get_initial_device(), omp_get_default_device())
template omp_target_memcpy_togpu*(dst: pointer, src: pointer; length: csize_t): cint =
  omp_target_memcpy(dst, src, length, 0, 0, omp_get_default_device(), omp_get_initial_device())
template omp_target_free*(device_ptr: pointer) =
  omp_target_free(device_ptr, omp_get_default_device())

template gpuMalloc*(size:csize_t):pointer = omp_target_alloc(size)
template gpuFree*(device_ptr:pointer) = omp_target_free(device_ptr)
template gpuMemCpyToCPU*(dst: pointer, src: pointer; length: csize_t): cint =
  omp_target_memcpy_tocpu(dst, src, length)
template gpuMemCpyToGPU*(dst: pointer, src: pointer; length: csize_t): cint =
  omp_target_memcpy_togpu(dst, src, length)

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

#template openmpDefs(body: untyped): untyped {.dirty.} =
template openmpDefs(body: untyped): untyped =
  # XXX check if GC matters
  #bind omp_get_num_teams, omp_get_team_num, omp_get_num_threads, omp_get_thread_num
  #bind ompBlock, inlineProcs
  let
    numTeams = omp_get_num_teams()
    teamNum = omp_get_team_num()
  ompBlock("parallel"):
    let
      numThreads = omp_get_num_threads()
      threadNum = omp_get_thread_num()
    template getThreadNum: untyped {.used.} = teamNum.int * numThreads.int + threadNum.int
    template getNumThreads: untyped {.used.} = numTeams.int * numThreads.int
    {.emit:"#define nimZeroMem(b,len) memset((b),0,(len))".}
    inlineProcs:
      body
    {.emit:"#undef nimZeroMem".}

proc prepareVars(n:NimNode):seq[NimNode] =
  # get a list of vars and new symbols to replace them, using let binding for now XXX
  #     <- [(id, varsym, letptrsym), ...]
  # the symbols in n is changed
  #echo "### prepareVars: ",n.treerepr
  var ignoreStack = newseq[NimNode]()
  var openvars = newseq[NimNode]()
  proc go(n:NimNode) =
    # ign is a stack for ignoring lexical bindings: [(outer,...), (inner,...), ...]
    #echo "go get: ",n.repr
    #block:
    #  var ignstr = ""
    #  for c in ignoreStack: ignstr &= ("\n" & c.repr)
    #  echo "ign has: ",ignstr
    var newscope = false
    if n.kind in {nnkBlockStmt, nnkBlockExpr, nnkIfExpr, nnkElifExpr, nnkElseExpr,
        nnkIfStmt, nnkElifBranch, nnkElse, nnkCaseStmt, nnkOfBranch,
        nnkWhileStmt, nnkForStmt} + RoutineNodes:
      # New lexical scope
      newscope = true
      ignoreStack.add newPar()
    for i in 0..<n.len:
      #echo "### ",n[i].lisprepr
      case n[i].kind
      of {nnkVarSection,nnkLetSection}:
        for cc in n[i]:
          for c in 0..cc.len-2:
            ignoreStack[^1].add cc[c]
      of nnkOpenSymChoice:
        if n.kind in Callnodes: continue
      of Callnodes:
        if n[i][0].kind in {nnkSym, nnkIdent}:
          var newid = true
          for c in ignoreStack[0]:
            if c == n[i][0]:
              newid = false
              break
          if newid:
            ignoreStack[0].add n[i][0]
      of {nnkSym, nnkIdent}:
        if n.kind == nnkDotExpr and i > 0: continue
        var ignore = false
        for cc in ignoreStack:
          for c in cc:
            if c.eqIdent n[i]:
              ignore = true
              break
          if ignore: break
        if not ignore:
          var newvar = true
          for c in openvars:
            if c[0].eqIdent n[i]:
              n[i] = newcall("gpuVarPtr",c[1],c[2])
              newvar = false
              break
          #echo "EXPR: ",n.lisprepr
          #echo "ID:   ",n[i].repr,"  newvar: ",newvar.repr
          #var rs = ""
          #for c in openvars:
          #  rs &= "  " & c.repr
          #echo "RES:  ",rs
          if newvar:
            let nv = gensym(nskvar, "gpu_var_" & $n[i])
            let np = gensym(nsklet, "gpu_ptr_" & $n[i])
            ignoreStack[0].add nv
            ignoreStack[0].add np
            openvars.add newpar(n[i], nv, np)
            n[i] = newcall("gpuVarPtr",nv,np)
      else:
        discard
      n[i].go
    if newscope: ignoreStack.setLen(ignoreStack.len-1)
  ignoreStack.add newPar(ident"gpuVarPtr")
  n.go
  openvars
type OffloadDummy*[T] = object
proc genCpuPrepare(n:seq[NimNode]):NimNode =
  template r(x,v,p:untyped):untyped =
    mixin offloadUsePtr, offloadUseVar, offloadPtr, offloadVar
    when offloadUsePtr(x):
      let p = offloadPtr(x)
    else:
      let p = cast[pointer](0)
    when offloadUseVar(x):
      var v = offloadVar(x,p)
    else:
      var v{.noinit.}:OffloadDummy[typeof(x)]
  result = newstmtlist()
  for c in n:
    result.add getast r(c[0],c[1],c[2])
proc genGpuPrepare(n:seq[NimNode]):NimNode =
  template r(x,v,p:untyped):untyped =
    mixin gpuPrepareOffload, rungpuPrepareOffload
    when rungpuPrepareOffload(x): gpuPrepareOffload(v,p)
  result = newstmtlist()
  for c in n:
    result.add getast r(c[0],c[1],c[2])
proc genCpuFinalize(n:seq[NimNode]):NimNode =
  template r(x,v,p:untyped):untyped =
    mixin cpuFinalizeOffload, runcpuFinalizeOffload
    when runcpuFinalizeOffload(x): cpuFinalizeOffload(x,v,p)
  result = newstmtlist()
  for c in n:
    result.add getast r(c[0],c[1],c[2])
proc declarePtrString(n:seq[NimNode]):NimNode =
  template res(ptrlist:untyped):untyped =
    const s = ptrlist
    when s.len == 0: "" else: "is_device_ptr(" & s[0..^2] & ")"
  template varname(x, xp:untyped):untyped =
    mixin offloadPtr
    when compiles(offloadPtr(x)): xp&"," else: ""
  var ps = newlit""
  for c in n:
    ps = infix(getast varname(c[0], $c[2]), "&", ps)
  result = getast res(ps)

macro onGpu*(body: untyped): untyped =
  # the architecture for cpugpuarray requires us replace body before it gets expanded, so we require untyped.
  template target(cpuPrepare, gpuPrepare, cpuFinalize, devicePtrDeclare, body: untyped): untyped =
    mixin hasGpuPtr, requireGpuMem
    {.push checks: off.}
    {.push stacktrace: off.}
    proc gpuProc {.gensym.} =
      cpuPrepare  # a let section declare and save device pointers
      const isDevicePtrList = devicePtrDeclare  # is_device_ptr(ptrList) in string
      ompBlock("target teams " & isDevicePtrList):
        openmpDefs:
          gpuPrepare
          body
      cpuFinalize
    gpuProc()
  let
    v = prepareVars(body)  # gather gpu pointers in symbols, body is changed accordingly
    cpuPrepare = genCpuPrepare v
    gpuPrepare = genGpuPrepare v
    cpuFinalize = genCpuFinalize v
    isDevicePtrs = declarePtrString v
  result = getast(target(cpuPrepare, gpuPrepare, cpuFinalize, isDevicePtrs, body))
  #echo result.repr

# XXX fix the following
template onGpu*(totalNumThreads, body: untyped): untyped = onGpu(body)
template onGpu*(totalNumThreads, numThreadsPerTeam, body: untyped): untyped = onGpu(body)

#[
template onGpu*(nn,tpb: untyped, body: untyped): untyped =
  block:
    var v = packVars(body, getGpuPtr)
    type ByCopy[T] {.bycopy.} = object
      d: T
    proc kern(xx: ByCopy[type(v)]) {.cudaGlobal.} =
      template deref(k: int): untyped = xx.d[k]
      substVars(body, deref)
    let ni = nn.int32
    let threadsPerBlock = tpb.int32
    let blocksPerGrid = (ni+threadsPerBlock-1) div threadsPerBlock
    #echo "launching kernel"
    cudaLaunch(kern, blocksPerGrid, threadsPerBlock, v)
    discard cudaDeviceSynchronize()
template onGpu*(nn: untyped, body: untyped): untyped = onGpu(nn, 64, body)
template onGpu*(body: untyped): untyped = onGpu(512*64, 64, body)
]#

template offloadUseVar*(x:SomeNumber):bool = true
template offloadUsePtr*(x:SomeNumber):bool = false
template rungpuPrepareOffload*(x:SomeNumber):bool = false
template runcpuFinalizeOffload*(x:SomeNumber):bool = false
template gpuVarPtr*(v:SomeNumber,p:untyped):untyped = v
template offloadVar*(x:SomeNumber,p:untyped):untyped = x

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
proc isIndex(n,i:NimNode):bool =
  result = n.eqident i
  if n.kind == nnkHiddenStdConv:
    result = n[1].eqident i
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
            ptrs.add newPar(v, n)
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

  proc test =
    var n = 50000.cint
    var
      a = newSeq[float32](n)
      b = newSeq[float32](n)
      c = newSeq[float32](n)

    template `[]`(x: FltArr, i: SomeInteger): untyped = x.a[][i]
    template `[]=`(x: FltArr, i: SomeInteger, y:untyped): untyped = x.a[][i] = y

    template offloadUseVar(x:seq):bool = true
    template offloadUsePtr(x:seq):bool = true
    template rungpuPrepareOffload(x:seq):bool = true
    template runcpuFinalizeOffload(x:seq):bool = true
    template gpuVarPtr(v:FltArr,p:untyped):untyped = v
    template offloadPtr(x:seq):untyped =
      let size = x.len * sizeof(x[0])
      let xp = omp_target_alloc(size)
      discard omp_target_memcpy_togpu(xp, x[0].addr, size)
      cast[ptr UncheckedArray[type(x[0])]](xp)
    template offloadVar(x:seq,p:untyped):untyped = FltArr(a:p)
    template gpuPrepareOffload(v:FltArr,p:untyped):untyped = v.a=p
    template cpuFinalizeOffload(x:seq,v,p:untyped):untyped = omp_target_free(p)

    macro dump(n:typed):typed =
      echo n.repr
      n
    dump:
      onGpu:
        let i = getThreadNum()
        if i < n:
          c[i] = a[i] + b[i]

  test()
