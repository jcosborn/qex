import macros, sets, typetraits

var ignore {.compileTime.}: seq[NimNode]
proc addIfNewSym(s: var seq[NimNode], x: NimNode): int =
  let sx = $x
  for i in 0..<ignore.len:
    if ignore[i].eqIdent sx: return -1
  for i in 0..<s.len:
    if s[i].eqIdent sx: return i
  result = s.len
  s.add x

proc cpNimNode(x: NimNode): NimNode =
  result = newNimNode(x.kind)
  case x.kind
  of nnkCharLit..nnkUInt64Lit:
    result.intVal = x.intVal
  of nnkFloatLit..nnkFloat64Lit:
    result.floatVal = x.floatVal
  of nnkStrLit..nnkTripleStrLit:
    result.strVal = x.strVal
  of nnkIdent:
    #result.ident = ident(x.repr)
    result = newIdentNode($x)
  of {nnkSym,nnkOpenSymChoice}:
    #echo "got sym"
    #quit -1
    #result = newIdentNode($x)
    result = x.copy
  else:
   discard

proc getVars*(v: var seq[NimNode], x,a: NimNode): NimNode =
  proc recurse(it: NimNode, vars: var seq[NimNode], a: NimNode): NimNode =
    var r0 = 0
    var r1 = it.len - 1
    case it.kind
    of {nnkSym, nnkIdent}:
      let i = vars.addIfNewSym(it)
      if i>=0:
        let ii = newLit(i)
        return newCall(a,it,ii)
    of nnkCallKinds: r0 = 1
    of nnkDotExpr: r1 = 0
    of {nnkVarSection,nnkLetSection}:
      result = it.cpNimNode
      for c in it:
        result.add c.cpNimNode
        for i in 0..(c.len-3):
          ignore.add c[i]
          result[^1].add c[i].cpNimNode
        result[^1].add c[^2].cpNimNode
        result[^1].add recurse(c[^1], vars, a)
      return
    else: discard
      #echo it.treerepr
    result = it.cpNimNode
    for i in 0..<r0:
      result.add it[i].cpNimNode
    for i in r0..r1:
      result.add recurse(it[i], vars, a)
    for i in (r1+1)..<it.len:
      result.add it[i].cpNimNode
  ignore.newSeq(0)
  result = recurse(x, v, a)

macro packVarsStmt*(x: untyped, f: untyped): auto =
  #echo x.treerepr
  var v = newSeq[NimNode](0)
  let a = ident("foo")
  discard getVars(v, x, a)
  var p = newStmtList()
  for vs in v:
    p.add newCall(f,vs)
  result = p
  #echo result.treerepr

macro packVars*(x: untyped, f: untyped): auto =
  #echo x.treerepr
  var v = newSeq[NimNode](0)
  let a = ident("foo")
  discard getVars(v, x, a)
  var p = newNimNode(nnkTupleConstr)
  if v.len==0:
    p.add newNimNode(nnkExprColonExpr).add(ident("Field0"),newLit(1))
  elif v.len==1:
    let vi = ident($v[0])
    p.add newNimNode(nnkExprColonExpr).add(ident("Field0"),newCall(f,vi))
  else:
    for vs in v:
      p.add newCall(f,vs)
  result = p
  #echo result.treerepr

macro substVars*(x: untyped, a: untyped): auto =
  #echo x.treerepr
  var v = newSeq[NimNode](0)
  let e = getVars(v, x, a)
  result = e
  #echo result.treerepr

# A kernel launch (onGpu, onGpuNowait) belongs in a proc: the device backends
# put OpenMP pragmas around it, and Nim emits those at file scope for code at
# the top level of a module, where the device compiler stops.
#
# The finalizer of a kernel launch gets its variables as arguments: for each one
# the address of the host variable, its value where it has no address (a call
# that gives a ref, a constant), or false for a type, and the device handle.
# It holds the addresses, so it has to run while the host variables are alive.

type HostVal*[T] = object
  ## a kernel variable without an address
  v*: T

template hostPtr*(x: typedesc): bool = false
template hostPtr*(x: typed): auto =
  when compiles(addr x): addr x else: HostVal[typeof(x)](v: x)

template hostFromGpu*(h: bool; g: untyped) = discard
template hostFromGpu*[T](h: ptr T; g: untyped) =
  mixin fromGpu
  fromGpu(h[], g)
template hostFromGpu*[T](h: HostVal[T]; g: untyped) =
  mixin fromGpu
  fromGpu(h.v, g)

template hostBack*(hp: untyped): untyped =
  ## A proc that calls fromGpu for the pairs of hp (genHostPairs).  It
  ## expands at the launch, so it finds the fromGpu declared there, also inside
  ## a proc; nimcall keeps it from capturing the locals of the launch, which
  ## the threads of a region share.
  proc back(ps: typeof(hp)) {.gensym, nimcall.} =
    for p in fields(ps): hostFromGpu(p[0], p[1])
  back

var gpuLaunchThreads* = 0
  ## The threads of the kernel being launched, for the toGpu of objects that
  ## hold a value per thread (GpuSum); written by the launching thread.

# Capture of the variables of a kernel body.
#
# Every outer variable that a kernel body names becomes a value of the
# launch: toGpu(x) where x has one, with fromGpu(x, g) owed afterwards, else
# a copy of x, which has to be plain data.  The body keeps its text: the
# kernel starts with `let x = getGpu(x, g)`, or the copy, under the name of
# the host variable.  getGpu gets a stand-in of the type of x there (gpuArg).
# Names that are types, constants, static values or routines stay as they
# are.  A name the body declares is its own, also where the body uses a
# host variable of that name elsewhere.  gpuObjects is declared after the
# device values when one of them is no number and no pointer.

macro gpuIsVar*(x: typed): untyped =
  ## Whether x is a variable or the value of an expression, which a kernel
  ## captures: no type, constant, static value or routine, and no variable
  ## of the threads or of the C compiler (threadvar, importc).
  var n = x
  while n.kind in {nnkHiddenDeref, nnkHiddenAddr, nnkHiddenStdConv, nnkHiddenSubConv}: n = n[^1]
  var r: bool
  if x.getType.typeKind == ntyTypeDesc: r = false
  else:
    case n.kind
    of nnkSym:
      r = n.symKind in {nskVar, nskLet, nskParam, nskResult, nskForVar, nskTemp}
      if r and n.symKind == nskVar:
        let d = n.getImpl
        if d.kind == nnkIdentDefs and d[0].kind == nnkPragmaExpr:
          for p in d[0][1]:
            let q = if p.kind == nnkExprColonExpr: p[0] else: p
            if q.kind in {nnkIdent, nnkSym} and (q.eqIdent"threadvar" or q.eqIdent"importc" or q.eqIdent"importcpp"): r = false
    of nnkCharLit..nnkNilLit: r = false
    else: r = true
  newLit(r)

proc gpuArg*[T](t: typedesc[T]): T {.importc: "QEX_HOST_VALUE_IN_KERNEL", nodecl.}
  ## Stands for the host variable in getGpu(x, g) inside a kernel, where the
  ## name of x is the device value: getGpu selects on its type and must not
  ## evaluate it (the C compiler stops on the name).

template gpuSel*(h, l: untyped): untyped =
  ## in a kernel: the device value l of the host variable h where the launch
  ## made one
  when declared(l): l else: h

proc declName(d: NimNode; ids: var HashSet[string]; syms: var seq[NimNode]) =
  case d.kind
  of nnkIdent: ids.incl d.strVal
  of nnkSym: syms.add d
  of nnkPragmaExpr: declName(d[0], ids, syms)
  of nnkPostfix: declName(d[1], ids, syms)
  of nnkAccQuoted:
    if d.len == 1: declName(d[0], ids, syms)
  else: discard

proc declNames(n: NimNode; ids: var HashSet[string]; syms: var seq[NimNode]) =
  ## the names n declares, at any depth
  case n.kind
  of nnkVarSection, nnkLetSection, nnkConstSection:
    for d in n:
      if d.kind in {nnkVarTuple, nnkIdentDefs, nnkConstDef}:
        for i in 0..<d.len-2: declName(d[i], ids, syms)
  of nnkTypeSection:
    for d in n: declName(d[0], ids, syms)
  of nnkForStmt:
    for i in 0..<n.len-2:
      if n[i].kind == nnkVarTuple:
        for c in n[i]: declName(c, ids, syms)
      else: declName(n[i], ids, syms)
  of RoutineNodes:
    declName(n[0], ids, syms)
    if n[3].kind == nnkFormalParams:
      for i in 1..<n[3].len:
        for j in 0..<n[3][i].len-2: declName(n[3][i][j], ids, syms)
  else: discard
  for c in n: declNames(c, ids, syms)

type KernelCap* = object
  body*: NimNode    ## the body, its outer variable symbols as gpuSel calls
  pre*: NimNode     ## before the launch: the device values
  shadow*: NimNode  ## first in the kernel: the device values under the names of the host variables
  pairs*: NimNode   ## ((hostPtr(x), device value), ...) for the finalizer (hostBack)

proc kernelCap*(body: NimNode): KernelCap =
  var ids: HashSet[string]
  var syms: seq[NimNode]
  declNames(body, ids, syms)
  var names: OrderedSet[string]  # identifiers, captured by name
  var caps: seq[NimNode]  # outer symbols (generic procs and templates bind them)
  proc walk(n: NimNode): NimNode =
    case n.kind
    of nnkIdent:
      if n.strVal notin ids: names.incl n.strVal
      result = n
    of nnkSym:
      if n.symKind in {nskVar, nskLet, nskParam, nskResult, nskForVar, nskTemp, nskTemplate} and n notin syms:
        var k = caps.find(n)
        if k < 0:
          k = caps.len
          caps.add n
        result = newCall(bindSym"gpuSel", n, ident("gpuLoc_" & n.strVal & "_" & $k))
      else: result = n
    of nnkDotExpr:  # x.f: f is a field or a routine
      result = n.copyNimNode
      result.add walk(n[0])
      result.add n[1]
    of nnkExprEqExpr, nnkExprColonExpr:  # named argument, field of a constructor
      result = n.copyNimNode
      result.add n[0]
      result.add walk(n[1])
    of nnkPragma, nnkAccQuoted: result = n
    of nnkCallKinds:
      result = n.copyNimNode
      if n[0].kind in {nnkIdent, nnkSym, nnkOpenSymChoice, nnkClosedSymChoice}: result.add n[0]  # a routine
      else: result.add walk(n[0])
      for i in 1..<n.len: result.add walk(n[i])
    else:
      result = n.copyNimNode
      for c in n: result.add walk(c)
  result.body = walk(body)
  result.pre = newStmtList()
  result.shadow = newStmtList()
  result.pairs = newNimNode(nnkTupleConstr)
  # g: the device value, a: the host type where toGpu gave g, l: g in the kernel
  template value(h, g, a: untyped) =
    when compiles(gpuIsVar(h)):
      when gpuIsVar(h):
        when compiles(toGpu(h)):
          let g {.used.} = toGpu(h)
          type a {.used.} = typeof(h)
        elif supportsCopyMem(typeof(h)):
          let g {.used.} = h
        else:
          {.error: "kernel variable '" & astToStr(h) & "' of type " & $typeof(h) & " has no toGpu and is not plain data".}
        when typeof(g) isnot (SomeNumber | bool | enum | ptr | pointer) and not declared(gpuObjects):
          const gpuObjects {.inject, used.} = true  # an object among the device values
  template pre(h, g, a: untyped) =
    when declared(h): value(h, g, a)
  template shadow(l, g, a: untyped) =
    when declared(a):
      let l {.used.} = getGpu(gpuArg(a), g)
    elif declared(g):
      let l {.used.} = g
  template pair(h, g, a: untyped): untyped =
    (when declared(a): (hostPtr(h), g) else: (false, false))
  for name in names:
    let h = ident(name)
    let g = ident("gpuDev_" & name)
    let a = ident("gpuHost_" & name)
    result.pre.add getAst(pre(h, g, a))
    result.shadow.add getAst(shadow(h, g, a))
    result.pairs.add getAst(pair(h, g, a))
  for k, h in caps:
    let s = h.strVal & "_" & $k
    let g = ident("gpuDev_" & s)
    let a = ident("gpuHost_" & s)
    result.pre.add getAst(value(h, g, a))
    result.shadow.add getAst(shadow(ident("gpuLoc_" & s), g, a))
    result.pairs.add getAst(pair(h, g, a))

when isMainModule:
  template test(x) =
    template getref(t: untyped): untyped = addr(t)
    let v = packVars(x,getref)
    proc foo(xx: type(v)) =
      template deref(i: int): untyped = xx[i][]
      substVars(x, deref)
    foo(v)

  macro dump(x: typed): auto =
    echo x.repr
    x

  var x,y,z: float

  dump:
    test:
      x = 1
      y = 2
      z = x + y
  echo x, y, z
