import macros, sets, typetraits

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
# a copy of x, which has to be plain data.  The body keeps its text: `let x =
# getGpu(x, g)`, or the copy, under the name of the host variable, comes
# first in the kernel (CPU backend) or last in the launch, read in place by
# the kernel (device backends).  getGpu gets a stand-in of the type of x
# there (gpuArg).  Names that are types, constants, static values or
# routines stay as they are.  A name the body declares is its own, also
# where the body uses a host variable of that name elsewhere.

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
  ## Stands for the host variable in getGpu(x, g) where the name of x becomes
  ## the device value: getGpu selects on its type and must not evaluate it
  ## (the C compiler stops on the name).

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

template gpuResident*(T: typedesc): bool =
  ## Whether a kernel leaves nothing to copy back to a host variable of type
  ## T: a number, an object that owns device memory.  onGpuQueue takes such
  ## variables only.
  false
template gpuResident*(T: typedesc[SomeNumber]): bool = true

template gpuQueued*(T: typedesc; name: static string) =
  mixin gpuResident
  when not gpuResident(T):
    {.error: "onGpuQueue: '" & name & "' of type " & $T & " is copied back after a kernel; use onGpu or onGpuNowait".}

const gpuLaunchTimers* {.booldefine.} = false
  ## true: profile timers around every kernel launch; at 16^4 on one GPU
  ## they cost 1 to 2% of the time of the staggered operator

type KernelCap* = object
  body*: NimNode    ## the body, its outer variable symbols as gpuSel calls
  pre*: NimNode     ## before the launch: the device values
  shadow*: NimNode  ## the device values under the names of the host variables: locals of the launch on a device, which the kernel captures and reads in place; the first statements of the kernel with the CPU backend
  pairs*: NimNode   ## ((hostPtr(x), device value), ...) for the finalizer (hostBack)
  queued*: NimNode  ## for a launch without finalizer: stops the compile for a variable that needs one
  copies*: NimNode  ## whether the finalizer copies something back: a variable got its device value from toGpu

proc names(n: NimNode; ids: HashSet[string]; syms: seq[NimNode]): bool =
  ## whether n names something that the body declares
  case n.kind
  of nnkIdent: n.strVal in ids
  of nnkSym: n in syms
  else:
    for c in n:
      if names(c, ids, syms): return true
    false

template gpuSites*(n: SomeInteger): int =
  ## The thread count of a launch given as a number.  Here, where the
  ## launches of every backend find it: a kernel in a generic proc is
  ## compiled with the names of the module that instantiates the proc.
  int(n)

iterator gpuOnce*(t: int): int =
  ## the share of thread t of n of the indices 0..<n
  yield t

proc kernelCap*(body0: NimNode; n: NimNode = nil): KernelCap =
  ## n: the launch count of a device kernel
  var ids: HashSet[string]
  var syms: seq[NimNode]
  declNames(body0, ids, syms)
  let nr = if n != nil: n.repr else: ""
  proc oneIndex(b: NimNode): NimNode =
    ## b with its loops `for i in gpuRange(n)` over the launch count n as
    ## loops over gpuOnce: thread t of n has the one index t, and a device
    ## kernel need not find that out.  The count of a loop is n when it reads
    ## as n (templates give one of them symbols, the other identifiers) and
    ## names nothing that the body declares.
    result = b.copyNimNode
    if b.kind == nnkForStmt and b.len == 3 and b[1].kind in nnkCallKinds and b[1].len == 2 and
       b[1][0].kind in {nnkIdent, nnkSym} and b[1][0].eqIdent("gpuRange") and
       b[1][1].repr == nr and not names(b[1][1], ids, syms):
      result.add b[0]
      result.add newCall(bindSym"gpuOnce", newCall(ident"int", ident"gpuTn"))  # gpuTn: the launch's name of the thread number
      result.add oneIndex(b[2])
    else:
      for c in b: result.add oneIndex(c)
  let body = if n != nil: oneIndex(body0) else: body0
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
    of nnkPragma, nnkAccQuoted, nnkBindStmt, nnkMixinStmt: result = n  # names, no values
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
  result.queued = newStmtList()
  result.copies = newLit(false)
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
  template pre(h, g, a: untyped) =
    when declared(h): value(h, g, a)
  template shadow(l, g, a: untyped) =
    when declared(a):
      let l {.used.} = getGpu(gpuArg(a), g)
    elif declared(g):
      let l {.used.} = g
  template pair(h, g, a: untyped): untyped =
    (when declared(a): (hostPtr(h), g) else: (false, false))
  template stays(a, name: untyped) =
    when declared(a): gpuQueued(a, name)
  var k = 0
  for name in names:
    # the counter keeps the names of x and X apart: Nim compares the letters
    # after the first one without case and underscores
    let h = ident(name)
    let g = ident("gpuDev_" & name & "_" & $k)
    let a = ident("gpuHost_" & name & "_" & $k)
    inc k
    result.pre.add getAst(pre(h, g, a))
    result.shadow.add getAst(shadow(h, g, a))
    result.pairs.add getAst(pair(h, g, a))
    result.queued.add getAst(stays(a, newLit(name)))
    result.copies = infix(result.copies, "or", newCall(ident"declared", a))
  for k, h in caps:
    let s = h.strVal & "_" & $k
    let g = ident("gpuDev_" & s)
    let a = ident("gpuHost_" & s)
    result.pre.add getAst(value(h, g, a))
    result.shadow.add getAst(shadow(ident("gpuLoc_" & s), g, a))
    result.pairs.add getAst(pair(h, g, a))
    result.queued.add getAst(stays(a, newLit(h.strVal)))
    result.copies = infix(result.copies, "or", newCall(ident"declared", a))
