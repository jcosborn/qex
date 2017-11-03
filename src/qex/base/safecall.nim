import macros

macro dumpSym*(x: typed): auto =
  result = newEmptyNode()
  echo x.treerepr
  if x.kind==nnkSym:
    echo x.symbol.getImpl.treerepr

var scArgCount {.compileTime.} = 0
proc checkarg(arg: NimNode): tuple[sl:NimNode,a:NimNode] =
  result.sl = newStmtList()
  case arg.kind
  of AtomicNodes:
    result.a = arg
  of {nnkStmtList,nnkStmtListExpr}:
    for i in 0..(arg.len-2):
      result.sl.add arg[i]
    let r = checkarg(arg[^1])
    r.sl.copyChildrenTo(result.sl)
    result.a = r.a
  of {nnkObjConstr}:
    let a = arg.copy
    for i in 1..<arg.len:
      let r = checkarg(arg[i][1])
      a[i][1] = r.a
      r.sl.copyChildrenTo(result.sl)
    result.a = a
  of {nnkBracketExpr,nnkDotExpr,nnkHiddenDeref}:
    let a = arg.copy
    for i in 0..<arg.len:
      let r = checkarg(arg[i])
      a[i] = r.a
      r.sl.copyChildrenTo(result.sl)
    result.a = a
  else:
    #echo arg.kind
    let n = genSym(nskLet, "scArg" & $scArgCount)
    inc scArgCount
    result.sl.add newLetStmt(n, arg)
    result.a = n

#macro safecall2*(dbg: static[bool], args: tuple, fn: untyped): auto =
macro safecall2*(dbg: static[bool], args: tuple, fn: varargs[untyped]): auto =
  #echo fn.treerepr
  if dbg: echo args.treerepr
  var sl = newStmtList()
  #var call = newCall(fn)
  var call = newCall(fn[0])
  for i in 0..<args.len:
    let (sl1, arg1) = checkarg(args[i][1])
    if sl1.len>0:
      sl1.copyChildrenTo(sl)
    call.add arg1
  if sl.len>0:
    result = sl
    result.add call
  else:
    result = call
  if dbg: echo result.treerepr
  #echo result.treerepr
macro safecallO*(fn: untyped, args: varargs[untyped]): auto =
  var t = newPar()
  for i in 0..<args.len:
    t.add newColonExpr(ident("Field"& $i), args[i])
  result = newCall(bindSym"safecall2", ident("false"))
  result.add t
  result.add fn
  #echo fn.treerepr
  #echo args.treerepr
macro safecall*(args: varargs[untyped]): auto =
  var t = newPar()
  for i in 1..<args.len:
    t.add newColonExpr(ident("Field"& $i), args[i])
  result = newCall(bindSym"safecall2", ident("false"))
  result.add t
  result.add args[0]
  #echo args[0].treerepr
macro safecallD*(fn: untyped, args: varargs[untyped]): auto =
  var t = newPar()
  for i in 0..<args.len:
    t.add newColonExpr(ident("Field"& $i), args[i])
  result = newCall(bindSym"safecall2", ident("true"))
  result.add t
  result.add fn


when isMainModule:
  macro dumpTyped*(x: typed): auto =
    result = newEmptyNode()
    echo x.treerepr

  template foo1(x: untyped): untyped =
    echo x+x

  block:
    dumpTyped:
      safecall(foo1):
        let x = 2
        x

  block:
    dumpTyped:
      safecall(foo1):
        let x = 2
        x+x

  block:
    dumpTyped:
      safecall(foo1):
        let x = 2
        float(x+x)

  proc `+`[N,T](x,y: array[N,T]): auto = x
  proc `$`[N,T](x: array[N,T]): string = $x[0]
  block:
    dumpTyped:
      safecall(foo1):
        let x = 2
        [x+1,x+2]
