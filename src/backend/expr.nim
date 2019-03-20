import macros

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
        return newCall(a,ii)
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
  let e = getVars(v, x, a)
  var p = newStmtList()
  for vs in v:
    p.add newCall(f,vs)
  result = p
  #echo result.treerepr

macro packVars*(x: untyped, f: untyped): auto =
  #echo x.treerepr
  var v = newSeq[NimNode](0)
  let a = ident("foo")
  let e = getVars(v, x, a)
  var p = newPar()
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
