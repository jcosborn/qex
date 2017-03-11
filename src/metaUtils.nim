import macros
import strUtils

proc symToIdent*(x: NimNode): NimNode =
  case x.kind:
    of nnkCharLit..nnkUInt64Lit:
      result = newNimNode(x.kind)
      result.intVal = x.intVal
    of nnkFloatLit..nnkFloat64Lit:
      result = newNimNode(x.kind)
      result.floatVal = x.floatVal
    of nnkStrLit..nnkTripleStrLit:
      result = newNimNode(x.kind)
      result.strVal = x.strVal
    of nnkIdent, nnkSym:
      result = newIdentNode($x)
    of nnkOpenSymChoice:
      result = newIdentNode($x[0])
    else:
      result = newNimNode(x.kind)
      for c in x:
        result.add symToIdent(c)

macro getConst*(x: static[int]): auto =
  return newLit(x)
#macro getConst*(x: typed): auto =
  #echo x.treerepr
  #result = newLit(3)
  #result = newLit(x.intVal)

macro delayExpansion*(x:untyped):auto = result = x

macro `$`*(t:typedesc):auto =
  result = newLit(t.getType[1].repr)

macro echoType*(x:typed):auto =
  result = newEmptyNode()
  let t1 = x.getType
  echo t1.treeRepr
  echo t1.getType.treeRepr
macro echoType*(x:typedesc):auto =
  result = newEmptyNode()
  let t1 = x.getType
  echo t1.treeRepr
  echo t1[1].getType.treeRepr

macro treerep*(x:typed):auto =
  echo x.treeRepr
  newEmptyNode()

macro echoAst*(x:untyped):untyped =
  echo x.treeRepr
  x

#template dump*(x:untyped):untyped =
#  echo $(x)
#  echo astToStr(x)
#  echo repr(x)
macro dump*(x:untyped):untyped =
  let s = x[0].strVal
  #echo s
  let v = parseExpr(s)
  #echo v.treeRepr
  #echo v.toStrLit.treeRepr
  result = quote do:
    echo `x`, ": ", `v`

macro toId*(s:static[string]):untyped =
  echo s
  newIdentNode(!s)

macro toId*(s:typed):untyped =
  echo s.treeRepr
  #newIdentNode(!s)

macro toString*(id:untyped):untyped =
  #echo id.repr
  echo id.treeRepr
  if id.kind==nnkSym:
    result = newLit($id)
  else:
    result = newLit($id[0])

macro catId*(x:varargs[untyped]):auto =
  #echo x.repr
  var s = ""
  for i in 0..<x.len:
     s &= x[i].repr
  result = ident(s)

macro setType*(x:untyped; s:static[string]):auto =
  let t = ident(s)
  result = quote do:
    type `x`* = `t`

macro map*(a:tuple; f:untyped; p:varargs[untyped]):untyped =
  #echo a.treeRepr
  #echo f.treeRepr
  #echo p.treeRepr
  let ti = a.getTypeImpl
  #echo ti.treeRepr
  let nargs = ti.len
  #echo nargs
  result = newPar()
  for i in 0..<nargs:
    #let c = newCall(f,newTree(nnkBracketExpr,a,newLit(i)))
    let c = newCall(f,newDotExpr(a,ti[i][0]))
    for pp in p: c.add(pp)
    #result.add(newColonExpr(ident("field" & $i),c))
    result.add(newColonExpr(ti[i][0],c))
  #echo result.repr

macro makeCall*(op:static[string],a:tuple):untyped =
  echo op
  echo a.repr
  #echo a[0].repr
  echo a.treeRepr
  result = newCall(!op)
  let nargs = a.getType.len - 1
  for i in 0..<nargs:
    result.add(a[i][1])
  echo result.repr
  #echo result.treeRepr

#macro makeCall*(op:static[string]; a:typed):untyped =
macro makeCall*(op:static[string],a:typed,idx:typed):untyped =
  #echo op
  #echo a.repr
  #echo a.treeRepr
  #echo a.getType.treeRepr
  #echo a.getType.len
  var opid = !op
  let nargs = a.getType.len - 1
  case nargs
    of 1:
      return quote do:
        `opid`(`a`[0][`idx`])
    of 2:
      return quote do:
        `opid`(`a`[0][`idx`],`a`[1][`idx`])
    else:
      quit("makeCall: unhandled number of arguments " & $nargs)

proc evalBackticR(body:NimNode):NimNode =
  #echo body.treeRepr
  if body.kind == nnkAccQuoted:
    var id = ""
    for c in body:
      id &= $c.repr
    result = newIdentNode(id)
  else:
    result = copyNimNode(body)
    for c in body.children:
      result.add(evalBackticR(c))
  #echo result.repr

macro evalBacktic*(body:untyped):untyped =
  result = evalBackticR(body)

proc replace(id,val,body:NimNode):NimNode =
  #echo(id.treeRepr)
  #echo(id.repr)
  #echo(" " & val.treeRepr)
  #echo(" " & val.repr)
  #echo(" " & body.treeRepr)
  if body == id:
    result = val
  else:
    result = copyNimNode(body)
    for c in body.children:
      result.add(replace(id, val, c))

macro makeTyped*(x:typed):auto = x
macro makeUntyped*(x:untyped):auto = x

macro echoTyped*(x:typed):auto =
  result = newEmptyNode()
  echo x.repr
macro echoTypedTree*(x:typed):auto =
  result = newEmptyNode()
  echo x.treeRepr

macro teeTyped*(x:typed):auto =
  result = x
  echo x.repr

macro teeTypedTree*(x:typed):auto =
  result = x
  echo x.treeRepr

proc dumpTyped(r:var NimNode; x:NimNode) =
  r = quote do:
    echoTyped:
      block:
        `x`
    `r`

var idNum{.compiletime.} = 1
macro makeUnique*(x:varargs[untyped]):auto =
  result = x[^1]
  #echo result.repr
  for i in 0..(x.len-2):
    echo x[i].repr
    let v = ident(($x[i])[0..^3] & $idNum & "_")
    idNum.inc
    result = replace(x[i], v, result)
  #echo result.repr
  result.dumpTyped(result)

macro subst*(x:varargs[untyped]):auto =
  let n = x.len
  result = x[n-1]
  #echo result.repr
  for i in countup(0, n-3, 2):
    #echo x[i].repr, " ", x[i+1].repr
    var t = x[i+1]
    #echo t.repr
    if t.repr == "_":
      t = ident($x[i] & "_" & $idNum & "_")
      idNum.inc
    result = replace(x[i], t, result)
  #echo result.repr
  #echo "subst: "
  #result.dumpTyped(result)

proc separateStmtListExpr(st: var NimNode, stex: NimNode): NimNode =
  if stex.kind == nnkStmtListExpr:
    for s in 0..(stex.len-2):
      st.add stex[s]
    result = separateStmtListExpr(st, stex[^1])
  else:
    result = copyNimNode(stex)
    for s in 0..<stex.len:
      result.add separateStmtListExpr(st, stex[s])

template newLet(a,b: untyped): untyped =
  let a = b
  #mixin simpleAssign
  #var a{.noInit.}:type(b)
  #simpleAssign(a, b)

macro lets*(x:varargs[untyped]):auto =
  var prestmts = newStmtList()
  result = x[^1]
  #echo "begin lets: ", result.repr
  for i in countup(0, x.len-3, 2):
    #echo x[i].repr, " ", x[i+1].repr
    var t = separateStmtListExpr(prestmts, x[i+1])
    #echo t.repr
    if t.kind == nnkInfix:
      echo "let: ", t.repr
      let v = genSym(nskLet, $x[i])
      prestmts.add getAst(newLet(v,t))
      t = v
    result = replace(x[i], t, result)
  result = newStmtList(prestmts, result)
  #echo result.repr
  #echo "lets: "
  #result.dumpTyped(result)

macro forStaticX2(a,b:static[int]; index,body:untyped):untyped =
  #echo(index.repr)
  #echo(index.treeRepr)
  #echo(body.repr)
  #echo(body.treeRepr)
  result = newStmtList()
  for i in a..b:
    #result.add(replace(index, newIntLitNode(i), body))
    result.add(newBlockStmt(replace(index, newIntLitNode(i), body)))
  #echo(result.repr)

macro forStaticX(slice: Slice[int]; index,body: untyped): untyped =
  #echo(index.repr)
  #echo(index.treeRepr)
  #echo(slice.repr)
  #echo(slice.treeRepr)
  #echo(body.repr)
  #echo(body.treeRepr)
  result = newStmtList()
  let a = slice[1][1].intVal
  let b = slice[1][2].intVal
  for i in a..b:
    #result.add(replace(index, newIntLitNode(i), body))
    result.add(newBlockStmt(replace(index, newIntLitNode(i), body)))
  #echo(result.repr)

#template forStatic*(index,slice,body:untyped):untyped =
#  bind forStaticX
#  forStaticX(slice, index, body)

template forStatic*(index,i0,i1,body:untyped):untyped =
  bind forStaticX2
  forStaticX2(i0, i1, index, body)

template forOpt*(i,r0,r1,b:untyped):untyped =
  when compiles((const x=r0;const y=r1;x)):
    forStatic i, r0, r1:
      b
  else:
    for i in r0..r1:
      b

template depthFirst*(body:untyped; action:untyped):untyped {.dirty.} =
  proc recurse(body:NimNode):NimNode =
    #echo body.treeRepr
    result = copyNimNode(body)
    for it in body:
      #echo "it: ", it.treeRepr
      action
      result.add recurse(it)
    #echo result.repr
  result = recurse(body)
  #echo result.treeRepr
  #echo result.repr
template depthFirst2*(body:untyped; action:untyped):untyped {.dirty.} =
  proc recurse(it:var NimNode):NimNode =
    action
    result = copyNimNode(it)
    for c in it:
      var cc = c
      result.add recurse(cc)
  var b{.genSym.} = body
  result = recurse(b)
template depthFirst3*(body:untyped; action:untyped):untyped {.dirty.} =
  proc recurse(it:NimNode) =
    action
    for c in it:
      recurse(c)
  recurse(body)

macro addImportC(prefix=""; body:untyped):auto =
  #echo body.treeRepr
  let p = prefix.strVal
  depthFirst(body):
    if it.kind==nnkProcDef:
      if it.pragma.kind == nnkEmpty:
        it.pragma = newNimNode(nnkPragma)
      it.pragma.add newColonExpr(ident("importC"), newLit(p & $it.name))
macro addPragma(prg:string; body:untyped):auto =
  #echo prg.repr
  let p = parseExpr(prg.strVal)
  #echo p.treerepr
  depthFirst(body):
    if it.kind==nnkProcDef:
      if it.pragma.kind == nnkEmpty:
        it.pragma = newNimNode(nnkPragma)
      p.copyChildrenTo it.pragma
macro addReturnType(t:untyped; body:untyped):auto =
  #echo t.repr
  #echo t.treerepr
  let tt = t
  depthFirst(body):
    if it.kind==nnkProcDef:
      it[3][0] = tt
macro addArgTypes(t:varargs[untyped]; body:untyped):auto =
  #echo t.repr
  #echo t.treerepr
  let tt = t
  var a = newSeq[NimNode]()
  for i in 0..<t.len:
    a.add newIdentDefs(ident($chr(ord('a')+i)),t[i])
  depthFirst(body):
    if it.kind==nnkProcDef:
      for s in a: it[3].add s

#nnkPostfix(nnkIdent(!"*"), nnkIdent(!"hello"))

macro neverInit*(p:untyped):auto =
  #echo p.treeRepr
  result = p
  template def = {.emit:"#define memset(a,b,c)".}
  template undef = {.emit:"#undef memset".}
  insert(result.body, 0, getAst(def()))
  add(result.body, getAst(undef()))
  #echo result.treeRepr
