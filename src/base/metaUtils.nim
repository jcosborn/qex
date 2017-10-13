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
macro getConst*(x: static[int64]): auto =
  return newLit(x)
#macro getConst*(x: typed): auto =
  #echo x.treerepr
  #result = newLit(3)
  #result = newLit(x.intVal)

macro makeIdent*(x: untyped): untyped =
  result = symToIdent(x)

macro delayExpansion*(x:untyped):auto = result = x

macro `$`*(t:typedesc):auto =
  result = newLit(t.getType[1].repr)

macro echoType*(x:typed):auto =
  result = newEmptyNode()
  echo x.getTypeInst.repr
  echo x.getTypeImpl.repr
macro echoType*(x:typedesc):auto =
  result = newEmptyNode()
  let t1 = x.getType
  echo t1.repr
  echo t1[1].getType.repr
macro echoTypeTree*(x:typed):auto =
  result = newEmptyNode()
  echo x.getTypeInst.treeRepr
  echo x.getTypeImpl.treeRepr
macro echoTypeTree*(x:typedesc):auto =
  result = newEmptyNode()
  let t1 = x.getType
  echo t1.treeRepr
  echo t1[1].getType.treeRepr

macro treerep*(x:typed):auto =
  echo x.lineinfo
  echo x.treeRepr
  newEmptyNode()

macro echoAst*(x:untyped):untyped =
  echo x.lineinfo
  echo x.treeRepr
  x

macro echoRepr*(x: untyped): untyped =
  echo x.lineinfo
  echo x.repr
  newEmptyNode()

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

when true:
  template staticTraceBegin*(x: untyped) = discard
  template staticTraceEnd*(x: untyped) = discard
  template staticTraceReturn*(x,body: untyped): untyped = body
else:
  var staticTraceLevel {.compiletime.} = 0
  macro staticTraceBegin*(x: untyped): untyped =
    #echo x.treerepr
    let s = x[0].repr
    echo repeat("  ", staticTraceLevel) & "-> " & s
    inc staticTraceLevel
    newEmptyNode()
  macro staticTraceEnd*(x: untyped): untyped =
    let s = x[0].repr
    dec staticTraceLevel
    echo repeat("  ", staticTraceLevel) & "<- " & s
    newEmptyNode()
  template staticTraceReturn*(x,body: untyped): untyped =
    staticTraceBegin: x
    let str = body
    staticTraceEnd: x
    str

macro toId*(s:static[string]):untyped =
  echo s
  newIdentNode(!s)

macro toId*(s:typed):untyped =
  echo s.treeRepr
  #newIdentNode(!s)

macro toString*(id: untyped): untyped =
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

#[
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
]#

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

proc replaceConv(id,val,body:NimNode):NimNode =
  if body == id:
    result = val
  elif body.kind == nnkHiddenStdConv and body[1] == id:
    result = val
  else:
    result = copyNimNode(body)
    for c in body.children:
      result.add(replaceConv(id, val, c))

macro makeTyped*(x:typed):auto = x
macro makeUntyped*(x:untyped):auto = x

macro echoUntyped*(x: untyped): auto =
  result = newEmptyNode()
  echo x.repr
macro echoUntypedTree*(x: untyped): auto =
  result = newEmptyNode()
  echo x.treeRepr
macro echoTyped*(x: typed): auto =
  result = newEmptyNode()
  echo x.repr
macro echoTypedTree*(x: typed): auto =
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

macro subst*(x: varargs[untyped]): untyped =
  let n = x.len
  result = x[n-1]
  #echo result.repr
  for i in countup(0, n-3, 2):
    #echo x[i].repr, " ", x[i+1].repr
    var t = x[i+1]
    #echo t.repr
    if t.repr == "_":
      t = ident($x[i] & "_" & $idNum & "_")
      inc idNum
    result = replace(x[i], t, result)
  #echo result.treerepr
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

#macro forStaticX2(a,b:static[int]; index,body:untyped):untyped =
macro forStaticX3(a0,b0: typed; index,body: untyped): untyped =
  #echo(index.repr)
  #echo(index.treeRepr)
  #echo(body.repr)
  #echo(body.treeRepr)
  #echo a0.treerepr
  #echo b0.treerepr
  let a = a0.intVal
  let b = b0.intVal
  result = newStmtList()
  for i in a..b:
    #result.add(replace(index, newIntLitNode(i), body))
    result.add(newBlockStmt(replace(index, newIntLitNode(i), body)))
  #echo(result.repr)

template forStaticX2(a0,b0: typed; index,body: untyped): untyped =
  forStaticX3(getConst(a0), getConst(b0), index, body)

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

macro forStaticU*(a0,b0: typed; fn: untyped): untyped =
  #echo a0.treerepr
  #echo b0.treerepr
  let a = a0.intVal
  let b = b0.intVal
  result = newStmtList()
  for i in a..b:
    result.add newBlockStmt(newCall(fn,newIntLitNode(i)))
  #echo result.repr

#template forStatic*(index,slice,body:untyped):untyped =
#  bind forStaticX
#  forStaticX(slice, index, body)

template forStaticUntyped*(index,i0,i1,body:untyped):untyped =
# template forStatic*(index,i0,i1,body: untyped): untyped {.dirty.}=
  bind forStaticX2
  forStaticX2(i0, i1, index, body)

proc unrollFor*(n:NimNode):NimNode =
  template must(p:bool) =
    if not p:
      echo "unrollFor can't handle it:"
      echo n.repr
      if n.len > 1:
        for c in n[1]:
          echo c.lineinfo," :: ",c.lisprepr
      error "Abort compilation."
  #echo n.treerepr
  must: n.kind == nnkForStmt
  must: n.len == 3
  must: n[1].kind == nnkInfix
  must: n[1].len == 3
  must: n[1][0].eqident ".."
  must: n[1][1].kind in nnkCharLit..nnkUInt64Lit
  must: n[1][2].kind in nnkCharLit..nnkUInt64Lit
  let
    a = n[1][1].intval
    b = n[1][2].intval
  result = newStmtList()
  for i in a..b:
    result.add newNimNode(nnkBlockStmt, n).add(
      ident("ITR: " & $i & " :: " & n.repr), replaceConv(n[0], newIntLitNode(i), n[2]))
  #echo result.treerepr
macro unrollFor*(n:typed):untyped =
  if n.kind == nnkForStmt: return n.unrollFor
  result = newstmtlist()
  for c in n:
    if c.kind == nnkForStmt: result.add c.unrollFor
    else: result.add c

template forStaticUnRollFor*(index,i0,i1,body:untyped):untyped =
  unrollFor:
    for index in i0..i1: body

template forStatic*(index,i0,i1,body:untyped):untyped =
  forStaticUntyped(index,i0,i1,body)
  # forStaticUnRollFor(index,i0,i1,body)

template forOpt*(i,r0,r1,b:untyped):untyped =
  when compiles((const x=r0;const y=r1;x;y)):
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
  #let tt = t
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

proc normalizeAstR(a: NimNode): NimNode =
  result = a
  case result.kind
  of {nnkStmtList,nnkStmtListExpr}:
    var nonempty,last = 0
    for i in 0..<result.len:
      result[i] = normalizeAstR(result[i])
      if result[i].kind notin {nnkEmpty,nnkDiscardStmt}:
        inc nonempty
        last = i
    case nonempty
    of 0: result = newEmptyNode()
    of 1: result = result[last]
    else: discard
  else:
    discard

macro normalizeAst*(a: typed): untyped =
  result = normalizeAstR(a)
  #echo "normalizeAst"
  #echo result.treerepr

proc optimizeAstR(a: NimNode): NimNode =
  result = a
  case result.kind
  of {nnkStmtList,nnkStmtListExpr}:
    var nonempty,last = 0
    for i in 0..<result.len:
      result[i] = normalizeAstR(result[i])
      if result[i].kind notin {nnkEmpty,nnkDiscardStmt}:
        inc nonempty
        last = i
    case nonempty
    of 0: result = newEmptyNode()
    of 1: result = result[last]
    else: discard
  else:
    discard

#proc optStmtList(x: NimNode, sym: var seq[NimSym],
#                 repl,stmts: var seq[NimNode]): NimNode =
#  result = x.copyNimNode
#  for i in 0..<x.len:
#    result.add optNimTree(x[i], sym,

var reccount{.compiletime.} = 0
proc inlineLetsR(x: NimNode, sym,repl,stmts: var seq[NimNode]): NimNode =
  #echo "new tree"
  #echo x.treeRepr
  case x.kind
  of nnkCommentStmt:
    result = newEmptyNode()
  of nnkStmtList:
    var bstmts = newSeq[NimNode](0)
    for i in 0..<x.len:
      let reccount0 = reccount
      inc reccount
      #echo "label", reccount0, ": stmtlistin"
      #echo x[i].repr
      let r = inlineLetsR(x[i], sym, repl, bstmts)
      #echo "label", reccount0, ": stmtlistout"
      #echo r.repr
      if r.kind != nnkEmpty:
        bstmts.add r
    case bstmts.len
    of 0: result = newEmptyNode()
    of 1: result = bstmts[0]
    else:
      result = x.copyNimNode
      #result = newNimNode(nnkStmtList)
      for c in bstmts:
        #echo "bstmts: ", c.repr
        result.add c
  of nnkStmtListExpr:
    for i in 0..(x.len-2):
      let r = inlineLetsR(x[i], sym, repl, stmts)
      if r.kind != nnkEmpty:
        stmts.add r
    result = inlineLetsR(x[^1], sym, repl, stmts)
  of nnkLetSection:
    result = x.copyNimNode
    #result = newNimNode(nnkLetSection)
    for i in 0..<x.len:
      if x[i].kind==nnkIdentDefs:
        let r = inlineLetsR(x[i][2], sym, repl, stmts)
        if r.kind in CallNodes:
          var id = x[i].copyNimNode
          #var id = newNimNode(nnkIdentDefs)
          id.add x[i][0]
          id.add x[i][1]
          id.add r
          result.add id
        else:
          #echo "let: ", x[i][0].repr, " = ", r.repr
          #echo c[id][2].treerepr
          sym.add x[i][0]
          #echo "sym: ", sym[^1]
          repl.add r
      else:
        #result.add x[i]
        echo "error: nnkLetSection expected nnkIdentDefs"
        echo x.treerepr
        quit -1
    if result.len==0: result = newEmptyNode()
  of nnkSym:
    var i = sym.len-1
    while i>=0:
      if sym[i].repr == x.repr: break
      dec i
    if i>=0:
      #echo "found: ", i, " : ", sym[i].repr, " -> ", repl[i].repr
      result = repl[i]
    else:
      result = x
  of nnkOpenSymChoice:
    result = x
  #of nnkVarSection:
  #  result = x
  #of nnkIdentDefs:
  #  result = x.copyNimNode
  #  #result = newNimNode(nnkIdentDefs)
  #  result.add x[0]
  #  result.add x[1]
  #  result.add inlineLetsR(x[2], sym, repl, stmts)
  of nnkBlockStmt:
    #echo "nnkBlockStmt"
    var bstmts = newSeq[NimNode](0)
    let nsym = sym.len
    for i in 1..<x.len:
      let reccount0 = reccount
      inc reccount
      #echo "label", reccount0, ": blockin"
      #echo x[i].repr
      let r = inlineLetsR(x[i], sym, repl, bstmts)
      #echo "label", reccount0, ": blockout"
      #echo r.repr
      if r.kind != nnkEmpty:
        bstmts.add r
    sym.setLen(nsym)
    repl.setLen(nsym)
    if bstmts.len>0 or x[0].kind!=nnkEmpty:
      result = x.copyNimNode
      #result = newNimNode(nnkBlockStmt)
      result.add x[0]
      for c in bstmts:
        #echo "bstmts: ", c.repr
        result.add c
    else:
      result = newEmptyNode()
  of nnkObjConstr:
    result = x.copyNimNode
    #result = newNimNode(nnkObjConstr)
    result.add inlineLetsR(x[0], sym, repl, stmts)
    for i in 1..<x.len:
      var t = x[i].copyNimNode
      #var t = newNimNode(nnkExprColonExpr)
      t.add x[i][0]
      t.add inlineLetsR(x[i][1], sym, repl, stmts)
      result.add t
  of nnkDotExpr:
    let o = inlineLetsR(x[0], sym, repl, stmts)
    if o.kind == nnkObjConstr:
    #if false:
      let ss = (if x[1].kind==nnkSym: x[1] else: x[1][0]).repr
      var i = o.len - 1
      while i>0:
        if o[i][0].repr == ss: break
        dec i
      if i==0:
        result = x.copyNimNode
        #result = newNimNode(nnkDotExpr)
        result.add o
        result.add x[1]
      else:
        result = o[i][1]
        #echo "objConstr:"
        #echo x.repr
        #echo result.repr
        #echo x.getTypeImpl.repr
        #echo result.getTypeImpl.repr
    else:
      result = x.copyNimNode
      #result = newNimNode(nnkDotExpr)
      result.add o
      result.add x[1]
  of nnkPragma:
    if x.len==1 and x[0].kind==nnkExprColonExpr and $x[0][0]=="emit":
      template emt(x): untyped =
        {.emit: x.}
      result = getAst(emt(x[0][1][0]))
    else:
      result = x.copyNimNode
      for i in 0..<x.len:
        result.add inlineLetsR(x[i], sym, repl, stmts)
  of {nnkNone, nnkEmpty, nnkIdent, nnkType}:
    result = x
  of nnkLiterals:
    result = x
  #of nnkTypeOfExpr:
  #  result = x
  #of nnkConv:
  #  result = x.copyNimNode
  #  result.add x[0]
  #  result.add inlineLetsR(x[1], sym, repl, stmts)
  #of nnkCall:
  #  echo "call: ", $x[0]
  #  if $x[0]=="type":
  #    result = x
  #  else:
  #    result = x.copyNimNode
  #    #result = newNimNode(x.kind)
  #    for i in 0..<x.len:
  #      #echo "Xelse"
  #      #echo x[i].repr
  #      result.add inlineLetsR(x[i], sym, repl, stmts)
  else:
    #if x.kind==nnkPragma:
      #echo x.treerepr
    result = x.copyNimNode
    #result = newNimNode(x.kind)
    for i in 0..<x.len:
      #echo "Xelse"
      #echo x[i].repr
      result.add inlineLetsR(x[i], sym, repl, stmts)
    #if x.kind==nnkStrLit:
    #  result = newLit(x.strval)

proc inlineLets(x: NimNode): NimNode =
  var sym = newSeq[NimNode](0)
  var repl = newSeq[NimNode](0)
  var stmts = newSeq[NimNode](0)
  let r = inlineLetsR(x, sym, repl, stmts)
  result = newStmtList()
  for c in stmts:
    result.add c
  result.add r

macro optimizeAst*(a: typed): untyped =
  #echo "optimizeAst in"
  #echo a.treerepr
  #echo a.repr
  #let ar = a.repr
  #result = a
  #result = optimizeAstR(a)
  result = inlineLets(a)
  #echo "optimizeAst out"
  #echo result.treerepr
  #echo result.repr
  #let rr = result.repr
  #echo "ar == rr: ", ar==rr

macro XoptimizeAst*(a: typed): untyped = a
