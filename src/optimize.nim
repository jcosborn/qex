import macros
import metaUtils

proc newStmtListForCall(n:NimNode):NimNode =
  if n[0].getTypeInst[0][0].kind == nnkEmpty:
    result = newStmtList()
  else:
    result = newNimNode(nnkStmtListExpr)

proc moveStmtExprCall(body:NimNode):NimNode =
  depthFirst2(body):
    if it.kind==nnkCall:
      var rewrite = false
      for s in it:
        if s.kind==nnkStmtListExpr:
          rewrite = true
      if rewrite:
        var sl = newStmtListForCall(it)
        for ia in 1..<it.len:
          let s = it[ia]
          if s.kind==nnkStmtListExpr:
            for i in 0..(s.len-2):
              sl.add s[i]
            it[ia] = s[s.len-1]
        sl.add it
        it = sl
        #echo sl.treerepr
  #echo result.treeRepr
  #echo result.repr

proc isInline(n:NimNode):bool =
  result = false
  let p = n[4]
  #echo p.treerepr
  if p.kind==nnkPragma:
    for c in p:
      if c.repr=="inline":
        #echo "inlined: ", n[0]
        result = true
        break

template formalParams(x:NimNode):untyped = x[3]

proc isVar(x:NimNode):bool =
  #echo x.treeRepr
  result = false
  case x.kind:
    of nnkVarTy:
      result = true
    of nnkSym:
      result = isVar(x.getTypeImpl)
    of nnkIdentDefs:
      result = isVar(x[1])
    of nnkBracketExpr:
      if x[0].repr=="var":
        result = true
      else:
        discard
        #for c in x:
        #  result = isVar(c)
        #  if result: return
    of nnkObjectTy: discard
    else:
      echo "unhandled node kind in isVar: ", x.kind

proc replaceSym(b,s,r:NimNode):NimNode =
  var ss:string
  if s.kind==nnkIdent: ss = $s
  else: ss = $(s.symbol)
  depthFirst2(b):
    #if it.kind==nnkSym and it.symbol==s.symbol: it = r
    #if it.kind==nnkSym:
    #  echo it.symbol.repr
    #  echo s.symbol.repr
    #  echo eqIdent(it, ss)
    if eqIdent(it, ss):
      return copyNimTree(r)

proc fixAsgn(b:NimNode):NimNode =
  depthFirst2(b):
    #if it.kind==nnkFastAsgn:
    #  echo "fastAsgn: ", it.repr
    #  var n = newNimNode(nnkAsgn)
    #  n.add it[0]
    #  n.add it[1][1]
    #  #return n
    #  return newEmptyNode()
    if it.kind==nnkVarSection and
        it[0].kind==nnkIdentDefs and
        it[0][1].kind==nnkEmpty and
        it[0][2].kind==nnkEmpty:
      it[0][1] = newIdentNode("int")

proc replaceSimpleLet(b:NimNode):NimNode =
  proc walkTree(x:var NimNode) =
    for ic in 0..<x.len:
      let c = x[ic]
      if c.kind==nnkLetSection:
        echo c.treerepr
        var id = 0
        while id<c.len:
          if c[id].kind==nnkIdentDefs and c[id][2].kind==nnkSym:
            x = replaceSym(x, c[id][0], c[id][2])
            if c.len==1:
              x.del(ic)
            else:
              x[ic].del(id)
            return
          inc id
      else:
        while true:
          var t = x[ic]
          walkTree(t)
          if t==x[ic]: break
          x[ic] = t
  var bb = b
  while true:
    let t = bb
    walkTree(bb)
    if bb==t: break
  return bb

proc inlineLets(b:NimNode):NimNode =
  proc walkTree(x:var NimNode) =
    for ic in 0..<x.len:
      let c = x[ic]
      if c.kind==nnkLetSection:
        echo c.treerepr
        var id = 0
        while id<c.len:
          if c[id].kind==nnkIdentDefs and c[id][2].kind!=nnkObjConstr:
            #echo c[id][2].treerepr
            x = replaceSym(x, c[id][0], c[id][2])
            if c.len==1:
              x.del(ic)
            else:
              x[ic].del(id)
            return
          inc id
      else:
        while true:
          var t = x[ic]
          walkTree(t)
          if t==x[ic]: break
          x[ic] = t
  var bb = b
  while true:
    let t = bb
    walkTree(bb)
    if bb==t: break
  return bb

var nInlineLet{.compileTime.} = 0
proc inlineProcs(body:NimNode):NimNode =
  depthFirst2(body):
    if it.kind==nnkCall and it[0].kind==nnkSym:
      var p = it[0].symbol.getImpl
      #echo p.treerepr
      if p.kind==nnkProcDef and p.isInline:
        #echo "inlined"
        template fp:expr = p.formalParams
        #let fp = p.formalParams
        #echo p.formalParams[1].treerepr
        #echo p.formalParams[1][1].treerepr
        #echo p.formalParams[1][1].getTypeImpl.treerepr
        #echo isVar(p.formalParams[1])
        var sl = newStmtListForCall(it)
        for ia in 1..<it.len:
          if isVar(fp[ia]):
            #echo "var: ", fp[ia].repr
            p = p.replaceSym(fp[ia][0], it[ia])
            #echo p.repr
          else:
            #echo "let: ", fp[ia].repr
            let t = genSym(nskLet, "arg" & $nInlineLet)
            inc nInlineLet
            template letX(x,y:untyped):untyped =
              let x = y
            let l = getAst(letX(t, it[ia]))
            sl.add l[0]
            p = p.replaceSym(fp[ia][0], t)
        #echo p.treeRepr
        sl.add p[6]
        return sl
  #xecho result.repr
  #echo result.treerepr

macro optimize(body:typed):auto =
  result = body
  echo result.repr
  result = moveStmtExprCall(result)
  #echo result.treerepr
  result = inlineProcs(result)
  #echo result.treerepr
  result = moveStmtExprCall(result)
  result = inlineProcs(result)
  #echo result.repr
  result = moveStmtExprCall(result)
  result = inlineProcs(result)
  #echo result.repr
  #echo result.treerepr
  #result = replaceSimpleLet(result)
  result = inlineLets(result)
  #echo result.repr
  echo result.repr

when isMainModule:
  import qex
  import qcdTypes
  import gaugeUtils

  let defaultLat = @[8,8,8,8]
  defaultSetup()

  var v1 = lo.ColorVector()
  var v2 = lo.ColorVector()

  threads:
    v1 := 1
    threadBarrier()

    proc foo(x:int) = discard
    optimize:
      v2 := g[0]*v1
      #foo((let x=1;x))
