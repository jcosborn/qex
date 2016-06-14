import macros
import metaUtils
import strUtils

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

proc cancelDerefAddr(body:NimNode):NimNode =
  depthFirst2(body):
    if it.kind in {nnkDerefExpr,nnkHiddenDeref}:
      if it[0].kind in {nnkAddr,nnkHiddenAddr}:
        it = it[0][0]

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
      case $x:
        of "int": discard
        else:
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

proc eqIdents(x,y:NimNode):bool =
  var s:string
  if y.kind==nnkIdent: s = $y
  else: s = $(y.symbol)
  result = eqIdent(x, s)

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
        #echo c.treerepr
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
  #echo b.treeRepr
  proc walkTree(x:var NimNode) =
    #echo "new tree"
    #echo x.treeRepr
    for ic in 0..<x.len:
      let c = x[ic]
      if c.kind==nnkLetSection:
        #echo c.treerepr
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
          #echo "tree changed"
          x[ic] = t
  var bb = b
  while true:
    let t = bb
    walkTree(bb)
    if bb==t: break
  return bb

proc unrollType(vl,ty:NimNode):seq[tuple[val:NimNode,typ:NimNode]] =
  result.newSeq(0)
  case ty.kind:
    of nnkObjectTy:
      let flds = ty[2]
      for i in 0..<flds.len:
        #echo flds[i].treerepr
        #echo flds[i][0].treerepr
        let vln = newDotExpr(vl, flds[i][0])
        #echo vln.treeRepr
        let r = unrollType(vln, flds[i][0].getTypeImpl)
        #for s in r:
        #  echo s.val.repr
        #  echo s.typ.repr
        result.add r
    of nnkTupleTy:
      for i in 0..<ty.len:
        #echo ty[i].treerepr
        #echo ty[i][0].treerepr
        let vln = newDotExpr(vl, ty[i][0])
        #echo vln.treeRepr
        let ti = ty[i][1].getTypeImpl
        let r = unrollType(vln, ty[i][1].getTypeImpl)
        if r.len==0:  # terminal type
          result.add ((vln,ty[i][1]))
        else:
          #for s in r:
          #  echo s.val.repr
          #  echo s.typ.repr
          result.add r
    of nnkBracketExpr:
      #  echo ty.treeRepr
      case $ty[0]:
        of "array":
          if ty[1].kind==nnkInfix and $ty[1][0]==".." and
             ty[1][1].kind==nnkIntLit and ty[1][2].kind==nnkIntLit:
            let i0 = ty[1][1].intVal
            let i1 = ty[1][2].intVal
            for i in i0..i1:
              let vln = newNimNode(nnkBracketExpr).add(vl).add(newLit(i))
              #echo vln.treeRepr
              let r = unrollType(vln, ty[2].getTypeImpl)
              #for s in r:
              #  echo s.val.repr
              #  echo s.typ.repr
              result.add r
        else:
          #echo "generic unrollType nnkBracketExpr: ", $ty[0]
          #echo ty.treeRepr
          let r = unrollType(vl, ty[0].getTypeImpl)
          #for s in r:
          #  echo s.val.repr
          #  echo s.typ.repr
          result.add r
    else:
      echo "unhandled unrollType: ", ty.kind
      echo ty.treeRepr

proc sameTree(x,y:NimNode):bool =
  result = true
  #echo "kind: ", x.kind
  case x.kind:
    of nnkCharLit..nnkUInt64Lit:
      if (y.kind notin nnkCharLit..nnkUInt64Lit) or (x.intVal!=y.intVal):
        result = false
    of nnkFloatLit..nnkFloat64Lit:
      if (y.kind notin nnkFloatLit..nnkFloat64Lit) or (x.floatVal!=y.floatVal):
        result = false
    of nnkStrLit..nnkTripleStrLit:
      if (y.kind notin nnkStrLit..nnkTripleStrLit) or (x.strVal==y.strVal):
        result = false
    of nnkIdent, nnkSym:
      if (y.kind!=nnkIdent and y.kind!=nnkSym) or not eqIdents(x, y):
        result = false
    else:
      if x.kind!=y.kind:
        result = false
      else:
        for i in 0..<x.len:
           if not sameTree(x[i], y[i]):
             result = false
             break
  #if result==true:
  #  echo "same:"
  #  echo " ", x.repr
  #  echo " ", y.repr

proc replaceTree(b,s,r:NimNode):NimNode =
  result = b
  if sameTree(s, b):
    result = r
  else:
    for i in 0..<result.len:
      result[i] = replaceTree(result[i], s, r)

proc unrollVars(b:NimNode):NimNode =
  proc walkTree(x:var NimNode) =
    var ic = 0
    while ic<x.len:
      let c = x[ic]
      if c.kind==nnkVarSection:
        var cc = c
        #echo c.treerepr
        var id = 0
        while id<c.len:
          let cid = c[id]
          if cid.kind==nnkIdentDefs and cid[2].kind==nnkEmpty:
            let t = cid[0].getTypeImpl
            #echo c[id][0].treerepr
            #echo c[id][0].getTypeImpl.treerepr
            #echo c[id][0].getTypeImpl.repr
            let r = unrollType(cid[0], t)
            #echo r.len
            if r.len>0:
              template vardef(v,t:untyped) =
                var v{.noInit.}:t
              var vs = newSeq[NimNode](r.len)
              for i in 0..<r.len:
                let vn = r[i].val.repr.replace("[","").replace("]","")
                  .replace(".","")
                let vr = genSym(nskVar, vn)
                let ga = getAst(vardef(vr, r[i].typ))
                vs[i] = vr
                echo vs[i].repr
                echo ga[0][0].treerepr
                if i==0:
                  c[id] = ga[0][0]
                else:
                  c.add ga[0][0]
                  #for j in id..<c.len:c.insert(ga[0][0], id)
                inc id
              for ie in (ic+1)..(<x.len):
                for i in 0..<r.len:
                  x[ie] = replaceTree(x[ie], r[i].val, vs[i])
          inc id
      else:
        while true:
          var t = x[ic]
          walkTree(t)
          if t==x[ic]: break
          x[ic] = t
      inc ic
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
        #if sl.kind==nnkStmtList:
        #  sl = newBlockStmt(sl)
        return sl
  #xecho result.repr
  #echo result.treerepr

macro checkOpt(body:typed):auto = body

macro optimize*(body:typed):auto =
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
  result = cancelDerefAddr(result)
  #echo result.repr
  #echo result.treerepr
  #result = replaceSimpleLet(result)
  result = inlineLets(result)
  echo result.treerepr
  result = unrollVars(result)
  #result = symToIdent(result)
  #result = newCall("checkOpt", result)
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
