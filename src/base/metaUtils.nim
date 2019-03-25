import macros
import strUtils

proc isMagic(x: NimNode): bool =
  # echo x.treerepr
  let pragmas = x[4]
  if pragmas.kind==nnkPragma and pragmas[0].kind==nnkExprColonExpr and
     $pragmas[0][0]=="magic": result = true

template isNotMagic(x: NimNode): bool = not isMagic(x)

proc getParam(fp: NimNode, n: int): auto =
  # n counts from 1
  var n = n-1
  for i in 1..<fp.len:
    let c = fp[i]
    if n >= c.len-2: n -= c.len-2
    else: return (c[n].copyNimTree, c[^2].copyNimTree)

proc has(n:NimNode, k:NimNodeKind):bool =
  for c in n:
    if c.kind == k: return true
  return false

proc replace(n,x,y:NimNode):NimNode =
  if n == x:
    result = y.copyNimTree
  else:
    result = n.copyNimNode
    for c in n:
      result.add c.replace(x,y)

proc replaceId(n,x,y:NimNode):NimNode =
  # Same as replace but only replace eqIdent identifier.
  if n.kind == nnkIdent and n.eqIdent(x):
    result = y.copyNimTree
  else:
    result = n.copyNimNode
    for c in n:
      result.add c.replaceId(x,y)

proc replaceAlt(n,x,y:NimNode, k:NimNodeKind):NimNode =
  # Same as replace but the optional parent node kind k is included in the replacement.
  if n.kind == k and n.len==1 and n[0] == x:
    result = y.copyNimTree
  elif n == x:
    result = y.copyNimTree
  else:
    result = n.copyNimNode
    for c in n:
      result.add c.replaceAlt(x,y,k)

proc replaceExcl(n,x,y:NimNode, k:NimNodeKind):NimNode =
  # Same as replace but the optional parent node kind k excludes the replacement.
  if n.kind == k and n.len==1 and n[0] == x:
    result = n.copyNimTree
  elif n == x:
    result = y.copyNimTree
  else:
    result = n.copyNimNode
    for c in n:
      result.add c.replaceExcl(x,y,k)

proc replaceNonDeclSym(b,s,r: NimNode, extra:NimNodeKind = nnkEmpty): NimNode =
  # Replace a symbol `s` that's not declared in the body `b` with `r`.
  # Assuming a unique symbol exists.  Only works with trees of symbols.
  var ss = s.strVal
  # echo "replacing ",ss
  var
    declSyms = newPar()
    theSym = newEmptyNode()
  proc checkSym(n:NimNode) =
    if theSym == n: return
    var f = false
    for c in declSyms:
      if c == n:
        f = true
        break
    if f: return
    elif theSym.kind == nnkEmpty: theSym = n
    else:
      echo "Internal ERROR: replaceNonDeclSym: multiple ",s.repr," found in:"
      echo b.treerepr
      echo "found: ",theSym.lineinfo," :: ",theSym.lisprepr
      echo "found: ",n.lineinfo," :: ",n.lisprepr
      quit 1
  proc find(n:NimNode) =
    # echo "declSyms: ",declSyms.repr
    # echo "theSym: ",theSym.repr
    case n.kind:
    of nnkSym:
      if n.eqIdent ss: checkSym n
    of nnkIdentDefs, nnkConstDef:
      for i in 0..<n.len-2:
        if n[i].eqIdent ss: declSyms.add n[i]
      find n[^1]
    of nnkExprColonExpr:
      find n[1]
    of nnkDotExpr:
      find n[0]
    of nnkProcDef, nnkMethodDef, nnkDo, nnkLambda, nnkIteratorDef,
       nnkTemplateDef, nnkConverterDef:
      echo "Internal ERROR: replaceNonDeclSym: unhandled cases."
      quit 1
    else:
      for c in n: find c
  find b
  # echo "declSyms: ",declSyms.repr
  # echo "theSym: ",theSym.repr
  if theSym.kind != nnkEmpty:
    result = b.replaceAlt(theSym, r, extra)
  else:
    result = b

proc rebuild*(n:NimNode):NimNode =
  # Typed AST has extra information in its nodes.
  # Replacing nodes in a typed AST can break its consistencies,
  # for which the compiler is not well prepared.
  # Here we simply rebuild some offensive nodes from scratch,
  # and force the compiler to rebuild its type information.

  # Note that the compiler currently (v0.17) only retype the AST
  # after a macro returns, so to preserve type information while
  # traversing the AST, call this proc on the result
  #    result = rebuild result
  # just before macro returns.

  # Special node kinds have to be taken care of.
  if n.kind == nnkAddr:
    # We process typed NimNode, so addr is already checked.
    result = newCall(bindsym"unsafeAddr", rebuild n[0])
  elif n.kind == nnkConv:
    result = newNimNode(nnkCall, n).add(rebuild n[0], rebuild n[1])
  elif n.kind in nnkCallKinds and n[^1].kind == nnkBracket and
       n[^1].len>0 and n[^1].has(nnkHiddenCallConv):
    # special case of varargs
    result = newCall(rebuild n[0])
    for i in 1..<n.len-1: result.add rebuild n[i]
    for c in n[^1]:
      if c.kind == nnkHiddenCallConv: result.add rebuild newcall(c[0], c[1])
      else: result.add rebuild c
#  elif n.kind in nnkCallKinds and n[0] == bindsym"echo" and n.len>0 and n[1].kind == nnkBracket:
#    # One dirty hack for the builtin echo, with no nnkHiddenCallConv (this is been caught above)
#    result = newCall(rebuild n[0])
#    for c in n[1]: result.add rebuild c
#    echo "In rebuild: echo: ",result.treerepr
#  elif n.kind in nnkCallKinds and n[^1].kind == nnkHiddenStdConv and n[^1][1].kind == nnkBracket and
#       n[^1][1].len>0 and n[^1][1].has(nnkHiddenCallConv):
  elif n.kind in nnkCallKinds and n[^1].kind == nnkHiddenStdConv and n[^1][1].kind == nnkBracket and
       n[^1][1].len>0:
    # Deals with varargs
    result = newCall(rebuild n[0])
    for i in 1..<n.len-1: result.add rebuild n[i]
    for c in n[^1][1]:
      if c.kind == nnkHiddenCallConv: result.add rebuild newcall(c[0], c[1])
      else: result.add rebuild c
  elif n.kind in nnkCallKinds:
    result = newCall(rebuild n[0])
    for i in 1..<n.len: result.add rebuild n[i]
  elif n.kind == nnkHiddenStdConv:
    # Generic HiddenStdConv
    result = rebuild n[1]
  elif n.kind == nnkHiddenAddr and n[0].kind == nnkHiddenDeref:
    result = rebuild n[0][0]
  elif n.kind == nnkHiddenDeref and n[0].kind == nnkHiddenAddr:
    result = rebuild n[0][0]
  elif n.kind in {nnkHiddenAddr,nnkHiddenDeref}:
    result = rebuild n[0]
  elif n.kind == nnkTypeSection:
    # Type section is special.  Once the type is instantiated, it exists, and we don't want duplicates.
    result = newNimNode(nnkDiscardStmt,n).add(newStrLitNode(n.repr))

  # Strip information from other kinds
  else:
    if n.kind in AtomicNodes:
      result = n.copyNimNode
    else:
      result = newNimNode(n.kind, n)
#[
    # If something breaks, try adding the offensive node here.
    #if n.kind in nnkCallKinds + {nnkBracketExpr,nnkBracket,nnkDotExpr}:
    if n.kind in nnkCallKinds + exprNodes + {nnkAsgn}:
      result = newNimNode(n.kind, n)
    # Copy other kinds of node.
    else:
      result = n.copyNimNode
]#
    for c in n:
      result.add rebuild c
  # echo result.treerepr
  # echo "### leave rebuild"

proc append(x,y:NimNode) =
  for c in y: x.add c

proc inlineLets(n:NimNode):NimNode =
  proc get(n:NimNode):NimNode =
    result = newPar()
    if n.kind == nnkLetSection:
      for d in n:
        if d.kind != nnkIdentDefs or d.len<3:
          echo "Internal ERROR: regenSym: get: can't handle:"
          echo n.treerepr
          quit 1
        if d[^1].kind in AtomicNodes:
          for i in 0..<d.len-2:   # Last 2 is type and value.
            if d[i].kind == nnkSym:
              result.add newPar(d[i],d[^1])
            else: error("inlineLets can't handel: " & n.treerepr)
          for c in d[^1]: result.append get c
    else:
      for c in n: result.append get c
  proc rep(n,x,y:NimNode):NimNode =
    if n == x: result = y.copy
    elif n.kind == nnkLetSection:
      var ll = n.copyNimNode
      for d in n:
        var dd = d.copyNimNode
        for i in 0..<d.len-2:
          if d[i] != x: dd.add d[i].copy
        if dd.len > 0:
          dd.add(d[^2].copy, d[^1].rep(x,y))
          ll.add dd
      if ll.len > 0: result = ll
      else: result = newNimNode(nnkDiscardStmt,n).add(newStrLitNode(n.repr))
    else:
      result = n.copyNimNode
      for c in n:
        result.add c.rep(x,y)
  result = n.copy
  for x in get n: result = result.rep(x[0],x[1])

proc regenSym(n:NimNode):NimNode =
  # Only regen nskVar and nskLet symbols.

  # We need to regenerate symbols for multiple inlined procs,
  # because cpp backend put variables on top level, although
  # the c backend works without this.
  proc get(n:NimNode,k:NimNodeKind):NimNode =
    result = newPar()
    if n.kind == k:
      for d in n:
        if d.kind != nnkIdentDefs or d.len<3:
          echo "Internal ERROR: regenSym: get: can't handle:"
          echo n.treerepr
          quit 1
        for i in 0..<d.len-2:   # Last 2 is type and value.
          if d[i].kind == nnkSym: result.add d[i]
        for c in d[^1]: result.append c.get k
    else:
      for c in n: result.append c.get k
  result = n.copyNimTree
  # We ignore anything inside a typeOfExpr, because we need the
  # type information in there, but our new symbols wouldn't have
  # any type info.
  for x in result.get nnkLetSection:
    #echo "Regen Let: ",x.repr
    let y = genSym(nskLet, x.strVal)
    result = result.replaceExcl(x,y,nnkTypeOfExpr)
  for x in result.get nnkVarSection:
    #echo "Regen Var: ",x.repr
    let y = genSym(nskVar, x.strVal)
    result = result.replaceExcl(x,y,nnkTypeOfExpr)
macro regenSym*(n: typed): untyped = regenSym(n)

proc matchGeneric(n,ty,g:NimNode):NimNode =
  ## Match generic type `ty`, with the instantiation node `n`, and return
  ## the instantiation type of the generic identifier `g`.
  # echo "MG:I: ",n.lisprepr
  # echo "MG:T: ",ty.lisprepr
  # echo "MG:G: ",g.lisprepr
  proc isG(n:NimNode):bool =
    n.kind == nnkIdent and n.eqIdent($g)
  proc typeof(n:NimNode):NimNode =
    newNimNode(nnkTypeOfExpr,g).add n
  proc getGParams(ti:NimNode):NimNode =
    # ti ~ type[G0,G1,...], is from gettypeinst
    # We go through the implementation to find the correct generic names.
    ti.expectKind nnkBracketExpr
    let tn = ti[0]
    tn.expectKind nnkSym
    let td = tn.getImpl
    td.expectKind nnkTypeDef
    result = td[1]
    result.expectKind nnkGenericParams
  proc matchT(ti,ty,g:NimNode):NimNode =
    # match instantiation type `ti`, with generic type `ty`
    # recursively find the chain of generic type variables
    # correponding to `g` in `ty`.
    result = newPar()
    var i = 0
    let tg = getGParams ti
    while i<ty.len:
      if ty[i].isG: break
      inc i
    if i == 0: return
    elif i < ty.len:
      if tg.len < i: return
      else: result.add tg[i-1]
    else:
      for i in 1..<ty.len:
        if ty[i].kind == nnkBracketExpr:
          if i < ti.len:
            ti[i].expectKind nnkBracketExpr
            result = matchT(ti[i],ty[i],g)
            if result.len > 0:
              result.add tg[i-1]
              return
  if ty.isG: return typeof n
  elif ty.kind == nnkBracketExpr:
    let ts = matchT(n.gettypeinst,ty,g)
    result = n
    if ts.len > 0:
      for i in countdown(ts.len-1,0): result = result.newDotExpr ts[i]
      return
  echo "Internal WARNING: matchGeneric: Unsupported"
  echo "MG:I: ",n.lisprepr
  echo "MG:T: ",ty.lisprepr
  echo "MG:G: ",g.lisprepr

proc cleanIterator(n:NimNode):NimNode =
  var fa = newPar()
  proc replaceFastAsgn(n:NimNode):NimNode =
    if n.kind == nnkFastAsgn:
      let n0 = genSym(nskLet, n[0].strVal)
      fa.add newPar(n[0],n[1],n0)
      template asgn(x,y:untyped):untyped =
        let x = y
      let n1 = replaceFastAsgn n[1]
      result = getAst(asgn(n0,n1))
    else:
      result = n.copyNimNode
      for c in n: result.add replaceFastAsgn c
  proc removeDeclare(n:NimNode):NimNode =
    if n.kind == nnkVarSection:
      var keep = newseq[int](0)
      for c in 0..<n.len:
        var i = 0
        while i < fa.len:
          var j = 0
          while j < n[c].len-2:
            if fa[i][0] == n[c][j]: break
            inc j
          if j < n[c].len-2:
            if n[c].len > 3:
              echo "Internal ERROR: cleanIterator: removeDeclare: unhandled situation"
              echo n.treerepr
              quit 1
            break
          inc i
        if i < fa.len and (n[c][^2].kind != nnkEmpty or n[c][^1].kind != nnkEmpty):
          echo "Internal ERROR: cleanIterator: removeDeclare: unhandled situation"
          echo n.treerepr
          quit 1
        elif i >= fa.len: keep.add c
      # echo keep," ",n.repr
      if keep.len == 0:
        # echo "Removing declaration: ",n.lisprepr
        result = newNimNode(nnkDiscardStmt,n).add newStrLitNode(n.repr)
      else:
        result = n.copyNimNode
        for i in keep: result.add removeDeclare n[i]
    else:
      result = n.copyNimNode
      for c in n: result.add removeDeclare c
  result = replaceFastAsgn n
  if fa.len > 0:
    result = result.removeDeclare
    for x in fa:
      result = result.replace(x[0],x[2])
      # echo x[0].lisprepr,"\n  :: ",x[0].gettypeinst.lisprepr
      # echo x[1].lisprepr,"\n  :: ",x[1].gettypeinst.lisprepr
  proc fixDeclare(n:NimNode):NimNode =
    # Inlined iterators have var sections that are not clearly typed.
    # We try to find inconsistencies from the type of the actual symbol being declared.
    result = n.copyNimNode
    if n.kind == nnkVarSection:
      for i in 0..<n.len:
        result.add n[i].copyNimTree
        if n[i][^2].kind == nnkEmpty and n[i][^1].kind != nnkEmpty:
          for j in 0..<n[i].len-2:
            # echo n.treerepr
            # echo "sym ",i," ",j," : ",n[i][j].repr
            # echo "    :- ",n[i][^1].repr
            let
              t = n[i][j].gettypeinst
              r = n[i][^1].gettypeinst
            # echo "    ty: ",t.lisprepr
            # echo "    <-: ",r.lisprepr
            # echo "    ??: ",t==r
            if result[i][^2].kind != nnkEmpty and result[i][^2] != r:
              echo "Internal ERROR: cleanIterator: fixDeclare: unhandled situation"
              echo n.treerepr
              quit 1
            # echo "Fixing declaration: ",n[i].lisprepr
            if t != r: result[i][^2] = newNimNode(nnkTypeOfExpr,n[i][j]).add n[i][j]
        result[i][^1] = fixDeclare result[i][^1]
      # echo result.repr
    else:
      for c in n: result.add fixDeclare c
  result = fixDeclare result
  # echo "<<<<<< cleanIterator"
  # echo result.treerepr

proc inlineProcsY(call: NimNode, procImpl: NimNode): NimNode =
  # echo ">>>>>> inlineProcsY"
  # echo "call:\n", call.lisprepr
  # echo "procImpl:\n", procImpl.treerepr
  let fp = procImpl[3]  # formal params
  proc removeRoutines(n:NimNode):NimNode =
    # We are inlining, so we don't need RoutineNodes anymore.
    if n.kind in RoutineNodes:
      result = newNimNode(nnkDiscardStmt,n).add(newStrLitNode(n.repr))
    else:
      result = n.copyNimNode
      for c in n: result.add removeRoutines c
  proc removeTypeSections(n:NimNode):NimNode =
    # Type section is special.  Once the type is instantiated, it exists, and we don't want duplicates.
    if n.kind == nnkTypeSection:
      result = newNimNode(nnkDiscardStmt,n).add(newStrLitNode(n.repr))
    else:
      result = n.copyNimNode
      for c in n: result.add removeTypeSections c
  var
    pre = newStmtList()
    body = procImpl.body.copyNimTree.removeRoutines.removeTypeSections
  # echo "### body w/o routines:"
  # echo body.repr
  body = cleanIterator body
  # echo "### body after clean up iterator:"
  # echo body.repr
  for i in 1..<call.len:  # loop over call arguments
    # We need to take care of the case when one argument use the same identifier
    # as one formal parameter.  Reusing the formal parameter identifiers is OK.
    let
      (sym,typ) = getParam(fp, i)
      t = genSym(nskLet, $sym)
      c = genSym(nskConst, $sym)
    template letX(x,y: untyped): untyped =
      let x = y
    template constX(x,y: untyped): untyped =
      const x = y
    # let p = if call[i].kind in {nnkHiddenAddr,nnkHiddenDeref}: call[i][0] else: call[i]
    let p = call[i]
    # echo "parameter: ",i," : ",p.treerepr
    # echo "sym: ",sym.lineinfo," :: ",sym.lisprepr
    # echo "typ: ",typ.lineinfo," :: ",typ.lisprepr
    if typ.kind == nnkStaticTy or (typ.kind == nnkBracketExpr and typ[0] == bindsym"static"):
      # echo typ.lisprepr
      # echo p.lisprepr
      if p.kind notin nnkLiterals:
        echo "ERROR: inlineProcsY: param type: ",typ.lisprepr
        echo "    received a non-literal node: ",p.lisprepr
        quit 1
      # We do nothing, assuming the compiler has finished constant unfolding.
    elif p.kind in nnkLiterals:
      pre.add getAst(constX(c, p))
      body = body.replaceNonDeclSym(sym, c)
    elif p.kind == nnkHiddenStdConv and p[^1].kind in nnkLiterals:
      pre.add getAst(constX(c, p[^1]))
      body = body.replaceNonDeclSym(sym, c)
    elif typ.kind == nnkVarTy:
      pre.add getAst(letX(t, newNimNode(nnkAddr,p).add p))
      body = body.replaceNonDeclSym(sym, newNimNode(nnkDerefExpr,p).add(t), nnkHiddenDeref)
    else:
      pre.add getAst(letX(t, p))
      body = body.replaceNonDeclSym(sym, t)
  # echo "### body with fp replaced:"
  # echo body.repr
  proc resolveGeneric(n:NimNode):NimNode =
    proc find(n:NimNode, s:string):bool =
      if n.kind == nnkDotExpr:
        # ignore n[1]
        return n[0].find s
      elif n.kind in RoutineNodes:
        return false
      elif n.kind == nnkIdent and n.eqIdent s:
        return true
      else:
        for c in n:
          if c.find s: return true
      return false
    var gs = newPar()
    if procImpl[5].kind == nnkBracket and procImpl[5].len>=2 and procImpl[5][1].kind == nnkGenericParams:
      let gp = procImpl[5][1]
      for c in gp:
        c.expectKind nnkIdentDefs
        for i in 0..<c.len-2:
          c[i].expectKind nnkIdent
          if n.find($c[i]): gs.add c[i]
    result = n
    # echo gs.lisprepr
    for g in gs:
      var j = 1
      var sym,typ:NimNode
      while j < call.len:
        (sym,typ) = fp.getParam j
        if typ.find($g): break
        inc j
      if j < call.len:
        # echo sym.treerepr
        # echo typ.treerepr
        # let tyi = call[j].gettypeinst
        # echo "timpl: ",call[j].gettypeimpl.lisprepr
        # echo "tinst: ",tyi.lisprepr
        # echo "impl: ",tyi[0].symbol.getimpl.lisprepr
        let inst = matchGeneric(call[j], typ, g)
        # echo inst.treerepr
        result = result.replaceId(g, inst)
      else:
        echo "Internal WARNING: resolveGeneric: couldn't find ",g.lisprepr
  body = resolveGeneric body
  # echo "### body after resolve generics:"
  # echo body.repr
  let blockname = genSym(nskLabel, $call[0])
  proc breakReturn(n:NimNode):NimNode =
    if n.kind == nnkReturnStmt:
      result = newStmtList()
      for c in n: result.add c
      result.add newNimNode(nnkBreakStmt, n).add blockname
      if result.len == 1: result = result[0]
    else:
      result = n.copyNimNode
      for c in n: result.add breakReturn c
  body = breakReturn body
  # echo "### body after replace return with break:"
  # echo body.repr
  var sl:NimNode
  if procImpl.len == 7:
    pre.add body
    sl = newBlockStmt(blockname, pre)
  elif procImpl.len == 8:
    # echo "TYPEof call: ",call.lisprepr," ",call.gettypeinst.treerepr
    # echo "TYPEof call: ",call.lisprepr," ",call.gettypeimpl.treerepr
    # echo "TYPEof call[0]: ",call[0].lisprepr," ",call[0].gettypeinst.treerepr
    # echo "TYPEof call[0]: ",call[0].lisprepr," ",call[0].gettypeimpl.treerepr
    # echo "TYPEof fp[0]: ",fp[0].lisprepr," ",fp[0].gettype.treerepr
    # echo "TYPEof pi[7]: ",procImpl[7].lisprepr," ",procImpl[7].gettype.treerepr
    template varX(x,t:untyped):untyped =
      var x: t
    template varXNI(x,t:untyped):untyped =
      var x {.noinit.} : t
    var calln:NimNode
    if call.kind == nnkHiddenCallConv:
      calln = newNimNode nnkCall
      call.copyChildrenTo calln
    else:
      calln = call.copyNimTree
    let
      #ty = call.gettypeinst
      #ty = call[0].gettypeinst[0][0]
      #ty = call[0].gettypeimpl[0][0]
      ty = newNimNode(nnkTypeOfExpr,call).add(calln)
      r = procImpl[7]
      z = genSym(nskVar, r.strVal)
      p = procImpl[4]
    #echo "    ty: ",ty.lisprepr
    var noinit = false
    if p.kind != nnkEmpty:
      # echo "pragmas: ", p.lisprepr
      p.expectKind nnkPragma
      for c in p:
        if c.eqIdent "noinit":
          noinit = true
          break
    let d = if noinit: getAst(varXNI(z,ty)) else: getAst(varX(z,ty))
    # if noinit: echo "noinit: ", d.lisprepr
    pre.add body.replace(r,z)
    sl = newBlockStmt(newNimNode(nnkStmtListExpr,call).add(d, newBlockStmt(blockname, pre), z))
  else:
    echo "Internal ERROR: inlineProcsY: unforeseen length of the proc implementation: ", procImpl.len
    quit 1
  # echo "====== sl"
  # echo sl.repr
  # echo "^^^^^^"
  # result = sl
  result = regenSym inlineLets sl
  # echo "<<<<<< inlineProcsY"
  # echo result.treerepr

proc callName(x: NimNode): NimNode =
  if x.kind in CallNodes: result = x[0]
  else: quit "callName: unknown kind (" & treeRepr(x) & ")\n" & repr(x)

proc inlineProcsX(body: NimNode): NimNode =
  # echo ">>>>>> inlineProcsX"
  # echo body.repr
  proc recurse(it: NimNode): NimNode =
    if it.kind == nnkTypeOfExpr: return it.copyNimTree
    if it.kind in CallNodes and it.callName.kind==nnkSym:
      let procImpl = it.callName.getImpl
      # echo "inspecting call"
      # echo it.lisprepr
      # echo procImpl.repr
      # echo isMagic(procImpl)
      if procImpl.kind == nnkTypeDef: return it.copyNimTree
      if procImpl.isNotMagic and
          procImpl.body.kind!=nnkEmpty and
          procImpl.kind != nnkIteratorDef:
        return recurse inlineProcsY(it, procImpl)
    result = copyNimNode(it)
    for c in it: result.add recurse c
  result = recurse(body)
  # echo "<<<<<< inlineProcsX"
  # echo result.repr

macro inlineProcs*(body: typed): auto =
  # echo ">>>>>> inlineProcs:"
  # echo body.repr
  # echo body.treerepr
  #result = body
  result = rebuild inlineProcsX body
  # echo "<<<<<< inlineProcs:"
  # echo result.repr
  # echo result.treerepr

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

macro `$`*(t: type): untyped =
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

macro echoRep*(x: typed): typed =
  echo x.lineinfo
  echo x.repr
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

macro toId*(s: static[string]):untyped =
  echo s
  ident(s)

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

proc peelStmt*(n:NimNode): NimNode =
  if n.len == 1 and n.kind == nnkStmtList: n[0]
  else: n

macro makeCall*(op:static[string],a:tuple):untyped =
  echo op
  echo a.repr
  #echo a[0].repr
  echo a.treeRepr
  result = newCall(ident(op))
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
  var opid = ident(op)
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

proc replaceComments*(n: NimNode; subs: varargs[(string,string)]): NimNode =
  #echo n.treeRepr
  if n.kind == nnkCommentStmt:
    echo n.strval
    let s = multiReplace(n.strVal, subs)
    echo s
    #echo subs
    result = newCommentStmtNode(s)
  #elif n.kind == nnkTypeDef:
  #  echo n.treerepr
  #  echo n[2].getType.treerepr
  else:
    result = copyNimNode(n)
    for c in n.children:
      result.add(replaceComments(c, subs))
  #echo result.repr

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
    result = result.replace(x[i], v)
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
    result = result.replace(x[i], t)
  #echo result.treerepr
  echo "subst: ", lineInfo(x[0])
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
    result = result.replace(x[i], t)
  result = newStmtList(prestmts, result)
  #echo result.repr
  echo "lets: ", lineInfo(x[0])
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
    #result.add(body.replace(index, newIntLitNode(i)))
    result.add(newBlockStmt(body.replace(index, newIntLitNode(i))))
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
    #result.add(body.replace(index, newIntLitNode(i)))
    result.add(newBlockStmt(body.replace(index, newIntLitNode(i))))
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
  #echo "### enter unrollFor"
  #echo n.repr
  #echo n.treerepr
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
  if n[1][1].kind == nnkHiddenStdConv:
    n[1][1] = n[1][1][1]
  must: n[1][1].kind in nnkCharLit..nnkUInt64Lit
  must: n[1][2].kind in nnkCharLit..nnkUInt64Lit
  let
    a = n[1][1].intval
    b = n[1][2].intval
  result = newStmtList()
  for i in a..b:
    result.add newNimNode(nnkBlockStmt, n).add(
        ident("ITR: " & $i & " :: \n#[" & n.repr & "]#\n"), n[2].replace(n[0], newIntLitNode(i)))
  #echo result.treerepr
  #echo result.repr
  #echo "### leave unrollFor"
macro unrollFor*(n:typed):untyped =
  if n.kind == nnkForStmt:
    result = n.unrollFor.rebuild
  else:
    result = newstmtlist()
    for c in n:
      if c.kind == nnkForStmt: result.add c.unrollFor.rebuild
      else: result.add c
  #echo n[^1].lineinfo
  #echo n.treerepr
  #echo result.treerepr

template forStaticUnRollFor*(index,i0,i1,body:untyped):untyped =
  const
    a = i0
    b = i1
  unrollFor:
    for index in a..b: body

template forStatic*(index,i0,i1,body:untyped):untyped =
  # forStaticUntyped(index,i0,i1,body)
  forStaticUnRollFor(index,i0,i1,body)

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

# transformations:
#  let x = y:StmtListExpr -> R( y[0..^2]; let x = y[^1] )
#  let x = y:ObjConstr -> R( let t_i = y[i] i=1..^1; let x=ObjConstr(t) )
# BlockStmt
# BracketExpr
# DotExpr
proc inlineLetsR(x: NimNode; sym,repl,stmts: var seq[NimNode]): NimNode

#[
proc optimizeObjConstr(x: NimNode; sym,repl: var seq[NimNode]): NimNode =
  x.expectKind(nnkObjConstr)
  result = x.copyNimNode
  result.add inlineLetsR(x[0], sym, repl)
  var sle = newNimNode(nnkStmtListExpr)
  for i in 1..<x.len:
    var t = x[i].copyNimNode
    #var t = newNimNode(nnkExprColonExpr)
    t.add x[i][0]
    var r = inlineLetsR(x[i][1], sym, repl, stmts)
    if r.kind == nnkStmtListExpr:
      for j in 0..(r.len-2):
        sle.add r[j]
      r = r[^1]
    t.add r
    result.add t
  if sle.len>0:
    sle.add result
    result = sle
]#

var reccount{.compiletime.} = 0
proc inlineLetsR(x: NimNode, sym,repl,stmts: var seq[NimNode]): NimNode =
  #echo "new tree"
  #echo x.treeRepr
  case x.kind
  of nnkCommentStmt:
    result = newEmptyNode()
  of nnkLiterals+{nnkNone, nnkEmpty, nnkIdent, nnkOpenSymChoice}:
    result = x
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

  of nnkLetSection:
    result = newStmtList()
    for i in 0..<x.len:
      if x[i].kind==nnkIdentDefs:
        var r = inlineLetsR(x[i][2], sym, repl, stmts)
        if r.kind == nnkStmtListExpr:
          for j in 0..(r.len-2):
            result.add r[j]
          r = r[^1]
        case r.kind
        of nnkSym:
          #echo "let: ", x[i][0].repr, " = ", r.repr
          #echo c[id][2].treerepr
          sym.add x[i][0]
          #echo "sym: ", sym[^1]
          repl.add r
        else:
          result.add newLetStmt(x[i][0], r)
      else:
        echo "error: nnkLetSection expected nnkIdentDefs"
        echo x.treerepr
        error "error"
    if result.len==0: result = newEmptyNode()

  of nnkStmtListExpr:
    result = x.copyNimNode
    for i in 0..(x.len-2):
      let r = inlineLetsR(x[i], sym, repl, stmts)
      if r.kind != nnkEmpty:
        result.add r
    let r = inlineLetsR(x[^1], sym, repl, stmts)
    if r.kind == nnkStmtListExpr:
      for i in 0..<r.len:
        result.add r[i]
    else:
      result.add r

  of nnkObjConstr:
    result = x.copyNimNode
    #result = newNimNode(nnkObjConstr)
    result.add inlineLetsR(x[0], sym, repl, stmts)
    var sle = newNimNode(nnkStmtListExpr)
    for i in 1..<x.len:
      var t = x[i].copyNimNode
      #var t = newNimNode(nnkExprColonExpr)
      t.add x[i][0]
      var r = inlineLetsR(x[i][1], sym, repl, stmts)
      if r.kind == nnkStmtListExpr:
        for j in 0..(r.len-2):
          sle.add r[j]
        r = r[^1]
      t.add r
      result.add t
    if sle.len>0:
      sle.add result
      result = sle

  of nnkBlockStmt:
    #echo "nnkBlockStmt"
    let nsym = sym.len
    result = x.copyNimNode
    result.add x[0]
    for i in 1..<x.len:
      result.add inlineLetsR(x[i], sym, repl, stmts)
    sym.setLen(nsym)
    repl.setLen(nsym)
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
        #echo "dotExpr objConstr:"
        #echo " ", x.repr
        #echo " ", o.repr
        #echo " ", x[1].repr
        result = x.copyNimNode
        #result = newNimNode(nnkDotExpr)
        result.add o
        result.add x[1]
      else:
        result = o[i][1]
        #echo "dotExpr objConstr:"
        #echo " ", x.repr
        #echo " ", result.repr
        #echo x.getTypeImpl.repr
        #echo result.getTypeImpl.repr
    else:
      result = x.copyNimNode
      #result = newNimNode(nnkDotExpr)
      result.add o
      result.add x[1]
  of nnkBracketExpr:
    let a = inlineLetsR(x[0], sym, repl, stmts)
    let k = inlineLetsR(x[1], sym, repl, stmts)
    var canindex = false
    var idx: int
    #echo "BracketExpr:"
    #echo " ", x.repr
    #echo " ", a.treerepr
    #echo " ", k.repr
    if a.kind==nnkBracket:
      var n = k
      if n.kind==nnkHiddenStdConv: n = n[1]
      if n.kind==nnkIntLit:
        idx = n.intval.int
        canindex = true
      #echo "BracketExpr Bracket:"
      #echo " ", x.repr
      #echo " ", a.repr
      #echo " ", k.repr
      #if not canindex: echo " failed"
    if canindex:
      result = a[idx]
      #echo " ", result.repr
    else:
      result = x.copyNimNode
      result.add a
      result.add k
      for i in 2..<x.len:
        result.add inlineLetsR(x[i], sym, repl, stmts)
  of nnkPragma:
    if x.len==1 and x[0].kind==nnkExprColonExpr and $x[0][0]=="emit":
      template emt(x): untyped =
        {.emit: x.}
      result = getAst(emt(x[0][1][0]))
    else:
      result = x.copyNimNode
      for i in 0..<x.len:
        result.add inlineLetsR(x[i], sym, repl, stmts)
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

import optlet

macro optimizeAst*(a: typed): untyped =
  #echo "optimizeAst in:"
  #echo a.treerepr
  #echo a.repr
  #let ar = a.repr
  #result = a
  #result = optimizeAstR(a)
  #result = inlineLets(a)
  result = optLets(a)
  #echo "optimizeAst out:"
  #echo result.treerepr
  #echo result.repr
  #let rr = result.repr
  #echo "ar == rr: ", ar==rr

macro XoptimizeAst*(a: typed): untyped = a

var flattenArgDebug {.compiletime.} = false
template debugFlattenCallArgs*(b: bool) =
  #bind flattenArgDebug
  static: flattenArgDebug = b

var flattenArgLetCount {.compileTime.} = 0
proc flattenArgP(arg: NimNode): tuple[sl:NimNode,a:NimNode] =
  result.sl = newStmtList()
  if flattenArgDebug:
    echo "flattenArgP: ", arg.kind
    echo arg.repr
  case arg.kind
  of AtomicNodes:
    result.a = arg
  of {nnkStmtList,nnkStmtListExpr}:
    for i in 0..(arg.len-2):
      result.sl.add arg[i]
      #let t = flattenArgP(arg[i])
      #t.sl.copyChildrenTo result.sl
      #result.sl.add t.a
    let r = flattenArgP(arg[^1])
    r.sl.copyChildrenTo result.sl
    result.a = r.a
  of {nnkObjConstr}:
    var a = arg.copy
    for i in 1..<arg.len:
      let r = flattenArgP(arg[i][1])
      r.sl.copyChildrenTo result.sl
      a[i][1] = r.a
    result.a = a
  # FIXME: should transform nnkCallKinds to pass in result var
  of nnkCallKinds:
    var a = arg.copy
    for i in 0..<arg.len:
      a[i] = arg[i].copy
      #let r = flattenArgP(arg[i])
      #r.sl.copyChildrenTo result.sl
      #a[i] = r.a
    let v = genSym(nskLet, "flattenTmp")
    result.sl.add getAst(newLet(v,a))
    result.a = v
  of nnkBracketExpr:
    var a = arg.copy
    let r = flattenArgP(arg[0])
    r.sl.copyChildrenTo result.sl
    a[0] = r.a
    for i in 1..<arg.len:
      #let r = flattenArgP(arg[i])
      #r.sl.copyChildrenTo result.sl
      #a[i] = r.a
      a[i] = arg[i].copy
    result.a = a
  of {nnkDotExpr,nnkDerefExpr,nnkHiddenDeref,nnkHiddenStdConv}:
    var a = arg.copy
    for i in 0..<arg.len:
      let r = flattenArgP(arg[i])
      r.sl.copyChildrenTo result.sl
      a[i] = r.a
    result.a = a
  else:
    #echo "flattenArgP else: ", arg.kind
    #if arg.kind in nnkCallKinds: echo "  ", arg[0].repr
    #let n = genSym(nskLet, "flattenArgP" & $flattenArgLetCount)
    #inc flattenArgLetCount
    #result.sl.add newLetStmt(n, arg)
    #result.a = n
    result.a = arg

macro flattenCallArgs2*(dbg: static[bool], fn: static[string],
                        args: varargs[untyped]): untyped =
  if dbg:
    echo "flattenCallArgs in: ", fn
    echo args.repr
  var sl = newStmtList()
  var call = newCall(ident(fn))
  for i in 0..<args.len:
    var t = flattenArgP(args[i])
    t.sl.copyChildrenTo sl
    call.add t.a
  if sl.len>0:
    result = sl
    result.add call
  else:
    result = call
  if dbg:
    echo "flattenCallArgs out: ", fn
    echo result.repr
macro flattenCallArgs*(args: varargs[untyped]): auto =
  result = newCall(bindSym"flattenCallArgs2")
  result.add ident("false")
  result.add newLit(args[0].repr)
  for i in 1..<args.len:
    result.add args[i]
macro flattenCallArgsD*(args: varargs[untyped]): auto =
  result = newCall(bindSym"flattenCallArgs2")
  result.add ident("true")
  result.add newLit(args[0].repr)
  for i in 1..<args.len:
    result.add args[i]
