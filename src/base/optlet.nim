import macros

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

type OptState = object
    scopes: seq[tuple[a:int,b:int,c:int]]
    letsyms: seq[NimNode]
    ssym: seq[NimNode]
    srepl: seq[NimNode]
    sym: seq[NimNode]
    repl: seq[NimNode]

proc init(os: var OptState) =
  os.scopes.newSeq(0)
  os.letsyms.newSeq(0)
  os.ssym.newSeq(0)
  os.srepl.newSeq(0)
  os.sym.newSeq(0)
  os.repl.newSeq(0)

proc pushScope(os: var OptState) =
  os.scopes.add((os.letsyms.len, os.ssym.len, os.sym.len))

proc popScope(os: var OptState) =
  let (n1,n2,n3) = os.scopes.pop
  os.letsyms.setLen(n1)
  os.ssym.setLen(n2)
  os.srepl.setLen(n2)
  os.sym.setLen(n3)
  os.repl.setLen(n3)

proc isLetSym(os: OptState, x: NimNode): bool =
  let xr = x.repr
  var i = os.letsyms.len-1
  while i>=0:
    if os.letsyms[i].repr == xr: break
    dec i
  if i>=0:
    #echo "found: ", i, " : ", sym[i].repr, " -> ", repl[i].repr
    result = true

proc replSyms(x: NimNode; syms,repl: seq[NimNode]): NimNode =
  let xr = x.repr
  var i = syms.len-1
  while i>=0:
    if syms[i].repr == xr: break
    dec i
  if i>=0:
    #echo "found: ", i, " : ", sym[i].repr, " -> ", repl[i].repr
    result = repl[i]
  else:
    result = x.copyNimNode

proc replSyms(x: NimNode; os: OptState): NimNode =
  replSyms(x, os.sym, os.repl)
proc replSSyms(os: OptState, x: NimNode): NimNode =
  replSyms(x, os.ssym, os.srepl)
#proc replAllSyms(x: NimNode; os: OptState): NimNode =
#  let y = replSyms(x, os.ssym, os.srepl)
#  replSyms(y, os.sym, os.repl)


proc optLetsR(x: NimNode; os: var OptState): NimNode

const letReplKinds = nnkLiterals+{nnkObjConstr,nnkBracket}

proc optimizeObjConstr(x: NimNode; os: var OptState): NimNode =
  x.expectKind(nnkObjConstr)
  result = x.copyNimNode
  result.add x[0]
  var sle = newNimNode(nnkStmtListExpr)
  for i in 1..<x.len:
    let c = x[i]
    c.expectKind(nnkExprColonExpr)
    var t = c.copyNimNode
    t.add c[0]
    var r = optLetsR(c[1], os)
    while r.kind == nnkStmtListExpr:
      for j in 0..(r.len-2):
        sle.add r[j]
      r = r[^1]
    if r.kind==nnkSym: r = replSyms(r, os)
    if r.kind in letReplKinds or (r.kind==nnkSym and os.isLetSym(r)):
      t.add r
    else:
      let s = genSym(nskLet, "optObjConstr" & $c[0])
      sle.add newLetStmt(s, r)
      t.add s
      os.letsyms.add s
    result.add t
  if sle.len>0:
    sle.add result
    result = sle

proc optimizeBracket(x: NimNode; os: var OptState): NimNode =
  x.expectKind(nnkBracket)
  result = x.copyNimNode
  var sle = newNimNode(nnkStmtListExpr)
  for i in 0..<x.len:
    var r = optLetsR(x[i], os)
    while r.kind == nnkStmtListExpr:
      for j in 0..(r.len-2):
        sle.add r[j]
      r = r[^1]
    if r.kind==nnkSym: r = replSyms(r, os)
    if r.kind in letReplKinds or (r.kind==nnkSym and os.isLetSym(r)):
      result.add r
    else:
      let s = genSym(nskLet, "optBracket" & $i)
      sle.add newLetStmt(s, r)
      result.add s
      os.letsyms.add s
  if sle.len>0:
    sle.add result
    result = sle

# main goal is to reduce DotExpr and BracketExpr
proc optLetsR(x: NimNode, os: var OptState): NimNode =
  #echo "new tree"
  #echo x.treeRepr
  case x.kind
  of nnkCommentStmt:
    result = newEmptyNode()
    #result = x
  of nnkLiterals+{nnkNone, nnkEmpty, nnkIdent, nnkOpenSymChoice}:
    result = x
  of nnkSym:
    #result = x.copyNimNode
    result = os.replSSyms(x)

  of nnkIdentDefs:
    result = x.copyNimNode
    #result.add x[0]
    #result.add x[1]
    #result.add x[2]
    for i in 0..<x.len:
      result.add x[i]

  of nnkLetSection:
    result = newStmtList()
    for i in 0..<x.len:
      let c = x[i]
      if c.kind==nnkIdentDefs:
        var r = optLetsR(c[2], os)
        while r.kind == nnkStmtListExpr:
          for j in 0..(r.len-2):
            result.add r[j]
          r = r[^1]
        if r.kind==nnkSym: r = replSyms(r, os)
        #echo "*let: ", c[0].repr, " = ", r.repr, " : ", r.kind
        if r.kind in letReplKinds:
          #echo "let: ", c[0].repr, " = ", r.repr
          #echo c[id][2].treerepr
          os.sym.add c[0]
          #echo "sym: ", sym[^1]
          os.repl.add r
        else:
          if r.kind==nnkSym and os.isLetSym(r):
            os.ssym.add c[0]
            os.srepl.add r
          else:
            os.letsyms.add c[0]
        result.add newLetStmt(c[0], r)
      else:
        echo "error: nnkLetSection expected nnkIdentDefs"
        echo x.treerepr
        error "error"
    if result.len==0: result = newEmptyNode()

  of nnkObjConstr:
    result = optimizeObjConstr(x, os)
    #result = x

  of nnkBracket:
    result = optimizeBracket(x, os)
    #result = x

  of nnkBlockStmt,nnkElifBranch,nnkElse:
    #echo "nnkBlockStmt"
    os.pushScope
    result = x.copyNimNode
    #result.add x[0]
    for i in 0..<x.len:
      result.add optLetsR(x[i], os)
    os.popScope

  of nnkDotExpr:
    var sle = newNimNode(nnkStmtListExpr)
    var o = optLetsR(x[0], os)
    while o.kind == nnkStmtListExpr:
      for j in 0..(o.len-2):
        sle.add o[j]
      o = o[^1]
    if o.kind==nnkSym: o = replSyms(o, os)
    var canindex = false
    var idx: int
    if o.kind == nnkObjConstr:
      # handle Sym or OpenSymChoice
      let ss = (if x[1].kind==nnkSym: x[1] else: x[1][0]).repr
      var i = o.len - 1
      while i>0:
        if o[i][0].repr == ss:
          idx = i
          canindex = true
          break
        dec i
    if canindex:
      result = o[idx][1]
    else:
      if o.kind notin {nnkSym,nnkCall,nnkDotExpr,nnkBracketExpr,nnkHiddenDeref}:
        echo "Opt: DotExpr failed: ", o.kind, ".", x[1]
        echo "  ", o.repr
      result = x.copyNimNode
      result.add o
      result.add x[1]
    if sle.len>0:
      sle.add result
      result = sle

  of nnkBracketExpr:
    var sle = newNimNode(nnkStmtListExpr)
    var a0 = optLetsR(x[0], os)
    var a = a0
    while a.kind == nnkStmtListExpr:
      for j in 0..(a.len-2):
        sle.add a[j]
      a = a[^1]
    if a.kind==nnkSym: a = replSyms(a, os)
    var k0 = optLetsR(x[1], os)
    var k = k0
    while k.kind == nnkStmtListExpr:
      for j in 0..(k.len-2):
        sle.add k[j]
      k = k[^1]
    var n = k
    var canindex = false
    var idx: int
    #echo "BracketExpr:"
    #echo " ", x.repr
    #echo " ", a.treerepr
    #echo " ", k.repr
    if a.kind==nnkBracket:
      if n.kind==nnkHiddenStdConv: n = n[1]
      if n.kind==nnkIntLit:
        idx = n.intval.int
        canindex = true
      #echo "BracketExpr Bracket:"
      #echo " ", x.repr
      #echo " ", a.repr
      #echo " ", k.repr
    if canindex:
      result = a[idx]
      #echo " ", result.repr
    else:
      if a.kind notin {nnkSym,nnkCall,nnkDotExpr,nnkBracketExpr,nnkHiddenDeref}:
        echo "Opt: BracketExpr failed: ", a.kind, "[", n.kind, "]"
      result = x.copyNimNode
      result.add a
      result.add k
    if sle.len>0:
      sle.add result
      result = sle

  of nnkCallKinds:
    result = x.copyNimNode
    result.add x[0]
    var sle = newNimNode(nnkStmtListExpr)
    for i in 1..<x.len:
      var o = optLetsR(x[i], os)
      while o.kind == nnkStmtListExpr:
        for j in 0..(o.len-2):
          sle.add o[j]
        o = o[^1]
      result.add o
    #echo sle.repr
    #echo result.repr
    if sle.len>0:
      sle.add result
      result = sle

  of nnkForStmt:
    result = x.copyNimNode
    result.add x[0].copyNimTree
    result.add x[1].copyNimTree
    result.add optLetsR(x[2], os)  # should be scoped

  of nnkPragma:
    result = x
    #[
    if x.len==1 and x[0].kind==nnkExprColonExpr and $x[0][0]=="emit":
      template emt(x): untyped =
        {.emit: x.}
      result = getAst(emt(x[0][1][0]))
    else:
      result = x.copyNimNode
      for i in 0..<x.len:
        result.add optLetsR(x[i], sym, repl)
    ]#

  else:
    result = x.copyNimNode
    for i in 0..<x.len:
      #echo "Xelse"
      #echo x[i].repr
      result.add optLetsR(x[i], os)

proc optLets*(x: NimNode): NimNode =
  var os: OptState
  os.init
  result = optLetsR(x, os)

macro optimizeAst*(a: typed): untyped =
  echo "optimizeAst in:"
  #echo a.treerepr
  echo a.repr
  #let ar = a.repr
  #result = a
  #result = optimizeAstR(a)
  result = optLets(a)
  echo "optimizeAst out:"
  #echo result.treerepr
  echo result.repr
  #let rr = result.repr
  #echo "ar == rr: ", ar==rr
