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

const Keywords = ["addr"]

proc prepareVars*(n:NimNode, deref:proc): seq[NimNode] =
  # get a list of vars and new symbols to replace them, using let binding for now XXX
  #     <- [(id, varsym, letptrsym), ...]
  # the symbols in n are changed
  #echo "### prepareVars: ",n.treerepr
  var ignoreStack = newseq[NimNode]()
  var openvars = newseq[NimNode]()
  proc go(n:NimNode) =
    # ign is a stack for ignoring lexical bindings: [(outer,...), (inner,...), ...]
    #echo "go get: ",n.repr
    #block:
    #  var ignstr = ""
    #  for c in ignoreStack: ignstr &= ("\n" & c.repr)
    #  echo "ign has: ",ignstr
    var newscope = false
    if n.kind in {nnkBlockStmt, nnkBlockExpr, nnkIfExpr, nnkElifExpr, nnkElseExpr,
        nnkIfStmt, nnkElifBranch, nnkElse, nnkCaseStmt, nnkOfBranch,
        nnkWhileStmt, nnkForStmt} + RoutineNodes:
      # New lexical scope
      newscope = true
      #ignoreStack.add newPar()
      ignoreStack.add newNimNode(nnkTupleConstr)
    for i in 0..<n.len:
      #echo "### ",n[i].lisprepr
      case n[i].kind
      of {nnkVarSection,nnkLetSection}:
        for cc in n[i]:
          for c in 0..cc.len-2:
            case cc[c].kind
            of nnkPragmaExpr:
              ignoreStack[^1].add cc[c][0]
            of nnkBracketExpr:
              for ccc in cc[c]:
                ignoreStack[^1].add ccc
            else:
              ignoreStack[^1].add cc[c]
      of nnkOpenSymChoice:
        if n.kind in Callnodes: continue
      of Callnodes:
        if n[i][0].kind in {nnkSym, nnkIdent}:
          var newid = true
          for c in ignoreStack[0]:
            if c == n[i][0]:
              newid = false
              break
          if newid:
            ignoreStack[0].add n[i][0]
      of nnkForStmt:
        #echo n[i].treerepr
        ignoreStack[^1].add n[i][0]
      of nnkWhenStmt:
        #echo n[i].treerepr
        if n[i][0][0].kind in {nnkIdent,nnkSym}:
          ignoreStack[^1].add n[i][0][0]
      of nnkPragma:
        continue
      of nnkCast:
        n[i][1].go
        continue
      of {nnkSym, nnkIdent}:
        if n.kind == nnkDotExpr and i > 0: continue
        if n[i].repr in Keywords: continue
        var ignore = false
        for cc in ignoreStack:
          for c in cc:
            if c.eqIdent n[i]:
              ignore = true
              break
          if ignore: break
        if not ignore:
          var newvar = true
          for c in openvars:
            if c[0].eqIdent n[i]:
              #n[i] = newcall("gpuVarPtr",c[1],c[2])
              n[i] = deref(c[0],c[1],c[2])
              newvar = false
              break
          #echo "EXPR: ",n.lisprepr
          #echo "ID:   ",n[i].repr,"  newvar: ",newvar.repr
          #var rs = ""
          #for c in openvars:
          #  rs &= "  " & c.repr
          #echo "RES:  ",rs
          if newvar:
            let nv = gensym(nskvar, "gpu_" & $n[i])
            #let np = gensym(nsklet, "gpu_ptr_" & $n[i])
            ignoreStack[0].add nv
            #ignoreStack[0].add np
            #openvars.add newpar(n[i], nv, np)
            let k = newLit openvars.len
            openvars.add newNimNode(nnkTupleConstr).add(n[i], nv, k)
            #n[i] = newcall("gpuVarPtr",nv,np)
            #n[i] = newcall("getGpu",n[i],nv)
            n[i] = deref(n[i],nv,k)
        continue
      else:
        discard
      n[i].go
    if newscope: ignoreStack.setLen(ignoreStack.len-1)
  #ignoreStack.add newPar(ident"gpuVarPtr")
  ignoreStack.add newNimNode(nnkTupleConstr)
  n.go
  openvars

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
