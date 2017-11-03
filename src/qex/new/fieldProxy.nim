import macros
import base
import strUtils

# f := x + y
# f.even = x + y
# f[e] = x + y

type
  FieldProxy*[T] = object
    field*: T
  FieldProxy2*[T] = FieldProxy[T]
  FieldProxy3*[T] = FieldProxy[T]
  FieldProxy4*[T] = FieldProxy[T]
  #SiteIndex* = distinct int
  #VSiteIndex* = distinct int

template `[]`*(x: FieldProxy): untyped =
  x.field

proc isFieldProxy(x: NimNode): bool =
  let t = x.getTypeInst
  echo t.treerepr
  if t.kind==nnkSym:
    if $t.symbol=="FieldProxy":
      result = true
  else:
    for i in 0..<t.len:
      result = isFieldProxy(t[i])
      if result: break
  if not result: echo t.treerepr

template indexFieldProxy*(x: FieldProxy, y: typed): untyped =
  x[][y]

const opNameOps = ['+','-','*','/']
const opNameNames = ["plus","minus","times","divd"]
proc opToName(op: string): string =
  result = ""
  for c in op:
    var i = opNameOps.indexOf(c)
    if i>=opNameOps.len:
      echo "unknown op: ", c, " in ", op
      quit 1
    result &= opNameNames[i]
proc nameToOp(op: string): string =
  result = ""
  var pos = 0
  while pos<op.len:
    var i = 0
    while i<opNameNames.len and not op.continuesWith(opNameNames[i],pos): inc i
    if i>=opNameNames.len:
      echo "unknown name: ", op[pos..^1], " in ", op
      quit 1
    result &= opNameOps[i]
    pos += opNameNames[i].len
proc toSiteLocalName(fn: string): string =
  var f = fn
  if f[0] notin IdentStartChars:
    f = "op_" & opToName(f)
  f &= "_SiteLocal"
  result = f
proc fromSiteLocalName(fn: string): string =
  var f = fn[0..^11]
  if f.startsWith("op_"):
    f = nameToOp(f[3..^1])
  result = f

# rename site local calls (f_SiteLocal)
# split equations
# - rename nonlocal calls (f_Impl)
#
proc isSiteLocal(x: string): bool =
  x.endsWith("_SiteLocal")
proc splitFieldExpression(x: NimNode): tuple[defs:NimNode,expr:NimNode] =
  var defs = newStmtList()
  var expr: NimNode
  if x.kind in CallNodes:
    if isSiteLocal($x.name):
      expr = newNimNode(x.kind)
      expr.add x.name
      for i in 1..<x.len:
        let (df,ex) = splitFieldExpression(x[i])
        defs.add df
        expr.add ex
    else: # FIXME: make recursive & rename nonlocal calls
      #template myLet(a,b) =
      #  let a = b
      #let t = genSym(nskLet, "tmpField")
      template myLet(a,b) =
        var a{.noInit.}: type(b)
        a := b
      let t = genSym(nskVar, "tmpField")
      defs.add getAst(myLet(t,x))
      expr = t
  else:
    expr = x
  result = (defs, expr)

macro indexFieldProxy*(x: FieldProxy{call}, y: typed): untyped =
  let f = $x[0]
  echo "site index: ", f
  if f.endsWith("_SiteLocal"):
    let f2 = fromSiteLocalName(f)
    result = newCall(ident(f2))
    for i in 1..<x.len:
      let xi = x[i]
      result.add newCall(!"indexFieldProxy", xi, y)
  else:
    #echo "not found: ", f
    result = newCall(ident"[]",newCall(ident"[]",x),y)

template `[]`*(x: FieldProxy, y: typed): untyped =
  indexFieldProxy(x, y)

template `[]=`*(x: FieldProxy, y: typed, z: typed): untyped =
  x[][y] = z

proc newFieldProxy*(x: FieldProxy, y: typedesc): auto =
  mixin newFieldImpl
  newFieldImpl(x[], y)

template `len`*(x: FieldProxy): untyped = x[].len

#macro `=`*[T](r: var FieldProxy[T], x: FieldProxy[T]): auto =
#  template doAssign(r: FieldProxy, x: FieldProxy): untyped =
#    mixin indices
#    for i in indices(r[]):
#      op(r[i], x[i])
#  let (df,ex) = splitFieldExpression(x)
#  result = df
#  result.add getAst(doAssign(r, ex))
#  echo result.repr

template assignOverloads(op: untyped) {.dirty.} =
  macro op*(r: FieldProxy, x: FieldProxy2): auto =
    template doAssign(r: FieldProxy, x: FieldProxy2): untyped =
      mixin indices
      for i in indices(r[]):
        op(r[i], x[i])
    let (df,ex) = splitFieldExpression(x)
    result = df
    result.add getAst(doAssign(r, ex))
    echo result.repr

assignOverloads(`:=`)
assignOverloads(`+=`)
assignOverloads(`-=`)
assignOverloads(`*=`)
assignOverloads(`/=`)

template unaryOverloads(fn: untyped) {.dirty.} =
  proc fn*(x: FieldProxy): auto {.noInit.} =
    mixin newFieldImpl, indices
    result = newFieldImpl(x[], type(fn(x[0])))
    for i in indices(result[]):
      result[i] = fn(x[i])

unaryOverloads(`-`)
unaryOverloads(adj)
unaryOverloads(toSingle)
unaryOverloads(toDouble)

template binaryOverloads2(fn,slfn: untyped) {.dirty.} =
  proc slfn*(x: FieldProxy, y: FieldProxy2): auto {.noInit.} =
    echo "unoptimized siteLocal"
    mixin newFieldImpl, indices
    result = newFieldImpl(x[], type(fn(x[0],y[0])))
    for i in indices(result[]):
      result[i] = fn(x[i], y[i])
  template fn*(x: FieldProxy, y: FieldProxy2): untyped = slfn(x,y)
macro binaryOverloads(fn: untyped): auto =
  let f = $fn.name
  let sln = toSiteLocalName(f)
  result = getAst(binaryOverloads2(fn, ident(sln)))

binaryOverloads(`+`)
binaryOverloads(`-`)
binaryOverloads(`*`)
binaryOverloads(`/`)

template assignOverloadScalar(op,typ: untyped) {.dirty.} =
  template op*(r: FieldProxy, x: typ): untyped =
    for i in indices(r[]):
      op(r[i], x)

template binaryOverloadScalar(fn,typ: untyped) {.dirty.} =
  proc fn*(x: FieldProxy, y: typ): auto {.noInit.} =
    mixin newFieldImpl, indices
    result = newFieldImpl(x[], type(fn(x[0],y)))
    for i in indices(result[]):
      result[i] = fn(x[i], y)
  proc fn*(x: typ, y: FieldProxy2): auto {.noInit.} =
    mixin newFieldImpl, indices
    result = newFieldImpl(y[], type(fn(x,y[0])))
    for i in indices(result[]):
      result[i] = fn(x, y[i])

template fieldScalarOverloads*(typ: untyped) {.dirty.} =
  bind assignOverloadScalar, binaryOverloadScalar
  template indexFieldProxy*(x: typ, y: typed): untyped = x
  assignOverloadScalar(`:=`, typ)
  assignOverloadScalar(`+=`, typ)
  assignOverloadScalar(`-=`, typ)
  assignOverloadScalar(`*=`, typ)
  assignOverloadScalar(`/=`, typ)
  binaryOverloadScalar(`+`, typ)
  binaryOverloadScalar(`-`, typ)
  binaryOverloadScalar(`*`, typ)
  binaryOverloadScalar(`/`, typ)

proc `$`*(x: FieldProxy): string =
  #let n = x.len
  #result =  "FieldProxy.len: " & $(x.len) & "\n"
  #result &= "FieldProxy.first = " & $x[0] & "\n"
  #result &= "FieldProxy[" & $n & "] = " & $x[n]
  $(x[])


#[
template `[]`*(x: FieldProxy, s: string): untyped =
  Subsetted[type(x),type(st)](field:x,subset:st)
template `[]`*(x:SomeField; st:Subset):untyped =
  Subsetted[type(x),type(st)](field:x,subset:st)

template numberType*(x: FieldProxy): untyped = numberType(x[])

template even*(x:Field):untyped = x["even"]
template odd*(x:Field):untyped = x["odd"]
template all*(x:Field):untyped = x["all"]
template `even=`*(x:Field; y:any):untyped = assign(x["even"], y)
template `odd=`*(x:Field; y:any):untyped =
  mixin assign
  assign(x["odd"], y)
template `all=`*(x:Field; y:any):untyped = assign(x["all"], y)

#proc `$`*(x:Field):string =
#  $(x[0])
proc `$`*(x:Field):string =
  mixin `$`
  let l = x.l
  result = ""
  for i in 0..<l.nSites:
    result.add $l.coords[0][i]
    for j in 1..<l.nDim:
      result.add " " & $l.coords[j][i]
    result.add ": "
    result.add $(x{i})
    if i<l.nSites-1: result.add "\n"

proc norm2P*(f:SomeField):auto =
  tic()
  mixin norm2, inorm2, simdSum, items, toDouble
  #var n2:type(norm2(f[0]))
  var n2:type(toDouble(norm2(f[0])))
  #echo n2
  #let t = f
  for x in items(f):
    inorm2(n2, f[x])
  toc("norm2 local")
  #echoAll n2
  result = simdSum(n2)
  toc("norm2 simd sum")
  #echoAll result
  #threadSum(result)
  #toc("norm2 thread sum")
  #rankSum(result)
  #toc("norm2 rank sum")
  threadRankSum(result)
  toc("norm2 thread rank sum")
template norm2*(f:SomeAllField):untyped =
  when declared(subsetObject):
    #echo "subsetObj" & s
    norm2P(f[subsetObject])
  elif declared(subsetString):
    #echo "subset norm2"
    norm2P(f[subsetString])
  else:
    norm2P(f)
template norm2*(f:Subsetted):untyped = norm2P(f)

proc dotP*(f1:SomeField; f2:SomeField2):auto =
  mixin dot, idot, simdSum, items, toDouble, eval
  #var d:type(dot(f1[0],f2[0]))
  var d:type(eval(toDouble(dot(f1[0],f2[0]))))
  let t1 = f1
  let t2 = f2
  for x in items(t1):
    idot(d, t1[x], t2[x])
  result = simdSum(d)
  #threadSum(result)
  #rankSum(result)
  threadRankSum(result)
template dot*(f1:SomeAllField; f2:SomeAllField2):untyped =
  when declared(subsetObject):
    #echo "subsetObj" & s
    dotP(f1[subsetObject], f2)
  elif declared(subsetString):
    dotP(f1[subsetString], f2)
  else:
    dotP(f1, f2)
template dot*(f1:Subsetted; f2:SomeAllField2):untyped = dotP(f1, f2)

proc redotP*(f1:SomeField; f2:SomeField2):auto =
  tic()
  mixin redot, iredot, simdSum, items, toDouble, eval
  #var d:type(redot(f1[0],f2[0]))
  var d: type(eval(toDouble(redot(f1[0],f2[0]))))
  let t1 = f1
  let t2 = f2
  for x in items(t1):
    #iredot(d, t1[x], t2[x])
    d += redot(t1[x], t2[x])
  toc("redot local")
  result = simdSum(d)
  toc("redot simd sum")
  #threadBarrier()
  #toc("thread barrier")
  #threadSum(result)
  #toc("thread sum")
  #rankSum(result)
  #toc("rank sum")
  threadRankSum(result)
  toc("redot thread rank sum")
template redot*(f1:SomeAllField; f2:SomeAllField2):untyped =
  when declared(subsetObject):
    #echo "subsetObj redot"
    redotP(f1[subsetObject], f2)
  elif declared(subsetString):
    redotP(f1[subsetString], f2)
    #echo "subset redot"
  else:
    redotP(f1, f2)
template redot*(f1:Subsetted; f2:SomeAllField2):untyped = redotP(f1, f2)

proc trace*(m:SomeField):auto =
  mixin trace, simdSum
  var tr:type(trace(m[0]))
  for x in m:
    tr += m[x].trace
  #echo tr
  result = simdSum(tr)
  #threadSum(result)
  #rankSum(result)
  threadRankSum(result)

proc sumP*(f:SomeField):auto =
  mixin inc, simdSum, items
  var s:type(f[0])
  let t = f
  for x in items(t):
    iadd(s, t[x])
  result = simdSum(s)
  threadSum(result)
  rankSum(result)
template sum*(f:SomeAllField):untyped =
  when declared(subsetString):
    sumP(f[subsetString])
  else:
    sumP(f)
template sum*(f:Subsetted):untyped = sumP(f)

template onSubset*(s:string; body:untyped):untyped =
  block:
    let subsetString{.inject.} = s
    threadBarrier()
    body
    threadBarrier()
template onNoSync*(s:Subset; body:untyped):untyped =
  block:
    let subsetObject{.inject.} = s
    body

when isMainModule:
  import qex/base/qexInternal
  qexInit()
  echo "rank ", myRank, "/", nRanks
  echo threadNum, "/", numThreads
  var lat = [8,8,8,8]
  #var lat = [16,16,8,8]
  #var lo = newLayout(lat)
  var lo = newLayout(lat, 1)

  var x = lo.newField(float)
  var y = lo.newField(float)
  var z = lo.newField(float)
  let vol = lo.physVol.float

  threads:
    x := 1
    y := -x
    z := x - y
    y := 2*x
    z := x * y
    z := x + 2*y
    z += 3*x + 4*y
    onSubset "even":
      z := 1
      z += 2*x + 3*y
    onSubset "odd":
      z := 2
    #threadBarrier()
    z.even := 0
    z.odd = 0
    threadBarrier()
    z.all += 3*y
    #z := fieldShift(y, 0, 1)
    #z := x * fieldShift(y, 0, 1)
    echo z.norm2/vol
    echo z.even.norm2/vol
    onSubset "even":
      echo z.norm2/vol
    echo z.odd.norm2/vol
    onSubset "odd":
      echo z.norm2/vol
      echo z.even.norm2/vol
      echo "dot: ", x.dot(y)/vol
    echo "sum: ", z.sum/vol
    echo "sum: ", z.even.sum/vol
    echo "sum: ", z.odd.sum/vol
    #echo x
    #echo y
    #echo z
    echo z[lo.nEven]

  qexFinalize()

# TODO
# future
# fma
# aggregate AddSub, Mul

]#


when isMainModule:
  type
    FieldObj[T] = object
      v: array[100,T]
    Field[T] = FieldProxy[FieldObj[T]]
    FieldIndexObj = object
    FieldIndexType = FieldProxy[FieldIndexObj]

  proc newFieldImpl[T](x: FieldObj, y: T): Field[T] = discard
  proc newFieldImpl[T](x: FieldObj, y: typedesc[T]): Field[T] = discard
  proc newField[T](): Field[T] = discard

  template `[]`(x: FieldObj, y: untyped): untyped =
    x.v[y]
  template `[]=`(x: FieldObj, y: untyped, z: untyped): untyped =
    x.v[y] = z
  template `len`*(x: FieldObj): untyped = x.v.len
  proc `$`*(x: FieldObj): string =
    let n = x.len - 1
    result  = "FieldObj[0]: " & $x[0] & "\n"
    result &= "FieldObj[" & $n & "]: " & $x[n]
  iterator indices(x: FieldObj): int =
    var i = 0
    var n = x.len
    while i<n:
      yield i
      inc i

  fieldScalarOverloads(SomeNumber)

  proc newFieldImpl[T](x: FieldIndexObj, y: T): Field[T] = discard
  proc newFieldImpl[T](x: FieldIndexObj, y: typedesc[T]): Field[T] = discard
  template `[]`(x: FieldIndexObj, y: untyped): untyped = y
  template fieldIndexVal(): untyped = FieldIndexType()
  template `:=`*(x: float, y: float) = x = y
  template `+`*(x: float, y: int): float = x + float(y)

  var x,y,z: Field[float]

  x[0] = (y+z)[0]
  x[0] = (y+z+y+z)[0]
  x[0] = (y+z+y+1)[0]

  x := y + fieldIndexVal()
  echo $x
  x := y + z + 1 + fieldIndexVal()
  echo $x
