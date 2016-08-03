import globals
import alignedMem
import threading
import comms
import layout
import macros
import stdUtils
import basicOps
import profile
import types

type
  FieldObj*[V:static[int],T] = object
    s*:alignedMem[T]
    l*:Layout[V]
    elemSize*:int
  Field*[V:static[int],T] = ref FieldObj[V,T]
  Subsetted*[F,S] = object
    field*:F
    subset*:S
  FieldBinop*[O:static[string],T1,T2] = object
    f1:T1
    f2:T2
  FieldAddSub*[S:static[tuple],T:tuple] = object
    field*:T
  FieldMul*[T:tuple] = object
    field*:T
  Shifted*[T] = object
    field*:T
    dir*:int
    ln*:int
  Field2* = distinct Field
  Field3* = distinct Field
  SomeField* = Field | Subsetted | FieldBinop | FieldAddSub | FieldMul|Shifted
  #SomeField2* = Field | Subsetted | FieldBinop | FieldAddSub|FieldMul|Shifted
  SomeField2* = concept x
    x is SomeField
  SomeAllField* = Field | FieldBinop | FieldAddSub | FieldMul
  SomeAllField2* = Field | FieldBinop | FieldAddSub | FieldMul
  notSomeField* = concept x
    #x isnot SomeField
    x isnot Field
    x isnot Subsetted
    x isnot FieldAddSub
    x isnot FieldMul
    x isnot Shifted
  notSomeField2* = concept x
    x isnot SomeField2

template numberType*(x:Field):untyped = numberType(x[0])

macro fieldAddSub*(sx:static[int],x:Field|FieldMul):auto =
  result = quote do:
    FieldAddSub[(a:`sx`),tuple[a:type(`x`)]](field:(a:`x`))
  #echo result.repr
macro fieldAddSub*(sx:static[int],x:Field|FieldMul,
                   sy:static[int],y:Field|FieldMul):auto =
  result = quote do:
    FieldAddSub[(a:`sx`,b:`sy`),tuple[a:type(`x`),b:type(`y`)]](
      field:(a:`x`,b:`y`) )
  #echo result.repr
macro fieldMul*(x:Field,y:Field2):auto =
  result = quote do:
    FieldMul[tuple[a:type(`x`),b:type(`y`)]](field:(a:`x`,b:`y`))
  #echo result.repr
macro fieldMul*(x:notSomeField,y:SomeField):auto =
  result = quote do:
    FieldMul[tuple[a:type(`x`),b:type(`y`)]](field:(a:`x`,b:`y`))
  #echo result.repr
macro fieldMul*(x:SomeField;y:notSomeField):auto =
  result = quote do:
    FieldMul[tuple[a:type(`x`),b:type(`y`)]](field:(a:`x`,b:`y`))
  #echo result.repr
macro fieldMul*(x:Field,y:Shifted):auto =
  result = quote do:
    FieldMul[tuple[a:type(`x`),b:type(`y`)]](field:(a:`x`,b:`y`))
  #echo result.repr
macro fieldBinop*(o:static[string];x,y:typed):auto =
  result = quote do:
    FieldBinop[`o`,type(`x`),type(`y`)](f1:`x`,f2:`y`)
  #echo result.repr
macro fieldShift*(x:SomeField, d,l:int):auto =
  result = quote do:
    Shifted[type(`x`)](field:`x`,dir:`d`,ln:`l`)
  #echo result.repr

proc new*[V:static[int],T](x:var FieldObj[V,T]; l:Layout[V]) =
  x.l = l
  x.s.new(l.nSitesOuter)
  #fence()
  x.elemSize = sizeOf(T)
proc new*[V:static[int],T](x:var Field[V,T]; l:Layout[V]) =
  x.new()
  new(x[], l)
proc new*[V:static[int],T](x:var FieldObj[V,T]; y:Field) = x.new(y.l)
proc new*[V:static[int],T](x:var Field[V,T]; y:Field) = x.new(y.l)
proc newField*[V:static[int],T](l:Layout[V]; t:typedesc[T]):Field[V,T] =
  result.new(l)
proc newOneOf*(x:Field):auto =
  var r:type(x)
  r.new(x.l)
  r

template `[]`*(x:Field; i:int):expr = x.s[i]
template `[]`*(x:Subsetted; i:int):expr = x.field[i]
template `[]`*(x:SomeField; st:string):untyped =
  Subsetted[type(x),type(st)](field:x,subset:st)
template `[]`*(x:SomeField; st:Subset):untyped =
  Subsetted[type(x),type(st)](field:x,subset:st)
template even*(x:Field):untyped = x["even"]
template odd*(x:Field):untyped = x["odd"]
template all*(x:Field):untyped = x["all"]
template `even=`*(x:Field; y:any):untyped = assign(x["even"], y)
template `odd=`*(x:Field; y:any):untyped =
  mixin assign
  assign(x["odd"], y)
template `all=`*(x:Field; y:any):untyped = assign(x["all"], y)
template indexField*(x:notSomeField; y:int):untyped = x
template indexField*(x:Field; y:int):untyped = x[y]
template indexField*(x:Adjointed[SomeField]; y:int):untyped = x[][y].adj
#macro indexField*(x:FieldBinop; y:int):auto =
#  echo x.treeRepr
#  var xx = x
#  echo xx.kind
#  if xx.kind == nnkStmtListExpr: xx = x[2]
#  let op = ident($xx[1])
#  let f1 = xx[2]
#  let f2 = xx[3]
#  result = quote do:
#    `op`(indexField(`f1`,`y`),indexField(`f2`,`y`))
#  echo result.repr
macro indexFieldM*(x:FieldAddSub, sx:tuple, y:int):auto =
  #echo x.repr
  #echo x.treeRepr
  #echo sx.treeRepr
  #var p = x[1][1]
  #while p.kind != nnkPar:
  #  echo p.kind
  #  p = p[1]
  #var t = p[0][1]
  result = quote do:
    indexField(`x`.field[0],`y`)
  if sx[0][1].intVal == -1:
    result = quote do:
      -`result`
  for i in 1..<sx.len:
    #t = p[i][1]
    if sx[i][1].intVal == 1:
      result = quote do:
        `result` + indexField(`x`.field[`i`],`y`)
    else:
      result = quote do:
        `result` - indexField(`x`.field[`i`],`y`)
  #echo result.repr
template indexField*(x:FieldAddSub, y:int):expr = indexFieldM(x, x.S, y)
template `[]`*(x:FieldAddSub, y:int):expr = indexField(x, y)
macro indexFieldM*(x:FieldMul; tx:typedesc; y:int):auto =
  #echo x.treeRepr
  #echo tx.getType.treeRepr
  let nt = tx.getType[1].len - 1
  #var p = x[1][1]
  #while p.kind != nnkPar:
  #  echo p.kind
  #  p = p[1]
  #var t = p[0][1]
  result = quote do:
    indexField(`x`.field[0],`y`)
  for i in 1..<nt:
    #t = p[i][1]
    result = quote do:
      `result` * indexField(`x`.field[`i`],`y`)
  #echo result.repr
template indexField*(x:FieldMul, y:int):expr = indexFieldM(x, x.T, y)
template `[]`*(x:FieldMul, y:int):expr = indexField(x, y)

template l*(x:FieldAddSub):untyped = x.field[0].l

template itemsI*(n0,n1:int):untyped =
  let n = n1 - n0
  var ti0 = n0 + ((threadNum*n) div numThreads)
  var ti1 = n0 + (((threadNum+1)*n) div numThreads)
  #echo "ti0: ", ti0, "  ti1: ", ti1
  var i = ti0
  while i < ti1:
    yield i
    inc(i)
iterator items*(l:Layout):int {.inline.} =
  let n = l.nSitesOuter
  itemsI(0, n)
iterator sites*(l:Layout):int {.inline.} =
  let n = l.nSites
  itemsI(0, n)
iterator items*(s:Subset):int {.inline.} =
  let n0 = s.lowOuter
  let n1 = s.highOuter
  itemsI(n0, n1)
iterator sites*(s:Subset):int {.inline.} =
  let n0 = s.low
  let n1 = s.high
  itemsI(n0, n1)
#iterator all*(x:Field):int {.inline.} =
#  let n = x.l.nSitesOuter
#  itemsI(0, n)
iterator items*(x:Field):int {.inline.} =
  let n = x.l.nSitesOuter
  #echo "n: ", n
  itemsI(0, n)
iterator items*(x:Subsetted):int {.inline.} =
  when x.subset is string:
    let s = getSubset(x.field.l, x.subset)
  else:
    let s = x.subset
  let n0 = s.lowOuter
  let n1 = s.highOuter
  #echo "n0: ", n0, " n1: ", n1
  itemsI(n0, n1)
iterator items*(x:FieldAddSub):int {.inline.} =
  let n = x.field[0].l.nSitesOuter
  itemsI(0, n)
iterator items*(x:FieldMul):int {.inline.} =
  let n = x.field[0].l.nSitesOuter
  itemsI(0, n)
#macro filter*(x:SomeField; s:string; pred:expr):Mask {.inline.} =
#iterator filter*(x:SomeField; s:string; pred:expr):Mask {.inline.} =
  ###

import types
export types

#proc `{}`*[V:static[int],T](f:Field[V,T]; i:int):Masked[T] =
proc `{}`*(f:Field; i:int):auto =
  let e = i div f.l.V
  let l = i mod f.l.V
  let mask = 1 shl l
  #echo i, " ", e, " ", r, " ", mask
  #result.pobj = f[e].addr
  #result.mask = mask
  #result = Masked[f.T](pobj:f[e], mask:mask)
  #var r:Masked[f.T]
  #r.pobj = f[e].addr
  #r.mask = mask
  #r
  #echoImm: "{}"
  result = masked(f[e], mask)
  
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


template indexField(x:Shifted, y:int):untyped = 0
proc applyOp1(x,y:NimNode; op:string):auto =
  let o = ident(op)
  result = quote do:
    let t = `y`
    #echoImm `x`[0] is VMconcept1
    #echoImm t isnot VMconcept2
    for e in `x`:
      mixin `o`
      #mixin isMatrix
      #echoAll isMatrix(`x`[e])
      `o`(`x`[e], t)
proc applyOp2(x,y:NimNode; ty:typedesc; op:string):auto =
  #echo ty.getType.treeRepr
  #echo ty.getType.getImpl.treeRepr
  let o = ident(op)
  result = quote do:
    let xx = `x`
    let yy = `y`
    for e in xx:
      `o`(xx[e], indexField(yy, e))
template makeOps(op,f,fM,s:untyped):untyped =
  macro f*(x:Subsetted; y:notSomeField2):auto = applyOp1(x,y,s)
  macro f*(x:Subsetted; y:SomeField2):auto = applyOp2(x,y,int,s)
  macro fM*(x:Field; y:notSomeField; ty:typedesc):auto = applyOp1(x,y,s)
  macro fM*(x:Field; y:SomeField; ty:typedesc):auto = applyOp2(x,y,ty,s)
  template f(x:Field; y:any):untyped =
    #when declaredInScope(subsetObject):
    when declared(subsetObject):
      #echo "subsetObj" & s
      f(x[subsetObject], y)
    elif declared(subsetString):
      #echo "subsetString" & s
      f(x[subsetString], y)
    else:
      fM(x, y, y.type)
  when profileEqns:
    template op*(x:Field; y:any):untyped =
      block:
        tic(-2)
        f(x, y)
        toc(-2)
    template op*(x:Subsetted; y:any):untyped =
      block:
        tic(-2)
        f(x, y)
        toc(-2)
  else:
    template op*(x:Field; y:any):untyped = f(x, y)
    #template op*(x:var Field; y:any):untyped = f(x, y)
    template op*(x:Subsetted; y:any):untyped = f(x, y)
makeOps(`:=`, assign, assignM, "assign")
makeOps(`+=`, iadd, iaddM, "iadd")
makeOps(`-=`, isub, isubM, "isub")

proc mul*(r:Field; x:Field2; y:Field3) =
  mixin mul
  for e in r:
    mul(r[e], x[e], y[e])

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
template norm2*(f:SomeAllField):expr =
  when declared(subsetObject):
    #echo "subsetObj" & s
    norm2P(f[subsetObject])
  elif declared(subsetString):
    #echo "subset norm2"
    norm2P(f[subsetString])
  else:
    norm2P(f)
template norm2*(f:Subsetted):expr = norm2P(f)

proc dotP*(f1:SomeField; f2:SomeField2):auto =
  mixin dot, idot, simdSum, items
  #var d:type(dot(f1[0],f2[0]))
  var d:type(toDouble(dot(f1[0],f2[0])))
  let t1 = f1
  let t2 = f2
  for x in items(t1):
    idot(d, t1[x], t2[x])
  result = simdSum(d)
  #threadSum(result)
  #rankSum(result)
  threadRankSum(result)
template dot*(f1:SomeAllField; f2:SomeAllField2):expr =
  when declared(subsetObject):
    #echo "subsetObj" & s
    dotP(f1[subsetObject], f2)
  elif declared(subsetString):
    dotP(f1[subsetString], f2)
  else:
    dotP(f1, f2)
template dot*(f1:Subsetted; f2:SomeAllField2):expr = dotP(f1, f2)
proc redotP*(f1:SomeField; f2:SomeField2):auto =
  tic()
  mixin redot, iredot, simdSum, items, toDouble
  #var d:type(redot(f1[0],f2[0]))
  var d:type(toDouble(redot(f1[0],f2[0])))
  let t1 = f1
  let t2 = f2
  for x in items(t1):
    iredot(d, t1[x], t2[x])
  toc("local")
  result = simdSum(d)
  toc("simd sum")
  #threadBarrier()
  #toc("thread barrier")
  #threadSum(result)
  #toc("thread sum")
  #rankSum(result)
  #toc("rank sum")
  threadRankSum(result)
  toc("thread rank sum")
template redot*(f1:SomeAllField; f2:SomeAllField2):expr =
  when declared(subsetObject):
    #echo "subsetObj redot"
    redotP(f1[subsetObject], f2)
  elif declared(subsetString):
    redotP(f1[subsetString], f2)
    #echo "subset redot"
  else:
    redotP(f1, f2)
template redot*(f1:Subsetted; f2:SomeAllField2):expr = redotP(f1, f2)

proc trace*(m:Field):auto =
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
template sum*(f:SomeAllField):expr =
  when declared(subsetString):
    sumP(f[subsetString])
  else:
    sumP(f)
template sum*(f:Subsetted):expr = sumP(f)

template `-`*(x:SomeField):expr = fieldAddSub(-1,x)
template `+`*(x:SomeField,y:SomeField2):expr = fieldAddSub(1,x,1,y)
template `+`*(x:SomeField,y:notSomeField2):expr = fieldAddSub(1,x,1,y)
template `+`*(x:notSomeField,y:SomeField2):expr = fieldAddSub(1,x,1,y)
template `-`*(x:SomeField,y:SomeField2):expr = fieldAddSub(1,x,-1,y)
template `-`*(x:SomeField,y:notSomeField2):expr = fieldAddSub(1,x,-1,y)
template `-`*(x:notSomeField,y:SomeField2):expr = fieldAddSub(1,x,-1,y)
template `*`*(x:SomeField,y:SomeField2):expr = fieldMul(x,y)
template `*`*(x:SomeField,y:notSomeField2):expr = fieldMul(x,y)
template `*`*(x:notSomeField,y:SomeField2):expr = fieldMul(x,y)
template `/`*(x:SomeField,y:SomeNumber):expr =
  let t = 1.0/y
  fieldMul(x,t)
#template makeBinop(f:untyped,s:string):untyped =
#  template `f`*(x:SomeField; y:SomeField2):expr =
#    echoImm x
#    echoImm y
#    fieldBinop(`s`, x, y)
#makeBinop(`+`,"+")
#makeBinop(`-`,"-")
#makeBinop(`*`,"*")
#makeBinop(`/`,"/")
#template `*`*(x:notSomeField; y:SomeField2):untyped =
#  fieldBinop("*", x, y)

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

import shifts

#template assign*(x:var Field; y:Shifted):untyped =
#  shift(x, y.dir, y.len, y.field)

when isMainModule:
  import qexInternal
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
