import base/globals
import base/alignedMem
import base
#import threading
#import comms
import layout/layoutX
import macros
#import stdUtils
#import basicOps
#import profile
import maths
import maths/types

type FieldOps* = enum
  foNeg, foAdd, foSub, foMul, foDiv, foAdj, foToSingle, foToDouble

type
  FieldObj*[V:static[int],T] = object
    s*:alignedMem[T]
    l*:Layout[V]
    elemSize*:int
  Field*[V:static[int],T] = ref FieldObj[V,T]
  FieldArray*[V:static[int],T] = object  ## share a single alignedMem
    shape*:seq[int]
    arr:seq[Field[V,T]]
  Subsetted*[F,S] = object
    field*:F
    subset*:S
  FieldUnop*[Op: static[FieldOps], T1] = object
    f1*: T1
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
  #Field2* = distinct Field
  #Field3* = distinct Field
  Field2*[V:static[int],T] = Field[V,T]
  Field3*[V:static[int],T] = Field[V,T]
  SomeField* = Field | Subsetted | FieldBinop | FieldAddSub | FieldMul|Shifted | FieldUnop
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
    x isnot FieldUnop
  notSomeField2* = concept x
    x isnot SomeField2

template elemType*(x:Field):typedesc = evalType(x[0])
template elemType*[V:static[int],T](x:typedesc[Field[V,T]]):typedesc = evalType(T)
template numberType*(x:Field):untyped = numberType(x[0])

#template fieldUnop*(o: static[FieldOps], x: SomeField): auto =
#  FieldUnop[o,type(x)](f1: x)
#macro fieldUnop*(o: static[FieldOps], x: SomeField): auto =
#  result = quote do:
#    FieldUnop[FieldOps(`o`),type(`x`)](f1: `x`)
template fieldUnop*(o: FieldOps, x: SomeField): untyped =
  FieldUnop[FieldOps(o),type(x)](f1: x)
macro fieldAddSub*(sx:static[int],x:auto):auto =
  result = quote do:
    FieldAddSub[(a:`sx`),tuple[a:type(`x`)]](field:(a:`x`))
  #echo result.repr
macro fieldAddSub*(sx:static[int],x:auto,
                   sy:static[int],y:auto):auto =
  result = quote do:
    FieldAddSub[(a:`sx`,b:`sy`),tuple[a:type(`x`),b:type(`y`)]](
      field:(a:`x`,b:`y`) )
  #echo result.repr
macro fieldMul*(x:SomeField,y:SomeField2):auto =
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
template adjImpl*(x: SomeField): untyped =
  fieldUnop(foAdj, x)
template toSingleImpl*(x: SomeField): untyped =
  fieldUnop(foToSingle, x)
template toDoubleImpl*(x: SomeField): untyped =
  fieldUnop(foToDouble, x)

template eval*[F:Field](x: typedesc[F]): typedesc =
  Field[F.V,eval(type F.T)]
template eval*[F:FieldObj](x: typedesc[F]): typedesc =
  FieldObj[F.V,eval(type F.T)]
template evalType*[F:Field](x: typedesc[FieldUnop[foToSingle,F]]): typedesc =
  mixin toSingle
  Field[F.V,eval(toSingle(type F.T))]

proc new*[V:static[int],T](x:var FieldObj[V,T]; l:Layout[V]) =
  # remember to change newFieldArray if the following changes
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
proc newOneOf*(x: Field): auto =
  var r: type(x)
  r.new(x.l)
  r
template l*(x: FieldUnop): untyped = x.f1.l
proc newOneOf*(x: FieldUnop): auto =
  var r: evalType(x)
  r.new(x.l)
  r
template new*(x: typedesc[Field], l: Layout): untyped =
  newField(l, x.T)

proc newFarrElem[V:static[int],T](f:var Field[V,T]; l:Layout[V]; s:alignedMem[T]; offset:int) =
  f.new()
  f.l = l
  f.s = s
  f.s.data = cast[typeof(s.data)](cast[int](s.data) + offset*l.nSitesOuter*s.stride)
  f.elemSize = sizeOf(T)

proc newFieldArray*[V:static[int],T](l:Layout[V]; t:typedesc[Field[V,T]]; n: int):FieldArray[V,T] {.noinit.} =
  result.shape = @[n]
  result.arr = newseq[t](n)
  var s:typeof(result.arr[0].s)
  s.new(l.nSitesOuter*n)
  for i in 0..<n:
    newFarrElem(result.arr[i], l, s, i)

template newFieldArray2*[V:static[int],T](l:Layout[V]; ty:typedesc[Field[V,T]];
    ns: array[2,int]; constraint: untyped):untyped =
  let n = ns[0]
  let m = ns[1]
  var r {.noinit.} :FieldArray[V,T]
  r.shape = @[n,m]
  r.arr = newseq[ty](n*m)
  var t = 0
  for i in 0..<n:
    for j in 0..<m:
      let mu {.inject.} = i
      let nu {.inject.} = j
      if constraint: inc t
  var s:typeof(r.arr[0].s)
  s.new(l.nSitesOuter*t)
  var k = 0
  t = 0
  for i in 0..<n:
    for j in 0..<m:
      let mu {.inject.} = i
      let nu {.inject.} = j
      if constraint:
        newFarrElem(r.arr[k], l, s, t)
        inc t
      inc k
  r

proc newOneOf*[V:static[int],T](fa:FieldArray[V,T]):FieldArray[V,T] {.noinit.} =
  result.shape = fa.shape
  result.arr = newseq[Field[V,T]](fa.arr.len)
  var s:typeof(fa.arr[0].s)
  let (l,n) = block:
    var n = -1
    for i in 0..<fa.arr.len:
      if fa.arr[i] != nil:
        n = i
        break
    if n < 0: return
    (fa.arr[n].l, fa.arr[n].s.len)
  s.new(n)
  var t = 0
  for i in 0..<fa.arr.len:
    if fa.arr[i] == nil: continue
    newFarrElem(result.arr[i], l, s, t)
    inc t

template dataPtr*[V:static[int],T](x: Field[V,T]): auto = x.s.data
template isWrapper*(x: SomeField): bool = false
template isWrapper*(x: typedesc[SomeField]): bool = false
template getT[V:static[int],T](x: Field[V,T]): typedesc = T
template getT[V:static[int],T](x: typedesc[Field[V,T]]): typedesc = T
template has*[F:Field](x: typedesc[F], y: typedesc): bool =
  mixin has, isWrapper
  #static: echo $F.T.type
  when y is Field: true
  else:
    when isWrapper(getT F):
      has(getT F, y)
    else: false

#template `[]`*[F:Field](x:typedesc[F]; i:int):typedesc = F.T
template index*[F:Field](x:typedesc[F]; i:typedesc[int]):typedesc = F.T
template `[]`*(x:Field; i:int):untyped = x.s[i]
#template `[]=`*(x:Field; i:int; y:typed) =
proc `[]=`*(x:Field; i:int; y:auto) =
  x.s[i] := y
template `[]`*(x:Subsetted; i:int):untyped = x.field[i]
template l*(x:Subsetted):untyped = x.field.l
template `[]`*(x:SomeField; st:string):untyped =
  Subsetted[type(x),type(st)](field:x,subset:st)
template `[]`*(x:SomeField; st:Subset):untyped =
  Subsetted[type(x),type(st)](field:x,subset:st)

template `[]`*(x: FieldUnop; i: int): untyped =
  when x.Op == foAdj: adj(x.f1[i])
  elif x.Op == foToSingle: toSingle(x.f1[i])
  elif x.Op == foToDouble: toDouble(x.f1[i])
  else: {.error.}

template `[]`*(x: Field; c: openarray): untyped =
  let ri = x.l.rankIndex(c)
  x{ri.index}
  #let t = eval(x{ri.index})
  #t

template `[]`*(x:FieldArray, i:int):untyped = x.arr[i]
template `[]`*(x:FieldArray, i,j:int):untyped = x.arr[i*x.shape[1]+j]

template even*(x:Field):untyped = x["even"]
template odd*(x:Field):untyped = x["odd"]
template all*(x:Field):untyped = x["all"]
template `even=`*(x:Field; y:auto):untyped = assign(x["even"], y)
template `odd=`*(x:Field; y:auto):untyped =
  mixin assign
  assign(x["odd"], y)
template `all=`*(x:Field; y:auto):untyped = assign(x["all"], y)
template indexField*(x:notSomeField; y:int):untyped = x
template indexField*(x:Field; y:int):untyped = x[y]
template indexField*(x:Subsetted; y:int):untyped = x.field[y]
template indexField*(x:Adjointed[SomeField]; y:int):untyped = x[][y].adj
template indexField*(x:FieldUnop; y:int):untyped = x[y]
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
    mixin indexField
    indexField(`x`.field[0],`y`)
  if sx[0][1].intVal == -1:
    result = quote do:
      -`result`
  for i in 1..<sx.len:
    #t = p[i][1]
    if sx[i][1].intVal == 1:
      template addResult(r,a,b,c: untyped): untyped =
        mixin indexField
        r + indexField(a.field[b],c)
      #result = quote do:
      #  `result` + indexField(`x`.field[`i`],`y`)
      result = getAst(addResult(result,x,i,y))
    else:
      result = quote do:
        `result` - indexField(`x`.field[`i`],`y`)
  #echo result.repr
template indexField*(x:FieldAddSub, y:int):untyped = indexFieldM(x, x.S, y)
template `[]`*(x:FieldAddSub, y:int):untyped = indexField(x, y)
template l*(x:FieldAddSub):untyped = x.field[0].l
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
template indexField*[T](x: FieldMul[T], y: int): untyped =
  indexFieldM(x, x.T, y)
template `[]`*(x:FieldMul, y:int):untyped = indexField(x, y)

template itemsI*(n0,n1:int):untyped =
  let n = n1 - n0
  var ti0 = n0 + ((threadNum*n) div numThreads)
  var ti1 = n0 + (((threadNum+1)*n) div numThreads)
  #echo "ti0: ", ti0, "  ti1: ", ti1
  var i = ti0
  while i < ti1:
    yield i
    inc(i)
template itemsI*(n0,n1,b:int):untyped =
  ## The same as itemsI(n0,n1), with the extra parameter `b` specifying
  ## the minimum block size which will not be separated by threads.
  let n = n1 - n0
  let nb = (n div b) + int(n mod b > 0)
  var ti0 = n0 + b*((threadNum*nb) div numThreads)
  var ti1 = n0 + b*(((threadNum+1)*nb) div numThreads)
  if ti1 > n1: ti1 = n1
  var i = ti0
  while i < ti1:
    yield i
    inc(i)
template itemsI2*(n0,n1:int):untyped =
  #let n = n1 - n0
  #var ti0 = n0 + ((threadNum*n) div numThreads)
  #var ti1 = n0 + (((threadNum+1)*n) div numThreads)
  #echo "ti0: ", ti0, "  ti1: ", ti1
  let s = 64
  var i = n0 + s*threadNum
  var j = s
  while i < n1:
    yield i
    inc i
    dec j
    if j==0:
      i += s*(numThreads-1)
      j = s
iterator items*(l:Layout):int {.inline.} =
  let n = l.nSitesOuter
  itemsI(0, n)
iterator sites*(l:Layout):int {.inline.} =
  let n = l.nSites
  itemsI(0, n, VLEN)
iterator sites*(f:Field):int {.inline.} =
  let n = f.l.nSites
  itemsI(0, n, VLEN)
iterator items*(s:Subset):int {.inline.} =
  let n0 = s.lowOuter
  let n1 = s.highOuter
  itemsI(n0, n1)
iterator sites*(s:Subset):int {.inline.} =
  let n0 = s.low
  let n1 = s.high
  itemsI(n0, n1, VLEN)
#iterator all*(x:Field):int {.inline.} =
#  let n = x.l.nSitesOuter
#  itemsI(0, n)
iterator items*(x:Field):int {.inline.} =
  #let n = x.l.nSitesOuter
  let l = x.l
  let n = l.nSitesOuter
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
iterator sites*(x: Subsetted): int {.inline.} =
  when x.subset is string:
    let s = getSubset(x.field.l, x.subset)
  else:
    let s = x.subset
  let n0 = s.low
  let n1 = s.high
  #echo "n0: ", n0, " n1: ", n1
  itemsI(n0, n1, VLEN)
iterator items*(x:FieldAddSub):int {.inline.} =
  let n = x.field[0].l.nSitesOuter
  itemsI(0, n)
iterator items*(x:FieldMul):int {.inline.} =
  let n = x.field[0].l.nSitesOuter
  itemsI(0, n)
#macro filter*(x:SomeField; s:string; pred:untyped):Mask {.inline.} =
#iterator filter*(x:SomeField; s:string; pred:untyped):Mask {.inline.} =
  ###

#import types
#export types

template fmask*(f: Field; i: int): untyped =
  when f.l.V == 1:
    when has(type(f),Simd):
      f[i][asSimd(0)]
    else:
      f[i]
  else:
    #when true:
    when false:
      mixin varMasked
      let e = i div f.l.V
      let l = i mod f.l.V
      let mask = 1 shl l
      #let fe = f[e]  # workaround for Nim codegen bug
      varMasked(f[e], mask)
    else:
      let e = i div f.l.V
      let l = i mod f.l.V
      #indexed(f[e], l)
      #static: echo "fmask: type f[e]: ", $type(f[e])
      f[e][asSimd(l)]

template `{}`*(f: Field; i: int): untyped =
  fmask(f, i)

template `{}`*(f: Subsetted; i: int): untyped =
  fmask(f.field, i)

#template mindex*(f: Field; i: int): untyped =
#  fmask(f, i)

#[
template `{}=`*(f: Field; i: int, y: typed): untyped =
    let e = i div f.l.V
    let l = i mod f.l.V
    f[e][Simd[l]] = y
]#

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
#[
proc applyOp1x(x,y:NimNode; op:string):auto =
  let o = ident(op)
  result = quote do:
    let tx = `x`
    let ty = `y`
    #echoImm `x`[0] is VMconcept1
    #echoImm t isnot VMconcept2
    for e in tx:
      mixin `o`
      #mixin isMatrix
      #echoAll isMatrix(`x`[e])
      `o`(tx[e], ty)
template applyOp1Impl(x,y,o:untyped) =
  let tx = `x`
  let ty = `y`
  #echoImm `x`[0] is VMconcept1
  #echoImm t isnot VMconcept2
  for e in tx:
    mixin `o`
    #mixin isMatrix
    #echoAll isMatrix(`x`[e])
    `o`(tx[e], ty)
proc applyOp1(x,y:NimNode; op:string):auto =
  let o = ident(op)
  result = getAst(applyOp1Impl(x,y,o))
]#

#[
var exprInstInfo {.compiletime.}: type(instantiationInfo())
macro debugExpr(body: typed): untyped =
  #let ii = instantiationInfo(1)
  #let ii = lineInfoObj(body)
  #let ii = lineInfo(body)
  let ii = $exprInstInfo
  let br = body.repr
  result = newStmtList()
  result.add quote do:
    {.emit: ["\n/* debugExpr\n", `ii`, "\n", `br`, "\n*/\n"] .}
  result.add body
  echo ii
  echo br
]#

#proc applyOp2(x,y:NimNode; ty:typedesc; op:string):auto =
proc applyOp2(x,y:NimNode; op:string):auto =
  #echo ty.getType.treeRepr
  #echo ty.getType.getImpl.treeRepr
  let o = ident(op)
  result = quote do:
    #debugExpr:
      let xx = `x`
      let yy = `y`
      for e in xx:
        when noAlias:
          staticTraceBegin: `o Field2`
          type Fpx = object
            v: type(xx[e])
          var xp = cast[ptr carray[Fpx]](xx[0].addr)
          `o`(xp[][e].v, indexField(yy, e))
          staticTraceEnd: `o Field2`
        else:
          staticTraceBegin: `o Field2`
          `o`(xx[e], indexField(yy, e))
          staticTraceEnd: `o Field2`
  #echo result.treerepr
#macro id(op:static string):auto =
#  result = newIdentNode(op)
#macro id(op:string):auto =
#  result = newIdentNode(op.strVal)
#template makeOps2(o,f,fM: untyped) {.dirty.} =
#  proc f*(x:Subsetted; y:notSomeField2) =
#    for e in x:
#      mixin o
#      o(tx[e], ty)
template makeOps(op,f,fM,s: untyped) {.dirty.} =
  #makeOps2(id(s),f,fM)
  #macro f*(x:Subsetted; y:notSomeField2):auto = applyOp1(x,y,s)
  proc f*(x:Subsetted; y:notSomeField2) =
    for e in x:
      mixin f
      f(x[e], y)
  macro f*(x:Subsetted; y:SomeField2):auto = applyOp2(x,y,s)
  #macro fM*(x:Field; y:notSomeField; ty:typedesc):auto = applyOp1(x,y,s)
  proc fM*(x:Field; y:notSomeField) =
    for e in x:
      mixin f
      f(x[e], y)
  macro fM*(x:Field; y:SomeField):auto = applyOp2(x,y,s)
  template f*(x:Field; y:auto):untyped =
    #when declaredInScope(subsetObject):
    when declared(subsetObject):
      #echo "subsetObj" & s
      f(x[subsetObject], y)
    elif declared(subsetString):
      #echo "subsetString" & s
      f(x[subsetString], y)
    else:
      #fM(x, y, y.type)
      staticTraceBegin: `f FieldAuto`
      fM(x, y)
      staticTraceEnd: `f FieldAuto`
  when profileEqns:
    template op*(x:Field; y:auto):untyped =
      #static: exprInstInfo = instantiationInfo(-1)
      block:
        tic(-2)
        f(x, y)
        toc(asttostr(op), -2)
    template op*(x:Subsetted; y:auto):untyped =
      #static: exprInstInfo = instantiationInfo(-1)
      block:
        tic(-2)
        f(x, y)
        toc(asttostr(op), -2)
  else:
    template op*(x:Field; y:auto):untyped =
      #static: exprInstInfo = instantiationInfo(1)
      f(x, y)
    #template op*(x:var Field; y:auto):untyped = f(x, y)
    template op*(x:Subsetted; y:auto):untyped =
      #static: exprInstInfo = instantiationInfo(1)
      f(x, y)
makeOps(`:=`, assign, assignM, "assign")
makeOps(`+=`, iadd, iaddM, "iadd")
makeOps(`-=`, isub, isubM, "isub")
makeOps(`*=`, imul, imulM, "imul")

proc mul*(r:Field; x:Field2; y:Field3) =
  mixin mul
  for e in r:
    mul(r[e], x[e], y[e])

proc norm2P*(f:SomeField):auto =
  tic()
  mixin norm2, inorm2, simdSum, items, toDouble
  #var n2:type(norm2(f[0]))
  var n2: evalType(norm2(toDouble(f[0])))
  #echo n2
  #let t = f
  for x in items(f):
    inorm2(n2, toDouble(f[x]))
  toc("norm2 local")
  #echoAll n2
  result = simdSum(n2)
  toc("norm2 simd sum")
  #echoAll myRank, ",", threadNum, ": ", result
  #threadSum(result)
  #toc("norm2 thread sum")
  #rankSum(result)
  #toc("norm2 rank sum")
  f.l.threadRankSum(result)
  #echo result
  toc("norm2 thread rank sum")
template norm2*(f:SomeAllField):auto =
  when declared(subsetObject):
    #echo "subsetObj" & s
    norm2P(f[subsetObject])
  elif declared(subsetString):
    #echo "subset norm2"
    norm2P(f[subsetString])
  else:
    norm2P(f)
template norm2*(f:Subsetted):auto = norm2P(f)

proc norm2subtract*(x: Field, y: float): float =
  var s: evalType(norm2(toDouble(x[0])))
  for i in x:
    s += x[i].toDouble.norm2 - y
  result = s.simdReduce
  x.l.threadRankSum(result)

proc norm2diffP*(f,g:SomeField):auto =
  tic()
  mixin norm2, inorm2, simdSum, items, toDouble
  #var n2:type(norm2(f[0]))
  var n2: evalType(norm2(toDouble(f[0])))
  #echo n2
  #let t = f
  for x in items(f):
    let t = toDouble(f[x]) - toDouble(g[x])
    inorm2(n2, t)
  toc("norm2 local")
  #echoAll n2
  result = simdSum(n2)
  toc("norm2 simd sum")
  #echoAll myRank, ",", threadNum, ": ", result
  #threadSum(result)
  #toc("norm2 thread sum")
  #rankSum(result)
  #toc("norm2 rank sum")
  f.l.threadRankSum(result)
  #echo result
  toc("norm2 thread rank sum")
template norm2diff*(f,g:SomeAllField):auto =
  when declared(subsetObject):
    #echo "subsetObj" & s
    norm2diffP(f[subsetObject], g[subsetObject])
  elif declared(subsetString):
    #echo "subset norm2"
    norm2diffP(f[subsetString], g[subsetString])
  else:
    norm2diffP(f, g)
template norm2diff*(f,g:Subsetted):auto = norm2diffP(f,g)

proc dotP*(f1:SomeField; f2:SomeField2):auto =
  tic()
  mixin dot, idot, simdSum, items, toDouble, eval
  #var d:type(dot(f1[0],f2[0]))
  var d: evalType(toDouble(dot(f1[0],f2[0])))
  let t1 = f1
  let t2 = f2
  for x in items(t1):
    #idot(d, t1[x], t2[x])
    d += dot(t1[x], t2[x])
  toc("dot local")
  result = simdSum(d)
  toc("dot simd sum")
  #threadSum(result)
  #rankSum(result)
  f1.l.threadRankSum(result)
  toc("dot thread rank sum")
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
  var d: evalType(toDouble(redot(f1[0],f2[0])))
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
  f1.l.threadRankSum(result)
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
  var tr: evalType(trace(m[0]))
  for x in m:
    tr += m[x].trace
  #echo tr
  result = simdSum(tr)
  #threadSum(result)
  #rankSum(result)
  m.l.threadRankSum(result)

proc sumP*(f:SomeField):auto =
  mixin inc, simdSum, items
  var s: evalType(f[0])
  let t = f
  for x in items(t):
    iadd(s, t[x])
  result = simdSum(s)
  #threadSum(result)
  #rankSum(result)
  f.l.threadRankSum(result)
template sum*(f:SomeAllField):untyped =
  when declared(subsetString):
    sumP(f[subsetString])
  else:
    sumP(f)
template sum*(f:Subsetted):untyped = sumP(f)

template `-`*(x:SomeField):untyped = fieldAddSub(-1,x)
template `+`*(x:SomeField,y:SomeField2):untyped = fieldAddSub(1,x,1,y)
template `+`*(x:SomeField,y:notSomeField2):untyped = fieldAddSub(1,x,1,y)
template `+`*(x:notSomeField,y:SomeField2):untyped = fieldAddSub(1,x,1,y)
template `-`*(x:SomeField,y:SomeField2):untyped = fieldAddSub(1,x,-1,y)
template `-`*(x:SomeField,y:notSomeField2):untyped = fieldAddSub(1,x,-1,y)
template `-`*(x:notSomeField,y:SomeField2):untyped = fieldAddSub(1,x,-1,y)
template `*`*(x:SomeField,y:SomeField2):untyped = fieldMul(x,y)
template `*`*(x:SomeField,y:notSomeField2):untyped = fieldMul(x,y)
template `*`*(x:notSomeField,y:SomeField2):untyped = fieldMul(x,y)
template `/`*(x:SomeField,y:SomeNumber):untyped =
  let t = 1.0/y
  fieldMul(x,t)
#template makeBinop(f:untyped,s:string):untyped =
#  template `f`*(x:SomeField; y:SomeField2):untyped =
#    echoImm x
#    echoImm y
#    fieldBinop(`s`, x, y)
#makeBinop(`+`,"+")
#makeBinop(`-`,"-")
#makeBinop(`*`,"*")
#makeBinop(`/`,"/")
#template `*`*(x:notSomeField; y:SomeField2):untyped =
#  fieldBinop("*", x, y)

#template toSingle*(x: SomeField): untyped =
#  fieldUnop(foToSingle, x)

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

#import layout/shifts

#template assign*(x:var Field; y:Shifted):untyped =
#  shift(x, y.dir, y.len, y.field)

when isMainModule:
  import qex
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
