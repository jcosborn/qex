import macros
export macros
import metaUtils
import strutils

template has*[T](x: T, y: typedesc): bool =
  mixin has
  has(type T, y)

template numberType*[T](x: T): untyped =
  mixin numberType
  numberType(type T)

template simdLength*[T](x: T): untyped =
  mixin numberType
  simdLength(type T)

type
  WrapperFlags* = enum
    wfDeref,wfEval,wfIndex,  # set if type will provide own implementation
    wfStorage,        # wrapper types
    wfUnwrap,wfWrapApply,wfDollar, wfGet, wfTest,
    wfPtr
  WFSet* = set[WrapperFlags]
const wfNone* = WFSet({})
type
  Wrapper* = concept C
    C.flags is WFSet
  WrapperNone* = concept C
    C.flags.card == 0
  WrapperNotNone* = concept C
    C.flags.card != 0
  WrapperDeref* = concept C
    wfDeref in C.flags
  WrapperNotDeref* = concept C
    wfDeref notin C.flags
  WrapperStorage* = concept C
    wfStorage in C.flags
  WrapperNotStorage* = concept C
    wfStorage notin C.flags

#[
macro makeCall(f: string, a: varargs[untyped]): untyped =
  result = newCall(ident(f.strVal))
  for i in 0..<a.len:
    result.add a[i]
]#

template makeFieldGetter(t,f,f0: untyped) {.dirty.} =
  template f*(x: t): untyped =
    #static: echo "fMulMV0"
    #echoKind: x
    x.f0
  macro f*(x: t{nkObjConstr}): untyped =
    #echo x.repr
    for i in countdown(x.len-1, 1):
      #echo x[i][0].repr
      #echo astToStr(f0)
      if x[i][0].repr == astToStr(f0):
        result = x[i][1]
        break
      if i==1:
        echo "error: field name not found: ", astToStr(f0)
        quit -1
  macro f*(x: t{nkStmtListExpr}): untyped =
    result = newNimNode(nnkStmtListExpr)
    for i in 0..(x.len-2):
      result.add x[i]
    result.add newCall(ident(astToStr(f)), x[^1])

#[
macro makeDerefOp(t: untyped): untyped =
  newCall(bindSym("makeFieldGetter"), t, ident("[]"), ident("f"&t.repr))
]#

template makeDeref*(t,u:untyped):untyped {.dirty.} =
  #[
  bind ctrace
  template `[]`*(x:t):untyped =
    when compiles(unsafeAddr(x)):
       #ctrace()
      cast[ptr u](unsafeAddr(x))[]
    else:
      #ctrace()
      (u)x
  template `[]=`*(x:t; y:auto):untyped =
    when compiles(unsafeAddr(x)):
      cast[ptr u](unsafeAddr(x))[] = y
    else:
      (u)x = y
  ]#
  #template callBracket*(x: t): untyped = x.v
  #macro callBracket*(x: t{nkObjConstr}): untyped =
  #  #echo x.treerepr
  #  result = x[1][1]
  #template `[]`*(x: t): untyped = callBracket(normalizeAst(x))
  #template `[]`*(x: t): untyped = callBracket(x)
  template `[]=`*(x:t; y:auto):untyped =
    x.v = y
  template `[]`*(x:t):untyped = x.v
  #template `[]=`*(x:t; y:auto):untyped =
  #  x.v[] = y

#[
template makeWrapper2(t,s,u: untyped): untyped {.dirty.} =
  #type t*[T] = distinct T
  type t*[T] = object
    v*: T

  template s*(xx: typed): untyped =
    let u = xx
    t[type(u)](v: u)

  makeDeref(t, 0)

#template makeWrapper*(t,s: untyped): untyped {.dirty.} =
#  bind makeWrapper2, x_id
#  makeWrapper2(t,s,x_id(s))

macro makeWrapper*(t,s: untyped): untyped =
  var xid = s
  case s.kind
  of nnkIdent:
    xid = newIdentNode("t_" & $s)
  of nnkAccQuoted:
    var a = "t_"
    for c in s:
      a &= $c
    #echo s.treerepr
    #echo a
    xid = newIdentNode(a)
  else:
    echo s.treerepr
  result = getAst(makeWrapper2(t,s,xid))
]#

template getPtr*[T](x: T): untyped =
  when T is RefWrap:
    x
  elif T is ptr or T is ref:
    x
  elif compiles(unsafeAddr x):
    unsafeAddr x
  else:
    var tGetPtr = x
    addr tGetPtr

template getPtr*[T](x: typedesc[T]): typedesc =
  when T is RefWrap:
    T
  elif T is ptr or T is ref:
    T
  else:
    ptr T

template derefPtr*[T](x: T): untyped =
  when T is ptr:
    x[]
  else:
    x

template derefPtr*[T](x: typedesc[T]): typedesc = T
template derefPtr*[T](x: typedesc[ptr T]): typedesc = T

#template makeWrapperTypeX(wf:WFSet, name,fName,asName,tasName: untyped) =
template makeWrapperTypeX(wf:typed, name,fName,asName,tasName: untyped) =
  #when wfPtr in wf:
  when false:
    type
      name*[T] {.byref.} = object # ## wrapper type
        fName*: ptr T
  else:
    type
      name*[T] {.byref.} = object # ## wrapper type
        fName*: T
  #template tasName*(x: typed): untyped =
  #  name[type(x)](fName: x)
  #proc tasName*[T](x: T): name[T] {.inline,noInit.} =
  #  result.fName = x
  #when wfPtr in wf:
  when true:
    template `asName Obj`*[T](x: T): untyped =
      name[type T](fName: x)
    template `asName View`*[T](x: T): untyped =
      `asName Obj`(getPtr x)
    template asName*[T](x: T): untyped =
      `asName Obj`(x)
  else:
    template asName*[T](x: T): untyped =
      name[type T](fName: x)
    # ## wrap an object, x, as a $NAME type
    #lets(x,xx):
    #static: echo "asColor typed"
    #dumpTree: xx
    #name[type(x)](fName: x)
    #let tasName = x
    #name[type(tasName)](fName: tasName)
    #Color[type(x_asColor)](x_asColor)
    #Color(x_asColor)
  #template asName*[T](x: typedesc[T]): untyped =
  #  name[type T]
  template `asName Obj`*[T](x: typedesc[T]): typedesc =
    name[typeof T]
  template `asName View`*[T](x: typedesc[T]): typedesc =
    `asName Obj`(getPtr typeof T)
  template asName*[T](x: typedesc[T]): typedesc =
    `asName Obj`(typeof T)
  #  flattenCallArgs(tasName, x)
  #proc asName*[T](x: T): name[T] {.inline,noInit.} =
  #  result.fName = x
  # ## dereference a $NAME object
  #template `[]`*(x: name): untyped =
  #  x.fName
  #makeDerefOp(name)
  #template derefXX*(x: name): untyped =
  #  x.fName
  makeFieldGetter(name, derefXX, fName)
  when wfDeref notin wf:
    template `[]`*[T](x: typedesc[name[T]]): typedesc =
      derefPtr(T)
    #when wfPtr in wf:
    when true:
      template `[]`*(x: name): untyped =
        #static: echo "wrapper []"
        #debugType: x
        derefPtr flattenCallArgs(derefXX, x)
        #derefPtr derefXX(x)
    else:
      template `[]`*(x: name): untyped =
        #static: echo "wrapper []"
        #debugType: x
        flattenCallArgs(derefXX, x)
  template isWrapper*(x: name): untyped = true
  template isWrapper*(x: typedesc[name]): untyped = true
  template asWrapper*(x: name, y: typed): untyped =
    #static: echo "asWrapper Color"
    #dumpTree: y
    asName(y)
  template asWrapper*(x: typedesc[name], y: typed): untyped =
    asName(y)

proc makeWrapperTypeP*(name: NimNode; docs: string, wf: WFSet): NimNode =
  let Name = capitalizeAscii(name.repr)
  let aName = if Name[0..1]=="As": "a"&Name[1..^1]
              else: "as"&Name
  let fName = ident("f" & Name)
  let asName = ident(aName)
  let tasName = ident("t_" & aName)
  let w = newLit(wf)
  result = getAst(makeWrapperTypeX(w,name,fName,asName,tasName))
  #result = result.replaceComments(("$DOCS",docs),("$NAME",Name))
  #echo result.repr

macro makeWrapperType*(name,docs: untyped): untyped =
  #echo $docs[0]
  #echo name.treerepr
  var d: string
  when docs.type is string: d = docs
  else: d = $docs[0]
  makeWrapperTypeP(name, d, wfNone)

macro makeWrapperType*(name: untyped): untyped =
  let d = "wrapper type for " & $name & " objects"
  makeWrapperTypeP(name, d, wfNone)

macro makeWrapperF*(wf: static[WFSet], name: untyped): auto =
  let d = "wrapper type for " & $name & " objects"
  makeWrapperTypeP(name, d, wf)
