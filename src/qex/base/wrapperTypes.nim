import macros
export macros
import stdUtils
import metaUtils
import strutils

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
  template `[]=`*(x:t; y:any):untyped =
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
  template `[]=`*(x:t; y:any):untyped =
    x.v = y
  template `[]`*(x:t):untyped = x.v
  #template `[]=`*(x:t; y:any):untyped =
  #  x.v[] = y

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
    xid = newIdentNode(!("t_" & $s))
  of nnkAccQuoted:
    var a = "t_"
    for c in s:
      a &= $c
    #echo s.treerepr
    #echo a
    xid = newIdentNode(!a)
  else:
    echo s.treerepr
  result = getAst(makeWrapper2(t,s,xid))

template makeWrapperTypeX(name,fName,asName,tasName: untyped) =
  type
    name*[T] = object # ## wrapper type
      fName*: T
  template asName*(x: typed): untyped =
    # ## wrap an object, x, as a $NAME type
    #lets(x,xx):
    #static: echo "asColor typed"
    #dumpTree: xx
    let tasName = x
    name[type(tasName)](fName: tasName)
    #Color[type(x_asColor)](x_asColor)
    #Color(x_asColor)
  # ## dereference a $NAME object
  template `[]`*(x: name): untyped =
    x.fName
  template isWrapper*(x: name): untyped = true
  template asWrapper*(x: name, y: typed): untyped =
    #static: echo "asWrapper Color"
    #dumpTree: y
    asName(y)

proc makeWrapperTypeP*(name: NimNode; docs: string): NimNode =
  let Name = capitalizeAscii($name)
  let fName = ident("f" & $Name)
  let asName = ident("as" & $Name)
  let tasName = ident("t_as" & $Name)
  result = getAst(makeWrapperTypeX(name,fName,asName,tasName))
  #result = result.replaceComments(("$DOCS",docs),("$NAME",Name))
  #echo result.repr

macro makeWrapperType*(name,docs: untyped): untyped =
  #echo $docs[0]
  #echo name.treerepr
  var d: string
  when docs.type is string: d = docs
  else: d = $docs[0]
  makeWrapperTypeP(name, d)

macro makeWrapperType*(name: untyped): untyped =
  let d = "wrapper type for " & $name & " objects"
  makeWrapperTypeP(name, d)
