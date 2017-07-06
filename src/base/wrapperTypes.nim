import stdUtils
#import metaUtils

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
  template callBracket*(x: t): untyped = x.v
  macro callBracket*(x: t{nkObjConstr}): untyped =
    #echo x.treerepr
    result = x[1][1]
  template `[]`*(x: t): untyped = callBracket(normalizeAst(x))
  template `[]=`*(x:t; y:any):untyped =
    x.v = y
  #template `[]`*(x:t):untyped = x.v[]
  #template `[]=`*(x:t; y:any):untyped =
  #  x.v[] = y

template makeWrapper*(t,s: untyped): untyped {.dirty.} =
  #type t*[T] = distinct T
  type t*[T] = object
    v*: T
  #type t*[T] = object
  #  v*:ptr T
  template s*(xx: typed): untyped =
    t[type(xx)](v: xx)
  #proc s*(xx:any):auto {.inline.} =
  #  lets(x,xx):
  #    when compiles(addr(x)):
  #    #when compiles(unsafeAddr(x)):
  #      #ctrace()
  #      #cast[ptr t[type(x)]](addr(x))[]
  #      cast[ptr t[type(x)]](unsafeAddr(x))[]
  #      #cast[t[type(x)]](x)
  #    else:
  #      #dumptree(x)
  #      #ctrace()
  #      #(t[type(x)])x
  #      cast[t[type(x)]](x)
  #      #var y = x
  #      #cast[t[type(x)]](addr(y))
  #makeDeref(t, x.T)
  makeDeref(t, 0)
