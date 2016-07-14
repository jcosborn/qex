import stdUtils
#import metaUtils

template makeDeref*(t,u:untyped):untyped {.dirty.} =
  #[
  bind ctrace
  template `[]`*(x:t):expr =
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
  template `[]`*(x:t):expr = x.v
  template `[]=`*(x:t; y:any):untyped =
    x.v = y
  #template `[]`*(x:t):expr = x.v[]
  #template `[]=`*(x:t; y:any):untyped =
  #  x.v[] = y

template makeWrapper*(t,s:untyped):untyped =
  #type t*[T] = distinct T
  type t*[T] = object
    v*:T
  #type t*[T] = object
  #  v*:ptr T
  template s*(xx:typed):expr =
  #proc s*(xx:any):auto {.inline.} =
    subst(x,xx):
      when compiles(addr(x)):
      #when compiles(unsafeAddr(x)):
        #ctrace()
        cast[ptr t[type(x)]](addr(x))[]
        #cast[t[type(x)]](x)
      else:
        #dumptree(x)
        #ctrace()
        #(t[type(x)])x
        cast[t[type(x)]](x)
        #var y = x
        #cast[t[type(x)]](addr(y))
  #makeDeref(t, x.T)
  makeDeref(t, 0)
