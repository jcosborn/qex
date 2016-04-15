import stdUtils
#import metaUtils

template makeDeref*(t,u:untyped):untyped {.dirty.} =
  bind ctrace
  template `[]`*(x:t):expr =
    when compiles(addr(x)):
      #ctrace()
      cast[ptr u](addr(x))[]
    else:
      #ctrace()
      (u)x
  template `[]=`*(x:t; y:any):untyped =
    when compiles(addr(x)):
      cast[ptr u](addr(x))[] = y
    else:
      (u)x = y

template makeWrapper*(t,s:untyped):untyped =
  type t*[T] = distinct T
  template s*(x:typed):expr =
    when compiles(addr(x)):
      #ctrace()
      cast[ptr t[type(x)]](addr(x))[]
    else:
      #ctrace()
      (t[type(x)])x
  makeDeref(t, x.T)
