import base/metaUtils

type
  M[N:static[int],T] = object
    d:array[N,T]
echo "* Object wrappers of generic types"
type
  W[T] = object
    o:T
  Walt[T] = object
    o:W[T]
proc toAlt[S](x:W[S]):auto = Walt[S](o:x)
proc toAlt2[S](x:W[W[S]]):auto = Walt[S](o:x.o)
block:
  var A {.noinit.} :M[3,float]
  for i in 0..<A.d.len: A.d[i] = i.float
  var w = W[type(A)](o:A)
  inlineProcs:
    var walt = w.toAlt
  for i in 0..<walt.o.o.d.len: echo walt.o.o.d[i]
  var w2 = W[type(w)](o:w)
  inlineProcs:
    var walt2 = w2.toAlt2
  for i in 0..<walt2.o.o.d.len: echo walt2.o.o.d[i]
