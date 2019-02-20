import base/metaUtils

echo "* Proc return an auto generic type"
proc rg[T](x:T):auto = x
type
  M[N:static[int],T] = object
    d:array[N,T]
proc rgt =
  var x,y {.noinit.}:M[3,float]
  for i in 0..<x.N: x.d[i] = 0.1+i.float
  inlineProcs:
    y = x.rg
  for i in 0..<y.N: echo i," ",y.d[i]
rgt()
