import base/metaUtils

echo "* multiple iterators"
type T = array[3,float]
proc loop(x:var T, y:T) =
  echo "loop"
  let n = 3.0
  for k in 0..<x.len:
    x[k] = n * y[k]
proc loop2(x:T,y:T):T =
  echo "loop2"
  let n = 0.1
  for k in 0..<x.len:
    result[k] = n * x[k] + y[k]
  for k in 0..<x.len:
    result[k] = n * result[k]
proc loop3(x:var T,y:T) =
  echo "loop3"
  x.loop y
  x = y.loop2 y
proc cl =
  var x {.noinit.}: T
  var z {.noinit.}: T
  for i in 0..<x.len: x[i] = i.float
  inlineProcs: z.loop x
  for i in 0..<x.len: echo z[i]
  inlineProcs: z = x.loop2 x
  for i in 0..<x.len: echo z[i]
  inlineProcs:
    z.loop x
    z = x.loop2 x
  for i in 0..<x.len: echo z[i]
  inlineProcs: z.loop3 x
  for i in 0..<x.len: echo z[i]
cl()
