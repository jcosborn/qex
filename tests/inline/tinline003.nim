import base/metaUtils

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

echo "* avoid duplicate computations"
proc inplace(x:var float, y:float) =
  x = x + y
  x = x + 1000*y
inlineProcs:
  var x {.noinit.}: T
  var z {.noinit.}: T
  for i in 0..<x.len: x[i] = i.float
  z.loop(x.loop2 x)
  for i in 0..<x.len:
    z[i].inplace(0.1*i.float)
    echo i," ",z[i]
