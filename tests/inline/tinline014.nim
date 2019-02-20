import base/metaUtils

type T = array[3,float]

echo "* avoid duplicate computations"
proc inplace(x:var float, y:float) =
  x = x + y
  x = x + 1000*y
inlineProcs:
  var x {.noinit.}: T
  for i in 0..<x.len: x[i] = i.float
  var s = 0.0
  for m in mitems(x):
    s += m
    m.inplace(1000 * s)
  for i in 0..<x.len: echo x[i]
