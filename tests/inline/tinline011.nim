import base/metaUtils

echo "* varargs"
proc square[T](x:T):float =
  let y = float(x)
  y*y
proc fv(z:var float, xs:varargs[float, square]) =
  for x in xs: z += x
block:
  inlineProcs:
    var
      s = 0.0
      x = 1
      y = 2.2
      z:float32 = 3.3
    s.fv(x,y,z)
    echo s
