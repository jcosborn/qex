import base/metaUtils

proc f1(r: var any; x: any) = r = 2*x
proc f2(x: any): auto = 2*x

proc a1(x: float) =
  inlineProcs:
    var r: float
    var s: type(r)
    f1(r, x)
proc a2(x: float) =
  inlineProcs:
    var r = f2(x)

echo "* Basics"
a1(1.0)
a2(1.0)
