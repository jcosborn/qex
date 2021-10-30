import base/metaUtils

proc f1(r: var auto; x: auto) = r = 2*x
proc f2(x: auto): auto = 2*x

proc a1(x: float) =
  inlineProcs:
    var r: float
    var s {.used.}: type(r)
    f1(r, x)
proc a2(x: float) =
  inlineProcs:
    var r {.used.} = f2(x)

echo "* Basics"
a1(1.0)
a2(1.0)
