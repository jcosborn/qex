import base/metaUtils

echo "* object construction"
proc oc(x:int):auto =
  type A = object
    x:int
  return A(x:x)
block:
  inlineProcs:
    var x = 3
    echo oc(x).x
