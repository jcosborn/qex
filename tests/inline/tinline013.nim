import base/metaUtils

echo "* static[T]"
proc fs(x:int, y:static[int]):int = x*y
inlineProcs:
  var x = 2
  let y = x.fs 3
  echo y
