import base/metaUtils

echo "* redeclaration of formal params"
proc redecl(x:var float, y:float) =
  block:
    var x = x
    x += y
    var y = y
    y += 1
    #echo x," ",y
  x += 3
  let x = x
  var y = y
  y += x
  #echo x," ",y
block:
  echo "Without inlining:"
  var x = 1.0
  var y = 0.1
  x.redecl(y+0.01)
  echo x," ",y
block:
  echo "With inlining:"
  inlineProcs:
    var x = 1.0
    var y = 0.1
    x.redecl(y+0.01)
  echo x," ",y
