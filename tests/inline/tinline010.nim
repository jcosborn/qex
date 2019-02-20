import base/metaUtils

echo "* Proc with local proc/template"
type Mt[F] = object
  m:array[3,F]
proc len[F](m:Mt[F]):int = m.m.len
template `[]`[F](x:Mt[F],i:int):F = x.m[i]
iterator items[F](m:Mt[F]):F =
  var i = 0
  while i < m.len:
    yield m[i]
    inc i
proc lp =
  proc `$`[F](m:Mt[F]):string =
    result = "Mt["
    for x in m: result &= " " & $x
    result &= " ]"
  template go[F](x:Mt[F],y:untyped) =
    for i in 0..<x.len: x[i] += y[i]
  var x = Mt[float](m:[1.0,2.0,3.0])
  var y = [0.1,0.2,0.3]
  for i in 0..<y.len: y[i] *= 0.1
  x.go y
  echo x
inlineProcs:
  lp()
