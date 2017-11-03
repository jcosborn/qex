import macros

proc volume*[T](x: openarray[T]): auto {.inline,noInit.} =
  result = x[0]
  for i in 1..<x.len: result *= x[i]

proc lexIndex*[T1,T2](coord: openarray[T1],
                      dims: openarray[T2]): auto {.inline,noInit.} =
  result = coord[^1]
  for i in 2..coord.len:
    result = result * dims[^i] + coord[^i]

proc lexIndexR*[T1,T2](coord: openarray[T1],
                       dims: openarray[T2]): auto {.inline,noInit.} =
  result = coord[0]
  for i in 1..<coord.len:
    result = result * dims[i] + coord[i]

proc lexCoord*[T1,T2](coord: var openarray[T1], l: SomeNumber,
                      dims: openarray[T2]) {.inline.} =
  var k = l
  coord[0] = k mod dims[0]
  for i in 1..<coord.len:
    k = k div dims[i-1]
    coord[i] = k mod dims[i]

proc lexCoordR*[T1,T2](coord: var openarray[T1], l: SomeNumber,
                       dims: openarray[T2]) {.inline.} =
  var k = l
  coord[^1] = k mod dims[^1]
  for i in 2..coord.len:
    k = k div dims[^(i-1)]
    coord[^i] = k mod dims[^i]

macro eqSum*[T](x: var openarray[T], y: varargs[untyped]): auto =
  #echo x.treerepr
  #echo y.treerepr
  template sumX(i,r,v: untyped) =
    for i in r.low..r.high:
      r[i] = v
  var i = ident"i"
  var v = newCall(ident"[]", y[0], i)
  for k in 1..<y.len:
    v = newCall(ident"+", v, newCall(ident"[]", y[k], i))
  result = getAst(sumX(i, x, v))
  #echo result.treerepr
  #echo result.repr

proc coordMod*(x: SomeNumber, y: any): auto {.inline,noInit.} =
  result = (x + abs(x)*y) mod y

proc eqMod*[T1,T2,T3](r: var openarray[T1], x: openarray[T2],
                      dims: openarray[T3]) {.inline,noInit.} =
  for i in r.low..r.high:
    r[i] = coordMod(x[i], dims[i])

proc `$`*[N,T](x: array[N,T]): string =
  result = "[" & $x[0]
  for i in 1..<x.len:
    result &= "," & $x[i]
  result &= "]"

when isMainModule:
  var dims = [10,11,12,13]
  var x = [1,2,3,4]
  var y = [2,3,4,5]

  echo x
  let ix = lexIndex(x, dims)
  x.lexCoord(ix, dims)
  echo x

  echo y
  let iy = lexIndexR(y, dims)
  y.lexCoordR(iy, dims)
  echo y

  x.eqSum(y,x)
  echo x

  x.eqSum(x,dims)
  echo x

  y.eqMod(x,dims)
  echo y

  x.eqSum(x,[-100,0,0,0])
  echo x

  y.eqMod(x,dims)
  echo y
