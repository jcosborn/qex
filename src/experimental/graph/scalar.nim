import core

type
  Gscalar* {.final.} = ref object of Gvalue
    ## Wrap a float as the graph value
    ## `getfloat=` changes the float value, useful during the graph construction
    ## `update` calls `getfloat=` and use `updated` to signal re-evaluation of the graph after graph construction
    sval: float

proc getfloat*(x: Gvalue): float = Gscalar(x).sval

proc `getfloat=`*(x: Gvalue, y: float) =
  let xs = Gscalar(x)
  xs.sval = y

proc update*(x: Gvalue, y: float) =
  x.getfloat = y
  x.updated

converter toGvalue*(x: float): Gvalue =
  result = Gscalar(sval: x)
  result.updated

method newOneOf*(x: Gscalar): Gvalue = Gscalar()
method valCopy*(z: Gscalar, x: Gscalar) = z.sval = x.sval

method `$`*(x: Gscalar): string = $x.sval

proc negsf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let z = Gscalar(v)
  z.sval = - x.sval

proc negsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  case i
  of 0:
    if zb == nil:
      return toGvalue(-1.0)
    else:
      return -zb
  else:
    raiseValueError("i must be 0, got: " & $i)

let gsneg = newGfunc(forward = negsf, backward = negsb, name = "-")

method `-`(x: Gscalar): Gvalue = Gscalar(inputs: @[Gvalue(x)], gfunc: gsneg)

proc addsf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let y = Gscalar(v.inputs[1])
  let z = Gscalar(v)
  z.sval = x.sval + y.sval

proc addsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  if zb == nil:
    return toGvalue(1.0)
  else:
    return zb

let gsadd = newGfunc(forward = addsf, backward = addsb, name = "+")

method `+`(x: Gscalar, y: Gscalar): Gvalue = Gscalar(inputs: @[Gvalue(x), y], gfunc: gsadd)

proc mulsf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let y = Gscalar(v.inputs[1])
  let z = Gscalar(v)
  z.sval = x.sval * y.sval

proc mulsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  if zb == nil:
    return z.inputs[1-i]
  else:
    return zb*z.inputs[1-i]

let gsmul = newGfunc(forward = mulsf, backward = mulsb, name = "*")

method `*`(x: Gscalar, y: Gscalar): Gvalue = Gscalar(inputs: @[Gvalue(x), y], gfunc: gsmul)

proc subsf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let y = Gscalar(v.inputs[1])
  let z = Gscalar(v)
  z.sval = x.sval - y.sval

proc subsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  case i
  of 0:
    if zb == nil:
      return toGvalue(1.0)
    else:
      return zb
  of 1:
    if zb == nil:
      return toGvalue(-1.0)
    else:
      return -zb
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

let gssub = newGfunc(forward = subsf, backward = subsb, name = "-")

method `-`(x: Gscalar, y: Gscalar): Gvalue = Gscalar(inputs: @[Gvalue(x), y], gfunc: gssub)

proc divsf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let y = Gscalar(v.inputs[1])
  let z = Gscalar(v)
  z.sval = x.sval / y.sval

proc divsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  # s = f(g(x))
  # ds = df/dg dg/dx dx
  # dz = dx/y - x/y^2 dy
  case i
  of 0:
    if zb == nil:
      return 1.0 / z.inputs[1]
    else:
      return zb / z.inputs[1]
  of 1:
    if zb == nil:
      return - z / z.inputs[1]
    else:
      return - zb * z / z.inputs[1]
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

let gsdiv = newGfunc(forward = divsf, backward = divsb, name = "/")

method `/`(x: Gscalar, y: Gscalar): Gvalue = Gscalar(inputs: @[Gvalue(x), y], gfunc: gsdiv)

when isMainModule:
  import math
  import std/assertions

  graphDebug = true

  let x = Gscalar()
  let y = Gscalar()
  let w = x-2.0
  let v = w+y
  let z = v*(-v)/w
  let dzdy = z.grad y

  func f(a, b: float): float = (a+b-2.0)*(2.0-a-b)/(a-2.0)
  func dfdb(a, b: float): float = -2.0*(a+b-2.0)/(a-2.0)

  let a = 1.1
  let b = 3.7
  let c = 1.3

  x.update a
  y.update b
  echo z.treeRepr
  echo dzdy.treeRepr
  z.eval
  dzdy.eval
  echo "z = ",z
  echo z.treeRepr
  echo "dzdy = ",dzdy
  echo dzdy.treeRepr

  dumpGradientList()

  doAssert almostEqual(z.getfloat, f(a,b))
  doAssert almostEqual(dzdy.getfloat, dfdb(a,b))

  y.update c
  z.eval
  dzdy.eval
  echo "z = ",z
  echo z.treeRepr
  echo "dzdy = ",dzdy
  echo dzdy.treeRepr
  doAssert almostEqual(z.getfloat, f(a,c))
  doAssert almostEqual(dzdy.getfloat, dfdb(a,c))

  # may need to change the following after we implement optimization passes
  doAssert gsneg.runCount == 4
  doAssert gsadd.runCount == 4
  doAssert gsmul.runCount == 6
  doAssert gssub.runCount == 1
  doAssert gsdiv.runCount == 3
