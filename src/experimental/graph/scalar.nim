import core
from math import exp

type
  Gscalar* {.final.} = ref object of Gvalue
    ## Wrap a float as the graph value
    ## `getfloat=` changes the float value, useful during the graph construction
    ## `update` calls `getfloat=` and use `updated` to signal re-evaluation of the graph after graph construction
    sval: float
  Gint* {.final.} = ref object of Gvalue
    ival: int

proc getfloat*(x: Gvalue): float = Gscalar(x).sval

proc `getfloat=`*(x: Gvalue, y: float) =
  let xs = Gscalar(x)
  xs.sval = y

method update*(x: Gscalar, y: float) =
  x.getfloat = y
  x.updated

converter toGvalue*(x: float): Gvalue =
  result = Gscalar(sval: x)
  result.updated

method newOneOf*(x: Gscalar): Gvalue = Gscalar()
method valCopy*(z: Gscalar, x: Gscalar) = z.sval = x.sval

method `$`*(x: Gscalar): string = $x.sval

method isZero*(x: Gscalar): bool = x.sval == 0.0

proc `getfloat=`*(x: Gvalue, y: int) =
  let xs = Gscalar(x)
  xs.sval = float(y)

method update*(x: Gscalar, y: int) =
  x.getfloat = y
  x.updated

proc getint*(x: Gvalue): int = Gint(x).ival

proc `getint=`*(x: Gvalue, y: int) =
  let xs = Gint(x)
  xs.ival = y

method update*(x: Gint, y: int) =
  x.getint = y
  x.updated

converter toGvalue*(x: int): Gvalue =
  result = Gint(ival: x)
  result.updated

method newOneOf*(x: Gint): Gvalue = Gint()
method valCopy*(z: Gint, x: Gint) = z.ival = x.ival

method `$`*(x: Gint): string = $x.ival

method isZero*(x: Gint): bool = x.ival == 0

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

method `-`*(x: Gscalar): Gvalue = Gscalar(inputs: @[Gvalue(x)], gfunc: gsneg)

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

method `+`*(x: Gscalar, y: Gscalar): Gvalue = Gscalar(inputs: @[Gvalue(x), y], gfunc: gsadd)

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

method `*`*(x: Gscalar, y: Gscalar): Gvalue = Gscalar(inputs: @[Gvalue(x), y], gfunc: gsmul)

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

method `-`*(x: Gscalar, y: Gscalar): Gvalue = Gscalar(inputs: @[Gvalue(x), y], gfunc: gssub)

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

method `/`*(x: Gscalar, y: Gscalar): Gvalue = Gscalar(inputs: @[Gvalue(x), y], gfunc: gsdiv)

proc expsf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let z = Gscalar(v)
  z.sval = exp(x.sval)

proc expsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  case i
  of 0:
    if zb == nil:
      return z
    else:
      return zb*z
  else:
    raiseValueError("i must be 0, got: " & $i)

let exps = newGfunc(forward = expsf, backward = expsb, name = "exps")

method exp*(x: Gscalar): Gvalue = Gscalar(inputs: @[Gvalue(x)], gfunc: exps)

method `<`*(x: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("`<`(" & $x & ", " & $y & ")")
method equal*(x: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("equal(" & $x & ", " & $y & ")")

proc newFalse(x: Gvalue): Gvalue =
  result = x.newOneOf
  result.update 0
proc newTrue(x: Gvalue): Gvalue =
  result = x.newOneOf
  result.update 1

proc `not`*(x: Gvalue): Gvalue = cond(x, x.newFalse, x.newTrue)
proc `and`*(x: Gvalue, y: Gvalue): Gvalue = cond(x, y, y.newFalse)  ## return type follows the second operand, as does cond
proc `or`*(x: Gvalue, y: Gvalue): Gvalue = cond(x, y.newTrue, y)  ## return type follows the second operand, as does cond
proc `xor`*(x: Gvalue, y: Gvalue): Gvalue = cond(x, not(y), y)  ## return type follows the second operand, as does cond

proc `>`*(x, y: Gvalue): Gvalue = not(x < y)
proc `>=`*(x, y: Gvalue): Gvalue = x > y or equal(x,y)
proc `<=`*(x, y: Gvalue): Gvalue = x < y or equal(x,y)

proc ltsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  case i
  of 0, 1:
    return toGvalue(0.0)
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc ltsf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let y = Gscalar(v.inputs[1])
  let z = Gscalar(v)
  z.sval = if x.sval < y.sval: 1.0 else: 0.0

let lts = newGfunc(forward = ltsf, backward = ltsb, name = "lts")

method `<`*(x: Gscalar, y: Gscalar): Gvalue = Gscalar(inputs: @[Gvalue(x), y], gfunc: lts)

proc equalsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  case i
  of 0, 1:
    return toGvalue(0.0)
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc equalsf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let y = Gscalar(v.inputs[1])
  let z = Gscalar(v)
  z.sval = if x.sval == y.sval: 1.0 else: 0.0

let equals = newGfunc(forward = equalsf, backward = equalsb, name = "equals")

method equal*(x: Gscalar, y: Gscalar): Gvalue = Gscalar(inputs: @[Gvalue(x), y], gfunc: equals)

proc ltib(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  case i
  of 0, 1:
    return toGvalue(0)
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc ltif(v: Gvalue) =
  let x = Gint(v.inputs[0])
  let y = Gint(v.inputs[1])
  let z = Gint(v)
  z.ival = if x.ival < y.ival: 1 else: 0

let lti = newGfunc(forward = ltif, backward = ltib, name = "lti")

method `<`*(x: Gint, y: Gint): Gvalue = Gint(inputs: @[Gvalue(x), y], gfunc: lti)

proc equalib(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  case i
  of 0, 1:
    return toGvalue(0)
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc equalif(v: Gvalue) =
  let x = Gint(v.inputs[0])
  let y = Gint(v.inputs[1])
  let z = Gint(v)
  z.ival = if x.ival == y.ival: 1 else: 0

let equali = newGfunc(forward = equalif, backward = equalib, name = "equali")

method equal*(x: Gint, y: Gint): Gvalue = Gint(inputs: @[Gvalue(x), y], gfunc: equali)

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
