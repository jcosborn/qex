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

proc update*(x: Gvalue, y: float) =
  x.getfloat = y
  x.updated

converter toGvalue*(x: float): Gvalue =
  result = Gscalar(sval: x)
  result.updated

method newOneOf*(x: Gscalar): Gvalue = Gscalar()
method valCopy*(z: Gscalar, x: Gscalar) = z.sval = x.sval

method `$`*(x: Gscalar): string = $x.sval

proc getint*(x: Gvalue): int = Gint(x).ival

proc `getint=`*(x: Gvalue, y: int) =
  let xs = Gint(x)
  xs.ival = y

proc update*(x: Gvalue, y: int) =
  x.getint = y
  x.updated

converter toGvalue*(x: int): Gvalue =
  result = Gint(ival: x)
  result.updated

method newOneOf*(x: Gint): Gvalue = Gint()
method valCopy*(z: Gint, x: Gint) = z.ival = x.ival

method `$`*(x: Gint): string = $x.ival

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

method `not`*(x: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("not(" & $x & ")")
method `and`*(x: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("and(" & $x & ", " & $y & ")")
method `or`*(x: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("or(" & $x & ", " & $y & ")")
method `<`*(x: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("`<`(" & $x & ", " & $y & ")")
method equal*(x: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("equal(" & $x & ", " & $y & ")")

method cond*(c: Gvalue, x: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("cond(" & $c & "," & $x & ", " & $y & ")")

proc `>`*(x, y: Gvalue): Gvalue = not(x < y)
proc `>=`*(x, y: Gvalue): Gvalue = x > y or equal(x,y)
proc `<=`*(x, y: Gvalue): Gvalue = x < y or equal(x,y)
proc `!=`*(x, y: Gvalue): Gvalue = not equal(x,y)
proc `xor`*(x, y: Gvalue): Gvalue = not equal(not(x), not(y))  # uses `not` to convert to 0/1

proc notsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  case i
  of 0:
    return toGvalue(0.0)
  else:
    raiseValueError("i must be 0, got: " & $i)

proc notsf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let z = Gscalar(v)
  z.sval = if x.sval == 0.0: 1.0 else: 0.0

let nots = newGfunc(forward = notsf, backward = notsb, name = "nots")

method `not`*(x: Gscalar): Gvalue = Gscalar(inputs: @[Gvalue(x)], gfunc: nots)

proc andsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  case i
  of 0, 1:
    return toGvalue(0.0)
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc andsf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let y = Gscalar(v.inputs[1])
  let z = Gscalar(v)
  z.sval = if x.sval == 0.0 or y.sval == 0.0: 0.0 else: 1.0

let ands = newGfunc(forward = andsf, backward = andsb, name = "ands")

method `and`*(x: Gscalar, y: Gscalar): Gvalue = Gscalar(inputs: @[Gvalue(x), y], gfunc: ands)

proc orsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  case i
  of 0, 1:
    return toGvalue(0.0)
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc orsf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let y = Gscalar(v.inputs[1])
  let z = Gscalar(v)
  z.sval = if x.sval == 0.0 and y.sval == 0.0: 0.0 else: 1.0

let ors = newGfunc(forward = orsf, backward = orsb, name = "ors")

method `or`*(x: Gscalar, y: Gscalar): Gvalue = Gscalar(inputs: @[Gvalue(x), y], gfunc: ors)

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

proc condsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  case i
  of 0:
    return toGvalue(0.0)
  of 1:
    if zb == nil:
      # the output must be a scalar, otherwise crash later
      return cond(z.inputs[0], toGvalue(1.0), toGvalue(0.0))
    else:
      return cond(z.inputs[0], zb, zb.newOneOf)
  of 2:
    if zb == nil:
      # the output must be a scalar, otherwise crash later
      return cond(z.inputs[0], toGvalue(0.0), toGvalue(1.0))
    else:
      return cond(z.inputs[0], zb.newOneOf, zb)
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc condsf(v: Gvalue) =
  let c = Gscalar(v.inputs[0])
  if c.sval == 0.0:
    v.valCopy v.inputs[2]
  else:
    v.valCopy v.inputs[1]

let conds = newGfunc(forward = condsf, backward = condsb, name = "conds")

method cond*(c: Gscalar, x: Gvalue, y: Gvalue): Gvalue =
  result = x.newOneOf
  result.inputs = @[Gvalue(c), x, y]
  result.gfunc = conds

proc notib(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  case i
  of 0:
    return toGvalue(0)
  else:
    raiseValueError("i must be 0, got: " & $i)

proc notif(v: Gvalue) =
  let x = Gint(v.inputs[0])
  let z = Gint(v)
  z.ival = if x.ival == 0: 1 else: 0

let noti = newGfunc(forward = notif, backward = notib, name = "noti")

method `not`*(x: Gint): Gvalue = Gint(inputs: @[Gvalue(x)], gfunc: noti)

proc andib(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  case i
  of 0, 1:
    return toGvalue(0)
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc andif(v: Gvalue) =
  let x = Gint(v.inputs[0])
  let y = Gint(v.inputs[1])
  let z = Gint(v)
  z.ival = if x.ival == 0 or y.ival == 0: 0 else: 1

let andi = newGfunc(forward = andif, backward = andib, name = "andi")

method `and`*(x: Gint, y: Gint): Gvalue = Gint(inputs: @[Gvalue(x), y], gfunc: andi)

proc orib(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  case i
  of 0, 1:
    return toGvalue(0)
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc orif(v: Gvalue) =
  let x = Gint(v.inputs[0])
  let y = Gint(v.inputs[1])
  let z = Gint(v)
  z.ival = if x.ival == 0 and y.ival == 0: 0 else: 1

let ori = newGfunc(forward = orif, backward = orib, name = "ori")

method `or`*(x: Gint, y: Gint): Gvalue = Gint(inputs: @[Gvalue(x), y], gfunc: ori)

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

proc condib(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  case i
  of 0:
    return toGvalue(0)
  of 1:
    if zb == nil:
      # the output must be a scalar, otherwise crash later
      return cond(z.inputs[0], toGvalue(1.0), toGvalue(0.0))
    else:
      return cond(z.inputs[0], zb, zb.newOneOf)
  of 2:
    if zb == nil:
      # the output must be a scalar, otherwise crash later
      return cond(z.inputs[0], toGvalue(0.0), toGvalue(1.0))
    else:
      return cond(z.inputs[0], zb.newOneOf, zb)
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc condif(v: Gvalue) =
  let c = Gint(v.inputs[0])
  if c.ival == 0:
    v.valCopy v.inputs[2]
  else:
    v.valCopy v.inputs[1]

let condi = newGfunc(forward = condif, backward = condib, name = "condi")

method cond*(c: Gint, x: Gvalue, y: Gvalue): Gvalue =
  result = x.newOneOf
  result.inputs = @[Gvalue(c), x, y]
  result.gfunc = condi

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
