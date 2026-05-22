import ../core
import types
import ../support/op
import math

type CondSelector = Gscalar | Gint
proc `-`*(x: Gscalar): Gscalar
proc `*`*(x: Gscalar, y: Gscalar): Gscalar
proc `/`*(x: Gscalar, y: Gscalar): Gscalar

proc zeroComparisonGrad[T: Gvalue](zero: T, i: int): T =
  case i
  of 0, 1:
    zero
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc affineUpstream(zb: Gvalue,
                    scale: float,
                    anchor: Gvalue): Gscalar =
  if zb == nil:
    return scalarLeafLike(anchor, scale)
  let upstream = requireUpstream(zb, Gscalar, "scalar backward")
  if scale == 1.0:
    return upstream
  if scale == -1.0:
    return -upstream
  scalarLeafLike(anchor, scale) * upstream

proc negsf(v: Gvalue) =
  let view = v.requireUnaryNodeView(Gscalar, "- forward")
  let x = view.x
  let z = v.requireScalar("- forward result")
  z.sval = -x.sval

proc negsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard dep
  unaryBackwardCase("- backward", i):
    return affineUpstream(zb, -1.0, z)

let gsneg = newGfunc(forward = negsf, backward = negsb, name = "-")

proc `-`*(x: Gscalar): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x)], gsneg, "-")

proc addsf(v: Gvalue) =
  let view = v.requireBinaryNodeView(Gscalar, Gscalar, "+ forward")
  let x = view.x
  let y = view.y
  let z = v.requireScalar("+ forward result")
  z.sval = x.sval + y.sval

proc addsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard dep
  binaryBackwardCase("+ backward", i,
    block:
      return affineUpstream(zb, 1.0, z),
    block:
      return affineUpstream(zb, 1.0, z))

let gsadd = newGfunc(forward = addsf, backward = addsb, name = "+")

proc `+`*(x: Gscalar, y: Gscalar): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x), Gvalue(y)], gsadd, "+")

proc mulsf(v: Gvalue) =
  let view = v.requireBinaryNodeView(Gscalar, Gscalar, "* forward")
  let x = view.x
  let y = view.y
  let z = v.requireScalar("* forward result")
  z.sval = x.sval * y.sval

proc mulsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView(Gscalar, Gscalar, "* backward")
  discard dep
  binaryBackwardCase("* backward", i,
    block:
      return scaledUpstreamOr(zb, Gscalar, view.y, "* backward"),
    block:
      return scaledUpstreamOr(zb, Gscalar, view.x, "* backward"))

let gsmul = newGfunc(forward = mulsf, backward = mulsb, name = "*")

proc `*`*(x: Gscalar, y: Gscalar): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x), Gvalue(y)], gsmul, "*")

proc subsf(v: Gvalue) =
  let view = v.requireBinaryNodeView(Gscalar, Gscalar, "- forward")
  let x = view.x
  let y = view.y
  let z = v.requireScalar("- forward result")
  z.sval = x.sval - y.sval

proc subsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard dep
  binaryBackwardCase("- backward", i,
    block:
      return affineUpstream(zb, 1.0, z),
    block:
      return affineUpstream(zb, -1.0, z))

let gssub = newGfunc(forward = subsf, backward = subsb, name = "-")

proc `-`*(x: Gscalar, y: Gscalar): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x), Gvalue(y)], gssub, "-")

proc divsf(v: Gvalue) =
  let view = v.requireBinaryNodeView(Gscalar, Gscalar, "/ forward")
  let x = view.x
  let y = view.y
  let z = v.requireScalar("/ forward result")
  if y.isZero:
    raiseValueError("division by zero:\n" & z.nodeRepr)
  z.sval = x.sval / y.sval

proc divsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView(Gscalar, Gscalar, "/ backward")
  discard dep
  # d(x / y) = dx / y - x * dy / y^2
  binaryBackwardCase("/ backward", i,
    block:
      return scaledUpstreamOr(
        zb,
        Gscalar,
        scalarLeafLike(view.y, 1.0) / view.y,
        "/ backward"),
    block:
      return scaledUpstreamOr(zb, Gscalar, -Gscalar(z) / view.y, "/ backward"))

let gsdiv = newGfunc(forward = divsf, backward = divsb, name = "/")

proc `/`*(x: Gscalar, y: Gscalar): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x), Gvalue(y)], gsdiv, "/")

proc expsf(v: Gvalue) =
  let view = v.requireUnaryNodeView(Gscalar, "exp forward")
  let x = view.x
  let z = v.requireScalar("exp forward result")
  z.sval = math.exp(x.sval)

proc expsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard dep
  unaryBackwardCase("exp backward", i):
    return scaledUpstreamOr(zb, Gscalar, Gscalar(z), "exp backward")

let exps = newGfunc(forward = expsf, backward = expsb, name = "exps")

proc exp*(x: Gscalar): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x)], exps, "exps")

proc `+`*(x: Gscalar, y: float): Gscalar =
  x + scalarLeafLike(x, y)

proc `+`*(x: float, y: Gscalar): Gscalar =
  scalarLeafLike(y, x) + y

proc `+`*(x: Gscalar, y: int): Gscalar =
  x + scalarLeafLike(x, y)

proc `+`*(x: int, y: Gscalar): Gscalar =
  scalarLeafLike(y, x) + y

proc `-`*(x: Gscalar, y: float): Gscalar =
  x - scalarLeafLike(x, y)

proc `-`*(x: float, y: Gscalar): Gscalar =
  scalarLeafLike(y, x) - y

proc `-`*(x: Gscalar, y: int): Gscalar =
  x - scalarLeafLike(x, y)

proc `-`*(x: int, y: Gscalar): Gscalar =
  scalarLeafLike(y, x) - y

proc `*`*(x: Gscalar, y: float): Gscalar =
  x * scalarLeafLike(x, y)

proc `*`*(x: float, y: Gscalar): Gscalar =
  scalarLeafLike(y, x) * y

proc `*`*(x: Gscalar, y: int): Gscalar =
  x * scalarLeafLike(x, y)

proc `*`*(x: int, y: Gscalar): Gscalar =
  scalarLeafLike(y, x) * y

proc `/`*(x: Gscalar, y: float): Gscalar =
  x / scalarLeafLike(x, y)

proc `/`*(x: float, y: Gscalar): Gscalar =
  scalarLeafLike(y, x) / y

proc `/`*(x: Gscalar, y: int): Gscalar =
  x / scalarLeafLike(x, y)

proc `/`*(x: int, y: Gscalar): Gscalar =
  scalarLeafLike(y, x) / y

method addLike*(prototype: Gscalar, x: Gvalue, y: Gvalue): Gvalue =
  discard prototype
  x.requireScalar("scalar gradient add left") + y.requireScalar("scalar gradient add right")

method scaleLike*(contribution: Gscalar, upstream: Gvalue): Gvalue =
  upstream.requireScalar("scalar scale upstream") * contribution

method addLike*(prototype: Gint, x: Gvalue, y: Gvalue): Gvalue =
  discard prototype
  let left = x.requireInt("int gradient add left")
  let right = y.requireInt("int gradient add right")
  if left.isZero:
    return right
  if right.isZero:
    return left
  raiseValueError("cannot accumulate non-zero int graph gradients")

proc cond*[C: CondSelector, T: Gvalue](c: C, x: T, y: T): T =
  newCondNode(c, x, y)

proc cond*(c: CondSelector, x: Gscalar, y: float): Gscalar =
  cond(c, x, scalarLeafLike(x, y))

proc cond*(c: CondSelector, x: float, y: Gscalar): Gscalar =
  cond(c, scalarLeafLike(y, x), y)

proc cond*(c: CondSelector, x: Gscalar, y: int): Gscalar =
  cond(c, x, scalarLeafLike(x, y))

proc cond*(c: CondSelector, x: int, y: Gscalar): Gscalar =
  cond(c, scalarLeafLike(y, x), y)

proc cond*(c: CondSelector, x: Gint, y: int): Gint =
  cond(c, x, intLeafLike(x, y))

proc cond*(c: CondSelector, x: int, y: Gint): Gint =
  cond(c, intLeafLike(y, x), y)

proc falseValue(x: Gscalar): Gscalar =
  scalarLeafLike(x, 0.0)

proc trueValue(x: Gscalar): Gscalar =
  scalarLeafLike(x, 1.0)

proc falseValue(x: Gint): Gint =
  intLeafLike(x, 0)

proc trueValue(x: Gint): Gint =
  intLeafLike(x, 1)

proc `not`*(x: Gscalar): Gscalar = cond(x, x.falseValue, x.trueValue)
proc `not`*(x: Gint): Gint = cond(x, x.falseValue, x.trueValue)

proc `and`*(x: Gscalar, y: Gscalar): Gscalar = cond(x, y, y.falseValue)
proc `and`*(x: Gscalar, y: Gint): Gint = cond(x, y, y.falseValue)
proc `and`*(x: Gint, y: Gscalar): Gscalar = cond(x, y, y.falseValue)
proc `and`*(x: Gint, y: Gint): Gint = cond(x, y, y.falseValue)

proc `or`*(x: Gscalar, y: Gscalar): Gscalar = cond(x, y.trueValue, y)
proc `or`*(x: Gscalar, y: Gint): Gint = cond(x, y.trueValue, y)
proc `or`*(x: Gint, y: Gscalar): Gscalar = cond(x, y.trueValue, y)
proc `or`*(x: Gint, y: Gint): Gint = cond(x, y.trueValue, y)

proc `xor`*(x: Gscalar, y: Gscalar): Gscalar = cond(x, not(y), y)
proc `xor`*(x: Gscalar, y: Gint): Gint = cond(x, not(y), y)
proc `xor`*(x: Gint, y: Gscalar): Gscalar = cond(x, not(y), y)
proc `xor`*(x: Gint, y: Gint): Gint = cond(x, not(y), y)

proc ltsf(v: Gvalue) =
  let view = v.requireBinaryNodeView(Gscalar, Gscalar, "< forward")
  let x = view.x
  let y = view.y
  let z = v.requireScalar("< forward result")
  z.sval = if x.sval < y.sval: 1.0 else: 0.0

proc ltsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard zb
  discard dep
  zeroComparisonGrad(numericLeafLike(z, 0.0), i)

let lts = newGfunc(forward = ltsf, backward = ltsb, name = "lts")

proc `<`*(x: Gscalar, y: Gscalar): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x), Gvalue(y)], lts, "lts")

proc equalsf(v: Gvalue) =
  let view = v.requireBinaryNodeView(Gscalar, Gscalar, "equal forward")
  let x = view.x
  let y = view.y
  let z = v.requireScalar("equal forward result")
  z.sval = if x.sval == y.sval: 1.0 else: 0.0

proc equalsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard zb
  discard dep
  zeroComparisonGrad(numericLeafLike(z, 0.0), i)

let equals = newGfunc(forward = equalsf, backward = equalsb, name = "equals")

proc equal*(x: Gscalar, y: Gscalar): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x), Gvalue(y)], equals, "equals")

proc ltif(v: Gvalue) =
  let view = v.requireBinaryNodeView(Gint, Gint, "int < forward")
  let x = view.x
  let y = view.y
  let z = v.requireInt("int < forward result")
  z.ival = if x.ival < y.ival: 1 else: 0

proc ltib(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard zb
  discard dep
  zeroComparisonGrad(numericLeafLike(z, 0), i)

let lti = newGfunc(forward = ltif, backward = ltib, name = "lti")

proc `<`*(x: Gint, y: Gint): Gint =
  graphNode(intNodeLike(x), @[Gvalue(x), Gvalue(y)], lti, "lti")

proc equalif(v: Gvalue) =
  let view = v.requireBinaryNodeView(Gint, Gint, "int equal forward")
  let x = view.x
  let y = view.y
  let z = v.requireInt("int equal forward result")
  z.ival = if x.ival == y.ival: 1 else: 0

proc equalib(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard zb
  discard dep
  zeroComparisonGrad(numericLeafLike(z, 0), i)

let equali = newGfunc(forward = equalif, backward = equalib, name = "equali")

proc equal*(x: Gint, y: Gint): Gint =
  graphNode(intNodeLike(x), @[Gvalue(x), Gvalue(y)], equali, "equali")

proc `>`*(x, y: Gscalar): auto = y < x
proc `>=`*(x, y: Gscalar): auto = (y < x) or equal(x, y)
proc `<=`*(x, y: Gscalar): auto = (x < y) or equal(x, y)

proc `>`*(x, y: Gint): auto = y < x
proc `>=`*(x, y: Gint): auto = (y < x) or equal(x, y)
proc `<=`*(x, y: Gint): auto = (x < y) or equal(x, y)

proc `<`*(x: Gscalar, y: float): Gscalar = x < scalarLeafLike(x, y)
proc `<`*(x: float, y: Gscalar): Gscalar = scalarLeafLike(y, x) < y
proc `<`*(x: Gscalar, y: int): Gscalar = x < scalarLeafLike(x, y)
proc `<`*(x: int, y: Gscalar): Gscalar = scalarLeafLike(y, x) < y

proc equal*(x: Gscalar, y: float): Gscalar = equal(x, scalarLeafLike(x, y))
proc equal*(x: float, y: Gscalar): Gscalar = equal(scalarLeafLike(y, x), y)
proc equal*(x: Gscalar, y: int): Gscalar = equal(x, scalarLeafLike(x, y))
proc equal*(x: int, y: Gscalar): Gscalar = equal(scalarLeafLike(y, x), y)

proc `>`*(x: Gscalar, y: float): Gscalar = y < x
proc `>`*(x: float, y: Gscalar): Gscalar = y < x
proc `>`*(x: Gscalar, y: int): Gscalar = y < x
proc `>`*(x: int, y: Gscalar): Gscalar = y < x

proc `>=`*(x: Gscalar, y: float): Gscalar = (y < x) or equal(x, y)
proc `>=`*(x: float, y: Gscalar): Gscalar = (y < x) or equal(x, y)
proc `>=`*(x: Gscalar, y: int): Gscalar = (y < x) or equal(x, y)
proc `>=`*(x: int, y: Gscalar): Gscalar = (y < x) or equal(x, y)

proc `<=`*(x: Gscalar, y: float): Gscalar = (x < y) or equal(x, y)
proc `<=`*(x: float, y: Gscalar): Gscalar = (x < y) or equal(x, y)
proc `<=`*(x: Gscalar, y: int): Gscalar = (x < y) or equal(x, y)
proc `<=`*(x: int, y: Gscalar): Gscalar = (x < y) or equal(x, y)

proc `<`*(x: Gint, y: int): Gint = x < intLeafLike(x, y)
proc `<`*(x: int, y: Gint): Gint = intLeafLike(y, x) < y

proc equal*(x: Gint, y: int): Gint = equal(x, intLeafLike(x, y))
proc equal*(x: int, y: Gint): Gint = equal(intLeafLike(y, x), y)

proc `>`*(x: Gint, y: int): Gint = y < x
proc `>`*(x: int, y: Gint): Gint = y < x

proc `>=`*(x: Gint, y: int): Gint = (y < x) or equal(x, y)
proc `>=`*(x: int, y: Gint): Gint = (y < x) or equal(x, y)

proc `<=`*(x: Gint, y: int): Gint = (x < y) or equal(x, y)
proc `<=`*(x: int, y: Gint): Gint = (x < y) or equal(x, y)
