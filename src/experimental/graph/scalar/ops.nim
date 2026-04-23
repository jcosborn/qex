import ../core
import types
import ../support/op
import math

template defineScalarUnaryForward(forwardName: untyped,
                                  label: static[string],
                                  forwardBody: untyped) =
  proc forwardName(v: Gvalue) =
    let view {.inject.} = v.requireUnaryNodeView(label)
    let x {.inject.} = view.x
    let z {.inject.} = v
    forwardBody

template defineScalarBinaryForward(forwardName: untyped,
                                   label: static[string],
                                   forwardBody: untyped) =
  proc forwardName(v: Gvalue) =
    let view {.inject.} = v.requireBinaryNodeView(label)
    let x {.inject.} = view.x
    let y {.inject.} = view.y
    let z {.inject.} = v
    forwardBody

template defineScalarComparisonForward(forwardName: untyped,
                                       label: static[string],
                                       predicate: untyped) =
  proc forwardName(v: Gvalue) =
    let view {.inject.} = v.requireBinaryNodeView(label)
    let x {.inject.} = view.x
    let y {.inject.} = view.y
    let z {.inject.} = v
    z.getfloat = if predicate: 1.0 else: 0.0

template defineIntComparisonForward(forwardName: untyped,
                                    label: static[string],
                                    predicate: untyped) =
  proc forwardName(v: Gvalue) =
    let view {.inject.} = v.requireBinaryNodeView(label)
    let x {.inject.} = view.x
    let y {.inject.} = view.y
    let z {.inject.} = v
    z.getint = if predicate: 1 else: 0

proc `-`*(x: Gscalar): Gscalar
proc `+`*(x: Gscalar, y: Gscalar): Gscalar
proc `*`*(x: Gscalar, y: Gscalar): Gscalar
proc `-`*(x: Gscalar, y: Gscalar): Gscalar
proc `/`*(x: Gscalar, y: Gscalar): Gscalar
proc exp*(x: Gscalar): Gscalar
proc `<`*(x: Gscalar, y: Gscalar): Gscalar
proc equal*(x: Gscalar, y: Gscalar): Gscalar
proc `<`*(x: Gint, y: Gint): Gint
proc equal*(x: Gint, y: Gint): Gint
type CondSelector = Gscalar | Gint
proc cond*[C: CondSelector, T: Gvalue](c: C, x: T, y: T): T
proc cond*(c: CondSelector, x: Gscalar, y: float): Gscalar
proc cond*(c: CondSelector, x: float, y: Gscalar): Gscalar
proc cond*(c: CondSelector, x: Gscalar, y: int): Gscalar
proc cond*(c: CondSelector, x: int, y: Gscalar): Gscalar
proc cond*(c: CondSelector, x: Gint, y: int): Gint
proc cond*(c: CondSelector, x: int, y: Gint): Gint

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

proc scalarLiteral(anchor: Gvalue,
                   value: float): Gscalar =
  scalarLeafLike(anchor, value)

proc scalarLiteral(anchor: Gvalue,
                   value: int): Gscalar =
  scalarLeafLike(anchor, value)

proc numericLiteral(anchor: Gvalue,
                    value: int): Gvalue =
  numericLeafLike(anchor, value)

proc numericLiteral(anchor: Gvalue,
                    value: float): Gvalue =
  numericLeafLike(anchor, value)

proc `+`*[T: Gvalue](x: T, y: float): auto =
  x + scalarLiteral(x, y)

proc `+`*[T: Gvalue](x: float, y: T): auto =
  scalarLiteral(y, x) + y

proc `+`*[T: Gvalue](x: T, y: int): auto =
  x + scalarLiteral(x, y)

proc `+`*[T: Gvalue](x: int, y: T): auto =
  scalarLiteral(y, x) + y

proc `-`*[T: Gvalue](x: T, y: float): auto =
  x - scalarLiteral(x, y)

proc `-`*[T: Gvalue](x: float, y: T): auto =
  scalarLiteral(y, x) - y

proc `-`*[T: Gvalue](x: T, y: int): auto =
  x - scalarLiteral(x, y)

proc `-`*[T: Gvalue](x: int, y: T): auto =
  scalarLiteral(y, x) - y

proc `*`*[T: Gvalue](x: T, y: float): auto =
  x * scalarLiteral(x, y)

proc `*`*[T: Gvalue](x: float, y: T): auto =
  scalarLiteral(y, x) * y

proc `*`*[T: Gvalue](x: T, y: int): auto =
  x * scalarLiteral(x, y)

proc `*`*[T: Gvalue](x: int, y: T): auto =
  scalarLiteral(y, x) * y

proc `/`*[T: Gvalue](x: T, y: float): auto =
  x / scalarLiteral(x, y)

proc `/`*[T: Gvalue](x: float, y: T): auto =
  scalarLiteral(y, x) / y

proc `/`*[T: Gvalue](x: T, y: int): auto =
  x / scalarLiteral(x, y)

proc `/`*[T: Gvalue](x: int, y: T): auto =
  scalarLiteral(y, x) / y

proc requireNonZeroDivisor(y: Gvalue, z: Gvalue) =
  if y.isZero:
    raiseValueError("division by zero:\n" & z.nodeRepr)

defineScalarUnaryForward(negsf, "- forward"):
  z.getfloat = -x.getfloat
proc negsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z.requireUnaryNodeView("- backward")
  discard dep
  unaryBackwardCase("- backward", i):
    return affineUpstream(zb, -1.0, z)
defineUnaryGraphOp(gsneg, `-`, Gscalar, x, scalarNodeLike(x), negsf, negsb, "-")

defineScalarBinaryForward(addsf, "+ forward"):
  z.getfloat = x.getfloat + y.getfloat
proc addsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z.requireBinaryNodeView("+ backward")
  discard dep
  binaryBackwardCase("+ backward", i,
    block:
      return affineUpstream(zb, 1.0, z),
    block:
      return affineUpstream(zb, 1.0, z))
defineBinaryGraphOp(gsadd, `+`, Gscalar, Gscalar, x, y, scalarNodeLike(x), addsf, addsb, "+")

defineScalarBinaryForward(mulsf, "* forward"):
  z.getfloat = x.getfloat * y.getfloat
proc mulsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView(Gscalar, Gscalar, "* backward")
  discard dep
  binaryBackwardCase("* backward", i,
    block:
      return scaledUpstreamOr(zb, Gscalar, view.y, "* backward"),
    block:
      return scaledUpstreamOr(zb, Gscalar, view.x, "* backward"))
defineBinaryGraphOp(gsmul, `*`, Gscalar, Gscalar, x, y, scalarNodeLike(x), mulsf, mulsb, "*")

defineScalarBinaryForward(subsf, "- forward"):
  z.getfloat = x.getfloat - y.getfloat
proc subsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z.requireBinaryNodeView("- backward")
  discard dep
  binaryBackwardCase("- backward", i,
    block:
      return affineUpstream(zb, 1.0, z),
    block:
      return affineUpstream(zb, -1.0, z))
defineBinaryGraphOp(gssub, `-`, Gscalar, Gscalar, x, y, scalarNodeLike(x), subsf, subsb, "-")

defineScalarBinaryForward(divsf, "/ forward"):
  requireNonZeroDivisor(y, z)
  z.getfloat = x.getfloat / y.getfloat
proc divsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView(Gscalar, Gscalar, "/ backward")
  discard dep
  # d(x / y) = dx / y - x * dy / y^2
  binaryBackwardCase("/ backward", i,
    block:
      return scaledUpstreamOr(zb, Gscalar, 1.0 / view.y, "/ backward"),
    block:
      return scaledUpstreamOr(zb, Gscalar, -Gscalar(z) / view.y, "/ backward"))
defineBinaryGraphOp(gsdiv, `/`, Gscalar, Gscalar, x, y, scalarNodeLike(x), divsf, divsb, "/")

defineScalarUnaryForward(expsf, "exp forward"):
  z.getfloat = math.exp(x.getfloat)
proc expsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z.requireUnaryNodeView("exp backward")
  discard dep
  unaryBackwardCase("exp backward", i):
    return scaledUpstreamOr(zb, Gscalar, Gscalar(z), "exp backward")
defineUnaryGraphOp(exps, exp, Gscalar, x, scalarNodeLike(x), expsf, expsb, "exps")

method addLike*(prototype: Gscalar, x: Gvalue, y: Gvalue): Gvalue =
  discard prototype
  x.requireScalar("scalar gradient add left") + y.requireScalar("scalar gradient add right")

proc `+`*(x: Gscalar, y: Gvalue): Gscalar =
  x + y.requireScalar("scalar + right")

method scaleLike*(contribution: Gscalar, upstream: Gvalue): Gvalue =
  upstream.requireScalar("scalar scale upstream") * contribution

proc `*`*(x: Gscalar, y: Gvalue): Gscalar =
  x * y.requireScalar("scalar * right")

proc `-`*(x: Gscalar, y: Gvalue): Gscalar =
  x - y.requireScalar("scalar - right")

proc `/`*(x: Gscalar, y: Gvalue): Gscalar =
  x / y.requireScalar("scalar / right")

method addLike*(prototype: Gint, x: Gvalue, y: Gvalue): Gvalue =
  discard prototype
  let left = x.requireInt("int gradient add left")
  let right = y.requireInt("int gradient add right")
  if left.isZero:
    return right
  if right.isZero:
    return left
  raiseValueError("cannot accumulate non-zero int graph gradients")

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

defineScalarComparisonForward(ltsf, "< forward", x.getfloat < y.getfloat)
proc ltsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard zb
  discard z.requireBinaryNodeView("< backward")
  discard dep
  zeroComparisonGrad(numericLiteral(z, 0.0), i)
defineBinaryGraphOp(lts, `<`, Gscalar, Gscalar, x, y, scalarNodeLike(x), ltsf, ltsb, "lts")

defineScalarComparisonForward(equalsf, "equal forward", x.getfloat == y.getfloat)
proc equalsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard zb
  discard z.requireBinaryNodeView("equal backward")
  discard dep
  zeroComparisonGrad(numericLiteral(z, 0.0), i)
defineBinaryGraphOp(equals, equal, Gscalar, Gscalar, x, y, scalarNodeLike(x), equalsf, equalsb, "equals")

defineIntComparisonForward(ltif, "int < forward", x.getint < y.getint)
proc ltib(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard zb
  discard z.requireBinaryNodeView("int < backward")
  discard dep
  zeroComparisonGrad(numericLiteral(z, 0), i)
defineBinaryGraphOp(lti, `<`, Gint, Gint, x, y, intNodeLike(x), ltif, ltib, "lti")

defineIntComparisonForward(equalif, "int equal forward", x.getint == y.getint)
proc equalib(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard zb
  discard z.requireBinaryNodeView("int equal backward")
  discard dep
  zeroComparisonGrad(numericLiteral(z, 0), i)
defineBinaryGraphOp(equali, equal, Gint, Gint, x, y, intNodeLike(x), equalif, equalib, "equali")

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
