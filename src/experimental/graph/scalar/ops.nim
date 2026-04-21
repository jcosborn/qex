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

proc zeroComparisonGrad(zero: Gvalue, i: int): Gvalue =
  case i
  of 0, 1:
    zero
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc affineUpstream(zb: Gvalue,
                    scale: float,
                    anchor: Gvalue): Gvalue =
  if zb == nil:
    return scalarLeafLike(anchor, scale)
  if scale == 1.0:
    return zb
  if scale == -1.0:
    return -zb
  scalarLeafLike(anchor, scale) * zb

proc scalarLiteral(anchor: Gvalue,
                   value: float): Gvalue =
  scalarLeafLike(anchor, value)

proc scalarLiteral(anchor: Gvalue,
                   value: int): Gvalue =
  scalarLeafLike(anchor, value)

proc numericLiteral(anchor: Gvalue,
                    value: int): Gvalue =
  numericLeafLike(anchor, value)

proc numericLiteral(anchor: Gvalue,
                    value: float): Gvalue =
  numericLeafLike(anchor, value)

proc `+`*(x: Gvalue, y: float): Gvalue =
  x + scalarLiteral(x, y)

proc `+`*(x: float, y: Gvalue): Gvalue =
  scalarLiteral(y, x) + y

proc `+`*(x: Gvalue, y: int): Gvalue =
  x + scalarLiteral(x, y)

proc `+`*(x: int, y: Gvalue): Gvalue =
  scalarLiteral(y, x) + y

proc `-`*(x: Gvalue, y: float): Gvalue =
  x - scalarLiteral(x, y)

proc `-`*(x: float, y: Gvalue): Gvalue =
  scalarLiteral(y, x) - y

proc `-`*(x: Gvalue, y: int): Gvalue =
  x - scalarLiteral(x, y)

proc `-`*(x: int, y: Gvalue): Gvalue =
  scalarLiteral(y, x) - y

proc `*`*(x: Gvalue, y: float): Gvalue =
  x * scalarLiteral(x, y)

proc `*`*(x: float, y: Gvalue): Gvalue =
  scalarLiteral(y, x) * y

proc `*`*(x: Gvalue, y: int): Gvalue =
  x * scalarLiteral(x, y)

proc `*`*(x: int, y: Gvalue): Gvalue =
  scalarLiteral(y, x) * y

proc `/`*(x: Gvalue, y: float): Gvalue =
  x / scalarLiteral(x, y)

proc `/`*(x: float, y: Gvalue): Gvalue =
  scalarLiteral(y, x) / y

proc `/`*(x: Gvalue, y: int): Gvalue =
  x / scalarLiteral(x, y)

proc `/`*(x: int, y: Gvalue): Gvalue =
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
  let view = z.requireBinaryNodeView("* backward")
  discard dep
  binaryBackwardCase("* backward", i,
    block:
      return scaledUpstreamOr(zb, view.y),
    block:
      return scaledUpstreamOr(zb, view.x))
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
  let view = z.requireBinaryNodeView("/ backward")
  discard dep
  # d(x / y) = dx / y - x * dy / y^2
  binaryBackwardCase("/ backward", i,
    block:
      return scaledUpstreamOr(zb, 1.0 / view.y),
    block:
      return scaledUpstreamOr(zb, -z / view.y))
defineBinaryGraphOp(gsdiv, `/`, Gscalar, Gscalar, x, y, scalarNodeLike(x), divsf, divsb, "/")

defineScalarUnaryForward(expsf, "exp forward"):
  z.getfloat = math.exp(x.getfloat)
proc expsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z.requireUnaryNodeView("exp backward")
  discard dep
  unaryBackwardCase("exp backward", i):
    return scaledUpstreamOr(zb, z)
defineUnaryGraphOp(exps, exp, Gscalar, x, scalarNodeLike(x), expsf, expsb, "exps")

method `<`*(x: Gvalue, y: Gvalue): Gvalue {.base.} =
  raiseErrorBaseMethod("`<`(" & $x & ", " & $y & ")")
method equal*(x: Gvalue, y: Gvalue): Gvalue {.base.} =
  raiseErrorBaseMethod("equal(" & $x & ", " & $y & ")")

proc `not`*(x: Gvalue): Gvalue = cond(x, x.constLike(0), x.constLike(1))
proc `and`*(x: Gvalue, y: Gvalue): Gvalue = cond(x, y, y.constLike(0))
proc `or`*(x: Gvalue, y: Gvalue): Gvalue = cond(x, y.constLike(1), y)
proc `xor`*(x: Gvalue, y: Gvalue): Gvalue = cond(x, not(y), y)

proc cond*(c: Gvalue,
           x: Gvalue,
           y: float): Gvalue =
  cond(c, x, numericLiteral(if x != nil: x else: c, y))

proc cond*(c: Gvalue,
           x: float,
           y: Gvalue): Gvalue =
  cond(c, numericLiteral(if y != nil: y else: c, x), y)

proc cond*(c: Gvalue,
           x: Gvalue,
           y: int): Gvalue =
  cond(c, x, numericLiteral(if x != nil: x else: c, y))

proc cond*(c: Gvalue,
           x: int,
           y: Gvalue): Gvalue =
  cond(c, numericLiteral(if y != nil: y else: c, x), y)

proc `>`*(x, y: Gvalue): Gvalue = y < x
proc `>=`*(x, y: Gvalue): Gvalue = (y < x) or equal(x, y)
proc `<=`*(x, y: Gvalue): Gvalue = (x < y) or equal(x, y)

proc `>`*(x, y: Gscalar): Gvalue = `>`(Gvalue(x), Gvalue(y))
proc `>=`*(x, y: Gscalar): Gvalue = `>=`(Gvalue(x), Gvalue(y))
proc `<=`*(x, y: Gscalar): Gvalue = `<=`(Gvalue(x), Gvalue(y))

proc `>`*(x, y: Gint): Gvalue = `>`(Gvalue(x), Gvalue(y))
proc `>=`*(x, y: Gint): Gvalue = `>=`(Gvalue(x), Gvalue(y))
proc `<=`*(x, y: Gint): Gvalue = `<=`(Gvalue(x), Gvalue(y))

proc `<`*(x: Gvalue, y: float): Gvalue = x < numericLiteral(x, y)
proc `<`*(x: float, y: Gvalue): Gvalue = numericLiteral(y, x) < y
proc `<`*(x: Gvalue, y: int): Gvalue = x < numericLiteral(x, y)
proc `<`*(x: int, y: Gvalue): Gvalue = numericLiteral(y, x) < y

proc equal*(x: Gvalue, y: float): Gvalue = equal(x, numericLiteral(x, y))
proc equal*(x: float, y: Gvalue): Gvalue = equal(numericLiteral(y, x), y)
proc equal*(x: Gvalue, y: int): Gvalue = equal(x, numericLiteral(x, y))
proc equal*(x: int, y: Gvalue): Gvalue = equal(numericLiteral(y, x), y)

proc `>`*(x: Gvalue, y: float): Gvalue = y < x
proc `>`*(x: float, y: Gvalue): Gvalue = y < x
proc `>`*(x: Gvalue, y: int): Gvalue = y < x
proc `>`*(x: int, y: Gvalue): Gvalue = y < x

proc `>=`*(x: Gvalue, y: float): Gvalue = (y < x) or equal(x, y)
proc `>=`*(x: float, y: Gvalue): Gvalue = (y < x) or equal(x, y)
proc `>=`*(x: Gvalue, y: int): Gvalue = (y < x) or equal(x, y)
proc `>=`*(x: int, y: Gvalue): Gvalue = (y < x) or equal(x, y)

proc `<=`*(x: Gvalue, y: float): Gvalue = (x < y) or equal(x, y)
proc `<=`*(x: float, y: Gvalue): Gvalue = (x < y) or equal(x, y)
proc `<=`*(x: Gvalue, y: int): Gvalue = (x < y) or equal(x, y)
proc `<=`*(x: int, y: Gvalue): Gvalue = (x < y) or equal(x, y)

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
