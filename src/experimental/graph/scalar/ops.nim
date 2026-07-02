import ../core
import types
import ../support/op
import math

type
  CondSelector = Gscalar | Gint
  ScalarLiteral = int | float
proc `-`*(x: Gscalar): Gscalar
proc `*`*(x: Gscalar, y: Gscalar): Gscalar
proc `/`*(x: Gscalar, y: Gscalar): Gscalar

proc affineUpstream(zb: Gvalue,
                    scale: float,
                    anchor: Gvalue): Gscalar =
  if zb == nil:
    return toGvalue(anchor.runtime, scale)
  let zb = Gscalar(zb)
  if scale == 1.0:
    return zb
  if scale == -1.0:
    return -zb
  toGvalue(anchor.runtime, scale) * zb

proc negsf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let z = Gscalar(v)
  z.sval = -x.sval

proc negsb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  affineUpstream(zb, -1.0, z)

let gsneg = Gfunc(forward: negsf, backward: negsb, name: "-")

proc `-`*(x: Gscalar): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x)], gsneg, "-")

proc addsf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let y = Gscalar(v.inputs[1])
  let z = Gscalar(v)
  z.sval = x.sval + y.sval

proc addsb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  affineUpstream(zb, 1.0, z)

let gsadd = Gfunc(forward: addsf, backward: addsb, name: "+")

proc `+`*(x: Gscalar, y: Gscalar): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x), Gvalue(y)], gsadd, "+")

proc mulsf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let y = Gscalar(v.inputs[1])
  let z = Gscalar(v)
  z.sval = x.sval * y.sval

proc mulsb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  bilinearBackward(zb, z, i, Gscalar)

let gsmul = Gfunc(forward: mulsf, backward: mulsb, name: "*")

proc `*`*(x: Gscalar, y: Gscalar): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x), Gvalue(y)], gsmul, "*")

# Real dot product over scalars is just the product. Reached through generic
# contraction/test helpers (e.g. `redot(grad(ff, x), a)`), not a direct caller.
proc redot*(x: Gscalar, y: Gscalar): Gscalar = x*y

proc subsf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let y = Gscalar(v.inputs[1])
  let z = Gscalar(v)
  z.sval = x.sval - y.sval

proc subsb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  if i == 0:
    return affineUpstream(zb, 1.0, z)
  affineUpstream(zb, -1.0, z)

let gssub = Gfunc(forward: subsf, backward: subsb, name: "-")

proc `-`*(x: Gscalar, y: Gscalar): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x), Gvalue(y)], gssub, "-")

proc divsf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let y = Gscalar(v.inputs[1])
  let z = Gscalar(v)
  z.sval = x.sval / y.sval

proc divsb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let y = Gscalar(z.inputs[1])
  # d(x / y) = dx / y - x * dy / y^2
  if i == 0:
    return scaledUpstreamOr(
      zb,
      Gscalar,
      toGvalue(y.runtime, 1.0) / y)
  scaledUpstreamOr(
    zb,
    Gscalar,
    -Gscalar(z) / y)

let gsdiv = Gfunc(forward: divsf, backward: divsb, name: "/")

proc `/`*(x: Gscalar, y: Gscalar): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x), Gvalue(y)], gsdiv, "/")

proc expsf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let z = Gscalar(v)
  z.sval = math.exp(x.sval)

proc expsb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  scaledUpstreamOr(zb, Gscalar, Gscalar(z))

let exps = Gfunc(forward: expsf, backward: expsb, name: "exps")

proc exp*(x: Gscalar): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x)], exps, "exps")

proc lazyScalef(v: Gvalue) =
  let upstream = Gscalar(v.inputs[0])
  let contribution = Gscalar(v.inputs[1])
  let z = Gscalar(v)
  if upstream.isZero:
    z.sval = 0.0
  else:
    z.sval = upstream.sval * contribution.sval

proc lazyScaleInputView(v: Gvalue,
                        mode: InputWalkMode,
                        visit: GnodeVisit) =
  let upstream = Gscalar(v.inputs[0])
  visit upstream
  if mode != iwmEval or not upstream.isZero:
    visit v.inputs[1]

proc lazyScaleb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let upstream = Gscalar(z.inputs[0])
  let contribution = Gscalar(z.inputs[1])
  if i == 0:
    return scaledUpstreamOr(zb, Gscalar, contribution)
  scaledUpstreamOr(zb, Gscalar, upstream)

let glazyScale = Gfunc(
  forward: lazyScalef,
  inputView: lazyScaleInputView,
  backward: lazyScaleb,
  name: "scalarScale")

proc lazyScale(upstream: Gscalar, contribution: Gscalar): Gscalar =
  graphNode(
    scalarNodeLike(contribution),
    @[Gvalue(upstream), Gvalue(contribution)],
    glazyScale,
    "scalarScale")

method addLike*(prototype: Gscalar, x: Gvalue, y: Gvalue): Gvalue =
  Gscalar(x) + Gscalar(y)

method scaleLike*(contribution: Gscalar, upstream: Gvalue): Gvalue =
  lazyScale(Gscalar(upstream), contribution)

method addLike*(prototype: Gint, x: Gvalue, y: Gvalue): Gvalue =
  let left = Gint(x)
  let right = Gint(y)
  if left.isZero:
    return right
  if right.isZero:
    return left
  raiseValueError("cannot accumulate non-zero int graph gradients")

proc cond*[C: CondSelector, T: Gvalue](c: C, x: T, y: T): T =
  newCondNode(c, x, y)

proc cond*[T: ScalarLiteral](c: CondSelector, x: Gscalar, y: T): Gscalar =
  cond(c, x, toGvalue(x.runtime, float(y)))

proc cond*[T: ScalarLiteral](c: CondSelector, x: T, y: Gscalar): Gscalar =
  cond(c, toGvalue(y.runtime, float(x)), y)

proc cond*(c: CondSelector, x: Gint, y: int): Gint =
  cond(c, x, toGvalue(x.runtime, y))

proc cond*(c: CondSelector, x: int, y: Gint): Gint =
  cond(c, toGvalue(y.runtime, x), y)

proc falseValue(x: Gscalar): Gscalar =
  toGvalue(x.runtime, 0.0)

proc trueValue(x: Gscalar): Gscalar =
  toGvalue(x.runtime, 1.0)

proc falseValue(x: Gint): Gint =
  toGvalue(x.runtime, 0)

proc trueValue(x: Gint): Gint =
  toGvalue(x.runtime, 1)

proc `not`*[T: CondSelector](x: T): T =
  cond(x, x.falseValue, x.trueValue)

proc `and`*[C: CondSelector, T: CondSelector](x: C, y: T): T =
  cond(x, y, y.falseValue)

proc `or`*[C: CondSelector, T: CondSelector](x: C, y: T): T =
  cond(x, y.trueValue, y)

proc `xor`*[C: CondSelector, T: CondSelector](x: C, y: T): T =
  cond(x, not(y), y)

# Comparisons are piecewise constant: the gradient w.r.t. either operand is zero.
# One shared backward serves every scalar/int comparison; `numericLeafLike` picks
# the zero-leaf type (Gscalar or Gint) from the result node `z`.
proc comparisonZeroBackward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  numericLeafLike(z, 0)

proc ltsf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let y = Gscalar(v.inputs[1])
  let z = Gscalar(v)
  z.sval = if x.sval < y.sval: 1.0 else: 0.0

let lts = Gfunc(forward: ltsf, backward: comparisonZeroBackward, name: "lts")

proc `<`*(x: Gscalar, y: Gscalar): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x), Gvalue(y)], lts, "lts")

proc equalsf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let y = Gscalar(v.inputs[1])
  let z = Gscalar(v)
  z.sval = if x.sval == y.sval: 1.0 else: 0.0

let equals = Gfunc(forward: equalsf, backward: comparisonZeroBackward, name: "equals")

proc equal*(x: Gscalar, y: Gscalar): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x), Gvalue(y)], equals, "equals")

proc ltif(v: Gvalue) =
  let x = Gint(v.inputs[0])
  let y = Gint(v.inputs[1])
  let z = Gint(v)
  z.ival = if x.ival < y.ival: 1 else: 0

let lti = Gfunc(forward: ltif, backward: comparisonZeroBackward, name: "lti")

proc `<`*(x: Gint, y: Gint): Gint =
  graphNode(intNodeLike(x), @[Gvalue(x), Gvalue(y)], lti, "lti")

proc equalif(v: Gvalue) =
  let x = Gint(v.inputs[0])
  let y = Gint(v.inputs[1])
  let z = Gint(v)
  z.ival = if x.ival == y.ival: 1 else: 0

let equali = Gfunc(forward: equalif, backward: comparisonZeroBackward, name: "equali")

proc equal*(x: Gint, y: Gint): Gint =
  graphNode(intNodeLike(x), @[Gvalue(x), Gvalue(y)], equali, "equali")

# Note: nim system defines `>` and `>=` as templates.

# Concrete per-type overloads (not a `CondSelector` generic): a generic here is
# ambiguous with system's `<=`/`>=`/`>` over `ref T`, since Gscalar/Gint are refs.
proc `<=`*(x, y: Gscalar): auto = (x < y) or equal(x, y)
proc `<=`*(x, y: Gint): auto = (x < y) or equal(x, y)

# Literal RHS/LHS sugar for binary ops: each op already exists value-vs-value
# (above), so one template pair emits both literal overloads uniformly. Covers
# arithmetic (+, -, *, /) and comparisons; the literal is anchored to the operand.
template litCmp(op) {.dirty.} =
  proc op*[T: ScalarLiteral](x: Gscalar, y: T): Gscalar =
    op(x, toGvalue(x.runtime, float(y)))
  proc op*[T: ScalarLiteral](x: T, y: Gscalar): Gscalar =
    op(toGvalue(y.runtime, float(x)), y)
template litCmpInt(op) {.dirty.} =
  proc op*(x: Gint, y: int): Gint = op(x, toGvalue(x.runtime, y))
  proc op*(x: int, y: Gint): Gint = op(toGvalue(y.runtime, x), y)

litCmp(`+`)
litCmp(`-`)
litCmp(`*`)
litCmp(`/`)
litCmp(`<`)
litCmp(equal)
litCmp(`<=`)
litCmpInt(`<`)
litCmpInt(equal)
litCmpInt(`<=`)
