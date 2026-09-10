import ../[core, scalar]
import ../support/op
import layout, gauge, physics/qcdTypes
from maths/matexp import newExpParam, ekPoly
import shared

# Section: Basic Gauge Ops

type GaugeLiteral = int | float

proc retr*(x: Ggauge): Gscalar
proc adj*(x: Ggauge): Ggauge
proc norm2*(x: Ggauge): Gscalar
proc redot*(x: Ggauge, y: Ggauge): Gscalar
proc exp*(x: Ggauge): Ggauge
proc expDeriv*(b: Ggauge, x: Ggauge, parity = -1, dir = 0): Ggauge
proc projTAH*(x: Ggauge): Ggauge
proc `-`*(x: Ggauge): Ggauge
proc `+`*(x: Gscalar, y: Ggauge): Ggauge
proc `+`*(x: Ggauge, y: Ggauge): Ggauge
proc `*`*(x: Gscalar, y: Ggauge): Ggauge
proc `*`*(x: Ggauge, y: Ggauge): Ggauge
proc `-`*(x: Ggauge, y: Gscalar): Ggauge
proc `-`*(x: Ggauge, y: Ggauge): Ggauge

proc `+`*[T: GaugeLiteral](x: T, y: Ggauge): Ggauge =
  toGvalue(y.runtime, float(x)) + y

proc `*`*[T: GaugeLiteral](x: T, y: Ggauge): Ggauge =
  toGvalue(y.runtime, float(x)) * y

proc `-`*[T: GaugeLiteral](x: Ggauge, y: T): Ggauge =
  x - toGvalue(x.runtime, float(y))

proc retrgb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  scaledUpstreamOr(zb, Gscalar, Ggauge(z.inputs[0]).unitGaugeLike)

proc retrgf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let z = Gscalar(v)
  threads:
    var t = 0.0
    for mu in 0..<x.gval.len:
      t += x.gval[mu].trace.re
    threadMaster: z.sval = t

let retrg = Gfunc(forward: retrgf, backward: retrgb, name: "retrg")

proc retr*(x: Ggauge): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x)], retrg, "retrg")

proc adjgb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  requireUpstream(zb, "adjg backward", Ggauge).adj

proc adjgf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let z = Ggauge(v)
  z.mapGaugeSites(x.gval[mu].adj)

let adjg = Gfunc(forward: adjgf, backward: adjgb, name: "adjg")

proc adj*(x: Ggauge): Ggauge =
  graphNode(x.gaugeNodeLike, @[Gvalue(x)], adjg, "adjg")

proc norm2gb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let x = Ggauge(z.inputs[0])
  scaledUpstreamOr(zb, Gscalar, toGvalue(x.runtime, 2.0) * x)

proc norm2gf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let z = Gscalar(v)
  threads:
    var t = 0.0
    for mu in 0..<x.gval.len:
      t += x.gval[mu].norm2
    threadMaster: z.sval = t

let norm2g = Gfunc(forward: norm2gf, backward: norm2gb, name: "norm2g")

proc norm2*(x: Ggauge): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x)], norm2g, "norm2g")

proc neggb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  -requireUpstream(zb, "-g backward", Ggauge)

proc neggf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let z = Ggauge(v)
  z.mapGaugeSites(-x.gval[mu])

let negg = Gfunc(forward: neggf, backward: neggb, name: "-g")

proc `-`*(x: Ggauge): Ggauge =
  graphNode(x.gaugeNodeLike, @[Gvalue(x)], negg, "-g")

proc addsgb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let upstream = requireUpstream(zb, "s+g backward", Ggauge)
  if i == 0:
    return retr(upstream)
  upstream

proc addsgf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let y = Ggauge(v.inputs[1])
  let z = Ggauge(v)
  z.mapGaugeSites(x.sval + y.gval[mu])

let addsg = Gfunc(forward: addsgf, backward: addsgb, name: "s+g")

proc `+`*(x: Gscalar, y: Ggauge): Ggauge =
  graphNode(y.gaugeNodeLike, @[Gvalue(x), Gvalue(y)], addsg, "s+g")

proc addggb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  requireUpstream(zb, "g+g backward", Ggauge)

proc addggf(v: Gvalue) =
  let z = Ggauge(v)
  threads:
    for mu in 0..<z.gval.len:
      for e in z.gval[mu]:
        var t {.noinit.}: evalType(z.gval[mu][e])
        t := Ggauge(v.inputs[0]).gval[mu][e]
        for i in 1..<v.inputs.len:
          t += Ggauge(v.inputs[i]).gval[mu][e]
        z.gval[mu][e] := t

let addgg = Gfunc(forward: addggf, backward: addggb, name: "g+g")

proc gaugeAddTerms*(x: Ggauge): seq[Ggauge] =
  ## Return the leaves of a gauge sum in evaluation order.
  if x.gfunc == addgg:
    for input in x.inputs:
      result.add gaugeAddTerms(Ggauge(input))
  else:
    result.add x

proc `+`*(x: Ggauge, y: Ggauge): Ggauge =
  x.requireSameGaugeShape(y, "g+g")
  var inputs: seq[Gvalue]
  if x.gfunc == addgg:
    inputs.add x.inputs
  else:
    inputs.add Gvalue(x)
  inputs.add Gvalue(y)
  graphNode(x.gaugeNodeLike, inputs, addgg, "g+g")

method addLike*(prototype: Ggauge, x: Gvalue, y: Gvalue): Gvalue =
  Ggauge(x) + Ggauge(y)

proc mulsgb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let x = Gscalar(z.inputs[0])
  let y = Ggauge(z.inputs[1])
  let upstream = requireUpstream(zb, "s*g backward", Ggauge)
  if i == 0:
    return redot(upstream, y)
  x * upstream

proc mulsgf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let y = Ggauge(v.inputs[1])
  let z = Ggauge(v)
  z.mapGaugeSites(x.sval * y.gval[mu])

let mulsg = Gfunc(forward: mulsgf, backward: mulsgb, name: "s*g")

proc `*`*(x: Gscalar, y: Ggauge): Ggauge =
  graphNode(y.gaugeNodeLike, @[Gvalue(x), Gvalue(y)], mulsg, "s*g")

proc mulggb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let x = Ggauge(z.inputs[0])
  let y = Ggauge(z.inputs[1])
  let upstream = requireUpstream(zb, "g*g backward", Ggauge)
  if i == 0:
    return upstream * y.adj
  x.adj * upstream

proc mulggf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let y = Ggauge(v.inputs[1])
  let z = Ggauge(v)
  z.mapGaugeSites(x.gval[mu] * y.gval[mu])

let mulgg = Gfunc(forward: mulggf, backward: mulggb, name: "g*g")

proc `*`*(x: Ggauge, y: Ggauge): Ggauge =
  graphNode(sameShapeGaugeNodeLike(x, y, "g*g"), @[Gvalue(x), Gvalue(y)], mulgg, "g*g")

method scaleLike*(contribution: Ggauge, upstream: Gvalue): Gvalue =
  if upstream of Gscalar:
    return Gscalar(upstream) * contribution
  if upstream of Ggauge:
    return Ggauge(upstream) * contribution
  raiseValueError("gauge scale upstream expects scalar or gauge value, got:\n" & upstream.nodeRepr)

proc redotggb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  bilinearBackward(zb, z, i, Ggauge)

proc redotggf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let y = Ggauge(v.inputs[1])
  let z = Gscalar(v)
  threads:
    var t = 0.0
    for mu in 0..<x.gval.len:
      t += redot(x.gval[mu], y.gval[mu])
    threadMaster: z.sval = t

let redotgg = Gfunc(forward: redotggf, backward: redotggb, name: "redotgg")

proc redot*(x: Ggauge, y: Ggauge): Gscalar =
  x.requireSameGaugeShape(y, "redotgg")
  graphNode(scalarNodeLike(x), @[Gvalue(x), Gvalue(y)], redotgg, "redotgg")

proc subgsb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let upstream = requireUpstream(zb, "g-s backward", Ggauge)
  if i == 0:
    return upstream
  -retr(upstream)

proc subgsf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let y = Gscalar(v.inputs[1])
  let z = Ggauge(v)
  z.mapGaugeSites(x.gval[mu] - y.sval)

let subgs = Gfunc(forward: subgsf, backward: subgsb, name: "g-s")

proc `-`*(x: Ggauge, y: Gscalar): Ggauge =
  graphNode(x.gaugeNodeLike, @[Gvalue(x), Gvalue(y)], subgs, "g-s")

proc subggb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let upstream = requireUpstream(zb, "g-g backward", Ggauge)
  if i == 0:
    return upstream
  -upstream

proc subggf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let y = Ggauge(v.inputs[1])
  let z = Ggauge(v)
  z.mapGaugeSites(x.gval[mu] - y.gval[mu])

let subgg = Gfunc(forward: subggf, backward: subggb, name: "g-g")

proc `-`*(x: Ggauge, y: Ggauge): Ggauge =
  graphNode(sameShapeGaugeNodeLike(x, y, "g-g"), @[Gvalue(x), Gvalue(y)], subgg, "g-g")

proc blendSubset*(parity, dir: int, cand, x: Ggauge): Ggauge =
  ## Use `cand` on one parity/direction subset and `x` elsewhere.
  requireParityDir(parity, dir, x.gval.len, "blendSubset")
  let sub = x.gval.paritySubset(parity)

  proc forward(v: Gvalue) =
    let
      cand = Ggauge(v.inputs[0])
      x = Ggauge(v.inputs[1])
      z = Ggauge(v)
    threads:
      for mu in 0..<z.gval.len:
        z.gval[mu] := x.gval[mu]
      threadBarrier()
      for e in sub:
        z.gval[dir][e] := cand.gval[dir][e]

  proc backward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
    let
      up = requireUpstream(zb, "blendSubset backward", Ggauge)
      zero = Ggauge(up.zeroLike)
    if i == 0:
      Gvalue(blendSubset(parity, dir, up, zero))
    else:
      Gvalue(blendSubset(parity, dir, zero, up))

  graphNode(sameShapeGaugeNodeLike(cand, x, "blendSubset"), @[Gvalue(cand), Gvalue(x)], Gfunc(forward: forward, backward: backward, name: "blendSubset"), "blendSubset")

proc maskSubset*(parity, dir: int, x: Ggauge): Ggauge =
  ## x on the (parity, dir) subset, zero elsewhere.
  blendSubset(parity, dir, x, Ggauge(x.zeroLike))

proc expgb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let x = Ggauge(z.inputs[0])
  expDeriv(requireUpstream(zb, "expg backward", Ggauge), x)

proc expgf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let z = Ggauge(v)
  z.mapGaugeElements:
    z.gval[mu][e] := exp(x.gval[mu][e])

let expg = Gfunc(forward: expgf, backward: expgb, name: "expg")

proc exp*(x: Ggauge): Ggauge =
  graphNode(x.gaugeNodeLike, @[Gvalue(x)], expg, "expg")

proc expPolyGraph*(x: Ggauge): Ggauge =
  ## Grad-complete graph replica of the site-local matrix exponential kernel
  ## (matexp ekPoly order-4 with 2^-scale scaling and repeated squaring):
  ##   y = x/2^scale
  ##   e = (y^2/24 + y/6 + 1/2)*y^2 + y      # expm1Poly4(y)
  ##   e <- e*(e+2), scale times             # (1+e)^2 = 1 + e*(e+2)
  ##   exp(x) ~ 1 + e
  ## For Nc > 1 this is the same scheme as the optimized `exp` kernel, so the
  ## two agree to roundoff. Built only from grad-complete ops, so it can be
  ## differentiated to any order.
  const params = newExpParam()
  static:
    doAssert params.kind == ekPoly and params.order == 4,
      "expPolyGraph must match the optimized matrix exponential"
  const scale = params.scale
  let y = (1.0 / float(1 shl scale)) * x
  let y2 = y * y
  let two = toGvalue(x.runtime, 2.0)
  var e = (0.5 + ((1.0/24.0) * y2 + (1.0/6.0) * y)) * y2 + y
  for _ in 1..scale:
    e = e * (two + e)
  result = 1.0 + e

proc expDerivContribution(u: Ggauge, z: Gvalue, i: int,
                          parity = -1, dir = 0): Gvalue =
  ## Backward of z = expDeriv(b, x) (the pullback of exp at x with cotangent
  ## b) under the raw upstream u; (parity, dir) select the subset variant.
  let b = Ggauge(z.inputs[0])
  let x = Ggauge(z.inputs[1])
  if i == 0:
    # z is linear in b with adjoint the pushforward of exp at x. For Nc>1,
    # the kernel is a real-coefficient matrix polynomial P, so
    # (dP_x)^adj = dP_{x^dag}; scalar exp has the same adjoint identity.
    # Since exp is site-local, the subset kernel needs no masked upstream.
    return expDeriv(u, x.adj, parity, dir)
  # Second derivative of exp: differentiate the polynomial replica.
  # <u, mask(w)> == <mask(u), w> masks the replica upstream for the subset.
  let useed = if parity < 0: u else: maskSubset(parity, dir, u)
  const nc = x.gval[0][0].nrows
  when nc == 1:
    # The optimized Nc=1 kernel is scalar exp, whose pullback is
    # exp(x^dag)*b. Differentiate that exact expression instead of the
    # polynomial replica used by the matrix kernel.
    # This expression is already partial in x and keeps b live.
    useed.adj * exp(x.adj) * b
  else:
    # secondPullback keeps this the partial x contribution even when b is x;
    # see its doc.
    secondPullback(x, b, useed, proc(slot: Ggauge): Gvalue = expPolyGraph(slot))

proc expDerivgb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  expDerivContribution(
    requireUpstream(zb, "expDeriv backward", Ggauge), z, i)

proc expDerivgf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let y = Ggauge(v.inputs[1])
  let z = Ggauge(v)
  z.mapGaugeElements:
    z.gval[mu][e] := expDeriv(
      y.gval[mu][e],
      x.gval[mu][e])

let expDerivg = Gfunc(forward: expDerivgf, backward: expDerivgb, name: "expDerivg")

proc expDeriv*(b: Ggauge, x: Ggauge, parity = -1, dir = 0): Ggauge =
  ## D exp(x)^*[b], on the whole field or only (parity,dir); zero elsewhere.
  if parity != -1:
    requireParityDir(parity, dir, b.gval.len, "expDeriv")
  let node = sameShapeGaugeNodeLike(b, x, "expDerivg")
  if parity < 0:
    return graphNode(node, @[Gvalue(b), Gvalue(x)], expDerivg, "expDerivg")
  let sub = b.gval.paritySubset(parity)
  proc fwd(v: Gvalue) =
    let b = Ggauge(v.inputs[0])
    let x = Ggauge(v.inputs[1])
    let z = Ggauge(v)
    forGaugeSubset(sub):
      # element order matches whole-field expDerivgf: expDeriv(inputs[1], inputs[0])
      z.gval[dir][e] := expDeriv(x.gval[dir][e], b.gval[dir][e])
  proc bwd(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
    expDerivContribution(
      requireUpstream(zb, "expDeriv subset backward", Ggauge), z, i, parity, dir)
  result = graphNode(node, @[Gvalue(b), Gvalue(x)],
    Gfunc(forward: fwd, backward: bwd, name: "expDerivg"),
    "expDerivg")
  result.zeroGaugeStorage

proc projTAHb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  projTAH(requireUpstream(zb, "projTAH backward", Ggauge))

proc projTAHf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let z = Ggauge(v)
  z.mapGaugeElements:
    z.gval[mu][e].projectTAH(x.gval[mu][e])

let projTAHg = Gfunc(forward: projTAHf, backward: projTAHb, name: "projTAH")

proc projTAH*(x: Ggauge): Ggauge =
  graphNode(x.gaugeNodeLike, @[Gvalue(x)], projTAHg, "projTAH")
