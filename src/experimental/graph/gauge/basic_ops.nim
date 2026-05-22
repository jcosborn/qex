import ../[core, scalar]
import ../scalar/types
import ../support/op
import layout, gauge, physics/qcdTypes
import shared

# Section: Basic Gauge Ops

proc retr*(x: Ggauge): Gscalar
proc adj*(x: Ggauge): Ggauge
proc norm2*(x: Ggauge): Gscalar
proc redot*(x: Gscalar, y: Gscalar): Gscalar
proc redot*(x: Ggauge, y: Ggauge): Gscalar
proc exp*(x: Ggauge): Ggauge
proc expDeriv*(b: Ggauge, x: Ggauge): Ggauge
proc projTAH*(x: Ggauge): Ggauge
proc `-`*(x: Ggauge): Ggauge
proc `+`*(x: Gscalar, y: Ggauge): Ggauge
proc `+`*(x: Ggauge, y: Ggauge): Ggauge
proc `*`*(x: Gscalar, y: Ggauge): Ggauge
proc `*`*(x: Ggauge, y: Ggauge): Ggauge
proc `-`*(x: Ggauge, y: Gscalar): Ggauge
proc `-`*(x: Ggauge, y: Ggauge): Ggauge

proc redot*(x: Gscalar, y: Gscalar): Gscalar = x*y

proc `+`*(x: float, y: Ggauge): Ggauge =
  scalarLeafLike(y, x) + y

proc `+`*(x: int, y: Ggauge): Ggauge =
  scalarLeafLike(y, x) + y

proc `*`*(x: float, y: Ggauge): Ggauge =
  scalarLeafLike(y, x) * y

proc `*`*(x: int, y: Ggauge): Ggauge =
  scalarLeafLike(y, x) * y

proc `-`*(x: Ggauge, y: float): Ggauge =
  x - scalarLeafLike(x, y)

proc `-`*(x: Ggauge, y: int): Ggauge =
  x - scalarLeafLike(x, y)

proc retrgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireUnaryNodeView(Ggauge, "retr backward")
  discard dep
  unaryBackwardCase("retr backward", i):
    let g = view.x.unsafeGaugeStorage.newOneOf
    threads:
      for f in g:
        f := 1.0
    return scaledUpstreamOr(zb, Gscalar, toGvalue(view.x.runtime, g), "retr backward")

proc retrgf(v: Gvalue) =
  let view = v.requireUnaryNodeView(Ggauge, "retr forward")
  let x = view.x
  let z = v.requireScalar("retr forward result")
  threads:
    var t = 0.0
    for mu in 0..<x.unsafeGaugeStorage.len:
      t += x.unsafeGaugeStorage[mu].trace.re
    threadMaster: z.sval = t

let retrg = newGfunc(forward = retrgf, backward = retrgb, name = "retrg")

proc retr*(x: Ggauge): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x)], retrg, "retrg")

proc adjgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z
  discard dep
  unaryBackwardCase("adj backward", i):
    return requireUpstream(zb, Ggauge, "adj backward").adj

proc adjgf(v: Gvalue) =
  let view = v.requireUnaryNodeView(Ggauge, "adj forward")
  let x = view.x
  let z = v.requireGauge("adj forward result")
  z.mapGaugeSites(x.unsafeGaugeStorage[mu].adj)

let adjg = newGfunc(forward = adjgf, backward = adjgb, name = "adjg")

proc adj*(x: Ggauge): Ggauge =
  graphNode(x.gaugeNodeLike, @[Gvalue(x)], adjg, "adjg")

proc norm2gb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireUnaryNodeView(Ggauge, "norm2 backward")
  discard dep
  unaryBackwardCase("norm2 backward", i):
    return scaledUpstreamOr(zb, Gscalar, 2.0 * view.x, "norm2 backward")

proc norm2gf(v: Gvalue) =
  let view = v.requireUnaryNodeView(Ggauge, "norm2 forward")
  let x = view.x
  let z = v.requireScalar("norm2 forward result")
  threads:
    var t = 0.0
    for mu in 0..<x.unsafeGaugeStorage.len:
      t += x.unsafeGaugeStorage[mu].norm2
    threadMaster: z.sval = t

let norm2g = newGfunc(forward = norm2gf, backward = norm2gb, name = "norm2g")

proc norm2*(x: Ggauge): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x)], norm2g, "norm2g")

proc neggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z
  discard dep
  unaryBackwardCase("neg backward", i):
    return -requireUpstream(zb, Ggauge, "neg backward")

proc neggf(v: Gvalue) =
  let view = v.requireUnaryNodeView(Ggauge, "neg forward")
  let x = view.x
  let z = v.requireGauge("neg forward result")
  z.mapGaugeSites(-x.unsafeGaugeStorage[mu])

let negg = newGfunc(forward = neggf, backward = neggb, name = "-g")

proc `-`*(x: Ggauge): Ggauge =
  graphNode(x.gaugeNodeLike, @[Gvalue(x)], negg, "-g")

proc addsgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z
  discard dep
  binaryBackwardCase("s+g backward", i,
    block:
      return retr(requireUpstream(zb, Ggauge, "s+g backward")),
    block:
      return requireUpstream(zb, Ggauge, "s+g backward"))

proc addsgf(v: Gvalue) =
  let view = v.requireBinaryNodeView(Gscalar, Ggauge, "s+g forward")
  let x = view.x
  let y = view.y
  let z = v.requireGauge("s+g forward result")
  z.mapGaugeSites(x.sval + y.unsafeGaugeStorage[mu])

let addsg = newGfunc(forward = addsgf, backward = addsgb, name = "s+g")

proc `+`*(x: Gscalar, y: Ggauge): Ggauge =
  graphNode(y.gaugeNodeLike, @[Gvalue(x), Gvalue(y)], addsg, "s+g")

proc addggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z
  discard dep
  binaryBackwardCase("g+g backward", i,
    block:
      return requireUpstream(zb, Ggauge, "g+g backward"),
    block:
      return requireUpstream(zb, Ggauge, "g+g backward"))

proc addggf(v: Gvalue) =
  let view = v.requireBinaryNodeView(Ggauge, Ggauge, "g+g forward")
  let x = view.x
  let y = view.y
  let z = v.requireGauge("g+g forward result")
  z.mapGaugeSites(x.unsafeGaugeStorage[mu] + y.unsafeGaugeStorage[mu])

let addgg = newGfunc(forward = addggf, backward = addggb, name = "g+g")

proc `+`*(x: Ggauge, y: Ggauge): Ggauge =
  graphNode(sameShapeGaugeNodeLike(x, y, "g+g"), @[Gvalue(x), Gvalue(y)], addgg, "g+g")

method addLike*(prototype: Ggauge, x: Gvalue, y: Gvalue): Gvalue =
  discard prototype
  x.requireGauge("gauge gradient add left") +
    y.requireGauge("gauge gradient add right")

proc mulsgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView(Gscalar, Ggauge, "s*g backward")
  discard dep
  binaryBackwardCase("s*g backward", i,
    block:
      return redot(requireUpstream(zb, Ggauge, "s*g backward"), view.y),
    block:
      return view.x * requireUpstream(zb, Ggauge, "s*g backward"))

proc mulsgf(v: Gvalue) =
  let view = v.requireBinaryNodeView(Gscalar, Ggauge, "s*g forward")
  let x = view.x
  let y = view.y
  let z = v.requireGauge("s*g forward result")
  z.mapGaugeSites(x.sval * y.unsafeGaugeStorage[mu])

let mulsg = newGfunc(forward = mulsgf, backward = mulsgb, name = "s*g")

proc `*`*(x: Gscalar, y: Ggauge): Ggauge =
  graphNode(y.gaugeNodeLike, @[Gvalue(x), Gvalue(y)], mulsg, "s*g")

proc mulggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView(Ggauge, Ggauge, "g*g backward")
  discard dep
  binaryBackwardCase("g*g backward", i,
    block:
      return requireUpstream(zb, Ggauge, "g*g backward") * view.y.adj,
    block:
      return view.x.adj * requireUpstream(zb, Ggauge, "g*g backward"))

proc mulggf(v: Gvalue) =
  let view = v.requireBinaryNodeView(Ggauge, Ggauge, "g*g forward")
  let x = view.x
  let y = view.y
  let z = v.requireGauge("g*g forward result")
  z.mapGaugeSites(x.unsafeGaugeStorage[mu] * y.unsafeGaugeStorage[mu])

let mulgg = newGfunc(forward = mulggf, backward = mulggb, name = "g*g")

proc `*`*(x: Ggauge, y: Ggauge): Ggauge =
  graphNode(sameShapeGaugeNodeLike(x, y, "g*g"), @[Gvalue(x), Gvalue(y)], mulgg, "g*g")

method scaleLike*(contribution: Ggauge, upstream: Gvalue): Gvalue =
  if upstream of Gscalar:
    return Gscalar(upstream) * contribution
  if upstream of Ggauge:
    return Ggauge(upstream) * contribution
  raiseValueError("gauge scale upstream expects scalar or gauge value, got:\n" & upstream.nodeRepr)

proc redotggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard dep
  let view = z.requireBinaryNodeView(Ggauge, Ggauge, "redotgg backward")
  binaryBackwardCase("redotgg backward", i,
    block:
      return scaledUpstreamOr(zb, Gscalar, view.y, "redotgg backward"),
    block:
      return scaledUpstreamOr(zb, Gscalar, view.x, "redotgg backward"))

proc redotggf(v: Gvalue) =
  let view = v.requireBinaryNodeView(Ggauge, Ggauge, "redot forward")
  let x = view.x
  let y = view.y
  let z = v.requireScalar("redot forward result")
  threads:
    var t = 0.0
    for mu in 0..<x.unsafeGaugeStorage.len:
      t += redot(x.unsafeGaugeStorage[mu], y.unsafeGaugeStorage[mu])
    threadMaster: z.sval = t

let redotgg = newGfunc(forward = redotggf, backward = redotggb, name = "redotgg")

proc redot*(x: Ggauge, y: Ggauge): Gscalar =
  x.requireSameGaugeShape(y, "redotgg")
  graphNode(scalarNodeLike(x), @[Gvalue(x), Gvalue(y)], redotgg, "redotgg")

proc subgsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z
  discard dep
  binaryBackwardCase("g-s backward", i,
    block:
      return requireUpstream(zb, Ggauge, "g-s backward"),
    block:
      return -retr(requireUpstream(zb, Ggauge, "g-s backward")))

proc subgsf(v: Gvalue) =
  let view = v.requireBinaryNodeView(Ggauge, Gscalar, "g-s forward")
  let x = view.x
  let y = view.y
  let z = v.requireGauge("g-s forward result")
  z.mapGaugeSites(x.unsafeGaugeStorage[mu] - y.sval)

let subgs = newGfunc(forward = subgsf, backward = subgsb, name = "g-s")

proc `-`*(x: Ggauge, y: Gscalar): Ggauge =
  graphNode(x.gaugeNodeLike, @[Gvalue(x), Gvalue(y)], subgs, "g-s")

proc subggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z
  discard dep
  binaryBackwardCase("g-g backward", i,
    block:
      return requireUpstream(zb, Ggauge, "g-g backward"),
    block:
      return -requireUpstream(zb, Ggauge, "g-g backward"))

proc subggf(v: Gvalue) =
  let view = v.requireBinaryNodeView(Ggauge, Ggauge, "g-g forward")
  let x = view.x
  let y = view.y
  let z = v.requireGauge("g-g forward result")
  z.mapGaugeSites(x.unsafeGaugeStorage[mu] - y.unsafeGaugeStorage[mu])

let subgg = newGfunc(forward = subggf, backward = subggb, name = "g-g")

proc `-`*(x: Ggauge, y: Ggauge): Ggauge =
  graphNode(sameShapeGaugeNodeLike(x, y, "g-g"), @[Gvalue(x), Gvalue(y)], subgg, "g-g")

proc expgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireUnaryNodeView(Ggauge, "exp backward")
  discard dep
  unaryBackwardCase("exp backward", i):
    return expDeriv(requireUpstream(zb, Ggauge, "exp backward"), view.x)

proc expgf(v: Gvalue) =
  let view = v.requireUnaryNodeView(Ggauge, "exp forward")
  let x = view.x
  let z = v.requireGauge("exp forward result")
  z.mapGaugeElements:
    z.unsafeGaugeStorage[mu][e] := exp(x.unsafeGaugeStorage[mu][e])

let expg = newGfunc(forward = expgf, backward = expgb, name = "expg")

proc exp*(x: Ggauge): Ggauge =
  graphNode(x.gaugeNodeLike, @[Gvalue(x)], expg, "expg")

proc expDerivgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard zb
  discard z
  discard dep
  binaryBackwardCase("expDeriv backward", i,
    block:
      raiseUnsupportedPath("expDeriv backward", "derivative with respect to force-direction input"),
    block:
      raiseUnsupportedPath("expDeriv backward", "second derivatives of matrix exponential"))

proc expDerivgf(v: Gvalue) =
  let view = v.requireBinaryNodeView(Ggauge, Ggauge, "expDeriv forward")
  let x = view.x
  let y = view.y
  let z = v.requireGauge("expDeriv forward result")
  z.mapGaugeElements:
    z.unsafeGaugeStorage[mu][e] := expDeriv(
      y.unsafeGaugeStorage[mu][e],
      x.unsafeGaugeStorage[mu][e])

let expDerivg = newGfunc(forward = expDerivgf, backward = expDerivgb, name = "expDerivg")

proc expDeriv*(b: Ggauge, x: Ggauge): Ggauge =
  graphNode(
    sameShapeGaugeNodeLike(b, x, "expDerivg"),
    @[Gvalue(b), Gvalue(x)],
    expDerivg,
    "expDerivg")

proc projTAHb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z
  discard dep
  unaryBackwardCase("projTAH backward", i):
    return projTAH(requireUpstream(zb, Ggauge, "projTAH backward"))

proc projTAHf(v: Gvalue) =
  let view = v.requireUnaryNodeView(Ggauge, "projTAH forward")
  let x = view.x
  let z = v.requireGauge("projTAH forward result")
  z.mapGaugeElements:
    z.unsafeGaugeStorage[mu][e].projectTAH(x.unsafeGaugeStorage[mu][e])

let projTAHg = newGfunc(forward = projTAHf, backward = projTAHb, name = "projTAH")

proc projTAH*(x: Ggauge): Ggauge =
  graphNode(x.gaugeNodeLike, @[Gvalue(x)], projTAHg, "projTAH")
