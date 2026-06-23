import ../[core, scalar]
import ../support/op
import layout, gauge, physics/qcdTypes
import shared

# Section: Basic Gauge Ops

type GaugeLiteral = int | float

proc retr*(x: Ggauge): Gscalar
proc adj*(x: Ggauge): Ggauge
proc norm2*(x: Ggauge): Gscalar
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

proc `+`*[T: GaugeLiteral](x: T, y: Ggauge): Ggauge =
  toGvalue(y.runtime, float(x)) + y

proc `*`*[T: GaugeLiteral](x: T, y: Ggauge): Ggauge =
  toGvalue(y.runtime, float(x)) * y

proc `-`*[T: GaugeLiteral](x: Ggauge, y: T): Ggauge =
  x - toGvalue(x.runtime, float(y))

proc retrgb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let x = Ggauge(z.inputs[0])
  let g = x.gval.newOneOf
  threads:
    for f in g:
      f := 1.0
  scaledUpstreamOr(zb, Gscalar, toGvalue(x.runtime, g))

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
  let x = Ggauge(v.inputs[0])
  let y = Ggauge(v.inputs[1])
  let z = Ggauge(v)
  z.mapGaugeSites(x.gval[mu] + y.gval[mu])

let addgg = Gfunc(forward: addggf, backward: addggb, name: "g+g")

proc `+`*(x: Ggauge, y: Ggauge): Ggauge =
  graphNode(sameShapeGaugeNodeLike(x, y, "g+g"), @[Gvalue(x), Gvalue(y)], addgg, "g+g")

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

proc expDerivgb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  if i == 0:
    raiseUnsupportedPath("expDeriv backward", "derivative with respect to force-direction input")
  raiseUnsupportedPath("expDeriv backward", "second derivatives of matrix exponential")

proc expDerivgf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let y = Ggauge(v.inputs[1])
  let z = Ggauge(v)
  z.mapGaugeElements:
    z.gval[mu][e] := expDeriv(
      y.gval[mu][e],
      x.gval[mu][e])

let expDerivg = Gfunc(forward: expDerivgf, backward: expDerivgb, name: "expDerivg")

proc expDeriv*(b: Ggauge, x: Ggauge): Ggauge =
  graphNode(
    sameShapeGaugeNodeLike(b, x, "expDerivg"),
    @[Gvalue(b), Gvalue(x)],
    expDerivg,
    "expDerivg")

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
