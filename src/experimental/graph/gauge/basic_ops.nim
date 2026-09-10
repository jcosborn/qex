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
