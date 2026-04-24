import ../[core, scalar]
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

proc gaugeScalarLiteral(anchor: Ggauge,
                        value: float): Gscalar =
  scalarLeafLike(anchor, value)

proc gaugeScalarLiteral(anchor: Ggauge,
                        value: int): Gscalar =
  scalarLeafLike(anchor, value)

proc `+`*(x: float, y: Ggauge): Ggauge =
  gaugeScalarLiteral(y, x) + y

proc `+`*(x: int, y: Ggauge): Ggauge =
  gaugeScalarLiteral(y, x) + y

proc `*`*(x: float, y: Ggauge): Ggauge =
  gaugeScalarLiteral(y, x) * y

proc `*`*(x: int, y: Ggauge): Ggauge =
  gaugeScalarLiteral(y, x) * y

proc `-`*(x: Ggauge, y: float): Ggauge =
  x - gaugeScalarLiteral(x, y)

proc `-`*(x: Ggauge, y: int): Ggauge =
  x - gaugeScalarLiteral(x, y)

proc localGaugeUpstreamRetr(zb: Gvalue, label: string): Gscalar {.inline.} =
  retr(gaugeUpstreamValue(zb, label))

proc localNegatedGaugeUpstream(zb: Gvalue, label: string): Ggauge {.inline.} =
  -gaugeUpstreamValue(zb, label)

proc localGaugeUpstreamProjTAH(zb: Gvalue, label: string): Ggauge {.inline.} =
  projTAH(gaugeUpstreamValue(zb, label))

proc localScalarGaugeAddBackward(zb: Gvalue,
                                 label: string,
                                 i: int): Gvalue =
  case i
  of 0:
    localGaugeUpstreamRetr(zb, label)
  of 1:
    gaugeUpstreamValue(zb, label)
  else:
    raiseInputIndexError(label, "0 or 1", i)

proc localSignedGaugeBinaryBackward(zb: Gvalue,
                                    label: string,
                                    i: int): Gvalue =
  case i
  of 0:
    gaugeUpstreamValue(zb, label)
  of 1:
    localNegatedGaugeUpstream(zb, label)
  else:
    raiseInputIndexError(label, "0 or 1", i)

proc localGaugeScalarSubBackward(zb: Gvalue,
                                 label: string,
                                 i: int): Gvalue =
  case i
  of 0:
    gaugeUpstreamValue(zb, label)
  of 1:
    -localGaugeUpstreamRetr(zb, label)
  else:
    raiseInputIndexError(label, "0 or 1", i)

proc retrgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireUnaryNodeView(Ggauge, "retr backward")
  discard dep
  unaryBackwardCase("retr backward", i):
    let g = view.x.getgauge.newOneOf
    threads:
      for f in g:
        f := 1.0
    return scaledUpstreamOr(zb, Gscalar, toGvalue(view.x.runtime, g), "retr backward")

defineUnaryForward(retrgf, Ggauge, Gscalar, "retr forward"):
  threads:
    var t = 0.0
    for mu in 0..<x.gval.len:
      t += x.gval[mu].trace.re
    threadMaster: z.getfloat = t

defineUnaryGraphOp(retrg, retr, Ggauge, x, scalarNodeLike(x), retrgf, retrgb, "retrg")

proc adjgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z.requireUnaryNodeView(Ggauge, "adj backward")
  discard dep
  unaryBackwardCase("adj backward", i):
    return gaugeUpstreamValue(zb, "adj backward").adj

defineUnaryForward(adjgf, Ggauge, Ggauge, "adj forward"):
  z.mapGaugeSites(x.gval[mu].adj)

defineUnaryGraphOp(adjg, adj, Ggauge, x, x.gaugeNodeLike, adjgf, adjgb, "adjg")

proc norm2gb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireUnaryNodeView(Ggauge, "norm2 backward")
  discard dep
  unaryBackwardCase("norm2 backward", i):
    return scaledUpstreamOr(zb, Gscalar, 2.0 * view.x, "norm2 backward")

defineUnaryForward(norm2gf, Ggauge, Gscalar, "norm2 forward"):
  threads:
    var t = 0.0
    for mu in 0..<x.gval.len:
      t += x.gval[mu].norm2
    threadMaster: z.getfloat = t

defineUnaryGraphOp(norm2g, norm2, Ggauge, x, scalarNodeLike(x), norm2gf, norm2gb, "norm2g")

proc neggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z.requireUnaryNodeView(Ggauge, "neg backward")
  discard dep
  unaryBackwardCase("neg backward", i):
    return localNegatedGaugeUpstream(zb, "neg backward")

defineUnaryForward(neggf, Ggauge, Ggauge, "neg forward"):
  z.mapGaugeSites(-x.gval[mu])

defineUnaryGraphOp(negg, `-`, Ggauge, x, x.gaugeNodeLike, neggf, neggb, "-g")

proc addsgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z.requireBinaryNodeView(Gscalar, Ggauge, "s+g backward")
  discard dep
  localScalarGaugeAddBackward(zb, "s+g backward", i)

defineBinaryForward(addsgf, Gscalar, Ggauge, Ggauge, "s+g forward"):
  z.mapGaugeSites(x.getfloat + y.gval[mu])

defineBinaryGraphOp(addsg, `+`, Gscalar, Ggauge, x, y, y.gaugeNodeLike, addsgf, addsgb, "s+g")

proc addggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z.requireBinaryNodeView(Ggauge, Ggauge, "g+g backward")
  discard dep
  sameGaugeBinaryBackward(zb, "g+g backward", i)

defineBinaryForward(addggf, Ggauge, Ggauge, Ggauge, "g+g forward"):
  z.mapGaugeSites(x.gval[mu] + y.gval[mu])

defineBinaryGraphOp(addgg, `+`, Ggauge, Ggauge, x, y, x.gaugeNodeLike, addggf, addggb, "g+g")

proc `+`*(x: Ggauge, y: Gvalue): Ggauge =
  x + y.requireGauge("gauge + right")

method addLike*(prototype: Ggauge, x: Gvalue, y: Gvalue): Gvalue =
  discard prototype
  if x == nil or not (x of Ggauge):
    raiseValueError("gauge gradient add left expects gauge value")
  if y == nil or not (y of Ggauge):
    raiseValueError("gauge gradient add right expects gauge value")
  Ggauge(x) + Ggauge(y)

proc mulsgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView(Gscalar, Ggauge, "s*g backward")
  discard dep
  binaryBackwardCase("s*g backward", i,
    block:
      return redot(gaugeUpstreamValue(zb, "s*g backward"), view.y),
    block:
      return view.x * gaugeUpstreamValue(zb, "s*g backward"))

defineBinaryForward(mulsgf, Gscalar, Ggauge, Ggauge, "s*g forward"):
  z.mapGaugeSites(x.getfloat * y.gval[mu])

defineBinaryGraphOp(mulsg, `*`, Gscalar, Ggauge, x, y, y.gaugeNodeLike, mulsgf, mulsgb, "s*g")

proc mulggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView(Ggauge, Ggauge, "g*g backward")
  discard dep
  binaryBackwardCase("g*g backward", i,
    block:
      return gaugeUpstreamValue(zb, "g*g backward") * view.y.adj,
    block:
      return view.x.adj * gaugeUpstreamValue(zb, "g*g backward"))

defineBinaryForward(mulggf, Ggauge, Ggauge, Ggauge, "g*g forward"):
  z.mapGaugeSites(x.gval[mu] * y.gval[mu])

defineBinaryGraphOp(mulgg, `*`, Ggauge, Ggauge, x, y, x.gaugeNodeLike, mulggf, mulggb, "g*g")

method scaleLike*(contribution: Ggauge, upstream: Gvalue): Gvalue =
  if upstream == nil:
    raiseValueError("gauge scale upstream requires non-nil value")
  if upstream of Gscalar:
    return Gscalar(upstream) * contribution
  if upstream of Ggauge:
    return Ggauge(upstream) * contribution
  raiseValueError("gauge scale upstream expects scalar or gauge value, got:\n" & upstream.nodeRepr)

proc `*`*(x: Ggauge, y: Gvalue): Ggauge =
  x * y.requireGauge("gauge * right")

proc redotggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard dep
  let view = z.requireBinaryNodeView(Ggauge, Ggauge, "redotgg backward")
  binaryBackwardCase("redotgg backward", i,
    block:
      return scaledUpstreamOr(zb, Gscalar, view.y, "redotgg backward"),
    block:
      return scaledUpstreamOr(zb, Gscalar, view.x, "redotgg backward"))

defineBinaryForward(redotggf, Ggauge, Ggauge, Gscalar, "redot forward"):
  threads:
    var t = 0.0
    for mu in 0..<x.gval.len:
      t += redot(x.gval[mu], y.gval[mu])
    threadMaster: z.getfloat = t

defineBinaryGraphOp(redotgg, redot, Ggauge, Ggauge, x, y, scalarNodeLike(x), redotggf, redotggb, "redotgg")

proc redot*(x: Ggauge, y: Gvalue): Gscalar =
  redot(x, y.requireGauge("redot right"))

proc subgsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z.requireBinaryNodeView(Ggauge, Gscalar, "g-s backward")
  discard dep
  localGaugeScalarSubBackward(zb, "g-s backward", i)

defineBinaryForward(subgsf, Ggauge, Gscalar, Ggauge, "g-s forward"):
  z.mapGaugeSites(x.gval[mu] - y.getfloat)

defineBinaryGraphOp(subgs, `-`, Ggauge, Gscalar, x, y, x.gaugeNodeLike, subgsf, subgsb, "g-s")

proc subggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z.requireBinaryNodeView(Ggauge, Ggauge, "g-g backward")
  discard dep
  localSignedGaugeBinaryBackward(zb, "g-g backward", i)

defineBinaryForward(subggf, Ggauge, Ggauge, Ggauge, "g-g forward"):
  z.mapGaugeSites(x.gval[mu] - y.gval[mu])

defineBinaryGraphOp(subgg, `-`, Ggauge, Ggauge, x, y, x.gaugeNodeLike, subggf, subggb, "g-g")

proc `-`*(x: Ggauge, y: Gvalue): Ggauge =
  if y of Gscalar:
    return x - Gscalar(y)
  x - y.requireGauge("gauge - right")

proc expgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireUnaryNodeView(Ggauge, "exp backward")
  discard dep
  unaryBackwardCase("exp backward", i):
    return expDeriv(gaugeUpstreamValue(zb, "exp backward"), view.x)

defineUnaryForward(expgf, Ggauge, Ggauge, "exp forward"):
  z.mapGaugeElements:
    z.gval[mu][e] := exp(x.gval[mu][e])

defineUnaryGraphOp(expg, exp, Ggauge, x, x.gaugeNodeLike, expgf, expgb, "expg")

proc expDerivgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard zb
  discard z.requireBinaryNodeView(Ggauge, Ggauge, "expDeriv backward")
  discard dep
  unaryBackwardCase("expDeriv backward", i):
    raiseUnsupportedPath("expDeriv backward", "second derivatives of matrix exponential")

defineBinaryForward(expDerivgf, Ggauge, Ggauge, Ggauge, "expDeriv forward"):
  z.mapGaugeElements:
    z.gval[mu][e] := expDeriv(y.gval[mu][e], x.gval[mu][e])

defineBinaryGraphOp(expDerivg, expDeriv, Ggauge, Ggauge, b, x, x.gaugeNodeLike, expDerivgf, expDerivgb, "expDerivg")

proc expDeriv*(b: Ggauge, x: Gvalue): Ggauge =
  expDeriv(b, x.requireGauge("expDeriv exponent"))

proc projTAHb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z.requireUnaryNodeView(Ggauge, "projTAH backward")
  discard dep
  unaryBackwardCase("projTAH backward", i):
    return localGaugeUpstreamProjTAH(zb, "projTAH backward")

defineUnaryForward(projTAHf, Ggauge, Ggauge, "projTAH forward"):
  z.mapGaugeElements:
    z.gval[mu][e].projectTAH(x.gval[mu][e])

defineUnaryGraphOp(projTAHg, projTAH, Ggauge, x, x.gaugeNodeLike, projTAHf, projTAHb, "projTAH")
