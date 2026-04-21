import ../[core, scalar]
import ../support/op
import layout, gauge, physics/qcdTypes
import shared

# Section: Basic Gauge Ops

method redot*(x: Gscalar, y: Gscalar): Gvalue = x*y

proc retrgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireUnaryNodeView("retr backward")
  discard dep
  unaryBackwardCase("retr backward", i):
    let g = view.x.getgauge.newOneOf
    threads:
      for f in g:
        f := 1.0
    return scaledUpstreamOr(zb, toGvalue(view.x.runtime, g))

defineUnaryForward(retrgf, Ggauge, Gscalar, "retr forward"):
  threads:
    var t = 0.0
    for mu in 0..<x.gval.len:
      t += x.gval[mu].trace.re
    threadMaster: z.getfloat = t

defineUnaryGraphOp(retrg, retr, Ggauge, x, scalarNodeLike(x), retrgf, retrgb, "retrg")

proc adjgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z.requireUnaryNodeView("adj backward")
  discard dep
  unaryBackwardCase("adj backward", i):
    return gaugeUpstreamValue(zb, "adj backward").adj

defineUnaryForward(adjgf, Ggauge, Ggauge, "adj forward"):
  z.mapGaugeSites(x.gval[mu].adj)

defineUnaryGraphOp(adjg, adj, Ggauge, x, x.newOneOf, adjgf, adjgb, "adjg")

proc norm2gb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireUnaryNodeView("norm2 backward")
  discard dep
  unaryBackwardCase("norm2 backward", i):
    return scaledUpstreamOr(zb, 2.0 * view.x)

defineUnaryForward(norm2gf, Ggauge, Gscalar, "norm2 forward"):
  threads:
    var t = 0.0
    for mu in 0..<x.gval.len:
      t += x.gval[mu].norm2
    threadMaster: z.getfloat = t

defineUnaryGraphOp(norm2g, norm2, Ggauge, x, scalarNodeLike(x), norm2gf, norm2gb, "norm2g")

proc neggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z.requireUnaryNodeView("neg backward")
  discard dep
  unaryBackwardCase("neg backward", i):
    return negatedGaugeUpstream(zb, "neg backward")

defineUnaryForward(neggf, Ggauge, Ggauge, "neg forward"):
  z.mapGaugeSites(-x.gval[mu])

defineUnaryGraphOp(negg, `-`, Ggauge, x, x.newOneOf, neggf, neggb, "-g")

proc addsgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z.requireBinaryNodeView("s+g backward")
  discard dep
  scalarGaugeAddBackward(zb, "s+g backward", i)

defineBinaryForward(addsgf, Gscalar, Ggauge, Ggauge, "s+g forward"):
  z.mapGaugeSites(x.getfloat + y.gval[mu])

defineBinaryGraphOp(addsg, `+`, Gscalar, Ggauge, x, y, y.newOneOf, addsgf, addsgb, "s+g")

proc addggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z.requireBinaryNodeView("g+g backward")
  discard dep
  sameGaugeBinaryBackward(zb, "g+g backward", i)

defineBinaryForward(addggf, Ggauge, Ggauge, Ggauge, "g+g forward"):
  z.mapGaugeSites(x.gval[mu] + y.gval[mu])

defineBinaryGraphOp(addgg, `+`, Ggauge, Ggauge, x, y, x.newOneOf, addggf, addggb, "g+g")

proc mulsgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView("s*g backward")
  discard dep
  binaryBackwardCase("s*g backward", i,
    block:
      return redot(gaugeUpstreamValue(zb, "s*g backward"), view.y),
    block:
      return view.x * gaugeUpstreamValue(zb, "s*g backward"))

defineBinaryForward(mulsgf, Gscalar, Ggauge, Ggauge, "s*g forward"):
  z.mapGaugeSites(x.getfloat * y.gval[mu])

defineBinaryGraphOp(mulsg, `*`, Gscalar, Ggauge, x, y, y.newOneOf, mulsgf, mulsgb, "s*g")

proc mulggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView("g*g backward")
  discard dep
  binaryBackwardCase("g*g backward", i,
    block:
      return gaugeUpstreamValue(zb, "g*g backward").muladj view.y,
    block:
      return view.x.adjmul gaugeUpstreamValue(zb, "g*g backward"))

defineBinaryForward(mulggf, Ggauge, Ggauge, Ggauge, "g*g forward"):
  z.mapGaugeSites(x.gval[mu] * y.gval[mu])

defineBinaryGraphOp(mulgg, `*`, Ggauge, Ggauge, x, y, x.newOneOf, mulggf, mulggb, "g*g")

proc redotggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard dep
  swappedScaledInputBackward(zb, z, "redotgg backward", i)

defineBinaryForward(redotggf, Ggauge, Ggauge, Gscalar, "redot forward"):
  threads:
    var t = 0.0
    for mu in 0..<x.gval.len:
      t += redot(x.gval[mu], y.gval[mu])
    threadMaster: z.getfloat = t

defineBinaryGraphOp(redotgg, redot, Ggauge, Ggauge, x, y, scalarNodeLike(x), redotggf, redotggb, "redotgg")

proc subgsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z.requireBinaryNodeView("g-s backward")
  discard dep
  gaugeScalarSubBackward(zb, "g-s backward", i)

defineBinaryForward(subgsf, Ggauge, Gscalar, Ggauge, "g-s forward"):
  z.mapGaugeSites(x.gval[mu] - y.getfloat)

defineBinaryGraphOp(subgs, `-`, Ggauge, Gscalar, x, y, x.newOneOf, subgsf, subgsb, "g-s")

proc subggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z.requireBinaryNodeView("g-g backward")
  discard dep
  signedGaugeBinaryBackward(zb, "g-g backward", i)

defineBinaryForward(subggf, Ggauge, Ggauge, Ggauge, "g-g forward"):
  z.mapGaugeSites(x.gval[mu] - y.gval[mu])

defineBinaryGraphOp(subgg, `-`, Ggauge, Ggauge, x, y, x.newOneOf, subggf, subggb, "g-g")

proc expgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireUnaryNodeView("exp backward")
  discard dep
  unaryBackwardCase("exp backward", i):
    return expDeriv(gaugeUpstreamValue(zb, "exp backward"), view.x)

defineUnaryForward(expgf, Ggauge, Ggauge, "exp forward"):
  z.mapGaugeElements:
    z.gval[mu][e] := exp(x.gval[mu][e])

defineUnaryGraphOp(expg, exp, Ggauge, x, x.newOneOf, expgf, expgb, "expg")

proc expDerivgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard zb
  discard z.requireBinaryNodeView("expDeriv backward")
  discard dep
  unaryBackwardCase("expDeriv backward", i):
    raiseUnsupportedPath("expDeriv backward", "second derivatives of matrix exponential")

defineBinaryForward(expDerivgf, Ggauge, Ggauge, Ggauge, "expDeriv forward"):
  z.mapGaugeElements:
    z.gval[mu][e] := expDeriv(y.gval[mu][e], x.gval[mu][e])

defineBinaryGraphOp(expDerivg, expDeriv, Ggauge, Ggauge, b, x, x.newOneOf, expDerivgf, expDerivgb, "expDerivg")

proc projTAHb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z.requireUnaryNodeView("projTAH backward")
  discard dep
  unaryBackwardCase("projTAH backward", i):
    return gaugeUpstreamProjTAH(zb, "projTAH backward")

defineUnaryForward(projTAHf, Ggauge, Ggauge, "projTAH forward"):
  z.mapGaugeElements:
    z.gval[mu][e].projectTAH(x.gval[mu][e])

defineUnaryGraphOp(projTAHg, projTAH, Ggauge, x, x.newOneOf, projTAHf, projTAHb, "projTAH")
