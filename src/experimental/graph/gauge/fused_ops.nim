import ../[core, scalar, multi]
import ../support/op
import layout, physics/qcdTypes
import shared

# Section: Fused Gauge Ops

proc adjmulggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView("adjmul backward")
  discard dep
  # adjmul(x, y) == x.adj * y
  binaryBackwardCase("adjmul backward", i,
    block:
      return view.y.muladj gaugeUpstreamValue(zb, "adjmul backward"),
    block:
      return view.x * gaugeUpstreamValue(zb, "adjmul backward"))

defineBinaryForward(adjmulggf, Ggauge, Ggauge, Ggauge, "adjmul forward"):
  z.mapGaugeSites(x.gval[mu].adj * y.gval[mu])

defineBinaryGraphOp(adjmulgg, adjmul, Ggauge, Ggauge, x, y, x.newOneOf, adjmulggf, adjmulggb, "g.adj*g")

proc muladjggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView("muladj backward")
  discard dep
  # muladj(x, y) == x * y.adj
  binaryBackwardCase("muladj backward", i,
    block:
      return gaugeUpstreamValue(zb, "muladj backward") * view.y,
    block:
      return gaugeUpstreamValue(zb, "muladj backward").adjmul view.x)

defineBinaryForward(muladjggf, Ggauge, Ggauge, Ggauge, "muladj forward"):
  z.mapGaugeSites(x.gval[mu] * y.gval[mu].adj)

defineBinaryGraphOp(muladjgg, muladj, Ggauge, Ggauge, x, y, x.newOneOf, muladjggf, muladjggb, "g*g.adj")

proc contractProjTAHb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView("contractProjTAH backward")
  discard dep
  # contractProjTAH(x, y) == projTAH(x * y.adj)
  let proj = gaugeUpstreamProjTAH(zb, "contractProjTAH backward")
  binaryBackwardCase("contractProjTAH backward", i,
    block:
      return proj * view.y,
    block:
      return proj.adjmul view.x)

defineBinaryForward(contractProjTAHf, Ggauge, Ggauge, Ggauge, "contractProjTAH forward"):
  z.mapGaugeElements:
    let s = x.gval[mu][e] * y.gval[mu][e].adj
    z.gval[mu][e].projectTAH s

defineBinaryGraphOp(contractProjTAHg, contractProjTAH, Ggauge, Ggauge, x, y, x.newOneOf, contractProjTAHf, contractProjTAHb, "projTAHg")

proc axexpb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView("axexp backward")
  discard dep
  # axexp(a, x) == exp(a * x)
  let deriv = expDeriv(gaugeUpstreamValue(zb, "axexp backward"), view.x * view.y)
  binaryBackwardCase("axexp backward", i,
    block:
      return deriv.redot view.y,
    block:
      return view.x * deriv)

defineBinaryForward(axexpf, Gscalar, Ggauge, Ggauge, "axexp forward"):
  let f = x.getfloat
  z.mapGaugeElements:
    z.gval[mu][e] := exp(f * y.gval[mu][e])

defineBinaryGraphOp(axexpg, axexp, Gscalar, Ggauge, a, x, x.newOneOf, axexpf, axexpb, "axexp")

const
  axexpmulyPackExpSlot = 0
  axexpmulyPackResultSlot = 1

proc packExpValue(pack: Gmulti): Gvalue =
  pack.slotValue(axexpmulyPackExpSlot)

proc packResultValue(pack: Gmulti): Gvalue =
  pack.slotValue(axexpmulyPackResultSlot)

proc packExpNode(pack: Gmulti): Gvalue =
  pack[axexpmulyPackExpSlot]

proc packResultNode(pack: Gmulti): Gvalue =
  pack[axexpmulyPackResultSlot]

proc packExpGauge(pack: Gmulti): Ggauge =
  Ggauge(pack.packExpValue)

proc packResultGauge(pack: Gmulti): Ggauge =
  Ggauge(pack.packResultValue)

proc axexpmulyPackDeriv(expUp: Gvalue,
                        resultUp: Gvalue,
                        ax: Gvalue,
                        y: Gvalue): Gvalue =
  expDeriv(expUp + resultUp.muladj y, ax)

proc axexpmulyPackb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireTernaryNodeView(
    "axexpmulyPack backward",
    "scale",
    "exponent",
    "value")
  discard dep
  # axexpmulyPack(a, x, y) returns:
  #   slot 0 = exp(a * x)
  #   slot 1 = exp(a * x) * y
  let upstream = requireMultiUpstream(zb, "axexpmulyPack backward")
  let expUp = upstream.packExpNode
  let resultUp = upstream.packResultNode
  let ax = view.x * view.y
  let deriv = axexpmulyPackDeriv(expUp, resultUp, ax, view.z)
  ternaryBackwardCase("axexpmulyPack backward", i,
    block:
      return deriv.redot view.y,
    block:
      return view.x * deriv,
    block:
      return z.requireMultiValue("axexpmulyPack backward result").packExpNode.adjmul resultUp)

proc axexpmulyPackf(v: Gvalue) =
  let view = v.requireTernaryNodeView(
    "axexpmulyPack forward",
    "scale",
    "exponent",
    "value")
  let a = Gscalar(view.x)
  let x = Ggauge(view.y)
  let y = Ggauge(view.z)
  let pack = v.requireMultiValue("axexpmulyPack forward result")
  let expax = pack.packExpGauge
  let result = pack.packResultGauge
  let f = a.getfloat
  threads:
    for mu in 0..<result.gval.len:
      for e in result.gval[mu]:
        var t{.noinit.}: evalType(x.gval[mu][e])
        t := exp(f * x.gval[mu][e])
        expax.gval[mu][e] := t
        result.gval[mu][e] := t * y.gval[mu][e]

let axexpmulyPackg = newGfunc(forward = axexpmulyPackf, backward = axexpmulyPackb, name = "axexpmulyPack")

method axexpmulyPack*(a: Gscalar, x: Ggauge, y: Ggauge): Gmulti =
  newMultiOutputNode(@[Gvalue(x), x], @[Gvalue(a), x, y], axexpmulyPackg, "axexpmulyPack")

method axexpmuly*(a: Gscalar, x: Ggauge, y: Ggauge): Gvalue =
  axexpmulyPack(a, x, y).packResultNode
