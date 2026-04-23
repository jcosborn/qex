import ../[core, scalar, multi]
import ../support/op
import layout, physics/qcdTypes
import shared, basic_ops

# Section: Fused Gauge Ops

proc adjmul*(x: Ggauge, y: Ggauge): Ggauge
proc muladj*(x: Ggauge, y: Ggauge): Ggauge
proc contractProjTAH*(x: Ggauge, y: Ggauge): Ggauge
proc axexp*(a: Gscalar, x: Ggauge): Ggauge
proc axexpmulyPack*(a: Gscalar, x: Ggauge, y: Ggauge): Gmulti
proc axexpmuly*(a: Gscalar, x: Ggauge, y: Ggauge): Ggauge

proc adjmulggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView(Ggauge, Ggauge, "adjmul backward")
  discard dep
  # adjmul(x, y) == x.adj * y
  binaryBackwardCase("adjmul backward", i,
    block:
      return view.y.muladj gaugeUpstreamValue(zb, "adjmul backward"),
    block:
      return view.x * gaugeUpstreamValue(zb, "adjmul backward"))

defineBinaryForward(adjmulggf, Ggauge, Ggauge, Ggauge, "adjmul forward"):
  z.mapGaugeSites(x.gval[mu].adj * y.gval[mu])

defineBinaryGraphOp(adjmulgg, adjmul, Ggauge, Ggauge, x, y, x.gaugeNodeLike, adjmulggf, adjmulggb, "g.adj*g")

proc adjmul*(x: Ggauge, y: Gvalue): Ggauge =
  adjmul(x, y.requireGauge("adjmul right"))

proc muladjggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView(Ggauge, Ggauge, "muladj backward")
  discard dep
  # muladj(x, y) == x * y.adj
  binaryBackwardCase("muladj backward", i,
    block:
      return gaugeUpstreamValue(zb, "muladj backward") * view.y,
    block:
      return gaugeUpstreamValue(zb, "muladj backward").adjmul view.x)

defineBinaryForward(muladjggf, Ggauge, Ggauge, Ggauge, "muladj forward"):
  z.mapGaugeSites(x.gval[mu] * y.gval[mu].adj)

defineBinaryGraphOp(muladjgg, muladj, Ggauge, Ggauge, x, y, x.gaugeNodeLike, muladjggf, muladjggb, "g*g.adj")

proc muladj*(x: Ggauge, y: Gvalue): Ggauge =
  muladj(x, y.requireGauge("muladj right"))

proc contractProjTAHb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView(Ggauge, Ggauge, "contractProjTAH backward")
  discard dep
  # contractProjTAH(x, y) == projTAH(x * y.adj)
  let proj = projTAH(gaugeUpstreamValue(zb, "contractProjTAH backward"))
  binaryBackwardCase("contractProjTAH backward", i,
    block:
      return proj * view.y,
    block:
      return proj.adjmul view.x)

defineBinaryForward(contractProjTAHf, Ggauge, Ggauge, Ggauge, "contractProjTAH forward"):
  z.mapGaugeElements:
    let s = x.gval[mu][e] * y.gval[mu][e].adj
    z.gval[mu][e].projectTAH s

defineBinaryGraphOp(contractProjTAHg, contractProjTAH, Ggauge, Ggauge, x, y, x.gaugeNodeLike, contractProjTAHf, contractProjTAHb, "projTAHg")

proc contractProjTAH*(x: Ggauge, y: Gvalue): Ggauge =
  contractProjTAH(x, y.requireGauge("contractProjTAH right"))

proc axexpb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView(Gscalar, Ggauge, "axexp backward")
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

defineBinaryGraphOp(axexpg, axexp, Gscalar, Ggauge, a, x, x.gaugeNodeLike, axexpf, axexpb, "axexp")

const
  axexpmulyPackExpSlot = 0
  axexpmulyPackResultSlot = 1

proc packExpValue(pack: Gmulti): Ggauge =
  Ggauge(pack.slotValue(axexpmulyPackExpSlot))

proc packResultValue(pack: Gmulti): Ggauge =
  Ggauge(pack.slotValue(axexpmulyPackResultSlot))

proc packExpNode(pack: Gmulti): Ggauge =
  Ggauge(pack[axexpmulyPackExpSlot])

proc packResultNode(pack: Gmulti): Ggauge =
  Ggauge(pack[axexpmulyPackResultSlot])

proc packExpGauge(pack: Gmulti): Ggauge =
  pack.packExpValue

proc packResultGauge(pack: Gmulti): Ggauge =
  pack.packResultValue

proc axexpmulyPackDeriv(expUp: Ggauge,
                        resultUp: Ggauge,
                        ax: Ggauge,
                        y: Ggauge): Ggauge =
  expDeriv(expUp + resultUp.muladj y, ax)

proc axexpmulyPackb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireTernaryNodeView(
    Gscalar,
    Ggauge,
    Ggauge,
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

proc axexpmulyPack*(a: Gscalar, x: Ggauge, y: Ggauge): Gmulti =
  newMultiOutputNode(@[Gvalue(x), x], @[Gvalue(a), x, y], axexpmulyPackg, "axexpmulyPack")

proc axexpmuly*(a: Gscalar, x: Ggauge, y: Ggauge): Ggauge =
  axexpmulyPack(a, x, y).packResultNode
