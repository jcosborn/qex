import ../[core, scalar, multi]
import ../support/op
import layout, physics/qcdTypes
import shared, basic_ops

# Section: Fused Gauge Ops

proc adjmul*(x: Ggauge, y: Ggauge): Ggauge
proc muladj*(x: Ggauge, y: Ggauge): Ggauge
proc contractProjTAH*(x: Ggauge, y: Ggauge): Ggauge
proc axexp*(a: Gscalar, x: Ggauge): Ggauge
proc axexpmuly*(a: Gscalar, x: Ggauge, y: Ggauge): Ggauge

proc requirePackedArity(pack: Gmulti,
                        expected: int,
                        label: string) =
  let actual = pack.len
  if actual != expected:
    raiseValueError(
      label & " expects " & $expected &
      " packed slots, got " & $actual)

proc packedSlotValue[T: Gvalue](pack: Gmulti,
                                index: int,
                                slotType: typedesc[T],
                                label: string,
                                slotLabel: string): T =
  let value = pack.slotValue(index)
  if not (value of T):
    raiseValueError(
      label & " expects " & slotLabel &
      " slot " & $index & " of type " & $slotType &
      ", got:\n" & value.nodeRepr)
  T(value)

proc packedSlotNode[T: Gvalue](pack: Gmulti,
                               index: int,
                               slotType: typedesc[T],
                               label: string,
                               slotLabel: string): T =
  discard pack.packedSlotValue(index, slotType, label, slotLabel)
  T(pack[index])

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

proc contractProjTAHPackedInputb(zb: Gvalue,
                                 z: Gvalue,
                                 i: int,
                                 dep: Gvalue): Gvalue =
  let view = z.requireUnaryNodeView(Gmulti, "contractProjTAH packed backward", "args")
  discard dep
  unaryBackwardCase("contractProjTAH packed backward", i):
    let args = view.x
    args.requirePackedArity(2, "contractProjTAH packed backward")
    let x = args.packedSlotNode(0, Ggauge, "contractProjTAH packed backward", "left")
    let y = args.packedSlotNode(1, Ggauge, "contractProjTAH packed backward", "right")
    let proj = projTAH(gaugeUpstreamValue(zb, "contractProjTAH packed backward"))
    return multiValues(
      [Gvalue(proj * y), Gvalue(proj.adjmul x)],
      "contractProjTAH packed input gradients")

proc contractProjTAHPackedInputf(v: Gvalue) =
  let view = v.requireUnaryNodeView(Gmulti, "contractProjTAH packed forward", "args")
  let args = view.x
  args.requirePackedArity(2, "contractProjTAH packed forward")
  let x = args.packedSlotValue(0, Ggauge, "contractProjTAH packed forward", "left")
  let y = args.packedSlotValue(1, Ggauge, "contractProjTAH packed forward", "right")
  let z = v.requireGauge("contractProjTAH packed forward result")
  z.mapGaugeElements:
    let s = x.gval[mu][e] * y.gval[mu][e].adj
    z.gval[mu][e].projectTAH s

let contractProjTAHPackedInputg = newGfunc(
  forward = contractProjTAHPackedInputf,
  backward = contractProjTAHPackedInputb,
  name = "contractProjTAH packed")

proc contractProjTAHPacked(args: Gmulti): Ggauge =
  args.requirePackedArity(2, "contractProjTAH packed")
  let x = args.packedSlotValue(0, Ggauge, "contractProjTAH packed", "left")
  discard args.packedSlotValue(1, Ggauge, "contractProjTAH packed", "right")
  graphNode(x.gaugeNodeLike, @[Gvalue(args)], contractProjTAHPackedInputg, "contractProjTAH packed")

proc contractProjTAH*(x: Ggauge, y: Ggauge): Ggauge =
  let args = multiValues([Gvalue(x), Gvalue(y)], "contractProjTAH args")
  contractProjTAHPacked(args)

proc contractProjTAH*(x: Ggauge, y: Gvalue): Ggauge =
  contractProjTAH(x, y.requireGauge("contractProjTAH right"))

proc axexpPackedInputb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireUnaryNodeView(Gmulti, "axexp packed backward", "args")
  discard dep
  unaryBackwardCase("axexp packed backward", i):
    let args = view.x
    args.requirePackedArity(2, "axexp packed backward")
    let a = args.packedSlotNode(0, Gscalar, "axexp packed backward", "scale")
    let x = args.packedSlotNode(1, Ggauge, "axexp packed backward", "exponent")
    let ax = a * x
    let deriv = expDeriv(gaugeUpstreamValue(zb, "axexp packed backward"), ax)
    return multiValues(
      [Gvalue(deriv.redot x), Gvalue(a * deriv)],
      "axexp packed input gradients")

proc axexpPackedInputf(v: Gvalue) =
  let view = v.requireUnaryNodeView(Gmulti, "axexp packed forward", "args")
  let args = view.x
  args.requirePackedArity(2, "axexp packed forward")
  let a = args.packedSlotValue(0, Gscalar, "axexp packed forward", "scale")
  let x = args.packedSlotValue(1, Ggauge, "axexp packed forward", "exponent")
  let z = v.requireGauge("axexp packed forward result")
  let f = a.getfloat
  z.mapGaugeElements:
    z.gval[mu][e] := exp(f * x.gval[mu][e])

let axexpPackedInputg = newGfunc(
  forward = axexpPackedInputf,
  backward = axexpPackedInputb,
  name = "axexp packed")

proc axexpPacked(args: Gmulti): Ggauge =
  args.requirePackedArity(2, "axexp packed")
  discard args.packedSlotValue(0, Gscalar, "axexp packed", "scale")
  let x = args.packedSlotValue(1, Ggauge, "axexp packed", "exponent")
  graphNode(x.gaugeNodeLike, @[Gvalue(args)], axexpPackedInputg, "axexp packed")

proc axexp*(a: Gscalar, x: Ggauge): Ggauge =
  let args = multiValues([Gvalue(a), Gvalue(x)], "axexp args")
  axexpPacked(args)

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

proc fillAxexpmulyPack(pack: Gmulti,
                       a: Gscalar,
                       x: Ggauge,
                       y: Ggauge) =
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

proc axexpmulyPackedInputPackb(zb: Gvalue,
                               z: Gvalue,
                               i: int,
                               dep: Gvalue): Gvalue =
  let view = z.requireUnaryNodeView(Gmulti, "axexpmulyPack packed backward", "args")
  discard dep
  unaryBackwardCase("axexpmulyPack packed backward", i):
    let args = view.x
    args.requirePackedArity(3, "axexpmulyPack packed backward")
    let a = args.packedSlotNode(0, Gscalar, "axexpmulyPack packed backward", "scale")
    let x = args.packedSlotNode(1, Ggauge, "axexpmulyPack packed backward", "exponent")
    let y = args.packedSlotNode(2, Ggauge, "axexpmulyPack packed backward", "value")
    let upstream = requireMultiUpstream(zb, "axexpmulyPack packed backward")
    let expUp = upstream.packExpNode
    let resultUp = upstream.packResultNode
    let ax = a * x
    let deriv = axexpmulyPackDeriv(expUp, resultUp, ax, y)
    let expNode = z.requireMultiValue("axexpmulyPack packed backward result").packExpNode
    return multiValues(
      [
        Gvalue(deriv.redot x),
        Gvalue(a * deriv),
        Gvalue(expNode.adjmul resultUp)
      ],
      "axexpmulyPack packed input gradients")

proc axexpmulyPackedInputPackf(v: Gvalue) =
  let view = v.requireUnaryNodeView(Gmulti, "axexpmulyPack packed forward", "args")
  let args = view.x
  args.requirePackedArity(3, "axexpmulyPack packed forward")
  let a = args.packedSlotValue(0, Gscalar, "axexpmulyPack packed forward", "scale")
  let x = args.packedSlotValue(1, Ggauge, "axexpmulyPack packed forward", "exponent")
  let y = args.packedSlotValue(2, Ggauge, "axexpmulyPack packed forward", "value")
  let pack = v.requireMultiValue("axexpmulyPack packed forward result")
  pack.fillAxexpmulyPack(a, x, y)

let axexpmulyPackedInputPackg = newGfunc(
  forward = axexpmulyPackedInputPackf,
  backward = axexpmulyPackedInputPackb,
  name = "axexpmulyPack packed")

proc axexpmulyPack(args: Gmulti): Gmulti =
  args.requirePackedArity(3, "axexpmulyPack packed")
  discard args.packedSlotValue(0, Gscalar, "axexpmulyPack packed", "scale")
  let x = args.packedSlotValue(1, Ggauge, "axexpmulyPack packed", "exponent")
  discard args.packedSlotValue(2, Ggauge, "axexpmulyPack packed", "value")
  newMultiOutputNode(@[Gvalue(x), x], @[Gvalue(args)], axexpmulyPackedInputPackg, "axexpmulyPack packed")

proc axexpmuly*(a: Gscalar, x: Ggauge, y: Ggauge): Ggauge =
  let args = multiValues([Gvalue(a), Gvalue(x), Gvalue(y)], "axexpmuly args")
  axexpmulyPack(args).packResultNode
