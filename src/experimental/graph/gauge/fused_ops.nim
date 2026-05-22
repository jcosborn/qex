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

proc adjmulggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView(Ggauge, Ggauge, "adjmul backward")
  discard dep
  # adjmul(x, y) == x.adj * y
  binaryBackwardCase("adjmul backward", i,
    block:
      return view.y.muladj(requireUpstream(zb, Ggauge, "adjmul backward")),
    block:
      return view.x * requireUpstream(zb, Ggauge, "adjmul backward"))

proc adjmulggf(v: Gvalue) =
  let view = v.requireBinaryNodeView(Ggauge, Ggauge, "adjmul forward")
  let x = view.x
  let y = view.y
  let z = v.requireGauge("adjmul forward result")
  z.mapGaugeSites(x.unsafeGaugeStorage[mu].adj * y.unsafeGaugeStorage[mu])

let adjmulgg = newGfunc(forward = adjmulggf, backward = adjmulggb, name = "g.adj*g")

proc adjmul*(x: Ggauge, y: Ggauge): Ggauge =
  graphNode(sameShapeGaugeNodeLike(x, y, "g.adj*g"), @[Gvalue(x), Gvalue(y)], adjmulgg, "g.adj*g")

proc muladjggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView(Ggauge, Ggauge, "muladj backward")
  discard dep
  # muladj(x, y) == x * y.adj
  binaryBackwardCase("muladj backward", i,
    block:
      return requireUpstream(zb, Ggauge, "muladj backward") * view.y,
    block:
      return requireUpstream(zb, Ggauge, "muladj backward").adjmul(view.x))

proc muladjggf(v: Gvalue) =
  let view = v.requireBinaryNodeView(Ggauge, Ggauge, "muladj forward")
  let x = view.x
  let y = view.y
  let z = v.requireGauge("muladj forward result")
  z.mapGaugeSites(x.unsafeGaugeStorage[mu] * y.unsafeGaugeStorage[mu].adj)

let muladjgg = newGfunc(forward = muladjggf, backward = muladjggb, name = "g*g.adj")

proc muladj*(x: Ggauge, y: Ggauge): Ggauge =
  graphNode(sameShapeGaugeNodeLike(x, y, "g*g.adj"), @[Gvalue(x), Gvalue(y)], muladjgg, "g*g.adj")

proc contractProjTAHPackedInputb(zb: Gvalue,
                                 z: Gvalue,
                                 i: int,
                                 dep: Gvalue): Gvalue =
  let view = z.requireUnaryNodeView(Gmulti, "contractProjTAH packed backward", "args")
  discard dep
  unaryBackwardCase("contractProjTAH packed backward", i):
    # Args: [x, y]
    let (x, y) = symbolicSlots[
      tuple[x: Ggauge, y: Ggauge]](view.x, "contractProjTAH packed backward")
    let proj = projTAH(
      requireUpstream(zb, Ggauge, "contractProjTAH packed backward"))
    return multiValues(
      "contractProjTAH packed input gradients",
      proj * y,
      proj.adjmul x)

proc contractProjTAHPackedInputf(v: Gvalue) =
  let view = v.requireUnaryNodeView(Gmulti, "contractProjTAH packed forward", "args")
  # Args: [x, y]
  let (x, y) = storedSlots[
    tuple[x: Ggauge, y: Ggauge]](view.x, "contractProjTAH packed forward")
  let z = v.requireGauge("contractProjTAH packed forward result")
  z.mapGaugeElements:
    let s = x.unsafeGaugeStorage[mu][e] * y.unsafeGaugeStorage[mu][e].adj
    z.unsafeGaugeStorage[mu][e].projectTAH s

let contractProjTAHPackedInputg = newGfunc(
  forward = contractProjTAHPackedInputf,
  backward = contractProjTAHPackedInputb,
  name = "contractProjTAH packed")

proc contractProjTAH*(x: Ggauge, y: Ggauge): Ggauge =
  x.requireSameGaugeShape(y, "contractProjTAH packed")
  let args = multiValues("contractProjTAH args", x, y)
  graphNode(
    x.gaugeNodeLike,
    @[Gvalue(args)],
    contractProjTAHPackedInputg,
    "contractProjTAH packed")

proc axexpPackedInputb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireUnaryNodeView(Gmulti, "axexp packed backward", "args")
  discard dep
  unaryBackwardCase("axexp packed backward", i):
    # Args: [a, x]
    let (a, x) = symbolicSlots[
      tuple[a: Gscalar, x: Ggauge]](view.x, "axexp packed backward")
    let ax = a * x
    let deriv = expDeriv(requireUpstream(zb, Ggauge, "axexp packed backward"), ax)
    return multiValues(
      "axexp packed input gradients",
      deriv.redot x,
      a * deriv)

proc axexpPackedInputf(v: Gvalue) =
  let view = v.requireUnaryNodeView(Gmulti, "axexp packed forward", "args")
  # Args: [a, x]
  let (a, x) = storedSlots[
    tuple[a: Gscalar, x: Ggauge]](view.x, "axexp packed forward")
  let z = v.requireGauge("axexp packed forward result")
  let f = a.sval
  z.mapGaugeElements:
    z.unsafeGaugeStorage[mu][e] := exp(f * x.unsafeGaugeStorage[mu][e])

let axexpPackedInputg = newGfunc(
  forward = axexpPackedInputf,
  backward = axexpPackedInputb,
  name = "axexp packed")

proc axexp*(a: Gscalar, x: Ggauge): Ggauge =
  let args = multiValues("axexp args", a, x)
  graphNode(x.gaugeNodeLike, @[Gvalue(args)], axexpPackedInputg, "axexp packed")

proc axexpmulyPackedInputPackb(zb: Gvalue,
                               z: Gvalue,
                               i: int,
                               dep: Gvalue): Gvalue =
  let view = z.requireUnaryNodeView(Gmulti, "axexpmulyPack packed backward", "args")
  discard dep
  unaryBackwardCase("axexpmulyPack packed backward", i):
    # Args: [a, x, y]
    let (a, x, y) = symbolicSlots[
      tuple[a: Gscalar, x: Ggauge, y: Ggauge]](
        view.x,
        "axexpmulyPack packed backward")
    let upstream = requireMultiUpstream(zb, "axexpmulyPack packed backward")
    # Upstream/result slots: [exp(a*x), exp(a*x)*y]
    let (upstreamExpax, upstreamValue) = symbolicSlots[
      tuple[expax: Ggauge, value: Ggauge]](upstream, "axexpmulyPack upstream")
    let ax = a * x
    let deriv = expDeriv(upstreamExpax + upstreamValue.muladj y, ax)
    let resultPack = z.requireMultiValue("axexpmulyPack packed backward result")
    let (resultExpax, _) = symbolicSlots[
      tuple[expax: Ggauge, value: Ggauge]](resultPack, "axexpmulyPack node")
    return multiValues(
      "axexpmulyPack packed input gradients",
      deriv.redot x,
      a * deriv,
      resultExpax.adjmul upstreamValue)

proc axexpmulyPackedInputPackf(v: Gvalue) =
  let view = v.requireUnaryNodeView(Gmulti, "axexpmulyPack packed forward", "args")
  # Args: [a, x, y]
  let (a, x, y) = storedSlots[
    tuple[a: Gscalar, x: Ggauge, y: Ggauge]](
      view.x,
      "axexpmulyPack packed forward")
  let pack = v.requireMultiValue("axexpmulyPack packed forward result")
  # Result slots: [exp(a*x), exp(a*x)*y]
  let (expax, value) = storedSlots[
    tuple[expax: Ggauge, value: Ggauge]](pack, "axexpmulyPack value")
  let f = a.sval
  threads:
    for mu in 0..<value.unsafeGaugeStorage.len:
      for e in value.unsafeGaugeStorage[mu]:
        var t{.noinit.}: evalType(x.unsafeGaugeStorage[mu][e])
        t := exp(f * x.unsafeGaugeStorage[mu][e])
        expax.unsafeGaugeStorage[mu][e] := t
        value.unsafeGaugeStorage[mu][e] := t * y.unsafeGaugeStorage[mu][e]

let axexpmulyPackedInputPackg = newGfunc(
  forward = axexpmulyPackedInputPackf,
  backward = axexpmulyPackedInputPackb,
  name = "axexpmulyPack packed")

proc axexpmulyPack(a: Gscalar, x: Ggauge, y: Ggauge): Gmulti =
  ## Input slots: [a, x, y]. Result slots: [exp(a*x), exp(a*x)*y].
  x.requireSameGaugeShape(y, "axexpmulyPack packed")
  let args = multiValues("axexpmuly args", a, x, y)
  newMultiOutputNode(
    @[Gvalue(x), Gvalue(x)],
    @[Gvalue(args)],
    axexpmulyPackedInputPackg,
    "axexpmulyPack packed")

proc axexpmuly*(a: Gscalar, x: Ggauge, y: Ggauge): Ggauge =
  let pack = axexpmulyPack(a, x, y)
  let (_, value) = symbolicSlots[
    tuple[expax: Ggauge, value: Ggauge]](pack, "axexpmulyPack node")
  result = value
