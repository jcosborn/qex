import ../[core, scalar, multi]
import ../support/op
import layout, physics/qcdTypes
import shared, basic_ops

# Section: Fused Gauge Ops

proc muladj*(x: Ggauge, y: Ggauge): Ggauge

proc adjmulggb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let x = Ggauge(z.inputs[0])
  let y = Ggauge(z.inputs[1])
  let upstream = requireUpstream(zb, "g.adj*g backward", Ggauge)
  # adjmul(x, y) == x.adj * y
  if i == 0:
    return y.muladj(upstream)
  x * upstream

proc adjmulggf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let y = Ggauge(v.inputs[1])
  let z = Ggauge(v)
  z.mapGaugeSites(x.gval[mu].adj * y.gval[mu])

let adjmulgg = Gfunc(forward: adjmulggf, backward: adjmulggb, name: "g.adj*g")

proc adjmul*(x: Ggauge, y: Ggauge): Ggauge =
  graphNode(sameShapeGaugeNodeLike(x, y, "g.adj*g"), @[Gvalue(x), Gvalue(y)], adjmulgg, "g.adj*g")

proc muladjggb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let x = Ggauge(z.inputs[0])
  let y = Ggauge(z.inputs[1])
  let upstream = requireUpstream(zb, "g*g.adj backward", Ggauge)
  # muladj(x, y) == x * y.adj
  if i == 0:
    return upstream * y
  upstream.adjmul(x)

proc muladjggf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let y = Ggauge(v.inputs[1])
  let z = Ggauge(v)
  z.mapGaugeSites(x.gval[mu] * y.gval[mu].adj)

let muladjgg = Gfunc(forward: muladjggf, backward: muladjggb, name: "g*g.adj")

proc muladj*(x: Ggauge, y: Ggauge): Ggauge =
  graphNode(sameShapeGaugeNodeLike(x, y, "g*g.adj"), @[Gvalue(x), Gvalue(y)], muladjgg, "g*g.adj")

proc contractProjTAHPackedInputb(zb: Gvalue,
                                 z: Gvalue,
                                 i: int,
                                 input: Gvalue): Gvalue =
  let args = Gmulti(z.inputs[0])
  # Args: [x, y]
  let x = Ggauge(args[0])
  let y = Ggauge(args[1])
  let proj = projTAH(
    requireUpstream(zb, "contractProjTAH packed backward", Ggauge))
  multiValues(
    "contractProjTAH packed input gradients",
    proj * y,
    proj.adjmul x)

proc contractProjTAHPackedInputf(v: Gvalue) =
  let args = Gmulti(v.inputs[0])
  # Args: [x, y]
  let x = Ggauge(args.storedSlot(0))
  let y = Ggauge(args.storedSlot(1))
  let z = Ggauge(v)
  z.mapGaugeElements:
    let s = x.gval[mu][e] * y.gval[mu][e].adj
    z.gval[mu][e].projectTAH s

let contractProjTAHPackedInputg = Gfunc(
  forward: contractProjTAHPackedInputf,
  backward: contractProjTAHPackedInputb,
  name: "contractProjTAH packed")

proc contractProjTAH*(x: Ggauge, y: Ggauge): Ggauge =
  x.requireSameGaugeShape(y, "contractProjTAH packed")
  let args = multiValues("contractProjTAH args", x, y)
  graphNode(
    x.gaugeNodeLike,
    @[Gvalue(args)],
    contractProjTAHPackedInputg,
    "contractProjTAH packed")

proc axexpPackedInputb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let args = Gmulti(z.inputs[0])
  # Args: [a, x]
  let a = Gscalar(args[0])
  let x = Ggauge(args[1])
  let ax = a * x
  let deriv = expDeriv(
    requireUpstream(zb, "axexp packed backward", Ggauge),
    ax)
  multiValues(
    "axexp packed input gradients",
    deriv.redot x,
    a * deriv)

proc axexpPackedInputf(v: Gvalue) =
  let args = Gmulti(v.inputs[0])
  # Args: [a, x]
  let a = Gscalar(args.storedSlot(0))
  let x = Ggauge(args.storedSlot(1))
  let z = Ggauge(v)
  let f = a.sval
  z.mapGaugeElements:
    z.gval[mu][e] := exp(f * x.gval[mu][e])

let axexpPackedInputg = Gfunc(
  forward: axexpPackedInputf,
  backward: axexpPackedInputb,
  name: "axexp packed")

proc axexp*(a: Gscalar, x: Ggauge): Ggauge =
  let args = multiValues("axexp args", a, x)
  graphNode(x.gaugeNodeLike, @[Gvalue(args)], axexpPackedInputg, "axexp packed")

proc axexpmulyPackedInputPackb(zb: Gvalue,
                               z: Gvalue,
                               i: int,
                               input: Gvalue): Gvalue =
  let args = Gmulti(z.inputs[0])
  # Args: [a, x, y]
  let a = Gscalar(args[0])
  let x = Ggauge(args[1])
  let y = Ggauge(args[2])
  let zb = requireUpstream(zb, "axexpmulyPack packed backward", Gmulti)
  # Upstream/result slots: [exp(a*x), exp(a*x)*y]
  let upstreamExpax = Ggauge(zb[0])
  let upstreamValue = Ggauge(zb[1])
  let ax = a * x
  let deriv = expDeriv(upstreamExpax + upstreamValue.muladj y, ax)
  let resultPack = Gmulti(z)
  let resultExpax = Ggauge(resultPack[0])
  multiValues(
    "axexpmulyPack packed input gradients",
    deriv.redot x,
    a * deriv,
    resultExpax.adjmul upstreamValue)

proc axexpmulyPackedInputPackf(v: Gvalue) =
  let args = Gmulti(v.inputs[0])
  # Args: [a, x, y]
  let a = Gscalar(args.storedSlot(0))
  let x = Ggauge(args.storedSlot(1))
  let y = Ggauge(args.storedSlot(2))
  let pack = Gmulti(v)
  # Result slots: [exp(a*x), exp(a*x)*y]
  let expax = Ggauge(pack.storedSlot(0))
  let value = Ggauge(pack.storedSlot(1))
  let f = a.sval
  threads:
    for mu in 0..<value.gval.len:
      for e in value.gval[mu]:
        var t{.noinit.}: evalType(x.gval[mu][e])
        t := exp(f * x.gval[mu][e])
        expax.gval[mu][e] := t
        value.gval[mu][e] := t * y.gval[mu][e]

let axexpmulyPackedInputPackg = Gfunc(
  forward: axexpmulyPackedInputPackf,
  backward: axexpmulyPackedInputPackb,
  name: "axexpmulyPack packed")

proc axexpmuly*(a: Gscalar, x: Ggauge, y: Ggauge): Ggauge =
  ## Input slots: [a, x, y]; packed result slots: [exp(a*x), exp(a*x)*y]. Returns
  ## slot 1; the fused carrier shares the exp(a*x) work across both outputs.
  x.requireSameGaugeShape(y, "axexpmulyPack packed")
  let args = multiValues("axexpmuly args", a, x, y)
  let pack = newMultiOutputNode(
    @[Gvalue(x), Gvalue(x)],
    @[Gvalue(args)],
    axexpmulyPackedInputPackg,
    "axexpmulyPack packed")
  Ggauge(pack[1])
