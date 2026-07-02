import ../[core, scalar, multi]
import ../support/op
import layout, physics/qcdTypes
import shared, basic_ops

# Section: Fused Gauge Ops

proc muladj*(x: Ggauge, y: Ggauge): Ggauge

proc axpygb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let a = Gscalar(z.inputs[0])
  let x = Ggauge(z.inputs[1])
  let upstream = requireUpstream(zb, "axpy backward", Ggauge)
  case i
  of 0:
    upstream.redot x
  of 1:
    a * upstream
  else:
    upstream

proc axpygf(v: Gvalue) =
  let a = Gscalar(v.inputs[0])
  let x = Ggauge(v.inputs[1])
  let y = Ggauge(v.inputs[2])
  let z = Ggauge(v)
  z.mapGaugeSites(a.sval * x.gval[mu] + y.gval[mu])

let axpyg = Gfunc(forward: axpygf, backward: axpygb, name: "axpy")

proc axpy*(a: Gscalar, x, y: Ggauge): Ggauge =
  ## Fused `a*x + y`, evaluated in one whole-gauge pass.
  graphNode(sameShapeGaugeNodeLike(x, y, "axpy"), @[Gvalue(a), Gvalue(x), Gvalue(y)], axpyg, "axpy")

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

proc blendSubset*(parity, dir: int, cand, x: Ggauge): Ggauge =
  ## Use `cand` on one parity/direction subset and `x` elsewhere.
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

proc contractProjTAHSumInputView(v: Gvalue, mode: InputWalkMode, visit: GnodeVisit) =
  case mode
  of iwmEval:
    for i in 1..<v.inputs.len:
      visit v.inputs[i]
  of iwmReachable, iwmBackward:
    visit v.inputs[0]

proc contractProjTAHSumInputf(v: Gvalue) =
  let
    y = Ggauge(v.inputs[1])
    z = Ggauge(v)
  z.mapGaugeElements:
    var x {.noinit.}: evalType(y.gval[mu][e])
    x := Ggauge(v.inputs[2]).gval[mu][e]
    for i in 3..<v.inputs.len:
      x += Ggauge(v.inputs[i]).gval[mu][e]
    let s = x * y.gval[mu][e].adj
    z.gval[mu][e].projectTAH s

let contractProjTAHSumInputg = Gfunc(
  forward: contractProjTAHSumInputf,
  backward: contractProjTAHPackedInputb,
  inputView: contractProjTAHSumInputView,
  name: "contractProjTAH packed")

proc contractProjTAH*(x: Ggauge, y: Ggauge, parity = -1, dir = 0): Ggauge =
  ## projectTAH(x*y†), on the whole field or only (parity,dir).
  x.requireSameGaugeShape(y, "contractProjTAH packed")
  let args = multiValues("contractProjTAH args", x, y)
  if parity < 0:
    let terms = x.gaugeAddTerms
    if terms.len > 1:
      var inputs = @[Gvalue(args), Gvalue(y)]
      for term in terms:
        inputs.add Gvalue(term)
      return graphNode(x.gaugeNodeLike, inputs, contractProjTAHSumInputg, "contractProjTAH packed")
    return graphNode(
      x.gaugeNodeLike, @[Gvalue(args)], contractProjTAHPackedInputg,
      "contractProjTAH packed")
  let sub = x.gval.paritySubset(parity)
  let terms = x.gaugeAddTerms
  if terms.len > 1:
    var inputs = @[Gvalue(args), Gvalue(y)]
    for term in terms:
      inputs.add Gvalue(term)
    proc sumf(v: Gvalue) =
      let
        y = Ggauge(v.inputs[1])
        z = Ggauge(v)
      forGaugeSubset(sub):
        var x {.noinit.}: evalType(y.gval[dir][e])
        x := Ggauge(v.inputs[2]).gval[dir][e]
        for i in 3..<v.inputs.len:
          x += Ggauge(v.inputs[i]).gval[dir][e]
        let s = x * y.gval[dir][e].adj
        z.gval[dir][e].projectTAH s
    proc sumb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
      let args = Gmulti(z.inputs[0])
      let
        x = Ggauge(args[0])
        y = Ggauge(args[1])
        up = requireUpstream(zb, "contractProjTAH subset backward", Ggauge)
        zero = Ggauge(up.zeroLike)
        proj = projTAH(blendSubset(parity, dir, up, zero))
      multiValues("contractProjTAH subset input gradients", proj * y, proj.adjmul x)
    result = graphNode(
      x.gaugeNodeLike, inputs,
      Gfunc(forward: sumf, backward: sumb, inputView: contractProjTAHSumInputView,
            name: "contractProjTAH packed"),
      "contractProjTAH packed")
    result.zeroGaugeStorage
    return
  proc fwd(v: Gvalue) =
    let args = Gmulti(v.inputs[0])
    let x = Ggauge(args.storedSlot(0))
    let y = Ggauge(args.storedSlot(1))
    let z = Ggauge(v)
    forGaugeSubset(sub):
      let s = x.gval[dir][e] * y.gval[dir][e].adj
      z.gval[dir][e].projectTAH s
  proc bwd(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
    let args = Gmulti(z.inputs[0])
    let
      x = Ggauge(args[0])
      y = Ggauge(args[1])
      up = requireUpstream(zb, "contractProjTAH subset backward", Ggauge)
      zero = Ggauge(up.zeroLike)
      proj = projTAH(blendSubset(parity, dir, up, zero))
    multiValues("contractProjTAH subset input gradients", proj * y, proj.adjmul x)
  result = graphNode(
    x.gaugeNodeLike, @[Gvalue(args)],
    Gfunc(forward: fwd, backward: bwd, name: "contractProjTAH packed"),
    "contractProjTAH packed")
  result.zeroGaugeStorage

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
        t[] := expAH(f * x.gval[mu][e][])
        expax.gval[mu][e] := t
        value.gval[mu][e] := t * y.gval[mu][e]

let axexpmulyPackedInputPackg = Gfunc(
  forward: axexpmulyPackedInputPackf,
  backward: axexpmulyPackedInputPackb,
  name: "axexpmulyPack packed")

proc axexpmuly*(a: Gscalar, x: Ggauge, y: Ggauge, parity = -1, dir = 0): Ggauge =
  ## exp(a*x)*y; cache exp(a*x) in the same carrier.
  ## For parity >= 0, compute exp and its pullback only on (parity,dir).
  ## x is anti-Hermitian, and traceless for SU(N).
  x.requireSameGaugeShape(y, "axexpmulyPack packed")
  let args = multiValues("axexpmuly args", a, x, y)
  if parity < 0:
    let pack = newMultiOutputNode(
      @[Gvalue(x), Gvalue(x)], @[Gvalue(args)],
      axexpmulyPackedInputPackg, "axexpmulyPack packed")
    return Ggauge(pack[1])
  let sub = x.gval.paritySubset(parity)
  proc fwd(v: Gvalue) =
    let args = Gmulti(v.inputs[0])
    let a = Gscalar(args.storedSlot(0))
    let x = Ggauge(args.storedSlot(1))
    let y = Ggauge(args.storedSlot(2))
    let pack = Gmulti(v)
    let expax = Ggauge(pack.storedSlot(0))
    let value = Ggauge(pack.storedSlot(1))
    let f = a.sval
    forGaugeSubset(sub):
      var t{.noinit.}: evalType(x.gval[dir][e])
      t[] := expAH(f * x.gval[dir][e][])
      expax.gval[dir][e] := t
      value.gval[dir][e] := t * y.gval[dir][e]
  proc bwd(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
    let args = Gmulti(z.inputs[0])
    let a = Gscalar(args[0])
    let x = Ggauge(args[1])
    let y = Ggauge(args[2])
    let zbm = requireUpstream(zb, "axexpmulyPack packed backward", Gmulti)
    let upstreamExpax = Ggauge(zbm[0])
    let upstreamValue = Ggauge(zbm[1])
    let ax = a * x
    let deriv = expDeriv(upstreamExpax + upstreamValue.muladj y, ax, parity, dir)
    let resultExpax = Ggauge(Gmulti(z)[0])
    multiValues(
      "axexpmulyPack packed input gradients",
      deriv.redot x,
      a * deriv,
      resultExpax.adjmul upstreamValue)
  let pack = newMultiOutputNode(
    @[Gvalue(x), Gvalue(x)], @[Gvalue(args)],
    Gfunc(forward: fwd, backward: bwd, name: "axexpmulyPack packed"),
    "axexpmulyPack packed")
  Ggauge(pack[0]).zeroGaugeStorage     # off-subset stays zero across evals
  Ggauge(pack[1]).zeroGaugeStorage
  Ggauge(pack[1])
