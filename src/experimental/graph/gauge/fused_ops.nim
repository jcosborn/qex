import ../[core, scalar, multi]
import ../support/op
import layout, physics/qcdTypes
import gauge/gaugeUtils
import types, basic_ops, matfun

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

let axpyg = Gfunc(bufferMode: bmFull, forward: axpygf, backward: axpygb, name: "axpy")

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

let adjmulgg = Gfunc(bufferMode: bmFull, forward: adjmulggf, backward: adjmulggb, name: "g.adj*g")

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

let muladjgg = Gfunc(bufferMode: bmFull, forward: muladjggf, backward: muladjggb, name: "g*g.adj")

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
  threads:
    contractProjectTAH(z.gval, x.gval, y.gval)

let contractProjTAHPackedInputg = Gfunc(bufferMode: bmFull,
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
    let x = gaugeTermSum(v.inputs.len-2,i,Ggauge(v.inputs[2+i]).gval[mu][e])
    let s = x * y.gval[mu][e].adj
    z.gval[mu][e].projectTAH s

let contractProjTAHSumInputg = Gfunc(bufferMode: bmFull,
  forward: contractProjTAHSumInputf,
  backward: contractProjTAHPackedInputb,
  inputView: contractProjTAHSumInputView,
  name: "contractProjTAH packed")

proc contractProjTAH*(x: Ggauge, y: Ggauge, parity = -1, dir = 0): Ggauge =
  ## projectTAH(x*y†), on the whole field or only (parity,dir).
  x.requireSameGaugeShape(y, "contractProjTAH packed")
  let args = multiValues("contractProjTAH args", x, y)
  if parity < 0:
    let terms = x.addTerms
    if terms.len > 1:
      var inputs = @[Gvalue(args), Gvalue(y)]
      for term in terms:
        inputs.add Gvalue(term)
      return graphNode(x.gaugeNodeLike, inputs, contractProjTAHSumInputg, "contractProjTAH packed")
    return graphNode(
      x.gaugeNodeLike, @[Gvalue(args)], contractProjTAHPackedInputg,
      "contractProjTAH packed")
  let sub = x.gval.paritySubset(parity)
  let terms = x.addTerms
  if terms.len > 1:
    var inputs = @[Gvalue(args), Gvalue(y)]
    for term in terms:
      inputs.add Gvalue(term)
    proc sumf(v: Gvalue) =
      let
        y = Ggauge(v.inputs[1])
        z = Ggauge(v)
      forGaugeSubset(sub):
        let x = gaugeTermSum(v.inputs.len-2,i,Ggauge(v.inputs[2+i]).gval[dir][e])
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
      Gfunc(bufferMode: bmZero, forward: sumf, backward: sumb, inputView: contractProjTAHSumInputView,
            name: "contractProjTAH packed"),
      "contractProjTAH packed")
    return
  proc fwd(v: Gvalue) =
    let args = Gmulti(v.inputs[0])
    let x = Ggauge(args.storedSlot(0))
    let y = Ggauge(args.storedSlot(1))
    let z = Ggauge(v)
    threads:
      contractProjectTAH(z.gval[dir], x.gval[dir], y.gval[dir], sub)
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
    Gfunc(bufferMode: bmZero, forward: fwd, backward: bwd, name: "contractProjTAH packed"),
    "contractProjTAH packed")

template axexpNodes(axexp, axexpmuly, kexp, kexpmuly: untyped; E: static string) =
  ## Graph nodes E(a x) and E(a x) y on the field kernels kexp and kexpmuly.
  ## Both families pull back through the exp tower's expDeriv.
  proc packedb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue {.gensym.} =
    let args = Gmulti(z.inputs[0])
    # Args: [a, x]
    let a = Gscalar(args[0])
    let x = Ggauge(args[1])
    let ax = a * x
    let deriv = expDeriv(
      requireUpstream(zb, E & " packed backward", Ggauge),
      ax)
    multiValues(
      E & " packed input gradients",
      deriv.redot x,
      a * deriv)

  proc packedf(v: Gvalue) {.gensym.} =
    let args = Gmulti(v.inputs[0])
    # Args: [a, x]
    let a = Gscalar(args.storedSlot(0))
    let x = Ggauge(args.storedSlot(1))
    let z = Ggauge(v)
    threads:
      kexp(z.gval, a.sval, x.gval)

  let packedg = Gfunc(bufferMode: bmFull,
    forward: packedf,
    backward: packedb,
    name: E & " packed")

  proc axexp*(a: Gscalar, x: Ggauge): Ggauge =
    let args = multiValues(E & " args", a, x)
    graphNode(x.gaugeNodeLike, @[Gvalue(args)], packedg, E & " packed")

  proc packb(zb: Gvalue,
             z: Gvalue,
             i: int,
             input: Gvalue): Gvalue {.gensym.} =
    let args = Gmulti(z.inputs[0])
    # Args: [a, x, y]
    let a = Gscalar(args[0])
    let x = Ggauge(args[1])
    let y = Ggauge(args[2])
    let zb = requireUpstream(zb, E & "mulyPack packed backward", Gmulti)
    # Upstream/result slots: [E(a*x), E(a*x)*y]
    let upstreamExpax = Ggauge(zb[0])
    let upstreamValue = Ggauge(zb[1])
    let ax = a * x
    let deriv = expDeriv(upstreamExpax + upstreamValue.muladj y, ax)
    let resultPack = Gmulti(z)
    let resultExpax = Ggauge(resultPack[0])
    multiValues(
      E & "mulyPack packed input gradients",
      deriv.redot x,
      a * deriv,
      resultExpax.adjmul upstreamValue)

  proc packf(v: Gvalue) {.gensym.} =
    let args = Gmulti(v.inputs[0])
    # Args: [a, x, y]
    let a = Gscalar(args.storedSlot(0))
    let x = Ggauge(args.storedSlot(1))
    let y = Ggauge(args.storedSlot(2))
    let pack = Gmulti(v)
    # Result slots: [E(a*x), E(a*x)*y]
    let expax = Ggauge(pack.storedSlot(0))
    let value = Ggauge(pack.storedSlot(1))
    threads:
      kexpmuly(value.gval, a.sval, x.gval, y.gval, expax.gval)

  let packg = Gfunc(bufferMode: bmFull,
    forward: packf,
    backward: packb,
    name: E & "mulyPack packed")

  proc axexpmuly*(a: Gscalar, x: Ggauge, y: Ggauge, parity = -1, dir = 0): Ggauge =
    ## E(a*x)*y; cache E(a*x) in the same carrier.
    ## For parity >= 0, compute E and its pullback only on (parity,dir).
    x.requireSameGaugeShape(y, E & "mulyPack packed")
    let args = multiValues(E & "muly args", a, x, y)
    if parity < 0:
      let pack = newMultiOutputNode(
        @[Gvalue(x), Gvalue(x)], @[Gvalue(args)],
        packg, E & "mulyPack packed")
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
      threads:
        kexpmuly(value.gval[dir], a.sval, x.gval[dir], y.gval[dir], sub, expax.gval[dir])
    proc bwd(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
      let args = Gmulti(z.inputs[0])
      let a = Gscalar(args[0])
      let x = Ggauge(args[1])
      let y = Ggauge(args[2])
      let zbm = requireUpstream(zb, E & "mulyPack packed backward", Gmulti)
      let upstreamExpax = Ggauge(zbm[0])
      let upstreamValue = Ggauge(zbm[1])
      let ax = a * x
      let deriv = expDeriv(upstreamExpax + upstreamValue.muladj y, ax, parity, dir)
      let resultExpax = Ggauge(Gmulti(z)[0])
      multiValues(
        E & "mulyPack packed input gradients",
        deriv.redot x,
        a * deriv,
        resultExpax.adjmul upstreamValue)
    let pack = newMultiOutputNode(
      @[Gvalue(x), Gvalue(x)], @[Gvalue(args)],
      Gfunc(bufferMode: bmZero, forward: fwd, backward: bwd, name: E & "mulyPack packed"),
      E & "mulyPack packed")
    Ggauge(pack[1])

# axexp/axexpmuly: any x, with the general field kernels (degree 4, scale 20).
axexpNodes(axexp, axexpmuly, axexp, axexpmuly, "axexp")
# axexpah/axexpahmuly: x anti-Hermitian, and traceless for SU(N), with expAH.
axexpNodes(axexpah, axexpahmuly, axexpah, axexpahmuly, "axexpah")
