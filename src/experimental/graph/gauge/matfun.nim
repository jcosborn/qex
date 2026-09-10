## Matrix exponential family: fused kernels at every order of the exp tower.
##
## Under redot(a, b) = Re tr(a^dag b), with P the real-coefficient kernel
## polynomial (matexp poly4 with scale-20 squaring), every derivative of exp
## is a coefficient of P at a jet point:
##   expTop(y; d_1..d_m) = d/de_1 ... d/de_m P(y + sum_i e_i d_i)   symmetric in d
##   exp(x)          = expTop(x)
##   expDeriv(b, x)  = dP_x^*[b] = dP_{x^dag}[b] = expTop(x^dag; b)
##   z = expTop(y; d), cotangent v:
##     d_i-bar = expTop(y^dag; v, d_j^dag for j != i)       order m
##     y-bar   = expTop(y^dag; v, d_1^dag, ..., d_m^dag)     order m + 1
## The kernel maths/matrixFunctions.expTop covers m <= 3 (66, 198, 594
## matrix products per site); past that expTopReplica nests the basic-op
## polynomial replica, so no order raises. Nc = 1 uses the exact scalar
## identities instead.

import ../[core, scalar]
import ../support/op
import layout, gauge, physics/qcdTypes
from maths/matexp import newExpParam, ekPoly
import shared, basic_ops

proc exp*(x: Ggauge): Ggauge
proc expDeriv*(b: Ggauge, x: Ggauge, parity = -1, dir = 0): Ggauge
proc expJet*(y: Ggauge, d: openArray[Ggauge]): Ggauge

proc expgb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let x = Ggauge(z.inputs[0])
  expDeriv(requireUpstream(zb, "expg backward", Ggauge), x)

proc expgf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let z = Ggauge(v)
  z.mapGaugeElements:
    z.gval[mu][e] := exp(x.gval[mu][e])

let expg = Gfunc(forward: expgf, backward: expgb, name: "expg")

proc exp*(x: Ggauge): Ggauge =
  graphNode(x.gaugeNodeLike, @[Gvalue(x)], expg, "expg")

proc expPolyGraph*(x: Ggauge): Ggauge =
  ## Grad-complete graph replica of the site-local matrix exponential kernel
  ## (matexp ekPoly order-4 with 2^-scale scaling and repeated squaring):
  ##   y = x/2^scale
  ##   e = (y^2/24 + y/6 + 1/2)*y^2 + y      # expm1Poly4(y)
  ##   e <- e*(e+2), scale times             # (1+e)^2 = 1 + e*(e+2)
  ##   exp(x) ~ 1 + e
  ## For Nc > 1 this is the same scheme as the optimized `exp` kernel, so the
  ## two agree to roundoff. Built only from grad-complete ops, so it can be
  ## differentiated to any order.
  const params = newExpParam()
  static:
    doAssert params.kind == ekPoly and params.order == 4,
      "expPolyGraph must match the optimized matrix exponential"
  const scale = params.scale
  let y = (1.0 / float(1 shl scale)) * x
  let y2 = y * y
  let two = toGvalue(x.runtime, 2.0)
  var e = (0.5 + ((1.0/24.0) * y2 + (1.0/6.0) * y)) * y2 + y
  for _ in 1..scale:
    e = e * (two + e)
  result = 1.0 + e

proc expTopReplica*(y: Ggauge, d: openArray[Ggauge]): Ggauge =
  ## expTop(y; d) as a basic-op graph. For m = 0 it is the polynomial replica;
  ## otherwise the y-bar rule read backwards, with one fewer direction:
  ##   expTop(y; d_1..d_m) = gradSeeded(expTop(slot; d_2^dag..d_m^dag), slot, d_1),
  ##   slot = slotVar(y^dag).
  ## The fallback past m = 3, and the test oracle for the fused nodes.
  if d.len == 0:
    return expPolyGraph(y)
  let slot = slotVar(y.adj)
  var rest: seq[Ggauge]
  for j in 1 ..< d.len:
    rest.add d[j].adj
  Ggauge(gradSeeded(expTopReplica(slot, rest), slot, d[0]))

template expJetKernel(v: Gvalue, M: static int) =
  let y = Ggauge(v.inputs[0])
  let z = Ggauge(v)
  z.mapGaugeElements:
    var yy, r {.noinit.}: evalType(y.gval[mu][e])
    var dd {.noinit.}: array[M, evalType(y.gval[mu][e])]
    yy := y.gval[mu][e]
    for k in 0 ..< M:
      dd[k] := Ggauge(v.inputs[k + 1]).gval[mu][e]
    when M == 1:
      r[] := expTop(yy[], [dd[0][]])
    elif M == 2:
      r[] := expTop(yy[], [dd[0][], dd[1][]])
    else:
      r[] := expTop(yy[], [dd[0][], dd[1][], dd[2][]])
    z.gval[mu][e] := r

proc expJetf1(v: Gvalue) = expJetKernel(v, 1)
proc expJetf2(v: Gvalue) = expJetKernel(v, 2)
proc expJetf3(v: Gvalue) = expJetKernel(v, 3)

proc expJetb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  # One rule for every slot: the cotangent joins the directions and the
  # slot's own direction leaves; slot 0 (y) keeps them all, so its order is m+1.
  let y = Ggauge(z.inputs[0])
  var ds = @[requireUpstream(zb, "expJet backward", Ggauge)]
  for j in 1 ..< z.inputs.len:
    if j != i:
      ds.add Ggauge(z.inputs[j]).adj
  expJet(y.adj, ds)

let expJetg = [
  Gfunc(forward: expJetf1, backward: expJetb, name: "expJet1"),
  Gfunc(forward: expJetf2, backward: expJetb, name: "expJet2"),
  Gfunc(forward: expJetf3, backward: expJetb, name: "expJet3")]

proc expJet*(y: Ggauge, d: openArray[Ggauge]): Ggauge =
  ## expTop(y; d_1..d_m) as a graph node; fused kernel for m <= 3.
  var inputs = @[Gvalue(y)]
  for x in d:
    y.requireSameGaugeShape(x, "expJet")
    inputs.add Gvalue(x)
  if d.len < 1 or d.len > 3:
    return expTopReplica(y, d)
  graphNode(y.gaugeNodeLike, inputs, expJetg[d.len - 1], "expJet")

proc expDerivContribution(u: Ggauge, z: Gvalue, i: int,
                          parity = -1, dir = 0): Gvalue =
  ## Backward of z = expDeriv(b, x) (the pullback of exp at x with cotangent
  ## b) under the raw upstream u; (parity, dir) select the subset variant.
  let b = Ggauge(z.inputs[0])
  let x = Ggauge(z.inputs[1])
  if i == 0:
    # z is linear in b with adjoint the pushforward of exp at x. For Nc>1,
    # the kernel is a real-coefficient matrix polynomial P, so
    # (dP_x)^adj = dP_{x^dag}; scalar exp has the same adjoint identity.
    # Since exp is site-local, the subset kernel needs no masked upstream.
    return expDeriv(u, x.adj, parity, dir)
  # Second derivative of exp. <u, mask(w)> == <mask(u), w> masks the upstream
  # for the subset variant.
  let useed = if parity < 0: u else: maskSubset(parity, dir, u)
  const nc = x.gval[0][0].nrows
  when nc == 1:
    # The Nc=1 kernel is scalar exp with pullback exp(x^dag)*b; differentiate
    # that exact expression. It is already partial in x and keeps b live.
    useed.adj * exp(x.adj) * b
  else:
    # z = expTop(x^dag; b): the y-bar rule gives expTop(x; u, b^dag), and the
    # adjoint of that is the x slot: x-bar = expTop(x^dag; u^dag, b).
    expJet(x.adj, [useed.adj, b])

proc expDerivgb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  expDerivContribution(
    requireUpstream(zb, "expDeriv backward", Ggauge), z, i)

proc expDerivgf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let y = Ggauge(v.inputs[1])
  let z = Ggauge(v)
  z.mapGaugeElements:
    z.gval[mu][e] := expDeriv(
      y.gval[mu][e],
      x.gval[mu][e])

let expDerivg = Gfunc(forward: expDerivgf, backward: expDerivgb, name: "expDerivg")

proc expDeriv*(b: Ggauge, x: Ggauge, parity = -1, dir = 0): Ggauge =
  ## D exp(x)^*[b], on the whole field or only (parity,dir); zero elsewhere.
  if parity != -1:
    requireParityDir(parity, dir, b.gval.len, "expDeriv")
  let node = sameShapeGaugeNodeLike(b, x, "expDerivg")
  if parity < 0:
    return graphNode(node, @[Gvalue(b), Gvalue(x)], expDerivg, "expDerivg")
  let sub = b.gval.paritySubset(parity)
  proc fwd(v: Gvalue) =
    let b = Ggauge(v.inputs[0])
    let x = Ggauge(v.inputs[1])
    let z = Ggauge(v)
    forGaugeSubset(sub):
      # element order matches whole-field expDerivgf: expDeriv(inputs[1], inputs[0])
      z.gval[dir][e] := expDeriv(x.gval[dir][e], b.gval[dir][e])
  proc bwd(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
    expDerivContribution(
      requireUpstream(zb, "expDeriv subset backward", Ggauge), z, i, parity, dir)
  result = graphNode(node, @[Gvalue(b), Gvalue(x)],
    Gfunc(forward: fwd, backward: bwd, name: "expDerivg"),
    "expDerivg")
  result.zeroGaugeStorage
