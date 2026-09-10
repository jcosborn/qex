import ../[core, scalar]
import ../support/op
import layout, gauge, physics/qcdTypes
import types

# Site-local algebra, stamped once for the gauge bundle and once for the
# single field. Backward hooks use the pairing dL = Re tr(G^dag dF): for
# z = x*y, x-bar = u*y^dag and y-bar = x^dag*u; retr, norm2, redot and
# projTAH are self-adjoint or have the unit/2x/partner adjoints below.

type Literal = int | float

template field(x: Ggauge, mu: int): untyped = x.gval[mu]
template field(x: GfieldOf, mu: int): untyped = x.fval
template nfields(x: Ggauge): int = x.gval.len
template nfields(x: GfieldOf): int = 1
template nodeLike(x: Ggauge): Ggauge = x.gaugeNodeLike
template nodeLike(x: GfieldOf): untyped = x.fieldNodeLike
template requireSameShape(x, y: Ggauge, label: string) = requireSameGaugeShape(x, y, label)
template requireSameShape(x, y: GfieldOf, label: string) = requireSameFieldShape(x, y, label)
template unitLike(x: Ggauge): Ggauge = x.unitGaugeLike
template unitLike(x: GfieldOf): untyped = unitField(x.runtime, x.fval)

template forFields(z: typed, body: untyped) =
  ## One threads region over the fields of z; mu injected.
  threads:
    for mu {.inject.} in 0 ..< z.nfields:
      body

template forElems(z: typed, body: untyped) =
  ## One threads region over every site of every field of z; mu, e injected.
  threads:
    for mu {.inject.} in 0 ..< z.nfields:
      for e {.inject.} in z.field(mu):
        body

# Hooks are generic in the node type; the stamp below binds them into
# operator records and constructors for each concrete type.

template ops =
  mixin field, nfields, nodeLike, unitLike, adj, retr, redot, projTAH, `*`, `+`, `-`

proc retrb[T](zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  ops
  scaledUpstreamOr(zb, Gscalar, T(z.inputs[0]).unitLike)

proc retrf[T](v: Gvalue) =
  ops
  let x = T(v.inputs[0])
  let z = Gscalar(v)
  threads:
    var t = 0.0
    for mu in 0 ..< x.nfields:
      t += x.field(mu).trace.re
    threadMaster: z.sval = t

proc adjb[T](zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  ops
  requireUpstream(zb, "adj backward", T).adj

proc adjf[T](v: Gvalue) =
  ops
  let x = T(v.inputs[0])
  let z = T(v)
  z.forFields:
    z.field(mu) := x.field(mu).adj

proc norm2b[T](zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  ops
  scaledUpstreamOr(zb, Gscalar, 2.0 * T(z.inputs[0]))

proc norm2f[T](v: Gvalue) =
  ops
  let x = T(v.inputs[0])
  let z = Gscalar(v)
  threads:
    var t = 0.0
    for mu in 0 ..< x.nfields:
      t += x.field(mu).norm2
    threadMaster: z.sval = t

proc negb[T](zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  ops
  -requireUpstream(zb, "neg backward", T)

proc negf[T](v: Gvalue) =
  ops
  let x = T(v.inputs[0])
  let z = T(v)
  z.forFields:
    z.field(mu) := -x.field(mu)

proc addsb[T](zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  ops
  let u = requireUpstream(zb, "s+ backward", T)
  if i == 0:
    return retr(u)
  u

proc addsf[T](v: Gvalue) =
  ops
  let x = Gscalar(v.inputs[0])
  let y = T(v.inputs[1])
  let z = T(v)
  z.forFields:
    z.field(mu) := x.sval + y.field(mu)

proc addb[T](zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  requireUpstream(zb, "add backward", T)

proc addf[T](v: Gvalue) =
  ops
  let z = T(v)
  z.forElems:
    var t {.noinit.}: evalType(z.field(mu)[e])
    t := T(v.inputs[0]).field(mu)[e]
    for i in 1 ..< v.inputs.len:
      t += T(v.inputs[i]).field(mu)[e]
    z.field(mu)[e] := t

proc mulsb[T](zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  ops
  let x = Gscalar(z.inputs[0])
  let y = T(z.inputs[1])
  let u = requireUpstream(zb, "s* backward", T)
  if i == 0:
    return redot(u, y)
  x * u

proc mulsf[T](v: Gvalue) =
  ops
  let x = Gscalar(v.inputs[0])
  let y = T(v.inputs[1])
  let z = T(v)
  z.forFields:
    z.field(mu) := x.sval * y.field(mu)

proc mulb[T](zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  ops
  let x = T(z.inputs[0])
  let y = T(z.inputs[1])
  let u = requireUpstream(zb, "mul backward", T)
  if i == 0:
    return u * y.adj
  x.adj * u

proc mulf[T](v: Gvalue) =
  ops
  let x = T(v.inputs[0])
  let y = T(v.inputs[1])
  let z = T(v)
  z.forFields:
    z.field(mu) := x.field(mu) * y.field(mu)

proc redotb[T](zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  bilinearBackward(zb, z, i, T)

proc redotf[T](v: Gvalue) =
  ops
  let x = T(v.inputs[0])
  let y = T(v.inputs[1])
  let z = Gscalar(v)
  threads:
    var t = 0.0
    for mu in 0 ..< x.nfields:
      t += redot(x.field(mu), y.field(mu))
    threadMaster: z.sval = t

proc subsb[T](zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  ops
  let u = requireUpstream(zb, "-s backward", T)
  if i == 0:
    return u
  -retr(u)

proc subsf[T](v: Gvalue) =
  ops
  let x = T(v.inputs[0])
  let y = Gscalar(v.inputs[1])
  let z = T(v)
  z.forFields:
    z.field(mu) := x.field(mu) - y.sval

proc subb[T](zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  ops
  let u = requireUpstream(zb, "sub backward", T)
  if i == 0:
    return u
  -u

proc subf[T](v: Gvalue) =
  ops
  let x = T(v.inputs[0])
  let y = T(v.inputs[1])
  let z = T(v)
  z.forFields:
    z.field(mu) := x.field(mu) - y.field(mu)

proc projTAHb[T](zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  ops
  projTAH(requireUpstream(zb, "projTAH backward", T))

proc projTAHf[T](v: Gvalue) =
  ops
  let x = T(v.inputs[0])
  let z = T(v)
  z.forElems:
    z.field(mu)[e].projectTAH(x.field(mu)[e])

template stampSiteOps(T: typedesc, tag: static string) =
  proc retr*(x: T): Gscalar
  proc adj*(x: T): T
  proc norm2*(x: T): Gscalar
  proc redot*(x: T, y: T): Gscalar
  proc projTAH*(x: T): T
  proc `-`*(x: T): T
  proc `+`*(x: Gscalar, y: T): T
  proc `+`*(x: T, y: T): T
  proc `*`*(x: Gscalar, y: T): T
  proc `*`*(x: T, y: T): T
  proc `-`*(x: T, y: Gscalar): T
  proc `-`*(x: T, y: T): T

  proc `+`*[L: Literal](x: L, y: T): T =
    toGvalue(y.runtime, float(x)) + y

  proc `*`*[L: Literal](x: L, y: T): T =
    toGvalue(y.runtime, float(x)) * y

  proc `-`*[L: Literal](x: T, y: L): T =
    x - toGvalue(x.runtime, float(y))

  let retrg = Gfunc(forward: retrf[T], backward: retrb[T], name: "retr" & tag)
  let adjg = Gfunc(forward: adjf[T], backward: adjb[T], name: "adj" & tag)
  let norm2g = Gfunc(forward: norm2f[T], backward: norm2b[T], name: "norm2" & tag)
  let negg = Gfunc(forward: negf[T], backward: negb[T], name: "-" & tag)
  let addsg = Gfunc(forward: addsf[T], backward: addsb[T], name: "s+" & tag)
  let addg = Gfunc(forward: addf[T], backward: addb[T], name: tag & "+" & tag)
  let mulsg = Gfunc(forward: mulsf[T], backward: mulsb[T], name: "s*" & tag)
  let mulg = Gfunc(forward: mulf[T], backward: mulb[T], name: tag & "*" & tag)
  let redotg = Gfunc(forward: redotf[T], backward: redotb[T], name: "redot" & tag)
  let subsg = Gfunc(forward: subsf[T], backward: subsb[T], name: tag & "-s")
  let subg = Gfunc(forward: subf[T], backward: subb[T], name: tag & "-" & tag)
  let projTAHg = Gfunc(forward: projTAHf[T], backward: projTAHb[T], name: "projTAH" & tag)

  proc retr*(x: T): Gscalar =
    graphNode(scalarNodeLike(x), @[Gvalue(x)], retrg, "retr" & tag)

  proc adj*(x: T): T =
    graphNode(x.nodeLike, @[Gvalue(x)], adjg, "adj" & tag)

  proc norm2*(x: T): Gscalar =
    graphNode(scalarNodeLike(x), @[Gvalue(x)], norm2g, "norm2" & tag)

  proc `-`*(x: T): T =
    graphNode(x.nodeLike, @[Gvalue(x)], negg, "-" & tag)

  proc `+`*(x: Gscalar, y: T): T =
    graphNode(y.nodeLike, @[Gvalue(x), Gvalue(y)], addsg, "s+" & tag)

  proc addTerms*(x: T): seq[T] =
    ## The leaves of a sum in evaluation order.
    if x.gfunc == addg:
      for input in x.inputs:
        result.add addTerms(T(input))
    else:
      result.add x

  proc `+`*(x: T, y: T): T =
    ## n-ary: a sum of sums flattens into one node.
    x.requireSameShape(y, tag & "+" & tag)
    var inputs: seq[Gvalue]
    if x.gfunc == addg:
      inputs.add x.inputs
    else:
      inputs.add Gvalue(x)
    inputs.add Gvalue(y)
    graphNode(x.nodeLike, inputs, addg, tag & "+" & tag)

  method addLike*(prototype: T, x: Gvalue, y: Gvalue): Gvalue =
    T(x) + T(y)

  proc `*`*(x: Gscalar, y: T): T =
    graphNode(y.nodeLike, @[Gvalue(x), Gvalue(y)], mulsg, "s*" & tag)

  proc `*`*(x: T, y: T): T =
    x.requireSameShape(y, tag & "*" & tag)
    graphNode(x.nodeLike, @[Gvalue(x), Gvalue(y)], mulg, tag & "*" & tag)

  method scaleLike*(contribution: T, upstream: Gvalue): Gvalue =
    if upstream of Gscalar:
      return Gscalar(upstream) * contribution
    if upstream of T:
      return T(upstream) * contribution
    raiseValueError(tag & " scale upstream expects a scalar or same-type value, got:\n" & upstream.nodeRepr)

  proc redot*(x: T, y: T): Gscalar =
    x.requireSameShape(y, "redot" & tag)
    graphNode(scalarNodeLike(x), @[Gvalue(x), Gvalue(y)], redotg, "redot" & tag)

  proc `-`*(x: T, y: Gscalar): T =
    graphNode(x.nodeLike, @[Gvalue(x), Gvalue(y)], subsg, tag & "-s")

  proc `-`*(x: T, y: T): T =
    x.requireSameShape(y, tag & "-" & tag)
    graphNode(x.nodeLike, @[Gvalue(x), Gvalue(y)], subg, tag & "-" & tag)

  proc projTAH*(x: T): T =
    graphNode(x.nodeLike, @[Gvalue(x)], projTAHg, "projTAH" & tag)

stampSiteOps(Ggauge, "g")
stampSiteOps(Gfield, "f")
when not cfieldIsGfield:
  stampSiteOps(Gcfield, "c")

proc blendSubset*(parity, dir: int, cand, x: Ggauge): Ggauge =
  ## Use `cand` on one parity/direction subset and `x` elsewhere.
  requireParityDir(parity, dir, x.gval.len, "blendSubset")
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

proc maskSubset*(parity, dir: int, x: Ggauge): Ggauge =
  ## x on the (parity, dir) subset, zero elsewhere.
  blendSubset(parity, dir, x, Ggauge(x.zeroLike))
