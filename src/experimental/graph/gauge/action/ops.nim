## Gauge-action graph values and operations.
import ../../[core, scalar]
import ../../scalar/types
import ../../support/op
import layout, gauge, physics/qcdTypes
import ../types, ../basic_ops, ../fused_ops, ../field_ops, ../transport, ../cfield, domain

# --- gauge-action coefficient value type and coefficient algebra ---

type
  Gactcoeff* {.final.} = ref object of Gvalue
    cval*: GaugeActionCoeffs

proc update*(x: Gactcoeff, c: GaugeActionCoeffs) =
  x.cval = c
  x.updated

proc toGvalue*(grt: GraphRuntime,
               x: GaugeActionCoeffs): Gactcoeff =
  result = Gactcoeff(
    runtime: grt,
    cval: x).assignStableNodeId
  result.updated

method newOneOf*(x: Gactcoeff): Gvalue =
  result = Gactcoeff(runtime: x.runtime).assignStableNodeId
method valCopy*(z: Gactcoeff, x: Gvalue) =
  z.cval = Gactcoeff(x).cval
method copyCompatible*(prototype: Gactcoeff, value: Gvalue): bool =
  value of Gactcoeff
method `$`*(x: Gactcoeff): string = $x.cval

proc raiseCoeffBackwardUnsupported(label: string) {.noreturn.} =
  raiseUnsupportedPath(label, "derivative with respect to gauge-action coefficients")

proc `*`*(x: Gscalar, y: Gactcoeff): Gactcoeff
proc `+`*(x, y: Gactcoeff): Gactcoeff
proc redot*(x: Gactcoeff, y: Gactcoeff): Gscalar

method addLike*(prototype: Gactcoeff, x: Gvalue, y: Gvalue): Gvalue =
  Gactcoeff(x) + Gactcoeff(y)

method scaleLike*(contribution: Gactcoeff, upstream: Gvalue): Gvalue =
  Gscalar(upstream) * contribution

proc adjCoeffb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let upstream = requireUpstream(zb, "adjCoeff backward", Gactcoeff)
  redot(upstream,
        toGvalue(z.runtime, GaugeActionCoeffs(adjplaq: 1.0)))

proc adjCoefff(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let z = Gactcoeff(v)
  z.cval = GaugeActionCoeffs(plaq: 1.0, adjplaq: x.sval)

let adjCoeffg = Gfunc(forward: adjCoefff, backward: adjCoeffb, name: "adjCoeff")

proc adjCoeff(adjFac: Gscalar): Gactcoeff =
  graphNode(
    Gactcoeff(runtime: adjFac.runtime),
    @[Gvalue(adjFac)],
    adjCoeffg,
    "adjCoeff")

proc mulscb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let x = Gscalar(z.inputs[0])
  let y = Gactcoeff(z.inputs[1])
  let upstream = requireUpstream(zb, "s*c backward", Gactcoeff)
  if i == 0:
    return redot(upstream, y)
  x * upstream

proc mulscf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let y = Gactcoeff(v.inputs[1])
  let z = Gactcoeff(v)
  z.cval = x.sval * y.cval

let mulsc = Gfunc(forward: mulscf, backward: mulscb, name: "s*c")

proc `*`*(x: Gscalar, y: Gactcoeff): Gactcoeff =
  graphNode(
    Gactcoeff(runtime: x.runtime),
    @[Gvalue(x), Gvalue(y)],
    mulsc,
    "s*c")

proc addccb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  requireUpstream(zb, "c+c backward", Gactcoeff)

proc addccf(v: Gvalue) =
  let x = Gactcoeff(v.inputs[0])
  let y = Gactcoeff(v.inputs[1])
  let z = Gactcoeff(v)
  z.cval = x.cval
  for a, b in fields(z.cval, y.cval):
    a += b

let addcc = Gfunc(forward: addccf, backward: addccb, name: "c+c")

proc `+`*(x, y: Gactcoeff): Gactcoeff =
  graphNode(
    Gactcoeff(runtime: x.runtime),
    @[Gvalue(x), Gvalue(y)],
    addcc,
    "c+c")

proc redotccb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  bilinearBackward(zb, z, i, Gactcoeff)

proc redotccf(v: Gvalue) =
  let x = Gactcoeff(v.inputs[0])
  let y = Gactcoeff(v.inputs[1])
  let z = Gscalar(v)
  var t = 0.0
  for a, b in fields(x.cval, y.cval):
    t += a * b
  z.sval = t

let redotcc = Gfunc(forward: redotccf, backward: redotccb, name: "redotcc")

proc redot*(x: Gactcoeff, y: Gactcoeff): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x), Gvalue(y)], redotcc, "redotcc")

proc actionCoeffLike(beta: Gscalar, coeffs: GaugeActionCoeffs): Gactcoeff =
  beta * toGvalue(beta.runtime, coeffs)

proc actWilson*(beta: Gscalar): Gactcoeff =
  beta.actionCoeffLike(GaugeActionCoeffs(plaq: 1.0))
proc actSymanzik*(beta: Gscalar): Gactcoeff =
  beta.actionCoeffLike(
    GaugeActionCoeffs(plaq: 1.0 - 8.0 * C1Symanzik, rect: C1Symanzik),
  )
proc actIwasaki*(beta: Gscalar): Gactcoeff =
  beta.actionCoeffLike(
    GaugeActionCoeffs(plaq: 1.0 - 8.0 * C1Iwasaki, rect: C1Iwasaki),
  )
proc actDBW2*(beta: Gscalar): Gactcoeff =
  beta.actionCoeffLike(
    GaugeActionCoeffs(plaq: 1.0 - 8.0 * C1DBW2, rect: C1DBW2),
  )
proc actAdj*(beta: Gscalar, adjFac: Gscalar): Gactcoeff = beta * adjCoeff(adjFac)

type Gcoeff = ref object of Gfunc
  basis: GaugeActionCoeffs
  family: proc(gc: GaugeActionCoeffs): bool {.nimcall.}

proc coeffb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  scaledUpstreamOr(zb, Gscalar, toGvalue(z.runtime, Gcoeff(z.gfunc).basis))

proc coefff(v: Gvalue) =
  # Guard the family on every evaluation: the reference graphs built through
  # this node cover only their kernel family, and the coefficient value can be
  # updated after the graph is built.
  let f = Gcoeff(v.gfunc)
  let gc = Gactcoeff(v.inputs[0]).cval
  if not f.family(gc):
    raiseUnsupportedPath(f.name, "coefficient set outside its family")
  Gscalar(v).sval =
    if f.basis.rect != 0: gc.rect
    elif f.basis.adjplaq != 0: gc.adjplaq
    else: gc.plaq

let coeffPlaqg = Gcoeff(forward: coefff, backward: coeffb, name: "coeffPlaq",
                        basis: GaugeActionCoeffs(plaq: 1.0), family: isPlaqRect)
let coeffPlaqOnlyg = Gcoeff(forward: coefff, backward: coeffb, name: "coeffPlaqOnly",
                            basis: GaugeActionCoeffs(plaq: 1.0), family: isPlaqOnly)
let coeffRectg = Gcoeff(forward: coefff, backward: coeffb, name: "coeffRect",
                        basis: GaugeActionCoeffs(rect: 1.0), family: isPlaqRect)
let coeffPlaqAg = Gcoeff(forward: coefff, backward: coeffb, name: "coeffPlaqA",
                         basis: GaugeActionCoeffs(plaq: 1.0), family: isAdjPlaq)
let coeffAdjg = Gcoeff(forward: coefff, backward: coeffb, name: "coeffAdj",
                       basis: GaugeActionCoeffs(adjplaq: 1.0), family: isAdjPlaq)

proc coeff(c: Gactcoeff, f: Gcoeff): Gscalar =
  ## One coefficient as a scalar graph value, with f's family guard.
  graphNode(scalarNodeLike(c), @[Gvalue(c)], f, f.name)

proc plaqSum(g: Ggauge): Gscalar =
  ## sum_{mu>nu} sum_x retr P_munu(x) over hop chains.
  let nd = g.gval.len
  let unit = g.unitFieldLike
  result = toGvalue(g.runtime, 0.0)
  for mu in 1..<nd:
    for nu in 0..<mu:
      result = result + retr(transport(g, unit, plaqPath(mu, nu)))

proc plaqActionGraph(c: Gactcoeff, g: Ggauge): Gscalar =
  ## Plaquette-only reference, S = -(c_plaq/nc) plaqSum; the replica behind
  ## gaugeActionDeriv2. Its coefficient node rejects other families.
  const nc = g.gval[0][0].nrows
  (toGvalue(g.runtime, -1.0/float(nc)) * coeff(c, coeffPlaqOnlyg)) * plaqSum(g)

proc gaugeActionGraph*(c: Gactcoeff, g: Ggauge): Gscalar =
  ## Reference action over the plaquette and rectangle families,
  ##   S = -(1/nc) [c_plaq sum_{mu>nu} sum_x retr P_munu
  ##                + c_rect sum_{mu>nu} sum_x (retr R_{2x1} + retr R_{1x2})],
  ## differentiable to any order in the field and in the coefficients. It
  ## matches `gaugeAction` for those families; the coefficient nodes reject
  ## the others at evaluation. The rectangles are shared path products
  ## (lineProducts) and are evaluated only when their coefficient is nonzero.
  discard sharedGraphRuntime(
    [Gvalue(c), Gvalue(g)], "gaugeActionGraph")
  let nd = g.gval.len
  const nc = g.gval[0][0].nrows
  var paths: seq[seq[int]]
  for mu in 1..<nd:
    for nu in 0..<mu:
      let a = mu + 1
      let b = nu + 1
      paths.add @[a, a, b, -a, -a, -b]
      paths.add @[a, b, b, -a, -b, -b]
  let zero = toGvalue(g.runtime, 0.0)
  let cr = coeff(c, coeffRectg)
  var sr = zero
  for p in lineProducts(g, paths, origin = false):
    sr = sr + retr(p)
  toGvalue(g.runtime, -1.0/float(nc)) *
    (coeff(c, coeffPlaqg) * plaqSum(g) + cond(equal(cr, zero), zero, cr * sr))

proc adjPlaqAction*(c: Gactcoeff, g: Ggauge): Gscalar =
  ## Reference for the adjoint-plaquette kernel family (actionA),
  ##   S = c_plaq (a0 - sum retr P / nc) + c_adj (a0 - sum |tr P|^2 / nc^2),
  ##   a0 = nd (nd-1)/2 * volume,
  ## differentiable to any order in the field and the coefficients. Matches
  ## `gaugeAction` when adjplaq != 0 (the kernel switches family on it); the
  ## coefficient nodes reject other families at evaluation.
  discard sharedGraphRuntime(
    [Gvalue(c), Gvalue(g)], "adjPlaqAction")
  let nd = g.gval.len
  const nc = g.gval[0][0].nrows
  let unit = g.unitFieldLike
  var sp = toGvalue(g.runtime, 0.0)
  var sa = toGvalue(g.runtime, 0.0)
  for mu in 1..<nd:
    for nu in 0..<mu:
      let p = transport(g, unit, plaqPath(mu, nu))
      sp = sp + retr(p)
      sa = sa + norm2(trace(p))
  let a0 = toGvalue(g.runtime, 0.5 * float(nd * (nd - 1) * g.gval[0].l.physVol))
  coeff(c, coeffPlaqAg) * (a0 - toGvalue(g.runtime, 1.0/float(nc)) * sp) +
    coeff(c, coeffAdjg) * (a0 - toGvalue(g.runtime, 1.0/float(nc*nc)) * sa)

proc gaugeActionDerivGraph(c: Gactcoeff, g: Ggauge): Ggauge =
  ## Grad-complete replica of gaugeActionDeriv for plaquette-only c.
  Ggauge(gradSeeded(plaqActionGraph(c, g), g, toGvalue(g.runtime, 1.0)))

# --- gauge-action graph operations ---

proc gaugeActionDeriv*(c: Gactcoeff, g: Ggauge): Ggauge
proc gaugeActionDeriv2*(b: Ggauge, c: Gactcoeff, g: Ggauge): Ggauge
proc gaugeActionDeriv2Subset(b: Ggauge, c: Gactcoeff, g: Ggauge, parity, dir: int): Ggauge

proc gaugeForce*(c: Gactcoeff, g: Ggauge): Ggauge

proc gaugeActionb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let c = Gactcoeff(z.inputs[0])
  let g = Ggauge(z.inputs[1])
  if i == 0:
    # This layer does not differentiate learned coefficients.
    raiseCoeffBackwardUnsupported("gaugeAction backward")
  scaledUpstreamOr(
    zb,
    Gscalar,
    gaugeActionDeriv(c, g))

proc gaugeActionf(v: Gvalue) =
  let c = Gactcoeff(v.inputs[0])
  let g = Ggauge(v.inputs[1])
  let gc = c.cval
  let z = Gscalar(v)
  z.sval = evalGaugeActionValue(gc, g.gval)

let gaugeActiong = Gfunc(
  forward: gaugeActionf,
  backward: gaugeActionb,
  name: "gaugeAction")

proc gaugeAction*(c: Gactcoeff, g: Ggauge): Gscalar =
  graphNode(scalarNodeLike(c), @[Gvalue(c), Gvalue(g)], gaugeActiong, "gaugeAction")

proc gaugeActionDerivb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let c = Gactcoeff(z.inputs[0])
  let g = Ggauge(z.inputs[1])
  if i == 0:
    # This layer does not differentiate learned coefficients.
    raiseCoeffBackwardUnsupported("gaugeActionDeriv backward")
  gaugeActionDeriv2(
    requireUpstream(zb, "gaugeActionDeriv backward", Ggauge),
    c,
    g)

proc gaugeActionDerivf(v: Gvalue) =
  let c = Gactcoeff(v.inputs[0])
  let g = Ggauge(v.inputs[1])
  let gc = c.cval
  let z = Ggauge(v)
  evalGaugeForceValue(gc, g.gval, z.gval)

let gaugeActionDerivg = Gfunc(
  forward: gaugeActionDerivf,
  backward: gaugeActionDerivb,
  name: "gaugeActionDeriv")

proc gaugeActionDeriv*(c: Gactcoeff, g: Ggauge): Ggauge =
  graphNode(g.gaugeNodeLike, @[Gvalue(c), Gvalue(g)], gaugeActionDerivg, "gaugeActionDeriv")

proc gaugeForceb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let
    c = Gactcoeff(z.inputs[0])
    g = Ggauge(z.inputs[1])
  if i == 0:
    raiseCoeffBackwardUnsupported("gaugeForce backward")
  let proj = projTAH(requireUpstream(zb, "gaugeForce backward", Ggauge))
  gaugeActionDeriv2(proj * g, c, g) + proj.adjmul(gaugeActionDeriv(c, g))

proc gaugeForcef(v: Gvalue) =
  let
    c = Gactcoeff(v.inputs[0])
    g = Ggauge(v.inputs[1])
    z = Ggauge(v)
  evalProjectedGaugeForceValue(c.cval, g.gval, z.gval)

let gaugeForceg = Gfunc(
  forward: gaugeForcef,
  backward: gaugeForceb,
  name: "gaugeForce")

proc gaugeForce*(c: Gactcoeff, g: Ggauge): Ggauge =
  ## Project the action derivative directly into the force output.
  graphNode(g.gaugeNodeLike, @[Gvalue(c), Gvalue(g)], gaugeForceg, "gaugeForce")

type GsubsetDeriv = ref object of Ggauge
  sd: Shifter[DLatticeColorMatrixV, DColorMatrixV]
  sf, sb: seq[ShiftB[DColorMatrixV]]
  parity, dir: int

proc subsetDerivNodeLike(x: Ggauge, parity, dir: int): GsubsetDeriv =
  let
    g = x.gval.newOneOf
    ps = if parity == 0: "even" else: "odd"
  g.zeroGaugeStorage
  GsubsetDeriv(
    runtime: x.runtime,
    gval: g,
    sd: newShifter(g[0], dir, 1),
    sf: createShiftBufs(g[0], 1, ps),
    sb: createShiftBufs(g[0], -1, ps),
    parity: parity,
    dir: dir).assignStableNodeId

method newOneOf(x: GsubsetDeriv): Gvalue =
  x.subsetDerivNodeLike(x.parity, x.dir)

proc gaugeActionDeriv*(c: Gactcoeff, g: Ggauge, parity, dir: int): Ggauge =
  ## Subset Wilson derivative; its pullback scatters to all staple neighbours.
  requireParityDir(parity, dir, g.gval.len, "gaugeActionDerivSubset")
  proc fwd(v: Gvalue) =
    let c = Gactcoeff(v.inputs[0])
    let g = Ggauge(v.inputs[1])
    let z = GsubsetDeriv(v)
    evalGaugeForceSubset(c.cval, g.gval, z.gval, z.sd, z.sf, z.sb, parity, dir)
  proc bwd(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
    let c = Gactcoeff(z.inputs[0])
    let g = Ggauge(z.inputs[1])
    if i == 0:
      raiseCoeffBackwardUnsupported("gaugeActionDerivSubset backward")
    gaugeActionDeriv2Subset(requireUpstream(zb, "gaugeActionDerivSubset backward", Ggauge), c, g, parity, dir)
  graphNode(g.subsetDerivNodeLike(parity, dir), @[Gvalue(c), Gvalue(g)], Gfunc(forward: fwd, backward: bwd, name: "gaugeActionDerivSubset"), "gaugeActionDerivSubset")

proc gaugeActionDeriv2b(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  ## z = gaugeActionDeriv2(b, c, g) is the action Hessian at g applied to b.
  if i == 1:
    raiseCoeffBackwardUnsupported("gaugeActionDeriv2 backward")
  let c = Gactcoeff(z.inputs[1])
  let g = Ggauge(z.inputs[2])
  let u = requireUpstream(zb, "gaugeActionDeriv2 backward", Ggauge)
  if i == 0:
    # The Hessian is self-adjoint, so the b cotangent is H[u].
    return gaugeActionDeriv2(u, c, g)
  # Third derivative: differentiate a grad-complete replica of z, the seeded
  # pullback of the reference action gradient.
  # secondPullback keeps this the partial g contribution even when b is g;
  # see its doc.
  # coeffPlaq guards the coefficient family at evaluation; graph building
  # must not read coefficient values.
  let b = Ggauge(z.inputs[0])
  secondPullback(g, b, u, proc(slot: Ggauge): Gvalue =
    gaugeActionDerivGraph(c, slot))

proc gaugeActionDeriv2f(v: Gvalue) =
  let b = Ggauge(v.inputs[0])
  let c = Gactcoeff(v.inputs[1])
  let g = Ggauge(v.inputs[2])
  let gc = c.cval
  let z = Ggauge(v)
  evalGaugeForceJacobian(b.gval, gc, g.gval, z.gval)

let gaugeActionDeriv2g = Gfunc(
  forward: gaugeActionDeriv2f,
  backward: gaugeActionDeriv2b,
  name: "gaugeActionDeriv2")

proc gaugeActionDeriv2*(b: Ggauge, c: Gactcoeff, g: Ggauge): Ggauge =
  graphNode(g.gaugeNodeLike, @[Gvalue(b), c, g], gaugeActionDeriv2g, "gaugeActionDeriv2")

type GsubsetHess = ref object of Ggauge
  hdir: DLatticeColorMatrixV

method newOneOf(x: GsubsetHess): Gvalue =
  let g = x.gval.newOneOf
  g.zeroGaugeStorage
  GsubsetHess(
    runtime: x.runtime,
    gval: g,
    hdir: x.hdir.newOneOf).assignStableNodeId

proc gaugeActionDeriv2Subset(b: Ggauge, c: Gactcoeff, g: Ggauge, parity, dir: int): Ggauge =
  let terms = b.addTerms
  let nterms = terms.len
  proc bwd(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
    ## z = H_g[mask(sum_j b_j)] where mask keeps only (parity, dir).
    if i == nterms:
      raiseCoeffBackwardUnsupported("gaugeActionDeriv2Subset backward")
    let
      c = Gactcoeff(z.inputs[nterms])
      g = Ggauge(z.inputs[nterms + 1])
      u = requireUpstream(zb, "gaugeActionDeriv2Subset backward", Ggauge)
    if i < nterms:
      # Linear in each b term with self-adjoint H: cotangent is mask(H[u]).
      return maskSubset(parity, dir, gaugeActionDeriv2(u, c, g))
    # The g cotangent differentiates a grad-complete replica of z.
    # secondPullback keeps this the partial g contribution even when b is g;
    # see its doc.
    # coeffPlaq guards the coefficient family at evaluation.
    var bsum = Ggauge(z.inputs[0])
    for j in 1..<nterms:
      bsum = bsum + Ggauge(z.inputs[j])
    secondPullback(g, maskSubset(parity, dir, bsum), u,
      proc(slot: Ggauge): Gvalue = gaugeActionDerivGraph(c, slot))
  if nterms == 1:
    proc fwd(v: Gvalue) =
      let
        b = Ggauge(v.inputs[0])
        c = Gactcoeff(v.inputs[1])
        g = Ggauge(v.inputs[2])
        z = Ggauge(v)
      evalGaugeForceJacobianSubset(b.gval, c.cval, g.gval, z.gval, parity, dir)
    return graphNode(g.gaugeNodeLike, @[Gvalue(terms[0]), Gvalue(c), Gvalue(g)], Gfunc(forward: fwd, backward: bwd, name: "gaugeActionDeriv2Subset"), "gaugeActionDeriv2Subset")
  var inputs = newSeq[Gvalue](nterms + 2)
  for i, term in terms:
    inputs[i] = Gvalue(term)
  inputs[nterms] = Gvalue(c)
  inputs[nterms + 1] = Gvalue(g)
  proc fwd(v: Gvalue) =
    let c = Gactcoeff(v.inputs[nterms])
    let g = Ggauge(v.inputs[nterms + 1])
    let z = GsubsetHess(v)
    var h = newSeq[DLatticeColorMatrixV](nterms)
    for i in 0..<nterms:
      h[i] = Ggauge(v.inputs[i]).gval[dir]
    evalGaugeForceJacobianSubsetSum(h, z.hdir, c.cval, g.gval, z.gval, parity, dir)
  let z = g.gaugeNodeLike
  graphNode(
    GsubsetHess(
      runtime: z.runtime,
      gval: z.gval,
      hdir: g.gval[dir].newOneOf),
    inputs,
    Gfunc(forward: fwd, backward: bwd, name: "gaugeActionDeriv2Subset"),
    "gaugeActionDeriv2Subset")
