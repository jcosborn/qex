## Gauge-action graph values and operations.
import ../../[core, scalar]
import ../../scalar/types
import ../../support/op
import layout, gauge, physics/qcdTypes
import ../types, ../basic_ops, ../fused_ops, ../field_ops, ../transport, ../cfield, domain
from ../stencil import stapleSum

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
method isZero*(x: Gactcoeff): bool =
  for f in x.cval.fields:
    if f != 0.0:
      return false
  true
method zeroLike*(x: Gactcoeff): Gvalue =
  result = x.newOneOf
  result.markStaticZeroLeaf
method valCopy*(z: Gactcoeff, x: Gvalue) =
  z.cval = Gactcoeff(x).cval
method copyCompatible*(prototype: Gactcoeff, value: Gvalue): bool =
  value of Gactcoeff
method `$`*(x: Gactcoeff): string = $x.cval

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

let adjCoeffg = Gfunc(bufferMode: bmFull, forward: adjCoefff, backward: adjCoeffb, name: "adjCoeff")

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

let mulsc = Gfunc(bufferMode: bmFull, forward: mulscf, backward: mulscb, name: "s*c")

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

let addcc = Gfunc(bufferMode: bmFull, forward: addccf, backward: addccb, name: "c+c")

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

let redotcc = Gfunc(bufferMode: bmFull, forward: redotccf, backward: redotccb, name: "redotcc")

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

proc coeffBasis(c: Gactcoeff, f: Gcoeff, basis: GaugeActionCoeffs): Gactcoeff

proc coeffBasisf(v: Gvalue) =
  let f = Gcoeff(v.gfunc)
  if not f.family(Gactcoeff(v.inputs[0]).cval):
    raiseUnsupportedPath(f.name, "coefficient set outside its family")
  Gactcoeff(v).cval = f.basis

proc coeffBasisb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  # Constant on the supported family; keep its domain visible at later orders.
  coeffBasis(Gactcoeff(z.inputs[0]), Gcoeff(z.gfunc), GaugeActionCoeffs())

proc coeffBasis(c: Gactcoeff, f: Gcoeff, basis: GaugeActionCoeffs): Gactcoeff =
  graphNode(Gactcoeff(runtime: c.runtime), @[Gvalue(c)],
    Gcoeff(bufferMode: bmFull, forward: coeffBasisf, backward: coeffBasisb, name: "coeffBasis",
           basis: basis, family: f.family), "coeffBasis")

proc coeffb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let f = Gcoeff(z.gfunc)
  scaledUpstreamOr(zb, Gscalar, coeffBasis(Gactcoeff(z.inputs[0]), f, f.basis))

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
    elif f.basis.pgm != 0: gc.pgm
    elif f.basis.adjplaq != 0: gc.adjplaq
    else: gc.plaq

let coeffPlaqg = Gcoeff(bufferMode: bmFull, forward: coefff, backward: coeffb, name: "coeffPlaq",
                        basis: GaugeActionCoeffs(plaq: 1.0), family: isFundamental)
let coeffRectg = Gcoeff(bufferMode: bmFull, forward: coefff, backward: coeffb, name: "coeffRect",
                        basis: GaugeActionCoeffs(rect: 1.0), family: isFundamental)
let coeffPgmg = Gcoeff(bufferMode: bmFull, forward: coefff, backward: coeffb, name: "coeffPgm",
                       basis: GaugeActionCoeffs(pgm: 1.0), family: isFundamental)
# Return adjplaq: zero selects the fundamental family, nonzero the adjoint family.
let coeffFamilyg = Gcoeff(bufferMode: bmFull, forward: coefff, backward: coeffb, name: "coeffFamily",
                          basis: GaugeActionCoeffs(adjplaq: 1.0), family: isGaugeAction)
let coeffPlaqAg = Gcoeff(bufferMode: bmFull, forward: coefff, backward: coeffb, name: "coeffPlaqA",
                         basis: GaugeActionCoeffs(plaq: 1.0), family: isAdjPlaq)
let coeffAdjg = Gcoeff(bufferMode: bmFull, forward: coefff, backward: coeffb, name: "coeffAdj",
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

proc loopSum(g: Ggauge, kind: int): Gscalar =
  ## Basis 0 = plaquette, 1 = rectangle, 2 = parallelogram.
  if kind == 0:
    return plaqSum(g)
  var paths: seq[seq[int]]
  for mu in 1..<g.gval.len:
    for nu in 0..<mu:
      let a = mu + 1
      let b = nu + 1
      if kind == 1:
        paths.add @[a, a, b, -a, -a, -b]
        paths.add @[a, b, b, -a, -b, -b]
      else:
        for sg in 0..<nu:
          let d = sg + 1
          # gaugeAction2's ts1, ts2, ts3, and ts7.
          paths.add @[a, b, d, -a, -b, -d]
          paths.add @[a, d, b, -a, -d, -b]
          paths.add @[b, a, d, -b, -a, -d]
          paths.add @[a, -b, d, -a, b, -d]
  result = toGvalue(g.runtime, 0.0)
  if paths.len == 0:
    return
  for p in lineProducts(g, paths, origin = false):
    result = result + retr(p)

proc gaugeActionGraph*(c: Gactcoeff, g: Ggauge): Gscalar =
  ## S = -(c_plaq sum Re tr P + c_rect sum Re tr R + c_pgm sum Re tr C)/Nc.
  ## The fundamental family is adjplaq == 0. scalarScale skips an inactive
  ## sum while preserving its coefficient derivative, including at zero.
  discard sharedGraphRuntime([Gvalue(c), Gvalue(g)], "gaugeActionGraph")
  const nc = g.gval[0][0].nrows
  toGvalue(g.runtime, -1.0/float(nc)) *
    (Gscalar(loopSum(g, 0).scaleLike(coeff(c, coeffPlaqg))) +
     Gscalar(loopSum(g, 1).scaleLike(coeff(c, coeffRectg))) +
     Gscalar(loopSum(g, 2).scaleLike(coeff(c, coeffPgmg))))

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
  ## Fundamental-family replica behind the numerical Hessian pullbacks.
  Ggauge(gradSeeded(gaugeActionGraph(c, g), g, toGvalue(g.runtime, 1.0)))

proc actionCoeffPullback(c: Gactcoeff, g: Ggauge, u: Gscalar): Gactcoeff =
  let s = slotVar(c)
  # Differentiate each family before selecting it, so no inactive family's
  # guarded coefficient basis enters a partial derivative of the other one.
  cond(equal(coeff(c, coeffFamilyg), toGvalue(c.runtime, 0.0)),
       Gactcoeff(gradSeeded(gaugeActionGraph(s, g), s, u)),
       Gactcoeff(gradSeeded(adjPlaqAction(s, g), s, u)))

proc derivCoeffPullback(c: Gactcoeff, g, u: Ggauge,
                        fundamental = false): Gactcoeff =
  let s = slotVar(c)
  let x = slotVar(g)
  let df = gaugeActionDerivGraph(s, x)
  let cf = Gactcoeff(gradSeeded(df, s, u))
  if fundamental:
    return cf
  let da = gradSeeded(adjPlaqAction(s, x), x, toGvalue(g.runtime, 1.0))
  cond(equal(coeff(c, coeffFamilyg), toGvalue(c.runtime, 0.0)),
       cf, Gactcoeff(gradSeeded(da, s, u)))

proc hessCoeffPullback(c: Gactcoeff, g, b, u: Ggauge): Gactcoeff =
  let s = slotVar(c)
  let x = slotVar(g)
  let d = gaugeActionDerivGraph(s, x)
  let h = gradSeeded(d, x, b)
  Gactcoeff(gradSeeded(h, s, u))

type GactionJet = ref object of Gfunc
  nseeds: int

proc actionJet(c: Gactcoeff, g: Ggauge, seeds: openArray[Ggauge]): Ggauge

proc actionJetInputs(v: Gvalue, mode: InputWalkMode, visit: GnodeVisit) =
  let off = GactionJet(v.gfunc).nseeds + 2
  for i in 0..<off:
    visit v.inputs[i]
  if mode == iwmBackward:
    return
  let c = Gactcoeff(v.inputs[0]).cval
  for k, a in [c.plaq, c.rect, c.pgm]:
    if mode != iwmEval or a != 0.0:
      visit v.inputs[off + k]

proc actionJetf(v: Gvalue) =
  let c = Gactcoeff(v.inputs[0]).cval
  if not c.isFundamental:
    raiseUnsupportedPath("actionJet", "coefficient set outside its fundamental family")
  let off = GactionJet(v.gfunc).nseeds + 2
  let z = Ggauge(v)
  threads:
    for mu in 0..<z.gval.len:
      z.gval[mu] := 0
      if c.plaq != 0:
        z.gval[mu] += c.plaq * Ggauge(v.inputs[off]).gval[mu]
      if c.rect != 0:
        z.gval[mu] += c.rect * Ggauge(v.inputs[off + 1]).gval[mu]
      if c.pgm != 0:
        z.gval[mu] += c.pgm * Ggauge(v.inputs[off + 2]).gval[mu]

proc actionJetb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let c = Gactcoeff(z.inputs[0])
  let g = Ggauge(z.inputs[1])
  let n = GactionJet(z.gfunc).nseeds
  let u = requireUpstream(zb, "actionJet backward", Ggauge)
  if i == 0:
    let off = n + 2
    return redot(u, Ggauge(z.inputs[off])) * coeffBasis(c, coeffPlaqg, GaugeActionCoeffs(plaq: 1.0)) +
           redot(u, Ggauge(z.inputs[off + 1])) * coeffBasis(c, coeffRectg, GaugeActionCoeffs(rect: 1.0)) +
           redot(u, Ggauge(z.inputs[off + 2])) * coeffBasis(c, coeffPgmg, GaugeActionCoeffs(pgm: 1.0))
  var seeds: seq[Ggauge]
  for j in 0..<n:
    seeds.add Ggauge(z.inputs[j + 2])
  if i == 1:
    seeds.add u
  else:
    seeds[i - 2] = u
  actionJet(c, g, seeds)

proc actionJet(c: Gactcoeff, g: Ggauge, seeds: openArray[Ggauge]): Ggauge =
  ## D^(seeds.len+1) S[g] contracted with every seed, leaving one gauge slot.
  ## The hidden basis replicas are eval inputs only. Public backward inputs
  ## remain c, g, and the live seeds, including aliases and dependent seeds.
  ## Each derivative order skips evaluation of zero-weight loop families.
  const nc = g.gval[0][0].nrows
  var inputs = @[Gvalue(c), Gvalue(g)]
  for seed in seeds:
    inputs.add Gvalue(seed)
  discard sharedGraphRuntime(inputs, "actionJet")
  for kind in 0..2:
    if kind == 0:
      inputs.add toGvalue(g.runtime, -1.0/float(nc)) * stapleSum(g, seeds)
      continue
    let slot = slotVar(g)
    let a = toGvalue(g.runtime, -1.0/float(nc)) * loopSum(slot, kind)
    var d = gradSeeded(a, slot, toGvalue(g.runtime, 1.0))
    for seed in seeds:
      d = gradSeeded(d, slot, seed)
    inputs.add d
  graphNode(g.gaugeNodeLike, inputs,
    GactionJet(bufferMode: bmFull, forward: actionJetf, backward: actionJetb, inputView: actionJetInputs,
               name: "actionJet", nseeds: seeds.len), "actionJet")

# --- gauge-action graph operations ---

proc gaugeActionDeriv*(c: Gactcoeff, g: Ggauge): Ggauge
proc gaugeActionDeriv2*(b: Ggauge, c: Gactcoeff, g: Ggauge): Ggauge
proc gaugeActionDeriv2Subset(b: Ggauge, c: Gactcoeff, g: Ggauge, parity, dir: int): Ggauge

proc gaugeForce*(c: Gactcoeff, g: Ggauge): Ggauge

proc gaugeActionb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let c = Gactcoeff(z.inputs[0])
  let g = Ggauge(z.inputs[1])
  if i == 0:
    return actionCoeffPullback(c, g, Gscalar(rootedUpstream(zb, z)))
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

let gaugeActiong = Gfunc(bufferMode: bmFull,
  forward: gaugeActionf,
  backward: gaugeActionb,
  name: "gaugeAction")

proc gaugeAction*(c: Gactcoeff, g: Ggauge): Gscalar =
  ## Coefficient pullbacks vary plaq/rect/pgm when adjplaq == 0, or plaq/adjplaq
  ## in ActionA (rect == pgm == 0, adjplaq != 0); other slots are zero.
  ## The family switch itself is not a differentiable coefficient direction.
  ## For actAdj(beta, adjFac) at beta == 0, this gives dS/dbeta = -sum Re tr P/Nc
  ## even with nonzero adjFac, since the zero coefficient value is fundamental.
  graphNode(scalarNodeLike(c), @[Gvalue(c), Gvalue(g)], gaugeActiong, "gaugeAction")

proc gaugeActionDerivb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let c = Gactcoeff(z.inputs[0])
  let g = Ggauge(z.inputs[1])
  if i == 0:
    return derivCoeffPullback(c, g, requireUpstream(zb, "gaugeActionDeriv backward", Ggauge))
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

let gaugeActionDerivg = Gfunc(bufferMode: bmFull,
  forward: gaugeActionDerivf,
  backward: gaugeActionDerivb,
  name: "gaugeActionDeriv")

proc gaugeActionDeriv*(c: Gactcoeff, g: Ggauge): Ggauge =
  graphNode(g.gaugeNodeLike, @[Gvalue(c), Gvalue(g)], gaugeActionDerivg, "gaugeActionDeriv")

proc gaugeForceb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let
    c = Gactcoeff(z.inputs[0])
    g = Ggauge(z.inputs[1])
  let proj = projTAH(requireUpstream(zb, "gaugeForce backward", Ggauge))
  if i == 0:
    return derivCoeffPullback(c, g, proj * g)
  gaugeActionDeriv2(proj * g, c, g) + proj.adjmul(gaugeActionDeriv(c, g))

proc gaugeForcef(v: Gvalue) =
  let
    c = Gactcoeff(v.inputs[0])
    g = Ggauge(v.inputs[1])
    z = Ggauge(v)
  evalProjectedGaugeForceValue(c.cval, g.gval, z.gval)

let gaugeForceg = Gfunc(bufferMode: bmFull,
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

method ensureStorage*(x: GsubsetDeriv) =
  procCall Ggauge(x).ensureStorage
  if x.sd.field.isNil:
    let ps = if x.parity == 0: "even" else: "odd"
    x.sd = newShifter(x.gval[0], x.dir, 1)
    x.sf = createShiftBufs(x.gval[0], 1, ps)
    x.sb = createShiftBufs(x.gval[0], -1, ps)

method releaseWork*(x: GsubsetDeriv) =
  x.sd = default(typeof(x.sd))
  x.sf = @[]
  x.sb = @[]

method releaseStorage*(x: GsubsetDeriv) =
  x.releaseWork
  procCall Ggauge(x).releaseStorage

proc subsetDerivNodeLike(x: Ggauge, parity, dir: int): GsubsetDeriv =
  result = GsubsetDeriv(runtime: x.runtime, gval: x.gaugeNodeLike.gval, parity: parity, dir: dir)
  result.assignStableNodeId

method newOneOf(x: GsubsetDeriv): Gvalue =
  x.subsetDerivNodeLike(x.parity, x.dir)

proc gaugeActionDeriv*(c: Gactcoeff, g: Ggauge, parity, dir: int): Ggauge =
  ## Masked fundamental derivative; its pullback reaches all affected links.
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
      return derivCoeffPullback(c, g, maskSubset(parity, dir,
        requireUpstream(zb, "gaugeActionDerivSubset backward", Ggauge)), true)
    gaugeActionDeriv2Subset(requireUpstream(zb, "gaugeActionDerivSubset backward", Ggauge), c, g, parity, dir)
  graphNode(g.subsetDerivNodeLike(parity, dir), @[Gvalue(c), Gvalue(g)], Gfunc(bufferMode: bmZero, forward: fwd, backward: bwd, name: "gaugeActionDerivSubset"), "gaugeActionDerivSubset")

proc gaugeActionDeriv2b(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  ## z = gaugeActionDeriv2(b, c, g) is the action Hessian at g applied to b.
  let c = Gactcoeff(z.inputs[1])
  let g = Ggauge(z.inputs[2])
  let u = requireUpstream(zb, "gaugeActionDeriv2 backward", Ggauge)
  if i == 1:
    return hessCoeffPullback(c, g, Ggauge(z.inputs[0]), u)
  if i == 0:
    # The Hessian is self-adjoint, so the b cotangent is H[u].
    return gaugeActionDeriv2(u, c, g)
  actionJet(c, g, [Ggauge(z.inputs[0]), u])

proc gaugeActionDeriv2f(v: Gvalue) =
  let b = Ggauge(v.inputs[0])
  let c = Gactcoeff(v.inputs[1])
  let g = Ggauge(v.inputs[2])
  let gc = c.cval
  let z = Ggauge(v)
  evalGaugeForceJacobian(b.gval, gc, g.gval, z.gval)

let gaugeActionDeriv2g = Gfunc(bufferMode: bmFull,
  forward: gaugeActionDeriv2f,
  backward: gaugeActionDeriv2b,
  name: "gaugeActionDeriv2")

proc gaugeActionDeriv2*(b: Ggauge, c: Gactcoeff, g: Ggauge): Ggauge =
  graphNode(g.gaugeNodeLike, @[Gvalue(b), c, g], gaugeActionDeriv2g, "gaugeActionDeriv2")

type GsubsetHess = ref object of Ggauge
  hdir: DLatticeColorMatrixV

method ensureStorage*(x: GsubsetHess) =
  procCall Ggauge(x).ensureStorage
  x.hdir.ensureFieldStorage

method releaseWork*(x: GsubsetHess) =
  x.hdir.releaseFieldStorage

method releaseStorage*(x: GsubsetHess) =
  x.releaseWork
  procCall Ggauge(x).releaseStorage

method newOneOf(x: GsubsetHess): Gvalue =
  GsubsetHess(runtime: x.runtime, gval: x.gaugeNodeLike.gval,
              hdir: x.hdir.newShape).assignStableNodeId

proc gaugeActionDeriv2Subset(b: Ggauge, c: Gactcoeff, g: Ggauge, parity, dir: int): Ggauge =
  let terms = b.addTerms
  let nterms = terms.len
  proc bwd(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
    ## z = H_g[mask(sum_j b_j)] where mask keeps only (parity, dir).
    let
      c = Gactcoeff(z.inputs[nterms])
      g = Ggauge(z.inputs[nterms + 1])
      u = requireUpstream(zb, "gaugeActionDeriv2Subset backward", Ggauge)
    if i < nterms:
      # Linear in each b term with self-adjoint H: cotangent is mask(H[u]).
      return maskSubset(parity, dir, gaugeActionDeriv2(u, c, g))
    var bsum = Ggauge(z.inputs[0])
    for j in 1..<nterms:
      bsum = bsum + Ggauge(z.inputs[j])
    if i == nterms:
      return hessCoeffPullback(c, g, maskSubset(parity, dir, bsum), u)
    actionJet(c, g, [maskSubset(parity, dir, bsum), u])
  if nterms == 1:
    proc fwd(v: Gvalue) =
      let
        b = Ggauge(v.inputs[0])
        c = Gactcoeff(v.inputs[1])
        g = Ggauge(v.inputs[2])
        z = Ggauge(v)
      evalGaugeForceJacobianSubset(b.gval, c.cval, g.gval, z.gval, parity, dir)
    return graphNode(g.gaugeNodeLike, @[Gvalue(terms[0]), Gvalue(c), Gvalue(g)], Gfunc(bufferMode: bmFull, forward: fwd, backward: bwd, name: "gaugeActionDeriv2Subset"), "gaugeActionDeriv2Subset")
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
      hdir: g.gval[dir].newShape),
    inputs,
    Gfunc(bufferMode: bmFull, forward: fwd, backward: bwd, name: "gaugeActionDeriv2Subset"),
    "gaugeActionDeriv2Subset")
