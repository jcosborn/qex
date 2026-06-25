## Gauge-action coefficient value type and all gauge-action graph operations.
## Merged from the former coeffs.nim; the physics-evaluation bridge (raw QEX
## action routines) stays separate in domain.nim.
import ../../[core, scalar]
import ../../scalar/types
import ../../support/op
import layout, gauge, physics/qcdTypes
import ../shared, ../basic_ops, ../fused_ops, domain

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
proc redot*(x: Gactcoeff, y: Gactcoeff): Gscalar

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

# --- gauge-action graph operations ---

proc gaugeActionDeriv*(c: Gactcoeff, g: Ggauge): Ggauge
proc gaugeActionDeriv2*(b: Ggauge, c: Gactcoeff, g: Ggauge): Ggauge

proc gaugeForce*(c: Gactcoeff, g: Ggauge): Ggauge =
  ## Gauge force is the projected derivative consumed by the HMC integrators.
  contractProjTAH(gaugeActionDeriv(c, g), g)

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

proc gaugeActionDeriv2b(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  if i == 0:
    raiseUnsupportedPath(
      "gaugeActionDeriv2 backward",
      "derivative with respect to force-direction input")
  if i == 1:
    raiseCoeffBackwardUnsupported("gaugeActionDeriv2 backward")
  raiseUnsupportedPath(
    "gaugeActionDeriv2 backward",
    "higher derivatives with respect to the gauge field")

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
