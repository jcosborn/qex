import ../../[core, scalar]
import ../../support/op
import ../shared, ../basic_ops, ../fused_ops, domain, coeffs

proc gaugeAction*(c: Gactcoeff, g: Ggauge): Gscalar
proc gaugeActionDeriv*(c: Gactcoeff, g: Ggauge): Ggauge
proc gaugeActionDeriv2*(b: Ggauge, c: Gactcoeff, g: Ggauge): Ggauge

proc gaugeForce*(c: Gactcoeff, g: Ggauge): Ggauge =
  ## Gauge force is the projected derivative consumed by the HMC integrators.
  contractProjTAH(gaugeActionDeriv(c, g), g)

proc gaugeActionb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView(Gactcoeff, Ggauge, "gaugeAction backward", "coefficients", "gauge")
  discard dep
  binaryBackwardCase("gaugeAction backward", i,
    block:
      # This layer does not differentiate learned coefficients.
      raiseCoeffBackwardUnsupported("gaugeAction backward"),
    block:
      return scaledUpstreamOr(zb, Gscalar, gaugeActionDeriv(view.x, view.y), "gaugeAction backward"))

proc gaugeActionf(v: Gvalue) =
  let view = v.requireBinaryNodeView(
    Gactcoeff,
    Ggauge,
    "gaugeAction forward",
    "coefficients",
    "gauge")
  let gc = view.x.getactcoeff
  let g = view.y
  let z = v.requireScalar("gaugeAction forward result")
  z.getfloat = evalGaugeActionValue(gc, g.gval)

defineBinaryGraphOp(gaugeActiong, gaugeAction, Gactcoeff, Ggauge, c, g, scalarNodeLike(c), gaugeActionf, gaugeActionb, "gaugeAction")

proc gaugeActionDerivb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView(Gactcoeff, Ggauge, "gaugeActionDeriv backward", "coefficients", "gauge")
  discard dep
  binaryBackwardCase("gaugeActionDeriv backward", i,
    block:
      # This layer does not differentiate learned coefficients.
      raiseCoeffBackwardUnsupported("gaugeActionDeriv backward"),
    block:
      return gaugeActionDeriv2(requireUpstream(zb, Ggauge, "gaugeActionDeriv backward"),
                               view.x,
                               view.y))

proc gaugeActionDerivf(v: Gvalue) =
  let view = v.requireBinaryNodeView(
    Gactcoeff,
    Ggauge,
    "gaugeActionDeriv forward",
    "coefficients",
    "gauge")
  let gc = view.x.getactcoeff
  let g = view.y
  let z = v.requireGauge("gaugeActionDeriv forward result")
  evalGaugeForceValue(gc, g.gval, z.gval)

defineBinaryGraphOp(gaugeActionDerivg, gaugeActionDeriv, Gactcoeff, Ggauge, c, g, g.gaugeNodeLike, gaugeActionDerivf, gaugeActionDerivb, "gaugeActionDeriv")

proc gaugeActionDeriv2b(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard zb
  discard z.requireTernaryNodeView(
    Ggauge,
    Gactcoeff,
    Ggauge,
    "gaugeActionDeriv2 backward",
    "direction",
    "coefficients",
    "gauge")
  discard dep
  ternaryBackwardCase("gaugeActionDeriv2 backward", i,
    block:
      raiseUnsupportedPath("gaugeActionDeriv2 backward", "derivative with respect to force-direction input"),
    block:
      raiseCoeffBackwardUnsupported("gaugeActionDeriv2 backward"),
    block:
      raiseUnsupportedPath("gaugeActionDeriv2 backward", "higher derivatives with respect to the gauge field"))

proc gaugeActionDeriv2f(v: Gvalue) =
  let view = v.requireTernaryNodeView(
    Ggauge,
    Gactcoeff,
    Ggauge,
    "gaugeActionDeriv2 forward",
    "direction",
    "coefficients",
    "gauge")
  let b = view.x
  let gc = view.y.getactcoeff
  let g = view.z
  let z = v.requireGauge("gaugeActionDeriv2 forward result")
  evalGaugeForceJacobian(b.gval, gc, g.gval, z.gval)

let gaugeActionDeriv2g = newGfunc(
  forward = gaugeActionDeriv2f,
  backward = gaugeActionDeriv2b,
  name = "gaugeActionDeriv2")

proc gaugeActionDeriv2*(b: Ggauge, c: Gactcoeff, g: Ggauge): Ggauge =
  graphNode(g.gaugeNodeLike, @[Gvalue(b), c, g], gaugeActionDeriv2g, "gaugeActionDeriv2")
