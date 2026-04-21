import ../../[core, scalar]
import ../../support/op
import ../shared, domain, coeffs

method gaugeAction*(c: Gvalue, g: Gvalue): Gvalue {.base.} =
  raiseErrorBaseMethod("gaugeAction(" & $c & "," & $g & ")")
method gaugeActionDeriv*(c: Gvalue, g: Gvalue): Gvalue {.base.} =
  raiseErrorBaseMethod("gaugeActionDeriv(" & $c & "," & $g & ")")
method gaugeActionDeriv2*(b: Gvalue, c: Gvalue, g: Gvalue): Gvalue {.base.} =
  raiseErrorBaseMethod("gaugeActionDeriv2(" & $b & "," & $c & "," & $g & ")")

proc gaugeForce*(c: Gvalue, g: Gvalue): Gvalue =
  ## Gauge force is the projected derivative consumed by the HMC integrators.
  contractProjTAH(gaugeActionDeriv(c, g), g)

proc gaugeActionb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView("gaugeAction backward", "coefficients", "gauge")
  discard dep
  binaryBackwardCase("gaugeAction backward", i,
    block:
      # This layer does not differentiate learned coefficients.
      raiseCoeffBackwardUnsupported("gaugeAction backward"),
    block:
      return scaledUpstreamOr(zb, gaugeActionDeriv(view.x, view.y)))

proc gaugeActionf(v: Gvalue) =
  let view = v.requireBinaryNodeView("gaugeAction forward", "coefficients", "gauge")
  let gc = view.x.getactcoeff
  let g = Ggauge(view.y)
  let z = Gscalar(v)
  z.getfloat = evalGaugeActionValue(gc, g.gval)

defineBinaryGraphOp(gaugeActiong, gaugeAction, Gactcoeff, Ggauge, c, g, scalarNodeLike(c), gaugeActionf, gaugeActionb, "gaugeAction")

proc gaugeActionDerivb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView("gaugeActionDeriv backward", "coefficients", "gauge")
  discard dep
  binaryBackwardCase("gaugeActionDeriv backward", i,
    block:
      # This layer does not differentiate learned coefficients.
      raiseCoeffBackwardUnsupported("gaugeActionDeriv backward"),
    block:
      return gaugeActionDeriv2(requireUpstream(zb, "gaugeActionDeriv backward"),
                               view.x,
                               view.y))

proc gaugeActionDerivf(v: Gvalue) =
  let view = v.requireBinaryNodeView("gaugeActionDeriv forward", "coefficients", "gauge")
  let gc = view.x.getactcoeff
  let g = Ggauge(view.y)
  let z = Ggauge(v)
  evalGaugeForceValue(gc, g.gval, z.gval)

defineBinaryGraphOp(gaugeActionDerivg, gaugeActionDeriv, Gactcoeff, Ggauge, c, g, g.newOneOf, gaugeActionDerivf, gaugeActionDerivb, "gaugeActionDeriv")

proc gaugeActionDeriv2b(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard zb
  discard z.requireTernaryNodeView(
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
    "gaugeActionDeriv2 forward",
    "direction",
    "coefficients",
    "gauge")
  let b = Ggauge(view.x)
  let gc = view.y.getactcoeff
  let g = Ggauge(view.z)
  let z = Ggauge(v)
  evalGaugeForceJacobian(b.gval, gc, g.gval, z.gval)

let gaugeActionDeriv2g = newGfunc(
  forward = gaugeActionDeriv2f,
  backward = gaugeActionDeriv2b,
  name = "gaugeActionDeriv2")

method gaugeActionDeriv2*(b: Ggauge, c: Gactcoeff, g: Ggauge): Gvalue =
  graphNode(g.newOneOf, @[Gvalue(b), c, g], gaugeActionDeriv2g, "gaugeActionDeriv2")
