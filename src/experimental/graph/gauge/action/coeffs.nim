import ../../[core, scalar]
import ../../support/op
import layout, gauge, physics/qcdTypes
import ../shared, domain

type
  Gactcoeff* {.final.} = ref object of Gvalue
    cval: GaugeActionCoeffs

proc coeffNodeIn(grt: GraphRuntime): Gactcoeff =
  Gactcoeff(runtime: grt)

proc coeffNodeLike(anchor: Gvalue): Gactcoeff =
  coeffNodeIn(anchor.runtime)

proc getactcoeff*(x: Gvalue): GaugeActionCoeffs = Gactcoeff(x).cval

proc `getactcoeff=`*(x: Gvalue, c: GaugeActionCoeffs) =
  let gc = Gactcoeff(x)
  gc.cval = c

proc update*(x: Gvalue, c: GaugeActionCoeffs) =
  x.getactcoeff = c
  x.updated

proc toGvalue*(grt: GraphRuntime,
               x: GaugeActionCoeffs): Gactcoeff =
  result = Gactcoeff(
    cval: x,
    runtime: grt)
  result.updated

method newOneOf*(x: Gactcoeff): Gvalue =
  result = Gactcoeff(runtime: x.runtime)
method valCopy*(z: Gactcoeff, x: Gactcoeff) = z.cval = x.cval
method copyCompatible*(prototype: Gactcoeff, value: Gactcoeff): bool =
  prototype != nil and value != nil
method `$`*(x: Gactcoeff): string = $x.cval

proc initCoeffLeaf(grt: GraphRuntime,
                   coeffs: GaugeActionCoeffs): Gactcoeff =
  result = Gactcoeff(
    cval: coeffs,
    runtime: grt)

proc raiseCoeffBackwardUnsupported*(label: string) {.noreturn.} =
  raiseUnsupportedPath(label, "derivative with respect to gauge-action coefficients")

proc adjCoeffb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z.requireUnaryNodeView("adjCoeff backward")
  discard dep
  unaryBackwardCase("adjCoeff backward", i,
    block:
      return redot(requireUpstream(zb, "adjCoeff backward"),
                   toGvalue(z.runtime, GaugeActionCoeffs(adjplaq: 1.0))))

defineUnaryForward(adjCoefff, Gscalar, Gactcoeff, "adjCoeff forward"):
  z.cval = GaugeActionCoeffs(plaq: 1.0, adjplaq: x.getfloat)

let adjCoeffg = newGfunc(forward = adjCoefff, backward = adjCoeffb, name = "adjCoeff")

proc adjCoeff(adjFac: Gscalar): Gvalue =
  graphNode(coeffNodeLike(adjFac), @[Gvalue(adjFac)], adjCoeffg, "adjCoeff")

proc mulscb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView("s*c backward")
  discard dep
  binaryBackwardCase("s*c backward", i,
    block:
      return redot(requireUpstream(zb, "s*c backward"), view.y),
    block:
      return view.x * requireUpstream(zb, "s*c backward"))

defineBinaryForward(mulscf, Gscalar, Gactcoeff, Gactcoeff, "s*c forward"):
  z.cval = x.getfloat * y.cval

defineBinaryGraphOp(mulsc, `*`, Gscalar, Gactcoeff, x, y, coeffNodeLike(x), mulscf, mulscb, "s*c")

proc redotccb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView("redotcc backward")
  discard dep
  binaryBackwardCase("redotcc backward", i,
    block:
      return scaledUpstreamOr(zb, view.y),
    block:
      return scaledUpstreamOr(zb, view.x))

defineBinaryForward(redotccf, Gactcoeff, Gactcoeff, Gscalar, "redotcc forward"):
  var t = 0.0
  for a, b in fields(x.cval, y.cval):
    t += a * b
  z.getfloat = t

defineBinaryGraphOp(redotcc, redot, Gactcoeff, Gactcoeff, x, y, scalarNodeLike(x), redotccf, redotccb, "redotcc")

method actWilson*(beta: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("actWilson(" & $beta & ")")
method actSymanzik*(beta: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("actSymanzik(" & $beta & ")")
method actIwasaki*(beta: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("actIwasaki(" & $beta & ")")
method actDBW2*(beta: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("actDBW2(" & $beta & ")")
method actAdj*(beta: Gvalue, adjFac: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("actAdj(" & $beta & "," & $adjFac & ")")

method actWilson*(beta: Gscalar): Gvalue =
  beta * initCoeffLeaf(beta.runtime, GaugeActionCoeffs(plaq: 1.0))
method actSymanzik*(beta: Gscalar): Gvalue =
  beta * initCoeffLeaf(beta.runtime,
    GaugeActionCoeffs(plaq: 1.0 - 8.0 * C1Symanzik, rect: C1Symanzik),
  )
method actIwasaki*(beta: Gscalar): Gvalue =
  beta * initCoeffLeaf(beta.runtime,
    GaugeActionCoeffs(plaq: 1.0 - 8.0 * C1Iwasaki, rect: C1Iwasaki),
  )
method actDBW2*(beta: Gscalar): Gvalue =
  beta * initCoeffLeaf(beta.runtime,
    GaugeActionCoeffs(plaq: 1.0 - 8.0 * C1DBW2, rect: C1DBW2),
  )
method actAdj*(beta: Gscalar, adjFac: Gscalar): Gvalue = beta * adjCoeff(adjFac)

proc actWilson*(grt: GraphRuntime,
                beta: float): Gvalue =
  actWilson(toGvalue(grt, beta))

proc actSymanzik*(grt: GraphRuntime,
                  beta: float): Gvalue =
  actSymanzik(toGvalue(grt, beta))

proc actIwasaki*(grt: GraphRuntime,
                 beta: float): Gvalue =
  actIwasaki(toGvalue(grt, beta))

proc actDBW2*(grt: GraphRuntime,
              beta: float): Gvalue =
  actDBW2(toGvalue(grt, beta))

proc actAdj*(grt: GraphRuntime,
             beta: float,
             adjFac: float): Gvalue =
  actAdj(toGvalue(grt, beta), toGvalue(grt, adjFac))
