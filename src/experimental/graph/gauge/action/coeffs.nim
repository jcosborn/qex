import ../../[core, scalar]
import ../../support/op
import layout, gauge, physics/qcdTypes
import ../shared, domain

type
  Gactcoeff* {.final.} = ref object of Gvalue
    cval: GaugeActionCoeffs

proc coeffNodeIn(grt: GraphRuntime): Gactcoeff =
  Gactcoeff().attachRuntime(grt)

proc coeffNodeLike(anchor: Gvalue): Gactcoeff =
  coeffNodeIn(anchor.runtime)

proc requireActCoeff*(value: Gvalue,
                      label: string): Gactcoeff =
  if value == nil:
    raiseValueError(label & " requires non-nil gauge-action coefficient value")
  if not (value of Gactcoeff):
    raiseValueError(
      label & " expects gauge-action coefficient value, got:\n" &
      value.nodeRepr)
  Gactcoeff(value)

proc getactcoeff*(x: Gactcoeff): GaugeActionCoeffs =
  x.cval

proc `getactcoeff=`*(x: Gactcoeff, c: GaugeActionCoeffs) =
  x.cval = c

proc update*(x: Gactcoeff, c: GaugeActionCoeffs) =
  x.getactcoeff = c
  x.updated

proc toGvalue*(grt: GraphRuntime,
               x: GaugeActionCoeffs): Gactcoeff =
  result = Gactcoeff(
    cval: x).attachRuntime(grt)
  result.updated

method newOneOf*(x: Gactcoeff): Gvalue =
  result = Gactcoeff().attachRuntime(x.runtime)
proc valCopy*(z: Gactcoeff, x: Gactcoeff) = z.cval = x.cval
method valCopy*(z: Gactcoeff, x: Gvalue) =
  z.valCopy(x.requireActCoeff("gauge-action coefficient copy"))
proc copyCompatible*(prototype: Gactcoeff, value: Gactcoeff): bool =
  prototype != nil and value != nil
method copyCompatible*(prototype: Gactcoeff, value: Gvalue): bool =
  prototype != nil and value != nil and value of Gactcoeff
method `$`*(x: Gactcoeff): string = $x.cval

proc initCoeffLeaf(grt: GraphRuntime,
                   coeffs: GaugeActionCoeffs): Gactcoeff =
  result = Gactcoeff(
    cval: coeffs).attachRuntime(grt)
  result.updated

proc raiseCoeffBackwardUnsupported*(label: string) {.noreturn.} =
  raiseUnsupportedPath(label, "derivative with respect to gauge-action coefficients")

proc `*`*(x: Gscalar, y: Gactcoeff): Gactcoeff
proc redot*(x: Gactcoeff, y: Gactcoeff): Gscalar

method scaleLike*(contribution: Gactcoeff, upstream: Gvalue): Gvalue =
  upstream.requireScalar("gauge-action coefficient scale upstream") * contribution

proc adjCoeffb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard z.requireUnaryNodeView(Gscalar, "adjCoeff backward")
  discard dep
  unaryBackwardCase("adjCoeff backward", i,
    block:
      return redot(requireUpstream(zb, Gactcoeff, "adjCoeff backward"),
                   toGvalue(z.runtime, GaugeActionCoeffs(adjplaq: 1.0))))

defineUnaryForward(adjCoefff, Gscalar, Gactcoeff, "adjCoeff forward"):
  z.cval = GaugeActionCoeffs(plaq: 1.0, adjplaq: x.getfloat)

let adjCoeffg = newGfunc(forward = adjCoefff, backward = adjCoeffb, name = "adjCoeff")

proc adjCoeff(adjFac: Gscalar): Gactcoeff =
  graphNode(coeffNodeLike(adjFac), @[Gvalue(adjFac)], adjCoeffg, "adjCoeff")

proc mulscb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView(Gscalar, Gactcoeff, "s*c backward")
  discard dep
  binaryBackwardCase("s*c backward", i,
    block:
      return redot(requireUpstream(zb, Gactcoeff, "s*c backward"), view.y),
    block:
      return view.x * requireUpstream(zb, Gactcoeff, "s*c backward"))

defineBinaryForward(mulscf, Gscalar, Gactcoeff, Gactcoeff, "s*c forward"):
  z.cval = x.getfloat * y.cval

defineBinaryGraphOp(mulsc, `*`, Gscalar, Gactcoeff, x, y, coeffNodeLike(x), mulscf, mulscb, "s*c")

proc redotccb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView(Gactcoeff, Gactcoeff, "redotcc backward")
  discard dep
  binaryBackwardCase("redotcc backward", i,
    block:
      return scaledUpstreamOr(zb, Gscalar, view.y, "redotcc backward"),
    block:
      return scaledUpstreamOr(zb, Gscalar, view.x, "redotcc backward"))

defineBinaryForward(redotccf, Gactcoeff, Gactcoeff, Gscalar, "redotcc forward"):
  var t = 0.0
  for a, b in fields(x.cval, y.cval):
    t += a * b
  z.getfloat = t

defineBinaryGraphOp(redotcc, redot, Gactcoeff, Gactcoeff, x, y, scalarNodeLike(x), redotccf, redotccb, "redotcc")

proc actWilson*(beta: Gscalar): auto =
  beta * initCoeffLeaf(beta.runtime, GaugeActionCoeffs(plaq: 1.0))
proc actSymanzik*(beta: Gscalar): auto =
  beta * initCoeffLeaf(beta.runtime,
    GaugeActionCoeffs(plaq: 1.0 - 8.0 * C1Symanzik, rect: C1Symanzik),
  )
proc actIwasaki*(beta: Gscalar): auto =
  beta * initCoeffLeaf(beta.runtime,
    GaugeActionCoeffs(plaq: 1.0 - 8.0 * C1Iwasaki, rect: C1Iwasaki),
  )
proc actDBW2*(beta: Gscalar): auto =
  beta * initCoeffLeaf(beta.runtime,
    GaugeActionCoeffs(plaq: 1.0 - 8.0 * C1DBW2, rect: C1DBW2),
  )
proc actAdj*(beta: Gscalar, adjFac: Gscalar): auto = beta * adjCoeff(adjFac)

proc actWilson*(grt: GraphRuntime,
                beta: float): Gactcoeff =
  actWilson(toGvalue(grt, beta))

proc actSymanzik*(grt: GraphRuntime,
                  beta: float): Gactcoeff =
  actSymanzik(toGvalue(grt, beta))

proc actIwasaki*(grt: GraphRuntime,
                 beta: float): Gactcoeff =
  actIwasaki(toGvalue(grt, beta))

proc actDBW2*(grt: GraphRuntime,
              beta: float): Gactcoeff =
  actDBW2(toGvalue(grt, beta))

proc actAdj*(grt: GraphRuntime,
             beta: float,
             adjFac: float): Gactcoeff =
  actAdj(toGvalue(grt, beta), toGvalue(grt, adjFac))
