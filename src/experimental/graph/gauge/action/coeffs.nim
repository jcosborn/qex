import ../../[core, scalar]
import ../../core/base
import ../../scalar/types
import ../../support/op
import layout, gauge, physics/qcdTypes
import domain

type
  Gactcoeff* {.final.} = ref object of Gvalue
    cval*: GaugeActionCoeffs

proc coeffNodeLike(anchor: Gvalue): Gactcoeff =
  Gactcoeff().attachRuntime(anchor.runtime)

proc requireActCoeff*(value: Gvalue,
                      label: string): Gactcoeff =
  if not (value of Gactcoeff):
    raiseValueError(
      label & " expects gauge-action coefficient value, got:\n" &
      value.nodeRepr)
  Gactcoeff(value)

proc update*(x: Gactcoeff, c: GaugeActionCoeffs) =
  x.cval = c
  x.updated

proc toGvalue*(grt: GraphRuntime,
               x: GaugeActionCoeffs): Gactcoeff =
  result = Gactcoeff(
    cval: x).attachRuntime(grt)
  result.updated

method newOneOf*(x: Gactcoeff): Gvalue =
  result = Gactcoeff().attachRuntime(x.runtime)
method valCopy*(z: Gactcoeff, x: Gvalue) =
  z.cval = x.requireActCoeff("gauge-action coefficient copy").cval
method copyCompatible*(prototype: Gactcoeff, value: Gvalue): bool =
  value of Gactcoeff
method `$`*(x: Gactcoeff): string = $x.cval

proc raiseCoeffBackwardUnsupported*(label: string) {.noreturn.} =
  raiseUnsupportedPath(label, "derivative with respect to gauge-action coefficients")

proc `*`*(x: Gscalar, y: Gactcoeff): Gactcoeff
proc redot*(x: Gactcoeff, y: Gactcoeff): Gscalar

method scaleLike*(contribution: Gactcoeff, upstream: Gvalue): Gvalue =
  upstream.requireScalar("gauge-action coefficient scale upstream") * contribution

proc adjCoeffb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard dep
  unaryBackwardCase("adjCoeff backward", i,
    block:
      return redot(requireUpstream(zb, Gactcoeff, "adjCoeff backward"),
                   toGvalue(z.runtime, GaugeActionCoeffs(adjplaq: 1.0))))

proc adjCoefff(v: Gvalue) =
  let view = v.requireUnaryNodeView(Gscalar, "adjCoeff forward")
  let x = view.x
  let z = v.requireActCoeff("adjCoeff forward result")
  z.cval = GaugeActionCoeffs(plaq: 1.0, adjplaq: x.sval)

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

proc mulscf(v: Gvalue) =
  let view = v.requireBinaryNodeView(Gscalar, Gactcoeff, "s*c forward")
  let x = view.x
  let y = view.y
  let z = v.requireActCoeff("s*c forward result")
  z.cval = x.sval * y.cval

let mulsc = newGfunc(forward = mulscf, backward = mulscb, name = "s*c")

proc `*`*(x: Gscalar, y: Gactcoeff): Gactcoeff =
  graphNode(coeffNodeLike(x), @[Gvalue(x), Gvalue(y)], mulsc, "s*c")

proc redotccb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let view = z.requireBinaryNodeView(Gactcoeff, Gactcoeff, "redotcc backward")
  discard dep
  binaryBackwardCase("redotcc backward", i,
    block:
      return scaledUpstreamOr(zb, Gscalar, view.y, "redotcc backward"),
    block:
      return scaledUpstreamOr(zb, Gscalar, view.x, "redotcc backward"))

proc redotccf(v: Gvalue) =
  let view = v.requireBinaryNodeView(Gactcoeff, Gactcoeff, "redotcc forward")
  let x = view.x
  let y = view.y
  let z = v.requireScalar("redotcc forward result")
  var t = 0.0
  for a, b in fields(x.cval, y.cval):
    t += a * b
  z.sval = t

let redotcc = newGfunc(forward = redotccf, backward = redotccb, name = "redotcc")

proc redot*(x: Gactcoeff, y: Gactcoeff): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x), Gvalue(y)], redotcc, "redotcc")

proc actWilson*(beta: Gscalar): Gactcoeff =
  beta * toGvalue(beta.runtime, GaugeActionCoeffs(plaq: 1.0))
proc actSymanzik*(beta: Gscalar): Gactcoeff =
  beta * toGvalue(beta.runtime,
    GaugeActionCoeffs(plaq: 1.0 - 8.0 * C1Symanzik, rect: C1Symanzik),
  )
proc actIwasaki*(beta: Gscalar): Gactcoeff =
  beta * toGvalue(beta.runtime,
    GaugeActionCoeffs(plaq: 1.0 - 8.0 * C1Iwasaki, rect: C1Iwasaki),
  )
proc actDBW2*(beta: Gscalar): Gactcoeff =
  beta * toGvalue(beta.runtime,
    GaugeActionCoeffs(plaq: 1.0 - 8.0 * C1DBW2, rect: C1DBW2),
  )
proc actAdj*(beta: Gscalar, adjFac: Gscalar): Gactcoeff = beta * adjCoeff(adjFac)
