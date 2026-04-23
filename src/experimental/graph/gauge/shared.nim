import ../[core, multi]
import ../support/op
import layout, gauge, physics/qcdTypes

# Shared Gauge Helpers

type
  Gauge* = seq[DLatticeColorMatrixV]

  Ggauge* {.final.} = ref object of Gvalue
    gval*: Gauge

template copyGaugeStorage(dst, src: untyped) =
  threads:
    for mu in 0..<dst.len:
      dst[mu] := src[mu]

template reunitGaugeImpl(g: untyped) =
  threads:
    g.projectSU
    threadBarrier()

proc reunitGauge*(g: Gauge) =
  reunitGaugeImpl(g)

proc reunitGauge*(g: auto) =
  reunitGaugeImpl(g)

proc getgauge*(x: Gvalue): Gauge = Ggauge(x).gval

proc requireGauge*(value: Gvalue,
                   label: string): Ggauge =
  if value == nil:
    raiseValueError(label & " requires non-nil gauge value")
  if not (value of Ggauge):
    raiseValueError(label & " expects gauge value, got:\n" & value.nodeRepr)
  Ggauge(value)

proc zeroGaugeStorage*(g: Gauge) =
  threads:
    for mu in 0..<g.len:
      g[mu] := 0.0

proc update*(x: Gvalue, g: Gauge, isZero = false) =
  let u = Ggauge(x)
  if isZero:
    u.gval.zeroGaugeStorage
  else:
    u.gval.copyGaugeStorage(g)
  x.updated

proc toGvalue*(grt: GraphRuntime,
               x: Gauge,
               isZero = false): Ggauge =
  # Use a proc instead of a converter so seq values are not converted implicitly.
  if isZero:
    let g = x.newOneOf
    g.zeroGaugeStorage
    result = Ggauge(gval: g).attachRuntime(grt)
  else:
    result = Ggauge(gval: x).attachRuntime(grt)
  result.updated

proc gaugeNodeLike*(x: Ggauge): Ggauge =
  let g = x.gval.newOneOf
  g.zeroGaugeStorage
  Ggauge(gval: g).attachRuntime(x.runtime)

method newOneOf*(x: Ggauge): Gvalue =
  x.gaugeNodeLike

proc valCopy*(z: Ggauge, x: Ggauge) =
  z.gval.copyGaugeStorage(x.gval)

method valCopy*(z: Ggauge, x: Gvalue) =
  if x == nil or not (x of Ggauge):
    raiseValueError("gauge copy expects gauge value")
  z.valCopy(Ggauge(x))

proc sameGaugeShape(a: Gauge, b: Gauge): bool =
  if a.len != b.len:
    return false
  for i in 0..<a.len:
    if a[i].l != b[i].l:
      return false
  true

proc copyCompatible*(prototype: Ggauge, value: Ggauge): bool =
  prototype != nil and value != nil and sameGaugeShape(prototype.gval, value.gval)

method copyCompatible*(prototype: Ggauge, value: Gvalue): bool =
  prototype != nil and value != nil and value of Ggauge and
    prototype.copyCompatible(Ggauge(value))

method `$`*(x: Ggauge): string =
  let v = x.gval[0][0][0,0]
  result = "Gauge (" & $v.re[0] & ", " & $v.im[0] & ")"

template defineUnaryForward*(forwardName, InputType, ResultType: untyped,
                             label: static[string],
                             forwardBody: untyped) =
  proc forwardName(v: Gvalue) =
    let view {.inject.} = v.requireUnaryNodeView(label)
    let x {.inject.} = InputType(view.x)
    let z {.inject.} = ResultType(v)
    forwardBody

template defineBinaryForward*(forwardName, LeftType, RightType, ResultType: untyped,
                              label: static[string],
                              forwardBody: untyped) =
  proc forwardName(v: Gvalue) =
    let view {.inject.} = v.requireBinaryNodeView(label)
    let x {.inject.} = LeftType(view.x)
    let y {.inject.} = RightType(view.y)
    let z {.inject.} = ResultType(v)
    forwardBody

template mapGaugeSites*(dst: Ggauge, valueExpr: untyped) =
  threads:
    for mu {.inject.} in 0..<dst.gval.len:
      dst.gval[mu] := valueExpr

template mapGaugeElements*(dst: Ggauge, body: untyped) =
  threads:
    for mu {.inject.} in 0..<dst.gval.len:
      for e {.inject.} in dst.gval[mu]:
        body

proc gaugeUpstreamValue*(zb: Gvalue, label: string): Ggauge {.inline.} =
  requireUpstream(zb, Ggauge, label)

proc sameGaugeBinaryBackward*(zb: Gvalue,
                              label: string,
                              i: int): Gvalue =
  case i
  of 0, 1:
    gaugeUpstreamValue(zb, label)
  else:
    raiseInputIndexError(label, "0 or 1", i)
