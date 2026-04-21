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
    result = Ggauge(gval: g, runtime: grt)
  else:
    result = Ggauge(gval: x, runtime: grt)
  result.updated

method newOneOf*(x: Ggauge): Gvalue =
  let g = x.gval.newOneOf
  g.zeroGaugeStorage
  result = Ggauge(gval: g, runtime: x.runtime)

method valCopy*(z: Ggauge, x: Ggauge) =
  z.gval.copyGaugeStorage(x.gval)

proc sameGaugeShape(a: Gauge, b: Gauge): bool =
  if a.len != b.len:
    return false
  for i in 0..<a.len:
    if a[i].l != b[i].l:
      return false
  true

method copyCompatible*(prototype: Ggauge, value: Ggauge): bool =
  prototype != nil and value != nil and sameGaugeShape(prototype.gval, value.gval)

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

method retr*(x: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("retr(" & $x & ")")
method adj*(x: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("adj(" & $x & ")")
method norm2*(x: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("norm2(" & $x & ")")
method redot*(x: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("redot(" & $x & "," & $y & ")")
method expDeriv*(b: Gvalue, x: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("expDeriv(" & $b & "," & $x & ")")
method projTAH*(x: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("projTAH(" & $x & ")")

method adjmul*(x: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("adjmul(" & $x & "," & $y & ")")
method muladj*(x: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("muladj(" & $x & "," & $y & ")")
method contractProjTAH*(x: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("contractProjTAH(" & $x & "," & $y & ")")
method axexp*(a: Gvalue, x: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("axexp(" & $a & "," & $x & ")")
method axexpmulyPack*(a: Gvalue, x: Gvalue, y: Gvalue): Gmulti {.base.} =
  raiseErrorBaseMethod("axexpmulyPack(" & $a & "," & $x & "," & $y & ")")
method axexpmuly*(a: Gvalue, x: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("axexpmuly(" & $a & "," & $x & "," & $y & ")")

proc gaugeUpstreamValue*(zb: Gvalue, label: string): Gvalue {.inline.} =
  requireUpstream(zb, label)

proc gaugeUpstreamRetr*(zb: Gvalue, label: string): Gvalue {.inline.} =
  retr(gaugeUpstreamValue(zb, label))

proc negatedGaugeUpstream*(zb: Gvalue, label: string): Gvalue {.inline.} =
  -gaugeUpstreamValue(zb, label)

proc negatedGaugeUpstreamRetr*(zb: Gvalue, label: string): Gvalue {.inline.} =
  -gaugeUpstreamRetr(zb, label)

proc gaugeUpstreamProjTAH*(zb: Gvalue, label: string): Gvalue {.inline.} =
  projTAH(gaugeUpstreamValue(zb, label))

proc sameGaugeBinaryBackward*(zb: Gvalue,
                              label: string,
                              i: int): Gvalue =
  case i
  of 0, 1:
    gaugeUpstreamValue(zb, label)
  else:
    raiseInputIndexError(label, "0 or 1", i)

proc signedGaugeBinaryBackward*(zb: Gvalue,
                                label: string,
                                i: int): Gvalue =
  case i
  of 0:
    gaugeUpstreamValue(zb, label)
  of 1:
    negatedGaugeUpstream(zb, label)
  else:
    raiseInputIndexError(label, "0 or 1", i)

proc scalarGaugeAddBackward*(zb: Gvalue,
                             label: string,
                             i: int): Gvalue =
  case i
  of 0:
    gaugeUpstreamRetr(zb, label)
  of 1:
    gaugeUpstreamValue(zb, label)
  else:
    raiseInputIndexError(label, "0 or 1", i)

proc gaugeScalarSubBackward*(zb: Gvalue,
                             label: string,
                             i: int): Gvalue =
  case i
  of 0:
    gaugeUpstreamValue(zb, label)
  of 1:
    negatedGaugeUpstreamRetr(zb, label)
  else:
    raiseInputIndexError(label, "0 or 1", i)

proc swappedScaledInputBackward*(zb: Gvalue,
                                 z: Gvalue,
                                 label: string,
                                 i: int): Gvalue =
  let view = z.requireBinaryNodeView(label)
  case i
  of 0:
    scaledUpstreamOr(zb, view.y)
  of 1:
    scaledUpstreamOr(zb, view.x)
  else:
    raiseInputIndexError(label, "0 or 1", i)
