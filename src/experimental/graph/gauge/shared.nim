import ../[core, multi]
import ../core/base
import layout, gauge, physics/qcdTypes

# Shared Gauge Helpers

type
  Gauge* = seq[DLatticeColorMatrixV]

  Ggauge* {.final.} = ref object of Gvalue
    gval: Gauge

proc requireGauge*(value: Gvalue,
                   label: string): Ggauge =
  if not (value of Ggauge):
    raiseValueError(label & " expects gauge value, got:\n" & value.nodeRepr)
  Ggauge(value)

template copyGaugeStorage(dst, src: untyped) =
  threads:
    for mu in 0..<dst.len:
      dst[mu] := src[mu]

template unsafeGaugeStorage*(x: Ggauge): untyped =
  ## Graph-owned mutable gauge storage. Prefer `gaugeSnapshot` for snapshots and
  ## `mutateGauge`/`update` for public writes so freshness is marked.
  x.gval

proc sameGaugeShape(a: Gauge, b: Gauge): bool =
  if a.len != b.len:
    return false
  for i in 0..<a.len:
    if a[i].l != b[i].l:
      return false
  true

proc requireSameGaugeShape(dst: Gauge,
                           src: Gauge,
                           label: string) =
  if not sameGaugeShape(dst, src):
    raiseValueError(label & " requires matching gauge shapes")

proc requireSameGaugeShape*(left: Ggauge,
                            right: Ggauge,
                            label: string) =
  left.gval.requireSameGaugeShape(right.gval, label)

proc reunitGauge*(g: Gauge) =
  threads:
    g.projectSU
    threadBarrier()

proc gaugeSnapshot*(x: Ggauge): Gauge =
  let storage = x.unsafeGaugeStorage
  let snapshot = storage.newOneOf
  snapshot.copyGaugeStorage(storage)
  result = snapshot

proc zeroGaugeStorage*(g: Gauge) =
  threads:
    for mu in 0..<g.len:
      g[mu] := 0.0

proc update*(x: Ggauge, g: Gauge) =
  x.gval.requireSameGaugeShape(g, "gauge update")
  x.gval.copyGaugeStorage(g)
  x.updated

template mutateGauge*(x: Ggauge, storageName: untyped, body: untyped) =
  block:
    let gaugeNode {.gensym.} = x
    let storageName {.inject.} = gaugeNode.unsafeGaugeStorage
    try:
      body
    finally:
      gaugeNode.updated

proc toGvalue*(grt: GraphRuntime,
               x: Gauge): Ggauge =
  # Use a proc instead of a converter so seq values are not converted implicitly.
  let g = x.newOneOf
  g.copyGaugeStorage(x)
  result = Ggauge(gval: g).attachRuntime(grt)
  result.updated

proc gaugeNodeLike*(x: Ggauge): Ggauge =
  let g = x.gval.newOneOf
  g.zeroGaugeStorage
  Ggauge(gval: g).attachRuntime(x.runtime)

proc sameShapeGaugeNodeLike*(x: Ggauge,
                             y: Ggauge,
                             label: string): Ggauge =
  x.requireSameGaugeShape(y, label)
  x.gaugeNodeLike

method newOneOf*(x: Ggauge): Gvalue =
  x.gaugeNodeLike

proc valCopy*(z: Ggauge, x: Ggauge) =
  z.gval.requireSameGaugeShape(x.gval, "gauge copy")
  z.gval.copyGaugeStorage(x.gval)

method valCopy*(z: Ggauge, x: Gvalue) =
  z.valCopy(x.requireGauge("gauge copy"))

proc copyCompatible*(prototype: Ggauge, value: Ggauge): bool =
  sameGaugeShape(prototype.gval, value.gval)

method copyCompatible*(prototype: Ggauge, value: Gvalue): bool =
  value of Ggauge and prototype.copyCompatible(Ggauge(value))

method `$`*(x: Ggauge): string =
  let v = x.gval[0][0][0,0]
  result = "Gauge (" & $v.re[0] & ", " & $v.im[0] & ")"

template mapGaugeSites*(dst: Ggauge, valueExpr: untyped) =
  threads:
    for mu {.inject.} in 0..<dst.gval.len:
      dst.gval[mu] := valueExpr

template mapGaugeElements*(dst: Ggauge, body: untyped) =
  threads:
    for mu {.inject.} in 0..<dst.gval.len:
      for e {.inject.} in dst.gval[mu]:
        body
