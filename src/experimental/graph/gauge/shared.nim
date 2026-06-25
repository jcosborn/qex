import ../core
import layout, gauge, physics/qcdTypes

# Shared Gauge Helpers

type
  Gauge* = seq[DLatticeColorMatrixV]

  Ggauge* {.final.} = ref object of Gvalue
    ## `gval` is graph-owned mutable gauge storage. Prefer `gaugeSnapshot` for
    ## snapshots and `mutateGauge`/`update` for public writes so freshness is marked.
    gval*: Gauge

template copyGaugeStorage(dst, src: untyped) =
  threads:
    for mu in 0..<dst.len:
      dst[mu] := src[mu]

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
  let storage = x.gval
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
    let storageName {.inject.} = gaugeNode.gval
    try:
      body
    finally:
      gaugeNode.updated

proc toGvalue*(grt: GraphRuntime,
               x: Gauge): Ggauge =
  # Use a proc instead of a converter so seq values are not converted implicitly.
  let g = x.newOneOf
  g.copyGaugeStorage(x)
  result = Ggauge(runtime: grt, gval: g).assignStableNodeId
  result.updated

proc gaugeNodeLike*(x: Ggauge): Ggauge =
  let g = x.gval.newOneOf
  g.zeroGaugeStorage
  Ggauge(runtime: x.runtime, gval: g).assignStableNodeId

proc sameShapeGaugeNodeLike*(x: Ggauge,
                             y: Ggauge,
                             label: string): Ggauge =
  x.requireSameGaugeShape(y, label)
  x.gaugeNodeLike

method newOneOf*(x: Ggauge): Gvalue =
  x.gaugeNodeLike

method valCopy*(z: Ggauge, x: Gvalue) =
  let src = Ggauge(x)
  z.gval.requireSameGaugeShape(src.gval, "gauge copy")
  z.gval.copyGaugeStorage(src.gval)

method copyCompatible*(prototype: Ggauge, value: Gvalue): bool =
  value of Ggauge and sameGaugeShape(prototype.gval, Ggauge(value).gval)

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
