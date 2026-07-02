import ../core
import layout, gauge, physics/qcdTypes

# Shared Gauge Helpers

type
  Gauge* = seq[DLatticeColorMatrixV]

  Ggauge* = ref object of Gvalue
    ## Graph-owned storage; public writes must mark freshness.
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
  # Project each link back onto its gauge group. SU(1) is trivial ({1}), so for
  # Nc==1 (U(1)) reunitize to the unit circle with projectU; otherwise projectSU.
  const nc = g[0][0].nrows
  threads:
    when nc == 1:
      g.projectU
    else:
      g.projectSU
    threadBarrier()

proc checkUnitary*(g: Gauge): tuple[avg, max: float] =
  ## Mean/max link distance from U(1) or SU(N); does not modify g.
  const nc = g[0][0].nrows
  var a, m: float
  threads:
    let d = when nc == 1: g.checkU else: g.checkSU
    threadMaster:
      a = d.avg
      m = d.max
  (avg: a, max: m)

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

method zeroLike*(x: Ggauge): Gvalue =
  result = x.gaugeNodeLike
  result.staticZeroLeaf = true

method isZero*(x: Ggauge): bool =
  ## Gauge zero leaves are marked when constructed; other gauges are not scanned.
  x.staticZeroLeaf

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

proc paritySubset*(g: Gauge, parity: int): Subset =
  ## Even (parity 0) / odd (parity 1) subset of the layout backing `g`.
  g[0].l.getSubset(if parity == 0: "even" else: "odd")

proc zeroGaugeStorage*(g: Ggauge) =
  ## Zero once; subset kernels leave off-subset entries zero across evaluations.
  threads:
    for mu in 0..<g.gval.len:
      g.gval[mu] := 0.0

template forGaugeSubset*(sub: Subset, body: untyped) =
  ## Run `body` (seeing injected outer index `e`) for each outer site in `sub`.
  ## Pair with a one-time `zeroGaugeStorage` on the output so off-subset stays 0.
  threads:
    for e {.inject.} in sub:
      body
