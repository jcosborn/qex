import ../core
import ../field/types
export types
import layout, gauge, physics/qcdTypes

# Gauge-layer value storage: the gauge bundle and the single-direction field.

type
  Gauge* = seq[DLatticeColorMatrixV]

  Ggauge* = ref object of Gvalue
    ## Graph-owned storage; public writes must mark freshness.
    gval*: Gauge

  Gfield* = GfieldOf[DLatticeColorMatrixV]
    ## One direction of a Gauge. Internal plumbing for cross-direction
    ## expressions and per-direction cotangents.

  Grmat*[n:static[int]] = GfieldOf[DLatticeRealMatrixV[n]]
    ## Graph storage bindings currently cover scalar (1) and local SU3 (8) matrices.
  Grfield* = Grmat[1]
  Grmat8* = Grmat[8]

const cfieldIsGfield* = DColorMatrixV is ColorMatrixN[1, DComplexV]
  ## With Nc = 1 the two field types coincide.

when cfieldIsGfield:
  type Gcfield* = Gfield
else:
  type Gcfield* = GfieldOf[DLatticeComplexMatrixV[1]]

method hasStorage*(x: Ggauge): bool = x.gval.hasFieldStorage

method ensureStorage*(x: Ggauge) =
  if not x.hasStorage:
    x.gval.ensureFieldStorage
    if x.restoreValue != nil: x.restoreValue(x)

method releaseStorage*(x: Ggauge) = x.gval.releaseFieldStorage

method valAlias*(z: Ggauge, x: Gvalue) =
  z.gval = Ggauge(x).gval

template copyGaugeStorage(dst, src: untyped) =
  threads:
    for mu in 0..<dst.len:
      dst[mu] := src[mu]

proc requireSameGaugeShape(dst: Gauge,
                           src: Gauge,
                           label: string) =
  if not sameFieldShape(dst, src):
    raiseValueError(label & " requires matching gauge shapes")

proc requireSameGaugeShape*(left: Ggauge,
                            right: Ggauge,
                            label: string) =
  left.gval.requireSameGaugeShape(right.gval, label)

proc gaugeNodeLike*(x: Ggauge): Ggauge =
  Ggauge(runtime: x.runtime, gval: x.gval.newShape).assignStableNodeId

method bufferProto*(x: Ggauge): Gvalue = x.gaugeNodeLike

method bufferCompatible*(x: Ggauge, y: Gvalue): bool =
  y of Ggauge and sameFieldShape(x.gval,Ggauge(y).gval)

method bindBuffer*(x: Ggauge, buffer: Gvalue) =
  x.gval = Ggauge(buffer).gval

method clearBuffer*(x: Ggauge) = x.gval.zeroFieldStorage

method bufferBytes*(x: Ggauge): int = x.gval.fieldBytes

proc gaugeSnapshot*(x: Ggauge): Gauge =
  if not x.hasStorage:
    discard x.eval
  let storage = x.gval
  let snapshot = storage.newOneOf
  snapshot.copyGaugeStorage(storage)
  result = snapshot

proc zeroGaugeStorage*(g: Gauge) = zeroFieldStorage(g)

proc update*(x: Ggauge, g: Gauge) =
  x.gval.requireSameGaugeShape(g, "gauge update")
  x.ensureStorage
  x.gval.copyGaugeStorage(g)
  x.updated

template mutateGauge*(x: Ggauge, storageName: untyped, body: untyped) =
  block:
    let gaugeNode {.gensym.} = x
    discard gaugeNode.eval
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

proc unitGaugeValue(v: Gvalue) =
  let x = Ggauge(v)
  threads:
    for f in x.gval:
      f := 1.0

proc unitGaugeLike*(x: Ggauge): Ggauge =
  ## Constant identity value restored from shape on demand.
  result = x.gaugeNodeLike
  result.updated
  # updated clears restoreValue; install the hook afterwards.
  result.restoreValue = unitGaugeValue
  result.valueOverride = false

proc sameShapeGaugeNodeLike*(x: Ggauge,
                             y: Ggauge,
                             label: string): Ggauge =
  x.requireSameGaugeShape(y, label)
  x.gaugeNodeLike

method newOneOf*(x: Ggauge): Gvalue =
  x.gaugeNodeLike

method valueLike*(x: Ggauge): Gvalue =
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
  z.ensureStorage
  z.gval.copyGaugeStorage(src.gval)

method copyCompatible*(prototype: Ggauge, value: Gvalue): bool =
  value of Ggauge and sameFieldShape(prototype.gval, Ggauge(value).gval)

method `$`*(x: Ggauge): string =
  if not x.hasStorage:
    return "Gauge (not resident)"
  let v = x.gval[0][0][0,0]
  result = "Gauge (" & $v.re[0] & ", " & $v.im[0] & ")"

proc requireLinkShape*(g: Ggauge, mu: int, f: DLatticeColorMatrixV, label: string) =
  if mu < 0 or mu >= g.gval.len:
    raiseValueError(label & " direction out of range")
  if f.l != g.gval[mu].l:
    raiseValueError(label & " requires matching field shapes")

proc unitFieldLike*(g: Ggauge): Gfield =
  ## Constant identity-matrix field leaf shaped like one direction of g.
  unitField(g.runtime, g.gval[0])

template matrixFieldMethods(T: typedesc, label: static string) =
  fieldMethods(T, label)

  method `$`*(x: T): string =
    if not x.hasStorage:
      return label & " (not resident)"
    let v = x.fval[0][0,0]
    result = label & " (" & $v.re[0] & ", " & $v.im[0] & ")"

matrixFieldMethods(Gfield, "GaugeField")
when not cfieldIsGfield:
  matrixFieldMethods(Gcfield, "ComplexField")
matrixFieldMethods(Grfield, "RealField")
matrixFieldMethods(Grmat8, "RealMatrix8")

type MatrixStorage = DLatticeColorMatrixV | DLatticeComplexMatrixV[1] | DLatticeRealMatrixV[1] | DLatticeRealMatrixV[8]

proc toGvalue*[F:MatrixStorage](grt: GraphRuntime, x: F): GfieldOf[F] =
  toGfield(grt, x)

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

proc requireParityDir*(parity, dir, nd: int, label: string) =
  if parity < 0 or parity > 1:
    raiseValueError(label & " parity must be 0 or 1")
  if dir < 0 or dir >= nd:
    raiseValueError(label & " direction out of range")

proc zeroGaugeStorage*(g: Ggauge) =
  ## Deferred outputs are zeroed when allocated, including subset complements.
  zeroGaugeStorage(g.gval)

template forGaugeSubset*(sub: Subset, body: untyped) =
  ## Run `body` (seeing injected outer index `e`) for each outer site in `sub`.
  ## Partial-output forwards declare bmZero so untouched sites stay zero.
  threads:
    for e {.inject.} in sub:
      body

template gaugeTermSum*(count: int, idx, term: untyped): untyped =
  ## count > 0. Evaluate terms in order, retaining the caller's storage access.
  block:
    var idx = 0
    var total {.noinit.}: evalType(term)
    total := term
    inc idx
    while idx < count:
      total += term
      inc idx
    total

template forGaugeBlend*(g: Gauge, sub, other: Subset, dir: int,
                        complement: bool, body: untyped) =
  ## Inside a threads region; inject mu/e and compile-time active for each site.
  ## The complement flag preserves fused kernels that supply their own base.
  if complement:
    for mu {.inject.} in 0..<g.len:
      if mu != dir:
        for e {.inject.} in g[mu]:
          const active {.inject.} = false
          body
    block:
      let mu {.inject.} = dir
      for e {.inject.} in other:
        const active {.inject.} = false
        body
  block:
    let mu {.inject.} = dir
    for e {.inject.} in sub:
      const active {.inject.} = true
      body
