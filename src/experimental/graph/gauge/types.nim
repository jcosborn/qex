import ../core
import layout, gauge, physics/qcdTypes

# Gauge-layer value storage: the gauge bundle and the single-direction field.

type
  Gauge* = seq[DLatticeColorMatrixV]

  Ggauge* = ref object of Gvalue
    ## Graph-owned storage; public writes must mark freshness.
    gval*: Gauge

  GfieldOf*[F] = ref object of Gvalue
    ## One lattice field of site matrices; graph-owned storage, public writes
    ## must mark freshness.
    fval*: F

  Gfield* = GfieldOf[DLatticeColorMatrixV]
    ## One direction of a Gauge. Internal plumbing for cross-direction
    ## expressions and per-direction cotangents.

  DLatticeCmatV* = Field[VLEN, ColorMatrixN[1, DComplexV]]
    ## Complex scalar field as 1x1 site matrices, so the matrix algebra and
    ## its kernels serve it unchanged (trace, |c|^2 = norm2, exp).

const cfieldIsGfield* = DColorMatrixV is ColorMatrixN[1, DComplexV]
  ## With Nc = 1 the two field types coincide.

when cfieldIsGfield:
  type Gcfield* = Gfield
else:
  type Gcfield* = GfieldOf[DLatticeCmatV]

proc zeroFieldStorage*[V:static[int],T](f: Field[V,T]) =
  if f.s.data.isNil:
    return
  threads:
    f := 0.0

proc ensureFieldStorage*[F](f: var F) =
  if f.s.data.isNil:
    f = f.newOneOf
    f.zeroFieldStorage

proc releaseFieldStorage*[F](f: var F) =
  # Replace the descriptor; aliases keep the old Field and RawMemRef alive.
  if not f.s.data.isNil:
    f = f.newShape

method hasStorage*(x: Ggauge): bool =
  for f in x.gval:
    if f.s.data.isNil:
      return false
  true

method ensureStorage*(x: Ggauge) =
  if x.hasStorage:
    return
  var fields = newSeq[DLatticeColorMatrixV](x.gval.len)
  for i, f in x.gval:
    fields[i] = f
    fields[i].ensureFieldStorage
  x.gval = fields
  if x.restoreValue != nil:
    x.restoreValue(x)

method releaseStorage*(x: Ggauge) =
  var resident = false
  for f in x.gval:
    if not f.s.data.isNil:
      resident = true
      break
  if not resident:
    return
  var fields = newSeq[DLatticeColorMatrixV](x.gval.len)
  for i, f in x.gval:
    fields[i] = f.newShape
  x.gval = fields

method valAlias*(z: Ggauge, x: Gvalue) =
  z.gval = Ggauge(x).gval

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

method bufferProto*(x: Ggauge): Gvalue =
  var g = newSeq[DLatticeColorMatrixV](x.gval.len)
  for i, f in x.gval:
    g[i] = f.newShape
  Ggauge(runtime:x.runtime,gval:g).assignStableNodeId

method bufferCompatible*(x: Ggauge, y: Gvalue): bool =
  y of Ggauge and sameGaugeShape(x.gval,Ggauge(y).gval)

method bindBuffer*(x: Ggauge, buffer: Gvalue) =
  x.gval = Ggauge(buffer).gval

method clearBuffer*(x: Ggauge) =
  threads:
    for f in x.gval:
      f := 0.0

method bufferBytes*(x: Ggauge): int =
  for f in x.gval:
    result += f.s.bytes

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
  if not x.hasStorage:
    discard x.eval
  let storage = x.gval
  let snapshot = storage.newOneOf
  snapshot.copyGaugeStorage(storage)
  result = snapshot

proc zeroGaugeStorage*(g: Gauge) =
  threads:
    for mu in 0..<g.len:
      if not g[mu].s.data.isNil:
        g[mu] := 0.0

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

proc gaugeNodeLike*(x: Ggauge): Ggauge =
  var g = newSeq[DLatticeColorMatrixV](x.gval.len)
  for i, f in x.gval:
    g[i] = f.newShape
  Ggauge(runtime: x.runtime, gval: g).assignStableNodeId

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
  value of Ggauge and sameGaugeShape(prototype.gval, Ggauge(value).gval)

method `$`*(x: Ggauge): string =
  if not x.hasStorage:
    return "Gauge (not resident)"
  let v = x.gval[0][0][0,0]
  result = "Gauge (" & $v.re[0] & ", " & $v.im[0] & ")"

proc requireSameFieldShape*[F](x, y: GfieldOf[F], label: string) =
  if x.fval.l != y.fval.l:
    raiseValueError(label & " requires matching field shapes")

proc requireLinkShape*(g: Ggauge, mu: int, f: DLatticeColorMatrixV, label: string) =
  if mu < 0 or mu >= g.gval.len:
    raiseValueError(label & " direction out of range")
  if f.l != g.gval[mu].l:
    raiseValueError(label & " requires matching field shapes")

proc fieldNodeLike*[F](x: GfieldOf[F]): GfieldOf[F] =
  let f = x.fval.newShape
  GfieldOf[F](runtime: x.runtime, fval: f).assignStableNodeId

proc sameShapeFieldNodeLike*[F](x, y: GfieldOf[F], label: string): GfieldOf[F] =
  x.requireSameFieldShape(y, label)
  x.fieldNodeLike

proc unitFieldValue[F](v: Gvalue) =
  let x = GfieldOf[F](v)
  threads:
    x.fval := 1.0

proc unitField*[F](grt: GraphRuntime, proto: F): GfieldOf[F] =
  ## Constant identity value restored from shape on demand.
  let f = proto.newShape
  result = GfieldOf[F](runtime: grt, fval: f).assignStableNodeId
  result.updated
  # updated clears restoreValue; install the hook afterwards.
  result.restoreValue = unitFieldValue[F]
  result.valueOverride = false

proc unitFieldLike*(g: Ggauge): Gfield =
  ## Constant identity-matrix field leaf shaped like one direction of g.
  unitField(g.runtime, g.gval[0])

template fieldMethods(T: typedesc, label: static string) =
  method bufferProto*(x: T): Gvalue =
    T(runtime:x.runtime,fval:x.fval.newShape).assignStableNodeId

  method bufferCompatible*(x: T, y: Gvalue): bool =
    y of T and x.fval.l == T(y).fval.l

  method bindBuffer*(x: T, buffer: Gvalue) =
    x.fval = T(buffer).fval

  method clearBuffer*(x: T) =
    threads:
      x.fval := 0.0

  method bufferBytes*(x: T): int = x.fval.s.bytes

  method hasStorage*(x: T): bool = not x.fval.s.data.isNil

  method ensureStorage*(x: T) =
    if not x.hasStorage:
      x.fval.ensureFieldStorage
      if x.restoreValue != nil:
        x.restoreValue(x)

  method releaseStorage*(x: T) =
    x.fval.releaseFieldStorage

  method valAlias*(z: T, x: Gvalue) =
    z.fval = T(x).fval

  method newOneOf*(x: T): Gvalue =
    x.fieldNodeLike

  method valueLike*(x: T): Gvalue =
    x.fieldNodeLike

  method zeroLike*(x: T): Gvalue =
    result = x.fieldNodeLike
    result.staticZeroLeaf = true

  method isZero*(x: T): bool =
    ## Zero leaves are marked when constructed; other fields are not scanned.
    x.staticZeroLeaf

  method valCopy*(z: T, x: Gvalue) =
    let src = T(x)
    if z.fval.l != src.fval.l:
      raiseValueError(label & " copy requires matching field shapes")
    z.ensureStorage
    threads:
      z.fval := src.fval

  method copyCompatible*(prototype: T, value: Gvalue): bool =
    value of T and prototype.fval.l == T(value).fval.l

  method `$`*(x: T): string =
    if not x.hasStorage:
      return label & " (not resident)"
    let v = x.fval[0][0,0]
    result = label & " (" & $v.re[0] & ", " & $v.im[0] & ")"

fieldMethods(Gfield, "GaugeField")
when not cfieldIsGfield:
  fieldMethods(Gcfield, "ComplexField")

type MatrixStorage = DLatticeColorMatrixV | DLatticeCmatV

proc toGvalue*[F:MatrixStorage](grt: GraphRuntime, x: F): GfieldOf[F] =
  let f = x.newOneOf
  threads:
    f := x
  result = GfieldOf[F](runtime: grt, fval: f).assignStableNodeId
  result.updated

proc update*[F](x: GfieldOf[F], f: F) =
  if x.fval.l != f.l:
    raiseValueError("field update requires matching field shapes")
  x.ensureStorage
  threads:
    x.fval := f
  x.updated

template mutateField*[F](x: GfieldOf[F], storageName: untyped, body: untyped) =
  block:
    let node {.gensym.} = x
    discard node.eval
    let storageName {.inject.} = node.fval
    try:
      body
    finally:
      node.updated

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
