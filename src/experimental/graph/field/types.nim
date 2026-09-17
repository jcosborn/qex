## Graph-owned QEX fields and collections of fields.
import qex
import ../core

type GfieldOf*[F] = ref object of Gvalue
  ## Public writes must mark freshness with `updated`.
  fval*: F

proc sameFieldShape*[F: Field](a, b: F): bool = a.l == b.l

proc sameFieldShape*[F: Field](a, b: seq[F]): bool =
  if a.len != b.len: return false
  for i in 0..<a.len:
    if not sameFieldShape(a[i], b[i]): return false
  true

proc requireSameFieldShape*[F](x, y: GfieldOf[F], label: string) =
  if not sameFieldShape(x.fval, y.fval):
    raiseValueError(label & " requires matching field shapes")

proc newShape*[F: Field](fs: seq[F]): seq[F] =
  result = newSeq[F](fs.len)
  for i, f in fs: result[i] = f.newShape

proc hasFieldStorage[F: Field](f: F): bool = not f.s.data.isNil

proc hasFieldStorage*[F: Field](fs: seq[F]): bool =
  for f in fs:
    if not f.hasFieldStorage: return false
  true

proc fieldBytes[F: Field](f: F): int = f.s.bytes

proc fieldBytes*[F: Field](fs: seq[F]): int =
  for f in fs: result += f.fieldBytes

proc zeroFieldStorage*[F: Field](f: F) =
  if not f.hasFieldStorage: return
  threads:
    f := 0

proc zeroFieldStorage*[F: Field](fs: seq[F]) =
  threads:
    for f in fs:
      if f.hasFieldStorage: f := 0

proc ensureFieldStorage*[F: Field](f: var F) =
  if not f.hasFieldStorage:
    f = f.newOneOf
    f.zeroFieldStorage

proc ensureFieldStorage*[F: Field](fs: var seq[F]) =
  if fs.hasFieldStorage: return
  var fields = newSeq[F](fs.len)
  for i, f in fs:
    fields[i] = f
    fields[i].ensureFieldStorage
  fs = fields

proc releaseFieldStorage*[F: Field](f: var F) =
  # Replace the descriptor; aliases keep the old Field and RawMemRef alive.
  if f.hasFieldStorage: f = f.newShape

proc releaseFieldStorage*[F: Field](fs: var seq[F]) =
  for f in fs:
    if f.hasFieldStorage:
      fs = fs.newShape
      return

proc copyFieldStorage*[F: Field](dst, src: F) =
  if not sameFieldShape(dst, src):
    raiseValueError("field copy requires matching field shapes")
  threads:
    dst := src

proc copyFieldStorage*[F: Field](dst, src: seq[F]) =
  if not sameFieldShape(dst, src):
    raiseValueError("field copy requires matching field shapes")
  threads:
    for i in 0..<dst.len: dst[i] := src[i]

proc fieldNodeLike*[F](x: GfieldOf[F]): GfieldOf[F] =
  let f = x.fval.newShape
  GfieldOf[F](runtime: x.runtime, fval: f).assignStableNodeId

proc sameShapeFieldNodeLike*[F](x, y: GfieldOf[F], label: string): GfieldOf[F] =
  x.requireSameFieldShape(y, label)
  x.fieldNodeLike

proc unitFieldValue[F](v: Gvalue) =
  let x = GfieldOf[F](v)
  threads:
    when F is seq:
      for f in x.fval: f := 1
    else:
      x.fval := 1

proc unitField*[F](rt: GraphRuntime, proto: F): GfieldOf[F] =
  ## Constant scalar one or matrix identity, restored from shape on demand.
  let f = proto.newShape
  result = GfieldOf[F](runtime: rt, fval: f).assignStableNodeId
  result.updated
  # updated clears restoreValue; install the hook afterwards.
  result.restoreValue = unitFieldValue[F]
  result.valueOverride = false

proc toGfield*[F](rt: GraphRuntime, f: F): GfieldOf[F] =
  ## Copy storage into a new leaf; collection channels share one layout.
  when F is seq:
    if f.len == 0: raiseValueError("field graph value requires channels")
    for c in f:
      if c.l != f[0].l: raiseValueError("field graph value layouts differ")
  let z = f.newOneOf
  z.copyFieldStorage(f)
  result = GfieldOf[F](runtime: rt, fval: z).assignStableNodeId
  result.updated

proc update*[F](x: GfieldOf[F], f: F) =
  if not sameFieldShape(x.fval, f):
    raiseValueError("field update requires matching field shapes")
  x.ensureStorage
  x.fval.copyFieldStorage(f)
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

template fieldMethods*(T: typedesc, label: static string) =
  ## Instantiate erased dispatch for each concrete storage specialization.
  method bufferProto*(x: T): Gvalue = x.fieldNodeLike
  method bufferCompatible*(x: T, y: Gvalue): bool =
    y of T and sameFieldShape(x.fval, T(y).fval)
  method bindBuffer*(x: T, buffer: Gvalue) = x.fval = T(buffer).fval
  method clearBuffer*(x: T) = x.fval.zeroFieldStorage
  method bufferBytes*(x: T): int = x.fval.fieldBytes
  method hasStorage*(x: T): bool = x.fval.hasFieldStorage
  method ensureStorage*(x: T) =
    if not x.hasStorage:
      x.fval.ensureFieldStorage
      if x.restoreValue != nil: x.restoreValue(x)
  method releaseStorage*(x: T) = x.fval.releaseFieldStorage
  method valAlias*(z: T, x: Gvalue) = z.fval = T(x).fval
  method newOneOf*(x: T): Gvalue = x.fieldNodeLike
  method valueLike*(x: T): Gvalue = x.fieldNodeLike
  method zeroLike*(x: T): Gvalue =
    result = x.fieldNodeLike
    result.staticZeroLeaf = true
  method oneLike*(x: T): Gvalue = unitField(x.runtime, x.fval)
  method isZero*(x: T): bool = x.staticZeroLeaf
  method valCopy*(z: T, x: Gvalue) =
    if not z.copyCompatible(x):
      raiseValueError(label & " copy requires matching field shapes")
    z.ensureStorage
    z.fval.copyFieldStorage(T(x).fval)
  method copyCompatible*(x: T, y: Gvalue): bool =
    y of T and sameFieldShape(x.fval, T(y).fval)
