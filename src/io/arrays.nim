## Dense array interchange with NumPy/JAX (ndarray.tofile / numpy.fromfile):
## a JSON descriptor {file, dtype, shape, byte_order, order} and a raw payload.
## `file` is nonempty and relative without ".."; dtype is float32, float64 or
## uint8; the last shape axis is fastest; the payload starts at byte zero and
## holds exactly arraySize(spec) elements with no header, framing or checksum.
## The application owns the manifest schema around these descriptors.
## SciDAC/LIME lattice records use io/reader and io/writer instead.
import std/[json, os, strutils]

type
  Dtype* = enum dtFloat32, dtFloat64, dtUint8
  ArrayMeta* = tuple[file: string, dtype: Dtype, size, bytes, width: int]

proc member(spec: JsonNode; key: string; kind: JsonNodeKind): JsonNode =
  if spec == nil or spec.kind != JObject:
    raise newException(ValueError, "array metadata must be an object")
  if not spec.hasKey(key) or spec[key].kind != kind:
    raise newException(ValueError, "array metadata requires " & $kind & " member '" & key & "'")
  spec[key]

proc arrayShape*(spec: JsonNode): seq[int] =
  for dim in spec.member("shape",JArray):
    if dim.kind != JInt or dim.getBiggestInt < 0 or dim.getBiggestInt > high(int):
      raise newException(ValueError, "array extents must be nonnegative integers representable by int")
    result.add dim.getInt

proc arraySize*(spec: JsonNode): int =
  ## An empty shape describes one scalar; any zero extent describes no elements.
  let shape = arrayShape(spec)
  for dim in shape:
    if dim == 0: return 0
  result = 1
  for dim in shape:
    if result > high(int) div dim:
      raise newException(ValueError, "array element count exceeds int")
    result *= dim

proc arrayMeta*(spec: JsonNode): ArrayMeta =
  ## The validated descriptor: relative file, element type and width, counts.
  result.file = spec.member("file",JString).getStr
  if result.file.len == 0:
    raise newException(ValueError, "array file name must not be empty")
  if result.file.isAbsolute or ".." in result.file.split({DirSep, AltSep}):
    raise newException(ValueError, "array file must be relative without '..' components: " & result.file)
  if spec.member("byte_order",JString).getStr != "little":
    raise newException(ValueError, "array byte_order must be little")
  if spec.member("order",JString).getStr != "C":
    raise newException(ValueError, "array order must be C")
  let dtype = spec.member("dtype",JString).getStr
  case dtype
  of "float32":
    result.dtype = dtFloat32
    result.width = 4
  of "float64":
    result.dtype = dtFloat64
    result.width = 8
  of "uint8":
    result.dtype = dtUint8
    result.width = 1
  else:
    raise newException(ValueError, "unsupported array dtype: " & dtype)
  result.size = arraySize(spec)
  if result.size > high(int) div result.width:
    raise newException(ValueError, "array byte count exceeds int")
  result.bytes = result.size*result.width

proc compatible[T: SomeFloat | uint8](dtype: Dtype) =
  when T is uint8:
    if dtype != dtUint8: raise newException(ValueError, "uint8 arrays require uint8 data")
  else:
    if dtype == dtUint8: raise newException(ValueError, "floating arrays require floating data")

proc load32(data: string; pos: int): float32 =
  var bits: uint32
  for b in 0..<4: bits = bits or (uint32(ord(data[pos+b])) shl (8*b))
  cast[float32](bits)

proc load64(data: string; pos: int): float64 =
  var bits: uint64
  for b in 0..<8: bits = bits or (uint64(ord(data[pos+b])) shl (8*b))
  cast[float64](bits)

proc readArray*[T: SomeFloat | uint8](dir: string; spec: JsonNode): seq[T] =
  ## Float dtypes convert to T. uint8 is read only as uint8.
  let meta = arrayMeta(spec)
  compatible[T](meta.dtype)
  let path = dir / meta.file
  let data = readFile(path)
  if data.len != meta.bytes:
    raise newException(ValueError, path & " has " & $data.len & " bytes; expected " & $meta.bytes)
  result = newSeq[T](meta.size)
  when T is uint8:
    for i in 0..<meta.size: result[i] = uint8(ord(data[i]))
  else:
    case meta.dtype
    of dtFloat32:
      for i in 0..<meta.size: result[i] = T(load32(data,4*i))
    of dtFloat64:
      for i in 0..<meta.size: result[i] = T(load64(data,8*i))
    of dtUint8: discard # Rejected before reading.

proc writeArray*[T: SomeFloat | uint8](dir: string; spec: JsonNode; data: openArray[T]; write = true) =
  ## The metadata dtype determines output precision; shape determines byte count.
  ## `write = false` validates the descriptor and the data length only.
  let meta = arrayMeta(spec)
  compatible[T](meta.dtype)
  if data.len != meta.size:
    raise newException(ValueError, "array data length " & $data.len & " differs from shape size " & $meta.size)
  if not write: return
  var bytes = newString(meta.bytes)
  when T is uint8:
    for i, value in data: bytes[i] = char(value)
  else:
    case meta.dtype
    of dtFloat32:
      for i, value in data:
        let bits = cast[uint32](float32(value))
        for b in 0..<4: bytes[4*i+b] = char((bits shr (8*b)) and 255)
    of dtFloat64:
      for i, value in data:
        let bits = cast[uint64](float64(value))
        for b in 0..<8: bytes[8*i+b] = char((bits shr (8*b)) and 255)
    of dtUint8: discard # Rejected before writing.
  let path = dir / meta.file
  if path.parentDir.len > 0: createDir(path.parentDir)
  writeFile(path,bytes)

proc readManifest*(dir: string): JsonNode =
  ## The caller owns schema and application-specific validation.
  result = parseFile(dir / "manifest.json")
  if result.kind != JObject:
    raise newException(ValueError, "array manifest must be an object")
