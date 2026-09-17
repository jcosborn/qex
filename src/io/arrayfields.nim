## Map C-order arrays to real QEX field channels using global coordinates.
## Axes follow `lo.physGeom` with the last lattice axis fastest; QIO site order
## has axis 0 fastest. Call outside `threads:`. Every operation is collective
## over `lo.comm`: each rank validates the descriptor, then transfers only its
## own sites through the parallel I/O maps, so no rank holds a global array.
## A failed transfer has no collective recovery.
import qex
import arrays, parallelIo
import std/[json, os, sequtils]

proc requireShape[F](fs: seq[F], spec: JsonNode): seq[int] =
  static: doAssert F is SLatticeRealV or F is DLatticeRealV
  result = arrayShape(spec)
  if fs.len == 0 or result != (@[fs.len] & fs[0].l.physGeom):
    raise newException(ValueError,"field array shape differs from its channels and lattice")
  for f in fs:
    if f.l != fs[0].l:
      raise newException(ValueError,"field array layouts differ")

proc fileMap(lo: Layout[VLEN]): WriteMap =
  ## Every rank transfers one contiguous chunk of each channel's C-order block.
  lo.setupWrite(lo.physGeom, newSeq[int](lo.nDim), toSeq(0..<lo.comm.commsize), corder = true)

template readChannels(lo: Layout[VLEN]; path: string; meta: ArrayMeta; B: typedesc; n: int; body: untyped) =
  ## Reads channel c into buf, in local scalar-site order, then runs body.
  if not fileExists(path) or getFileSize(path) != meta.bytes:
    raise newException(ValueError, path & " must hold exactly " & $meta.bytes & " bytes")
  let wm = fileMap(lo)
  var pr = lo.comm.openRead(path)
  if system.cpuEndian == bigEndian and sizeof(B) > 1: pr.setSwap(8*sizeof(B))
  var buf {.inject.} = newSeq[B](lo.nSites)
  for c {.inject.} in 0..<n:
    pr.read(wm, sizeof(B), addr buf[0])
    body
  pr.close

proc loadFields*[F](fs: seq[F], dir: string, spec: JsonNode) =
  ## File shape is [channels] followed by the global lattice dimensions.
  discard requireShape(fs,spec)
  let meta = arrayMeta(spec)
  let lo = fs[0].l
  template load(B: typedesc) =
    readChannels(lo, dir/meta.file, meta, B, fs.len):
      threads:
        for s in lo.sites: fs[c]{s} := buf[s]
  case meta.dtype
  of dtFloat32: load(float32)
  of dtFloat64: load(float64)
  of dtUint8: raise newException(ValueError, "floating fields require floating data")

proc saveFields*[F](fs: seq[F], dir: string, spec: JsonNode) =
  ## Each rank writes its sites of every channel into the pre-sized file.
  discard requireShape(fs,spec)
  let meta = arrayMeta(spec)
  if meta.dtype == dtUint8:
    raise newException(ValueError, "floating fields require floating data")
  let lo = fs[0].l
  let path = dir/meta.file
  if lo.comm.isMaster and path.parentDir.len > 0: createDir(path.parentDir)
  let wm = fileMap(lo)
  var pw = lo.comm.openCreate(path, meta.bytes)
  if system.cpuEndian == bigEndian: pw.setSwap(8*meta.width)
  template save(B: typedesc) =
    var buf = newSeq[B](lo.nSites)
    for c in 0..<fs.len:
      threads:
        for s in lo.sites: buf[s] := fs[c]{s}
      pw.write(wm, sizeof(B), addr buf[0])
  case meta.dtype
  of dtFloat32: save(float32)
  of dtFloat64: save(float64)
  of dtUint8: discard # Rejected above.
  pw.close

proc readMask*(lo: Layout[VLEN], dir: string, spec: JsonNode): SLatticeRealV =
  ## A byte array over the global lattice; nonzero samples select sites.
  if arrayShape(spec) != lo.physGeom:
    raise newException(ValueError,"mask array shape differs from its lattice")
  let meta = arrayMeta(spec)
  if meta.dtype != dtUint8:
    raise newException(ValueError, "mask arrays require uint8 data")
  result = lo.RealS()
  let mask = result
  readChannels(lo, dir/meta.file, meta, uint8, 1):
    threads:
      for s in lo.sites: mask{s} := (if buf[s] != 0: 1'f32 else: 0'f32)
