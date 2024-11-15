import parallelIo
import tables

proc modReadString(pr: var ParallelReader): string =
  var bytes = pr.readBigInt32()
  result = newString(bytes)
  pr.read(result)

proc modWrite(pw: var ParallelWriter, s: string) =
  var bytes = s.len
  pw.writeBigInt32 bytes
  pw.write s

type
  ModFileHeader* = object
    magic*: string
    version*: int32
    userdata*: string
    mapstart*: int

proc newModFileHeader*(ud = ""): ModFileHeader =
  result.magic = "XXXXQDPLazyDiskMapObjFileXXXX"
  result.version = 1
  result.userdata = ud
  result.mapstart = 0

proc modReadHeader*(pr: var ParallelReader): ModFileHeader =
  result.magic = pr.modReadString()
  result.version = pr.readBigInt32()
  result.userdata = pr.modReadString()
  discard pr.readBigInt64()
  result.mapstart = pr.readBigInt64().int
  #echo result

proc write*(pw: var ParallelWriter, h: ModFileHeader) =
  #echo h.userdata
  #echo h.mapstart
  let a = pw.active
  pw.setSingle()
  pw.modWrite h.magic
  pw.writeBigInt32 h.version
  pw.modWrite h.userdata
  pw.writeBigInt64 0
  pw.writeBigInt64 h.mapstart
  pw.setActive a


type
  ModFileMap* = OrderedTableRef[string,int]

proc newModFileMap*(): ModFileMap =
  result = newOrderedTable[string,int]()

proc modReadMap*(pr: var ParallelReader, mdstart: int): ModFileMap =
  result = newModFileMap()
  pr.seekSet(mdstart)
  pr.beginLocalChecksum()
  let num = pr.readBigInt32()
  #echo num
  for i in 0..<num:
    let k = pr.modReadString()
    discard pr.readBigInt64()
    let v = pr.readBigInt64().int
    #echo k.toHex()
    #echo i, ": ", v
    result[k] = v
  pr.endLocalChecksum()
  let cks = pr.readBigInt32().uint32
  #echo "cksum: ", cks, "  ", pr.crc32
  doAssert(cks == pr.crc32)

proc getPos*(m: ModFileMap, x: auto): int =
  let l = sizeof(x)
  var s = newString(l)
  copyMem(s[0].addr, unsafeAddr(x), l)
  #echo s.toHex
  m[s]

proc getPos*(m: ModFileMap, a,b: int): int =
  let aa = toBigEndian(a.int32)
  let bb = toBigEndian(b.int32)
  m.getPos((aa,bb))

proc packKey*(v: seq[int]): string =
  let n = v.len * sizeof(int32)
  result = newString(n)
  let a = cast[ptr UncheckedArray[int32]](addr result[0])
  for i in 0..<v.len:
    a[i] = toBigEndian(v[i].int32)

proc unpackKey*(k: string): seq[int] =
  let a = cast[ptr UncheckedArray[int32]](unsafeAddr k[0])
  let n = k.len div sizeof(int32)
  result.newSeq(n)
  for i in 0..<n:
    result[i] = fromBigEndian(a[i])

proc add*(m: var ModFileMap, s: string, pos: int) =
  m[s] = pos

proc write*(pw: var ParallelWriter, m: ModFileMap) =
  let num = m.len
  #echo num
  let a = pw.active
  pw.setSingle()
  pw.beginChecksum
  pw.writeBigInt32 num
  for k in m.keys:
    pw.modWrite k
    pw.writeBigInt64 0
    pw.writeBigInt64 m[k]
  pw.endChecksum
  pw.writeBigInt32 pw.crc32.int32
  pw.setActive a



type
  ModFileReader* = object
    hdr*: ModFileHeader
    map*: ModFileMap
    r*: ParallelReader

proc newModFileReader*(r: var ParallelReader): ModFileReader =
  result.hdr = r.modReadHeader()
  result.map = r.modReadMap(result.hdr.mapstart)
  result.r = r

proc newModFileReader*(fn: string): ModFileReader =
  var r = openRead(fn)
  r.newModFileReader()

proc close*(mr: var ModFileReader) =
  mr.r.close()


type
  ModFileWriter* = object
    hdr*: ModFileHeader
    map*: ModFileMap
    w*: ParallelWriter

proc newModFileWriter*(w: var ParallelWriter, ud = ""): ModFileWriter =
  result.hdr = newModFileHeader(ud)
  result.map = newModFileMap()
  w.write result.hdr
  result.w = w

proc newModFileWriter*(fn: string, ud = ""): ModFileWriter =
  var w = openCreate(fn)
  w.newModFileWriter(ud)

proc beginWrite*(mw: var ModFileWriter, key: string) =
  mw.map.add(key, mw.w.pos)

proc endWrite*(mw: var ModFileWriter) =
  mw.hdr.mapstart = max(mw.hdr.mapstart, mw.w.pos)

proc close*(mw: var ModFileWriter) =
  mw.w.seekSet(mw.hdr.mapstart)
  mw.w.write(mw.map)
  mw.w.seekSet(0)
  mw.w.write(mw.hdr)
  mw.w.close()

when isMainModule:
  import qex
  import ../comms/comms
  qexInit()
  var comm = getComm()
  echo "rank: ", comm.rank, "/", comm.size
  var lat = intSeqParam("lat", @[4,4,4,16])
  let sfn = "colorvec.mod"
  let fn = "test_file.mod"
  var nio = intParam("nio", sqrt(comm.size.float).int)
  var ioranks = @[0]
  for i in 1..<nio:
    ioranks.add( (i*(comm.size-1)) div (nio-1) )
  echo "ioranks: ", ioranks

  let lo1 = newLayout(lat, 1)
  var nt = lat[^1]
  var size = lat
  size[^1] = 1
  var offset = @[0,0,0,0]
  var wm = newSeq[WriteMap](nt)
  for t in 0..<nt:
    offset[^1] = t
    wm[t] = lo1.setupWrite(size, offset, ioranks)

  let bytes = 24 * size.prod
  echo "bytes: ", bytes
  var cv1 = lo1.ColorVectorS1()

  var mr = newModFileReader(sfn)
  let buf = cast[pointer](alloc(bytes))
  var mw = newModFileWriter(fn, mr.hdr.userdata)

  for k in mr.map.keys:
    #echo k.toHex()
    let ti = unpackKey(k)
    let t = ti[0]
    #echo "t: ", ti[0], "  i: ", ti[1]
    let p = mr.map[k]
    #echo p
    mr.r.seekSet(p)
    mr.r.beginChecksum()
    #mr.r.readSingle(buf, bytes)
    mr.r.read(cv1, wm[t])
    mr.r.endChecksum()
    let cks = mr.r.readBigInt32().uint32
    doAssert(cks == mr.r.crc32)
    #echo "cksum: ", cks, "  ", mr.r.crc32
    mw.w.seekSet(p)
    mw.beginWrite(k)
    mw.w.beginChecksum()
    #mw.w.writeSingle(buf, bytes)
    mw.w.write(cv1, wm[t])
    mw.w.endChecksum()
    doAssert(cks == mw.w.crc32)
    mw.w.setSingle()
    mw.w.writeBigInt32 cks
    mw.w.setActive true
    mw.endWrite()
    #echo mw.hdr.mapstart

  mr.close()
  mw.close()

  #mw.beginWrite("key1")
  #mw.w.write("test object 1")

  #var pr = openRead(fn)
  #var srchdr = pr.modReadHeader()
  #var srcmap = pr.modReadMap(srchdr.mapstart)
  #echo srchdr

  #commsFinalize()
  qexFinalize()
