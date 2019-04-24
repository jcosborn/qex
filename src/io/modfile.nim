import parallelIo
import endians, tables, strutils

proc modReadString(pr: var ParallelReader): string =
  var bytes = pr.readBigInt32()
  result = newString(bytes)
  pr.readAll(result)

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
  var dum = pr.readBigInt64()
  result.mapstart = pr.readBigInt64().int
  #echo result

proc write*(pw: var ParallelWriter, h: ModFileHeader) =
  #echo h.userdata
  #echo h.mapstart
  pw.setSingle()
  pw.modWrite h.magic
  pw.writeBigInt32 h.version
  pw.modWrite h.userdata
  pw.writeBigInt64 0
  pw.writeBigInt64 h.mapstart
  pw.setActive true


type
  ModFileMap* = OrderedTableRef[string,int]

proc newModFileMap*(): ModFileMap =
  result = newOrderedTable[string,int]()

proc modReadMap*(pr: var ParallelReader, mdstart: int): ModFileMap =
  result = newModFileMap()
  pr.seekSet(mdstart)
  let num = pr.readBigInt32()
  #echo num
  for i in 0..<num:
    let k = pr.modReadString()
    var dum = pr.readBigInt64()
    let v = pr.readBigInt64().int
    #echo k.toHex()
    #echo i, ": ", v
    result.add(k, v)

proc getPos*(m: ModFileMap, x: any): int =
  let l = sizeof(x)
  var s = newString(l)
  copyMem(s[0].addr, unsafeAddr(x), l)
  #echo s.toHex
  m[s]

proc getPos*(m: ModFileMap, a,b: int): int =
  let aa = toBigEndian(a.int32)
  let bb = toBigEndian(b.int32)
  m.getPos((aa,bb))

proc add*(m: var ModFileMap, s: string, pos: int) =
  m[s] = pos

proc write*(pw: var ParallelWriter, m: ModFileMap) =
  let num = m.len
  #echo num
  pw.beginChecksum
  pw.writeBigInt32 num
  for k in m.keys:
    pw.modWrite k
    pw.writeBigInt64 0
    pw.writeBigInt64 m[k]
  pw.endChecksum
  pw.writeBigInt32 pw.crc32.uint32



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

proc close*(mr: ModFileReader) =
  mr.r.close()


type
  ModFileWriter* = object
    hdr*: ModFileHeader
    map*: ModFileMap
    w*: ParallelWriter

proc newModFileWriter*(w: var ParallelWriter, ud = ""): ModFileWriter =
  result.hdr = newModFileHeader(ud)
  result.map = newModFileMap()
  result.w = w
  w.write result.hdr

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
  import ../comms/comms
  commsInit()
  echo "rank: ", myRank, "/", nRanks

  let sfn = "colorvec.mod"
  var mr = newModFileReader(sfn)
  let bytes = 1540
  let buf = cast[pointer](alloc(bytes))

  let fn = "test_file.mod"
  var mw = newModFileWriter(fn, mr.hdr.userdata)

  for k in mr.map.keys:
    let p = mr.map[k]
    #echo p
    mr.r.seekSet(p)
    mr.r.readSingle(buf, bytes)
    mw.w.seekSet(p)
    mw.beginWrite(k)
    mw.w.writeSingle(buf, bytes)
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

  commsFinalize()
