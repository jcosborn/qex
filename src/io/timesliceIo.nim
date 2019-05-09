import parallelIo
import modfile
import qex
import ../comms/comms

type TimesliceIo* = object
    wm*: seq[WriteMap]
    ioRanks*: seq[int]

proc getIoRanks*(ranks: int, nio0 = -1): seq[int] =
  var nio = nio0
  if nio<0: nio = intParam("nio", sqrt(ranks.float).int)
  var dior = @[0]
  for i in 1..<nio:
    dior.add( (i*(ranks-1)) div (nio-1) )
  result = intSeqParam("ior", dior)

proc newTimesliceIo*(lo: Layout, ioRanks: seq[int]): TimesliceIo =
  var size = lo.physGeom
  var nt = size[^1]
  size[^1] = 1
  var offset = @[0,0,0,0]
  result.wm = newSeq[WriteMap](nt)
  for t in 0..<nt:
    offset[^1] = t
    result.wm[t] = lo.setupWrite(size, offset, ioRanks)
  result.ioRanks = ioRanks

proc newTimesliceIo*(lo: Layout): TimesliceIo =
  var ioRanks = getIoRanks(lo.nranks)
  newTimesliceIo(lo, ioRanks)

proc read*(ti: TimesliceIo, pr: var ParallelReader, f: Field, t: int) =
  pr.beginChecksum()
  pr.read(f, ti.wm[t])
  pr.endChecksum()
  let cks = pr.readBigInt32().uint32
  if cks != pr.crc32:
    echo "timesliceIo read checksum error got: ", pr.crc32, " wanted: ", cks
    qexAbort(-1)

proc write*(ti: TimesliceIo, pw: var ParallelWriter, f: Field, t: int) =
  pw.beginChecksum()
  pw.write(f, ti.wm[t])
  pw.endChecksum()
  let a = pw.active
  pw.setSingle()
  pw.writeBigInt32 pw.crc32
  pw.setActive a

when isMainModule:
  qexInit()
  var comm = getComm()
  echo "rank: ", comm.rank, "/", comm.size
  var lat = intSeqParam("lat", @[4,4,4,16])
  echo "ioranks: ", getIoRanks(comm.size)

  let fn = "test_file.dat"
  let lo1 = newLayout(lat, 1)
  let tsio = newTimesliceIo(lo1)
  var nt = lat[^1]

  var cv1 = lo1.ColorVectorS1()
  var cv2 = lo1.ColorVectorS1()
  var r = newRngField(lo1, RngMilc6)
  cv1.gaussian(r)
  echo cv1.norm2

  var pw = openCreate(fn)
  for t in 0..<nt:
    tsio.write(pw, cv1, t)
  pw.close()

  var pr = openRead(fn)
  for t in 0..<nt:
    tsio.read(pr, cv2, t)
  pr.close()

  echo cv2.norm2
  cv2 -= cv1
  echo cv2.norm2

  qexFinalize()
