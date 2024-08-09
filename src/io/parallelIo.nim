import posix
import ../comms/comms
#import strformat
import ../comms/gather
import layout
import endians, crc32

#template `&`(x: char): untyped = cast[ptr UncheckedArray[char]](unsafeAddr(x))
#template `&`[T](x: seq[T]): untyped =
#  cast[ptr UncheckedArray[T]](unsafeAddr(x[0]))
#template `&&`(x: char): untyped = cast[pointer](unsafeAddr(x))
template `&&`(x: seq): untyped = cast[pointer](unsafeAddr(x[0]))
template `&&`(x: SomeNumber): untyped = cast[pointer](unsafeAddr(x))

proc toBigEndian*(x: int32): int32 =
  when system.cpuEndian == bigEndian:
    result = x
  else:
    swapEndian32(&&result, &&x)
template fromBigEndian*(x: int32): int32 = toBigEndian(x)

proc toBigEndian*(x: float32): float32 =
  when system.cpuEndian == bigEndian:
    result = x
  else:
    swapEndian32(&&result, &&x)

proc swapEndian32*(p: pointer, bytes: int) =
  let s = cast[ptr UncheckedArray[int32]](p)
  let n = bytes div 4
  for i in 0..<n:
    var t = s[][i]
    swapEndian32(&&s[][i], &&t)

proc swapEndian64*(p: pointer, bytes: int) =
  let s = cast[ptr UncheckedArray[int64]](p)
  let n = bytes div 8
  for i in 0..<n:
    var t = s[][i]
    swapEndian64(&&s[][i], &&t)


proc posixOpenWrite*(fn: string): cint =
  let flags = O_WRONLY
  let mode = Mode 0o666
  open(fn, flags, mode)

proc posixCreate*(fn: string, size=0): cint =
  let flags = O_CREAT or O_WRONLY
  let mode = Mode 0o666
  let fd = open(fn, flags, mode)
  discard ftruncate(fd, size)
  fd


type
  ParallelReader* = object
    comm*: Comm
    pos*: int
    retval*: int
    swap*: int
    fd*: cint
    crc32*: Crc32
    active*: bool
    doChecksum*: bool

proc openRead*(c: Comm, fn: string): ParallelReader =
  let flags = O_RDONLY
  result.fd = open(fn, flags)
  if result.fd < 0:
    echo "error opening file: ", fn
  result.comm = c
  result.pos = 0
  result.retval = 0
  result.swap = 0
  result.active = true
  result.doChecksum = false

proc openRead*(fn: string): ParallelReader =
  var c = getComm()
  c.openRead(fn)

proc close*(r: var ParallelReader) =
  r.comm.barrier()
  r.retval = close(r.fd)

proc setActive*(r: var ParallelReader, a: bool) =
  r.active = a

proc setSingle*(r: var ParallelReader) =
  if r.comm.isMaster:
    r.setActive true
  else:
    r.setActive false

proc setSwap*(r: var ParallelReader, sw: int) =
  r.swap = sw

proc setBig32*(r: var ParallelReader) =
  when system.cpuEndian == bigEndian:
    r.swap = 0
  else:
    r.swap = 32

proc setBig64*(r: var ParallelReader) =
  when system.cpuEndian == bigEndian:
    r.swap = 0
  else:
    r.swap = 64

proc beginChecksum*(r: var ParallelReader) =
  r.doChecksum = true
  if r.comm.isMaster:
    r.crc32 = InitCrc32
  else:
    r.crc32 = Crc32(0)

proc endChecksum*(r: var ParallelReader) =
  r.doChecksum = false
  var t = int(r.crc32)
  r.comm.allReduceXor(t)
  r.crc32 = finishCrc32(Crc32(t))

proc beginLocalChecksum*(r: var ParallelReader) =
  r.doChecksum = true
  r.crc32 = InitCrc32

proc endLocalChecksum*(r: var ParallelReader) =
  r.doChecksum = false
  r.crc32 = finishCrc32(r.crc32)

proc seekSet*(r: var ParallelReader, offset: int) =
  r.retval = lseek(r.fd, offset, SEEK_SET).int
  if r.doChecksum:
    let offs = offset - r.pos
    doAssert(offs >= 0)
    r.crc32 = zeroPadCrc32(r.crc32, offs)
  r.pos = offset

proc seekCur*(r: var ParallelReader, offset: int) =
  r.retval = lseek(r.fd, offset, SEEK_CUR).int
  if r.doChecksum:
    doAssert(offset >= 0)
    r.crc32 = zeroPadCrc32(r.crc32, offset)
  r.pos += offset

proc readRaw*(r: var ParallelReader, buf: pointer, bytes: int) =
  r.retval = read(r.fd, buf, bytes)
  if r.doChecksum:
    doAssert(bytes >= 0)
    r.crc32 = updateCrc32(r.crc32, buf, bytes)
  r.pos += bytes

proc read*(r: var ParallelReader, buf: pointer, bytes: int) =
  if r.active:
    readRaw(r, buf, bytes)
    case r.swap
    of 32: swapEndian32(buf, bytes)
    of 64: swapEndian64(buf, bytes)
    else: discard
  else:
    r.seekCur(bytes)

proc read*(r: var ParallelReader, s: var SomeNumber) =
  r.read(addr(s), sizeof(s))

proc read*(r: var ParallelReader, s: var string) =
  r.read(addr(s[0]), s.len)

proc readAll*(r: var ParallelReader, buf: pointer, bytes: int) =
  let a = r.active
  r.setActive true
  r.read(buf, bytes)
  r.setActive a

proc readSingle*(r: var ParallelReader, buf: pointer, bytes: int) =
  let a = r.active
  r.setActive r.comm.isMaster
  r.read(buf, bytes)
  r.setActive a

proc readSingle*(r: var ParallelReader, s: var string) =
  r.readSingle(addr(s[0]), s.len)

proc readBigInt32*(pr: var ParallelReader): int32 =
  let s = pr.swap
  pr.setBig32()
  pr.read(&&result, sizeof(result))
  pr.setSwap s

proc readBigInt64*(pr: var ParallelReader): int64 =
  let s = pr.swap
  pr.setBig64()
  pr.read(&&result, sizeof(result))
  pr.setSwap s


type
  ParallelWriter* = object
    comm*: Comm
    pos*: int
    retval*: int
    swap*: int
    fd*: cint
    crc32*: Crc32
    active*: bool
    doChecksum*: bool

proc openCreate*(c: Comm, fn: string, size=0): ParallelWriter =
  if c.isMaster:
    result.fd = posixCreate(fn, size)
    c.barrier
  else:
    c.barrier
    result.fd = posixOpenWrite(fn)
  c.barrier
  result.comm = c
  result.pos = 0
  result.retval = 0
  result.swap = 0
  result.active = true
  result.doChecksum = false

proc openCreate*(fn: string, size=0): ParallelWriter =
  var c = getComm()
  c.openCreate(fn, size)

proc close*(w: var ParallelWriter) =
  w.comm.barrier()
  w.retval = close(w.fd)

proc setActive*(w: var ParallelWriter, a: bool) =
  w.active = a

proc setSingle*(w: var ParallelWriter) =
  if w.comm.isMaster:
    w.setActive true
  else:
    w.setActive false

proc setSwap*(w: var ParallelWriter, sw: int) =
  w.swap = sw

proc setBig32*(w: var ParallelWriter) =
  when system.cpuEndian == bigEndian:
    w.swap = 0
  else:
    w.swap = 32

proc setBig64*(w: var ParallelWriter) =
  when system.cpuEndian == bigEndian:
    w.swap = 0
  else:
    w.swap = 64

proc beginChecksum*(w: var ParallelWriter) =
  w.doChecksum = true
  if w.comm.isMaster:
    w.crc32 = InitCrc32
  else:
    w.crc32 = Crc32(0)

proc endChecksum*(w: var ParallelWriter) =
  w.doChecksum = false
  var t = int(w.crc32)
  w.comm.allReduceXor(t)
  w.crc32 = finishCrc32(Crc32(t))

proc seekSet*(w: var ParallelWriter, offset: int) =
  w.retval = lseek(w.fd, offset, SEEK_SET).int
  if w.doChecksum:
    let offs = offset - w.pos
    doAssert(offs >= 0)
    w.crc32 = zeroPadCrc32(w.crc32, offs)
  w.pos = offset

proc seekCur*(w: var ParallelWriter, offset: int) =
  w.retval = lseek(w.fd, offset, SEEK_CUR).int
  if w.doChecksum:
    doAssert(offset >= 0)
    w.crc32 = zeroPadCrc32(w.crc32, offset)
  w.pos += offset

proc writeRaw*(w: var ParallelWriter, buf: pointer, bytes: int) =
  #echo "writing at: ", w.pos, "  bytes: ", bytes
  w.retval = write(w.fd, buf, bytes)
  if w.doChecksum:
    doAssert(bytes >= 0)
    w.crc32 = updateCrc32(w.crc32, buf, bytes)
  w.pos += bytes

proc write*(w: var ParallelWriter, buf: pointer, bytes: int) =
  if w.active:
    case w.swap
    of 32: swapEndian32(buf, bytes)
    of 64: swapEndian64(buf, bytes)
    else: discard
    writeRaw(w, buf, bytes)
    case w.swap
    of 32: swapEndian32(buf, bytes)
    of 64: swapEndian64(buf, bytes)
    else: discard
  else:
    w.seekCur(bytes)

proc write*(w: var ParallelWriter, s: SomeNumber) =
  w.write(unsafeAddr(s), sizeof(s))

proc write*(w: var ParallelWriter, s: string) =
  w.write(unsafeAddr(s[0]), s.len)

proc writeSingle*(w: var ParallelWriter, buf: pointer, bytes: int) =
  let a = w.active
  w.setActive w.comm.isMaster
  w.write(buf, bytes)
  w.setActive a

proc writeSingle*(w: var ParallelWriter, s: string) =
  w.writeSingle(unsafeAddr(s[0]), s.len)

proc writeBigInt32*(pw: var ParallelWriter, x: SomeNumber) =
  let s = pw.swap
  pw.setBig32()
  pw.write(x.int32)
  pw.setSwap s

proc writeBigInt64*(pw: var ParallelWriter, x: SomeNumber) =
  let s = pw.swap
  pw.setBig64()
  pw.write(x.int64)
  pw.setSwap s


type
  WriteMap* = object
    offset*: int       # offset of bytes to write
    endoffset*: int    # offset to end of data after write
    #nwrite*: int       # number of elements to write
    #nsend*: int        # number of elememnt to send
    #sidx*: seq[int32]  # indices in src buffer to send (size nsend)
    #smsginfo*: seq[MsgInfo]  # send ranks, position and length in send buffer
    #lidx*: seq[int32]  # indices in src buffer of local elements
    #nrecv*: int        # number of elements to receive
    #rmsginfo*: seq[MsgInfo]  # recv ranks, position and length in recv buffer
    #rdest*: seq[int32] # indices in write buffer of recv buf elements (nrecv)
    #ldest*: seq[int32] # indices in write buffer of src buffer elements
    gm*: GatherMap

proc hyperIndex[T,U,V](x: seq[T], o: seq[U], s: seq[V]): int =
  for i in countdown(x.len-1,0):
    let k = x[i] - o[i]
    if k<0 or k>=s[i]:
      return -1
    result = result*s[i] + k

iterator localIndicesInHyper(lo: Layout, size: seq[int], offset: seq[int]):
    tuple[lidx:int,didx:int] =
  var x = newSeq[int32](lo.nDim)
  for i in 0..<lo.nSites:
    lo.coord(x, i)
    let k = hyperIndex(x, offset, size)
    if k>=0:
      yield (i,k)

proc setupWrite*(lo: Layout, size: seq[int], offset: seq[int],
                 ioranks: seq[int]): WriteMap =
  ## setup hypercubic subset write
  let nsites = size.product
  let nwriters = ioranks.len
  var c = getComm()
  var rank = c.commrank()

  # find my sites to write
  var sl = newSeq[SendList](0)
  for i,k in lo.localIndicesInHyper(size, offset):
    # i: local, k: file
    #echo "i: ", i, "  k: ", k
    let ri = (k * nwriters) div nsites
    let d = k - ((ri*nsites+nwriters-1) div nwriters)
    let r = ioranks[ri]
    sl.add SendList(sidx: i.int32, drank: r.int32, didx: d.int32)
  result.gm = c.makeGatherMap(sl)
  let ri = ioranks.find(rank)
  if ri<0:
    result.offset = 0
    result.endoffset = nsites
  else:
    let o = (ri*nsites+nwriters-1) div nwriters
    let o2 = ((ri+1)*nsites+nwriters-1) div nwriters
    result.offset = o
    result.endoffset = nsites - o2

proc read*(pr: var ParallelReader, wm: WriteMap,
           elemsize: int, dbuf: pointer) =
  let rsize = elemsize * wm.gm.ndest
  #echo "elemsize: ", elemsize, "  rsize: ", rsize
  var rp: pointer
  # read data
  if rsize>0:
    var rbuf = newSeq[char](rsize)
    rp = &&rbuf
    pr.seekCur(elemsize*wm.offset)
    pr.read(rp, rsize)

  pr.comm.gatherReversed(wm.gm, elemsize, dbuf, rp)

  pr.seekCur(elemsize*wm.endoffset)

proc write*(pw: var ParallelWriter, wm: WriteMap,
            elemsize: int, sbuf: pointer) =
  let wsize = elemsize * wm.gm.ndest
  #echo "elemsize: ", elemsize, "  wsize: ", wsize
  var wp: pointer
  if wsize>0:
    var wbuf = newSeq[char](wsize)
    wp = &&wbuf

  pw.comm.gather(wm.gm, elemsize, wp, sbuf)

  # write data
  if wsize > 0:
    pw.seekCur(elemsize*wm.offset)
    pw.write(wp, wsize)
  pw.seekCur(elemsize*wm.endoffset)

import field

#proc getAddr[T](x: var T): ptr T = addr x

proc read*(pr: var ParallelReader, f: Field, wm: WriteMap) =
  when f.V==1:
    let elemsize = sizeof(f[0])
    pr.read(wm, elemsize, cast[pointer](addr f[0]))
    #pr.read(wm, elemsize, cast[pointer](getAddr f[0]))
  else:
    echo "parallel read V: ", f.V, " not supported yet!"
    qexAbort(-1)

proc write*(pw: var ParallelWriter, f: Field, wm: WriteMap) =
  when f.V==1:
    let elemsize = sizeof(f[0])
    pw.write(wm, elemsize, cast[pointer](addr f[0]))
    #pw.write(wm, elemsize, cast[pointer](getAddr f[0]))
  else:
    echo "parallel write V: ", f.V, " not supported yet!"
    qexAbort(-1)

#[
type
  IndexMapObj* = object
    srcidx*: int32
    destidx*: int32
proc createWriteMap*(wim: seq[IndexMapObj], c: Comm, nwriters=0): WriteMap =
  let rank = c.commrank
  let nranks = c.commsize

type
  ReadMap* = object
    offset*: int       # offset of bytes to read
    endoffset*: int    # offset to end of data after read
    nread*: int        # number of elements to read
    nsend*: int        # number of elememnt to send
    sidx*: seq[int32]  # indices in read buffer to send (size nsend)
    smsginfo*: seq[MsgInfo]  # send ranks, position and length in send buffer
    lidx*: seq[int32]  # indices in read buffer of local elements
    nrecv*: int        # number of elements to receive
    rmsginfo*: seq[MsgInfo]  # recv ranks, position and length in recv buffer
    rdest*: seq[int32] # indices in dest array of recv buffer elements (nrecv)
    ldest*: seq[int32] # indices in dest array of read buffer elements

proc read*(pr: var ParallelReader, rm: ReadMap, c: Comm,
           elemsize: int, buf: ptr UncheckedArray[char]) =
  # start recv messages
  let recvsize = elemsize * rm.nrecv
  var recvbuf = newSeq[char](recvsize)
  for i in 0..<rm.rmsginfo.len:
    c.pushRecv(rm.rmsginfo[i].rank, &&recvbuf[elemsize*rm.rmsginfo[i].start],
               elemsize*rm.rmsginfo[i].count)

  # read data
  let readsize = elemsize * rm.nread
  var readbuf = newSeq[char](readsize)
  if readsize > 0:
    pr.seekCur(rm.offset)
    pr.read(&&readbuf, readsize)
  pr.seekCur(rm.endoffset)

  # prepare send buffer and send
  let sendsize = elemsize * rm.nsend
  var sendbuf = newSeq[char](sendsize)
  var sendi = 0
  var nextsend = rm.smsginfo[0].count - 1
  for i in 0..<rm.sidx.len:
    copyMem(&&sendbuf[elemsize*i], &&readbuf[elemsize*rm.sidx[i]], elemsize)
    if i == nextsend:
      c.pushSend(rm.smsginfo[sendi].rank,
                 &&sendbuf[elemsize*rm.smsginfo[sendi].start],
                 elemsize*rm.smsginfo[sendi].count)
      nextsend = rm.smsginfo[sendi].start + rm.smsginfo[sendi].count - 1
      sendi += 1

  # copy local data
  for i in 0..<rm.ldest.len:
    copymem(addr buf[elemsize*rm.ldest[i]],
            addr readbuf[elemsize*rm.lidx[i]], elemsize)

  # wait recv and copy
  var recvi = 0
  var nextrecv = 0
  for i in 0..<rm.nrecv:
    if i == nextrecv:
      c.waitRecv(recvi, free=false)
      nextrecv = rm.rmsginfo[recvi].start
      recvi += 1
    copymem(addr buf[elemsize*rm.rdest[i]],
            addr recvbuf[elemsize*i], elemsize)
  c.freeRecvs(rm.rmsginfo.len)

  # wait sends
  c.waitSends(rm.smsginfo.len)
]#

when isMainModule:
  proc dowrite(fn: string; x: seq; srcoffset,srcstride: int;
               c: Comm; nwriters: int) =
    let nranks = c.commsize
    let rank = c.commrank
    let nsites = (x.len-srcoffset+srcstride-1) div srcstride
    var wim = newSeq[IndexMapObj](0)
    for i in 0..<nsites:
      let s = i*srcstride + srcoffset
      let d = rank*nsites + i
      #let d = i*nranks + rank
      wim.add IndexMapObj(srcidx: s.int32, destidx: d.int32)
    let wm = createWriteMap(wim, c)
    let wr = c.openCreate(fn)
    wr.write(wm, c, sizeof(x[0]), &&x)
    wr.write(wm, c, sizeof(x[0]), &&x)

  commsInit()
  echo "rank: ", myRank, "/", nRanks

  var fn = "testrw.out"
  var c = getComm()
  let nranks = c.commsize
  let rank = c.commrank
  let nsites = 2
  let srcstride = 1
  let srcoffset = 0
  let nwriters = nranks
  var x = newSeq[int](nsites)
  echo fmt"commrank: {rank}/{nranks}"
  dowrite(fn, x, srcoffset, srcstride, c, nwriters)

  commsFinalize()
