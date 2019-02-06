import posix
import ../comms/comms
import strformat
import ../comms/gather
import layout
import endians

template `&`(x: char): untyped = cast[ptr UncheckedArray[char]](unsafeAddr(x))
template `&`[T](x: seq[T]): untyped =
  cast[ptr UncheckedArray[T]](unsafeAddr(x[0]))
template `&&`(x: char): untyped = cast[pointer](unsafeAddr(x))
template `&&`(x: seq): untyped = cast[pointer](unsafeAddr(x[0]))

proc posixOpenWrite*(fn: string): cint =
  let flags = O_WRONLY
  let mode = 0o666
  open(fn, flags, mode)

proc posixCreate*(fn: string, size=0): cint =
  let flags = O_CREAT or O_WRONLY
  let mode = 0o666
  let fd = open(fn, flags, mode)
  discard ftruncate(fd, size)
  fd


type
  ParallelReader* = object
    comm*: Comm
    pos*: int
    fd*: cint

proc openRead*(c: Comm, fn: string): ParallelReader =
  result.comm = c
  let flags = O_RDONLY
  result.fd = open(fn, flags)

proc openRead*(fn: string): ParallelReader =
  var c = getComm()
  c.openRead(fn)

proc close*(r: ParallelReader) =
  r.comm.barrier()
  discard close(r.fd)

proc seekCur*(r: var ParallelReader, offset: int) =
  r.pos += offset
  discard lseek(r.fd, offset, SEEK_CUR)

proc read*(r: var ParallelReader, buf: pointer, bytes: int) =
  r.pos += bytes
  discard read(r.fd, buf, bytes)

proc readAll*(r: var ParallelReader, buf: pointer, bytes: int) =
  r.read(buf, bytes)

proc readAll*(r: var ParallelReader, s: var SomeNumber) =
  r.readAll(addr(s), sizeof(s))

proc readAll*(r: var ParallelReader, s: var string) =
  r.readAll(addr(s[0]), s.len)

proc readSingle*(r: var ParallelReader, buf: pointer, bytes: int) =
  if r.comm.isMaster:
    r.read(buf, bytes)
  else:
    r.seekCur(bytes)

proc readSingle*(r: var ParallelReader, s: var string) =
  r.readSingle(addr(s[0]), s.len)

proc readBigInt32(pr: var ParallelReader): int32 =
  when system.cpuEndian == bigEndian:
    pr.readAll(result)
  else:
    var t: int32
    pr.readAll(t)
    template `&&`(x: int32): untyped = cast[pointer](unsafeAddr(x))
    swapEndian32(&&result, &&t)


type
  ParallelWriter* = object
    comm*: Comm
    pos*: int
    fd*: cint

proc openCreate*(c: Comm, fn: string, size=0): ParallelWriter =
  result.comm = c
  if c.isMaster:
    result.fd = posixCreate(fn, size)
    c.barrier
  else:
    c.barrier
    result.fd = posixOpenWrite(fn)
  c.barrier

proc openCreate*(fn: string, size=0): ParallelWriter =
  var c = getComm()
  c.openCreate(fn, size)

proc close*(w: ParallelWriter) =
  w.comm.barrier()
  discard close(w.fd)

proc seekCur*(w: var ParallelWriter, offset: int) =
  w.pos += offset
  discard lseek(w.fd, offset, SEEK_CUR)

proc write*(w: var ParallelWriter, buf: pointer, bytes: int) =
  w.pos += bytes
  discard write(w.fd, buf, bytes)

proc writeSingle*(w: var ParallelWriter, buf: pointer, bytes: int) =
  if w.comm.isMaster:
    w.write(buf, bytes)
  else:
    w.seekCur(bytes)

proc writeSingle*(w: var ParallelWriter, s: string) =
  w.writeSingle(unsafeAddr(s[0]), s.len)



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

#[
type
  IndexMapObj* = object
    srcidx*: int32
    destidx*: int32
proc createWriteMap*(wim: seq[IndexMapObj], c: Comm, nwriters=0): WriteMap =
  let rank = c.commrank
  let nranks = c.commsize
]#

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

  pr.comm.gatherReversed(wm.gm, elemsize, rp, dbuf)

  pr.seekCur(elemsize*wm.endoffset)

proc write*(pw: var ParallelWriter, wm: WriteMap,
            elemsize: int, sbuf: pointer) =
  let wsize = elemsize * wm.gm.ndest
  #echo "elemsize: ", elemsize, "  wsize: ", wsize
  var wp: pointer
  if wsize>0:
    var wbuf = newSeq[char](wsize)
    wp = &&wbuf

  pw.comm.gather(wm.gm, elemsize, sbuf, wp)

  # write data
  if wsize > 0:
    pw.seekCur(elemsize*wm.offset)
    pw.write(wp, wsize)
  pw.seekCur(elemsize*wm.endoffset)

import field

proc read*(pr: var ParallelReader, f: Field, wm: WriteMap) =
  let elemsize = sizeof(f[0])
  when f.V==1:
    pr.read(wm, elemsize, cast[pointer](unsafeAddr f[0]))
  else:
    echo "parallel read V: ", f.V, " not supported yet!"
    qexAbort(-1)

proc write*(pw: var ParallelWriter, f: Field, wm: WriteMap) =
  let elemsize = sizeof(f[0])
  when f.V==1:
    pw.write(wm, elemsize, cast[pointer](unsafeAddr f[0]))
  else:
    echo "parallel write V: ", f.V, " not supported yet!"
    qexAbort(-1)


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
