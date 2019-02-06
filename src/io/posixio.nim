import base
import posix

type
  PosixFile* = object
    fd*: cint

proc posixOpenRead*(fn: string): PosixFile =
  let flags = O_RDONLY
  result.fd = open(fn, flags)

proc posixOpenWrite*(fn: string): PosixFile =
  let flags = O_WRONLY
  result.fd = open(fn, flags)

proc posixCreate*(fn: string, size=0): PosixFile =
  let flags = O_CREAT or O_WRONLY
  let fd = open(fn, flags)
  if size>0:
    discard ftruncate(fd, size)


type
  ReaderMap* = object
    nsites*: int
    nsends*: int
    sendRanks*: seq[int32]
    nsendSites*: seq[int32]
    sendSiteOrder*: seq[int32]
    rbuf*: seq[char]
    sbuf*: seq[char]

# int aio_read(struct aiocb *aiocbp);
#let rsize2 = read(rfd, rbuf[re0].voidaddr, ren)

# start recv
# read block
# inplace reorder
# send



proc read(rm: ReaderMap, rfd: cint, elemsize: int) =
  #for i in 0..<rm.nrecvs:

  let rsize = elemsize * rm.nsites
  if rm.buf.len<rsize:
    rm.buf.setlen(rsize)

  let br = rfd.read(cast[pointer](rm.buf[0].addr), rsize)
  if br<rsize:
    echo fmt"warning: bytes read({br})<request({rsize})"

  for i in 0..<rm.nsites:
    copyMem(sbuf[elemsize*i].addr, rbuf[elemsize*rm.sendSitesOrder[i]].addr)

  for i in 0..<rm.nsends:
    startsend(



