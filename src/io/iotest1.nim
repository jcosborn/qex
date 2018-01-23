import ../comms/mpi
import posix
import times
import os
import strutils

type TimerObj = tuple[s:string,d:float,n:int]
var timer: float
var deltas = newSeq[ptr TimerObj](0)
template tic =
  timer = epochTime()
template toc(st: string) =
  block:
    var first{.global.} = true
    var to{.global.}: TimerObj
    if first:
      first = false
      deltas.add to.addr
      to.s = st
      to.d = 0.0
      to.n = 0
    let t = epochTime()
    to.d += t - timer
    timer = t
    inc to.n

discard MPI_init()
let comm = MPI_COMM_WORLD
var rank, size: cint
discard MPI_Comm_rank(comm, rank.addr)
discard MPI_Comm_size(comm, size.addr)
var fn = "testout.dat"
if paramCount()>=1:
  fn = paramStr(1)
var nwriters = size.int
var nreaders = size.int
if paramCount()>=3:
  nwriters = parseInt(paramStr(3))
nreaders = nwriters
if paramCount()>=4:
  nreaders = parseInt(paramStr(4))
var maxBlock = 16*1024*1024*1024

proc testContig(n: int) =
  var wbuf = newSeq[cint](n)
  var rbuf = newSeq[cint](n)
  for i in 0..<n:
    wbuf[i] = cint(i+1)

  let wflags0 = O_CREAT or O_WRONLY
  let wflags1 = O_WRONLY
  let wmode = 0o666
  let wsize = wbuf.len*sizeof(type(wbuf[0]))
  var woffset = rank*wsize

  let rflags = O_RDONLY
  let rmode = 0o666
  let rsize = rbuf.len*sizeof(type(rbuf[0]))
  var roffset = rank*rsize

  tic()
  var wfd: cint
  if rank==0:
    wfd = open(fn, wflags0, wmode)
    let gsize = wsize*size
    discard ftruncate(wfd, gsize)
    discard MPI_barrier(comm)
    discard MPI_barrier(comm)
  else:
    discard MPI_barrier(comm)
    wfd = open(fn, wflags1, wmode)
    discard MPI_barrier(comm)
  toc("open write")
  #let woff2 = lseek(wfd, woffset, SEEK_SET)
  #discard MPI_barrier(comm)
  #toc("seek write")
  let wMaxElems = maxBlock div sizeof(type(wbuf[0]))
  let wb = (n+wMaxElems-1) div wMaxElems
  let wn = (size+nwriters-1) div nwriters
  let woffi1 = size div wn
  let woffi2 = rank div wn
  for b in 0..<wb:
    for i in 0..<wn:
      if i == rank mod wn:
        let we0 = b*wMaxElems
        let we1 = min((b+1)*wMaxElems,n)
        let wen = (we1-we0)*sizeof(type(wbuf[0]))
        woffset = wen * (((b)*wn+i)*woffi1+woffi2)
        let woff2 = lseek(wfd, woffset, SEEK_SET)
        let wsize2 = write(wfd, wbuf[we0].voidaddr, wen)
      discard MPI_barrier(comm)
  toc("write")
  discard close(wfd)
  discard MPI_barrier(comm)
  toc("write close")

  tic()
  let rfd = open(fn, rflags, rmode)
  discard MPI_barrier(comm)
  toc("open read")
  #let roff2 = lseek(rfd, roffset, SEEK_SET)
  #discard MPI_barrier(comm)
  #toc("seek read")
  let rMaxElems = maxBlock div sizeof(type(rbuf[0]))
  let rb = (n+rMaxElems-1) div rMaxElems
  let rn = (size+nreaders-1) div nreaders
  let roffi1 = size div wn
  let roffi2 = rank div wn
  for b in 0..<rb:
    for i in 0..<rn:
      if i == rank mod rn:
        let re0 = b*rMaxElems
        let re1 = min((b+1)*rMaxElems,n)
        let ren = (re1-re0)*sizeof(type(rbuf[0]))
        roffset = ren * (((b)*rn+i)*roffi1+roffi2)
        let roff2 = lseek(rfd, roffset, SEEK_SET)
        let rsize2 = read(rfd, rbuf[re0].voidaddr, ren)
      discard MPI_barrier(comm)
  toc("read")
  discard close(rfd)
  discard MPI_barrier(comm)
  toc("read close")

  var errors,errors2 = 0
  for i in 0..<n:
    if rbuf[i] != wbuf[i]:
      #echo rank, ": ", i, ": rbuf: ", rbuf[i], " wbuf: ", wbuf[i]
      inc errors
      #break
  discard MPI_Allreduce(errors.voidaddr, errors2.voidaddr,
                        1.cint, MPI_INT64_T, MPI_SUM, comm)
  if errors2>0:
    if rank==0:
      echo "errors: ", errors2 / size

  if rank==0: removeFile(fn)


var n = 1024
var nmax = 2*72*16*1024
if paramCount()>=2:
  nmax = parseInt(paramStr(2))*1024*1024

proc disp =
  if rank==0:
    echo "total ranks: ", size
    echo "using file: ", fn
    echo "nwriters: ", nwriters
    echo "nreaders: ", nreaders
    echo "maxBlock: ", maxBlock
    echo "sizeof(Off): ", sizeof(Off)
    {.emit:"""printf("sizeof(off_t): %i\n", sizeof(off_t));""".}
disp()

while n<=nmax:
  if rank==0:
    let b = formatSize(n*size*sizeof(cint))
    echo "total file size: ", b
  #testContig(n)
  #testContig(n)
  for i in 0..<deltas.len:
    deltas[i].d = 0.0
    deltas[i].n = 0
  testContig(n)
  testContig(n)
  for i in 0..<deltas.len:
    if rank==0:
      let sp = max(0, 11-deltas[i].s.len)
      let t = formatFloat(deltas[i].d/deltas[i].n.float, ffDecimal, 9)
      echo deltas[i].s, spaces(sp), ": ", align(t,15)
  n *= 2

discard MPI_Finalize()
