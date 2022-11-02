import math
import comms

template `&&`(x: char): untyped = cast[pointer](unsafeAddr(x))
#template `&&`(x: int): untyped = cast[pointer](unsafeAddr(x))

type
  RankScatterDescriptor*[T] = object
    rank*: int
    data*: T
  RankScatterSeq*[T] = seq[RankScatterDescriptor[T]]
  RankScatterDescriptor2*[T] = object
    srank*: int
    desc*: RankScatterDescriptor[T]
  RankScatterMem* = object
    rank*: int
    bytes*: int
    data*: pointer
  RankScatterMemSeq* = seq[RankScatterMem]

proc exchange[T](desc: var seq[RankScatterDescriptor2[T]], comm: auto,
                 rank,srank,rrank,nranks: int, bitwise=true) =
  template send(r,d) = comm.pushSend(r,d)
  template recv(r,d) = comm.pushRecv(r,d)
  template waitSend() = comm.waitSends
  template waitRecv() = comm.waitRecvs
  template waitAll() = comm.waitAll
  #let n2 = nranks div 2
  #template rankDiff(x,y: int): int = abs(((nranks+x-y+n2) mod nranks) - n2)
  template rankDiff(x,y: int): int = ((nranks+x-y) mod nranks)

  var nmsgs = 0
  recv(rrank, nmsgs)
  var smsg = newSeq[RankScatterDescriptor2[T]](0)
  var rmsg = newSeq[RankScatterDescriptor2[T]](0)
  var newdesc = newSeq[RankScatterDescriptor2[T]](0)
  for i in 0..<desc.len:
    let drank = desc[i].desc.rank
    if bitwise:
      if (drank xor srank) < (drank xor rank):
        #echo "bsend: ", drank, "  ", rank, "  ", srank
        smsg.add desc[i]
      else:
        #echo "bsave: ", drank, "  ", rank, "  ", srank
        newdesc.add desc[i]
    else:
      let rd = rankDiff(drank,rank)
      let rds = rankDiff(drank,srank)
      let d = rd xor rds
      #if rds <= rd:
      if rd!=0 and (rds==0 or (d and (d-1)) == 0):
        #echo "send: ", drank, "  ", rank, "  ", srank
        smsg.add desc[i]
      else:
        #echo "save: ", drank, "  ", rank, "  ", srank
        newdesc.add desc[i]
  send(srank, smsg.len)
  waitRecv()
  if nmsgs>0:
    rmsg.setLen(nmsgs)
    recv(rrank, rmsg)
  waitSend()
  if smsg.len>0:
    send(srank, smsg)
  waitAll()

  if rrank > rank:
    newdesc.add rmsg
    when defined(gcArc) or defined(gcOrc):
      desc = move newdesc
    else:
      shallowCopy desc, newdesc
  else:
    rmsg.add newdesc
    when defined(gcArc) or defined(gcOrc):
      desc = move rmsg
    else:
      shallowCopy desc, rmsg

proc scatter*[T](desc0: RankScatterSeq[T], comm: auto): RankScatterSeq[T] =
  let myrank = comm.commrank
  let nranks = comm.commsize
  #template `&`(x: var char): untyped = cast[pointer](addr x)

  var desc = newSeq[RankScatterDescriptor2[T]](desc0.len)
  for i in 0..<desc0.len:
    desc[i].srank = myrank
    desc[i].desc = desc0[i]

  var d = 1
  while (nranks mod (2*d)) == 0:
    let srank = myrank xor d
    let rrank = srank
    exchange(desc, comm, myrank, srank, rrank, nranks)
    d = 2*d

  #d = nextPowerOfTwo(nranks) div 2
  #while d>0 and (nranks mod (2*d)) != 0:
  while d<nranks:
    let srank = (myrank + d + nranks) mod nranks
    let rrank = (myrank - d + nranks) mod nranks
    exchange(desc, comm, myrank, srank, rrank, nranks, false)
    #d = d div 2
    d = 2*d

  result.newSeq(desc.len)
  for i in 0..<desc.len:
    result[i] = desc[i].desc
    result[i].rank = desc[i].srank

proc scatter*(mem: RankScatterMemSeq, c: Comm): tuple[buf:seq[char],mem:RankScatterMemSeq] =
  let rank = c.commrank

  # find recv ranks and sizes
  type Rsd = RankScatterDescriptor[int]
  var desc = newSeq[Rsd](0)
  for i in 0..<mem.len:
    desc.add Rsd(rank: mem[i].rank, data: mem[i].bytes)
  let b = desc.scatter(c)

  # start recvs
  var rbytes = 0
  for i in 0..<b.len:
    rbytes += b[i].data
  result.buf.newSeq(rbytes)
  result.mem.newSeq(b.len)
  var rpos = 0
  for i in 0..<b.len:
    if b[i].rank != rank:
      c.pushRecv(b[i].rank, &&result.buf[rpos], b[i].data)
    result.mem[i].rank = b[i].rank
    result.mem[i].bytes = b[i].data
    result.mem[i].data = &&result.buf[rpos]
    rpos += b[i].data

  # start sends
  for i in 0..<mem.len:
    if mem[i].rank != rank:
      c.pushSend(mem[i].rank, mem[i].data, mem[i].bytes)

  # copy local
  var j = 0
  for i in 0..<b.len:
    if b[i].rank == rank:
      while mem[j].rank != rank: inc j
      copyMem(result.mem[i].data, mem[j].data, mem[j].bytes)

  c.waitRecvs(b.len)
  c.waitSends(mem.len)

when isMainModule:
  import qex
  qexInit()
  echo "rank ", myRank, "/", nRanks

  var c = getComm()
  type Rsd = RankScatterDescriptor[int]
  var desc = newSeq[Rsd](0)
  let fac = 100
  for i in 0..<nRanks:
    desc.add Rsd(rank: i, data: myRank*fac + i)
  let r = desc.scatter(c)
  echoAll r
  for i in 0..<r.len:
    let d = r[i].data
    let mr = d mod fac
    let sr = d div fac
    if (myrank!=mr) or (r[i].rank!=sr):
      echoAll "error: ", myrank, " ", r[i]

  var mem = newSeq[RankScatterMem](0)
  for i in 0..<desc.len:
    mem.add RankScatterMem(rank: desc[i].rank, bytes: sizeof(desc[i].data), data: &&desc[i].data)
  let m = mem.scatter(c)
  for i in 0..<m.mem.len:
    if (m.mem[i].rank != r[i].rank) or (cast[ptr int](m.mem[i].data)[] != r[i].data):
      echoAll "error: ", myrank, " ", r[i]

  qexFinalize()
