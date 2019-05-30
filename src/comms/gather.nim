import comms
import rankScatter
import strformat
import algorithm
import base/profile

template `&&`(x: char): untyped = cast[pointer](unsafeAddr(x))
template `&&`(x: int32): untyped = cast[pointer](unsafeAddr(x))
template `&&`(x: seq): untyped = cast[pointer](unsafeAddr(x[0]))

type
  RecvList* = object
    didx*: int32  # destination index on this rank
    srank*: int32 # source rank
    sidx*: int32  # source index
  SendList* = object
    sidx*: int32  # source index on this rank
    drank*: int32 # dest rank
    didx*: int32  # dest index
  MsgInfo* = object
    rank*: int32
    start*: int32
    count*: int32
  GatherMap* = object
    nsrc*: int
    ndest*: int
    sidx*: seq[int32]  # indices in src buffer to send (size nsend)
    smsginfo*: seq[MsgInfo]  # send ranks, position and length in send buffer
    lidx*: seq[int32]  # indices in src buffer of local elements
    rmsginfo*: seq[MsgInfo]  # recv ranks, position and length in recv buffer
    rdest*: seq[int32] # indices in write buffer of recv buf elements (nrecv)
    ldest*: seq[int32] # indices in write buffer of src buffer elements
    sendbuf*: seq[char]
    recvbuf*: seq[char]

proc makeGatherMap*(c: Comm, rl: var seq[RecvList]): GatherMap =
  let rank = c.commrank
  rl.sort do (x, y: RecvList) -> int:
    result = cmp(x.srank, y.srank)
    if result == 0:
      result = cmp(x.didx, y.didx)

  result.lidx.newSeq(0)
  result.ldest.newSeq(0)
  var sidx = newSeq[int32](0)
  result.rdest.newSeq(0)
  result.rmsginfo.newSeq(0)
  var mem = newSeq[RankScatterMem](0)
  var lastrank,nsites: int32
  var i0,rst,i = 0
  while true:
    if i>=rl.len or rl[i].srank != lastrank:
      if nsites>0 and lastrank!=rank:
        mem.add RankScatterMem(rank: lastrank, bytes: nsites*sizeof(int32),
                               data: &&sidx[i0])
        result.rmsginfo.add MsgInfo(rank: lastrank, start: rst.int32,
                                    count: nsites.int32)
        rst += nsites
      if i>=rl.len:
        break
      lastrank = rl[i].srank
      i0 = sidx.len
      nsites = 0
    if rl[i].srank == rank:
      result.lidx.add rl[i].sidx
      result.ldest.add rl[i].didx
    else:
      sidx.add rl[i].sidx
      result.rdest.add rl[i].didx
    inc nsites
    inc i

  let r = mem.scatter(c)
  #echoAll r

  let nsend = r.buf.len div sizeof(int32)
  let sbuf = cast[ptr UncheckedArray[int32]](unsafeAddr r.buf[0])
  result.sidx.newSeq(nsend)
  for i in 0..<nsend:
    result.sidx[i] = sbuf[i]

  result.smsginfo.newSeq(r.mem.len)
  var st = 0
  for i in 0..<r.mem.len:
    let ln = r.mem[i].bytes div sizeof(int32)
    result.smsginfo[i] = MsgInfo(rank: r.mem[i].rank.int32, start: st.int32,
                                 count: ln.int32)
    st += ln

  result.nsrc = nsend + result.lidx.len
  result.ndest = rl.len

proc makeGatherMap*(c: Comm, sl: var seq[SendList]): GatherMap =
  let rank = c.commrank
  sl.sort do (x, y: SendList) -> int:
    result = cmp(x.drank, y.drank)
    if result == 0:
      result = cmp(x.didx, y.didx)
  #echoAll fmt"{rank}: {sl.len}"

  result.lidx.newSeq(0)
  result.ldest.newSeq(0)
  var didx = newSeq[int32](0)
  result.sidx.newSeq(0)
  result.smsginfo.newSeq(0)
  var mem = newSeq[RankScatterMem](0)
  var lastrank,nsites: int32
  var i0,sst,i = 0
  while true:
    if i>=sl.len or sl[i].drank != lastrank:
      if nsites>0 and lastrank!=rank:
        mem.add RankScatterMem(rank: lastrank, bytes: nsites*sizeof(int32),
                               data: &&didx[i0])
        result.smsginfo.add MsgInfo(rank: lastrank, start: sst.int32,
                                    count: nsites.int32)
        sst += nsites
      if i>=sl.len:
        break
      lastrank = sl[i].drank
      i0 = didx.len
      nsites = 0
    if sl[i].drank == rank:
      result.lidx.add sl[i].sidx
      result.ldest.add sl[i].didx
    else:
      didx.add sl[i].didx
      result.sidx.add sl[i].sidx
    inc nsites
    inc i

  let r = mem.scatter(c)
  #echoAll r

  let nrecv = r.buf.len div sizeof(int32)
  if nrecv>0:
    let rbuf = cast[ptr UncheckedArray[int32]](unsafeAddr r.buf[0])
    result.rdest.newSeq(nrecv)
    for i in 0..<nrecv:
      result.rdest[i] = rbuf[i]

  result.rmsginfo.newSeq(r.mem.len)
  var st = 0
  for i in 0..<r.mem.len:
    let ln = r.mem[i].bytes div sizeof(int32)
    result.rmsginfo[i] = MsgInfo(rank: r.mem[i].rank.int32, start: st.int32,
                                 count: ln.int32)
    st += ln

  result.nsrc = sl.len
  result.ndest = nrecv + result.lidx.len
  #echoAll fmt"{rank}: {result.smsginfo}, {result.rmsginfo}"

proc gather*(c: Comm; gm: GatherMap; elemsize: int; dest,src: pointer) =
  tic()
  let srcbuf = cast[ptr UncheckedArray[char]](src)
  let destbuf = cast[ptr UncheckedArray[char]](dest)

  # start recvs
  let rsize = elemsize*gm.rdest.len
  #echo "rsize: ", rsize
  var recvbuf = newSeqUninitialized[int8](rsize)
  for i in 0..<gm.rmsginfo.len:
    #echoAll "pushRecv: ", i, "  ", gm.rmsginfo[i].rank
    c.pushRecv(gm.rmsginfo[i].rank, &&recvbuf[elemsize*gm.rmsginfo[i].start],
               elemsize*gm.rmsginfo[i].count)
  toc("gather: start recvs")

  # fill send buffers and send
  let ssize = elemsize*gm.sidx.len
  #echo "ssize: ", ssize
  var sendbuf = newSeqUninitialized[int8](ssize)
  toc("gather: alloc sendbuf")
  for i in 0..<gm.smsginfo.len:
    #tic()
    #echoAll "pushSend: ", i
    threads:
      tfor j, 0..<gm.smsginfo[i].count:
        let k = gm.smsginfo[i].start + j
        copyMem(&&sendbuf[elemsize*k], &&srcbuf[elemsize*gm.sidx[k]], elemsize)
    #toc("gather: fill")
    c.pushSend(gm.smsginfo[i].rank,
               &&sendbuf[elemsize*gm.smsginfo[i].start],
               elemsize*gm.smsginfo[i].count)
    #toc("gather: send")
  toc("gather: fill and send")

  # copy local sites
  threads:
    tfor i, 0..<gm.ldest.len:
      copyMem(&&destbuf[elemsize*gm.ldest[i]], &&srcbuf[elemsize*gm.lidx[i]], elemsize)
  toc("gather: copy local")

  # wait recvs and copy
  for i in 0..<gm.rmsginfo.len:
    #tic()
    #echo "waitRecv: ", i
    c.waitRecv(i, free=false)
    #toc("gather: wait recv")
    threads:
      tfor j, 0..<gm.rmsginfo[i].count:
        let k = gm.rmsginfo[i].start + j
        copyMem(&&destbuf[elemsize*gm.rdest[k]], &&recvbuf[elemsize*k], elemsize)
    #toc("gather: copy")
  toc("gather: wait recvs and copy")

  # free recvs and wait sends
  c.freeRecvs(gm.rmsginfo.len)
  toc("gather: free recvs")
  c.waitSends(gm.smsginfo.len)
  toc("gather: wait sends")

proc gatherReversed*(c: Comm; gm: GatherMap; elemsize: int;
                     dest,src: pointer) =
  var r: GatherMap
  r.nsrc = gm.ndest
  r.ndest = gm.nsrc
  r.sidx = gm.rdest
  r.rdest = gm.sidx
  r.smsginfo = gm.rmsginfo
  r.rmsginfo = gm.smsginfo
  r.lidx = gm.ldest
  r.ldest = gm.lidx
  r.sendbuf = gm.recvbuf
  r.recvbuf = gm.sendbuf
  c.gather(r, elemsize, dest, src)

when isMainModule:
  commsInit()
  echo "rank: ", myRank, "/", nRanks
  let c = getComm()
  let nranks = c.commsize
  let rank = c.commrank
  echo &"commrank: {rank}/{nranks}"

  var rl = newSeq[RecvList](0)
  for i in 0..9:
    let d = i.int32
    let r = (i mod nranks).int32
    let s = (i xor 1).int32
    rl.add RecvList(didx: d, srank: r, sidx: s)

  let gm = c.makeGatherMap(rl)
  #echoAll gm

  var src1 = newSeq[int](10)
  var dst1 = newSeq[int](10)
  for i in 0..9:
    src1[i] = i xor 1
    dst1[i] = -1

  c.gather(gm, sizeof(int), &&dst1, &&src1)

  echo dst1

  commsFinalize()
