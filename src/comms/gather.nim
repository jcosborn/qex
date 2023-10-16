import comms, rankScatter, base/profile, field
import algorithm

template `&&`(x: char): untyped = cast[pointer](unsafeAddr(x))
template `&&`(x: int32): untyped = cast[pointer](unsafeAddr(x))
template `&&`(x: seq): untyped = cast[pointer](unsafeAddr(x[0]))
template newSeqOfCap[T](x: var seq[T], n: int) =
  x = newSeqOfCap[T](n)
template newSeqUninitialized[T](x: var seq[T], n: int) =
  x = newSeqUninitialized[T](n)

proc adjust[T](i: int, x: openarray[T], b: int): int =
  result = i
  if i > 0:
    let n = x.len
    while result<n:
      if x[result-1] div b != x[result] div b:
        break
      inc result

proc splitThreads[T](x: openarray[T], b: int, nt,myt: int): tuple[a:int,b:int] =
  let n = x.len
  var i0 = (n*myt) div nt;
  var i1 = (n*(myt+1)) div nt
  i0 = adjust(i0, x, b)
  i1 = adjust(i1, x, b)
  result = (i0,i1)

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
  GatherMap* = ref object
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
  result.new
  let rank = c.commrank
  rl.sort do (x, y: RecvList) -> int:
    result = cmp(x.srank, y.srank)
    if result == 0:
      result = cmp(x.didx, y.didx)
      if result == 0:
        result = cmp(x.sidx, y.sidx)

  result.lidx.newSeq(0)
  result.ldest.newSeq(0)
  var sidx = newSeqOfCap[int32](rl.len)
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
  if nsend > 0:
    let sbuf = cast[ptr UncheckedArray[int32]](unsafeAddr r.buf[0])
    result.sidx.newSeqUninitialized(nsend)
    for i in 0..<nsend:
      result.sidx[i] = sbuf[i]

    result.smsginfo.newSeqOfCap(r.mem.len)
    var st = 0
    for i in 0..<r.mem.len:
      let ln = r.mem[i].bytes div sizeof(int32)
      result.smsginfo.add MsgInfo(rank: r.mem[i].rank.int32, start: st.int32,
                                  count: ln.int32)
      st += ln

  result.nsrc = nsend + result.lidx.len
  result.ndest = rl.len

proc makeGatherMap*(c: Comm, sl: var seq[SendList]): GatherMap =
  result.new
  let rank = c.commrank
  sl.sort do (x, y: SendList) -> int:
    result = cmp(x.drank, y.drank)
    if result == 0:
      result = cmp(x.didx, y.didx)
      if result == 0:
        result = cmp(x.sidx, y.sidx)
  #echoAll fmt"{rank}: {sl.len}"

  result.lidx.newSeq(0)
  result.ldest.newSeq(0)
  var didx = newSeqOfCap[int32](sl.len)
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
    result.rdest.newSeqUninitialized(nrecv)
    for i in 0..<nrecv:
      result.rdest[i] = rbuf[i]

  result.rmsginfo.newSeqOfCap(r.mem.len)
  var st = 0
  for i in 0..<r.mem.len:
    let ln = r.mem[i].bytes div sizeof(int32)
    result.rmsginfo.add MsgInfo(rank: r.mem[i].rank.int32, start: st.int32,
                                count: ln.int32)
    st += ln

  result.nsrc = sl.len
  result.ndest = nrecv + result.lidx.len
  #echoAll fmt"{rank}: {result.smsginfo}, {result.rmsginfo}"

proc reverse*(gm: GatherMap): GatherMap =
  result.new
  result.nsrc = gm.ndest
  result.ndest = gm.nsrc
  result.sidx = gm.rdest
  result.rdest = gm.sidx
  result.smsginfo = gm.rmsginfo
  result.rmsginfo = gm.smsginfo
  result.lidx = gm.ldest
  result.ldest = gm.lidx
  result.sendbuf = gm.recvbuf
  result.recvbuf = gm.sendbuf

# TODO: startRecvs, startSends, doLocal, wait
proc gather*(c: Comm; gm: GatherMap; data: auto) =
  tic("gather")
  # start recvs
  let elemsize = data.elemsize
  let rsize = elemsize*gm.rdest.len
  #echo "rsize: ", rsize
  var recvbuf = newSeqUninitialized[int8](rsize)
  for i in 0..<gm.rmsginfo.len:
    #echoAll "pushRecv: ", i, "  ", gm.rmsginfo[i].rank
    c.pushRecv(gm.rmsginfo[i].rank, &&recvbuf[elemsize*gm.rmsginfo[i].start],
               elemsize*gm.rmsginfo[i].count)
  toc("start recvs")
  c.barrier  # make sure recvs posted before sends start
  toc("barrier")

  # fill send buffers and send
  let ssize = elemsize*gm.sidx.len
  #echo "ssize: ", ssize
  var sendbuf = newSeqUninitialized[int8](ssize)
  #toc("alloc sendbuf")
  var bytes = 0
  for i in 0..<gm.smsginfo.len:
    bytes += gm.smsginfo[i].count * elemsize
    #tic()
    #echoAll "pushSend: ", i, " count ", gm.smsginfo[i].count
    threads:
      tfor j, 0..<gm.smsginfo[i].count:
        let k = gm.smsginfo[i].start + j
        data.copy(&&sendbuf[elemsize*k], gm.sidx[k])
    #toc("gather: fill")
    c.pushSend(gm.smsginfo[i].rank,
               &&sendbuf[elemsize*gm.smsginfo[i].start],
               elemsize*gm.smsginfo[i].count)
    #toc("gather: send")
  toc("fill and send", flops=bytes)

  # copy local sites
  threads:
    #tfor i, 0..<gm.ldest.len:
    # make sure threads don't share inner sites
    let range = splitThreads(gm.ldest, data.vlen, numThreads, threadNum)
    when true:
      for i in range[0] ..< range[1]:
        data.copy(gm.ldest[i], gm.lidx[i])
    else:
      var desti: array[16,int]
      var idx: array[16,int]
      var ni = 0
      var vs0 = gm.ldest[range[0]] div data.vlen
      for i in range[0] ..< range[1]:
        let vs = gm.ldest[i] div data.vlen
        let vi = gm.ldest[i] mod data.vlen
        if vs != vs0:
          data.copy(vs0, desti, idx, ni)
          ni = 0
          vs0 = vs
        desti[ni] = vi
        idx[ni] = gm.lidx[i]
        inc ni
      if ni>0: data.copy(vs0, desti, idx, ni)
  bytes = gm.ldest.len * elemsize
  toc("copy local", flops=bytes)

  # wait recvs and copy
  bytes = 0
  for i in 0..<gm.rmsginfo.len:
    bytes += gm.rmsginfo[i].count * elemsize
    #tic()
    #echo "waitRecv: ", i
    c.waitRecv(i, free=false)
    #toc("gather: wait recv")
    #for j in 0..<gm.rmsginfo[i].count:
    #  let k = gm.rmsginfo[i].start + j
    #  data.copy(gm.rdest[k], &&recvbuf[elemsize*k])
    threads:
      let i0 = gm.rmsginfo[i].start
      let i1 = gm.rmsginfo[i].start + gm.rmsginfo[i].count;
      #tfor k, i0 ..< i1:
      # make sure threads don't share inner sites
      let range = splitThreads(toOpenArray(gm.rdest,i0,i1-1),
                               data.vlen, numThreads, threadNum)
      for j in range[0] ..< range[1]:
        let k = i0 + j
        data.copy(gm.rdest[k], &&recvbuf[elemsize*k])
    #toc("gather: copy")
  toc("wait recvs and copy", flops=bytes)

  # free recvs and wait sends
  c.freeRecvs(gm.rmsginfo.len)
  toc("free recvs")
  c.waitSends(gm.smsginfo.len)
  toc("wait sends")

type GatherPointer = object
  src: ptr UncheckedArray[char]
  dest: ptr UncheckedArray[char]
  elemsize: int
template vlen(gd: GatherPointer): int = 1
template copy(gd: GatherPointer, d: pointer, s: SomeInteger) =
  copyMem(d, cast[pointer](addr gd.src[gd.elemsize*s]), gd.elemsize)
template copy(gd: GatherPointer, d: SomeInteger, s: pointer) =
  copyMem(cast[pointer](addr gd.dest[gd.elemsize*d]), s, gd.elemsize)
template copy(gd: GatherPointer, d: SomeInteger, s: SomeInteger) =
  copyMem(cast[pointer](addr gd.dest[gd.elemsize*d]),
          cast[pointer](addr gd.src[gd.elemsize*s]), gd.elemsize)
proc gather*(c: Comm; gm: GatherMap; elemsize: int; d,s: pointer) =
  let dest = cast[ptr UncheckedArray[char]](d)
  let src = cast[ptr UncheckedArray[char]](s)
  var gd = GatherPointer(src: src, dest: dest, elemsize: elemsize)
  c.gather(gm, gd)

type GatherField[S,D] = object
  src: S
  dest: D
  elemsize: int
  vlen: int
template copy(gd: GatherField, d: pointer, s: SomeInteger) {.used.} =
  mixin evalType
  type T = evalType(gd.src{s})
  let p = cast[ptr T](d)
  p[] := gd.src{s}
template copy(gd: GatherField, d: SomeInteger, s: pointer) {.used.} =
  mixin evalType
  type T = evalType(gd.dest{d})
  let p = cast[ptr T](s)
  gd.dest{d} := p[]
template copy(gd: GatherField, d: SomeInteger, s: SomeInteger) {.used.} =
  mixin evalType
  gd.dest{d} := gd.src{s}
proc gather*[D,S:Field](c: Comm; gm: GatherMap; d: D, s: S) =
  mixin evalType
  type T = evalType(d{0})
  var gd = GatherField[S,D](src:s,dest:d,elemsize:sizeof(T),vlen:max(D.V,S.V))
  c.gather(gm, gd)

proc gatherReversed*(c: Comm; gm: GatherMap; elemsize: int;
                     dest,src: pointer) =
  var r = GatherMap.new
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
  import strformat
  import qex
  qexInit()

  echo "rank: ", myRank, "/", nRanks
  let c = getComm()
  let nranks = c.commsize
  let rank = c.commrank
  echo &"commrank: {rank}/{nranks}"

  proc test1 =
    let n = 20
    var rl = newSeq[RecvList](0)
    for i in 0..<n:
      let d = i.int32
      let r = (i mod nranks).int32
      let s = (i xor 1).int32
      rl.add RecvList(didx: d, srank: r, sidx: s)

    let gm = c.makeGatherMap(rl)
    #echoAll gm

    var src1 = newSeq[int](n)
    var dst1 = newSeq[int](n)
    for i in 0..<n:
      dst1[i] = -1
      src1[i xor 1] = -1
      if (i mod nranks) == rank:
        src1[i xor 1] = i

    c.gather(gm, sizeof(int), &&dst1, &&src1)

    echo dst1
    for i in 0..<n:
      if dst1[i] != i:
        echoAll "Failed ", rank, " ", i, " : ", dst1[i]

  test1()

  proc test2 =
    defaultSetup()
    var seed = 987654321'u
    var rng = newRngField(lo, MRG32k3a, seed)
    var cv0 = lo.ColorVector()
    var cv1 = lo.ColorVector()
    var cv2 = lo.ColorVector()
    cv0.gaussian rng
    let nd = lo.nDim
    var x = newSeq[int32](nd)
    for d in 0..<lo.ndim:
      var sl = newSeq[SendList](0)
      for i in lo.sites:
        lo.coord(x, i)
        x[d] = int32 (x[d]+1) mod lo.physGeom[d]
        let ri = lo.rankIndex(x)
        sl.add SendList(sidx: int32 i, drank: int32 ri.rank, didx: int32 ri.index)
      let gm = c.makeGatherMap(sl)
      c.gather(gm, cv1, cv0)
      let sh = newShifter(cv0, d, -1)
      cv2 := sh ^* cv0
      cv2 -= cv1
      echo cv2.norm2

  test2()

  qexFinalize()
