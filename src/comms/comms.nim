import commsQmp
export commsQmp
#import base/profile

import qmp

template `&&`(x: typed): untyped = cast[pointer](unsafeAddr x)

# commBuffer
# declareMessageGroup
# MessageGroup: startRecvs, startSends, start, wait
#   getRecvBuf(nrecv), getSendBuf(nsend)

type
  Comm* = ref object of RootObj

  CommQmp* = ref object of Comm
    comm: QMP_comm_t
    smem: seq[QMP_msgmem_t]
    smsg: seq[QMP_msghandle_t]
    rmem: seq[QMP_msgmem_t]
    rmsg: seq[QMP_msghandle_t]

method commrank*(c: Comm): int {.base.} = discard
method commsize*(c: Comm): int {.base.} = discard
method isMaster*(c: Comm): bool {.base.} = discard

template rank*(c: Comm): int = commrank(c)
template size*(c: Comm): int = commsize(c)

method barrier*(c: Comm) {.base.} = discard
method allReduce*(c: Comm, x: var int) {.base.} = discard
method allReduce*(c: Comm, x: ptr float32, n: int) {.base.} = discard
method allReduce*(c: Comm, x: var UncheckedArray[float32], n: int) {.base.} = discard
method allReduceXor*(c: Comm, x: var int) {.base.} = discard

method nsends*(c: Comm): int {.base.} = discard
method nrecvs*(c: Comm): int {.base.} = discard
method pushSend*(c: Comm, rank: int, p: pointer, bytes: int) {.base.} =
  discard
method pushRecv*(c: Comm, rank: int, p: pointer, bytes: int) {.base.} =
  discard
method waitSends*(c: Comm; i0,i1: int) {.base.} = discard
method freeRecvs*(c: Comm; i0,i1: int) {.base.} = discard
method waitRecvs*(c: Comm; i0,i1: int; free=true) {.base.} = discard

template pushSend*(c: Comm, rank: int, xx: SomeNumber) =
  var x = xx
  pushSend(c, rank, &&x, sizeof(x))

template pushSend*(c: Comm, rank: int, x: seq) =
  pushSend(c, rank, &&x[0], x.len*sizeof(x[0]))

template waitSends*(c: Comm) =
  c.waitSends(0, c.nsends-1)

template waitSends*(c: Comm, k: int) =
  let n = c.nsends
  let i0 = max(0, n-k)
  waitSends(c, i0, n-1)

template pushRecv*(c: Comm, rank: int, x: SomeNumber) =
  pushRecv(c, rank, &&x, sizeof(x))

template pushRecv*(c: Comm, rank: int, x: seq) =
  pushRecv(c, rank, &&x[0], x.len*sizeof(x[0]))

template freeRecvs*(c: Comm; k: int) =
  let n = c.nrecvs
  let i0 = max(0, n-k)
  freeRecvs(c, i0, n-1)

template waitRecv*(c: Comm, free=true) =
  let n = c.nrecvs - 1
  waitRecvs(c, n, n, free)

template waitRecv*(c: Comm, i: int, free=true) =
  waitRecvs(c, i, i, free)

template waitRecvs*(c: Comm, free=true) =
  waitRecvs(c, 0, c.nrecvs-1, free)

template waitRecvs*(c: Comm, k: int) =
  let n = c.nrecvs
  let i0 = max(0, n-k)
  waitRecvs(c, i0, n-1)

template waitAll*(c: Comm) =
  c.waitSends
  c.waitRecvs

#var commList: seq[Comm]
#proc push(c: Comm) = commList.add c
#proc pop(c: Comm) = commList.pop
#proc getComm(): Comm = commList[^1]
#proc getInitialComm(): Comm = commList[0]
#proc myRank(): int = myRank(getComm())
#proc myInitialRank(): int = myRank(getInitialComm())
#method myRank(c: Comm): int {.base.}
#method myRank(c: CommQmp): int = ...
#method declareSend
#method declareMultiple


proc getComm*(): CommQmp =
  result.new
  result.comm = QMP_comm_get_default()
  result.smem.newSeq(0)
  result.smsg.newSeq(0)
  result.rmem.newSeq(0)
  result.rmsg.newSeq(0)

method commrank*(c: CommQmp): int =
  QMP_comm_get_node_number(c.comm).int

method commsize*(c: CommQmp): int =
  QMP_comm_get_number_of_nodes(c.comm).int

method isMaster*(c: CommQmp): bool =
  c.commrank == 0

method barrier*(c: CommQmp) =
  QMP_comm_barrier(c.comm)

method allReduce*(c: CommQmp, x: var SomeNumber) =
  var t = float(x)
  QMP_comm_sum_double(c.comm, addr t)
  x = (type x)(t)

method allReduce*(c: CommQmp, x: ptr float32, n: int) =
  QMP_comm_sum_float_array(c.comm, x, n.cint)

method allReduce*(c: CommQmp, x: var UncheckedArray[float32], n: int) =
  QMP_comm_sum_float_array(c.comm, addr x[0], n.cint)

method allReduceXor*(c: CommQmp, x: var int) =
  var t = cast[ptr culong](addr x)
  QMP_comm_xor_ulong(c.comm, t)

method nsends*(c: CommQmp): int = c.smsg.len
method nrecvs*(c: CommQmp): int = c.rmsg.len

method pushSend*(c: CommQmp, rank: int, p: pointer, bytes: int) =
  let m = QMP_declare_msgmem(p, bytes)
  let h = QMP_declare_send_to(m, rank.cint, 0.cint)
  let stat = QMP_start(h)
  c.smem.add m
  c.smsg.add h

method pushRecv*(c: CommQmp, rank: int, p: pointer, bytes: int) =
  let m = QMP_declare_msgmem(p, bytes)
  let h = QMP_declare_receive_from(m, rank.cint, 0.cint)
  let stat = QMP_start(h)
  c.rmem.add m
  c.rmsg.add h

method waitSends*(c: CommQmp; i0,i1: int) =
  #tic()
  for i in i0..i1:
    #tic()
    discard QMP_wait(c.smsg[i])
    #toc("waitSends: wait")
    QMP_free_msghandle(c.smsg[i])
    #toc("waitSends: free msghandle")
    QMP_free_msgmem(c.smem[i])
    #toc("waitSends: free msgmem")
  #toc("waitSends: wait and free")
  for i in (i1+1)..<c.smsg.len:
    let k = i - i1 - 1 + i0
    c.smem[k] = c.smem[i]
    c.smsg[k] = c.smsg[i]
  #toc("waitSends: copy")
  let n = c.smsg.len - (i1-i0+1)
  c.smem.setLen(n)
  c.smsg.setLen(n)
  #toc("waitSends: setLen")

method freeRecvs*(c: CommQmp; i0,i1: int) =
  for i in i0..i1:
    QMP_free_msghandle(c.rmsg[i])
    QMP_free_msgmem(c.rmem[i])
  for i in (i1+1)..<c.rmsg.len:
    let k = i - i1 - 1 + i0
    c.rmem[k] = c.rmem[i]
    c.rmsg[k] = c.rmsg[i]
  let n = c.rmsg.len - (i1-i0+1)
  c.rmem.setLen(n)
  c.rmsg.setLen(n)

method waitRecvs*(c: CommQmp; i0,i1: int; free=true) =
  for i in i0..i1:
    discard QMP_wait(c.rmsg[i])
    if free:
      QMP_free_msghandle(c.rmsg[i])
      QMP_free_msgmem(c.rmem[i])
  if free:
    for i in (i1+1)..<c.rmsg.len:
      let k = i - i1 - 1 + i0
      c.rmem[k] = c.rmem[i]
      c.rmsg[k] = c.rmsg[i]
    let n = c.rmsg.len - (i1-i0+1)
    c.rmem.setLen(n)
    c.rmsg.setLen(n)

when isMainModule:
  commsInit()
  echo "rank ", myRank, "/", nRanks
  printf("rank %i/%i\n", myRank, nRanks)

  var c = getComm()
  let orank = 1-myRank

  var sx = 10
  var rx = 0
  c.pushRecv(orank, rx)
  c.pushSend(orank, sx)
  c.waitAll
  echo rx

  var sy = @[1.0, 2, 3, 4, 5, 6, 7]
  var ry = newSeq[float](sy.len)
  c.pushRecv(orank, ry)
  c.pushSend(orank, sy)
  c.waitAll
  echo ry

  commsFinalize()
