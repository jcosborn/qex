type
  Comm* = ref object of RootObj

var defaultComm*: Comm
template getDefaultComm*(): Comm = defaultComm

var myRank* = 0
var nRanks* = 1

# base methods

method commrank*(c: Comm): int {.base.} = discard
method commsize*(c: Comm): int {.base.} = discard
method isMaster*(c: Comm): bool {.base.} = discard

method barrier*(c: Comm) {.base.} = discard
method broadcast*(c: Comm, p: pointer, bytes: int) {.base.} = discard
method allReduce*(c: Comm, x: ptr float32, n: int) {.base.} = discard
method allReduce*(c: Comm, x: ptr float64, n: int) {.base.} = discard
method allReduceXor*(c: Comm, x: var int) {.base.} = discard
# max, min

method nsends*(c: Comm): int {.base.} = discard
method nrecvs*(c: Comm): int {.base.} = discard
method pushSend*(c: Comm, rank: int, p: pointer, bytes: int) {.base.} =
  discard
method pushRecv*(c: Comm, rank: int, p: pointer, bytes: int) {.base.} =
  discard
method waitSends*(c: Comm; i0,i1: int) {.base.} = discard
method freeRecvs*(c: Comm; i0,i1: int) {.base.} = discard
method waitRecvs*(c: Comm; i0,i1: int; free=true) {.base.} = discard
# freeSends
# declareSend, declareRecv, pair(nsends, nrecvs), start(nsends, nrecvs)
# dup, split

# convenience templates

template `&&`(x: typed): untyped = cast[pointer](unsafeAddr x)

template rank*(c: Comm): int = commrank(c)
template size*(c: Comm): int = commsize(c)

template broadcast*(c: Comm, x: var SomeNumber) =
  c.broadcast(cast[pointer](addr x), sizeof(x))
template broadcast*(c: Comm, x: var string) =
  var n = x.len
  c.broadcast(n)
  x.setLen(n)
  c.broadcast(cast[pointer](addr x[0]), n)

template allReduce*(c: Comm, x: var float32) = c.allReduce(addr x, 1)
template allReduce*(c: Comm, x: var float64) = c.allReduce(addr x, 1)
template allReduce*(c: Comm, x: var SomeInteger) =
  var t = float(x)
  c.allreduce(t)
  x = (type x)(t)
template allReduce*(c: Comm, x: var UncheckedArray[float32], n: int) =
  c.allReduce(addr x[0], n.cint)
template allReduce*(c: Comm, x: var UncheckedArray[float64], n: int) =
  c.allReduce(addr x[0], n.cint)

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
  ## frees last k receives
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

# commBuffer
# declareMessageGroup
# MessageGroup: startRecvs, startSends, start, wait
#   getRecvBuf(nrecv), getSendBuf(nsend)

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

import commsQmp
export commsQmp

import commsUtils
export commsUtils

proc commsInit*() =
  commsInitQmp()
  defaultComm = getComm()
  myRank = defaultComm.rank
  nRanks = defaultComm.size


when isMainModule:
  commsInit()
  #echoAll "rank ", myRank, "/", nRanks
  #printf("rank %i/%i\n", myRank, nRanks)

  var c = getComm()
  let orank = 1-myRank

  c.barrier

  var bcst = myRank + 10
  echoAll myRank, "/", nRanks, " bcst: ", bcst
  c.broadcast(bcst)
  echoAll myRank, "/", nRanks, " bcst: ", bcst

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
