import qmp

var inited = false

proc commsInitQmp* =
  var argc {.importc:"cmdCount", global.}:cint
  var argv {.importc:"cmdLine", global.}:ptr cstring
  if not inited:
    inited = true
    var prv = QMP_THREAD_FUNNELED
    #var prv = QMP_THREAD_SERIALIZED
    let err = QMP_init_msg_passing(argc.addr, argv.addr, prv, prv.addr)
    doAssert(err == QMP_SUCCESS)
    #discard err
    #myRank = int(QMP_get_node_number())
    #nRanks = int(QMP_get_number_of_nodes())
    #defaultComm = getComm()

proc commsFinalizeQmp* =
  if inited:
    GC_fullCollect()  # attempt to free any dangling message handles
    QMP_finalize_msg_passing()

proc commsAbortQmp*(status = -1) =
  QMP_abort(status.cint)

proc qmpSum*(v:var int) =
  var t = v.float
  QmpSumDouble(t.addr)
  v = t.int

template qmpSum*(v:float32):untyped = QmpSumFloat(v.addr)
template qmpSum*(v:float64):untyped = QmpSumDouble(v.addr)
template qmpSum*(v:ptr float32, n:int):untyped = QmpSumFloatArray(v,n.cint)
template qmpSum*(v:ptr float64, n:int):untyped = QmpSumDoubleArray(v,n.cint)
#template qmpSum*(v:ptr array, n:int):untyped =
#  qmpSum(v[][0].addr, n*v[].len)
#template qmpSum*(v:ptr tuple, n:int):untyped =
#  qmpSum(v[][0].addr, n*(sizeOf(v) div sizeOf(v[0])))
#template qmpSum*(v:ptr object, n:int):untyped =
#  qmpSum(v[][].addr, n)
#template qmpSum*(v: object) =
#  qmpSum(asNumberPtr(v), numNumbers(v))
#template qmpSum*(v:ptr typed, n:int):untyped =
#  qmpSum(v[][].addr, n)
#template QmpSum(v:array[int,int]):untyped =
#  var tQmpSumDoubleArray(v)
template qmpSum*[I,T](v:array[I,T]):untyped =
  qmpSum(v[0].addr, v.len)
#template qmpSum*(v:openArray[float64]):untyped =
#  QmpSumDoubleArray(v[0].addr,v.len.cint)
template qmpSum*[T](v:seq[T]):untyped =
  qmpSum(v[0].addr, v.len)
#template qmpSum*[I,T](v:seq[array[I,T]]):untyped =
#  qmpSum(v[0][0].addr, v.len.cint*sizeOf(v[0]))
#template qmpSum*(v:openArray[array]):untyped =
#  qmpSum(v[0][0].addr, v.len.cint*sizeOf(v[0]))
template qmpSum*(v:tuple):untyped =
  qmpSum(v[0].addr, sizeOf(v) div sizeOf(v[0]))
#template qmpSum*[T](v:T):untyped =
#template qmpSum*(v:typed):untyped =
#  qmpSum(v[])
#template qmpSum*[T](v:T):untyped =
#  qmpSum(v[])
template qmpSum*(v: typed): untyped =
  when numberType(v) is float64:
    qmpSum(cast[ptr float64](addr v), sizeof(v) div sizeof(float64))
  elif numberType(v) is float32:
    qmpSum(cast[ptr float32](addr v), sizeof(v) div sizeof(float32))
  else:
    qmpSum(v[])

template qmpMax*(v:float32):untyped = QmpMaxFloat(v.addr)
template qmpMax*(v:float64):untyped = QmpMaxDouble(v.addr)
template qmpMin*(v:float32):untyped = QmpMinFloat(v.addr)
template qmpMin*(v:float64):untyped = QmpMinDouble(v.addr)

# generic comms interface

import commsTypes

type
  CommQmp* = ref object of Comm
    comm: QMP_comm_t
    smem: seq[QMP_msgmem_t]
    smsg: seq[QMP_msghandle_t]
    rmem: seq[QMP_msgmem_t]
    rmsg: seq[QMP_msghandle_t]

proc getCommQmp*(): CommQmp =
  commsInitQmp()
  result.new
  result.comm = QMP_comm_get_default()
  result.smem.newSeq(0)
  result.smsg.newSeq(0)
  result.rmem.newSeq(0)
  result.rmsg.newSeq(0)

proc getQmpComm*(): Comm =
  result = getCommQmp()

method name*(c: CommQmp): string = "QMP"

method commrank*(c: CommQmp): int =
  QMP_comm_get_node_number(c.comm).int

method commsize*(c: CommQmp): int =
  QMP_comm_get_number_of_nodes(c.comm).int

method isMaster*(c: CommQmp): bool =
  c.commrank == 0

method abort*(c: CommQmp, status: int) = commsAbortQmp(status)

method barrier*(c: CommQmp) =
  QMP_comm_barrier(c.comm)

method broadcast*(c: CommQMP, p: pointer, bytes: int) =
  QMP_comm_broadcast(c.comm, p, bytes)

#method allReduce*(c: CommQmp, x: var float64) =
#  QMP_comm_sum_double(c.comm, addr x)

method allReduce*(c: CommQmp, x: ptr float32, n: int) =
  QMP_comm_sum_float_array(c.comm, x, n.cint)

method allReduce*(c: CommQmp, x: ptr float64, n: int) =
  QMP_comm_sum_double_array(c.comm, x, n.cint)

method allReduceXor*(c: CommQmp, x: var int) =
  var t = cast[ptr culong](addr x)
  QMP_comm_xor_ulong(c.comm, t)

method nsends*(c: CommQmp): int = c.smsg.len
method nrecvs*(c: CommQmp): int = c.rmsg.len

method pushSend*(c: CommQmp, rank: int, p: pointer, bytes: int) =
  let m = QMP_declare_msgmem(p, bytes.csize_t)
  let h = QMP_comm_declare_send_to(c.comm, m, rank.cint, 0.cint)
  let stat = QMP_start(h)
  discard stat
  c.smem.add m
  c.smsg.add h

method pushRecv*(c: CommQmp, rank: int, p: pointer, bytes: int) =
  let m = QMP_declare_msgmem(p, bytes.csize_t)
  let h = QMP_comm_declare_receive_from(c.comm, m, rank.cint, 0.cint)
  let stat = QMP_start(h)
  discard stat
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
  if c.smsg.len > 0:
    c.smem.setLen(n)
    c.smsg.setLen(n)
  else:
    c.smem.newSeq(n)
    c.smsg.newSeq(n)
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
  if c.rmsg.len > 0:
    c.rmem.setLen(n)
    c.rmsg.setLen(n)
  else:
    c.rmem.newSeq(n)
    c.rmsg.newSeq(n)

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
    if c.rmsg.len > 0:
      c.rmem.setLen(n)
      c.rmsg.setLen(n)
    else:
      c.rmem.newSeq(n)
      c.rmsg.newSeq(n)




