import os
import macros

when existsEnv("QMPDIR"):
  const qmpDir = getEnv("QMPDIR")
else:
  const homeDir = getHomeDir()
  const qmpDir = homeDir & "lqcd/install/qmp"
{. passC: "-I" & qmpDir & "/include" .}
{. passL: "-L" & qmpDir & "/lib -lqmp" .}
{. pragma: qmp, importc, header:"qmp.h" .}

type QMP_status_t{.qmp.} = enum
  test
type QMP_thread_level_t*{.qmp.} = enum
  QMP_THREAD_SINGLE,
  QMP_THREAD_FUNNELED,
  QMP_THREAD_SERIALIZED,
  QMP_THREAD_MULTIPLE
type
  QMP_comm_t*{.qmp.} = pointer
  QMP_msgmem_t*{.qmp.} = pointer
  QMP_msghandle_t*{.qmp.} = object

var emptyQmpMsgHandle: QMP_msghandle_t

proc clear*(x: var QMP_msghandle_t) =
  x = emptyQmpMsgHandle

proc isEmpty*(x: QMP_msghandle_t): bool =
  x == emptyQmpMsgHandle

proc QMP_init_msg_passing*(argc:ptr cint; argv:ptr ptr cstring;
                           required:QMP_thread_level_t;
                           provided:ptr QMP_thread_level_t):QMP_status_t{.qmp.}
proc QMP_finalize_msg_passing*() {.qmp.}
proc QMP_abort*(error_code:cint) {.qmp.}
proc QMP_get_number_of_nodes*():cint {.qmp.}
proc QMP_get_node_number*():cint {.qmp.}
proc QMP_barrier*() {.qmp.}
proc QMP_sum_float*(value:ptr cfloat) {.qmp.}
proc QMP_sum_double*(value:ptr cdouble) {.qmp.}
proc QMP_comm_sum_double*(comm: QMP_comm_t, value:ptr cdouble) {.qmp.}
proc QMP_sum_float_array*(value:ptr cfloat, length:cint) {.qmp.}
proc QMP_sum_double_array*(value:ptr cdouble, length:cint) {.qmp.}
proc QMP_comm_sum_float_array*(comm: QMP_comm_t, value:ptr cfloat, length:cint) {.qmp.}
proc QMP_comm_sum_double_array*(comm: QMP_comm_t, value:ptr cdouble, length:cint) {.qmp.}
proc QMP_max_float*(value:ptr cfloat) {.qmp.}
proc QMP_max_double*(value:ptr cdouble) {.qmp.}
proc QMP_min_float*(value:ptr cfloat) {.qmp.}
proc QMP_min_double*(value:ptr cdouble) {.qmp.}
proc QMP_comm_xor_ulong*(comm: QMP_comm_t, value: ptr culong) {.qmp.}
proc QMP_comm_get_default*(): QMP_comm_t {.qmp.}
proc QMP_comm_get_number_of_nodes*(comm: QMP_comm_t):cint {.qmp.}
proc QMP_comm_get_node_number*(comm: QMP_comm_t):cint {.qmp.}
proc QMP_comm_barrier*(comm: QMP_comm_t) {.qmp.}




proc qmpSum*(v:var int) =
  var t = v.float
  QmpSumDouble(t.addr)
  v = t.int

template qmpSum*(v:float32):untyped = QmpSumFloat(v.addr)
template qmpSum*(v:float64):untyped = QmpSumDouble(v.addr)
template qmpSum*(v:ptr float32, n:int):untyped = QmpSumFloatArray(v,n.cint)
template qmpSum*(v:ptr float64, n:int):untyped = QmpSumDoubleArray(v,n.cint)
template qmpSum*(v:ptr array, n:int):untyped =
  qmpSum(v[][0].addr, n*v[].len)
template qmpSum*(v:ptr tuple, n:int):untyped =
  qmpSum(v[][0].addr, n*(sizeOf(v) div sizeOf(v[0])))
template qmpSum*(v:ptr object, n:int):untyped =
  qmpSum(v[][].addr, n)
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

proc QMP_declare_msgmem*(mem: pointer; nbytes: csize): QMP_msgmem_t {.
    importc: "QMP_declare_msgmem", header: "qmp.h".}
proc QMP_declare_send_to*(m: QMP_msgmem_t; rem_node_rank: cint;
                          priority: cint): QMP_msghandle_t {.
    importc: "QMP_declare_send_to", header: "qmp.h".}
proc QMP_declare_receive_from*(m: QMP_msgmem_t; rem_node_rank: cint;
                               priority: cint): QMP_msghandle_t {.
    importc: "QMP_declare_receive_from", header: "qmp.h".}
proc QMP_declare_send_recv_pairs*(msgh: ptr QMP_msghandle_t;
                                  num: cint): QMP_msghandle_t {.
    importc: "QMP_declare_send_recv_pairs", header: "qmp.h".}
proc QMP_start*(h: QMP_msghandle_t): QMP_status_t {.importc: "QMP_start",
    header: "qmp.h".}
proc QMP_wait*(h: QMP_msghandle_t): QMP_status_t {.importc: "QMP_wait",
    header: "qmp.h".}
type QMP_clear_to_send_t* {.importC, header:"qmp.h".} = cint
var QMP_CTS_DISABLED*{.importC, header:"qmp.h".}: QMP_clear_to_send_t
var QMP_CTS_NOT_READY*{.importC, header:"qmp.h".}: QMP_clear_to_send_t
var QMP_CTS_READY*{.importC, header:"qmp.h".}: QMP_clear_to_send_t
proc QMP_clear_to_send*(mh: QMP_msghandle_t;
                        cts: QMP_clear_to_send_t): QMP_status_t {.
    importc: "QMP_clear_to_send", header: "qmp.h".}
proc QMP_free_msghandle*(h: QMP_msghandle_t) {.importc: "QMP_free_msghandle",
    header: "qmp.h".}
proc QMP_free_msgmem*(m: QMP_msgmem_t) {.importc: "QMP_free_msgmem", header: "qmp.h".}


when isMainModule:
  var argc {.importc:"cmdCount", global.}:cint
  var argv {.importc:"cmdLine", global.}:ptr cstring
  var prv = QMP_THREAD_SERIALIZED
  let err = QMP_init_msg_passing(argc.addr, argv.addr, prv, prv.addr)
  let rank = QMP_get_node_number()
  let size = QMP_get_number_of_nodes()
  echo "rank " & $rank & "/" & $size
  QMP_finalize_msg_passing()
