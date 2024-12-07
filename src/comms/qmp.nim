import os

when existsEnv("QMPDIR"):
  const qmpDir {.strDefine.} = getEnv("QMPDIR")
else:
  const qmpDir {.strDefine.} = getHomeDir() & "lqcd/install/qmp"
const qmpPassC = "-I" & qmpDir & "/include"
const qmpPassL* = "-L" & qmpDir & "/lib -lqmp -Wl,-rpath," & qmpDir & "/lib"
static:
  echo "Using QMP: ", qmpDir
  echo "QMP compile flags: ", qmpPassC
  echo "QMP link flags: ", qmpPassL
{. passC: qmpPassC .}
{. passL: qmpPassL .}
{. pragma: qmp, importc, header:"qmp.h" .}

type QMP_status_t*{.qmp.} = object
var QMP_SUCCESS*{.qmp.}:QMP_status_t

type QMP_thread_level_t*{.qmp.} = enum
  QMP_THREAD_SINGLE,
  QMP_THREAD_FUNNELED,
  QMP_THREAD_SERIALIZED,
  QMP_THREAD_MULTIPLE
type
  QMP_comm_t*{.qmp.} = pointer
  QMP_msgmem_t*{.qmp.} = pointer
  QMP_msghandle_t*{.qmp.} = object

proc clear*(x: var QMP_msghandle_t) =
  var p = cast[pointer](x)
  p = nil

proc isEmpty*(x: QMP_msghandle_t): bool =
  var p = cast[pointer](x)
  p == nil

proc QMP_init_msg_passing*(argc:ptr cint; argv:ptr ptr cstring;
                           required:QMP_thread_level_t;
                           provided:ptr QMP_thread_level_t):QMP_status_t{.qmp.}
proc QMP_finalize_msg_passing*() {.qmp.}
proc QMP_abort*(error_code:cint) {.qmp.}
proc QMP_get_number_of_nodes*():cint {.qmp.}
proc QMP_get_node_number*():cint {.qmp.}
proc QMP_barrier*() {.qmp.}
proc QMP_broadcast*(buffer: pointer, nbytes: csize_t) {.qmp.}
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
proc QMP_comm_broadcast*(comm: QMP_comm_t, buffer: pointer, nbytes: csize_t) {.qmp.}

template QMP_comm_broadcast*(c: QMP_comm_t, b: pointer, n: int) =
  QMP_comm_broadcast(c, b, (csize_t)n)

proc QMP_declare_msgmem*(mem: pointer; nbytes: csize_t): QMP_msgmem_t {.
    importc: "QMP_declare_msgmem", header: "qmp.h".}
proc QMP_declare_send_to*(m: QMP_msgmem_t; rem_node_rank: cint;
                          priority: cint): QMP_msghandle_t {.
    importc: "QMP_declare_send_to", header: "qmp.h".}
proc QMP_comm_declare_send_to*(comm: QMP_comm_t; m: QMP_msgmem_t;
                               rem_node_rank: cint; priority: cint):
                                 QMP_msghandle_t {.
    importc: "QMP_comm_declare_send_to", header: "qmp.h".}
proc QMP_declare_receive_from*(m: QMP_msgmem_t; rem_node_rank: cint;
                               priority: cint): QMP_msghandle_t {.
    importc: "QMP_declare_receive_from", header: "qmp.h".}
proc QMP_comm_declare_receive_from*(comm: QMP_comm_t; m: QMP_msgmem_t;
                                    rem_node_rank: cint; priority: cint):
                                      QMP_msghandle_t {.
    importc: "QMP_comm_declare_receive_from", header: "qmp.h".}
proc QMP_declare_send_recv_pairs*(msgh: ptr QMP_msghandle_t;
                                  num: cint): QMP_msghandle_t {.
    importc: "QMP_declare_send_recv_pairs", header: "qmp.h".}
proc QMP_declare_multiple*(msgh: ptr QMP_msghandle_t; num: cint): QMP_msghandle_t {.
    importc: "QMP_declare_multiple", header: "qmp.h".}
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
