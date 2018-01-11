import macros

macro applyCall(f: untyped, idents: untyped): untyped =
  result = newStmtList()
  for i in 0..<idents.len:
    result.add newCall(f, idents[i])
  #echo result.repr
macro applyCall2(f: untyped, idents: untyped): untyped =
  #echo idents.treerepr
  result = newStmtList()
  for i in 0..<idents.len:
    var t = idents[i][1]
    if t.kind == nnkStmtList: t = t[0]
    result.add newCall(f, idents[i][0], t)
  #echo result.treerepr

{. pragma: mpih, importc, header: "mpi.h" .}

# import opaque types
template tdef(i: untyped) {.dirty.} =
  type i* {.mpih.} = object
applyCall tdef:
  MPI_Aint
  MPI_Count
  MPI_Comm
  MPI_Datatype
  MPI_Errhandler
  MPI_File
  MPI_Group
  MPI_Info
  MPI_Op
  MPI_Request
  MPI_Message
  MPI_Status
  MPI_Win

type MPI_Offset* {.mpih.} = cint

# import constants
template objdef(i,t: untyped) {.dirty.} =
  var i* {.mpih.}: t
applyCall2 objdef:
  MPI_ANY_SOURCE: cint
  MPI_PROC_NULL: cint
  MPI_ROOT: MPI_Comm
  MPI_ANY_TAG: cint
  MPI_MAX_PROCESSOR_NAME: cint
  MPI_MAX_ERROR_STRING: cint
  MPI_MAX_OBJECT_NAME: cint
  MPI_MAX_LIBRARY_VERSION_STRING: cint
  MPI_UNDEFINED: cint
  MPI_DIST_GRAPH: cint
  MPI_CART: cint
  MPI_GRAPH: cint
  MPI_KEYVAL_INVALID: cint

  MPI_UNWEIGHTED: ptr cint
  MPI_WEIGHTS_EMPTY: ptr cint
  MPI_BOTTOM: pointer
  MPI_IN_PLACE: pointer
  MPI_BSEND_OVERHEAD: cint
  MPI_MAX_INFO_KEY: cint
  MPI_MAX_INFO_VAL: cint
  MPI_ARGV_NULL: ptr ptr cchar
  MPI_ARGVS_NULL: ptr ptr ptr cchar
  MPI_ERRCODES_IGNORE: ptr cint
  MPI_MAX_PORT_NAME: cint
  MPI_ORDER_C: cint
  MPI_ORDER_FORTRAN: cint
  MPI_DISTRIBUTE_BLOCK: cint
  MPI_DISTRIBUTE_CYCLIC: cint
  MPI_DISTRIBUTE_NONE: cint
  MPI_DISTRIBUTE_DFLT_DARG: cint

  MPI_MODE_CREATE: cint
  MPI_MODE_RDONLY: cint
  MPI_MODE_WRONLY: cint
  MPI_MODE_RDWR: cint
  MPI_MODE_DELETE_ON_CLOSE: cint
  MPI_MODE_UNIQUE_OPEN: cint
  MPI_MODE_EXCL: cint
  MPI_MODE_APPEND: cint
  MPI_MODE_SEQUENTIAL: cint
  MPI_DISPLACEMENT_CURRENT: cint
  MPI_SEEK_SET: cint
  MPI_SEEK_CUR: cint
  MPI_SEEK_END: cint
  MPI_MAX_DATAREP_STRING: cint

  MPI_MODE_NOCHECK: cint
  MPI_MODE_NOPRECEDE: cint
  MPI_MODE_NOPUT: cint
  MPI_MODE_NOSTORE: cint
  MPI_MODE_NOSUCCEED: cint
  MPI_LOCK_EXCLUSIVE: cint
  MPI_LOCK_SHARED: cint
  MPI_WIN_FLAVOR_CREATE: cint
  MPI_WIN_FLAVOR_ALLOCATE: cint
  MPI_WIN_FLAVOR_DYNAMIC: cint
  MPI_WIN_FLAVOR_SHARED: cint
  MPI_WIN_UNIFIED: cint
  MPI_WIN_SEPARATE: cint

  MPI_TAG_UB: cint
  MPI_HOST: cint
  MPI_IO: cint
  MPI_WTIME_IS_GLOBAL: cint
  MPI_APPNUM: cint
  MPI_LASTUSEDCODE: cint
  MPI_UNIVERSE_SIZE: cint
  MPI_WIN_BASE: cint
  MPI_WIN_SIZE: cint
  MPI_WIN_DISP_UNIT: cint
  MPI_WIN_CREATE_FLAVOR: cint
  MPI_WIN_MODEL: cint

  MPI_SUCCESS: cint
  MPI_ERR_BUFFER: cint
  MPI_ERR_COUNT: cint
  MPI_ERR_TYPE: cint
  MPI_ERR_TAG: cint
  MPI_ERR_COMM: cint
  MPI_ERR_RANK: cint
  MPI_ERR_REQUEST: cint
  MPI_ERR_ROOT: cint
  MPI_ERR_GROUP: cint
  MPI_ERR_OP: cint
  MPI_ERR_TOPOLOGY: cint
  MPI_ERR_DIMS: cint
  MPI_ERR_ARG: cint
  MPI_ERR_UNKNOWN: cint
  MPI_ERR_TRUNCATE: cint
  MPI_ERR_OTHER: cint
  MPI_ERR_INTERN: cint
  MPI_ERR_IN_STATUS: cint
  MPI_ERR_PENDING: cint
  MPI_ERR_ACCESS: cint
  MPI_ERR_AMODE: cint
  MPI_ERR_ASSERT: cint
  MPI_ERR_BAD_FILE: cint
  MPI_ERR_BASE: cint
  MPI_ERR_CONVERSION: cint
  MPI_ERR_DISP: cint
  MPI_ERR_DUP_DATAREP: cint
  MPI_ERR_FILE_EXISTS: cint
  MPI_ERR_FILE_IN_USE: cint
  MPI_ERR_FILE: cint
  MPI_ERR_INFO_KEY: cint
  MPI_ERR_INFO_NOKEY: cint
  MPI_ERR_INFO_VALUE: cint
  MPI_ERR_INFO: cint
  MPI_ERR_IO: cint
  MPI_ERR_KEYVAL: cint
  MPI_ERR_LOCKTYPE: cint
  MPI_ERR_NAME: cint
  MPI_ERR_NO_MEM: cint
  MPI_ERR_NOT_SAME: cint
  MPI_ERR_NO_SPACE: cint
  MPI_ERR_NO_SUCH_FILE: cint
  MPI_ERR_PORT: cint
  MPI_ERR_QUOTA: cint
  MPI_ERR_READ_ONLY: cint
  MPI_ERR_RMA_CONFLICT: cint
  MPI_ERR_RMA_SYNC: cint
  MPI_ERR_SERVICE: cint
  MPI_ERR_SIZE: cint
  MPI_ERR_SPAWN: cint
  MPI_ERR_UNSUPPORTED_DATAREP: cint
  MPI_ERR_UNSUPPORTED_OPERATION: cint
  MPI_ERR_WIN: cint
  MPI_ERR_RMA_RANGE: cint
  MPI_ERR_RMA_ATTACH: cint
  MPI_ERR_RMA_FLAVOR: cint
  MPI_ERR_RMA_SHARED: cint
  MPI_ERR_LASTCODE: cint

  MPI_IDENT: cint
  MPI_CONGRUENT: cint
  MPI_SIMILAR: cint
  MPI_UNEQUAL: cint

  MPI_THREAD_SINGLE: cint
  MPI_THREAD_FUNNELED: cint
  MPI_THREAD_SERIALIZED: cint
  MPI_THREAD_MULTIPLE: cint

  MPI_COMBINER_NAMED: cint
  MPI_COMBINER_DUP: cint
  MPI_COMBINER_CONTIGUOUS: cint
  MPI_COMBINER_VECTOR: cint
  MPI_COMBINER_HVECTOR_INTEGER: cint
  MPI_COMBINER_HVECTOR: cint
  MPI_COMBINER_INDEXED: cint
  MPI_COMBINER_HINDEXED_INTEGER: cint
  MPI_COMBINER_HINDEXED: cint
  MPI_COMBINER_INDEXED_BLOCK: cint
  MPI_COMBINER_STRUCT_INTEGER: cint
  MPI_COMBINER_STRUCT: cint
  MPI_COMBINER_SUBARRAY: cint
  MPI_COMBINER_DARRAY: cint
  MPI_COMBINER_F90_REAL: cint
  MPI_COMBINER_F90_COMPLEX: cint
  MPI_COMBINER_F90_INTEGER: cint
  MPI_COMBINER_RESIZED: cint
  MPI_COMBINER_HINDEXED_BLOCK: cint

  MPI_COMM_TYPE_SHARED: cint

  MPI_GROUP_NULL: MPI_Group
  MPI_COMM_NULL: MPI_Comm
  MPI_REQUEST_NULL: MPI_Request
  MPI_MESSAGE_NULL: MPI_Message
  MPI_OP_NULL: MPI_Op
  MPI_ERRHANDLER_NULL: MPI_Errhandler
  MPI_INFO_NULL: MPI_Info
  MPI_WIN_NULL: MPI_Win
  MPI_FILE_NULL: MPI_File

  MPI_INFO_ENV: MPI_Info

  MPI_STATUS_IGNORE: ptr MPI_Status
  MPI_STATUSES_IGNORE: ptr MPI_Status

  MPI_COMM_WORLD: MPI_Comm
  MPI_COMM_SELF: MPI_Comm

  MPI_GROUP_EMPTY: MPI_Group

  MPI_MESSAGE_NO_PROC: MPI_Message

  MPI_MAX: MPI_Op
  MPI_MIN: MPI_Op
  MPI_SUM: MPI_Op
  MPI_PROD: MPI_Op
  MPI_LAND: MPI_Op
  MPI_BAND: MPI_Op
  MPI_LOR: MPI_Op
  MPI_BOR: MPI_Op
  MPI_LXOR: MPI_Op
  MPI_BXOR: MPI_Op
  MPI_MAXLOC: MPI_Op
  MPI_MINLOC: MPI_Op
  MPI_REPLACE: MPI_Op
  MPI_NO_OP: MPI_Op

  MPI_DATATYPE_NULL: MPI_Datatype
  MPI_BYTE: MPI_Datatype
  MPI_PACKED: MPI_Datatype
  MPI_CHAR: MPI_Datatype
  MPI_SHORT: MPI_Datatype
  MPI_INT: MPI_Datatype
  MPI_LONG: MPI_Datatype
  MPI_FLOAT: MPI_Datatype
  MPI_DOUBLE: MPI_Datatype
  MPI_LONG_DOUBLE: MPI_Datatype
  MPI_UNSIGNED_CHAR: MPI_Datatype
  MPI_SIGNED_CHAR: MPI_Datatype
  MPI_UNSIGNED_SHORT: MPI_Datatype
  MPI_UNSIGNED_LONG: MPI_Datatype
  MPI_UNSIGNED: MPI_Datatype
  MPI_FLOAT_INT: MPI_Datatype
  MPI_DOUBLE_INT: MPI_Datatype
  MPI_LONG_DOUBLE_INT: MPI_Datatype
  MPI_LONG_INT: MPI_Datatype
  MPI_SHORT_INT: MPI_Datatype
  MPI_2INT: MPI_Datatype
  MPI_UB: MPI_Datatype
  MPI_LB: MPI_Datatype
  MPI_WCHAR: MPI_Datatype
  MPI_LONG_LONG_INT: MPI_Datatype
  MPI_LONG_LONG: MPI_Datatype
  MPI_UNSIGNED_LONG_LONG: MPI_Datatype
  MPI_2COMPLEX: MPI_Datatype
  MPI_2DOUBLE_COMPLEX: MPI_Datatype

  MPI_INT8_T: MPI_Datatype
  MPI_UINT8_T: MPI_Datatype
  MPI_INT16_T: MPI_Datatype
  MPI_UINT16_T: MPI_Datatype
  MPI_INT32_T: MPI_Datatype
  MPI_UINT32_T: MPI_Datatype
  MPI_INT64_T: MPI_Datatype
  MPI_UINT64_T: MPI_Datatype
  MPI_C_BOOL: MPI_Datatype
  MPI_C_COMPLEX: MPI_Datatype
  MPI_C_FLOAT_COMPLEX: MPI_Datatype
  MPI_C_DOUBLE_COMPLEX: MPI_Datatype
  MPI_C_LONG_DOUBLE_COMPLEX: MPI_Datatype
  MPI_CXX_BOOL: MPI_Datatype
  MPI_CXX_COMPLEX: MPI_Datatype
  MPI_CXX_FLOAT_COMPLEX: MPI_Datatype
  MPI_CXX_DOUBLE_COMPLEX: MPI_Datatype
  MPI_CXX_LONG_DOUBLE_COMPLEX: MPI_Datatype

  MPI_ERRORS_ARE_FATAL: MPI_Errhandler
  MPI_ERRORS_RETURN: MPI_Errhandler

  MPI_TYPECLASS_INTEGER: cint
  MPI_TYPECLASS_REAL: cint
  MPI_TYPECLASS_COMPLEX: cint

var MPI_AINT_Datatype* {.importc: "MPI_AINT", header: "mpi.h".}: MPI_Datatype
var MPI_OFFSET_Datatype* {.importc: "MPI_OFFSET", header: "mpi.h".}: MPI_Datatype
var MPI_COUNT_Datatype* {.importc: "MPI_COUNT", header: "mpi.h".}: MPI_Datatype
