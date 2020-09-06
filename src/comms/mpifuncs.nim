# mpifuncs.nim
# created with c2nim from mpifuncs.h
import mpitypes

type 
  MPI_Copy_function* = proc (a2: MPI_Comm; a3: cint; a4: pointer; a5: pointer; 
                             a6: pointer; a7: ptr cint): cint
  MPI_Delete_function* = proc (a2: MPI_Comm; a3: cint; a4: pointer; a5: pointer): cint
  MPI_Datarep_extent_function* = proc (a2: MPI_Datatype; a3: ptr MPI_Aint; 
                                       a4: pointer): cint
  MPI_Datarep_conversion_function* = proc (a2: pointer; a3: MPI_Datatype; 
      a4: cint; a5: pointer; a6: MPI_Offset; a7: pointer): cint
  MPI_Comm_errhandler_function* = proc (a2: ptr MPI_Comm; a3: ptr cint) {.
      varargs.}
  MPI_Comm_errhandler_fn* = MPI_Comm_errhandler_function
  MPI_File_errhandler_function* = proc (a2: ptr MPI_File; a3: ptr cint) {.
      varargs.}
  MPI_Win_errhandler_function* = proc (a2: ptr MPI_Win; a3: ptr cint) {.varargs.}
  MPI_Win_errhandler_fn* = MPI_Win_errhandler_function
  MPI_Handler_function* = proc (a2: ptr MPI_Comm; a3: ptr cint) {.varargs.}
  MPI_User_function* = proc (a2: pointer; a3: pointer; a4: ptr cint; 
                             a5: ptr MPI_Datatype)
  MPI_Comm_copy_attr_function* = proc (a2: MPI_Comm; a3: cint; a4: pointer; 
                                       a5: pointer; a6: pointer; a7: ptr cint): cint
  MPI_Comm_delete_attr_function* = proc (a2: MPI_Comm; a3: cint; a4: pointer; 
      a5: pointer): cint
  MPI_Type_copy_attr_function* = proc (a2: MPI_Datatype; a3: cint; a4: pointer; 
                                       a5: pointer; a6: pointer; a7: ptr cint): cint
  MPI_Type_delete_attr_function* = proc (a2: MPI_Datatype; a3: cint; 
      a4: pointer; a5: pointer): cint
  MPI_Win_copy_attr_function* = proc (a2: MPI_Win; a3: cint; a4: pointer; 
                                      a5: pointer; a6: pointer; a7: ptr cint): cint
  MPI_Win_delete_attr_function* = proc (a2: MPI_Win; a3: cint; a4: pointer; 
                                        a5: pointer): cint
  MPI_Grequest_query_function* = proc (a2: pointer; a3: ptr MPI_Status): cint
  MPI_Grequest_free_function* = proc (a2: pointer): cint
  MPI_Grequest_cancel_function* = proc (a2: pointer; a3: cint): cint

proc MPI_Abort*(comm: MPI_Comm; errorcode: cint): cint {.importc: "MPI_Abort", 
    header: "mpi.h".}
proc MPI_Accumulate*(origin_addr: pointer; origin_count: cint; 
                     origin_datatype: MPI_Datatype; target_rank: cint; 
                     target_disp: MPI_Aint; target_count: cint; 
                     target_datatype: MPI_Datatype; op: MPI_Op; win: MPI_Win): cint {.
    importc: "MPI_Accumulate", header: "mpi.h".}
proc MPI_Add_error_class*(errorclass: ptr cint): cint {.
    importc: "MPI_Add_error_class", header: "mpi.h".}
proc MPI_Add_error_code*(errorclass: cint; errorcode: ptr cint): cint {.
    importc: "MPI_Add_error_code", header: "mpi.h".}
proc MPI_Add_error_string*(errorcode: cint; string: cstring): cint {.
    importc: "MPI_Add_error_string", header: "mpi.h".}
proc MPI_Address*(location: pointer; address: ptr MPI_Aint): cint {.
    importc: "MPI_Address", header: "mpi.h".}
proc MPI_Aint_add*(base: MPI_Aint; disp: MPI_Aint): MPI_Aint {.
    importc: "MPI_Aint_add", header: "mpi.h".}
proc MPI_Aint_diff*(addr1: MPI_Aint; addr2: MPI_Aint): MPI_Aint {.
    importc: "MPI_Aint_diff", header: "mpi.h".}
proc MPI_Allgather*(sendbuf: pointer; sendcount: cint; sendtype: MPI_Datatype; 
                    recvbuf: pointer; recvcount: cint; recvtype: MPI_Datatype; 
                    comm: MPI_Comm): cint {.importc: "MPI_Allgather", 
    header: "mpi.h".}
proc MPI_Iallgather*(sendbuf: pointer; sendcount: cint; sendtype: MPI_Datatype; 
                     recvbuf: pointer; recvcount: cint; recvtype: MPI_Datatype; 
                     comm: MPI_Comm; request: ptr MPI_Request): cint {.
    importc: "MPI_Iallgather", header: "mpi.h".}
proc MPI_Allgatherv*(sendbuf: pointer; sendcount: cint; sendtype: MPI_Datatype; 
                     recvbuf: pointer; recvcounts: ptr cint; displs: ptr cint; 
                     recvtype: MPI_Datatype; comm: MPI_Comm): cint {.
    importc: "MPI_Allgatherv", header: "mpi.h".}
proc MPI_Iallgatherv*(sendbuf: pointer; sendcount: cint; sendtype: MPI_Datatype; 
                      recvbuf: pointer; recvcounts: ptr cint; displs: ptr cint; 
                      recvtype: MPI_Datatype; comm: MPI_Comm; 
                      request: ptr MPI_Request): cint {.
    importc: "MPI_Iallgatherv", header: "mpi.h".}
proc MPI_Alloc_mem*(size: MPI_Aint; info: MPI_Info; baseptr: pointer): cint {.
    importc: "MPI_Alloc_mem", header: "mpi.h".}
proc MPI_Allreduce*(sendbuf: pointer; recvbuf: pointer; count: cint; 
                    datatype: MPI_Datatype; op: MPI_Op; comm: MPI_Comm): cint {.
    importc: "MPI_Allreduce", header: "mpi.h".}
proc MPI_Iallreduce*(sendbuf: pointer; recvbuf: pointer; count: cint; 
                     datatype: MPI_Datatype; op: MPI_Op; comm: MPI_Comm; 
                     request: ptr MPI_Request): cint {.
    importc: "MPI_Iallreduce", header: "mpi.h".}
proc MPI_Alltoall*(sendbuf: pointer; sendcount: cint; sendtype: MPI_Datatype; 
                   recvbuf: pointer; recvcount: cint; recvtype: MPI_Datatype; 
                   comm: MPI_Comm): cint {.importc: "MPI_Alltoall", 
    header: "mpi.h".}
proc MPI_Ialltoall*(sendbuf: pointer; sendcount: cint; sendtype: MPI_Datatype; 
                    recvbuf: pointer; recvcount: cint; recvtype: MPI_Datatype; 
                    comm: MPI_Comm; request: ptr MPI_Request): cint {.
    importc: "MPI_Ialltoall", header: "mpi.h".}
proc MPI_Alltoallv*(sendbuf: pointer; sendcounts: ptr cint; sdispls: ptr cint; 
                    sendtype: MPI_Datatype; recvbuf: pointer; 
                    recvcounts: ptr cint; rdispls: ptr cint; 
                    recvtype: MPI_Datatype; comm: MPI_Comm): cint {.
    importc: "MPI_Alltoallv", header: "mpi.h".}
proc MPI_Ialltoallv*(sendbuf: pointer; sendcounts: ptr cint; sdispls: ptr cint; 
                     sendtype: MPI_Datatype; recvbuf: pointer; 
                     recvcounts: ptr cint; rdispls: ptr cint; 
                     recvtype: MPI_Datatype; comm: MPI_Comm; 
                     request: ptr MPI_Request): cint {.
    importc: "MPI_Ialltoallv", header: "mpi.h".}
proc MPI_Alltoallw*(sendbuf: pointer; sendcounts: ptr cint; sdispls: ptr cint; 
                    sendtypes: ptr MPI_Datatype; recvbuf: pointer; 
                    recvcounts: ptr cint; rdispls: ptr cint; 
                    recvtypes: ptr MPI_Datatype; comm: MPI_Comm): cint {.
    importc: "MPI_Alltoallw", header: "mpi.h".}
proc MPI_Ialltoallw*(sendbuf: pointer; sendcounts: ptr cint; sdispls: ptr cint; 
                     sendtypes: ptr MPI_Datatype; recvbuf: pointer; 
                     recvcounts: ptr cint; rdispls: ptr cint; 
                     recvtypes: ptr MPI_Datatype; comm: MPI_Comm; 
                     request: ptr MPI_Request): cint {.
    importc: "MPI_Ialltoallw", header: "mpi.h".}
proc MPI_Attr_delete*(comm: MPI_Comm; keyval: cint): cint {.
    importc: "MPI_Attr_delete", header: "mpi.h".}
proc MPI_Attr_get*(comm: MPI_Comm; keyval: cint; attribute_val: pointer; 
                   flag: ptr cint): cint {.importc: "MPI_Attr_get", 
    header: "mpi.h".}
proc MPI_Attr_put*(comm: MPI_Comm; keyval: cint; attribute_val: pointer): cint {.
    importc: "MPI_Attr_put", header: "mpi.h".}
proc MPI_Barrier*(comm: MPI_Comm): cint {.importc: "MPI_Barrier", 
    header: "mpi.h".}
proc MPI_Ibarrier*(comm: MPI_Comm; request: ptr MPI_Request): cint {.
    importc: "MPI_Ibarrier", header: "mpi.h".}
proc MPI_Bcast*(buffer: pointer; count: cint; datatype: MPI_Datatype; 
                root: cint; comm: MPI_Comm): cint {.importc: "MPI_Bcast", 
    header: "mpi.h".}
proc MPI_Bsend*(buf: pointer; count: cint; datatype: MPI_Datatype; dest: cint; 
                tag: cint; comm: MPI_Comm): cint {.importc: "MPI_Bsend", 
    header: "mpi.h".}
proc MPI_Ibcast*(buffer: pointer; count: cint; datatype: MPI_Datatype; 
                 root: cint; comm: MPI_Comm; request: ptr MPI_Request): cint {.
    importc: "MPI_Ibcast", header: "mpi.h".}
proc MPI_Bsend_init*(buf: pointer; count: cint; datatype: MPI_Datatype; 
                     dest: cint; tag: cint; comm: MPI_Comm; 
                     request: ptr MPI_Request): cint {.
    importc: "MPI_Bsend_init", header: "mpi.h".}
proc MPI_Buffer_attach*(buffer: pointer; size: cint): cint {.
    importc: "MPI_Buffer_attach", header: "mpi.h".}
proc MPI_Buffer_detach*(buffer: pointer; size: ptr cint): cint {.
    importc: "MPI_Buffer_detach", header: "mpi.h".}
proc MPI_Cancel*(request: ptr MPI_Request): cint {.importc: "MPI_Cancel", 
    header: "mpi.h".}
proc MPI_Cart_coords*(comm: MPI_Comm; rank: cint; maxdims: cint; 
                      coords: ptr cint): cint {.importc: "MPI_Cart_coords", 
    header: "mpi.h".}
proc MPI_Cart_create*(old_comm: MPI_Comm; ndims: cint; dims: ptr cint; 
                      periods: ptr cint; reorder: cint; comm_cart: ptr MPI_Comm): cint {.
    importc: "MPI_Cart_create", header: "mpi.h".}
proc MPI_Cart_get*(comm: MPI_Comm; maxdims: cint; dims: ptr cint; 
                   periods: ptr cint; coords: ptr cint): cint {.
    importc: "MPI_Cart_get", header: "mpi.h".}
proc MPI_Cart_map*(comm: MPI_Comm; ndims: cint; dims: ptr cint; 
                   periods: ptr cint; newrank: ptr cint): cint {.
    importc: "MPI_Cart_map", header: "mpi.h".}
proc MPI_Cart_rank*(comm: MPI_Comm; coords: ptr cint; rank: ptr cint): cint {.
    importc: "MPI_Cart_rank", header: "mpi.h".}
proc MPI_Cart_shift*(comm: MPI_Comm; direction: cint; disp: cint; 
                     rank_source: ptr cint; rank_dest: ptr cint): cint {.
    importc: "MPI_Cart_shift", header: "mpi.h".}
proc MPI_Cart_sub*(comm: MPI_Comm; remain_dims: ptr cint; new_comm: ptr MPI_Comm): cint {.
    importc: "MPI_Cart_sub", header: "mpi.h".}
proc MPI_Cartdim_get*(comm: MPI_Comm; ndims: ptr cint): cint {.
    importc: "MPI_Cartdim_get", header: "mpi.h".}
proc MPI_Close_port*(port_name: cstring): cint {.importc: "MPI_Close_port", 
    header: "mpi.h".}
proc MPI_Comm_accept*(port_name: cstring; info: MPI_Info; root: cint; 
                      comm: MPI_Comm; newcomm: ptr MPI_Comm): cint {.
    importc: "MPI_Comm_accept", header: "mpi.h".}
proc MPI_Comm_call_errhandler*(comm: MPI_Comm; errorcode: cint): cint {.
    importc: "MPI_Comm_call_errhandler", header: "mpi.h".}
proc MPI_Comm_compare*(comm1: MPI_Comm; comm2: MPI_Comm; result: ptr cint): cint {.
    importc: "MPI_Comm_compare", header: "mpi.h".}
proc MPI_Comm_connect*(port_name: cstring; info: MPI_Info; root: cint; 
                       comm: MPI_Comm; newcomm: ptr MPI_Comm): cint {.
    importc: "MPI_Comm_connect", header: "mpi.h".}
proc MPI_Comm_create_errhandler*(function: ptr MPI_Comm_errhandler_function; 
                                 errhandler: ptr MPI_Errhandler): cint {.
    importc: "MPI_Comm_create_errhandler", header: "mpi.h".}
proc MPI_Comm_create_keyval*(comm_copy_attr_fn: ptr MPI_Comm_copy_attr_function; 
    comm_delete_attr_fn: ptr MPI_Comm_delete_attr_function; 
                             comm_keyval: ptr cint; extra_state: pointer): cint {.
    importc: "MPI_Comm_create_keyval", header: "mpi.h".}
proc MPI_Comm_create_group*(comm: MPI_Comm; group: MPI_Group; tag: cint; 
                            newcomm: ptr MPI_Comm): cint {.
    importc: "MPI_Comm_create_group", header: "mpi.h".}
proc MPI_Comm_create*(comm: MPI_Comm; group: MPI_Group; newcomm: ptr MPI_Comm): cint {.
    importc: "MPI_Comm_create", header: "mpi.h".}
proc MPI_Comm_delete_attr*(comm: MPI_Comm; comm_keyval: cint): cint {.
    importc: "MPI_Comm_delete_attr", header: "mpi.h".}
proc MPI_Comm_disconnect*(comm: ptr MPI_Comm): cint {.
    importc: "MPI_Comm_disconnect", header: "mpi.h".}
proc MPI_Comm_dup*(comm: MPI_Comm; newcomm: ptr MPI_Comm): cint {.
    importc: "MPI_Comm_dup", header: "mpi.h".}
proc MPI_Comm_idup*(comm: MPI_Comm; newcomm: ptr MPI_Comm; 
                    request: ptr MPI_Request): cint {.importc: "MPI_Comm_idup", 
    header: "mpi.h".}
proc MPI_Comm_dup_with_info*(comm: MPI_Comm; info: MPI_Info; 
                             newcomm: ptr MPI_Comm): cint {.
    importc: "MPI_Comm_dup_with_info", header: "mpi.h".}
proc MPI_Comm_free_keyval*(comm_keyval: ptr cint): cint {.
    importc: "MPI_Comm_free_keyval", header: "mpi.h".}
proc MPI_Comm_free*(comm: ptr MPI_Comm): cint {.importc: "MPI_Comm_free", 
    header: "mpi.h".}
proc MPI_Comm_get_attr*(comm: MPI_Comm; comm_keyval: cint; 
                        attribute_val: pointer; flag: ptr cint): cint {.
    importc: "MPI_Comm_get_attr", header: "mpi.h".}
proc MPI_Dist_graph_create*(comm_old: MPI_Comm; n: cint; nodes: ptr cint; 
                            degrees: ptr cint; targets: ptr cint; 
                            weights: ptr cint; info: MPI_Info; reorder: cint; 
                            newcomm: ptr MPI_Comm): cint {.
    importc: "MPI_Dist_graph_create", header: "mpi.h".}
proc MPI_Dist_graph_create_adjacent*(comm_old: MPI_Comm; indegree: cint; 
                                     sources: ptr cint; sourceweights: ptr cint; 
                                     outdegree: cint; destinations: ptr cint; 
                                     destweights: ptr cint; info: MPI_Info; 
                                     reorder: cint; 
                                     comm_dist_graph: ptr MPI_Comm): cint {.
    importc: "MPI_Dist_graph_create_adjacent", header: "mpi.h".}
proc MPI_Dist_graph_neighbors*(comm: MPI_Comm; maxindegree: cint; 
                               sources: ptr cint; sourceweights: ptr cint; 
                               maxoutdegree: cint; destinations: ptr cint; 
                               destweights: ptr cint): cint {.
    importc: "MPI_Dist_graph_neighbors", header: "mpi.h".}
proc MPI_Dist_graph_neighbors_count*(comm: MPI_Comm; inneighbors: ptr cint; 
                                     outneighbors: ptr cint; weighted: ptr cint): cint {.
    importc: "MPI_Dist_graph_neighbors_count", header: "mpi.h".}
proc MPI_Comm_get_errhandler*(comm: MPI_Comm; erhandler: ptr MPI_Errhandler): cint {.
    importc: "MPI_Comm_get_errhandler", header: "mpi.h".}
proc MPI_Comm_get_info*(comm: MPI_Comm; info_used: ptr MPI_Info): cint {.
    importc: "MPI_Comm_get_info", header: "mpi.h".}
proc MPI_Comm_get_name*(comm: MPI_Comm; comm_name: cstring; resultlen: ptr cint): cint {.
    importc: "MPI_Comm_get_name", header: "mpi.h".}
proc MPI_Comm_get_parent*(parent: ptr MPI_Comm): cint {.
    importc: "MPI_Comm_get_parent", header: "mpi.h".}
proc MPI_Comm_group*(comm: MPI_Comm; group: ptr MPI_Group): cint {.
    importc: "MPI_Comm_group", header: "mpi.h".}
proc MPI_Comm_join*(fd: cint; intercomm: ptr MPI_Comm): cint {.
    importc: "MPI_Comm_join", header: "mpi.h".}
proc MPI_Comm_rank*(comm: MPI_Comm; rank: ptr cint): cint {.
    importc: "MPI_Comm_rank", header: "mpi.h".}
proc MPI_Comm_remote_group*(comm: MPI_Comm; group: ptr MPI_Group): cint {.
    importc: "MPI_Comm_remote_group", header: "mpi.h".}
proc MPI_Comm_remote_size*(comm: MPI_Comm; size: ptr cint): cint {.
    importc: "MPI_Comm_remote_size", header: "mpi.h".}
proc MPI_Comm_set_attr*(comm: MPI_Comm; comm_keyval: cint; 
                        attribute_val: pointer): cint {.
    importc: "MPI_Comm_set_attr", header: "mpi.h".}
proc MPI_Comm_set_errhandler*(comm: MPI_Comm; errhandler: MPI_Errhandler): cint {.
    importc: "MPI_Comm_set_errhandler", header: "mpi.h".}
proc MPI_Comm_set_info*(comm: MPI_Comm; info: MPI_Info): cint {.
    importc: "MPI_Comm_set_info", header: "mpi.h".}
proc MPI_Comm_set_name*(comm: MPI_Comm; comm_name: cstring): cint {.
    importc: "MPI_Comm_set_name", header: "mpi.h".}
proc MPI_Comm_size*(comm: MPI_Comm; size: ptr cint): cint {.
    importc: "MPI_Comm_size", header: "mpi.h".}
proc MPI_Comm_spawn*(command: cstring; argv: ptr cstring; maxprocs: cint; 
                     info: MPI_Info; root: cint; comm: MPI_Comm; 
                     intercomm: ptr MPI_Comm; array_of_errcodes: ptr cint): cint {.
    importc: "MPI_Comm_spawn", header: "mpi.h".}
proc MPI_Comm_spawn_multiple*(count: cint; array_of_commands: ptr cstring; 
                              array_of_argv: ptr cstringArray; 
                              array_of_maxprocs: ptr cint; 
                              array_of_info: ptr MPI_Info; root: cint; 
                              comm: MPI_Comm; intercomm: ptr MPI_Comm; 
                              array_of_errcodes: ptr cint): cint {.
    importc: "MPI_Comm_spawn_multiple", header: "mpi.h".}
proc MPI_Comm_split*(comm: MPI_Comm; color: cint; key: cint; 
                     newcomm: ptr MPI_Comm): cint {.importc: "MPI_Comm_split", 
    header: "mpi.h".}
proc MPI_Comm_split_type*(comm: MPI_Comm; split_type: cint; key: cint; 
                          info: MPI_Info; newcomm: ptr MPI_Comm): cint {.
    importc: "MPI_Comm_split_type", header: "mpi.h".}
proc MPI_Comm_test_inter*(comm: MPI_Comm; flag: ptr cint): cint {.
    importc: "MPI_Comm_test_inter", header: "mpi.h".}
proc MPI_Compare_and_swap*(origin_addr: pointer; compare_addr: pointer; 
                           result_addr: pointer; datatype: MPI_Datatype; 
                           target_rank: cint; target_disp: MPI_Aint; 
                           win: MPI_Win): cint {.
    importc: "MPI_Compare_and_swap", header: "mpi.h".}
proc MPI_Dims_create*(nnodes: cint; ndims: cint; dims: ptr cint): cint {.
    importc: "MPI_Dims_create", header: "mpi.h".}
proc MPI_Errhandler_create*(function: ptr MPI_Handler_function; 
                            errhandler: ptr MPI_Errhandler): cint {.
    importc: "MPI_Errhandler_create", header: "mpi.h".}
proc MPI_Errhandler_free*(errhandler: ptr MPI_Errhandler): cint {.
    importc: "MPI_Errhandler_free", header: "mpi.h".}
proc MPI_Errhandler_get*(comm: MPI_Comm; errhandler: ptr MPI_Errhandler): cint {.
    importc: "MPI_Errhandler_get", header: "mpi.h".}
proc MPI_Errhandler_set*(comm: MPI_Comm; errhandler: MPI_Errhandler): cint {.
    importc: "MPI_Errhandler_set", header: "mpi.h".}
proc MPI_Error_class*(errorcode: cint; errorclass: ptr cint): cint {.
    importc: "MPI_Error_class", header: "mpi.h".}
proc MPI_Error_string*(errorcode: cint; string: cstring; resultlen: ptr cint): cint {.
    importc: "MPI_Error_string", header: "mpi.h".}
proc MPI_Exscan*(sendbuf: pointer; recvbuf: pointer; count: cint; 
                 datatype: MPI_Datatype; op: MPI_Op; comm: MPI_Comm): cint {.
    importc: "MPI_Exscan", header: "mpi.h".}
proc MPI_Fetch_and_op*(origin_addr: pointer; result_addr: pointer; 
                       datatype: MPI_Datatype; target_rank: cint; 
                       target_disp: MPI_Aint; op: MPI_Op; win: MPI_Win): cint {.
    importc: "MPI_Fetch_and_op", header: "mpi.h".}
proc MPI_Iexscan*(sendbuf: pointer; recvbuf: pointer; count: cint; 
                  datatype: MPI_Datatype; op: MPI_Op; comm: MPI_Comm; 
                  request: ptr MPI_Request): cint {.importc: "MPI_Iexscan", 
    header: "mpi.h".}
proc MPI_File_call_errhandler*(fh: MPI_File; errorcode: cint): cint {.
    importc: "MPI_File_call_errhandler", header: "mpi.h".}
proc MPI_File_create_errhandler*(function: ptr MPI_File_errhandler_function; 
                                 errhandler: ptr MPI_Errhandler): cint {.
    importc: "MPI_File_create_errhandler", header: "mpi.h".}
proc MPI_File_set_errhandler*(file: MPI_File; errhandler: MPI_Errhandler): cint {.
    importc: "MPI_File_set_errhandler", header: "mpi.h".}
proc MPI_File_get_errhandler*(file: MPI_File; errhandler: ptr MPI_Errhandler): cint {.
    importc: "MPI_File_get_errhandler", header: "mpi.h".}
proc MPI_File_open*(comm: MPI_Comm; filename: cstring; amode: cint; 
                    info: MPI_Info; fh: ptr MPI_File): cint {.
    importc: "MPI_File_open", header: "mpi.h".}
proc MPI_File_close*(fh: ptr MPI_File): cint {.importc: "MPI_File_close", 
    header: "mpi.h".}
proc MPI_File_delete*(filename: cstring; info: MPI_Info): cint {.
    importc: "MPI_File_delete", header: "mpi.h".}
proc MPI_File_set_size*(fh: MPI_File; size: MPI_Offset): cint {.
    importc: "MPI_File_set_size", header: "mpi.h".}
proc MPI_File_preallocate*(fh: MPI_File; size: MPI_Offset): cint {.
    importc: "MPI_File_preallocate", header: "mpi.h".}
proc MPI_File_get_size*(fh: MPI_File; size: ptr MPI_Offset): cint {.
    importc: "MPI_File_get_size", header: "mpi.h".}
proc MPI_File_get_group*(fh: MPI_File; group: ptr MPI_Group): cint {.
    importc: "MPI_File_get_group", header: "mpi.h".}
proc MPI_File_get_amode*(fh: MPI_File; amode: ptr cint): cint {.
    importc: "MPI_File_get_amode", header: "mpi.h".}
proc MPI_File_set_info*(fh: MPI_File; info: MPI_Info): cint {.
    importc: "MPI_File_set_info", header: "mpi.h".}
proc MPI_File_get_info*(fh: MPI_File; info_used: ptr MPI_Info): cint {.
    importc: "MPI_File_get_info", header: "mpi.h".}
proc MPI_File_set_view*(fh: MPI_File; disp: MPI_Offset; etype: MPI_Datatype; 
                        filetype: MPI_Datatype; datarep: cstring; info: MPI_Info): cint {.
    importc: "MPI_File_set_view", header: "mpi.h".}
proc MPI_File_get_view*(fh: MPI_File; disp: ptr MPI_Offset; 
                        etype: ptr MPI_Datatype; filetype: ptr MPI_Datatype; 
                        datarep: cstring): cint {.importc: "MPI_File_get_view", 
    header: "mpi.h".}
proc MPI_File_read_at*(fh: MPI_File; offset: MPI_Offset; buf: pointer; 
                       count: cint; datatype: MPI_Datatype; 
                       status: ptr MPI_Status): cint {.
    importc: "MPI_File_read_at", header: "mpi.h".}
proc MPI_File_read_at_all*(fh: MPI_File; offset: MPI_Offset; buf: pointer; 
                           count: cint; datatype: MPI_Datatype; 
                           status: ptr MPI_Status): cint {.
    importc: "MPI_File_read_at_all", header: "mpi.h".}
proc MPI_File_write_at*(fh: MPI_File; offset: MPI_Offset; buf: pointer; 
                        count: cint; datatype: MPI_Datatype; 
                        status: ptr MPI_Status): cint {.
    importc: "MPI_File_write_at", header: "mpi.h".}
proc MPI_File_write_at_all*(fh: MPI_File; offset: MPI_Offset; buf: pointer; 
                            count: cint; datatype: MPI_Datatype; 
                            status: ptr MPI_Status): cint {.
    importc: "MPI_File_write_at_all", header: "mpi.h".}
proc MPI_File_iread_at*(fh: MPI_File; offset: MPI_Offset; buf: pointer; 
                        count: cint; datatype: MPI_Datatype; 
                        request: ptr MPI_Request): cint {.
    importc: "MPI_File_iread_at", header: "mpi.h".}
proc MPI_File_iwrite_at*(fh: MPI_File; offset: MPI_Offset; buf: pointer; 
                         count: cint; datatype: MPI_Datatype; 
                         request: ptr MPI_Request): cint {.
    importc: "MPI_File_iwrite_at", header: "mpi.h".}
proc MPI_File_iread_at_all*(fh: MPI_File; offset: MPI_Offset; buf: pointer; 
                            count: cint; datatype: MPI_Datatype; 
                            request: ptr MPI_Request): cint {.
    importc: "MPI_File_iread_at_all", header: "mpi.h".}
proc MPI_File_iwrite_at_all*(fh: MPI_File; offset: MPI_Offset; buf: pointer; 
                             count: cint; datatype: MPI_Datatype; 
                             request: ptr MPI_Request): cint {.
    importc: "MPI_File_iwrite_at_all", header: "mpi.h".}
proc MPI_File_read*(fh: MPI_File; buf: pointer; count: cint; 
                    datatype: MPI_Datatype; status: ptr MPI_Status): cint {.
    importc: "MPI_File_read", header: "mpi.h".}
proc MPI_File_read_all*(fh: MPI_File; buf: pointer; count: cint; 
                        datatype: MPI_Datatype; status: ptr MPI_Status): cint {.
    importc: "MPI_File_read_all", header: "mpi.h".}
proc MPI_File_write*(fh: MPI_File; buf: pointer; count: cint; 
                     datatype: MPI_Datatype; status: ptr MPI_Status): cint {.
    importc: "MPI_File_write", header: "mpi.h".}
proc MPI_File_write_all*(fh: MPI_File; buf: pointer; count: cint; 
                         datatype: MPI_Datatype; status: ptr MPI_Status): cint {.
    importc: "MPI_File_write_all", header: "mpi.h".}
proc MPI_File_iread*(fh: MPI_File; buf: pointer; count: cint; 
                     datatype: MPI_Datatype; request: ptr MPI_Request): cint {.
    importc: "MPI_File_iread", header: "mpi.h".}
proc MPI_File_iwrite*(fh: MPI_File; buf: pointer; count: cint; 
                      datatype: MPI_Datatype; request: ptr MPI_Request): cint {.
    importc: "MPI_File_iwrite", header: "mpi.h".}
proc MPI_File_iread_all*(fh: MPI_File; buf: pointer; count: cint; 
                         datatype: MPI_Datatype; request: ptr MPI_Request): cint {.
    importc: "MPI_File_iread_all", header: "mpi.h".}
proc MPI_File_iwrite_all*(fh: MPI_File; buf: pointer; count: cint; 
                          datatype: MPI_Datatype; request: ptr MPI_Request): cint {.
    importc: "MPI_File_iwrite_all", header: "mpi.h".}
proc MPI_File_seek*(fh: MPI_File; offset: MPI_Offset; whence: cint): cint {.
    importc: "MPI_File_seek", header: "mpi.h".}
proc MPI_File_get_position*(fh: MPI_File; offset: ptr MPI_Offset): cint {.
    importc: "MPI_File_get_position", header: "mpi.h".}
proc MPI_File_get_byte_offset*(fh: MPI_File; offset: MPI_Offset; 
                               disp: ptr MPI_Offset): cint {.
    importc: "MPI_File_get_byte_offset", header: "mpi.h".}
proc MPI_File_read_shared*(fh: MPI_File; buf: pointer; count: cint; 
                           datatype: MPI_Datatype; status: ptr MPI_Status): cint {.
    importc: "MPI_File_read_shared", header: "mpi.h".}
proc MPI_File_write_shared*(fh: MPI_File; buf: pointer; count: cint; 
                            datatype: MPI_Datatype; status: ptr MPI_Status): cint {.
    importc: "MPI_File_write_shared", header: "mpi.h".}
proc MPI_File_iread_shared*(fh: MPI_File; buf: pointer; count: cint; 
                            datatype: MPI_Datatype; request: ptr MPI_Request): cint {.
    importc: "MPI_File_iread_shared", header: "mpi.h".}
proc MPI_File_iwrite_shared*(fh: MPI_File; buf: pointer; count: cint; 
                             datatype: MPI_Datatype; request: ptr MPI_Request): cint {.
    importc: "MPI_File_iwrite_shared", header: "mpi.h".}
proc MPI_File_read_ordered*(fh: MPI_File; buf: pointer; count: cint; 
                            datatype: MPI_Datatype; status: ptr MPI_Status): cint {.
    importc: "MPI_File_read_ordered", header: "mpi.h".}
proc MPI_File_write_ordered*(fh: MPI_File; buf: pointer; count: cint; 
                             datatype: MPI_Datatype; status: ptr MPI_Status): cint {.
    importc: "MPI_File_write_ordered", header: "mpi.h".}
proc MPI_File_seek_shared*(fh: MPI_File; offset: MPI_Offset; whence: cint): cint {.
    importc: "MPI_File_seek_shared", header: "mpi.h".}
proc MPI_File_get_position_shared*(fh: MPI_File; offset: ptr MPI_Offset): cint {.
    importc: "MPI_File_get_position_shared", header: "mpi.h".}
proc MPI_File_read_at_all_begin*(fh: MPI_File; offset: MPI_Offset; buf: pointer; 
                                 count: cint; datatype: MPI_Datatype): cint {.
    importc: "MPI_File_read_at_all_begin", header: "mpi.h".}
proc MPI_File_read_at_all_end*(fh: MPI_File; buf: pointer; 
                               status: ptr MPI_Status): cint {.
    importc: "MPI_File_read_at_all_end", header: "mpi.h".}
proc MPI_File_write_at_all_begin*(fh: MPI_File; offset: MPI_Offset; 
                                  buf: pointer; count: cint; 
                                  datatype: MPI_Datatype): cint {.
    importc: "MPI_File_write_at_all_begin", header: "mpi.h".}
proc MPI_File_write_at_all_end*(fh: MPI_File; buf: pointer; 
                                status: ptr MPI_Status): cint {.
    importc: "MPI_File_write_at_all_end", header: "mpi.h".}
proc MPI_File_read_all_begin*(fh: MPI_File; buf: pointer; count: cint; 
                              datatype: MPI_Datatype): cint {.
    importc: "MPI_File_read_all_begin", header: "mpi.h".}
proc MPI_File_read_all_end*(fh: MPI_File; buf: pointer; status: ptr MPI_Status): cint {.
    importc: "MPI_File_read_all_end", header: "mpi.h".}
proc MPI_File_write_all_begin*(fh: MPI_File; buf: pointer; count: cint; 
                               datatype: MPI_Datatype): cint {.
    importc: "MPI_File_write_all_begin", header: "mpi.h".}
proc MPI_File_write_all_end*(fh: MPI_File; buf: pointer; status: ptr MPI_Status): cint {.
    importc: "MPI_File_write_all_end", header: "mpi.h".}
proc MPI_File_read_ordered_begin*(fh: MPI_File; buf: pointer; count: cint; 
                                  datatype: MPI_Datatype): cint {.
    importc: "MPI_File_read_ordered_begin", header: "mpi.h".}
proc MPI_File_read_ordered_end*(fh: MPI_File; buf: pointer; 
                                status: ptr MPI_Status): cint {.
    importc: "MPI_File_read_ordered_end", header: "mpi.h".}
proc MPI_File_write_ordered_begin*(fh: MPI_File; buf: pointer; count: cint; 
                                   datatype: MPI_Datatype): cint {.
    importc: "MPI_File_write_ordered_begin", header: "mpi.h".}
proc MPI_File_write_ordered_end*(fh: MPI_File; buf: pointer; 
                                 status: ptr MPI_Status): cint {.
    importc: "MPI_File_write_ordered_end", header: "mpi.h".}
proc MPI_File_get_type_extent*(fh: MPI_File; datatype: MPI_Datatype; 
                               extent: ptr MPI_Aint): cint {.
    importc: "MPI_File_get_type_extent", header: "mpi.h".}
proc MPI_File_set_atomicity*(fh: MPI_File; flag: cint): cint {.
    importc: "MPI_File_set_atomicity", header: "mpi.h".}
proc MPI_File_get_atomicity*(fh: MPI_File; flag: ptr cint): cint {.
    importc: "MPI_File_get_atomicity", header: "mpi.h".}
proc MPI_File_sync*(fh: MPI_File): cint {.importc: "MPI_File_sync", 
    header: "mpi.h".}
proc MPI_Finalize*(): cint {.importc: "MPI_Finalize", header: "mpi.h".}
proc MPI_Finalized*(flag: ptr cint): cint {.importc: "MPI_Finalized", 
    header: "mpi.h".}
proc MPI_Free_mem*(base: pointer): cint {.importc: "MPI_Free_mem", 
    header: "mpi.h".}
proc MPI_Gather*(sendbuf: pointer; sendcount: cint; sendtype: MPI_Datatype; 
                 recvbuf: pointer; recvcount: cint; recvtype: MPI_Datatype; 
                 root: cint; comm: MPI_Comm): cint {.importc: "MPI_Gather", 
    header: "mpi.h".}
proc MPI_Igather*(sendbuf: pointer; sendcount: cint; sendtype: MPI_Datatype; 
                  recvbuf: pointer; recvcount: cint; recvtype: MPI_Datatype; 
                  root: cint; comm: MPI_Comm; request: ptr MPI_Request): cint {.
    importc: "MPI_Igather", header: "mpi.h".}
proc MPI_Gatherv*(sendbuf: pointer; sendcount: cint; sendtype: MPI_Datatype; 
                  recvbuf: pointer; recvcounts: ptr cint; displs: ptr cint; 
                  recvtype: MPI_Datatype; root: cint; comm: MPI_Comm): cint {.
    importc: "MPI_Gatherv", header: "mpi.h".}
proc MPI_Igatherv*(sendbuf: pointer; sendcount: cint; sendtype: MPI_Datatype; 
                   recvbuf: pointer; recvcounts: ptr cint; displs: ptr cint; 
                   recvtype: MPI_Datatype; root: cint; comm: MPI_Comm; 
                   request: ptr MPI_Request): cint {.importc: "MPI_Igatherv", 
    header: "mpi.h".}
proc MPI_Get_address*(location: pointer; address: ptr MPI_Aint): cint {.
    importc: "MPI_Get_address", header: "mpi.h".}
proc MPI_Get_count*(status: ptr MPI_Status; datatype: MPI_Datatype; 
                    count: ptr cint): cint {.importc: "MPI_Get_count", 
    header: "mpi.h".}
proc MPI_Get_elements*(status: ptr MPI_Status; datatype: MPI_Datatype; 
                       count: ptr cint): cint {.importc: "MPI_Get_elements", 
    header: "mpi.h".}
proc MPI_Get_elements_x*(status: ptr MPI_Status; datatype: MPI_Datatype; 
                         count: ptr MPI_Count): cint {.
    importc: "MPI_Get_elements_x", header: "mpi.h".}
proc MPI_Get*(origin_addr: pointer; origin_count: cint; 
              origin_datatype: MPI_Datatype; target_rank: cint; 
              target_disp: MPI_Aint; target_count: cint; 
              target_datatype: MPI_Datatype; win: MPI_Win): cint {.
    importc: "MPI_Get", header: "mpi.h".}
proc MPI_Get_accumulate*(origin_addr: pointer; origin_count: cint; 
                         origin_datatype: MPI_Datatype; result_addr: pointer; 
                         result_count: cint; result_datatype: MPI_Datatype; 
                         target_rank: cint; target_disp: MPI_Aint; 
                         target_count: cint; target_datatype: MPI_Datatype; 
                         op: MPI_Op; win: MPI_Win): cint {.
    importc: "MPI_Get_accumulate", header: "mpi.h".}
proc MPI_Get_library_version*(version: cstring; resultlen: ptr cint): cint {.
    importc: "MPI_Get_library_version", header: "mpi.h".}
proc MPI_Get_processor_name*(name: cstring; resultlen: ptr cint): cint {.
    importc: "MPI_Get_processor_name", header: "mpi.h".}
proc MPI_Get_version*(version: ptr cint; subversion: ptr cint): cint {.
    importc: "MPI_Get_version", header: "mpi.h".}
proc MPI_Graph_create*(comm_old: MPI_Comm; nnodes: cint; index: ptr cint; 
                       edges: ptr cint; reorder: cint; comm_graph: ptr MPI_Comm): cint {.
    importc: "MPI_Graph_create", header: "mpi.h".}
proc MPI_Graph_get*(comm: MPI_Comm; maxindex: cint; maxedges: cint; 
                    index: ptr cint; edges: ptr cint): cint {.
    importc: "MPI_Graph_get", header: "mpi.h".}
proc MPI_Graph_map*(comm: MPI_Comm; nnodes: cint; index: ptr cint; 
                    edges: ptr cint; newrank: ptr cint): cint {.
    importc: "MPI_Graph_map", header: "mpi.h".}
proc MPI_Graph_neighbors_count*(comm: MPI_Comm; rank: cint; nneighbors: ptr cint): cint {.
    importc: "MPI_Graph_neighbors_count", header: "mpi.h".}
proc MPI_Graph_neighbors*(comm: MPI_Comm; rank: cint; maxneighbors: cint; 
                          neighbors: ptr cint): cint {.
    importc: "MPI_Graph_neighbors", header: "mpi.h".}
proc MPI_Graphdims_get*(comm: MPI_Comm; nnodes: ptr cint; nedges: ptr cint): cint {.
    importc: "MPI_Graphdims_get", header: "mpi.h".}
proc MPI_Grequest_complete*(request: MPI_Request): cint {.
    importc: "MPI_Grequest_complete", header: "mpi.h".}
proc MPI_Grequest_start*(query_fn: ptr MPI_Grequest_query_function; 
                         free_fn: ptr MPI_Grequest_free_function; 
                         cancel_fn: ptr MPI_Grequest_cancel_function; 
                         extra_state: pointer; request: ptr MPI_Request): cint {.
    importc: "MPI_Grequest_start", header: "mpi.h".}
proc MPI_Group_compare*(group1: MPI_Group; group2: MPI_Group; result: ptr cint): cint {.
    importc: "MPI_Group_compare", header: "mpi.h".}
proc MPI_Group_difference*(group1: MPI_Group; group2: MPI_Group; 
                           newgroup: ptr MPI_Group): cint {.
    importc: "MPI_Group_difference", header: "mpi.h".}
proc MPI_Group_excl*(group: MPI_Group; n: cint; ranks: ptr cint; 
                     newgroup: ptr MPI_Group): cint {.importc: "MPI_Group_excl", 
    header: "mpi.h".}
proc MPI_Group_free*(group: ptr MPI_Group): cint {.importc: "MPI_Group_free", 
    header: "mpi.h".}
proc MPI_Group_incl*(group: MPI_Group; n: cint; ranks: ptr cint; 
                     newgroup: ptr MPI_Group): cint {.importc: "MPI_Group_incl", 
    header: "mpi.h".}
proc MPI_Group_intersection*(group1: MPI_Group; group2: MPI_Group; 
                             newgroup: ptr MPI_Group): cint {.
    importc: "MPI_Group_intersection", header: "mpi.h".}
proc MPI_Group_range_excl*(group: MPI_Group; n: cint; 
                           ranges: ptr array[3, cint]; newgroup: ptr MPI_Group): cint {.
    importc: "MPI_Group_range_excl", header: "mpi.h".}
proc MPI_Group_range_incl*(group: MPI_Group; n: cint; 
                           ranges: ptr array[3, cint]; newgroup: ptr MPI_Group): cint {.
    importc: "MPI_Group_range_incl", header: "mpi.h".}
proc MPI_Group_rank*(group: MPI_Group; rank: ptr cint): cint {.
    importc: "MPI_Group_rank", header: "mpi.h".}
proc MPI_Group_size*(group: MPI_Group; size: ptr cint): cint {.
    importc: "MPI_Group_size", header: "mpi.h".}
proc MPI_Group_translate_ranks*(group1: MPI_Group; n: cint; ranks1: ptr cint; 
                                group2: MPI_Group; ranks2: ptr cint): cint {.
    importc: "MPI_Group_translate_ranks", header: "mpi.h".}
proc MPI_Group_union*(group1: MPI_Group; group2: MPI_Group; 
                      newgroup: ptr MPI_Group): cint {.
    importc: "MPI_Group_union", header: "mpi.h".}
proc MPI_Ibsend*(buf: pointer; count: cint; datatype: MPI_Datatype; dest: cint; 
                 tag: cint; comm: MPI_Comm; request: ptr MPI_Request): cint {.
    importc: "MPI_Ibsend", header: "mpi.h".}
proc MPI_Improbe*(source: cint; tag: cint; comm: MPI_Comm; flag: ptr cint; 
                  message: ptr MPI_Message; status: ptr MPI_Status): cint {.
    importc: "MPI_Improbe", header: "mpi.h".}
proc MPI_Imrecv*(buf: pointer; count: cint; `type`: MPI_Datatype; 
                 message: ptr MPI_Message; request: ptr MPI_Request): cint {.
    importc: "MPI_Imrecv", header: "mpi.h".}
proc MPI_Info_create*(info: ptr MPI_Info): cint {.importc: "MPI_Info_create", 
    header: "mpi.h".}
proc MPI_Info_delete*(info: MPI_Info; key: cstring): cint {.
    importc: "MPI_Info_delete", header: "mpi.h".}
proc MPI_Info_dup*(info: MPI_Info; newinfo: ptr MPI_Info): cint {.
    importc: "MPI_Info_dup", header: "mpi.h".}
proc MPI_Info_free*(info: ptr MPI_Info): cint {.importc: "MPI_Info_free", 
    header: "mpi.h".}
proc MPI_Info_get*(info: MPI_Info; key: cstring; valuelen: cint; value: cstring; 
                   flag: ptr cint): cint {.importc: "MPI_Info_get", 
    header: "mpi.h".}
proc MPI_Info_get_nkeys*(info: MPI_Info; nkeys: ptr cint): cint {.
    importc: "MPI_Info_get_nkeys", header: "mpi.h".}
proc MPI_Info_get_nthkey*(info: MPI_Info; n: cint; key: cstring): cint {.
    importc: "MPI_Info_get_nthkey", header: "mpi.h".}
proc MPI_Info_get_valuelen*(info: MPI_Info; key: cstring; valuelen: ptr cint; 
                            flag: ptr cint): cint {.
    importc: "MPI_Info_get_valuelen", header: "mpi.h".}
proc MPI_Info_set*(info: MPI_Info; key: cstring; value: cstring): cint {.
    importc: "MPI_Info_set", header: "mpi.h".}
proc MPI_Init*(argc: ptr cint; argv: ptr cstringArray): cint {.
    importc: "MPI_Init", header: "mpi.h".}
proc MPI_Initialized*(flag: ptr cint): cint {.importc: "MPI_Initialized", 
    header: "mpi.h".}
proc MPI_Init_thread*(argc: ptr cint; argv: ptr cstringArray; required: cint; 
                      provided: ptr cint): cint {.importc: "MPI_Init_thread", 
    header: "mpi.h".}
proc MPI_Intercomm_create*(local_comm: MPI_Comm; local_leader: cint; 
                           bridge_comm: MPI_Comm; remote_leader: cint; 
                           tag: cint; newintercomm: ptr MPI_Comm): cint {.
    importc: "MPI_Intercomm_create", header: "mpi.h".}
proc MPI_Intercomm_merge*(intercomm: MPI_Comm; high: cint; 
                          newintercomm: ptr MPI_Comm): cint {.
    importc: "MPI_Intercomm_merge", header: "mpi.h".}
proc MPI_Iprobe*(source: cint; tag: cint; comm: MPI_Comm; flag: ptr cint; 
                 status: ptr MPI_Status): cint {.importc: "MPI_Iprobe", 
    header: "mpi.h".}
proc MPI_Irecv*(buf: pointer; count: cint; datatype: MPI_Datatype; source: cint; 
                tag: cint; comm: MPI_Comm; request: ptr MPI_Request): cint {.
    importc: "MPI_Irecv", header: "mpi.h".}
proc MPI_Irsend*(buf: pointer; count: cint; datatype: MPI_Datatype; dest: cint; 
                 tag: cint; comm: MPI_Comm; request: ptr MPI_Request): cint {.
    importc: "MPI_Irsend", header: "mpi.h".}
proc MPI_Isend*(buf: pointer; count: cint; datatype: MPI_Datatype; dest: cint; 
                tag: cint; comm: MPI_Comm; request: ptr MPI_Request): cint {.
    importc: "MPI_Isend", header: "mpi.h".}
proc MPI_Issend*(buf: pointer; count: cint; datatype: MPI_Datatype; dest: cint; 
                 tag: cint; comm: MPI_Comm; request: ptr MPI_Request): cint {.
    importc: "MPI_Issend", header: "mpi.h".}
proc MPI_Is_thread_main*(flag: ptr cint): cint {.importc: "MPI_Is_thread_main", 
    header: "mpi.h".}
proc MPI_Keyval_create*(copy_fn: ptr MPI_Copy_function; 
                        delete_fn: ptr MPI_Delete_function; keyval: ptr cint; 
                        extra_state: pointer): cint {.
    importc: "MPI_Keyval_create", header: "mpi.h".}
proc MPI_Keyval_free*(keyval: ptr cint): cint {.importc: "MPI_Keyval_free", 
    header: "mpi.h".}
proc MPI_Lookup_name*(service_name: cstring; info: MPI_Info; port_name: cstring): cint {.
    importc: "MPI_Lookup_name", header: "mpi.h".}
proc MPI_Mprobe*(source: cint; tag: cint; comm: MPI_Comm; 
                 message: ptr MPI_Message; status: ptr MPI_Status): cint {.
    importc: "MPI_Mprobe", header: "mpi.h".}
proc MPI_Mrecv*(buf: pointer; count: cint; `type`: MPI_Datatype; 
                message: ptr MPI_Message; status: ptr MPI_Status): cint {.
    importc: "MPI_Mrecv", header: "mpi.h".}
proc MPI_Neighbor_allgather*(sendbuf: pointer; sendcount: cint; 
                             sendtype: MPI_Datatype; recvbuf: pointer; 
                             recvcount: cint; recvtype: MPI_Datatype; 
                             comm: MPI_Comm): cint {.
    importc: "MPI_Neighbor_allgather", header: "mpi.h".}
proc MPI_Ineighbor_allgather*(sendbuf: pointer; sendcount: cint; 
                              sendtype: MPI_Datatype; recvbuf: pointer; 
                              recvcount: cint; recvtype: MPI_Datatype; 
                              comm: MPI_Comm; request: ptr MPI_Request): cint {.
    importc: "MPI_Ineighbor_allgather", header: "mpi.h".}
proc MPI_Neighbor_allgatherv*(sendbuf: pointer; sendcount: cint; 
                              sendtype: MPI_Datatype; recvbuf: pointer; 
                              recvcounts: ptr cint; displs: ptr cint; 
                              recvtype: MPI_Datatype; comm: MPI_Comm): cint {.
    importc: "MPI_Neighbor_allgatherv", header: "mpi.h".}
proc MPI_Ineighbor_allgatherv*(sendbuf: pointer; sendcount: cint; 
                               sendtype: MPI_Datatype; recvbuf: pointer; 
                               recvcounts: ptr cint; displs: ptr cint; 
                               recvtype: MPI_Datatype; comm: MPI_Comm; 
                               request: ptr MPI_Request): cint {.
    importc: "MPI_Ineighbor_allgatherv", header: "mpi.h".}
proc MPI_Neighbor_alltoall*(sendbuf: pointer; sendcount: cint; 
                            sendtype: MPI_Datatype; recvbuf: pointer; 
                            recvcount: cint; recvtype: MPI_Datatype; 
                            comm: MPI_Comm): cint {.
    importc: "MPI_Neighbor_alltoall", header: "mpi.h".}
proc MPI_Ineighbor_alltoall*(sendbuf: pointer; sendcount: cint; 
                             sendtype: MPI_Datatype; recvbuf: pointer; 
                             recvcount: cint; recvtype: MPI_Datatype; 
                             comm: MPI_Comm; request: ptr MPI_Request): cint {.
    importc: "MPI_Ineighbor_alltoall", header: "mpi.h".}
proc MPI_Neighbor_alltoallv*(sendbuf: pointer; sendcounts: ptr cint; 
                             sdispls: ptr cint; sendtype: MPI_Datatype; 
                             recvbuf: pointer; recvcounts: ptr cint; 
                             rdispls: ptr cint; recvtype: MPI_Datatype; 
                             comm: MPI_Comm): cint {.
    importc: "MPI_Neighbor_alltoallv", header: "mpi.h".}
proc MPI_Ineighbor_alltoallv*(sendbuf: pointer; sendcounts: ptr cint; 
                              sdispls: ptr cint; sendtype: MPI_Datatype; 
                              recvbuf: pointer; recvcounts: ptr cint; 
                              rdispls: ptr cint; recvtype: MPI_Datatype; 
                              comm: MPI_Comm; request: ptr MPI_Request): cint {.
    importc: "MPI_Ineighbor_alltoallv", header: "mpi.h".}
proc MPI_Neighbor_alltoallw*(sendbuf: pointer; sendcounts: ptr cint; 
                             sdispls: ptr MPI_Aint; sendtypes: ptr MPI_Datatype; 
                             recvbuf: pointer; recvcounts: ptr cint; 
                             rdispls: ptr MPI_Aint; recvtypes: ptr MPI_Datatype; 
                             comm: MPI_Comm): cint {.
    importc: "MPI_Neighbor_alltoallw", header: "mpi.h".}
proc MPI_Ineighbor_alltoallw*(sendbuf: pointer; sendcounts: ptr cint; 
                              sdispls: ptr MPI_Aint; 
                              sendtypes: ptr MPI_Datatype; recvbuf: pointer; 
                              recvcounts: ptr cint; rdispls: ptr MPI_Aint; 
                              recvtypes: ptr MPI_Datatype; comm: MPI_Comm; 
                              request: ptr MPI_Request): cint {.
    importc: "MPI_Ineighbor_alltoallw", header: "mpi.h".}
proc MPI_Op_commutative*(op: MPI_Op; commute: ptr cint): cint {.
    importc: "MPI_Op_commutative", header: "mpi.h".}
proc MPI_Op_create*(function: ptr MPI_User_function; commute: cint; 
                    op: ptr MPI_Op): cint {.importc: "MPI_Op_create", 
    header: "mpi.h".}
proc MPI_Open_port*(info: MPI_Info; port_name: cstring): cint {.
    importc: "MPI_Open_port", header: "mpi.h".}
proc MPI_Op_free*(op: ptr MPI_Op): cint {.importc: "MPI_Op_free", 
    header: "mpi.h".}
proc MPI_Pack_external*(datarep: ptr char; inbuf: pointer; incount: cint; 
                        datatype: MPI_Datatype; outbuf: pointer; 
                        outsize: MPI_Aint; position: ptr MPI_Aint): cint {.
    importc: "MPI_Pack_external", header: "mpi.h".}
proc MPI_Pack_external_size*(datarep: ptr char; incount: cint; 
                             datatype: MPI_Datatype; size: ptr MPI_Aint): cint {.
    importc: "MPI_Pack_external_size", header: "mpi.h".}
proc MPI_Pack*(inbuf: pointer; incount: cint; datatype: MPI_Datatype; 
               outbuf: pointer; outsize: cint; position: ptr cint; 
               comm: MPI_Comm): cint {.importc: "MPI_Pack", header: "mpi.h".}
proc MPI_Pack_size*(incount: cint; datatype: MPI_Datatype; comm: MPI_Comm; 
                    size: ptr cint): cint {.importc: "MPI_Pack_size", 
    header: "mpi.h".}
proc MPI_Pcontrol*(level: cint): cint {.varargs, importc: "MPI_Pcontrol", 
                                        header: "mpi.h".}
proc MPI_Probe*(source: cint; tag: cint; comm: MPI_Comm; status: ptr MPI_Status): cint {.
    importc: "MPI_Probe", header: "mpi.h".}
proc MPI_Publish_name*(service_name: cstring; info: MPI_Info; port_name: cstring): cint {.
    importc: "MPI_Publish_name", header: "mpi.h".}
proc MPI_Put*(origin_addr: pointer; origin_count: cint; 
              origin_datatype: MPI_Datatype; target_rank: cint; 
              target_disp: MPI_Aint; target_count: cint; 
              target_datatype: MPI_Datatype; win: MPI_Win): cint {.
    importc: "MPI_Put", header: "mpi.h".}
proc MPI_Query_thread*(provided: ptr cint): cint {.importc: "MPI_Query_thread", 
    header: "mpi.h".}
proc MPI_Raccumulate*(origin_addr: pointer; origin_count: cint; 
                      origin_datatype: MPI_Datatype; target_rank: cint; 
                      target_disp: MPI_Aint; target_count: cint; 
                      target_datatype: MPI_Datatype; op: MPI_Op; win: MPI_Win; 
                      request: ptr MPI_Request): cint {.
    importc: "MPI_Raccumulate", header: "mpi.h".}
proc MPI_Recv_init*(buf: pointer; count: cint; datatype: MPI_Datatype; 
                    source: cint; tag: cint; comm: MPI_Comm; 
                    request: ptr MPI_Request): cint {.importc: "MPI_Recv_init", 
    header: "mpi.h".}
proc MPI_Recv*(buf: pointer; count: cint; datatype: MPI_Datatype; source: cint; 
               tag: cint; comm: MPI_Comm; status: ptr MPI_Status): cint {.
    importc: "MPI_Recv", header: "mpi.h".}
proc MPI_Reduce*(sendbuf: pointer; recvbuf: pointer; count: cint; 
                 datatype: MPI_Datatype; op: MPI_Op; root: cint; comm: MPI_Comm): cint {.
    importc: "MPI_Reduce", header: "mpi.h".}
proc MPI_Ireduce*(sendbuf: pointer; recvbuf: pointer; count: cint; 
                  datatype: MPI_Datatype; op: MPI_Op; root: cint; 
                  comm: MPI_Comm; request: ptr MPI_Request): cint {.
    importc: "MPI_Ireduce", header: "mpi.h".}
proc MPI_Reduce_local*(inbuf: pointer; inoutbuf: pointer; count: cint; 
                       datatype: MPI_Datatype; op: MPI_Op): cint {.
    importc: "MPI_Reduce_local", header: "mpi.h".}
proc MPI_Reduce_scatter*(sendbuf: pointer; recvbuf: pointer; 
                         recvcounts: ptr cint; datatype: MPI_Datatype; 
                         op: MPI_Op; comm: MPI_Comm): cint {.
    importc: "MPI_Reduce_scatter", header: "mpi.h".}
proc MPI_Ireduce_scatter*(sendbuf: pointer; recvbuf: pointer; 
                          recvcounts: ptr cint; datatype: MPI_Datatype; 
                          op: MPI_Op; comm: MPI_Comm; request: ptr MPI_Request): cint {.
    importc: "MPI_Ireduce_scatter", header: "mpi.h".}
proc MPI_Reduce_scatter_block*(sendbuf: pointer; recvbuf: pointer; 
                               recvcount: cint; datatype: MPI_Datatype; 
                               op: MPI_Op; comm: MPI_Comm): cint {.
    importc: "MPI_Reduce_scatter_block", header: "mpi.h".}
proc MPI_Ireduce_scatter_block*(sendbuf: pointer; recvbuf: pointer; 
                                recvcount: cint; datatype: MPI_Datatype; 
                                op: MPI_Op; comm: MPI_Comm; 
                                request: ptr MPI_Request): cint {.
    importc: "MPI_Ireduce_scatter_block", header: "mpi.h".}
proc MPI_Register_datarep*(datarep: cstring; read_conversion_fn: ptr MPI_Datarep_conversion_function; 
    write_conversion_fn: ptr MPI_Datarep_conversion_function; 
    dtype_file_extent_fn: ptr MPI_Datarep_extent_function; extra_state: pointer): cint {.
    importc: "MPI_Register_datarep", header: "mpi.h".}
proc MPI_Request_free*(request: ptr MPI_Request): cint {.
    importc: "MPI_Request_free", header: "mpi.h".}
proc MPI_Request_get_status*(request: MPI_Request; flag: ptr cint; 
                             status: ptr MPI_Status): cint {.
    importc: "MPI_Request_get_status", header: "mpi.h".}
proc MPI_Rget*(origin_addr: pointer; origin_count: cint; 
               origin_datatype: MPI_Datatype; target_rank: cint; 
               target_disp: MPI_Aint; target_count: cint; 
               target_datatype: MPI_Datatype; win: MPI_Win; 
               request: ptr MPI_Request): cint {.importc: "MPI_Rget", 
    header: "mpi.h".}
proc MPI_Rget_accumulate*(origin_addr: pointer; origin_count: cint; 
                          origin_datatype: MPI_Datatype; result_addr: pointer; 
                          result_count: cint; result_datatype: MPI_Datatype; 
                          target_rank: cint; target_disp: MPI_Aint; 
                          target_count: cint; target_datatype: MPI_Datatype; 
                          op: MPI_Op; win: MPI_Win; request: ptr MPI_Request): cint {.
    importc: "MPI_Rget_accumulate", header: "mpi.h".}
proc MPI_Rput*(origin_addr: pointer; origin_count: cint; 
               origin_datatype: MPI_Datatype; target_rank: cint; 
               target_disp: MPI_Aint; target_cout: cint; 
               target_datatype: MPI_Datatype; win: MPI_Win; 
               request: ptr MPI_Request): cint {.importc: "MPI_Rput", 
    header: "mpi.h".}
proc MPI_Rsend*(ibuf: pointer; count: cint; datatype: MPI_Datatype; dest: cint; 
                tag: cint; comm: MPI_Comm): cint {.importc: "MPI_Rsend", 
    header: "mpi.h".}
proc MPI_Rsend_init*(buf: pointer; count: cint; datatype: MPI_Datatype; 
                     dest: cint; tag: cint; comm: MPI_Comm; 
                     request: ptr MPI_Request): cint {.
    importc: "MPI_Rsend_init", header: "mpi.h".}
proc MPI_Scan*(sendbuf: pointer; recvbuf: pointer; count: cint; 
               datatype: MPI_Datatype; op: MPI_Op; comm: MPI_Comm): cint {.
    importc: "MPI_Scan", header: "mpi.h".}
proc MPI_Iscan*(sendbuf: pointer; recvbuf: pointer; count: cint; 
                datatype: MPI_Datatype; op: MPI_Op; comm: MPI_Comm; 
                request: ptr MPI_Request): cint {.importc: "MPI_Iscan", 
    header: "mpi.h".}
proc MPI_Scatter*(sendbuf: pointer; sendcount: cint; sendtype: MPI_Datatype; 
                  recvbuf: pointer; recvcount: cint; recvtype: MPI_Datatype; 
                  root: cint; comm: MPI_Comm): cint {.importc: "MPI_Scatter", 
    header: "mpi.h".}
proc MPI_Iscatter*(sendbuf: pointer; sendcount: cint; sendtype: MPI_Datatype; 
                   recvbuf: pointer; recvcount: cint; recvtype: MPI_Datatype; 
                   root: cint; comm: MPI_Comm; request: ptr MPI_Request): cint {.
    importc: "MPI_Iscatter", header: "mpi.h".}
proc MPI_Scatterv*(sendbuf: pointer; sendcounts: ptr cint; displs: ptr cint; 
                   sendtype: MPI_Datatype; recvbuf: pointer; recvcount: cint; 
                   recvtype: MPI_Datatype; root: cint; comm: MPI_Comm): cint {.
    importc: "MPI_Scatterv", header: "mpi.h".}
proc MPI_Iscatterv*(sendbuf: pointer; sendcounts: ptr cint; displs: ptr cint; 
                    sendtype: MPI_Datatype; recvbuf: pointer; recvcount: cint; 
                    recvtype: MPI_Datatype; root: cint; comm: MPI_Comm; 
                    request: ptr MPI_Request): cint {.importc: "MPI_Iscatterv", 
    header: "mpi.h".}
proc MPI_Send_init*(buf: pointer; count: cint; datatype: MPI_Datatype; 
                    dest: cint; tag: cint; comm: MPI_Comm; 
                    request: ptr MPI_Request): cint {.importc: "MPI_Send_init", 
    header: "mpi.h".}
proc MPI_Send*(buf: pointer; count: cint; datatype: MPI_Datatype; dest: cint; 
               tag: cint; comm: MPI_Comm): cint {.importc: "MPI_Send", 
    header: "mpi.h".}
proc MPI_Sendrecv*(sendbuf: pointer; sendcount: cint; sendtype: MPI_Datatype; 
                   dest: cint; sendtag: cint; recvbuf: pointer; recvcount: cint; 
                   recvtype: MPI_Datatype; source: cint; recvtag: cint; 
                   comm: MPI_Comm; status: ptr MPI_Status): cint {.
    importc: "MPI_Sendrecv", header: "mpi.h".}
proc MPI_Sendrecv_replace*(buf: pointer; count: cint; datatype: MPI_Datatype; 
                           dest: cint; sendtag: cint; source: cint; 
                           recvtag: cint; comm: MPI_Comm; status: ptr MPI_Status): cint {.
    importc: "MPI_Sendrecv_replace", header: "mpi.h".}
proc MPI_Ssend_init*(buf: pointer; count: cint; datatype: MPI_Datatype; 
                     dest: cint; tag: cint; comm: MPI_Comm; 
                     request: ptr MPI_Request): cint {.
    importc: "MPI_Ssend_init", header: "mpi.h".}
proc MPI_Ssend*(buf: pointer; count: cint; datatype: MPI_Datatype; dest: cint; 
                tag: cint; comm: MPI_Comm): cint {.importc: "MPI_Ssend", 
    header: "mpi.h".}
proc MPI_Start*(request: ptr MPI_Request): cint {.importc: "MPI_Start", 
    header: "mpi.h".}
proc MPI_Startall*(count: cint; array_of_requests: ptr MPI_Request): cint {.
    importc: "MPI_Startall", header: "mpi.h".}
proc MPI_Status_set_cancelled*(status: ptr MPI_Status; flag: cint): cint {.
    importc: "MPI_Status_set_cancelled", header: "mpi.h".}
proc MPI_Status_set_elements*(status: ptr MPI_Status; datatype: MPI_Datatype; 
                              count: cint): cint {.
    importc: "MPI_Status_set_elements", header: "mpi.h".}
proc MPI_Status_set_elements_x*(status: ptr MPI_Status; datatype: MPI_Datatype; 
                                count: MPI_Count): cint {.
    importc: "MPI_Status_set_elements_x", header: "mpi.h".}
proc MPI_Testall*(count: cint; array_of_requests: ptr MPI_Request; 
                  flag: ptr cint; array_of_statuses: ptr MPI_Status): cint {.
    importc: "MPI_Testall", header: "mpi.h".}
proc MPI_Testany*(count: cint; array_of_requests: ptr MPI_Request; 
                  index: ptr cint; flag: ptr cint; status: ptr MPI_Status): cint {.
    importc: "MPI_Testany", header: "mpi.h".}
proc MPI_Test*(request: ptr MPI_Request; flag: ptr cint; status: ptr MPI_Status): cint {.
    importc: "MPI_Test", header: "mpi.h".}
proc MPI_Test_cancelled*(status: ptr MPI_Status; flag: ptr cint): cint {.
    importc: "MPI_Test_cancelled", header: "mpi.h".}
proc MPI_Testsome*(incount: cint; array_of_requests: ptr MPI_Request; 
                   outcount: ptr cint; array_of_indices: ptr cint; 
                   array_of_statuses: ptr MPI_Status): cint {.
    importc: "MPI_Testsome", header: "mpi.h".}
proc MPI_Topo_test*(comm: MPI_Comm; status: ptr cint): cint {.
    importc: "MPI_Topo_test", header: "mpi.h".}
proc MPI_Type_commit*(`type`: ptr MPI_Datatype): cint {.
    importc: "MPI_Type_commit", header: "mpi.h".}
proc MPI_Type_contiguous*(count: cint; oldtype: MPI_Datatype; 
                          newtype: ptr MPI_Datatype): cint {.
    importc: "MPI_Type_contiguous", header: "mpi.h".}
proc MPI_Type_create_darray*(size: cint; rank: cint; ndims: cint; 
                             gsize_array: ptr cint; distrib_array: ptr cint; 
                             darg_array: ptr cint; psize_array: ptr cint; 
                             order: cint; oldtype: MPI_Datatype; 
                             newtype: ptr MPI_Datatype): cint {.
    importc: "MPI_Type_create_darray", header: "mpi.h".}
proc MPI_Type_create_f90_complex*(p: cint; r: cint; newtype: ptr MPI_Datatype): cint {.
    importc: "MPI_Type_create_f90_complex", header: "mpi.h".}
proc MPI_Type_create_f90_integer*(r: cint; newtype: ptr MPI_Datatype): cint {.
    importc: "MPI_Type_create_f90_integer", header: "mpi.h".}
proc MPI_Type_create_f90_real*(p: cint; r: cint; newtype: ptr MPI_Datatype): cint {.
    importc: "MPI_Type_create_f90_real", header: "mpi.h".}
proc MPI_Type_create_hindexed_block*(count: cint; blocklength: cint; 
                                     array_of_displacements: ptr MPI_Aint; 
                                     oldtype: MPI_Datatype; 
                                     newtype: ptr MPI_Datatype): cint {.
    importc: "MPI_Type_create_hindexed_block", header: "mpi.h".}
proc MPI_Type_create_hindexed*(count: cint; array_of_blocklengths: ptr cint; 
                               array_of_displacements: ptr MPI_Aint; 
                               oldtype: MPI_Datatype; newtype: ptr MPI_Datatype): cint {.
    importc: "MPI_Type_create_hindexed", header: "mpi.h".}
proc MPI_Type_create_hvector*(count: cint; blocklength: cint; stride: MPI_Aint; 
                              oldtype: MPI_Datatype; newtype: ptr MPI_Datatype): cint {.
    importc: "MPI_Type_create_hvector", header: "mpi.h".}
proc MPI_Type_create_keyval*(type_copy_attr_fn: ptr MPI_Type_copy_attr_function; 
    type_delete_attr_fn: ptr MPI_Type_delete_attr_function; 
                             type_keyval: ptr cint; extra_state: pointer): cint {.
    importc: "MPI_Type_create_keyval", header: "mpi.h".}
proc MPI_Type_create_indexed_block*(count: cint; blocklength: cint; 
                                    array_of_displacements: ptr cint; 
                                    oldtype: MPI_Datatype; 
                                    newtype: ptr MPI_Datatype): cint {.
    importc: "MPI_Type_create_indexed_block", header: "mpi.h".}
proc MPI_Type_create_struct*(count: cint; array_of_block_lengths: ptr cint; 
                             array_of_displacements: ptr MPI_Aint; 
                             array_of_types: ptr MPI_Datatype; 
                             newtype: ptr MPI_Datatype): cint {.
    importc: "MPI_Type_create_struct", header: "mpi.h".}
proc MPI_Type_create_subarray*(ndims: cint; size_array: ptr cint; 
                               subsize_array: ptr cint; start_array: ptr cint; 
                               order: cint; oldtype: MPI_Datatype; 
                               newtype: ptr MPI_Datatype): cint {.
    importc: "MPI_Type_create_subarray", header: "mpi.h".}
proc MPI_Type_create_resized*(oldtype: MPI_Datatype; lb: MPI_Aint; 
                              extent: MPI_Aint; newtype: ptr MPI_Datatype): cint {.
    importc: "MPI_Type_create_resized", header: "mpi.h".}
proc MPI_Type_delete_attr*(`type`: MPI_Datatype; type_keyval: cint): cint {.
    importc: "MPI_Type_delete_attr", header: "mpi.h".}
proc MPI_Type_dup*(`type`: MPI_Datatype; newtype: ptr MPI_Datatype): cint {.
    importc: "MPI_Type_dup", header: "mpi.h".}
proc MPI_Type_extent*(`type`: MPI_Datatype; extent: ptr MPI_Aint): cint {.
    importc: "MPI_Type_extent", header: "mpi.h".}
proc MPI_Type_free*(`type`: ptr MPI_Datatype): cint {.importc: "MPI_Type_free", 
    header: "mpi.h".}
proc MPI_Type_free_keyval*(type_keyval: ptr cint): cint {.
    importc: "MPI_Type_free_keyval", header: "mpi.h".}
proc MPI_Type_get_attr*(`type`: MPI_Datatype; type_keyval: cint; 
                        attribute_val: pointer; flag: ptr cint): cint {.
    importc: "MPI_Type_get_attr", header: "mpi.h".}
proc MPI_Type_get_contents*(mtype: MPI_Datatype; max_integers: cint; 
                            max_addresses: cint; max_datatypes: cint; 
                            array_of_integers: ptr cint; 
                            array_of_addresses: ptr MPI_Aint; 
                            array_of_datatypes: ptr MPI_Datatype): cint {.
    importc: "MPI_Type_get_contents", header: "mpi.h".}
proc MPI_Type_get_envelope*(`type`: MPI_Datatype; num_integers: ptr cint; 
                            num_addresses: ptr cint; num_datatypes: ptr cint; 
                            combiner: ptr cint): cint {.
    importc: "MPI_Type_get_envelope", header: "mpi.h".}
proc MPI_Type_get_extent*(`type`: MPI_Datatype; lb: ptr MPI_Aint; 
                          extent: ptr MPI_Aint): cint {.
    importc: "MPI_Type_get_extent", header: "mpi.h".}
proc MPI_Type_get_extent_x*(`type`: MPI_Datatype; lb: ptr MPI_Count; 
                            extent: ptr MPI_Count): cint {.
    importc: "MPI_Type_get_extent_x", header: "mpi.h".}
proc MPI_Type_get_name*(`type`: MPI_Datatype; type_name: cstring; 
                        resultlen: ptr cint): cint {.
    importc: "MPI_Type_get_name", header: "mpi.h".}
proc MPI_Type_get_true_extent*(datatype: MPI_Datatype; true_lb: ptr MPI_Aint; 
                               true_extent: ptr MPI_Aint): cint {.
    importc: "MPI_Type_get_true_extent", header: "mpi.h".}
proc MPI_Type_get_true_extent_x*(datatype: MPI_Datatype; true_lb: ptr MPI_Count; 
                                 true_extent: ptr MPI_Count): cint {.
    importc: "MPI_Type_get_true_extent_x", header: "mpi.h".}
proc MPI_Type_hindexed*(count: cint; array_of_blocklengths: ptr cint; 
                        array_of_displacements: ptr MPI_Aint; 
                        oldtype: MPI_Datatype; newtype: ptr MPI_Datatype): cint {.
    importc: "MPI_Type_hindexed", header: "mpi.h".}
proc MPI_Type_hvector*(count: cint; blocklength: cint; stride: MPI_Aint; 
                       oldtype: MPI_Datatype; newtype: ptr MPI_Datatype): cint {.
    importc: "MPI_Type_hvector", header: "mpi.h".}
proc MPI_Type_indexed*(count: cint; array_of_blocklengths: ptr cint; 
                       array_of_displacements: ptr cint; oldtype: MPI_Datatype; 
                       newtype: ptr MPI_Datatype): cint {.
    importc: "MPI_Type_indexed", header: "mpi.h".}
proc MPI_Type_lb*(`type`: MPI_Datatype; lb: ptr MPI_Aint): cint {.
    importc: "MPI_Type_lb", header: "mpi.h".}
proc MPI_Type_match_size*(typeclass: cint; size: cint; `type`: ptr MPI_Datatype): cint {.
    importc: "MPI_Type_match_size", header: "mpi.h".}
proc MPI_Type_set_attr*(`type`: MPI_Datatype; type_keyval: cint; 
                        attr_val: pointer): cint {.importc: "MPI_Type_set_attr", 
    header: "mpi.h".}
proc MPI_Type_set_name*(`type`: MPI_Datatype; type_name: cstring): cint {.
    importc: "MPI_Type_set_name", header: "mpi.h".}
proc MPI_Type_size*(`type`: MPI_Datatype; size: ptr cint): cint {.
    importc: "MPI_Type_size", header: "mpi.h".}
proc MPI_Type_size_x*(`type`: MPI_Datatype; size: ptr MPI_Count): cint {.
    importc: "MPI_Type_size_x", header: "mpi.h".}
proc MPI_Type_struct*(count: cint; array_of_blocklengths: ptr cint; 
                      array_of_displacements: ptr MPI_Aint; 
                      array_of_types: ptr MPI_Datatype; 
                      newtype: ptr MPI_Datatype): cint {.
    importc: "MPI_Type_struct", header: "mpi.h".}
proc MPI_Type_ub*(mtype: MPI_Datatype; ub: ptr MPI_Aint): cint {.
    importc: "MPI_Type_ub", header: "mpi.h".}
proc MPI_Type_vector*(count: cint; blocklength: cint; stride: cint; 
                      oldtype: MPI_Datatype; newtype: ptr MPI_Datatype): cint {.
    importc: "MPI_Type_vector", header: "mpi.h".}
proc MPI_Unpack*(inbuf: pointer; insize: cint; position: ptr cint; 
                 outbuf: pointer; outcount: cint; datatype: MPI_Datatype; 
                 comm: MPI_Comm): cint {.importc: "MPI_Unpack", header: "mpi.h".}
proc MPI_Unpublish_name*(service_name: cstring; info: MPI_Info; 
                         port_name: cstring): cint {.
    importc: "MPI_Unpublish_name", header: "mpi.h".}
proc MPI_Unpack_external*(datarep: ptr char; inbuf: pointer; insize: MPI_Aint; 
                          position: ptr MPI_Aint; outbuf: pointer; 
                          outcount: cint; datatype: MPI_Datatype): cint {.
    importc: "MPI_Unpack_external", header: "mpi.h".}
proc MPI_Waitall*(count: cint; array_of_requests: ptr MPI_Request; 
                  array_of_statuses: ptr MPI_Status): cint {.
    importc: "MPI_Waitall", header: "mpi.h".}
proc MPI_Waitany*(count: cint; array_of_requests: ptr MPI_Request; 
                  index: ptr cint; status: ptr MPI_Status): cint {.
    importc: "MPI_Waitany", header: "mpi.h".}
proc MPI_Wait*(request: ptr MPI_Request; status: ptr MPI_Status): cint {.
    importc: "MPI_Wait", header: "mpi.h".}
proc MPI_Waitsome*(incount: cint; array_of_requests: ptr MPI_Request; 
                   outcount: ptr cint; array_of_indices: ptr cint; 
                   array_of_statuses: ptr MPI_Status): cint {.
    importc: "MPI_Waitsome", header: "mpi.h".}
proc MPI_Win_allocate*(size: MPI_Aint; disp_unit: cint; info: MPI_Info; 
                       comm: MPI_Comm; baseptr: pointer; win: ptr MPI_Win): cint {.
    importc: "MPI_Win_allocate", header: "mpi.h".}
proc MPI_Win_allocate_shared*(size: MPI_Aint; disp_unit: cint; info: MPI_Info; 
                              comm: MPI_Comm; baseptr: pointer; win: ptr MPI_Win): cint {.
    importc: "MPI_Win_allocate_shared", header: "mpi.h".}
proc MPI_Win_attach*(win: MPI_Win; base: pointer; size: MPI_Aint): cint {.
    importc: "MPI_Win_attach", header: "mpi.h".}
proc MPI_Win_call_errhandler*(win: MPI_Win; errorcode: cint): cint {.
    importc: "MPI_Win_call_errhandler", header: "mpi.h".}
proc MPI_Win_complete*(win: MPI_Win): cint {.importc: "MPI_Win_complete", 
    header: "mpi.h".}
proc MPI_Win_create*(base: pointer; size: MPI_Aint; disp_unit: cint; 
                     info: MPI_Info; comm: MPI_Comm; win: ptr MPI_Win): cint {.
    importc: "MPI_Win_create", header: "mpi.h".}
proc MPI_Win_create_dynamic*(info: MPI_Info; comm: MPI_Comm; win: ptr MPI_Win): cint {.
    importc: "MPI_Win_create_dynamic", header: "mpi.h".}
proc MPI_Win_create_errhandler*(function: ptr MPI_Win_errhandler_function; 
                                errhandler: ptr MPI_Errhandler): cint {.
    importc: "MPI_Win_create_errhandler", header: "mpi.h".}
proc MPI_Win_create_keyval*(win_copy_attr_fn: ptr MPI_Win_copy_attr_function; 
    win_delete_attr_fn: ptr MPI_Win_delete_attr_function; win_keyval: ptr cint; 
                            extra_state: pointer): cint {.
    importc: "MPI_Win_create_keyval", header: "mpi.h".}
proc MPI_Win_delete_attr*(win: MPI_Win; win_keyval: cint): cint {.
    importc: "MPI_Win_delete_attr", header: "mpi.h".}
proc MPI_Win_detach*(win: MPI_Win; base: pointer): cint {.
    importc: "MPI_Win_detach", header: "mpi.h".}
proc MPI_Win_fence*(assert: cint; win: MPI_Win): cint {.
    importc: "MPI_Win_fence", header: "mpi.h".}
proc MPI_Win_flush*(rank: cint; win: MPI_Win): cint {.importc: "MPI_Win_flush", 
    header: "mpi.h".}
proc MPI_Win_flush_all*(win: MPI_Win): cint {.importc: "MPI_Win_flush_all", 
    header: "mpi.h".}
proc MPI_Win_flush_local*(rank: cint; win: MPI_Win): cint {.
    importc: "MPI_Win_flush_local", header: "mpi.h".}
proc MPI_Win_flush_local_all*(win: MPI_Win): cint {.
    importc: "MPI_Win_flush_local_all", header: "mpi.h".}
proc MPI_Win_free*(win: ptr MPI_Win): cint {.importc: "MPI_Win_free", 
    header: "mpi.h".}
proc MPI_Win_free_keyval*(win_keyval: ptr cint): cint {.
    importc: "MPI_Win_free_keyval", header: "mpi.h".}
proc MPI_Win_get_attr*(win: MPI_Win; win_keyval: cint; attribute_val: pointer; 
                       flag: ptr cint): cint {.importc: "MPI_Win_get_attr", 
    header: "mpi.h".}
proc MPI_Win_get_errhandler*(win: MPI_Win; errhandler: ptr MPI_Errhandler): cint {.
    importc: "MPI_Win_get_errhandler", header: "mpi.h".}
proc MPI_Win_get_group*(win: MPI_Win; group: ptr MPI_Group): cint {.
    importc: "MPI_Win_get_group", header: "mpi.h".}
proc MPI_Win_get_info*(win: MPI_Win; info_used: ptr MPI_Info): cint {.
    importc: "MPI_Win_get_info", header: "mpi.h".}
proc MPI_Win_get_name*(win: MPI_Win; win_name: cstring; resultlen: ptr cint): cint {.
    importc: "MPI_Win_get_name", header: "mpi.h".}
proc MPI_Win_lock*(lock_type: cint; rank: cint; assert: cint; win: MPI_Win): cint {.
    importc: "MPI_Win_lock", header: "mpi.h".}
proc MPI_Win_lock_all*(assert: cint; win: MPI_Win): cint {.
    importc: "MPI_Win_lock_all", header: "mpi.h".}
proc MPI_Win_post*(group: MPI_Group; assert: cint; win: MPI_Win): cint {.
    importc: "MPI_Win_post", header: "mpi.h".}
proc MPI_Win_set_attr*(win: MPI_Win; win_keyval: cint; attribute_val: pointer): cint {.
    importc: "MPI_Win_set_attr", header: "mpi.h".}
proc MPI_Win_set_errhandler*(win: MPI_Win; errhandler: MPI_Errhandler): cint {.
    importc: "MPI_Win_set_errhandler", header: "mpi.h".}
proc MPI_Win_set_info*(win: MPI_Win; info: MPI_Info): cint {.
    importc: "MPI_Win_set_info", header: "mpi.h".}
proc MPI_Win_set_name*(win: MPI_Win; win_name: cstring): cint {.
    importc: "MPI_Win_set_name", header: "mpi.h".}
proc MPI_Win_shared_query*(win: MPI_Win; rank: cint; size: ptr MPI_Aint; 
                           disp_unit: ptr cint; baseptr: pointer): cint {.
    importc: "MPI_Win_shared_query", header: "mpi.h".}
proc MPI_Win_start*(group: MPI_Group; assert: cint; win: MPI_Win): cint {.
    importc: "MPI_Win_start", header: "mpi.h".}
proc MPI_Win_sync*(win: MPI_Win): cint {.importc: "MPI_Win_sync", 
    header: "mpi.h".}
proc MPI_Win_test*(win: MPI_Win; flag: ptr cint): cint {.
    importc: "MPI_Win_test", header: "mpi.h".}
proc MPI_Win_unlock*(rank: cint; win: MPI_Win): cint {.
    importc: "MPI_Win_unlock", header: "mpi.h".}
proc MPI_Win_unlock_all*(win: MPI_Win): cint {.importc: "MPI_Win_unlock_all", 
    header: "mpi.h".}
proc MPI_Win_wait*(win: MPI_Win): cint {.importc: "MPI_Win_wait", 
    header: "mpi.h".}
proc MPI_Wtick*(): cdouble {.importc: "MPI_Wtick", header: "mpi.h".}
proc MPI_Wtime*(): cdouble {.importc: "MPI_Wtime", header: "mpi.h".}
