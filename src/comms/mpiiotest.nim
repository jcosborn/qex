import mpi
import times
import os
import strutils

type TimerObj = tuple[s:string,d:float,n:int]
var timer: float
var deltas = newSeq[ptr TimerObj](0)
template tic =
  timer = epochTime()
template toc(st: string) =
  block:
    var first{.global.} = true
    var to{.global.}: TimerObj
    if first:
      first = false
      deltas.add to.addr
      to.s = st
      to.d = 0.0
      to.n = 0
    let t = epochTime()
    to.d += t - timer
    timer = t
    inc to.n

discard MPI_init()
let comm = MPI_COMM_WORLD
var rank, size: cint
discard MPI_Comm_rank(comm, rank.addr)
discard MPI_Comm_size(comm, size.addr)
var fn = "testout.dat"
if paramCount()>=1:
  fn = paramStr(1)
var sinfo = true

proc showInfo(fh: MPI_File) =
  var info: MPI_Info
  discard MPI_File_get_info(fh, info.addr)
  var nkeys: cint
  discard MPI_Info_get_nkeys(info, nkeys.addr)
  var key = newSeq[char](MPI_MAX_INFO_KEY)
  var ckey = cast[cstring](key[0].addr)
  var val = newSeq[char](MPI_MAX_INFO_VAL)
  var cval = cast[cstring](val[0].addr)
  var lval = val.len.cint
  var flag: cint
  for i in 0..<nkeys:
    discard MPI_Info_get_nthkey(info, i, ckey)
    discard MPI_Info_get(info, ckey, lval, cval, flag.addr)
    echo ckey, ": ", cval

proc testContig(n: int) =
  var buf = newSeq[cint](n)
  var filetype: MPI_Datatype
  var disp = MPI_Offset(rank*sizeof(cint)*n)
  var etype = MPI_INT
  discard MPI_Type_contiguous(n.cint, etype, filetype.addr)
  discard MPI_Type_commit(filetype.addr)

  var info: MPI_Info
  info = MPI_INFO_NULL
  #discard MPI_Info_create(info.addr)
  #let cbnodes = size div 8
  #discard MPI_Info_set(info, "cb_nodes", $cbnodes);

  # MPI_File_open(comm, filename, MPI_MODE_RDONLY, info, &fh);
  #MPI_TYPE_CREATE_SUBARRAY(ndims, iarray_of_sizes, iarray_of_subsizes,
  #                         iarray_of_starts, MPI_ORDER_FORTRAN, MPI_REAL,
  #                         ifiletype, ierr)
  #MPI_TYPE_COMMIT(ifiletype#, ierr)

  var fh: MPI_File
  let wmode = MPI_MODE_CREATE or MPI_MODE_WRONLY
  tic()
  discard MPI_File_open(comm, fn, wmode, info, fh.addr)
  toc("MPI_File_open write")
  discard MPI_File_set_view(fh, disp, etype, filetype, "native", info)
  if rank==0 and sinfo: showInfo(fh)
  toc("MPI_File_set_view write")
  discard MPI_File_write_all(fh,buf[0].voidaddr,n.cint,etype,MPI_STATUS_IGNORE)
  discard MPI_barrier(comm)
  toc("MPI_File_write_all")
  discard MPI_File_close(fh.addr)
  discard MPI_barrier(comm)
  toc("MPI_File_close write")

  let rmode = MPI_MODE_RDONLY
  discard MPI_File_open(comm, fn, rmode, info, fh.addr)
  toc("MPI_File_open read")
  discard MPI_File_set_view(fh, disp, etype, filetype, "native", info)
  toc("MPI_File_set_view read")
  discard MPI_File_read_all(fh,buf[0].voidaddr,n.cint,etype,MPI_STATUS_IGNORE)
  discard MPI_barrier(comm)
  toc("MPI_File_read_all")
  discard MPI_File_close(fh.addr)
  discard MPI_barrier(comm)
  toc("MPI_File_close read")

  #discard MPI_Info_free(info.addr)
  if rank==0: removeFile(fn)


if rank==0:
  echo "total ranks: ", size
  echo "using file: ", fn

var n = 1024
var nmax = 2*72*16*1024
if paramCount()>=2:
  nmax = parseInt(paramStr(2))*1024*1024
while n<=nmax:
  if rank==0:
    let b = formatSize(n*size*sizeof(cint))
    echo "total file size: ", b
  testContig(n)
  sinfo = false
  testContig(n)
  for i in 0..<deltas.len:
    deltas[i].d = 0.0
    deltas[i].n = 0
  testContig(n)
  testContig(n)
  for i in 0..<deltas.len:
    if rank==0:
      let sp = max(0, 23-deltas[i].s.len)
      let t = formatFloat(1e6*deltas[i].d/deltas[i].n.float, ffDecimal, 3)
      echo deltas[i].s, spaces(sp), ": ", align(t,15)
  n *= 2

discard MPI_Finalize()
