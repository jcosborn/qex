include mpitypes
include mpifuncs

template voidAddr*(x: typed): pointer = cast[pointer](x.addr)

proc MPI_Init*(): cint =
  var argc {.importc: "cmdCount", global.}: cint
  var argv {.importc: "cmdLine", global.}: cstringArray
  result = MPI_Init(argc.addr, argv.addr)

proc MPI_Init_thread*(required: cint; provided: ptr cint): cint =
  var argc {.importc: "cmdCount", global.}: cint
  var argv {.importc: "cmdLine", global.}: cstringArray
  result = MPI_Init_thread(argc.addr, argv.addr, required, provided)

when isMainModule:
  var prv = MPI_THREAD_MULTIPLE
  let err = MPI_init_thread(prv, prv.addr)
  let comm = MPI_COMM_WORLD
  var rank, size: cint
  discard MPI_Comm_rank(MPI_COMM_WORLD, rank.addr)
  discard MPI_Comm_size(MPI_COMM_WORLD, size.addr)
  var ver, subver: cint
  discard MPI_Get_version(ver.addr, subver.addr)

  if rank==0:
    echo "rank: ", rank, "  of: ", size
    echo "err: ", err
    echo "provided: ", prv
    echo "ver: ", ver, "  subver: ", subver
    stdout.flushFile
  discard MPI_barrier(comm)

  var s = "message from " & $rank
  var slen = s.len.cint
  var srank = ((rank+1)mod size).cint
  var sreq: MPI_Request
  var sstat: MPI_Status
  var tag = 10.cint
  var rlen: cint
  var rrank = ((rank+size-1)mod size).cint
  var rreq: MPI_Request
  var rstat: MPI_Status
  template vaddr(x: typed): untyped = voidAddr(x)
  discard MPI_Irecv(rlen.vaddr, 1.cint, MPI_INT, rrank, tag, comm, rreq.addr)
  discard MPI_Isend(slen.vaddr, 1.cint, MPI_INT, srank, tag, comm, sreq.addr)
  discard MPI_Wait(sreq.addr, sstat.addr)
  discard MPI_Wait(rreq.addr, rstat.addr)
  echo "rank: ", rank, "  rlen: ", rlen
  echo "rank: ", rank, "  slen: ", slen
  discard MPI_barrier(comm)

  var buf = newString(rlen)
  discard MPI_Irecv(buf[0].vaddr, rlen, MPI_CHAR, rrank, tag, comm, rreq.addr)
  discard MPI_Isend(s[0].vaddr, slen, MPI_CHAR, srank, tag, comm, sreq.addr)
  discard MPI_Wait(sreq.addr, sstat.addr)
  discard MPI_Wait(rreq.addr, rstat.addr)
  echo "rank: ", rank, "  buf: ", buf
  discard MPI_barrier(comm)

  discard MPI_finalize()
