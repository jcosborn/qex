import mpi

var init = false
var shmcomm: MPI_Comm

proc getSharedMpiComm*: MpiComm =
  if not init:
    let err = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                                  MPI_INFO_NULL, addr shmcomm)
    init = true
  result = shmcomm

proc getRanksPerNode*: int =
  let c = getSharedMpiComm()
  var shmsize: int32
  let err = MPI_Comm_size(c, addr shmsize)
  result = shmsize

proc isSharedRanksFastest*: bool =
  let c = getSharedMpiComm()
  var shmsize,shmrank,rmax,rmin,allfastest: int32
  var err: int32
  err = MPI_Comm_size(c, addr shmsize)
  err = MPI_Comm_rank(c, addr shmrank)
  err = MpiAllreduce(addr shmrank, addr rmax, 1, MpiInt, MpiMax, c)
  err = MpiAllreduce(addr shmrank, addr rmin, 1, MpiInt, MpiMin, c)
  var fastest = ord( (rmax-rmin) == (shmsize-1) )
  err = MpiAllreduce(addr fastest, addr allfastest, 1, MpiInt, MpiMin, c)
  result = bool allfastest


when isMainModule:
  import comms
  commsInit()
  echo "rpn: ", getRanksPerNode()
  echo "is fastest: ", isSharedRanksFastest()
  commsFinalize()

