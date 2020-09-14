import mpi

var init = false
var shmcomm: MPI_Comm

proc getRanksPerNode*: int =
  if not init:
    let err = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                                  MPI_INFO_NULL, addr shmcomm)
    init = true
  var shmsize,shmrank: int32
  let err = MPI_Comm_size(shmcomm, addr shmsize)
  #MPI_Comm_rank(shmcomm, &shmrank)
  result = shmsize
