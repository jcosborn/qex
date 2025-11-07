import base/[threading,profile], comms/comms

type
  BenchResult* = object
    reps*: int
    secs*: float
    warms*: int

template benchSingle*(body: auto): BenchResult =
  proc benchSingleProc: BenchResult {.gensym.} =
    let comm = getDefaultComm()
    var nrep = 1
    var dt = 0.0
    while true:
      comm.barrier()
      let t0 = getTics()
      for rep in 1..nrep:
        body
      let t1 = getTics()
      comm.barrier()
      dt = ticDiffSecs(t1,t0)
      comm.max(dt)
      #echo "nrep: ", nrep, "  dt: ", dt
      if result.warms>0 and dt>1.0: break
      result.warms += nrep
      let nnrep = 1 + int(1.1*nrep.float/(dt+1e-9))
      nrep = min(10*nrep, nnrep)
    result.reps = nrep
    result.secs = dt
  benchSingleProc()

template benchThreaded*(body: auto): BenchResult =
  proc benchThreadedProc: BenchResult {.gensym.} =
    doAssert(getDefaultComm().size==1)
    var nrep = 1
    var dt = 0.0
    while true:
      threads:
        threadBarrier()
        let t0 = getTics()
        for rep in 1..nrep:
          body
        let t1 = getTics()
        threadBarrier()
        var dtt = ticDiffSecs(t1,t0)
        threadMax(dtt)
        threadSingle: dt = dtt
      if result.warms>0 and dt>1.0: break
      result.warms += nrep
      let nnrep = 1 + int(1.1*nrep.float/(dt+1e-9))
      nrep = min(10*nrep, nnrep)
    result.reps = nrep
    result.secs = dt
  benchThreadedProc()
