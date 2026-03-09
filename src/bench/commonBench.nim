import base/[threading,profile], comms/comms

type
  BenchResult* = object
    reps*: int
    secs*: float
    warms*: int

type
  Bench* = object
    br*: BenchResult
    comm*: Comm
    nrep*: int
    t0: TicType
    dt: float

proc newBench*: Bench =
  result.comm = getDefaultComm()
  result.nrep = 1
  result.dt = 0.0
  result.br.warms = 0
  result.br.reps = 0

proc start*(b: var Bench) =
  b.br.warms += b.br.reps
  if b.br.warms > 0:
    let nnrep = 1 + int(1.1*b.nrep.float/(b.dt+1e-9))
    b.nrep = min(10*b.nrep, nnrep)
  b.comm.barrier()
  b.t0 = getTics()

proc stop*(b: var Bench) =
  let t1 = getTics()
  b.comm.barrier()
  b.dt = ticDiffSecs(t1,b.t0)
  b.comm.max(b.dt)
  #echo "nrep: ", nrep, "  dt: ", dt
  b.br.reps = b.nrep
  b.br.secs = b.dt

template notDone*(b: Bench): bool =
  b.br.warms<=0 or b.dt<1.0

template perNs*(b: Bench): float =
  # 1/(iteration time in ns)
  1e-9 * b.br.reps.float / b.br.secs

template benchSingle*(bb: Bench, body: untyped) =
  proc benchSingleProc(b: var Bench) {.gensym.} =
    while b.notDone:
      b.start
      body
      b.stop
  benchSingleProc(bb)

#[
template benchSingle*(body: untyped): BenchResult =
  var b = newBench()
  benchSingle2(b):
    for rep in 1..b.nrep:
      body
  b.br
]#

#[
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
]#

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
