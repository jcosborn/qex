## Communication totals for -d:commStats, reported at communicator shutdown.
## Thread 0 counts exchanges, bytes sent, and time waiting for communication.

when defined(commStats):
  import std/[monotimes, times, strformat]
  import commsTypes

  # Templates instantiated in generic callers need public field names.
  type CommStats = object
    shifts*, msgs*, bytes*: int  # shift starts; those with remote messages; bytes sent
    wait*: float                 # waits for shift messages
    gathers*, gmsgs*, gbytes*: int
    gwait*: float                # gather barriers and message waits
    reduces*: int
    rwait*: float                # all-reduce calls

  var commStats: CommStats

  template commAdd*(f: untyped, n: int) = commStats.f += n

  template commTime*(f, body: untyped) =
    let t0 = getMonoTime()
    body
    commStats.f += float((getMonoTime() - t0).inNanoseconds) * 1e-9

  proc commReport*(label: string) =
    let v = commStats
    if myRank == 0:
      echo &"[commStats] {label}: shifts msgs MB wait/s gathers gmsgs gMB gwait/s reduces rwait/s"
      echo &"[commStats] {label}: {v.shifts} {v.msgs} {v.bytes.float*1e-6:.3f} {v.wait:.4f} " &
        &"{v.gathers} {v.gmsgs} {v.gbytes.float*1e-6:.3f} {v.gwait:.4f} {v.reduces} {v.rwait:.4f}"
    commStats = CommStats()

else:
  template commAdd*(f: untyped, n: int) = discard
  template commTime*(f, body: untyped) = body
  template commReport*(label: string) = discard

template commReduce*(body: untyped) =
  commAdd(reduces, 1)
  commTime(rwait, body)
