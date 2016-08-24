import globals
import threading
import comms
import profile


proc qexInit*() =
  threadsInit()
  commsInit()
  #echo "rank " & $rank & "/" & $size

proc qexFinalize*() =
  echo("mem: (used+free)/total: (", getOccupiedMem(), "+", getFreeMem(), ")/",
       getTotalMem())
  echo GC_getStatistics()
  GC_fullCollect()
  echo("mem: (used+free)/total: (", getOccupiedMem(), "+", getFreeMem(), ")/",
       getTotalMem())
  echo GC_getStatistics()
  commsFinalize()
  #when profileEqns:
  echoTimers()

proc qexExit*(status = 0) =
  commsFinalize()
  quit(status)

proc qexAbort*(status = -1) =
  commsAbort(status)
