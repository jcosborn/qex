import globals
export globals
import threading
export threading
import comms/comms
export comms
import profile
export profile

import algorithm, times

var
  qexGlobalInitializers* = newseq[proc()]()    ## Will be run in qexInit in forward order
  qexGlobalFinalizers* = newseq[proc()]()    ## Will be run in qexFinalize in backward order
  qexStartTime: float

proc qexTime*: float = epochTime() - qexStartTime

proc qexInit* =
  qexStartTime = epochTime()
  threadsInit()
  commsInit()
  for p in qexGlobalInitializers: p()
  when defined(FUELCompat):
    echo "FUEL compatibility mode: ON"
  #echo "rank " & $rank & "/" & $size

proc qexFinalize* =
  for p in qexGlobalFinalizers.reversed: p()
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
  echo "Total time (Init - Finalize): ",qexTime()," seconds."

proc qexExit*(status = 0) =
  for p in qexGlobalFinalizers.reversed: p()
  commsFinalize()
  echo "Total time (Init - Finalize): ",qexTime()," seconds."
  quit(status)

proc qexAbort*(status = -1) =
  commsAbort(status)
