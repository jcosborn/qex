import globals
import threading
import comms
import profile


proc qexInit*() =
  threadsInit()
  commsInit()
  #echo "rank " & $rank & "/" & $size

proc qexFinalize*() =
  commsFinalize()
  #when profileEqns:
  echoTimers()
