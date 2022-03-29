import globals
export globals
import threading
export threading
import comms/comms
export comms
import profile
export profile

import algorithm, strutils

var
  qexGlobalInitializers* = newseq[proc()]() ## Will be run in qexInit in forward order
  qexGlobalFinalizers* = newseq[proc()]() ## Will be run in qexFinalize in backward order
  qexGlobalPreInit* = newseq[proc()]() ## like qexGlobalInitializers but before QEX init
  qexGlobalPostFinal* = newseq[proc()]() ## like qexGlobalFinalizers but after QEX fini
  qexFinalizeComms = true
  qexStartTime: TicType

proc qexTime*: float = ticDiffSecs(getTics(), qexStartTime)

template qexLog*(s:varargs[string,`$`]) =
  let t = qexTime()
  if s.len > 0:
    echo "[", formatFloat(t,ffDecimal,3), " s] ", s.join

template qexWarn*(s:varargs[string,`$`]) =
  let ii = instantiationInfo()
  echo "Warning: ", ii.filename, ":", ii.line, ":"
  if s.len > 0:
    echo "  ", s.join

template qexError*(s:varargs[string,`$`]) =
  let ii = instantiationInfo()
  echo "Error: ", ii.filename, ":", ii.line, ":"
  if s.len > 0:
    echo "  ", s.join
  flushFile stdout
  flushFile stderr
  getComm().barrier
  qexAbort()

proc qexInit* =
  qexStartTime = getTics()
  for p in qexGlobalPreInit: p()
  threadsInit()
  commsInit()
  for p in qexGlobalInitializers: p()
  when defined(FUELCompat):
    echo "FUEL compatibility mode: ON"
  #echo "rank " & $rank & "/" & $size

proc qexSetFinalizeComms*(val: bool) =
  qexFinalizeComms = val

proc qexFinalize*() =
  flushFile stdout
  flushFile stderr
  GC_fullCollect()
  getComm().barrier
  for p in qexGlobalFinalizers.reversed: p()
  #echo("mem: (used+free)/total: (", getOccupiedMem(), "+", getFreeMem(), ")/",
  #     getTotalMem())
  #echo GC_getStatistics()
  #echo("mem: (used+free)/total: (", getOccupiedMem(), "+", getFreeMem(), ")/",
  #     getTotalMem())
  #echo GC_getStatistics()
  if qexFinalizeComms:
    commsFinalize()
  #when profileEqns:
  #echoTimers()
  for p in qexGlobalPostFinal.reversed: p()
  qexLog "Total time (Init - Finalize): ",qexTime()," seconds."

proc qexExit*(status = 0) =
  for p in qexGlobalFinalizers.reversed: p()
  commsFinalize()
  qexLog "Total time (Init - Finalize): ",qexTime()," seconds."
  quit(status)

proc qexAbort*(status = -1) =
  commsAbort(status)
