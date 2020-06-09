import ./qexInternal
import ./alignedMem

var VerboseGCStats* = false

template qexGC*(label:string = "") =
  if unlikely VerboseGCStats:
    if label.len > 0:
      echo "# "&label
    else:
      const
        ii = instantiationInfo()
        s = "# " & ii.filename & ":" & $ii.line & ":" & $ii.column
      echo s
    echo "[RAW] allocated memory: " & $getRawMemAllocated()
    echo "[RAW] used memory: " & $getRawMemUsed()
    echo "[RAW] max used memory: " & $getRawMemMaxUsed()
    echo GC_getStatistics()
  GC_fullCollect()
  if unlikely VerboseGCStats:
    echo "[RAW] allocated memory: " & $getRawMemAllocated()
    echo "[RAW] used memory: " & $getRawMemUsed()
    echo "[RAW] max used memory: " & $getRawMemMaxUsed()
    echo GC_getStatistics()
