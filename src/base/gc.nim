import ./qexInternal

var VerboseGCStats* = false

proc qexGC*(label:string = "") =
  if unlikely VerboseGCStats:
    if label.len > 0: echo "# ",label
    echo GC_getStatistics()
  GC_fullCollect()
  if unlikely VerboseGCStats:
    echo GC_getStatistics()
