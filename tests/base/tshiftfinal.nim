import qex

# Releases shift buffers in a collection after commsFinalize, as happens when a stale
# stack word kept them through the earlier collections. On several ranks they hold QMP
# message handles, which must not go through QMP once it is finalized.
qexInit()
let lo = newLayout([8,8,8,8])
var x = lo.ColorVector()
var s = newShifters(x, 1)
qexGlobalPostFinal.add proc() =
  reset s
  GC_fullCollect()
qexFinalize()
