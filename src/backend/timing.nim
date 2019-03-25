include system/timers

template timex*(rep:var int, nn:int, s:untyped): float =
  let n = nn
  var dt {.global.}:float
  let t = getTicks()
  for i in 0..<n: s             # repeats the expression, `s`, `n` times
  var dtt = 1e-9*float(getTicks()-t) # seconds elapsed
  threadSum(dtt)
  threadSingle:
    dt = dtt/getNumThreads().float
    rep += n
  dt
