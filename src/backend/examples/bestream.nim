## Device memory bandwidth of gpuFor kernels, as STREAM: copy c = a, scale
## b = s c, add c = a + b, triad a = b + s c, the best of reps timings, with
## 1, 4 and 8 independent elements i = t + j n/k, j < k, per thread t.
##   -n: reals per array (default 2^27, 1 GiB; arrays of 64 MiB stay in the
##       L2 of a PVC tile)  -reps: repetitions
import qex
import backend/accel
import base/metaUtils
import times

proc stream(n, reps: int) =
  let a = cast[ptr UncheckedArray[float]](gpuMalloc(n*sizeof(float)))
  let b = cast[ptr UncheckedArray[float]](gpuMalloc(n*sizeof(float)))
  let c = cast[ptr UncheckedArray[float]](gpuMalloc(n*sizeof(float)))
  gpuFor(i, n):
    a[i] = 1.0
    b[i] = 2.0
    c[i] = 0.0
  let s = 3.0
  template bench(name: string; arrays: int; body: untyped) =
    var best = 0.0
    for r in 0..<reps:
      let t0 = epochTime()
      body
      let dt = epochTime() - t0
      if r == 0 or dt < best: best = dt
    echo name, ": ", 1e-9*float(arrays*n*sizeof(float))/best, " GB/s"
  template ops(k: static int) =
    let m = n div k
    bench("copy x" & $k, 2):
      gpuFor(t, m):
        forStatic j, 0, k-1: c[t + j*m] = a[t + j*m]
    bench("scale x" & $k, 2):
      gpuFor(t, m):
        forStatic j, 0, k-1: b[t + j*m] = s*c[t + j*m]
    bench("add x" & $k, 3):
      gpuFor(t, m):
        forStatic j, 0, k-1: c[t + j*m] = a[t + j*m] + b[t + j*m]
    bench("triad x" & $k, 3):
      gpuFor(t, m):
        forStatic j, 0, k-1: a[t + j*m] = b[t + j*m] + s*c[t + j*m]
  ops(1)
  ops(4)
  ops(8)
  for p in [a, b, c]: gpuFree(p)

qexInit()
stream(intParam("n", 1 shl 27), intParam("reps", 20))
qexFinalize()
