import backend/timing, backend/cpugpuarray, base/metaUtils, math

proc test(vecLen, memLen: static[int]; N: int) =
  var
    x = newColorMatrixArray(vecLen,memLen,N) # array of N 3x3 single prec complex matrices
    y = newColorMatrixArray(vecLen,memLen,N)
    z = newColorMatrixArray(vecLen,memLen,N)
    rep = 0                     # accumulates the number of runs

  let
    mr = float(3 * 8 * x.T.N * x.T.N * N) / float(1024 * 1024 * 1024) # Resident memory in 2^30 bytes
    mt = 4 * mr / 3             # Memory transaction
    fp = float(8 * x.T.N * x.T.N * x.T.N * N) * 1e-9 # Floating point op / 10^9
  template timeit(label:string, s:untyped) =
    var
      R {.global.}:int
      T {.global.}:float
    threadSingle:
      R = 128                   # Base repeat
      T = 1.0                   # Time limit
    var t = timex(rep, R, s)    # Always warm up cache
    while true:
      threadSingle:
        R = min(2*R,max(R,int(R.float*0.4/t))) # set up to run for at least 0.4 sec or 2*R
      t = timex(rep, R, s)
      threadSingle: T -= t
      if T < 0: break
    threadSingle:               # Use the last R & t for performance measure
      printf("%8d %3d %d %-8s rep: %7d KB: %8.0f ms: %8.4f GF/s: %7.2f GB/s: %7.2f\n",
             N, vecLen, memLen, label, R, 1024*1024*mr, 1e3*t/R.float, fp*R.float/t, mt*R.float/t)

  threads:                      # CPU threads
    x := 0                      # set them to diagonal matrices on CPU
    y := 1
    z := 2
    timeit "CPU": x += y * z

  timeit "GPU5":                # includes kernel launching and synchronization
    onGpu(N, 32):               # Number of threads, threads per block
      x += y * z
  timeit "GPU6": onGpu(N, 64): x += y * z
  timeit "GPU7": onGpu(N, 128): x += y * z

  threads: timeit "CPU": x += y * z # back to CPU threads again

  let scale = 0.5 / (sqrt(3.0) * rep.float)
  threads:
    x *= scale
    var n = x.norm2
    threadSingle: echo "# Final scaled x.norm2: ",n,"  rep: ",rep
  x.free
  y.free
  z.free

for n in 10..26:
  forstaticUntyped v, 2, 7:
    when (1 shl v) >= (structsize(vectorizedElementType(float32)) div sizeof(float32)):
      forstaticUntyped ml, 1, 2:
        test(1 shl v, ml, 1 shl n)
