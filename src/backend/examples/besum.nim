## Times the ways a kernel can sum two dot products over n sites, each as
## kernel + read, in microseconds per cycle:
##   store   one value per site and no sum: the kernel and its wait alone
##   atomic  gpuAtomicAdd of doubles into 512 slots, whose result depends on
##           the order of the device threads from 3 terms per slot on
##   fixed   one value per site, then sumFixed
##   exact   SumX with s slots, the exact sum, and its kernel alone
## and the parts of a read: the store kernel followed by a copy of B bytes to
## the host after a wait, queued behind the kernel, and three queued copies.
##   -n: sites (list), -reps: cycles per timing, -slots: slot counts of SumX,
##   0 for the count newSumX takes by default
import qex
import backend/accel
import times, strformat

qexInit()
let ns = intSeqParam("n", @[4096, 32768, 165888, 524288])
let reps = intParam("reps", 200)
let slots = intSeqParam("slots", @[0, 64, 256, 1024, 4096])
const nSlot = 512

proc dev(n: int): ptr UncheckedArray[float] =
  cast[ptr UncheckedArray[float]](gpuMalloc(max(n, 1)*sizeof(float)))

template dotsOf(i, y, d, yy, dy: untyped) =
  ## y.y and d.y of the 6 reals of site i
  var yy = 0.0
  var dy = 0.0
  for c in 0..<6:
    let yc = y[6*i + c]
    yy += yc*yc
    dy += d[6*i + c]*yc

template timeIt(body: untyped): float =
  ## microseconds per cycle of body, after one cycle to warm up
  block:
    body
    gpuWaitAsync()
    let t0 = epochTime()
    for rep in 0..<reps: body
    gpuWaitAsync()
    1e6*(epochTime() - t0)/float(reps)

# In a proc: at the top level of a module Nim puts the statements that the
# device backends emit outside every function.
proc main =
  for n0 in ns:
    let n = n0
    let y = dev(6*n)
    let d = dev(6*n)
    onGpu(6*n):
      for k in gpuRange(6*n):
        y[k] = float((k*2654435761) and 1023)/512.0 - 1.0
        d[k] = float((k*40503 + 17) and 1023)/512.0 - 1.0
    let hb = cast[ptr UncheckedArray[float]](gpuMallocHost(2*sumHost*sizeof(float)))
    let rs = dev(2*max(n, nSlot))
    let rz0 = dev(2*nSlot)
    let ws = dev(2*(n div 15 + 16))
    var r: array[2, float]

    let tStore = timeIt:
      onGpuQueue(n):
        for i in gpuRange(n):
          dotsOf(i, y, d, yy, dy)
          rs[i] = yy
          rs[n + i] = dy
      gpuWaitAsync()

    # the two buffers of the atomic slots: a kernel zeroes those of the next one
    var ra = rs
    var rz = rz0
    onGpu(2*nSlot):
      for k in gpuRange(2*nSlot):
        ra[k] = 0.0
        rz[k] = 0.0
    var sAtomic: array[2, float]
    let tAtomic = timeIt:
      let a = ra
      let z = rz
      onGpuQueue(n):
        for i in gpuRange(n):
          dotsOf(i, y, d, yy, dy)
          if i < nSlot:
            z[i] = 0.0
            z[nSlot + i] = 0.0
          gpuAtomicAdd(a, i mod nSlot, yy)
          gpuAtomicAdd(a, nSlot + i mod nSlot, dy)
      gpuWaitAsync()
      gpuMemCpyToCpu(hb, a, 2*nSlot*sizeof(float))
      sAtomic = [0.0, 0.0]
      for k in 0..<nSlot:
        sAtomic[0] += hb[k]
        sAtomic[1] += hb[nSlot + k]
      swap(ra, rz)

    var sFixed: array[2, float]
    let tFixed = timeIt:
      onGpuQueue(n):
        for i in gpuRange(n):
          dotsOf(i, y, d, yy, dy)
          rs[i] = yy
          rs[n + i] = dy
      sumFixed(sFixed, rs, ws, hb, 2, n)

    let tKA = timeIt:  # the kernel of atomic alone
      let a = ra
      onGpuQueue(n):
        for i in gpuRange(n):
          dotsOf(i, y, d, yy, dy)
          gpuAtomicAdd(a, i mod nSlot, yy)
          gpuAtomicAdd(a, nSlot + i mod nSlot, dy)
      gpuWaitAsync()
    echo &"n {n:8}  store {tStore:8.2f}  atomic {tAtomic:8.2f} (+{tAtomic-tStore:6.2f})  kernel alone {tKA:8.2f} (+{tKA-tStore:6.2f})  fixed {tFixed:8.2f} (+{tFixed-tStore:6.2f})  us per cycle"

    const sizes = [4096, 8192, 12288, 16384, 16400, 24576, 49152, 98304, 393216]
    let hc = cast[ptr UncheckedArray[byte]](gpuMallocHost(3*sizes[^1]))
    let dc = cast[ptr UncheckedArray[byte]](gpuMalloc(3*sizes[^1]))
    template store =
      onGpuQueue(n):
        for i in gpuRange(n):
          dotsOf(i, y, d, yy, dy)
          rs[i] = yy
          rs[n + i] = dy
    for bytes in sizes:
      let tC = timeIt:
        store()
        gpuWaitAsync()
        gpuMemCpyToCpu(hc, dc, bytes)
      let tQ = timeIt:
        store()
        gpuMemCpyToCpuAsync(hc, dc, bytes)
        gpuWaitAsync()
      let tQ3 = timeIt:
        store()
        for j in 0..2: gpuMemCpyToCpuAsync(addr hc[j*bytes], addr dc[j*bytes], bytes)
        gpuWaitAsync()
      echo &"n {n:8}  store + copy of {bytes:6} B: after a wait {tC:8.2f} (+{tC-tStore:6.2f})  queued {tQ:8.2f} (+{tQ-tStore:6.2f})  three queued {tQ3:8.2f} (+{tQ3-tStore:6.2f})"
    gpuFreeHost(hc)
    gpuFree(dc)

    var seen: seq[int]
    for s0 in slots:
      var x = newSumX(2, n, s0)
      if s0 == 0 or x.s notin seen:
        seen.add x.s
        let tExact = timeIt:
          gpuSumX(x, i, n):
            dotsOf(i, y, d, yy, dy)
            [yy, dy]
          x.read(r)
        let ok = abs(r[0] - sFixed[0]) <= 1e-12*abs(sFixed[0]) and abs(r[1] - sFixed[1]) <= 1e-9*abs(sFixed[0])
        # the kernel alone: without reads the words wrap, which the read after it forgets
        let xa = x.a
        let xs = x.s
        let xst = x.st
        let xw = x.w
        let tKX = timeIt:
          onGpuQueue(n):
            for i in gpuRange(n):
              dotsOf(i, y, d, yy, dy)
              let k = 1 + (i and (xs-1))
              sumXAdd(xa, k, xst, xw, yy)
              sumXAdd(xa, k + xs, xst, xw, dy)
          gpuWaitAsync()
        var z: array[2, float]
        x.read(z)
        echo &"n {n:8}  exact, {x.s:5} slots {tExact:8.2f} (+{tExact-tStore:6.2f})  kernel alone {tKX:8.2f} (+{tKX-tStore:6.2f})  copy {sumWin*x.st*8:7} B  sums {r[0]:.17g} {r[1]:.17g}",
          (if s0 == 0: "  default" else: ""), (if ok: "" else: "  DIFFERS from fixed " & $sFixed)
      x.free
    for p in [y, d, rs, rz0, ws]: gpuFree(p)
    gpuFreeHost(hb)

main()
qexFinalize()
