import testutils
import qex, base/hyper
import comms/[halo, halogpu]
import backend/accel

# The kernels of backend/accel, which the CPU backend runs on the host:
# queued kernels run in order, sumFixed and gpuSum add in an order fixed by
# the number of terms, and GpuHaloEx brings the neighbors x+-mu of the sites
# of either parity, as stagGpu takes them: local, in another lane, or
# received from another rank.

qexInit()

proc dev(n: int): ptr UncheckedArray[float] =
  cast[ptr UncheckedArray[float]](gpuMalloc(max(n, 1)*sizeof(float)))

proc host(d: ptr UncheckedArray[float]; n: int): seq[float] =
  result = newSeq[float](n)
  if n > 0: gpuMemCpyToCpu(addr result[0], d, n*sizeof(float))

proc host(d: ptr UncheckedArray[int32]; n: int): seq[int32] =
  result = newSeq[int32](n)
  if n > 0: gpuMemCpyToCpu(addr result[0], d, n*sizeof(int32))

proc site(lo: auto; e, k: int; o: seq[int32]): int =
  ## global lex index of lane k of outer site e, shifted by o
  let xs = lo.vcoords(e)
  var x = newSeq[int](lo.nDim)
  for d in 0..<lo.nDim:
    x[d] = (int(xs[d][k]) + int(o[d]) + lo.physGeom[d]) mod lo.physGeom[d]
  lexIndex(x, lo.physGeom)

proc checkHalo(lo: auto; reps: int): float =
  ## Mismatches of the neighbors x+-mu of the sites x of each parity, local
  ## or through a GpuHaloEx of the other parity, over reps exchanges of
  ## fields with ne reals per site: real c of site y is 8 lex(y) + c + rep/2.
  const V = lo.V
  const ne = 3
  let nd = lo.nDim
  let no = lo.nSitesOuter
  let n = V*no
  let c = getDefaultComm()
  var w = newSeq[int](nd)
  for d in 0..<nd: w[d] = 1
  let hl = lo.makeHaloLayout(w, w)
  var offs: seq[seq[int32]]
  for mu in 0..<nd:
    for s in [1, -1]:
      var o = newSeq[int32](nd)
      o[mu] = int32 s
      offs.add o
  let z = newSeq[int32](nd)
  let f = dev(ne*n)
  var hf = newSeq[float](ne*n)
  for p in 0..1:  # parity of the sites sent
    let hm = hl.makeHaloMap(c, offs, p)
    let ex = newGpuHaloEx[float](hm.gather, ne, n, V, c)
    let src = hl.haloSource(hm.gather)
    let ro = host(ex.rofs, ex.nrecv)
    let rs = host(ex.rstr, ex.nrecv)
    let q = 1 - p
    let o0 = if q == 0: 0 else: lo.nEvenOuter
    let o1 = if q == 0: lo.nEvenOuter else: no
    for rep in 0..<reps:
      for e in 0..<no:
        for k in 0..<V:
          for a in 0..<ne:
            hf[(e*ne + a)*V + k] = float(8*lo.site(e, k, z) + a) + 0.5*float(rep)
      gpuMemCpyToGpu(f, addr hf[0], ne*n*sizeof(float))
      ex.pack(f)
      ex.start
      ex.wait
      let hr = host(ex.rbuf, ne*ex.nrecv)
      for e in o0..<o1:
        for mu in 0..<nd:
          for fb in 0..1:
            let x = int(if fb == 0: hl.neighborFwd[mu][e] else: hl.neighborBck[mu][e])
            for k in 0..<V:
              let j = if x < no: V*x + k else: int src[V*(x-no) + k]
              let y = lo.site(e, k, offs[2*mu + fb])
              for a in 0..<ne:
                let got = if j < 0: -1.0
                          elif j < n: hf[((j div V)*ne + a)*V + j mod V]
                          else: hr[int(ro[j-n]) + a*int(rs[j-n])]
                let want = float(8*y + a) + 0.5*float(rep)
                if got != want:
                  if result < 5:
                    echo "parity ", q, " site ", e, " lane ", k, " mu ", mu, " fb ", fb, ": ", got, " != ", want
                  result += 1
      c.barrier  # the peers read their receive buffers before the next exchange
    ex.free
  gpuFree(f)
  c.allReduce(result)

# The kernels run in procs: the GPU backends emit their code, which at the
# top level, as in a test block, lands outside any function.

proc queued(): int =
  ## mismatches after three queued kernels, each reading the output of the
  ## one before, and a synchronous one after the wait
  let n = 1000
  let a = dev(n)
  let b = dev(n)
  gpuForAsync(i, n): a[i] = float(i)
  gpuForAsync(i, n): b[i] = 2.0*a[i] + 1.0
  gpuForAsync(i, n): a[i] = b[i] - a[i]
  gpuWaitAsync()
  gpuFor(i, n): b[i] = a[i] - 1.0
  let h = host(b, n)
  for i in 0..<n:
    if h[i] != float(i): inc result
  gpuFree(a)
  gpuFree(b)

proc sums(): int =
  ## mismatches of sumFixed at the sizes where its levels change, with 1,
  ## 2 and 8 sums: exact sums of exactly representable terms, and rounded
  ## sums equal twice and near the host sum
  let hb = cast[ptr UncheckedArray[float]](gpuMallocHost(8*sumHost*sizeof(float)))
  for n in [0, 1, 15, 16, 17, 511, 512, 513, 8193]:
    for m in [1, 2, 8]:
      let a = dev(m*n)
      let w = dev(m*(n div 15 + 16))
      gpuFor(k, m*n):  # a[c n + i] = i + 1 + c/2
        let c = k div n
        a[k] = float(k - c*n + 1) + 0.5*float(c)
      var r = newSeq[float](m)
      sumFixed(r, a, w, hb, m, n)
      var bad = 0
      for c in 0..<m:
        if r[c] != float(n*(n+1) div 2) + 0.5*float(c*n): inc bad
      gpuFor(k, m*n): a[k] = 1.0/float(k + 3)
      var r1, r2 = newSeq[float](m)
      sumFixed(r1, a, w, hb, m, n)
      sumFixed(r2, a, w, hb, m, n)
      for c in 0..<m:
        var t = 0.0
        for i in 0..<n: t += 1.0/float(c*n + i + 3)
        if r1[c] != r2[c] or abs(r1[c] - t) > 1e-12*t: inc bad
      if bad != 0: echo "sumFixed n ", n, " m ", m, ": ", r, " ", r1, " ", r2
      result += bad
      gpuFree(a)
      gpuFree(w)
  gpuFreeHost(hb)

proc gsums(): int =
  ## mismatches of gpuSum of 1, 2 and 8 exact sums
  for n in [0, 1, 15, 16, 17, 511, 512, 513, 8193, 16*8193 + 5]:
    let s1 = gpuSum(i, n, 1, [float(i + 1)])
    let s2 = gpuSum(i, n, 2, [float(i), 1.0])
    let s8 = gpuSum(i, n, 8):
      var v {.noInit.}: array[8, float]
      for c in 0..<8: v[c] = float(i + c)
      v
    var bad = 0
    if s1[0] != float(n*(n+1) div 2): inc bad
    if s2[0] != float(n*(n-1) div 2) or s2[1] != float(n): inc bad
    for c in 0..<8:
      if s8[c] != float(n*(n-1) div 2 + c*n): inc bad
    if bad != 0: echo "gpuSum n ", n, ": ", s1, " ", s2, " ", s8
    result += bad

suite "backend/accel":
  test "queued kernels run in order":
    check queued() == 0

  test "sumFixed":
    check sums() == 0

  test "gpuSum":
    check gsums() == 0

  test "GpuHaloEx neighbors":
    for ll in [[4,4,4,8], [8,4,4,4], [4,8,8,8]]:
      let lo = latticeFromLocalLattice(ll, nRanks).newLayout
      echo "lattice ", lo.physGeom, "  ranks ", lo.rankGeom, "  lanes ", lo.innerGeom
      check checkHalo(lo, 3) == 0

qexFinalize()
