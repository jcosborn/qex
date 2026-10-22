import testutils
import qex
import backend/[accel, cpugpu, cgfield]
import comms/[halo, halogpu]
import sequtils

# backend/gpumem keeps the device copies of host objects in a table keyed by
# their host addresses.  The entries stay where the handles taken from them
# point while other entries come and go, and the entry of freed host memory
# does not pass to a new object at the same address.  tgpumem.nims gives the
# CPU backend device copies apart from host memory, a collection at every
# allocation and freed memory filled with 0xFF.

qexInit()

proc dev[T](p: pointer; n: int): seq[T] =
  ## the n values of type T at device address p
  result = newSeq[T](n)
  if n > 0: gpuMemCpyToCpu(addr result[0], p, n*sizeof(T))

proc host(p: pointer; n: int): seq[byte] =
  result = newSeq[byte](n)
  copyMem(addr result[0], p, n)

proc handles[T](cf: var T): int =
  ## Handles of 300 entries, the first ones taken in a nearly empty table, kept
  ## while every other entry is freed: the handles that moved, copies through
  ## them that moved the wrong bytes, and the CgField's entry.
  const n = 300
  var hs = newSeq[seq[float]](n)
  var ps = newSeq[ptr GpuMem](n)
  for i in 0..<n:
    hs[i] = newSeq[float](i+1)
    for j in 0..i: hs[i][j] = float(1000*i + j)
    ps[i] = getGpuMem(addr hs[i][0], (i+1)*sizeof(float))
  for i in countup(1, n-1, 2): freeGpuMem(addr hs[i][0])
  inc kernelCallCount
  for i in countup(0, n-1, 2):
    let p = getGpuMem(addr hs[i][0])
    if p != ps[i]:
      inc result
      continue
    ps[i].copyIn(addr hs[i][0])
    if dev[float](p.p, i+1) != hs[i]: inc result
    var t = newSeq[float](i+1)
    for j in 0..i: t[j] = -float(j)
    gpuMemCpyToGpu(p.p, addr t[0], (i+1)*sizeof(float))
    ps[i].copyOut(addr hs[i][0])
    if hs[i] != t: inc result
  if getGpuMem(addr cf.cpu[0]) != cf.gpu.pgm:
    inc result
  else:
    cf.copyToGpu
    if getGpuMem(addr cf.cpu[0]).copyInCount != 1: inc result
  for i in countup(0, n-1, 2): freeGpuMem(addr hs[i][0])

proc nested(): int =
  ## A seq of 1000 seqs, more than the entries the test before leaves room
  ## for, uploads its inner seqs, new entries, holding the handle of its
  ## own: the copy of the container counts in its entry, and the device
  ## container holds the inner device pointers, whose data match.
  const n = 1000
  var c = newSeq[seq[int32]](n)
  for i in 0..<n: c[i] = @[int32 i, int32 2*i, int32 3*i]
  inc kernelCallCount
  var g: GpuSeq[ptr UncheckedArray[int32]]
  g.toGpu(c)
  if getGpuMem(addr c[0]).copyInCount != 1: inc result
  let t = dev[ptr UncheckedArray[int32]](g.p, n)
  for i in 0..<n:
    let e = getGpuMem(addr c[i][0])
    if t[i] != e.p or dev[int32](e.p, 3) != c[i]: inc result
  for i in 0..<n: freeGpuMem(addr c[i][0])
  freeGpuMem(addr c[0])

proc storage(k: int): (pointer, int) {.noinline.} =
  ## new storage of 1000 reals 1000 k + i, the first one copied in once:
  ## its address and the mismatches of its copies in and out
  var a: alignedMem[float]
  a.new(1000)
  for i in 0..<a.len: a[i] = float(1000*k + i)
  if k == 0: gpuMemFlagsIncl(addr a[0], {gmCpuWriteOnce})
  inc kernelCallCount
  var g: GpuSeq[float]
  g.toGpu(a)
  var bad = 0
  let d = dev[float](g.p, a.len)
  for i in 0..<a.len:
    if d[i] != a[i]: inc bad
  var t = newSeq[float](a.len)  # as a kernel would change them
  for i in 0..<a.len: t[i] = -float(1000*k + i)
  gpuMemCpyToGpu(g.p, addr t[0], a.len*sizeof(float))
  a.fromGpu
  for i in 0..<a.len:
    if a[i] != t[i]: inc bad
  (pointer(addr a[0]), bad)

proc field(lo: Layout; k: int): (pointer, int) {.noinline.} =
  ## a new field of k+1, the first one copied in once and never out: its
  ## address and the copies in and, after the first, out that differ
  var f = lo.ColorVector()
  threads: f := float(k+1)
  if k == 0:
    f.gpuFlagsIncl {gmCpuWriteOnce}
    f.gpuFlagsExcl {gmGpuWrite}
  inc kernelCallCount
  let g = toGpu(f)
  let n = g.bytes
  var bad = int(dev[byte](g.p, n) != host(addr f[0], n))
  if k > 0:
    var t = newSeq[byte](n)  # as a kernel would change them
    for i in 0..<n: t[i] = byte((i + k) mod 256)
    gpuMemCpyToGpu(g.p, addr t[0], n)
    f.fromGpu
    bad += int(host(addr f[0], n) != t)
  (pointer(addr f[0]), bad)

proc layout(lo: Layout; fwd, bck: openArray[int]): (pointer, int) {.noinline.} =
  ## a new halo layout: the address of its seq of forward neighbor tables
  ## and the mismatches of the tables on the device
  let hl = lo.makeHaloLayout(fwd, bck)
  inc kernelCallCount
  let g = toGpu(hl)
  var bad = 0
  for (d, h) in [(g.neighborFwd, hl.neighborFwd), (g.neighborBck, hl.neighborBck)]:
    let t = dev[ptr UncheckedArray[int32]](d.p, h.len)
    for mu in 0..<h.len:
      if dev[int32](t[mu], h[mu].len) != h[mu]: inc bad
  (pointer(addr hl.neighborFwd[0]), bad)

proc halos(lo: Layout; hl: HaloLayout): int =
  ## A seq of two halos of fields 1 and 2, whose second halo becomes one of
  ## a field 3 after the first upload: the device fields of the device seq
  ## that differ from those of its halos, in both uploads.
  var fs: array[3, typeof(lo.ColorVector())]
  for i in 0..2:
    fs[i] = lo.ColorVector()
    threads: fs[i] := float(i+1)
  var h = @[makeHalo(hl, fs[0]), makeHalo(hl, fs[1])]
  for rep in 0..1:
    inc kernelCallCount
    let g = toGpu(h)
    let t = dev[typeof(g[0])](g.p, 2)
    for i in 0..1:
      let f = h[i].field
      if pointer(t[i].field.p) != getGpuMem(addr f[0]).p or
         dev[byte](t[i].field.p, t[i].field.bytes) != host(addr f[0], t[i].field.bytes): inc result
    h[1] = makeHalo(hl, fs[2])

proc reals(n, k: int): (pointer, int) {.noinline.} =
  ## a new seq of n reals 1000 k + i: its address and whether its device
  ## copy differs
  var s = newSeq[float](n)
  for i in 0..<n: s[i] = float(1000*k + i)
  inc kernelCallCount
  let g = toGpu(s)
  (pointer(addr s[0]), int(dev[float](g.p, n) != s))

template reuse(name: string; body: untyped) =
  ## runs body, a call giving an address and mismatches, four times with a
  ## full collection after each, and prints how many addresses were reused
  var ps: seq[pointer]
  var bs: seq[int]
  for k {.inject.} in 0..<4:
    let (p, b) = body
    ps.add p
    bs.add b
    GC_fullCollect()
  echo name, ": ", ps.len - ps.deduplicate.len, " of ", ps.len - 1, " addresses reused, mismatches ", bs
  check bs == @[0, 0, 0, 0]

let lo = latticeFromLocalLattice([4,4,4,8], nRanks).newLayout

suite "backend/gpumem":
  test "handles of entries":
    var cf = newCgField(lo.ColorVector())
    check handles(cf) == 0

  test "nested copies":
    check nested() == 0

  test "freed storage":
    reuse("alignedMem"): storage(k)

  test "freed fields":
    reuse("Field"): field(lo, k)

  test "dropped halo layouts":
    reuse("HaloLayout"):
      if k mod 2 == 0: layout(lo, [1,1,1,1], [0,0,0,0])
      else: layout(lo, [0,0,0,0], [1,1,1,1])

  test "seqs of halos":
    check halos(lo, lo.makeHaloLayout([1,1,1,1], [1,1,1,1])) == 0

  test "freed seqs of another length":
    reuse("seq"): reals(100 - k mod 2, k)

qexFinalize()
