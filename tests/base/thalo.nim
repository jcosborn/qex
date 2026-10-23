import testutils
import qex, base/hyper, comms/[gather, halo]

# One pending local message exercises gather buffers with any QMP backend.
type
  LoopbackComm = ref object of Comm
    send, recv: pointer
    bytes: int
  GatherInts = object
    src, dest: ptr UncheckedArray[int]
    elemsize, vlen: int

method nsends(c: LoopbackComm): int = ord(c.send != nil)
method nrecvs(c: LoopbackComm): int = ord(c.recv != nil)
method pushRecv(c: LoopbackComm, rank: int, p: pointer, bytes: int) =
  doAssert rank == 0 and c.recv == nil
  c.recv = p
  c.bytes = bytes
method pushSend(c: LoopbackComm, rank: int, p: pointer, bytes: int) =
  doAssert rank == 0 and c.send == nil and c.recv != nil and bytes == c.bytes
  c.send = p
method waitRecvs(c: LoopbackComm, i0, i1: int, free = true) =
  doAssert i0 == 0 and i1 == 0 and c.send != nil and c.recv != nil
  copyMem(c.recv, c.send, c.bytes)
  if free: c.recv = nil
method freeRecvs(c: LoopbackComm, i0, i1: int) =
  doAssert i0 == 0 and i1 == 0 and c.recv != nil
  c.recv = nil
method waitSends(c: LoopbackComm, i0, i1: int) =
  doAssert i0 == 0 and i1 == 0 and c.send != nil
  c.send = nil

template copy(data: GatherInts, d, s: pointer|SomeInteger) =
  when d is pointer: cast[ptr int](d)[] = data.src[s]
  elif s is pointer: data.dest[d] = cast[ptr int](s)[]
  else: data.dest[d] = data.src[s]

proc exerciseGatherBuffers(): int {.noinline.} =
  # Return before collecting: refc conservatively retains locals on live stacks.
  let before = getOccupiedSharedMem()
  let c = LoopbackComm()
  var send, recv: ptr UncheckedArray[char]
  var capacity = 0
  for count in [4, 512, 8]:
    var src, dest = newSeq[int](count)
    let gm = GatherMap(nsrc: count, ndest: count,
      sidx: newSeq[int32](count), rdest: newSeq[int32](count),
      smsginfo: @[MsgInfo(rank: 0, count: count.int32)],
      rmsginfo: @[MsgInfo(rank: 0, count: count.int32)])
    for i in 0..<count:
      src[i] = 1000*count + i
      gm.sidx[i] = int32(count-1-i)
      gm.rdest[i] = int32(i)
    let data = GatherInts(src: cast[ptr UncheckedArray[int]](src[0].addr),
      dest: cast[ptr UncheckedArray[int]](dest[0].addr), elemsize: sizeof(int), vlen: 1)
    c.gather(gm, data, keep = true)
    doAssert c.send == nil and c.recv == nil
    for i in 0..<count: doAssert dest[i] == src[count-1-i]
    capacity = max(capacity, count*sizeof(int))
    doAssert c.gatherBuffers.send.n == capacity
    doAssert c.gatherBuffers.recv.n == capacity
    if count == 8:
      doAssert c.gatherBuffers.send.p == send
      doAssert c.gatherBuffers.recv.p == recv
    send = c.gatherBuffers.send.p
    recv = c.gatherBuffers.recv.p
  getOccupiedSharedMem()-before

# global lex index of lane k of outer site e, shifted by o and wrapped
proc site(lo: auto, e, k: int, o: seq[int32]): int =
  let xs = lo.vcoords(e)
  var x = newSeq[int](lo.nDim)
  for d in 0..<lo.nDim:
    x[d] = (int(xs[d][k]) + int(o[d]) + lo.physGeom[d]) mod lo.physGeom[d]
  lexIndex(x, lo.physGeom)

# compose unit steps on the extended lattice, -1 once it leaves the box
proc nbr(hl: auto, i: int, o: seq[int32]): int =
  result = i
  for d in 0..<o.len:
    for s in 1..abs(o[d]):
      if result < 0: return
      result = if o[d] > 0: hl.neighborFwd[d][result] else: hl.neighborBck[d][result]

# [0,0].re of every lane holds the site's global lex index
proc fill(g: auto) =
  let lo = g[0].l
  const V = lo.V
  let z = newSeq[int32](lo.nDim)
  for mu in 0..<g.len:
    for e in 0..<lo.nSitesOuter:
      var t: array[V,float]
      for k in 0..<V: t[k] = float lo.site(e, k, z)
      g[mu][e] := 0
      g[mu][e][0,0].re := t

# number of (site, offset, lane) mismatches between halo values and coordinates
proc checkUpdate(g: auto, hl: auto, hm: auto, offsets: seq[seq[int32]]): int =
  let lo = g[0].l
  const V = lo.V
  let comm = getDefaultComm()
  fill g
  for mu in 0..<g.len:
    let h = makeHalo(hl, g[mu])
    h.update hm, comm
    for i in 0..<hl.nOut:
      for o in offsets:
        let j = hl.nbr(i, o)
        if j < 0:
          inc result
          continue
        for k in 0..<V:
          let want = float lo.site(i, k, o)
          let got = h[j][0,0].re[k]
          if got != want:
            if result < 5: echo "mu ", mu, " site ", i, " offset ", o, " lane ", k, ": ", got, " != ", want
            inc result

# number of lane mismatches after accumulating 1 from every mapped halo cell
proc checkRev(g: auto, hl: auto, hm: auto, offsets: seq[seq[int32]]): int =
  let lo = g[0].l
  const V = lo.V
  let comm = getDefaultComm()
  let z = newSeq[int32](lo.nDim)
  # expected count per global site: each halo cell the offsets reach, once per
  # lane, summed over ranks since every rank's halo sends back to the owner
  var tally = newSeq[float](lo.physVol)
  var seen = newSeq[bool](hl.nExt - hl.nOut)
  for i in 0..<hl.nOut:
    for o in offsets:
      let j = hl.nbr(i, o)
      if j < hl.nOut or seen[j-hl.nOut]: continue
      seen[j-hl.nOut] = true
      for k in 0..<V: tally[lo.site(i, k, o)] += 1.0
  rankSum tally
  for mu in 0..<g.len:
    g[mu] := 0
    let h = makeHalo(hl, g[mu])
    for c in 0..<h.halo.len: h.halo[c][0,0].re := 1.0
    h.updateRev hm, comm
    for e in 0..<lo.nSitesOuter:
      for k in 0..<V:
        let want = tally[lo.site(e, k, z)]
        let got = g[mu][e][0,0].re[k]
        if got != want:
          if result < 5: echo "mu ", mu, " site ", e, " lane ", k, ": ", got, " != ", want
          inc result

# per plane sum_x Re tr(U_mu(x) U_nu(x+mu) U_mu(x+nu)^+ U_nu(x)^+) / (nc physVol) = 6 plaq
proc haloPlaq(g: auto, hl: auto, hm: auto): seq[float] =
  let lo = g[0].l
  let nd = lo.nDim
  let comm = getDefaultComm()
  type H = type makeHalo(hl, g[0])
  var h = newSeq[H](nd)
  for mu in 0..<nd:
    h[mu] = makeHalo(hl, g[mu])
    h[mu].update hm, comm
  result = newSeq[float]((nd*(nd-1)) div 2)
  for i in 0..<hl.nOut:
    var k = 0
    for mu in 1..<nd:
      let n0 = hl.neighborFwd[mu][i]
      for nu in 0..<mu:
        let n1 = hl.neighborFwd[nu][i]
        let a = g[mu][i] * h[nu][n0]
        let b = g[nu][i] * h[mu][n1]
        result[k] += simdReduce redot(a, b)
        inc k
  rankSum result
  let vf = 1.0/(g[0][0].nrows*lo.physVol)
  for k in 0..<result.len: result[k] *= vf

qexInit()

suite "Halo":
  echo "rank ", myRank, "/", nRanks
  var (lo, g, r) = setupLattice([8,8,4,4])
  let nd = lo.nDim
  let comm = getDefaultComm()
  let hl = haloLayout(lo, [1,1,1,1], [1,1,1,1])
  echo "V: ", lo.V, " outerGeom: ", lo.outerGeom, " innerGeom: ", lo.innerGeom,
    " nOut: ", hl.nOut, " nExt: ", hl.nExt
  # offset sets: the 16 corners (+-1)^4 and the 8 unit axis steps
  var corners = newSeq[seq[int32]]()
  for c in 0..15:
    var t = newSeq[int32](nd)
    for d in 0..<nd: t[d] = if ((c shr d) and 1) == 1: -1 else: 1
    corners.add t
  var axes = newSeq[seq[int32]]()
  for d in 0..<nd:
    for s in [-1'i32, 1'i32]:
      var t = newSeq[int32](nd)
      t[d] = s
      axes.add t
  let hmC = haloMap(hl, comm, corners)
  let hmA = haloMap(hl, comm, axes)

  test "map indices":
    # gather: local sites -> halo cells; gatherRev: halo cells -> local sites
    proc inRange(x: seq[int32], n: int): bool =
      for v in x:
        if v < 0 or v >= n: return false
      true
    let nh = (hl.nExt - hl.nOut) * lo.V
    for hm in [hmC, hmA]:
      check(hm.sites.sidx.inRange(lo.nSitesOuter) and hm.sites.lidx.inRange(lo.nSitesOuter))
      check(hm.sites.rdest.inRange(hl.nExt-hl.nOut) and hm.sites.ldest.inRange(hl.nExt-hl.nOut))
      check(hm.gather.sidx.inRange(lo.nSites) and hm.gather.lidx.inRange(lo.nSites))
      check(hm.gather.rdest.inRange(nh) and hm.gather.ldest.inRange(nh))
      check(hm.gatherRev.sidx.inRange(nh) and hm.gatherRev.lidx.inRange(nh))
      check(hm.gatherRev.rdest.inRange(lo.nSites) and hm.gatherRev.ldest.inRange(lo.nSites))

  test "retained gather buffers grow, reuse storage, and release with their owner":
    discard exerciseGatherBuffers()  # Warm threading and profiling state.
    GC_fullCollect()
    GC_fullCollect()
    let before = getOccupiedSharedMem()
    let allocated = exerciseGatherBuffers()
    check allocated >= 2*512*sizeof(int)
    GC_fullCollect()
    GC_fullCollect()
    check getOccupiedSharedMem() == before

  test "retained gather buffers belong to their communicator":
    fill g
    let hm = haloMap(hl, comm, axes)
    let h = makeHalo(hl, g[0])
    h.update(hm, comm)
    let first = comm.gatherBuffers
    check first != nil
    let send = first.send.p
    let recv = first.recv.p
    h.update(hm, comm)
    check first.send.p == send
    check first.recv.p == recv
    let other = getCommQmp()
    h.update(hm, other)
    check other.gatherBuffers != first
    check other.gatherBuffers.send.p != send
    check other.gatherBuffers.recv.p != recv
    other.gatherBuffers.releaseGatherBuffers
    other.gatherBuffers.releaseGatherBuffers
    check other.gatherBuffers.send.p == nil
    check other.gatherBuffers.recv.p == nil
    check other.gatherBuffers.send.n == 0
    check other.gatherBuffers.recv.n == 0
    h.update(hm, other)
    check other.gatherBuffers.send.p != nil
    check other.gatherBuffers.recv.p != nil
    other.gatherBuffers.releaseGatherBuffers
    check first.send.p == send
    check first.recv.p == recv

  test "update corners":
    check(checkUpdate(g, hl, hmC, corners) == 0)
  test "update axes":
    check(checkUpdate(g, hl, hmA, axes) == 0)
  test "updateRev corners":
    check(checkRev(g, hl, hmC, corners) == 0)
  test "updateRev axes":
    check(checkRev(g, hl, hmA, axes) == 0)

  test "update at other lane counts":
    # V = 4 leaves directions without lanes, so split ranks send whole sites.
    template layoutCase(lat: openArray[int], v: static int) =
      block:
        let lo2 = newLayout(lat, v)
        let n = lo2.nDim
        var g2 = newSeq[type(lo2.ColorMatrix(3))](n)
        for mu in 0..<n: g2[mu] = lo2.ColorMatrix(3)
        var ones = newSeq[int](n)
        for d in 0..<n: ones[d] = 1
        let hl2 = haloLayout(lo2, ones, ones)
        var cs, ax: seq[seq[int32]]
        for c in 0..<(1 shl n):
          var t = newSeq[int32](n)
          for d in 0..<n: t[d] = if ((c shr d) and 1) == 1: -1 else: 1
          cs.add t
        for d in 0..<n:
          for s in [-1'i32, 1'i32]:
            var t = newSeq[int32](n)
            t[d] = s
            ax.add t
        echo "V: ", lo2.V, " innerGeom: ", lo2.innerGeom, " rankGeom: ", lo2.rankGeom
        check(checkUpdate(g2, hl2, haloMap(hl2, comm, cs), cs) == 0)
        check(checkUpdate(g2, hl2, haloMap(hl2, comm, ax), ax) == 0)
        check(checkRev(g2, hl2, haloMap(hl2, comm, cs), cs) == 0)
    layoutCase(lo.physGeom, 4)

  test "plaquette":
    threads: g.random r
    let p = g.plaq
    let hp = haloPlaq(g, hl, hmA)
    for k in 0..<p.len: check(hp[k] ~ 6.0*p[k])

  test "offset layout and index table":
    # Widths from the extreme offsets: forward 2 along d0, backward 1 along d1.
    let offs = @[@[2'i32,0,0,0], @[0'i32,-1,0,0], @[1'i32,-1,0,0], @[0'i32,0,0,0]]
    let hlo = haloLayout(lo, offs)
    check(hlo == haloLayout(lo, [2,0,0,0], [0,1,0,0]))
    let idx = haloIndex(hlo, offs)
    check(idx.len == hlo.nOut*offs.len)
    var bad = 0
    for i in 0..<hlo.nOut:
      for t, o in offs:
        if int(idx[i*offs.len+t]) != hlo.nbr(i, o): inc bad
    check(bad == 0)
    expect ValueError: discard haloLayout(lo, @[@[8'i32,0,0,0]])
    expect ValueError: discard haloIndex(hlo, @[@[3'i32,0,0,0]])
    expect ValueError: discard haloIndex(hlo, @[@[0'i32,-2,0,0]])

  test "cache":
    check(haloLayout(lo, [1,1,1,1], [1,1,1,1]) == hl)
    check(haloLayout(lo, [1'i32,1,1,1], [1'i32,1,1,1]) == hl)
    check(haloLayout(lo, [1,1,1,1], [0,0,0,0]) != hl)
    check(haloMap(hl, comm, corners) == hmC)
    check(haloMap(hl, comm, axes) == hmA)
    check(hmC != hmA)
    var ca = newSeq[array[4,int32]]()
    for o in corners: ca.add [o[0], o[1], o[2], o[3]]
    check(haloMap(hl, comm, ca) == hmC)

qexFinalize()
