import testutils
import qex, base/hyper, comms/halo

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

  test "update corners":
    check(checkUpdate(g, hl, hmC, corners) == 0)
  test "update axes":
    check(checkUpdate(g, hl, hmA, axes) == 0)
  test "updateRev corners":
    check(checkRev(g, hl, hmC, corners) == 0)
  test "updateRev axes":
    check(checkRev(g, hl, hmA, axes) == 0)

  test "plaquette":
    threads: g.random r
    let p = g.plaq
    let hp = haloPlaq(g, hl, hmA)
    for k in 0..<p.len: check(hp[k] ~ 6.0*p[k])

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
