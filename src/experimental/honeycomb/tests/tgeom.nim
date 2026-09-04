## Task G test: validate `hcgeom.nim` against doc/FORMULATION.md sections 1-2 and
## the acceptance list in doc/PLAN.md (task G).
##
##   export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
##   cd build_mac && make src/experimental/honeycomb/tests/tgeom.nim && ./bin/tgeom
##
## Pure combinatorics: no QEX Layout, no MPI, no threads.

import std/[math, tables, random, strformat, unittest]
import ../hcgeom

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

template ok(msg: string, cond: bool) =
  ## `check` plus a human readable PASS/FAIL line carrying the measured number
  let c = cond
  if c: echo "PASS: ", msg
  else: echo "FAIL: ", msg
  check c

# ---------------------------------------------------------------------------
# helpers that re-derive the geometry independently of hcgeom.walk/step
# ---------------------------------------------------------------------------

proc pos2(s: Site): Vec2 =
  ## doubled coordinates of a site:  A(y) -> 2y,  B(y) -> 2y + (1,1,1,1)
  for mu in 0..<nDim:
    result[mu] = 2*s.cell[mu] + (if s.sub == 1: 1 else: 0)

proc linkEnds(l: LinkRef): tuple[src, dst: Site] =
  ## the two endpoints of one link traversal, honouring `dag`.
  ## Derived straight from the table in FORMULATION 1.4, not from hcgeom.step.
  var a, b: Site
  case l.kind
  of lkA:
    a = Site(cell: l.cell, sub: 0)
    b = a
    b.cell[l.idx] += 1
  of lkB:
    a = Site(cell: l.cell, sub: 1)
    b = a
    b.cell[l.idx] += 1
  of lkD:
    a = Site(cell: l.cell, sub: 1)
    b = Site(cell: l.cell, sub: 0)
    for mu in 0..<nDim: b.cell[mu] += (l.idx shr mu) and 1
  result = if l.dag: (b, a) else: (a, b)

proc linkDisp(l: LinkRef): Vec2 =
  let (s, d) = linkEnds(l)
  sub2(pos2(d), pos2(s))

proc walkPath(p: Path, start: Site): tuple[closes, connected: bool, dest: Site] =
  ## follow the links of `p` one at a time and report whether they chain up
  var s = start
  var conn = true
  for l in p:
    let (a, b) = linkEnds(l)
    if a != s: conn = false
    s = b
  (s == start, conn, s)

proc isAxisVec(v: Vec2): bool =
  var nz = 0
  for mu in 0..<nDim:
    if v[mu] != 0:
      inc nz
      if abs(v[mu]) != 2: return false
  nz == 1

proc isDiagVec(v: Vec2): bool =
  for mu in 0..<nDim:
    if abs(v[mu]) != 1: return false
  true

proc wrap(c: Cell, L: int): Cell =
  for mu in 0..<nDim: result[mu] = ((c[mu] mod L) + L) mod L

proc key(l: LinkRef, L: int): (int, int, int, int, int, int) =
  let c = wrap(l.cell, L)
  (l.kind.ord, l.idx, c[0], c[1], c[2], c[3])

# ---------------------------------------------------------------------------

suite "hcgeom: 16-cell honeycomb geometry":

  test "1. 24 unit neighbour directions":
    ok(&"nDirs = {nDirs}", nDirs == 24)
    ok(&"nAxis = {nAxis}, nDiag = {nDiag}", nAxis == 8 and nDiag == 16)
    var allUnit = true
    for d in 0..<nDirs:
      if dot2(allDirs[d], allDirs[d]) != 4: allUnit = false
    ok("all 24 have |n| = 1  (dot2 = 4)", allUnit)
    # all distinct, and closed under negation
    var distinct1 = true
    for a in 0..<nDirs:
      for b in (a+1)..<nDirs:
        if allDirs[a] == allDirs[b]: distinct1 = false
    ok("all 24 distinct", distinct1)
    var negClosed = true
    for d in 0..<nDirs:
      if allDirs[opposite(d)] != neg2(allDirs[d]): negClosed = false
    ok("opposite() is consistent with -n", negClosed)
    var axisOk = true
    for d in 0..<nAxis:
      if not isAxisVec(allDirs[d]): axisOk = false
    for d in nAxis..<nDirs:
      if not isDiagVec(allDirs[d]): axisOk = false
    ok("first 8 are axis, last 16 are diagonal", axisOk)

  test "1b. step returns the directed link between neighbouring sites":
    for c in [Cell([0, 0, 0, 0]), Cell([2, -1, 3, -4])]:
      for sub in 0..1:
        let s = Site(cell: c, sub: sub)
        for d in 0..<nDirs:
          let (l, dst) = step(s, d)
          let (a, b) = linkEnds(l)
          checkpoint(&"cell={c} sub={sub} dir={d}")
          check a == s
          check b == dst
          check sub2(pos2(dst), pos2(s)) == dirVec(d)

  test "2. sum_i n_mu n_nu = 6 delta_munu  (2 delta_munu for the 8 axis)":
    # everything is doubled, so the doubled sum is 4x the true one
    var s24, s8: array[nDim, array[nDim, int]]
    for d in 0..<nDirs:
      for a in 0..<nDim:
        for b in 0..<nDim:
          s24[a][b] += allDirs[d][a]*allDirs[d][b]
          if d < nAxis: s8[a][b] += allDirs[d][a]*allDirs[d][b]
    var g24 = true
    var g8 = true
    for a in 0..<nDim:
      for b in 0..<nDim:
        let e = if a == b: 1 else: 0
        if s24[a][b] != 4*6*e: g24 = false
        if s8[a][b] != 4*2*e: g8 = false
    echo &"  sum_24 n_0 n_0 = {0.25*s24[0][0].float}, sum_8 n_0 n_0 = {0.25*s8[0][0].float}"
    ok("sum over 24 = 6 delta_munu", g24)
    ok("sum over 8 axis = 2 delta_munu", g8)

  test "3. 32 unordered zero-sum triples, each 1 axis + 2 diagonals":
    var nTrip = 0
    var shapeOk = true
    for i in 0..<nDirs:
      for j in (i+1)..<nDirs:
        for k in (j+1)..<nDirs:
          let s = add2(add2(allDirs[i], allDirs[j]), allDirs[k])
          if s == [0, 0, 0, 0]:
            inc nTrip
            var na = 0
            var nd = 0
            for t in [i, j, k]:
              if isAxis(t): inc na else: inc nd
            if na != 1 or nd != 2: shapeOk = false
    echo &"  zero-sum triples = {nTrip}"
    ok(&"32 unordered zero-sum triples (got {nTrip})", nTrip == 32)
    ok("every triple = 1 axis + 2 diagonal edges", shapeOk)

  test "4. 96 triangles through a site, 32 with the apex there":
    var thru = 0
    var apexPairs = initTable[(int, int), int]()
    for u in 0..<nDirs:
      for v in (u+1)..<nDirs:
        let w = sub2(allDirs[v], allDirs[u])
        if dirIndexOf(w) >= 0:
          inc thru
          if not isAxis(u) and not isAxis(v):
            apexPairs[(u, v)] = 0
    echo &"  triangles through a site = {thru} ({2*thru} oriented)"
    ok(&"96 triangles through a site (got {thru})", thru == 96)
    ok(&"32 of them have their apex there (got {apexPairs.len})",
       apexPairs.len == 32)
    # apexTris must enumerate exactly those 32, with no duplicates
    ok(&"apexTris.len = {apexTris.len}", apexTris.len == nTriPerSite)
    var seen = initTable[(int, int), int]()
    var inSet = true
    for t in apexTris:
      var a = diagIndex(t.delta)
      var b = diagIndex(t.deltaP)
      if a > b: swap(a, b)
      if not apexPairs.hasKey((a, b)): inSet = false
      seen[(a, b)] = seen.getOrDefault((a, b)) + 1
    var noDup = true
    for _, c in seen:
      if c != 1: noDup = false
    ok("apexTris are all apex pairs of the site", inSet)
    ok(&"apexTris has no duplicates ({seen.len} distinct)",
       noDup and seen.len == 32)
    # label sanity: bit mu of delta clear, deltaP = delta or 2^mu
    var lab = true
    for t in apexTris:
      if ((t.delta shr t.mu) and 1) != 0: lab = false
      if t.deltaP != (t.delta or (1 shl t.mu)): lab = false
    ok("apex labels (delta, mu) canonical: bit mu of delta clear", lab)

  test "5. every link is in exactly 8 triangles (4^4 cell torus)":
    const L = 4
    var tally = initTable[(int, int, int, int, int, int), int]()
    var nTri = 0
    for x0 in 0..<L:
      for x1 in 0..<L:
        for x2 in 0..<L:
          for x3 in 0..<L:
            let y: Cell = [x0, x1, x2, x3]
            for sub in 0..1:
              for t in apexTris:
                let p = triPath(Site(cell: y, sub: sub), t)
                inc nTri
                for l in p:
                  let k = key(l, L)
                  tally[k] = tally.getOrDefault(k) + 1
    let nCells = L*L*L*L
    echo &"  cells = {nCells}, triangles = {nTri} ({nTri/nCells} per cell)"
    echo &"  distinct links = {tally.len} (expected {24*nCells})"
    ok(&"64 triangles per cell (got {nTri div nCells})", nTri == 64*nCells)
    ok(&"all {24*nCells} links appear", tally.len == 24*nCells)
    var mult8 = true
    var minm = high(int)
    var maxm = 0
    for _, c in tally:
      if c != 8: mult8 = false
      minm = min(minm, c)
      maxm = max(maxm, c)
    echo &"  link multiplicity: min {minm}, max {maxm}"
    ok("every link is contained in exactly 8 triangles", mult8)

  test "6. triPath closes and matches FORMULATION 2.2":
    var closes = true
    var connected = true
    var dispOk = true
    var lenOk = true
    let y: Cell = [3, -1, 2, 0]
    for sub in 0..1:
      let apex = Site(cell: y, sub: sub)
      for t in apexTris:
        let p = triPath(apex, t)
        if p.len != 3: lenOk = false
        let w = walkPath(p, apex)
        if not w.closes: closes = false
        if not w.connected: connected = false
        var s: Vec2
        var na = 0
        var nd = 0
        for l in p:
          let d = linkDisp(l)
          s = add2(s, d)
          if isAxisVec(d): inc na
          elif isDiagVec(d): inc nd
        if s != [0, 0, 0, 0]: dispOk = false
        if na != 1 or nd != 2: dispOk = false
    ok("every triPath has 3 links", lenOk)
    ok("links chain head-to-tail", connected)
    ok("64 triPaths (2 sublattices x 32) all close", closes)
    ok("displacements sum to zero: 1 axis + 2 diagonals", dispOk)

    # explicit comparison with the hand-derived forms of FORMULATION 2.2
    var formA = true
    var formB = true
    for t in apexTris:
      let
        db = t.delta xor 15
        dbp = t.deltaP xor 15          # = db - 2^mu
        # apex A(y): uD[db](y-db)^dag . uB[mu](y-db) . uD[db'](y-db')
        cdb = block:
          var c = y
          for mu in 0..<nDim: c[mu] -= (db shr mu) and 1
          c
        cdbp = block:
          var c = y
          for mu in 0..<nDim: c[mu] -= (dbp shr mu) and 1
          c
        wantA: Path = @[
          LinkRef(kind: lkD, idx: db, cell: cdb, dag: true),
          LinkRef(kind: lkB, idx: t.mu, cell: cdb, dag: false),
          LinkRef(kind: lkD, idx: dbp, cell: cdbp, dag: false)]
        # apex B(y): uD[delta](y) . uA[mu](y+delta) . uD[delta'](y)^dag
        cd = block:
          var c = y
          for mu in 0..<nDim: c[mu] += (t.delta shr mu) and 1
          c
        wantB: Path = @[
          LinkRef(kind: lkD, idx: t.delta, cell: y, dag: false),
          LinkRef(kind: lkA, idx: t.mu, cell: cd, dag: false),
          LinkRef(kind: lkD, idx: t.deltaP, cell: y, dag: true)]
      if dbp != db - (1 shl t.mu): formA = false
      if triPath(Site(cell: y, sub: 0), t) != wantA: formA = false
      if triPath(Site(cell: y, sub: 1), t) != wantB: formB = false
    ok("apex A path = uD[db](y-db)^dag uB[mu](y-db) uD[db'](y-db')", formA)
    ok("apex B path = uD[d](y) uA[mu](y+d) uD[d'](y)^dag", formB)

  test "7. 16 hexagons, coplanar, 60 degrees, 2 apex triangles each":
    ok(&"hexagons.len = {hexagons.len}", hexagons.len == nHexPerSite)
    var unit = true
    var cons = true
    var opp = true
    var distinctRing = true
    for h in hexagons:
      var s = initTable[int, int]()
      for k in 0..<6:
        let
          a = allDirs[h.ring[k]]
          b = allDirs[h.ring[(k+1) mod 6]]
          c = allDirs[h.ring[(k+3) mod 6]]
        if dot2(a, a) != 4: unit = false
        if dot2(a, b) != 2: cons = false      # cos 60 = 1/2  -> doubled dot 2
        if dot2(a, c) != -4: opp = false      # opposite
        s[h.ring[k]] = 0
      if s.len != 6: distinctRing = false
    ok("all hexagon vectors are unit", unit)
    ok("consecutive dot products = 1/2 (60 degrees)", cons)
    ok("opposite dot products = -1", opp)
    ok("each hexagon has 6 distinct directions", distinctRing)

    # coplanarity: every ring vector lies in the plane spanned by omega
    var coplanar = true
    for h in hexagons:
      let om = omega(h)
      for k in 0..<6:
        let v = toFloat(allDirs[h.ring[k]])
        # om is antisymmetric, rank 2; P = -om*om projects onto the plane
        for a in 0..<nDim:
          var pv = 0.0
          for b in 0..<nDim:
            for c in 0..<nDim:
              pv -= om[a][b]*om[b][c]*v[c]
          if abs(pv - v[a]) > 1e-12: coplanar = false
    ok("all 6 vectors lie in the omega plane", coplanar)

    # 16 x 6 = 96 consecutive pairs = the 96 triangles through the site
    var pairs = initTable[(int, int), int]()
    var apexCount = 0
    for h in hexagons:
      var nApexHere = 0
      for k in 0..<6:
        var a = h.ring[k]
        var b = h.ring[(k+1) mod 6]
        if not isAxis(a) and not isAxis(b): inc nApexHere
        if a > b: swap(a, b)
        pairs[(a, b)] = pairs.getOrDefault((a, b)) + 1
      apexCount += nApexHere
      if nApexHere != 2: apexCount = -1000
    echo &"  distinct hexagon edge-pairs = {pairs.len}"
    ok(&"16 x 6 = 96 distinct triangles (got {pairs.len})", pairs.len == 96)
    var once = true
    for _, c in pairs:
      if c != 1: once = false
    ok("each of the 96 appears in exactly one hexagon", once)
    ok(&"16 x 2 = 32 apex triangles (got {apexCount})", apexCount == 32)

    # hexTriPaths: 6 closing loops, all the same orientation
    var hexCloses = true
    var hexConn = true
    var orient = true
    let centre = Site(cell: [1, 0, -2, 3], sub: 1)
    for h in hexagons:
      let ps = hexTriPaths(centre, h)
      var s0 = 0.0
      for k in 0..<6:
        let w = walkPath(ps[k], centre)
        if not w.closes: hexCloses = false
        if not w.connected: hexConn = false
        # signed area of the wedge w_k ^ w_{k+1} in the hexagon plane
        let
          a = toFloat(allDirs[h.ring[k]])
          b = toFloat(allDirs[h.ring[(k+1) mod 6]])
          om = omega(h)
        var s = 0.0
        for p in 0..<nDim:
          for q in 0..<nDim:
            s += om[p][q]*(a[p]*b[q] - a[q]*b[p])
        if k == 0: s0 = s
        if s0*s <= 0.0: orient = false
        if abs(abs(s) - abs(s0)) > 1e-12: orient = false
    ok("all 96 hexagon triangle loops close", hexCloses)
    ok("their links chain head-to-tail", hexConn)
    ok("all 6 triangles of a hexagon have the same orientation", orient)

    # also on the A sublattice
    var hexClosesA = true
    for h in hexagons:
      for p in hexTriPaths(Site(cell: [0, 0, 0, 0], sub: 0), h):
        if not walkPath(p, Site(cell: [0, 0, 0, 0], sub: 0)).closes:
          hexClosesA = false
    ok("same on the A sublattice", hexClosesA)

  test "8. clover reconstruction: fromOmega(projectOmega(F)) = F":
    var rng = initRand(20260821)
    var worst = 0.0
    for trial in 0..<20:
      var f: array[nDim, array[nDim, float]]
      for a in 0..<nDim:
        for b in (a+1)..<nDim:
          let v = 2.0*rng.rand(1.0) - 1.0
          f[a][b] = v
          f[b][a] = -v
      let g = fromOmega(projectOmega(f))
      for a in 0..<nDim:
        for b in 0..<nDim:
          worst = max(worst, abs(g[a][b] - f[a][b]))
    echo &"  max |fromOmega(projectOmega(F)) - F| = {worst:.3e}"
    ok(&"clover reconstruction exact to {worst:.2e} (validates the 3/8)",
       worst < 1e-12)
    # the 3/8 is the unique constant that works: perturb it and it must fail
    var scale = 0.0
    block:
      var f: array[nDim, array[nDim, float]]
      f[0][1] = 1.0
      f[1][0] = -1.0
      scale = fromOmega(projectOmega(f))[0][1]
    echo &"  reconstruction gain on F_01 = {scale:.15f}"
    ok("gain is exactly 1 (so 3/8 is the right prefactor)",
       abs(scale - 1.0) < 1e-14)

  test "9. point group orders 1152 / 384":
    let g16 = pointGroupOrder(honeycombDirs())
    let gc = pointGroupOrder(cubicDirs())
    echo &"  |point group| 16-cell = {g16}, cubic = {gc}"
    ok(&"16-cell point group = 1152 (got {g16})", g16 == 1152)
    ok(&"cubic point group = 384 (got {gc})", gc == 384)
