## Geometry of the 4D 16-cell honeycomb (D4* lattice).
##
## This module is pure combinatorics: no QEX `Layout`, no fields, no threads.
## It is the normative source of the direction / triangle / hexagon tables that
## `hcgauge`, `hcaction`, `hctopo` and `hcwilson` all build on.
##
## See doc/FORMULATION.md for the derivations; every number produced here is
## checked in tests/tgeom.nim.
##
## Representation
## --------------
## Sites are `Z^4 union (Z+1/2)^4`, encoded as a hypercubic lattice of *cells*
## with a two-site basis:
##   sublattice A (sub=0) at cell y,  sublattice B (sub=1) at y + (1/2,1/2,1/2,1/2).
##
## The 24 nearest-neighbour directions all have length 1:
##   8 axis     +- e_mu
##  16 diagonal  d(delta)_mu = delta_mu - 1/2,   delta in {0,1}^4
## To keep everything exact they are stored doubled, `Vec2 = 2*n`, so axis
## directions are `(+-2,0,0,0)` and diagonals are `(+-1,+-1,+-1,+-1)`.

import std/math

const
  nDim* = 4
  nDirs* = 24          ## nearest neighbours (the 4D kissing number)
  nAxis* = 8
  nDiag* = 16
  nTriPerSite* = 32    ## triangles with their apex at a given site
  nTriThruSite* = 96   ## triangles containing a given site
  nHexPerSite* = 16    ## hexagons (coplanar rings of 6 triangles)
  nTriPerLink* = 8
  nLinkPerSite* = 12
  nSubs* = 2

type
  Vec2* = array[nDim, int]
    ## a 4-vector stored doubled: the true vector is `Vec2 / 2`.

  Cell* = array[nDim, int]
    ## integer cell offset

  Site* = object
    ## a lattice site, as a cell offset plus sublattice index
    cell*: Cell
    sub*: int          ## 0 = A, 1 = B

  LinkKind* = enum
    lkA,               ## uA[mu]: A(y) -> A(y+e_mu)
    lkB,               ## uB[mu]: B(y) -> B(y+e_mu)
    lkD                ## uD[delta]: B(y) -> A(y+delta)

  LinkRef* = object
    ## one traversal of one link, relative to a base cell
    kind*: LinkKind
    idx*: int          ## mu in 0..3 for lkA/lkB, delta in 0..15 for lkD
    cell*: Cell        ## cell offset at which the link field lives
    dag*: bool         ## true = traverse the link backwards (use the adjoint)

  Path* = seq[LinkRef]
    ## a sequence of link traversals, in order; the leftmost is the first step

  ApexTri* = object
    ## a triangle labelled by its unique apex corner.
    ## Its two apex edges are the diagonals d(delta) and d(deltaP);
    ## the third (axis) edge runs along +e_mu between their far ends.
    delta*, deltaP*: int   ## deltaP = delta or (1 shl mu); delta has bit mu clear
    mu*: int

  Hexagon* = object
    ## six coplanar unit vectors at 60 degrees, hence six triangles sharing the
    ## centre site.  Labelled by the axis direction mu it contains and by the
    ## diagonal delta (bit mu clear) of one of its two apex triangles.
    mu*: int
    delta*: int
    ring*: array[6, int]   ## the 6 direction indices, in cyclic (60 degree) order

# ---------------------------------------------------------------------------
# direction tables
# ---------------------------------------------------------------------------

func axisIndex*(mu: int, neg: bool): int {.inline.} =
  ## global direction index of +-e_mu; 0..7
  2*mu + (if neg: 1 else: 0)

func diagIndex*(delta: int): int {.inline.} =
  ## global direction index of d(delta); 8..23
  nAxis + delta

func isAxis*(dir: int): bool {.inline.} = dir < nAxis
func dirDelta*(dir: int): int {.inline.} = dir - nAxis
func dirMu*(dir: int): int {.inline.} = dir div 2
func dirNeg*(dir: int): bool {.inline.} = (dir and 1) == 1

func opposite*(dir: int): int {.inline.} =
  ## the direction index of -n
  if dir < nAxis: dir xor 1
  else: nAxis + ((dir - nAxis) xor 15)

func deltaVec*(delta: int): Vec2 =
  ## 2*d(delta); components are +-1
  for mu in 0..<nDim:
    result[mu] = 2*((delta shr mu) and 1) - 1

func dirVec*(dir: int): Vec2 =
  ## 2*n for global direction index `dir`
  if dir < nAxis:
    let mu = dir div 2
    result[mu] = if (dir and 1) == 1: -2 else: 2
  else:
    result = deltaVec(dir - nAxis)

const
  allDirs* = block:
    var a: array[nDirs, Vec2]
    for d in 0..<nDirs: a[d] = dirVec(d)
    a

func dot2*(a, b: Vec2): int {.inline.} =
  ## 4 * (a/2).(b/2) -- i.e. the true dot product times 4
  for mu in 0..<nDim: result += a[mu]*b[mu]

func add2*(a, b: Vec2): Vec2 {.inline.} =
  for mu in 0..<nDim: result[mu] = a[mu] + b[mu]

func sub2*(a, b: Vec2): Vec2 {.inline.} =
  for mu in 0..<nDim: result[mu] = a[mu] - b[mu]

func neg2*(a: Vec2): Vec2 {.inline.} =
  for mu in 0..<nDim: result[mu] = -a[mu]

func dirIndexOf*(v: Vec2): int =
  ## index of `v` among the 24 directions, or -1
  for d in 0..<nDirs:
    if allDirs[d] == v: return d
  -1

func toFloat*(v: Vec2): array[nDim, float] =
  for mu in 0..<nDim: result[mu] = 0.5*v[mu].float

# ---------------------------------------------------------------------------
# walking the lattice
# ---------------------------------------------------------------------------

func shiftCell(c: Cell, mu, s: int): Cell {.inline.} =
  result = c
  result[mu] += s

func subCell(c: Cell, delta: int): Cell {.inline.} =
  ## c - delta, delta read as a 0/1 bit vector
  result = c
  for mu in 0..<nDim: result[mu] -= (delta shr mu) and 1

func addCell(c: Cell, delta: int): Cell {.inline.} =
  result = c
  for mu in 0..<nDim: result[mu] += (delta shr mu) and 1

proc step*(s: Site, dir: int): tuple[link: LinkRef, dest: Site] =
  ## Move from site `s` one step along direction `dir`, returning the link that
  ## must be traversed and the destination site.
  if dir < nAxis:
    let
      mu = dir div 2
      neg = (dir and 1) == 1
      kind = if s.sub == 0: lkA else: lkB
    if neg:
      let c = shiftCell(s.cell, mu, -1)
      result.link = LinkRef(kind: kind, idx: mu, cell: c, dag: true)
      result.dest = Site(cell: c, sub: s.sub)
    else:
      result.link = LinkRef(kind: kind, idx: mu, cell: s.cell, dag: false)
      result.dest = Site(cell: shiftCell(s.cell, mu, 1), sub: s.sub)
  else:
    let delta = dir - nAxis
    if s.sub == 1:
      # B(y) -> A(y+delta) is exactly uD[delta](y)
      result.link = LinkRef(kind: lkD, idx: delta, cell: s.cell, dag: false)
      result.dest = Site(cell: addCell(s.cell, delta), sub: 0)
    else:
      # A(y) + d(delta) = B(y - dbar), reached by reversing uD[dbar](y-dbar)
      let
        dbar = delta xor 15
        c = subCell(s.cell, dbar)
      result.link = LinkRef(kind: lkD, idx: dbar, cell: c, dag: true)
      result.dest = Site(cell: c, sub: 1)

proc walk*(start: Site, dirs: openArray[int]): tuple[path: Path, dest: Site] =
  ## follow a sequence of directions, collecting the links traversed
  var s = start
  for d in dirs:
    let (l, n) = step(s, d)
    result.path.add l
    s = n
  result.dest = s

proc loop*(start: Site, dirs: openArray[int]): Path =
  ## like `walk`, but assert that the path closes
  let (p, dest) = walk(start, dirs)
  doAssert dest == start, "loop does not close"
  p

# ---------------------------------------------------------------------------
# triangles
# ---------------------------------------------------------------------------

const
  apexTris* = block:
    ## The 32 triangles whose apex is the site in question, labelled (delta, mu)
    ## with bit mu of delta clear.  Ordering: mu slowest, then the 3 free bits.
    var a: array[nTriPerSite, ApexTri]
    var n = 0
    for mu in 0..<nDim:
      for k in 0..<8:
        # insert a 0 bit at position mu into the 3-bit number k
        let lo = k and ((1 shl mu) - 1)
        let hi = (k shr mu) shl (mu + 1)
        let delta = lo or hi
        a[n] = ApexTri(delta: delta, deltaP: delta or (1 shl mu), mu: mu)
        inc n
    doAssert n == nTriPerSite
    a

proc triDirs*(t: ApexTri): array[3, int] =
  ## the three direction indices traversed, starting from the apex
  [diagIndex(t.delta), axisIndex(t.mu, false), opposite(diagIndex(t.deltaP))]

proc triPath*(apex: Site, t: ApexTri): Path =
  ## The closed 3-link loop of triangle `t` based at its apex.
  ##
  ## apex = A(y):  uD[db](y-db)^dag . uB[mu](y-db) . uD[db'](y-db')
  ## apex = B(y):  uD[delta](y)     . uA[mu](y+delta) . uD[delta'](y)^dag
  ## with db = delta xor 15, db' = deltaP xor 15.
  loop(apex, triDirs(t))

# ---------------------------------------------------------------------------
# hexagons
# ---------------------------------------------------------------------------

const
  hexagons* = block:
    ## 16 hexagons per site.  Each apex triangle (delta,mu) sits in one hexagon;
    ## the map is 2:1, the partner being (delta bitcomplemented in the 3 free
    ## bits, same mu).  Canonical representative: the smaller delta.
    var a: array[nHexPerSite, Hexagon]
    var n = 0
    for mu in 0..<nDim:
      for k in 0..<8:
        let lo = k and ((1 shl mu) - 1)
        let hi = (k shr mu) shl (mu + 1)
        let delta = lo or hi                 # bit mu clear
        # partner apex triangle of the same hexagon
        let partner = (delta xor 15) xor (1 shl mu)
        if delta > partner: continue         # keep one representative
        let dm = diagIndex(delta)            # d^- , mu-component -1/2
        let dp = diagIndex(delta or (1 shl mu))   # d^+ , mu-component +1/2
        let em = axisIndex(mu, false)
        # cyclic order, 60 degrees apart:  d^-, d^+, e_mu, -d^-, -d^+, -e_mu
        a[n] = Hexagon(mu: mu, delta: delta,
                       ring: [dm, dp, em, opposite(dm), opposite(dp),
                              opposite(em)])
        inc n
    doAssert n == nHexPerSite
    a

proc hexTriPaths*(centre: Site, h: Hexagon): array[6, Path] =
  ## The six triangles of hexagon `h`, each as a closed loop based at `centre`
  ## and all traversed with the same (positive) orientation in the plane.
  ## Triangle k has vertices  centre, centre+w_k, centre+w_{k+1}.
  for k in 0..<6:
    let
      wk = h.ring[k]
      wk1 = h.ring[(k+1) mod 6]
      # edge from centre+w_k to centre+w_{k+1} is w_{k+2}
      wk2 = h.ring[(k+2) mod 6]
    result[k] = loop(centre, [wk, wk2, opposite(wk1)])

proc omega*(h: Hexagon): array[nDim, array[nDim, float]] =
  ## Unit area 2-form of the plane of hexagon `h`:
  ##   e = e_mu,  f = (1/sqrt 3) sum_{nu != mu} sigma_nu e_nu,
  ##   Omega_ab = e_a f_b - e_b f_a.
  ## sigma_nu is read off the diagonal d(delta): 2 d_nu = sigma_nu.
  let
    mu = h.mu
    dv = deltaVec(h.delta)
    s = 1.0/sqrt(3.0)
  for nu in 0..<nDim:
    if nu == mu: continue
    let f = s*dv[nu].float
    result[mu][nu] = f
    result[nu][mu] = -f

proc fromOmega*(fOmega: array[nHexPerSite, float]): array[nDim, array[nDim, float]] =
  ## Reconstruct F_munu from the 16 in-plane field strengths, FORMULATION (4.2):
  ##   F_munu = (3/8) sum_h Omega^(h)_munu F_Omega^(h)
  for h in 0..<nHexPerSite:
    let om = omega(hexagons[h])
    for a in 0..<nDim:
      for b in 0..<nDim:
        result[a][b] += 0.375*om[a][b]*fOmega[h]

proc projectOmega*(f: array[nDim, array[nDim, float]]): array[nHexPerSite, float] =
  ## F_Omega^(h) = (1/2) Omega^(h)_ab F_ab -- the inverse companion of `fromOmega`
  for h in 0..<nHexPerSite:
    let om = omega(hexagons[h])
    var s = 0.0
    for a in 0..<nDim:
      for b in 0..<nDim:
        s += om[a][b]*f[a][b]
    result[h] = 0.5*s

# ---------------------------------------------------------------------------
# point symmetry group
# ---------------------------------------------------------------------------

proc pointGroupOrder*(vecs: openArray[Vec2]): int =
  ## Number of 4x4 orthogonal matrices mapping the neighbour set onto itself.
  ## Expected: 1152 for the 24 honeycomb directions (Weyl group of F4),
  ##            384 for the 8 axis directions (Weyl group of B4).
  ##
  ## A candidate is fixed by the images of e_0..e_3, each of which must itself be
  ## a neighbour vector.  Work with M = 2A so that everything stays integral:
  ## A orthogonal  <=>  M^T M = 4 I;  and (M v2)/2 is the image of v2.
  let n = vecs.len
  var idx = newSeq[Vec2](n)
  for i in 0..<n: idx[i] = vecs[i]
  proc inSet(v: Vec2): bool =
    for w in idx:
      if w == v: return true
    false
  var m: array[nDim, Vec2]      # m[mu] = column mu of M = image of e_mu, doubled
  var cnt = 0
  for i0 in 0..<n:
    m[0] = idx[i0]
    for i1 in 0..<n:
      m[1] = idx[i1]
      if dot2(m[0], m[1]) != 0: continue
      for i2 in 0..<n:
        m[2] = idx[i2]
        if dot2(m[0], m[2]) != 0 or dot2(m[1], m[2]) != 0: continue
        for i3 in 0..<n:
          m[3] = idx[i3]
          if dot2(m[0], m[3]) != 0 or dot2(m[1], m[3]) != 0 or
             dot2(m[2], m[3]) != 0: continue
          # columns are orthogonal and each has (doubled) norm^2 = 4, so M^T M = 4I
          var ok = true
          for v in idx:
            var w: Vec2
            for a in 0..<nDim:
              var t = 0
              for b in 0..<nDim: t += m[b][a]*v[b]
              if (t and 1) != 0: ok = false; break
              w[a] = t div 2
            if not ok: break
            if not inSet(w): ok = false; break
          if ok: inc cnt
  cnt

# ---------------------------------------------------------------------------
# convenience: the direction sets
# ---------------------------------------------------------------------------

proc honeycombDirs*(): seq[Vec2] =
  for d in 0..<nDirs: result.add allDirs[d]

proc cubicDirs*(): seq[Vec2] =
  for d in 0..<nAxis: result.add allDirs[d]

when isMainModule:
  echo "16-cell honeycomb geometry"
  echo "  directions      : ", nDirs
  echo "  apex triangles  : ", apexTris.len
  echo "  hexagons        : ", hexagons.len
  echo "  point group     : ", pointGroupOrder(honeycombDirs())
  echo "  cubic pt group  : ", pointGroupOrder(cubicDirs())
