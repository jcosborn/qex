## Gauge fields on the 16-cell honeycomb (task **L**).
##
## 24 link fields per cell, indexed exactly as in doc/FORMULATION.md section 1.4:
##
##   uA[mu]    mu = 0..3    A(y) -> A(y+e_mu)
##   uB[mu]    mu = 0..3    B(y) -> B(y+e_mu)
##   uD[delta] delta = 0..15  B(y) -> A(y+delta)     (the 16 diagonals)
##
## The observable here is `triangleSum`, the average over the 64 triangles per
## cell of `Re Tr P / N`.  The triangle paths follow `hcgeom.apexTris` /
## `hcgeom.triPath`; this module never re-derives the geometry.
##
## Correctness first (PLAN.md section 1.2): the shifters are rebuilt on every
## call and nothing is cached.

import std/math
import base, layout, field, maths, rng
import physics/qcdTypes
import gauge
import hcgeom, hclayout

export hcgeom, hclayout

type
  HcGauge*[V: static[int], F] = object
    ## Two type parameters, not one: `V` is the SIMD length of the cell layout,
    ## which Nim cannot recover from the field type `F` in a field declaration.
    hl*: HcLayout[V]
    uA*, uB*: array[nDim, F]   ## A(y)->A(y+e_mu),  B(y)->B(y+e_mu)
    uD*: array[nDiag, F]       ## B(y)->A(y+delta)

template lo*(g: HcGauge): untyped = g.hl.lo

template eachLink*(g: HcGauge, u, body: untyped) =
  ## Run `body` once for each of the 24 link fields, bound to `u`, in the
  ## canonical order uA[0..3], uB[0..3], uD[0..15].  Allocation free, so it is
  ## safe inside a `threads:` block -- unlike `allLinks`.
  for hcMu in 0..<nDim:
    block:
      template u: untyped = g.uA[hcMu]
      body
  for hcMu in 0..<nDim:
    block:
      template u: untyped = g.uB[hcMu]
      body
  for hcDelta in 0..<nDiag:
    block:
      template u: untyped = g.uD[hcDelta]
      body

proc allLinks*[V: static[int], F](g: HcGauge[V, F]): seq[F] =
  ## All 24 link fields as a flat sequence: uA[0..3], uB[0..3], uD[0..15].
  ##
  ## **Allocates**, so it must be called *outside* a `threads:` block: QEX runs
  ## `threads:` bodies on raw OpenMP threads, which have no initialised Nim
  ## thread-local allocator, and heap traffic there segfaults at random.
  result = newSeq[F](nDirs)
  for mu in 0..<nDim:
    result[mu] = g.uA[mu]
    result[nDim+mu] = g.uB[mu]
  for d in 0..<nDiag:
    result[2*nDim+d] = g.uD[d]

proc link*[V: static[int], F](g: HcGauge[V, F], k: LinkKind, idx: int): F =
  ## uniform accessor matching `hcgeom.LinkRef`
  case k
  of lkA: g.uA[idx]
  of lkB: g.uB[idx]
  of lkD: g.uD[idx]

proc unit*(g: HcGauge) =
  ## all 24 link fields := 1
  eachLink(g, u):
    u := 1

proc newHcGauge*[V: static[int]](hl: HcLayout[V],
                                 nc: static[int] = getDefaultNc()): auto =
  ## 24 unit link fields on the cell layout of `hl`
  type F = type(hl.lo.ColorMatrix(nc))
  var g: HcGauge[V, F]
  g.hl = hl
  for mu in 0..<nDim:
    g.uA[mu] = hl.lo.ColorMatrix(nc)
    g.uB[mu] = hl.lo.ColorMatrix(nc)
  for d in 0..<nDiag:
    g.uD[d] = hl.lo.ColorMatrix(nc)
  g.unit
  g

proc newOneOf*[V: static[int], F](g: HcGauge[V, F]): HcGauge[V, F] =
  ## a second gauge object with the same shape (contents undefined)
  result.hl = g.hl
  for mu in 0..<nDim:
    result.uA[mu] = g.uA[mu].newOneOf
    result.uB[mu] = g.uB[mu].newOneOf
  for d in 0..<nDiag:
    result.uD[d] = g.uD[d].newOneOf

# The next five are allocation free and must be called inside a `threads:`
# block, as everywhere in QEX.

proc random*(g: HcGauge, r: var RNGField) =
  eachLink(g, u):
    when g.uA[0][0].nrows == 1:
      randomU(u, r)
    else:
      randomSU(u, r)

proc warm*(g: HcGauge, s: float, r: var RNGField) =
  eachLink(g, u):
    when g.uA[0][0].nrows == 1:
      u.gaussian r
      u := (1-s) + s*u
      u.projectU
    else:
      u.warmSU(s, r)

proc reunit*(g: HcGauge) =
  ## projectSU (projectU for Nc = 1) on all 24 link fields
  eachLink(g, u):
    when g.uA[0][0].nrows == 1:
      u.projectU
    else:
      u.projectSU

proc checkSU*(g: HcGauge): tuple[avg, max: float] =
  ## same measure as QEX's `checkSU`, averaged over all 24 link fields
  var a, b: float
  eachLink(g, u):
    for s in u:
      let dev = u[s].checkSU
      a += dev.simdSum
      let m = dev.simdMax
      if b < m: b = m
  threadRankSum a
  threadRankMax b
  const nc = g.uA[0][0].nrows
  let vol = g.uA[0].l.physVol
  let c = float(2*(nc*nc+1))
  (sqrt(a/(c*float(nDirs*vol))), sqrt(b/c))

proc `:=`*(a: HcGauge, b: HcGauge) =
  for mu in 0..<nDim:
    a.uA[mu] := b.uA[mu]
    a.uB[mu] := b.uB[mu]
  for d in 0..<nDiag:
    a.uD[d] := b.uD[d]

# ---------------------------------------------------------------------------
# triangle loops
# ---------------------------------------------------------------------------

proc triangleTrace*[V: static[int], F](g: HcGauge[V, F]): tuple[re, im: float] =
  ## `sum_x sum_{i=1..32} Tr P_i(x)` over all `2 N_cells` sites, i.e. over all
  ## `64 N_cells` triangles, each counted once at its apex.
  ##
  ## apex on B, cell y (FORMULATION 2.2):
  ##   `uD[d](y) . uA[mu](y+d) . uD[d'](y)^dag`,  d' = d or 2^mu
  ## apex on A, cell y:
  ##   `uD[db](y-db)^dag . uB[mu](y-db) . uD[db'](y-db')`,  db = d xor 15,
  ##   db' = db - 2^mu.
  ## The A loop is evaluated re-based at the corner `z = y-db` (the trace is
  ## cyclic and `y -> z` is a bijection of the torus), which turns three
  ## diagonal shifts into one single-axis shift:
  ##   `uB[mu](z) . uD[db'](z+e_mu) . uD[db](z)^dag`.
  # 16-way forward shifts of the A-sublattice axis links: uA[mu](y+delta)
  type SH = type(newHcShift16(g.uA[0], 1))
  var shA: array[nDim, SH]
  for mu in 0..<nDim:
    shA[mu] = newHcShift16(g.uA[mu], 1)
  # one reusable single-axis +e_mu shifter per direction, for the A apexes
  type SS = type(newShifter(g.uD[0], 0, 1))
  var sD: array[nDim, SS]
  for mu in 0..<nDim:
    sD[mu] = newShifter(g.uD[0], mu, 1)
  var m = g.uA[0].newOneOf
  var tr: type(trace(m))
  threads:
    m := 0
    for mu in 0..<nDim:
      shA[mu].run
    threadBarrier()
    for t in apexTris:                       # 32 apexes on B
      m += (g.uD[t.delta] * shA[t.mu].f[t.delta]) * g.uD[t.deltaP].adj
    for t in apexTris:                       # 32 apexes on A
      let
        db = t.delta xor 15
        dbp = t.deltaP xor 15                # = db - 2^mu
        s = sD[t.mu] ^* g.uD[dbp]
      m += (g.uB[t.mu] * s) * g.uD[db].adj
    threadBarrier()
    tr = trace(m)
  (tr.re, tr.im)

proc triangleSum*(g: HcGauge): float =
  ## `(1/(32 N_sites)) sum_x sum_{i=1..32} Re Tr P_i(x) / N`, with
  ## `N_sites = 2 N_cells`.  Unit gauge gives exactly 1.
  ## Only the real part enters the action (FORMULATION section 3).
  const nc = g.uA[0][0].nrows
  let tr = triangleTrace(g)
  tr.re/float(nTriPerSite*g.hl.nSites*nc)

proc triangleSumIm*(g: HcGauge): float =
  ## The imaginary companion of `triangleSum`.  It is exactly zero for unit and
  ## pure-gauge configurations, but **not** identically zero in general: the 32
  ## apex triangles are enumerated with a fixed orientation, so unlike the cubic
  ## plaquette sum there is no reversed partner to cancel against.  For a random
  ## configuration it is the same statistical size as the real part.
  ## Recomputes everything; use `triangleTrace` if you want both parts.
  const nc = g.uA[0][0].nrows
  let tr = triangleTrace(g)
  tr.im/float(nTriPerSite*g.hl.nSites*nc)

# ---------------------------------------------------------------------------
# gauge transformations
# ---------------------------------------------------------------------------

proc gaugeTransform*[V: static[int], F](g: HcGauge[V, F], vA, vB: F) =
  ## `u -> V(start) u V(end)^dag` on all 24 link fields:
  ##   uA[mu](y)    -> vA(y) uA[mu](y) vA(y+e_mu)^dag
  ##   uB[mu](y)    -> vB(y) uB[mu](y) vB(y+e_mu)^dag
  ##   uD[delta](y) -> vB(y) uD[delta](y) vA(y+delta)^dag
  ## `vA` lives on the A sites, `vB` on the B sites; both are cell fields.
  var shA = newHcShift16(vA, 1)              # vA(y+delta), 15 shifts
  type SS = type(newShifter(vB, 0, 1))
  var sB: array[nDim, SS]
  for mu in 0..<nDim:
    sB[mu] = newShifter(vB, mu, 1)
  var t = vA.newOneOf
  threads:
    shA.run
    for mu in 0..<nDim:
      discard sB[mu] ^* vB
    threadBarrier()
    for mu in 0..<nDim:
      # vA(y+e_mu) is just the 2^mu leaf of the shift tree
      t := vA * g.uA[mu]
      g.uA[mu] := t * shA.f[1 shl mu].adj
      t := vB * g.uB[mu]
      g.uB[mu] := t * sB[mu].field.adj
    for d in 0..<nDiag:
      t := vB * g.uD[d]
      g.uD[d] := t * shA.f[d].adj
