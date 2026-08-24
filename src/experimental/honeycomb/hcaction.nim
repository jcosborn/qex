## Triangle gauge action, staple derivative and MD force on the 16-cell
## honeycomb (task **A**).
##
## Action (doc/FORMULATION.md section 3):
##
##   S(U) = (beta/2) sum_x sum_{i=1}^{32} ( 1 - (1/N) Re Tr P_i(x) )
##
## with `x` over all `2 N_cells` sites and `i` over the 32 apex triangles at
## `x`, i.e. over all `64 N_cells` triangles, each counted exactly once
## (`hcgeom.apexTris`).
##
## Staple recipe -- derived from `hcgeom.triPath`, not from prose
## --------------------------------------------------------------
## Every link field sits in exactly 8 triangles.  For each occurrence write the
## triangle loop re-based (cyclically, reversing orientation where the link
## enters daggered -- Re Tr is orientation blind) so it starts with the link
## un-daggered, anchored at the cell where the link field lives:
##
##   P = U_l(x) . W(x)          W = product of the other two edges.
##
## QEX's staple convention keeps V = W^dag so that P = U_l V^dag, exactly like
## `makeStaples`/`gaugeActionDeriv`/`contractProjectTAH` on the cubic lattice.
## At module load this file *scans all 64 triangles per cell* via
## `hcgeom.triPath`, collects all 24x8 = 192 (link, triangle) incidences, and
## `doAssert`s them equal to the closed-form six families the executor below
## implements.  Action and force therefore cannot drift apart from the
## geometry or from each other.
##
## The six families, for each `apexTris` entry `(mu, delta)` [bit mu of delta
## clear, `deltaP = delta or 2^mu`, `db = delta xor 15` (bit mu set),
## `dbp = db xor 2^mu = deltaP xor 15`]:
##
##   apex-B triangle  uD[delta](z) . uA[mu](z+delta) . uD[deltaP](z)^dag :
##     V(uD[delta])(z)  = uD[deltaP](z) . uA[mu](z+delta)^dag
##     V(uD[deltaP])(z) = uD[delta](z)  . uA[mu](z+delta)
##     V(uA[mu])(y)     = uD[delta](y-delta)^dag . uD[deltaP](y-delta)
##   apex-A triangle (re-based at z = apex-db, see hcgauge.triangleTrace):
##                    uD[db](z)^dag . uB[mu](z) . uD[dbp](z+e_mu)  reversed ->
##     V(uB[mu])(z)     = uD[db](z)  . uD[dbp](z+e_mu)^dag
##     V(uD[db])(z)     = uB[mu](z)  . uD[dbp](z+e_mu)
##     V(uD[dbp])(y)    = uB[mu](y-e_mu)^dag . uD[db](y-e_mu)
##
## Force convention (pinned numerically in tests/taction.nim)
## ----------------------------------------------------------
##   hcActionDeriv:  D_l = (beta/(2N)) sum_{k=1}^{8} V_k     (NOT TAH projected)
##   hcForce:        f_l = projectTAH( U_l . D_l^dag )
##                       = (beta/(2N)) projectTAH( sum_{k=1}^{8} U_l W_k )
##
## and f_l satisfies, for arbitrary traceless anti-Hermitian momenta P_l (one
## per link field),
##
##   d/ds S( exp(s P_l) U_l ) |_{s=0}  =  sum_l redot( P_l, f_l )
##
## with redot(a,b) = Re tr(a^dag b) summed over all sites (QEX field `redot`).
## This is the same convention QEX's `gaugeForce` satisfies (both are verified
## against numerical derivatives in tests/taction.nim), so the standard
## leapfrog of examples/puregaugehmc.nim,
##
##   p_l -= dt * f_l;        U_l := exp(dt * p_l) * U_l
##
## conserves H = S + (1/2) sum_l redot(p_l, p_l).
##
## Performance: `HcActionWork` holds all shifters and scratch fields; create it
## once (outside `threads:`) and reuse it every MD step.  Nothing below
## allocates inside a `threads:` block.  The convenience overloads without a
## work object allocate a fresh one per call -- fine for tests only.

import std/algorithm
import base, layout, field, maths, rng
import physics/qcdTypes
import gauge
import hcgeom, hclayout, hcgauge

export hcgauge

# ---------------------------------------------------------------------------
# staple recipe: scan hcgeom.triPath and check the executor's six families
# ---------------------------------------------------------------------------

type
  HcStapleFactor = object
    kind: LinkKind
    idx: int
    cell: Cell     ## offset relative to the cell of the target link
    dag: bool

func relCell(a, b: Cell): Cell =
  for i in 0..<nDim: result[i] = a[i] - b[i]

func negCell(a: Cell): Cell =
  for i in 0..<nDim: result[i] = -a[i]

func deltaCell(delta: int): Cell =
  for i in 0..<nDim: result[i] = (delta shr i) and 1

func muCell(mu: int): Cell =
  result[mu] = 1

proc termKey(tk: LinkKind, ti: int, a, b: HcStapleFactor): string =
  ## canonical string form of one staple term:
  ##   V(x) = a(x + a.cell)^{a.dag} . b(x + b.cell)^{b.dag}  for link (tk, ti)
  $tk & " " & $ti & " | " & $a.kind & " " & $a.idx & " " & $a.cell & " " &
    $a.dag & " | " & $b.kind & " " & $b.idx & " " & $b.cell & " " & $b.dag

proc linkSlot(k: LinkKind, idx: int): int =
  case k
  of lkA: idx
  of lkB: nDim + idx
  of lkD: 2*nDim + idx

proc scanStapleTerms(): seq[string] =
  ## Every (link, triangle) incidence obtained by walking all 64 triangles of a
  ## cell with `hcgeom.triPath`, re-basing each loop at each of its three links
  ## and converting to the staple V = W^dag.  Independent of the executor.
  var count: array[nDirs, int]
  for sub in 0..1:
    for t in apexTris:
      let path = triPath(Site(cell: [0, 0, 0, 0], sub: sub), t)
      for pos in 0..2:
        var lp: array[3, LinkRef]
        if not path[pos].dag:
          for k in 0..<3: lp[k] = path[(pos+k) mod 3]
        else:
          # reverse the loop so our link enters un-daggered, then rotate it
          # to the front
          var rev: array[3, LinkRef]
          for k in 0..<3:
            rev[k] = path[2-k]
            rev[k].dag = not rev[k].dag
          let rpos = 2 - pos
          for k in 0..<3: lp[k] = rev[(rpos+k) mod 3]
        doAssert (not lp[0].dag) and lp[0].kind == path[pos].kind and
                 lp[0].idx == path[pos].idx and lp[0].cell == path[pos].cell
        # W = lp[1] . lp[2]  =>  V = W^dag = lp[2]^flip . lp[1]^flip
        let a = HcStapleFactor(kind: lp[2].kind, idx: lp[2].idx,
                               cell: relCell(lp[2].cell, lp[0].cell),
                               dag: not lp[2].dag)
        let b = HcStapleFactor(kind: lp[1].kind, idx: lp[1].idx,
                               cell: relCell(lp[1].cell, lp[0].cell),
                               dag: not lp[1].dag)
        result.add termKey(lp[0].kind, lp[0].idx, a, b)
        inc count[linkSlot(lp[0].kind, lp[0].idx)]
  for s in 0..<nDirs:
    doAssert count[s] == nTriPerLink,
      "hcaction: link slot " & $s & " sits in " & $count[s] & " triangles"

proc executorStapleTerms(): seq[string] =
  ## the six families hcActionDeriv implements, in the same canonical form
  template add(tk, ti, ak, ai, ac, ad, bk, bi, bc, bd: untyped) =
    result.add termKey(tk, ti,
      HcStapleFactor(kind: ak, idx: ai, cell: ac, dag: ad),
      HcStapleFactor(kind: bk, idx: bi, cell: bc, dag: bd))
  let zero: Cell = [0, 0, 0, 0]
  for t in apexTris:
    let
      mu = t.mu
      delta = t.delta
      deltaP = t.deltaP
      db = delta xor 15
      dbp = deltaP xor 15
      cD = deltaCell(delta)
      cDm = negCell(cD)
      cMu = muCell(mu)
      cMum = negCell(cMu)
    # apex-B triangle (mu, delta)
    add(lkD, delta,  lkD, deltaP, zero, false,  lkA, mu, cD, true)
    add(lkD, deltaP, lkD, delta,  zero, false,  lkA, mu, cD, false)
    add(lkA, mu,     lkD, delta,  cDm,  true,   lkD, deltaP, cDm, false)
    # apex-A triangle (mu, db)
    add(lkB, mu,     lkD, db,     zero, false,  lkD, dbp, cMu, true)
    add(lkD, db,     lkB, mu,     zero, false,  lkD, dbp, cMu, false)
    add(lkD, dbp,    lkB, mu,     cMum, true,   lkD, db,  cMum, false)

block:
  ## run the recipe check once at module load; pure combinatorics, no QEX state
  var a = scanStapleTerms()
  var b = executorStapleTerms()
  doAssert a.len == 3*2*nTriPerSite and b.len == a.len
  sort a
  sort b
  doAssert a == b,
    "hcaction: staple recipe disagrees with hcgeom.triPath -- refusing to run"

# ---------------------------------------------------------------------------
# persistent work space
# ---------------------------------------------------------------------------

type
  HcActionWork*[V: static[int], F, SH, SS] = ref object
    ## All shifters and scratch fields for `hcAction`/`hcActionDeriv`/
    ## `hcForce`.  Create once with `newHcActionWork` (allocates), then every
    ## call is allocation free.  The object is rebound to the gauge field
    ## passed to each call, so one work object can serve several `HcGauge`s of
    ## the same shape.
    shA*: array[nDim, SH]           ## forward 16-trees: uA[mu](x+delta)
    sF*: array[nDim, SS]            ## +e_mu shifters: uD[dbp](x+e_mu)
    sB*: array[3, array[nDim, SS]]  ## -e_mu shifter chains (3 levels)
    t*: F                           ## 2-link product scratch

proc newHcActionWork*[V: static[int], F](g: HcGauge[V, F]): auto =
  ## Allocate shifters and scratch for gauge fields shaped like `g`.
  ## Call outside `threads:`.
  type SH = type(newHcShift16(g.uA[0], 1))
  type SS = type(newShifter(g.uD[0], 0, 1))
  var w = HcActionWork[V, F, SH, SS]()
  for mu in 0..<nDim:
    w.shA[mu] = newHcShift16(g.uA[mu], 1)
    w.sF[mu] = newShifter(g.uD[0], mu, 1)
    for l in 0..<3:
      w.sB[l][mu] = newShifter(g.uD[0], mu, -1)
  w.t = g.uA[0].newOneOf
  w

proc rebind[V: static[int], F, SH, SS](w: HcActionWork[V, F, SH, SS],
                                       g: HcGauge[V, F]) =
  ## point the shift trees at g's axis links (outside `threads:` -- these are
  ## ref assignments)
  for mu in 0..<nDim:
    w.shA[mu].setSrc g.uA[mu]

# ---------------------------------------------------------------------------
# action
# ---------------------------------------------------------------------------

proc hcAction*[V: static[int], F, SH, SS](w: HcActionWork[V, F, SH, SS],
                                          beta: float, g: HcGauge[V, F]): float =
  ## S = (beta/2) sum_x sum_{i=1}^{32} (1 - Re Tr P_i(x)/N).
  ##
  ## Computed from the same staple recipe as the force (each triangle once,
  ## through redot(staple, link) at its uD[delta] / uB[mu] edge).  This is an
  ## independent code path from `hcgauge.triangleTrace` (which accumulates the
  ## loop matrices and traces once); the two are cross-checked in
  ## tests/taction.nim.
  const nc = g.uA[0][0].nrows
  rebind(w, g)
  var tsum = 0.0
  threads:
    for mu in 0..<nDim:
      w.shA[mu].run
    var a = 0.0
    for t in apexTris:
      let
        mu = t.mu
        delta = t.delta
        deltaP = t.deltaP
        db = delta xor 15
        dbp = deltaP xor 15
      # apex-B triangle: Re Tr[ uD[delta](z) uA[mu](z+delta) uD[deltaP](z)^dag ]
      #                = redot( V(uD[delta]), uD[delta] )
      w.t := g.uD[deltaP] * w.shA[mu].f[delta].adj
      a += redot(w.t, g.uD[delta])
      # apex-A triangle: Re Tr[ uB[mu](z) uD[dbp](z+e_mu) uD[db](z)^dag ]
      #                = redot( V(uB[mu]), uB[mu] )
      let s2 = w.sF[mu] ^* g.uD[dbp]
      w.t := g.uD[db] * s2.adj
      a += redot(w.t, g.uB[mu])
    threadMaster: tsum = a
  let nTri = float(nTriPerSite*g.hl.nSites)   # 64 per cell
  0.5*beta*(nTri - tsum/float(nc))

proc hcAction*[V: static[int], F](beta: float, g: HcGauge[V, F]): float =
  ## Convenience overload: allocates a fresh work object (tests only; HMC
  ## should create one `HcActionWork` and reuse it).
  var w = newHcActionWork(g)
  hcAction(w, beta, g)

# ---------------------------------------------------------------------------
# derivative (staple sums) and force
# ---------------------------------------------------------------------------

proc hcActionDeriv*[V: static[int], F, SH, SS](w: HcActionWork[V, F, SH, SS],
                                               beta: float, g: HcGauge[V, F],
                                               f: var HcGauge[V, F]) =
  ## f_l := (beta/(2N)) sum_{k=1}^{8} V_k(l), the staple sums in QEX's
  ## convention (P = U_l V^dag), NOT yet TAH projected.  `f` is overwritten.
  ## The 24 fields of `f` must be distinct from those of `g`.
  const nc = g.uA[0][0].nrows
  let ff = f              # value copy: the 24 fields are refs
  rebind(w, g)
  let cf = beta/(2.0*float(nc))
  threads:
    eachLink(ff, u):
      u := 0
    for mu in 0..<nDim:
      w.shA[mu].run
    for t in apexTris:
      let
        mu = t.mu
        delta = t.delta
        deltaP = t.deltaP
        db = delta xor 15
        dbp = deltaP xor 15
      # ---- apex-B triangle  uD[delta](z) uA[mu](z+delta) uD[deltaP](z)^dag
      ff.uD[delta]  += g.uD[deltaP] * w.shA[mu].f[delta].adj
      ff.uD[deltaP] += g.uD[delta] * w.shA[mu].f[delta]
      # staple of uA[mu] lives shifted: V(y) = [uD[delta]^dag uD[deltaP]](y-delta)
      w.t := g.uD[delta].adj * g.uD[deltaP]
      block:
        var cur = w.t
        var lev = 0
        for b in 0..<nDim:
          if ((delta shr b) and 1) != 0:
            cur = w.sB[lev][b] ^* cur
            inc lev
        ff.uA[mu] += cur
      # ---- apex-A triangle  uD[db](z)^dag uB[mu](z) uD[dbp](z+e_mu)
      let s2 = w.sF[mu] ^* g.uD[dbp]
      ff.uB[mu] += g.uD[db] * s2.adj
      ff.uD[db] += g.uB[mu] * s2
      # staple of uD[dbp]: V(y) = [uB[mu]^dag uD[db]](y-e_mu)
      w.t := g.uB[mu].adj * g.uD[db]
      let rs = w.sB[0][mu] ^* w.t
      ff.uD[dbp] += rs
    eachLink(ff, u):
      u := cf*u

proc hcActionDeriv*[V: static[int], F](beta: float, g: HcGauge[V, F],
                                       f: var HcGauge[V, F]) =
  ## Convenience overload: allocates a fresh work object per call.
  var w = newHcActionWork(g)
  hcActionDeriv(w, beta, g, f)

template contractPTAH(u, ff: untyped) =
  ## f := projectTAH( u f^dag ) elementwise; the body of QEX's
  ## `contractProjectTAH` (src/gauge/gaugeUtils.nim:406), inlined here to stay
  ## allocation free (no seq of the 24 fields needed).
  for e in ff:
    let s = u[e] * ff[e].adj
    ff[e].projectTAH s

proc hcForce*[V: static[int], F, SH, SS](w: HcActionWork[V, F, SH, SS],
                                         beta: float, g: HcGauge[V, F],
                                         f: var HcGauge[V, F]) =
  ## The Lie-algebra force: f_l = projectTAH( U_l D_l^dag ) with D_l the
  ## staple sums of `hcActionDeriv`.  Satisfies (see module docs)
  ##   d/ds S(exp(s P) U)|_0 = sum_l redot(P_l, f_l).
  hcActionDeriv(w, beta, g, f)
  let ff = f
  threads:
    for mu in 0..<nDim:
      contractPTAH(g.uA[mu], ff.uA[mu])
      contractPTAH(g.uB[mu], ff.uB[mu])
    for d in 0..<nDiag:
      contractPTAH(g.uD[d], ff.uD[d])

proc hcForce*[V: static[int], F](beta: float, g: HcGauge[V, F],
                                 f: var HcGauge[V, F]) =
  ## Convenience overload: allocates a fresh work object per call.
  var w = newHcActionWork(g)
  hcForce(w, beta, g, f)

# ---------------------------------------------------------------------------
# small helpers for HMC and tests
# ---------------------------------------------------------------------------

proc redot*[V: static[int], F](a, b: HcGauge[V, F]): float =
  ## sum over all 24 link fields and all sites of Re tr(a^dag b)
  var res = 0.0
  threads:
    var s = 0.0
    for mu in 0..<nDim:
      s += redot(a.uA[mu], b.uA[mu])
      s += redot(a.uB[mu], b.uB[mu])
    for d in 0..<nDiag:
      s += redot(a.uD[d], b.uD[d])
    threadMaster: res = s
  res

proc randomTAH*(g: HcGauge, r: var RNGField) =
  ## gaussian traceless anti-Hermitian matrices on all 24 link fields
  ## (HMC momenta); call inside a `threads:` block, like `hcgauge.random`.
  eachLink(g, u):
    u.randomTAH r

when isMainModule:
  import std/monotimes, std/times
  qexInit()
  echo "hcaction: staple recipe check passed (192 terms match hcgeom.triPath)"
  let hl = newHcLayout([4, 4, 4, 4])
  var g = newHcGauge(hl)
  var f = newOneOf(g)
  var w = newHcActionWork(g)
  echo "unit gauge: hcAction = ", hcAction(w, 6.0, g)
  hcForce(w, 6.0, g, f)
  echo "unit gauge: |force|^2 = ", redot(f, f)
  # rough per-call cost on 8^4 cells (for task M planning)
  block:
    let hl8 = newHcLayout([8, 8, 8, 8])
    var r = hl8.lo.newRNGField(RngMilc6, 7777'u64)
    var g8 = newHcGauge(hl8)
    var f8 = newOneOf(g8)
    threads:
      g8.warm(0.35, r)
    var w8 = newHcActionWork(g8)
    hcForce(w8, 6.0, g8, f8)          # warm up comms/buffers
    let n = 20
    let t0 = getMonoTime()
    for i in 0..<n:
      hcForce(w8, 6.0, g8, f8)
    let t1 = getMonoTime()
    for i in 0..<n:
      discard hcAction(w8, 6.0, g8)
    let t2 = getMonoTime()
    echo "8^4 cells: hcForce  ", (t1-t0).inMicroseconds.float/(1e6*n.float),
         " s/call"
    echo "8^4 cells: hcAction ", (t2-t1).inMicroseconds.float/(1e6*n.float),
         " s/call"
  qexFinalize()
