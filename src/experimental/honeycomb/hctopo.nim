## Hexagon-clover field strength, energy density and topological charge on the
## 16-cell honeycomb (task **W2**).  FORMULATION.md section 4.1-4.3.
##
## Construction
## ------------
## For every site x (both sublattices) and each of the 16 hexagons h at x,
## `C_h(x) = sum_{k=1}^{6} P_k(x)` is the sum of the 6 triangle loops of the
## hexagon, all based at x and all traversed in `hcgeom.hexTriPaths` ring
## order.  Then
##
##   Fhat_Omega^(h)(x) = sign * (4/sqrt 3) * TAH[ (1/6) C_h(x) ]        (4.1')
##   Fhat_munu(x)      = (3/8) sum_h Omega^(h)_munu Fhat_Omega^(h)(x)   (4.2)
##
## with TAH the traceless anti-Hermitian projection and Omega the unit area
## 2-form of hcgeom.omega.  Since TAH is linear, the whole thing collapses to
##
##   Fhat_munu(x) = sum_h  c_h,munu * TAH[C_h(x)],
##   c_h,munu = sign * (3/8)(4/sqrt3)(1/6) * Omega^(h)_munu
##            = sign * Omega^(h)_munu/(4 sqrt 3)  =  -+ 1/12 .
##
## **The sign.**  hcgeom's hexagon ring `d^-, d^+, e_mu, ...` runs at angles
## 120, 60, 0, -60, ... degrees in the (e_mu, f) plane of `omega`, i.e.
## CLOCKWISE with respect to Omega = e ^ f.  Each triangle loop therefore
## encloses flux  Phi_k = -(sqrt3/4) a^2 F_Omega  (area sqrt3/4, negative
## orientation), so TAH[(1/6)C_h] = -i (sqrt3/4) a^2 F_Omega T for a weak
## Abelian field U = exp(i phi T):  FORMULATION (4.1) needs `sign = -1` for
## Fhat_munu to equal +i a^2 F_munu T, the standard convention in which
## U_l = exp(+i int A.dl).  `hcCloverSign = -1` below; it is validated
## site-by-site against the exact continuum F of a weak plane wave and by the
## Atiyah-Singer constant-flux test (tests/ttopo.nim), both of which would
## come out with ratio -1 if the sign were wrong.
##
## Observables (FORMULATION 4.3)
## -----------------------------
##   E(x)   = -1/2 sum_{mu,nu} Tr[ Fhat_munu Fhat_munu ]
##   avgE   = (1/N_sites) sum_x E(x)                       (intensive, no a^4/2)
##   q(x)   = -(1/(32 pi^2)) eps_{mu nu rho sigma} Tr[ Fhat_munu Fhat_rhosigma ]
##   Q      = (1/2) sum_x q(x)                             (a^4/2 volume per site)
##
## Derivation of the q reduction and its prefactor (independent, do not copy):
##
## 1. Continuum:  Q = (1/(32 pi^2)) int d4x eps Tr[F_munu F_rhosigma] for
##    HERMITIAN F (check: F^a F^a = -2 Tr FF wait, with anti-Hermitian
##    generators T^a normalised Tr T^aT^b = -delta^ab/2 ... we work directly
##    with the matrix identity below).  Our clover Fhat is ANTI-Hermitian,
##    Fhat = i a^2 F_H with F_H Hermitian, so Tr[Fhat Fhat] = -a^4 Tr[F_H F_H]
##    and q picks up an overall minus:  q = -(1/32pi^2) eps Tr[Fhat Fhat].
##    (BPST check: a self-dual Q=1 instanton has int eps Tr[F_H F_H] = 32 pi^2.)
## 2. eps contraction: group the 24 permutations by the pair partition of
##    {0,1,2,3}.  For {01|23} all 8 index assignments give +Tr(F_01 F_23)
##    (each sign flip of eps is compensated by the antisymmetry of F, and the
##    trace is cyclic); eps_0213 = -1 and eps_0312 = +1 give the other two:
##      eps_{munurhosigma} Tr[F_mn F_rs]
##        = 8 [ Tr(F_01 F_23) - Tr(F_02 F_13) + Tr(F_03 F_12) ].
##    Hence, with t(..) = Re tr (Tr of a product of anti-Hermitians is real),
##      q(x) = -(1/(4 pi^2)) [ t(F_01F_23) - t(F_02F_13) + t(F_03F_12) ].
##    This is identical to QEX's cubic `topoQ` (gaugeUtils.nim:1274): its
##    a - b + c uses F_10F_32 etc., where both index swaps cancel.  The old
##    "factor 2 off" suspicion was wrong (settled exactly by task C's
##    Atiyah-Singer test; the mis-remembered pairing was 1/32pi^2 with
##    eps F^a F^a, which actually carries 1/64pi^2).
## 3. In the pair storage f[p] = Fhat_{ab}, a>b (p: (1,0)(2,0)(2,1)(3,0)(3,1)(3,2)):
##      t(F_01F_23) = t(f0 f5),  t(F_02F_13) = t(f1 f4),  t(F_03F_12) = t(f3 f2)
##    and for anti-Hermitian fields  sum_x Re tr(a b) = -redot(a, b), so
##      sum_x q(x) = +(1/(4 pi^2)) [redot(f0,f5) - redot(f1,f4) + redot(f2,f3)]
##    summed over both sublattices, and  Q = (1/2) sum_x q(x).
##
## Implementation notes
## --------------------
## The 2 x 16 x 6 = 192 loop recipes are generated at module load by walking
## `hcgeom.hexTriPaths` (never hand-derived), each as 3 factors
## (link slot, dagger, cell offset).  All offsets have components in {-1,0,1}
## (asserted).  Each loop is evaluated as a field product anchored at
## intermediate offsets, chain-shifting the partial product with single-axis
## shifters (SIMD- and MPI-safe), with the association order (left-to-right or
## right-to-left) chosen per loop to minimise the number of shifts.
## `HcTopoWork` holds all shifters and scratch; after `newHcTopoWork` nothing
## allocates inside a `threads:` block.

import std/math
import base, layout, field, maths
import physics/qcdTypes
import gauge
import hcgeom, hclayout, hcgauge

export hcgauge

const hcCloverSign* = -1.0
  ## overall sign of (4.1) required by hcgeom's (clockwise) hexagon ring
  ## orientation; see module docs.  Validated in tests/ttopo.nim.

# ---------------------------------------------------------------------------
# loop plans, generated from hcgeom.hexTriPaths at module load
# ---------------------------------------------------------------------------

type
  HcTopoFactor = object
    slot: int          ## 0..3 uA[mu], 4..7 uB[mu], 8..23 uD[delta]
    dag: bool
  HcTopoPlan = object
    ## one triangle loop  M0(x+o0)^d0 M1(x+o1)^d1 M2(x+o2)^d2, and its
    ## evaluation plan:
    ##  l2r: t1 = M0^d0 * sh(M1,s1)^d1;  t2 = t1 * sh(M2,s2)^d2;  c += sh(t2,sf)
    ##  r2l: t1 = M1^d1 * sh(M2,s1)^d2;  t2 = M0^d0 * sh(t1,s2);  c += sh(t2,sf)
    ## where sh(x,v)(z) = x(z+v) is a chain of single-axis shifts.
    m: array[3, HcTopoFactor]
    l2r: bool
    s1, s2, sf: Cell

func topoSlot(k: LinkKind, idx: int): int =
  case k
  of lkA: idx
  of lkB: nDim + idx
  of lkD: 2*nDim + idx

func subCell(a, b: Cell): Cell =
  for i in 0..<nDim: result[i] = a[i] - b[i]

func shiftCost(c: Cell): int =
  for i in 0..<nDim:
    doAssert c[i] >= -1 and c[i] <= 1,
      "hctopo: loop offset component out of range"
    result += abs(c[i])

proc buildTopoPlans(): array[nSubs, array[nHexPerSite, array[6, HcTopoPlan]]] =
  for sub in 0..<nSubs:
    for h in 0..<nHexPerSite:
      let paths = hexTriPaths(Site(cell: [0, 0, 0, 0], sub: sub), hexagons[h])
      for k in 0..<6:
        let p = paths[k]
        doAssert p.len == 3
        var o: array[3, Cell]
        var pl: HcTopoPlan
        for j in 0..<3:
          pl.m[j] = HcTopoFactor(slot: topoSlot(p[j].kind, p[j].idx),
                                 dag: p[j].dag)
          o[j] = p[j].cell
        let
          cL = shiftCost(subCell(o[1], o[0])) + shiftCost(subCell(o[2], o[0]))
          cR = shiftCost(subCell(o[2], o[1])) + shiftCost(subCell(o[1], o[0]))
        pl.l2r = cL < cR
        if pl.l2r:
          pl.s1 = subCell(o[1], o[0])
          pl.s2 = subCell(o[2], o[0])
        else:
          pl.s1 = subCell(o[2], o[1])
          pl.s2 = subCell(o[1], o[0])
        pl.sf = o[0]
        discard shiftCost(pl.sf)
        result[sub][h][k] = pl

let hcTopoPlans = buildTopoPlans()

# pair storage order for Fhat_{ab}, a > b
func pairIndex*(a, b: int): int =
  ## (1,0)->0 (2,0)->1 (2,1)->2 (3,0)->3 (3,1)->4 (3,2)->5
  doAssert a > b
  a*(a-1) div 2 + b

proc buildHexWeights(): array[nHexPerSite, array[3, tuple[p: int, w: float]]] =
  ## c_h,ab = hcCloverSign * Omega^(h)_ab / (4 sqrt 3)  for the 3 pairs (a>b)
  ## with a or b equal to the hexagon's axis direction; = -+1/12.
  for h in 0..<nHexPerSite:
    let
      hex = hexagons[h]
      om = omega(hex)
    var n = 0
    for nu in 0..<nDim:
      if nu == hex.mu: continue
      let
        a = max(hex.mu, nu)
        b = min(hex.mu, nu)
      result[h][n] = (p: pairIndex(a, b),
                      w: hcCloverSign*om[a][b]/(4.0*sqrt(3.0)))
      inc n
    doAssert n == 3

let hcHexWeights = buildHexWeights()

# ---------------------------------------------------------------------------
# work space and executor
# ---------------------------------------------------------------------------

type
  HcTopoWork*[V: static[int], F, SS] = ref object
    ## Shifters and scratch for `hcFmunu`/`hcEQ`.  Create once with
    ## `newHcTopoWork` (allocates); afterwards every call is allocation free.
    ## Rebound to the gauge field passed to each call.
    links: array[nDirs, F]          ## flat link table (topoSlot order)
    f*: array[nSubs, array[6, F]]   ## Fhat_{ab}, pairIndex order, per sublattice
    c, t1, t2: F                    ## scratch
    sp, sm: array[nDim, SS]         ## +e_mu / -e_mu shifters

proc newHcTopoWork*[V: static[int], F](g: HcGauge[V, F]): auto =
  ## Allocate shifters, scratch and the 12 Fhat fields.  Call outside `threads:`.
  type SS = type(newShifter(g.uA[0], 0, 1))
  var w = HcTopoWork[V, F, SS]()
  for sub in 0..<nSubs:
    for p in 0..<6:
      w.f[sub][p] = g.uA[0].newOneOf
  w.c = g.uA[0].newOneOf
  w.t1 = g.uA[0].newOneOf
  w.t2 = g.uA[0].newOneOf
  for mu in 0..<nDim:
    w.sp[mu] = newShifter(g.uA[0], mu, 1)
    w.sm[mu] = newShifter(g.uA[0], mu, -1)
  w

proc rebind[V: static[int], F, SS](w: HcTopoWork[V, F, SS], g: HcGauge[V, F]) =
  ## ref assignments: outside `threads:`
  for mu in 0..<nDim:
    w.links[mu] = g.uA[mu]
    w.links[nDim+mu] = g.uB[mu]
  for d in 0..<nDiag:
    w.links[2*nDim+d] = g.uD[d]

template hcShiftBy(w, xx, vv: untyped): untyped =
  ## field with value xx(z + vv); either xx itself or a shifter buffer
  block:
    var cur = xx
    let v = vv
    for smu in 0..<nDim:
      if v[smu] == 1:
        cur = w.sp[smu] ^* cur
      elif v[smu] == -1:
        cur = w.sm[smu] ^* cur
    cur

template mulDD(t, a, b: untyped; da, db: bool) =
  ## t := a^{da} * b^{db}
  if da:
    if db: t := a.adj * b.adj
    else: t := a.adj * b
  else:
    if db: t := a * b.adj
    else: t := a * b

proc hcFmunu*[V: static[int], F, SS](w: HcTopoWork[V, F, SS],
                                     g: HcGauge[V, F]) =
  ## Fill w.f[sub][pairIndex(a,b)] with the hexagon-clover Fhat_{ab}(x) for
  ## the A (sub=0) and B (sub=1) sites of every cell.  Traceless
  ## anti-Hermitian; Fhat ~ +i a^2 F_munu T for weak fields (see module docs).
  rebind(w, g)
  threads:
    for sub in 0..<nSubs:
      for p in 0..<6:
        w.f[sub][p] := 0
    for sub in 0..<nSubs:
      for h in 0..<nHexPerSite:
        w.c := 0
        for k in 0..<6:
          let pl = hcTopoPlans[sub][h][k]
          let m0 = w.links[pl.m[0].slot]
          let m1 = w.links[pl.m[1].slot]
          let m2 = w.links[pl.m[2].slot]
          if pl.l2r:
            let b1 = hcShiftBy(w, m1, pl.s1)
            mulDD(w.t1, m0, b1, pl.m[0].dag, pl.m[1].dag)
            let b2 = hcShiftBy(w, m2, pl.s2)
            mulDD(w.t2, w.t1, b2, false, pl.m[2].dag)
          else:
            let b1 = hcShiftBy(w, m2, pl.s1)
            mulDD(w.t1, m1, b1, pl.m[1].dag, pl.m[2].dag)
            let b2 = hcShiftBy(w, w.t1, pl.s2)
            mulDD(w.t2, m0, b2, pl.m[0].dag, false)
          let res = hcShiftBy(w, w.t2, pl.sf)
          w.c += res
        let hw = hcHexWeights[h]
        for e in w.c:
          var mm {.noinit.}: type(load1(w.c[0]))
          mm.projectTAH w.c[e]
          w.f[sub][hw[0].p][e] += hw[0].w*mm
          w.f[sub][hw[1].p][e] += hw[1].w*mm
          w.f[sub][hw[2].p][e] += hw[2].w*mm

proc hcEQ*[V: static[int], F, SS](w: HcTopoWork[V, F, SS],
                                  g: HcGauge[V, F]): tuple[e, q: float] =
  ## (avgE, Q) from the hexagon clover:
  ##   avgE = (1/N_sites) sum_x E(x),  E(x) = -1/2 sum_{munu} Tr[Fhat Fhat]
  ##   Q    = (1/2) sum_x q(x)         (the 1/2 is the a^4/2 volume per site)
  ## Sign conventions and the q reduction are derived in the module docs.
  hcFmunu(w, g)
  var ee, qq: float
  threads:
    var es = 0.0
    var qs = 0.0
    for sub in 0..<nSubs:
      for p in 0..<6:
        # -sum_x tr(Fhat_ab Fhat_ab) = +redot(f_p, f_p); E counts each (a,b)
        # ordered pair twice, so E(x) = sum_{a>b} |Fhat_ab|^2 exactly.
        es += redot(w.f[sub][p], w.f[sub][p])
      qs += redot(w.f[sub][0], w.f[sub][5]) -
            redot(w.f[sub][1], w.f[sub][4]) +
            redot(w.f[sub][2], w.f[sub][3])
    threadMaster:
      ee = es
      qq = qs
  result.e = ee/float(g.hl.nSites)
  result.q = 0.5*qq/(4.0*PI*PI)

when isMainModule:
  import std/monotimes, std/times
  import rng
  qexInit()
  echo "hctopo: ", 2*nHexPerSite*6, " loop plans built from hcgeom.hexTriPaths"
  block:
    let hl = newHcLayout([4, 4, 4, 4])
    var g = newHcGauge(hl)
    var w = newHcTopoWork(g)
    let (e0, q0) = hcEQ(w, g)
    echo "unit gauge: avgE = ", e0, "  Q = ", q0
    doAssert abs(e0) < 1e-24 and abs(q0) < 1e-24
  block:
    let hl = newHcLayout([8, 8, 8, 8])
    var r = hl.lo.newRNGField(RngMilc6, 13579'u64)
    var g = newHcGauge(hl)
    threads:
      g.warm(0.35, r)
    var w = newHcTopoWork(g)
    let (e1, q1) = hcEQ(w, g)
    echo "warm 8^4: avgE = ", e1, "  Q = ", q1
    let n = 10
    let t0 = getMonoTime()
    for i in 0..<n:
      discard hcEQ(w, g)
    let t1 = getMonoTime()
    echo "8^4 cells: hcEQ ", (t1-t0).inMicroseconds.float/(1e6*n.float),
         " s/call"
  qexFinalize()
