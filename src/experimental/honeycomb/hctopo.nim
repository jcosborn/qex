## Hexagon-clover field strength, energy density and topological charge.
##
## For site x and hexagon h (hcgeom.hexagons, ring order of hexTriPaths):
##   C_h(x)     = sum_{k=1}^{6} P_k(x)                (6 triangle loops, area sqrt3/4)
##   Fhat_ab(x) = sum_h c_h,ab TAH[C_h(x)],   c_h,ab = s Omega^(h)_ab/(4 sqrt3) = -+1/12
## from Fhat_Omega = s (4/sqrt3) TAH[C_h/6] and Fhat_ab = (3/8) sum_h Omega_ab Fhat_Omega.
## s = cloverSign = -1 because the ring runs clockwise w.r.t. omega; with it
## Fhat_ab = +i a^2 F_ab T for U = exp(i int A.dl).  The sign is pinned by the
## weak-field test of tests/ttopo.nim (E and Q are even in Fhat and cannot
## see it); it matters only for the clover term of hcwilson.
## With f_p = Fhat_ab, a > b, in pairIndex order (Fhat anti-Hermitian):
##   E(x) = -(1/2) sum_{ab} Tr Fhat_ab Fhat_ab = sum_{a>b} |Fhat_ab|^2,   avgE = (1/N_sites) sum_x E(x)
##   q(x) = -(1/32 pi^2) eps Tr[Fhat Fhat] = -(1/4 pi^2) [tr F01F23 - tr F02F13 + tr F03F12]
##   sum_x q = (1/4 pi^2) [redot(f0,f5) - redot(f1,f4) + redot(f2,f3)],   Q = (1/2) sum_x q(x)
## the 1/2 being the a^4/2 volume per site.  The same reductions as QEX
## densityE/topoQ, so E and Q are on the cubic footing.
##
## The 192 loop plans are generated from hexTriPaths at load; each loop is a
## 3-link product chained through single-axis shifts (offsets in {-1,0,1}^4),
## associated to minimise shifts.  TopoWork holds all shifters and scratch.

import std/math
import base, layout, field, maths
import physics/qcdTypes
import gauge
import hcgeom, hcgauge

export hcgauge

const cloverSign* = -1.0

type
  TopoFactor = object
    slot: int          ## 0..3 uA[mu], 4..7 uB[mu], 8..23 uD[delta]
    dag: bool
  TopoPlan = object
    ## loop M0(x+o0)^d0 M1(x+o1)^d1 M2(x+o2)^d2:
    ##  l2r: t1 = M0^d0 sh(M1,s1)^d1;  t2 = t1 sh(M2,s2)^d2;  c += sh(t2,sf)
    ##  r2l: t1 = M1^d1 sh(M2,s1)^d2;  t2 = M0^d0 sh(t1,s2);  c += sh(t2,sf)
    m: array[3, TopoFactor]
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
    doAssert c[i] >= -1 and c[i] <= 1, "loop offset component out of range"
    result += abs(c[i])

proc buildTopoPlans(): array[nSubs, array[nHexPerSite, array[6, TopoPlan]]] =
  for sub in 0..<nSubs:
    for h in 0..<nHexPerSite:
      let paths = hexTriPaths(Site(cell: [0, 0, 0, 0], sub: sub), hexagons[h])
      for k in 0..<6:
        let p = paths[k]
        doAssert p.len == 3
        var o: array[3, Cell]
        var pl: TopoPlan
        for j in 0..<3:
          pl.m[j] = TopoFactor(slot: topoSlot(p[j].kind, p[j].idx),
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

let topoPlans = buildTopoPlans()

func pairIndex*(a, b: int): int =
  ## storage of Fhat_ab, a > b: (1,0)->0 (2,0)->1 (2,1)->2 (3,0)->3 (3,1)->4 (3,2)->5
  doAssert a > b
  a*(a-1) div 2 + b

proc buildHexWeights(): array[nHexPerSite, array[3, tuple[p: int, w: float]]] =
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
      result[h][n] = (p: pairIndex(a, b), w: cloverSign*om[a][b]/(4.0*sqrt(3.0)))
      inc n
    doAssert n == 3

let hexWeights = buildHexWeights()

type
  TopoWork*[F, SS] = ref object
    ## shifters, scratch and the 12 Fhat fields; build once outside `threads:`
    links: array[nDirs, F]
    f*: array[nSubs, array[6, F]]   ## Fhat_ab, pairIndex order, per sublattice
    c, t1, t2: F
    sp, sm: array[nDim, SS]         ## +e_mu / -e_mu shifters

proc newTopoWork*[F](g: HcGauge[F]): auto =
  type SS = type(newShifter(g.uA[0], 0, 1))
  var w = TopoWork[F, SS]()
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

template shiftBy(w, xx, vv: untyped): untyped =
  ## field with value xx(z + vv): xx itself or a shifter buffer
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
  ## t := a^{da} b^{db}
  if da:
    if db: t := a.adj * b.adj
    else: t := a.adj * b
  else:
    if db: t := a * b.adj
    else: t := a * b

proc fmunu*[F, SS](w: TopoWork[F, SS], g: HcGauge[F]) =
  ## w.f[sub][pairIndex(a,b)] := Fhat_ab on the A (0) and B (1) sites
  for i in 0..<nDirs:
    w.links[i] = g.links[i]
  threads:
    for sub in 0..<nSubs:
      for p in 0..<6:
        w.f[sub][p] := 0
    for sub in 0..<nSubs:
      for h in 0..<nHexPerSite:
        w.c := 0
        for k in 0..<6:
          let pl = topoPlans[sub][h][k]
          let m0 = w.links[pl.m[0].slot]
          let m1 = w.links[pl.m[1].slot]
          let m2 = w.links[pl.m[2].slot]
          if pl.l2r:
            let b1 = shiftBy(w, m1, pl.s1)
            mulDD(w.t1, m0, b1, pl.m[0].dag, pl.m[1].dag)
            let b2 = shiftBy(w, m2, pl.s2)
            mulDD(w.t2, w.t1, b2, false, pl.m[2].dag)
          else:
            let b1 = shiftBy(w, m2, pl.s1)
            mulDD(w.t1, m1, b1, pl.m[1].dag, pl.m[2].dag)
            let b2 = shiftBy(w, w.t1, pl.s2)
            mulDD(w.t2, m0, b2, pl.m[0].dag, false)
          let res = shiftBy(w, w.t2, pl.sf)
          w.c += res
        let hw = hexWeights[h]
        for e in w.c:
          var mm {.noinit.}: type(load1(w.c[0]))
          mm.projectTAH w.c[e]
          w.f[sub][hw[0].p][e] += hw[0].w*mm
          w.f[sub][hw[1].p][e] += hw[1].w*mm
          w.f[sub][hw[2].p][e] += hw[2].w*mm

proc EQ*[F, SS](w: TopoWork[F, SS], g: HcGauge[F]): tuple[e, q: float] =
  ## (avgE, Q) from the hexagon clover (computes fmunu first)
  fmunu(w, g)
  var ee, qq: float
  threads:
    var es = 0.0
    var qs = 0.0
    for sub in 0..<nSubs:
      for p in 0..<6:
        es += redot(w.f[sub][p], w.f[sub][p])
      qs += redot(w.f[sub][0], w.f[sub][5]) -
            redot(w.f[sub][1], w.f[sub][4]) +
            redot(w.f[sub][2], w.f[sub][3])
    threadMaster:
      ee = es
      qq = qs
  result.e = ee/float(2*g.lo.physVol)
  result.q = 0.5*qq/(4.0*PI*PI)
