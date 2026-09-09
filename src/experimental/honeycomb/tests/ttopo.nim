## `hctopo.nim` (hexagon clover F_munu, E, Q).
##
## 1. brute-force reference: an independent single-site implementation of
##    FORMULATION 4.1-4.3 (hexTriPaths + scalar complex algebra, literal
##    (4/sqrt3), (3/8), Omega sums -- no code shared with hctopo's collapsed
##    executor) must agree with `fmunu`/`EQ` to 1e-12;
## 2. weak-field validation: for an exact-line-integral Abelian plane wave the
##    clover Fhat_munu(x) must equal +i a^2 F_munu(x) T site by site with
##    O(p^2 a^2) errors (factor ~4 between k=1 and k=2) -- this validates the
##    4/sqrt3, the 3/8, the hexagon orientations and cloverSign at once;
## 3. the same for E;
## 4. THE DECISIVE q-prefactor test (exact, Atiyah-Singer): a constant-field-
##    strength Abelian configuration in the Cartan direction T = diag(1,-1,0)
##    with integer fluxes (n1, n2) in the (0,1) and (2,3) planes has
##    Q_exact = sum_i q_i^2 n1 n2 = 2 n1 n2; the honeycomb clover must
##    reproduce it up to the O(phi^2) ~ 1/L^4 clover artifact.  A missing
##    a^4/2 site-volume factor would show up as exactly 2x.
##    (helpers.setFluxHc, transferred from task C's refCubicMeas -abeliantest.)
## 5. gauge invariance of E and Q;
## 6. integer-Q clustering (secondary, per task C's warning the plain clover
##    at coarse spacing is loose): warm rough 8^4 configs flowed deep must
##    give Q clustering at integers.

import std/[math, strformat, unittest]
import qex except epsilon
import physics/qcdTypes
import ../hcgeom
import ../hcgauge
import ../hcaction
import ../hcflow
import ../hctopo
import helpers

qexInit()

# ---------------------------------------------------------------------------
# brute-force clover per FORMULATION 4.1-4.3 (literal, uncollapsed)
# ---------------------------------------------------------------------------

proc fmunuRefSite(g: auto, y: Cell, sub: int): array[6, Mat] =
  ## Fhat_{ab} (a>b, hctopo.pairIndex order) at one site, straight from
  ## hexTriPaths / omega:  Fhat_Omega = cloverSign*(4/sqrt3)*TAH[(1/6)C_h],
  ## Fhat_ab = (3/8) sum_h Omega_ab Fhat_Omega.
  let
    lo = g.lo
    geom = g.lo.physGeom
  for h in 0..<nHexPerSite:
    var c: Mat
    for path in hexTriPaths(Site(cell: y, sub: sub), hexagons[h]):
      var p = mid()
      for l in path:
        var m = getMat(g.link(l.kind, l.idx), siteIndex(lo, l.cell, geom))
        if l.dag: m = mdag(m)
        p = mmul(p, m)
      c = madd(c, p)
    let fom = mscale(mTAH(mscale(c, 1.0/6.0)), cloverSign*4.0/sqrt(3.0))
    let om = omega(hexagons[h])
    for a in 1..<nDim:
      for b in 0..<a:
        if om[a][b] != 0.0:
          result[pairIndex(a, b)] =
            madd(result[pairIndex(a, b)], mscale(fom, 0.375*om[a][b]))

proc eqRef(g: auto): tuple[e, q: float, f: seq[array[6, Mat]]] =
  ## brute-force avgE, Q and the full Fhat table (index 2*site + sub)
  let
    lo = g.lo
    geom = g.lo.physGeom
  var fs = newSeq[array[6, Mat]](2*lo.nSites)
  var esum = 0.0
  var qsum = 0.0
  for y in lexCells(geom):
    let i = siteIndex(lo, y, geom)
    for sub in 0..1:
      let f = fmunuRefSite(g, y, sub)
      fs[2*i+sub] = f
      for p in 0..<6:
        esum -= reTr(mmul(f[p], f[p]))
      # q(x) = -(1/4pi^2)[t(F01 F23) - t(F02 F13) + t(F03 F12)]
      #      = -(1/4pi^2)[t(f0 f5) - t(f1 f4) + t(f3 f2)]   (f_p = Fhat_{ab}, a>b)
      qsum -= (reTr(mmul(f[0], f[5])) - reTr(mmul(f[1], f[4])) +
               reTr(mmul(f[3], f[2])))/(4.0*PI*PI)
  result.e = esum/float(g.lo.physVol*2)
  result.q = 0.5*qsum
  result.f = fs

# ---------------------------------------------------------------------------

suite "hctopo":

  test "1. field implementation == brute-force FORMULATION 4.1-4.3":
    let hl = newLayout([4, 4, 4, 6])
    let lo = hl
    var r = lo.newRNGField(RngMilc6, 97531'u64)
    var g = newHcGauge(hl)
    threads:
      g.random r
    var w = newTopoWork(g)
    let (e, q) = EQ(w, g)
    let (er, qr, fr) = eqRef(g)
    var fdev = 0.0
    for i in 0..<lo.nSites:
      for sub in 0..1:
        for p in 0..<6:
          fdev = max(fdev, mmaxdiff(getMat(w.f[sub][p], i), fr[2*i+sub][p]))
    echo &"  random config: avgE field {e:.15g}  ref {er:.15g}"
    echo &"                 Q    field {q:.15g}  ref {qr:.15g}"
    echo &"  max site-wise |Fhat(field) - Fhat(ref)| = {fdev:.3e}"
    ok(&"Fhat matches brute force ({fdev:.2e})", fdev < 1e-12)
    ok(&"avgE matches brute force (rel {abs(e-er)/er:.2e})",
       abs(e-er) < 1e-12*abs(er))
    ok(&"Q matches brute force ({abs(q-qr):.2e})",
       abs(q-qr) < 1e-12*(abs(qr)+1.0))
    # and once on a V=1 layout (SIMD-lane check)
    let hl1 = newLayout([4, 4, 4, 6], 1)
    var r1 = hl1.newRNGField(RngMilc6, 97531'u64)
    var g1 = newHcGauge(hl1)
    threads:
      g1.random r1
    var w1 = newTopoWork(g1)
    let (e1, q1) = EQ(w1, g1)
    let (er1, qr1, _) = eqRef(g1)
    ok(&"V=1 layout: avgE (rel {abs(e1-er1)/er1:.2e}), Q ({abs(q1-qr1):.2e})",
       abs(e1-er1) < 1e-12*abs(er1) and abs(q1-qr1) < 1e-12*(abs(qr1)+1.0))
    # and on [2,4,4,6]: L=2 makes the +-1 shift chains wrap
    let hl2 = newLayout([2, 4, 4, 6])
    var r2 = hl2.newRNGField(RngMilc6, 13570'u64)
    var g2 = newHcGauge(hl2)
    threads:
      g2.random r2
    var w2 = newTopoWork(g2)
    let (e2, q2) = EQ(w2, g2)
    let (er2, qr2, _) = eqRef(g2)
    ok(&"[2,4,4,6] (shift wraps): avgE (rel {abs(e2-er2)/er2:.2e}), " &
       &"Q ({abs(q2-qr2):.2e})",
       abs(e2-er2) < 1e-12*abs(er2) and abs(q2-qr2) < 1e-12*(abs(qr2)+1.0))

  test "2. weak-field Fhat and E vs exact continuum (O(p^2 a^2) scaling)":
    # A_mu = eps_mu cos(p.x) embedded via T: exact F_ab(x) = -(eps_b p_a -
    # eps_a p_b) sin(p.x); the clover must give Fhat_ab = +i F_ab T site by
    # site.  This validates 4/sqrt3, 3/8, ring orientation and cloverSign.
    const eps0 = 1e-4
    proc weakCase(ns, k: int): tuple[fdev, erat: float] =
      var pv, eps: array[4, float]
      pv[0] = 2.0*PI*float(k)/float(ns)
      for mu in 0..<4: eps[mu] = eps0*ehat[mu]
      let hl = newLayout([ns, ns, ns, ns])
      let lo = hl
      var g = newHcGauge(hl)
      setAbelianHc(g, eps, pv)
      var w = newTopoWork(g)
      let eMeas = EQ(w, g).e
      var fdev = 0.0
      var fscale = 0.0
      var eex = 0.0
      for i in lo.sites:
        var y: array[4, float]
        for mu in 0..<4: y[mu] = lo.coords[mu][i].float
        for sub in 0..1:
          var x: array[4, float]
          for mu in 0..<4: x[mu] = y[mu] + 0.5*sub.float
          var px = 0.0
          for mu in 0..<4: px += pv[mu]*x[mu]
          var esite = 0.0
          for a in 1..<4:
            for b in 0..<a:
              let fex = -(eps[b]*pv[a] - eps[a]*pv[b])*sin(px)
              esite += 2.0*fex*fex
              # measured phase: Fhat = i phi T  =>  phi = Im [0,0] element
              let m = getMat(w.f[sub][pairIndex(a, b)], i)
              fdev = max(fdev, abs(m[0][0][1] - fex))
              fscale = max(fscale, abs(fex))
              # structure: [1,1] = -i phi, off-diag and real parts ~ 0
              fdev = max(fdev, abs(m[1][1][1] + fex))
              fdev = max(fdev, abs(m[2][2][1]))
              fdev = max(fdev, abs(m[0][0][0]))
              fdev = max(fdev, abs(m[0][1][0]) + abs(m[0][1][1]))
          eex += esite
      eex = eex/float(2*hl.physVol)
      (fdev/fscale, eMeas/eex)
    let (f1, e1) = weakCase(12, 1)
    let (f2, e2) = weakCase(12, 2)
    echo &"  Ns=12 k=1: max site |Fhat - i F_exact T|/max|F| = {f1:.5f}   E/E_exact = {e1:.6f}"
    echo &"  Ns=12 k=2: max site |Fhat - i F_exact T|/max|F| = {f2:.5f}   E/E_exact = {e2:.6f}"
    echo &"  scaling k=1 -> k=2:  F dev ratio {f2/f1:.3f} ~ 4,  ",
         &"(1-E ratio): {(1.0-e2)/(1.0-e1):.3f} ~ 4"
    ok(&"Fhat correct incl. sign at k=1 (rel dev {f1:.4f} < 0.05)", f1 < 0.05)
    ok(&"F dev scales as O(p^2) (ratio {f2/f1:.2f} in [3,5.5])",
       f2/f1 > 3.0 and f2/f1 < 5.5)
    ok(&"E/E_exact -> 1 (k=1: {e1:.4f})", abs(e1 - 1.0) < 0.05)
    ok(&"E dev scales as O(p^2) (ratio {(1.0-e2)/(1.0-e1):.2f} in [3,5.5])",
       (1.0-e2)/(1.0-e1) > 3.0 and (1.0-e2)/(1.0-e1) < 5.5)

  test "3. DECISIVE: Atiyah-Singer constant-flux Q = 2 n1 n2 (exact)":
    echo "  L    n1 n2  Q_exact   Q_measured      Q/Q_exact     E/E_cont"
    proc fluxCase(lsize, n1, n2: int): tuple[qrat, erat: float] =
      let hl = newLayout([lsize, lsize, lsize, lsize])
      var g = newHcGauge(hl)
      setFluxHc(g, n1, n2)
      var w = newTopoWork(g)
      let (e, q) = EQ(w, g)
      let
        f1 = 2.0*PI*float(n1)/float(lsize*lsize)
        f2 = 2.0*PI*float(n2)/float(lsize*lsize)
        qex = 2.0*float(n1*n2)
        eex = 2.0*(f1*f1 + f2*f2)
      echo &"  {lsize:2d}  {n1:3d} {n2:2d}  {qex:8.4f}  {q:.10f}  {q/qex:.8f}  {e/eex:.8f}"
      (q/qex, e/eex)
    let (q44, _) = fluxCase(4, 1, 1)
    let (q88, e88) = fluxCase(8, 1, 1)
    let (q12, _) = fluxCase(12, 1, 1)
    let (qm, _) = fluxCase(8, 1, -1)
    let (q21, _) = fluxCase(8, 2, 1)
    ok(&"Q/2n1n2 = 1 to 2% at L=8 ({q88:.6f})", abs(q88 - 1.0) < 0.02)
    ok(&"Q/2n1n2 = 1 to 1% at L=12 ({q12:.6f})", abs(q12 - 1.0) < 0.01)
    ok(&"Q odd under n2 -> -n2 (ratio {qm:.6f})", abs(qm - 1.0) < 0.02)
    ok(&"fluxes (2,1): Q/4 = {q21:.6f}", abs(q21 - 1.0) < 0.03)
    ok(&"E/E_cont = 1 at L=8 ({e88:.6f})", abs(e88 - 1.0) < 0.02)
    let s48 = (1.0-q44)/(1.0-q88)
    let s812 = (1.0-q88)/(1.0-q12)
    echo &"  artifact scaling: dev(L=4)/dev(L=8) = {s48:.2f} (~16 for 1/L^4),",
         &"  dev(L=8)/dev(L=12) = {s812:.2f} (~5.1)"
    ok("clover artifact shrinks with L", s48 > 4.0 and s812 > 2.0)

  test "4. gauge invariance of E and Q":
    let hl = newLayout([4, 4, 4, 6])
    let lo = hl
    var r = lo.newRNGField(RngMilc6, 8642097'u64)
    var g = newHcGauge(hl)
    threads:
      g.warm(0.5, r)
    var w = newTopoWork(g)
    let (e0, q0) = EQ(w, g)
    var vA = lo.ColorMatrix(nc)
    var vB = lo.ColorMatrix(nc)
    threads:
      vA.randomSU r
      vB.randomSU r
    g.gaugeTransform(vA, vB)
    let (e1, q1) = EQ(w, g)
    echo &"  E before {e0:.15g}  after {e1:.15g}  (rel {abs(e1-e0)/e0:.2e})"
    echo &"  Q before {q0:.15g}  after {q1:.15g}  (abs {abs(q1-q0):.2e})"
    ok(&"E gauge invariant (rel {abs(e1-e0)/e0:.2e})", abs(e1-e0) < 1e-12*e0)
    ok(&"Q gauge invariant (abs {abs(q1-q0):.2e})",
       abs(q1-q0) < 1e-12*(abs(q0)+1.0))

  test "5. flow recovers the exact integer Q of a perturbed known sector":
    # Start from the exact constant-flux configuration of sector Q = 2 n1 n2,
    # kick every link with multiplicative random-algebra noise U -> exp(eta X) U
    # (strong enough that the raw clover Q is visibly off), flow, and demand
    # Q returns to the KNOWN integer.  Unlike generic rough configs this has
    # topological structure at scale L, so it must anneal cleanly -- and a
    # factor-2 (or 1/2) normalisation error would land on 2 Q_exact (Q_exact/2)
    # instead.
    let hl = newLayout([8, 8, 8, 8])
    var r = hl.newRNGField(RngMilc6, 555777999'u64)
    var g = newHcGauge(hl)
    var p = newOneOf(g)
    var wt = newTopoWork(g)
    const eta = 0.35
    let sectors = [(0, 0), (1, 1), (1, -1), (1, 1), (2, 1), (1, 1), (0, 0),
                   (1, -1)]
    echo &"  eta = {eta}   cfg: (n1,n2)  Q_exact   Q(t=0)     Q(t=1)     Q(t=2)"
    var maxdev = 0.0
    var maxplat = 0.0
    for icfg in 0..<sectors.len:
      let (n1, n2) = sectors[icfg]
      setFluxHc(g, n1, n2)
      threads:
        p.randomTAH r
      expMul(g, g, p, eta)
      let qex = 2.0*float(n1*n2)
      let q0 = EQ(wt, g).q
      var q1, q2: float
      var nstep = 0
      g.flow(0.02, cflow):
        inc nstep
        if nstep == 50:
          q1 = EQ(wt, g).q
        elif nstep == 100:
          q2 = EQ(wt, g).q
          break
      echo &"  {icfg:3d}  ({n1:2d},{n2:2d})  {qex:7.2f}  {q0:9.4f}  {q1:9.4f}  {q2:9.4f}"
      maxdev = max(maxdev, abs(q2 - qex))
      maxplat = max(maxplat, abs(q2 - q1))
    ok(&"flowed Q returns to the exact sector integer " &
       &"(max |Q - Q_exact| = {maxdev:.4f} < 0.05)", maxdev < 0.05)
    ok(&"Q plateaus in t (max |Q(2)-Q(1)| = {maxplat:.4f} < 0.02)",
       maxplat < 0.02)

  test "6. integer-Q clustering on flowed rough configs (secondary, loose)":
    # Generic rough (warm 0.85) configs on 8^4, flowed to t = 10 (200 RK3
    # steps of 0.05; larger steps change which dislocations annihilate).
    # Their lumps sit at the cutoff scale, so -- exactly as task C found on
    # the cubic side -- Q keeps drifting through dislocations and the
    # clustering is LOOSE; the sharp normalisation statements are tests 3
    # and 5.  Here we check the distribution favours integers over
    # half-odd-integers (the factor-1/2 failure mode).
    let hl = newLayout([8, 8, 8, 8])
    const ncfg = 2
    var r = hl.newRNGField(RngMilc6, 1122334455'u64)
    var g = newHcGauge(hl)
    var wt = newTopoWork(g)
    var qend: seq[float]
    echo "  cfg   Q(t=2)    Q(t=5)    Q(t=10)   dist(Q(10),int)"
    for icfg in 0..<ncfg:
      threads:
        g.warm(0.85, r)
      var q1, q2, q3: float
      var nstep = 0
      g.flow(0.05, cflow):
        inc nstep
        if nstep == 40:
          q1 = EQ(wt, g).q
        elif nstep == 100:
          q2 = EQ(wt, g).q
        elif nstep == 200:
          q3 = EQ(wt, g).q
          break
      qend.add q3
      echo &"  {icfg:3d}  {q1:8.4f}  {q2:8.4f}  {q3:8.4f}   {abs(q3-round(q3)):.4f}"
    var dInt = 0.0
    var dHalfOdd = 0.0
    var nnz = 0
    for i in 0..<ncfg:
      let di = abs(qend[i]-round(qend[i]))
      dInt += di
      dHalfOdd += 0.5 - di
      if abs(qend[i]) > 0.5: inc nnz
    dInt = dInt/ncfg.float
    dHalfOdd = dHalfOdd/ncfg.float
    echo &"  mean |Q - nearest integer|      = {dInt:.4f}   (uniform: 0.25)"
    echo &"  mean |Q - nearest half-odd int| = {dHalfOdd:.4f}"
    echo &"  {nnz}/{ncfg} configs with |Q| > 0.5"
    ok(&"nonzero sectors sampled ({nnz}/{ncfg})", nnz > 0)
    ok(&"Q clusters at integers (mean dist {dInt:.3f} < 0.2; uniform 0.25)",
       dInt < 0.2)
    ok(&"distribution favours integers over half-odd-integers " &
       &"({dInt:.3f} < {dHalfOdd:.3f})", dInt < dHalfOdd)

qexFinalize()
