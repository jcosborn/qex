## Task W test, part 2: `hctopo.nim` (hexagon clover F_munu, E, Q).
##
##   export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
##   cd build_mac && make src/experimental/honeycomb/tests/ttopo.nim && \
##     OMP_NUM_THREADS=4 ./bin/ttopo
##
## 1. brute-force reference: an independent single-site implementation of
##    FORMULATION 4.1-4.3 (hexTriPaths + hand-written 3x3 complex algebra,
##    literal (4/sqrt3), (3/8), Omega sums -- no code shared with hctopo's
##    collapsed executor) must agree with `hcFmunu`/`hcEQ` to 1e-12;
## 2. weak-field validation: for an exact-line-integral Abelian plane wave the
##    clover Fhat_munu(x) must equal +i a^2 F_munu(x) T site by site with
##    O(p^2 a^2) errors (factor ~4 between k=1 and k=2) -- this validates the
##    4/sqrt3, the 3/8, the hexagon orientations and hcCloverSign at once;
## 3. the same for E;
## 4. THE DECISIVE q-prefactor test (exact, Atiyah-Singer): a constant-field-
##    strength Abelian configuration in the Cartan direction T = diag(1,-1,0)
##    with integer fluxes (n1, n2) in the (0,1) and (2,3) planes has
##    Q_exact = sum_i q_i^2 n1 n2 = 2 n1 n2; the honeycomb clover must
##    reproduce it up to the O(phi^2) ~ 1/L^4 clover artifact.  A missing
##    a^4/2 site-volume factor would show up as exactly 2x.
##    (Construction transferred from task C's refCubicMeas -abeliantest:
##    U_l = exp(i phi_l T), phi_l the exact line integral of A_1 = f1 x0,
##    A_3 = f2 x2 in the fundamental domain, plus the transition-function
##    correction -f1 L0 x1(endpoint) / -f2 L2 x3(endpoint) for links crossing
##    the x0 / x2 boundary -- on the honeycomb the diagonal links cross with
##    nonzero transverse displacement, whence the endpoint rule; it reduces
##    to task C's construction on the cubic lattice.)
## 5. gauge invariance of E and Q;
## 6. integer-Q clustering (secondary, per task C's warning the plain clover
##    at coarse spacing is loose): warm rough 6^4 configs flowed deep must
##    give Q clustering at integers, stable in t.

import std/[math, strformat, unittest]
import qex except epsilon
import physics/qcdTypes
import ../hcgeom
import ../hclayout
import ../hcgauge
import ../hcaction
import ../hcflow
import ../hctopo

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

template ok(msg: string, cond: bool) =
  let c = cond
  if c: echo "PASS: ", msg
  else: echo "FAIL: ", msg
  check c

qexInit()

const nc = getDefaultNc()

# ---------------------------------------------------------------------------
# independent 3x3 complex matrix machinery (pattern of tests/tgauge.nim)
# ---------------------------------------------------------------------------

type Mat = array[nc, array[nc, array[2, float]]]   ## [row][col][re/im]

proc mid(): Mat =
  for i in 0..<nc: result[i][i][0] = 1.0

proc mmul(a, b: Mat): Mat =
  for i in 0..<nc:
    for j in 0..<nc:
      var re = 0.0
      var im = 0.0
      for k in 0..<nc:
        re += a[i][k][0]*b[k][j][0] - a[i][k][1]*b[k][j][1]
        im += a[i][k][0]*b[k][j][1] + a[i][k][1]*b[k][j][0]
      result[i][j][0] = re
      result[i][j][1] = im

proc mdag(a: Mat): Mat =
  for i in 0..<nc:
    for j in 0..<nc:
      result[i][j][0] = a[j][i][0]
      result[i][j][1] = -a[j][i][1]

proc madd(a, b: Mat): Mat =
  for i in 0..<nc:
    for j in 0..<nc:
      result[i][j][0] = a[i][j][0] + b[i][j][0]
      result[i][j][1] = a[i][j][1] + b[i][j][1]

proc mscale(a: Mat, s: float): Mat =
  for i in 0..<nc:
    for j in 0..<nc:
      result[i][j][0] = s*a[i][j][0]
      result[i][j][1] = s*a[i][j][1]

proc reTr(a: Mat): float =
  for i in 0..<nc: result += a[i][i][0]

proc mTAH(a: Mat): Mat =
  ## (a - a^dag)/2 - tr/nc
  let d = mdag(a)
  var tre = 0.0
  var tim = 0.0
  for i in 0..<nc:
    for j in 0..<nc:
      result[i][j][0] = 0.5*(a[i][j][0] - d[i][j][0])
      result[i][j][1] = 0.5*(a[i][j][1] - d[i][j][1])
  for i in 0..<nc:
    tre += result[i][i][0]
    tim += result[i][i][1]
  for i in 0..<nc:
    result[i][i][0] -= tre/nc.float
    result[i][i][1] -= tim/nc.float

proc mmaxdiff(a, b: Mat): float =
  for i in 0..<nc:
    for j in 0..<nc:
      result = max(result, abs(a[i][j][0]-b[i][j][0]))
      result = max(result, abs(a[i][j][1]-b[i][j][1]))

template toF(x: untyped): float =
  block:
    var v: float
    v := x
    v

proc getMat(f: auto, idx: int): Mat =
  for a in 0..<nc:
    for b in 0..<nc:
      result[a][b][0] = toF f{idx}[a, b].re
      result[a][b][1] = toF f{idx}[a, b].im

proc wrapCoord(c: Cell, geom: openArray[int]): array[4, cint] =
  for mu in 0..<nDim:
    result[mu] = cint(((c[mu] mod geom[mu]) + geom[mu]) mod geom[mu])

proc siteIndex(lo: Layout, c: Cell, geom: openArray[int]): int =
  var cc = wrapCoord(c, geom)
  let ri = lo.rankIndex(cc)
  doAssert ri.rank == lo.myRank, "ttopo is a single-rank test"
  ri.index

iterator lexCells(geom: openArray[int]): Cell =
  var c: Cell
  for x3 in 0..<geom[3]:
    for x2 in 0..<geom[2]:
      for x1 in 0..<geom[1]:
        for x0 in 0..<geom[0]:
          c = [x0, x1, x2, x3]
          yield c

# ---------------------------------------------------------------------------
# brute-force clover per FORMULATION 4.1-4.3 (literal, uncollapsed)
# ---------------------------------------------------------------------------

proc fmunuRefSite(g: auto, y: Cell, sub: int): array[6, Mat] =
  ## Fhat_{ab} (a>b, hctopo.pairIndex order) at one site, straight from
  ## hexTriPaths / omega:  Fhat_Omega = hcCloverSign*(4/sqrt3)*TAH[(1/6)C_h],
  ## Fhat_ab = (3/8) sum_h Omega_ab Fhat_Omega.
  let
    lo = g.hl.lo
    geom = g.hl.geom
  for h in 0..<nHexPerSite:
    var c: Mat
    for path in hexTriPaths(Site(cell: y, sub: sub), hexagons[h]):
      var p = mid()
      for l in path:
        var m = getMat(g.link(l.kind, l.idx), siteIndex(lo, l.cell, geom))
        if l.dag: m = mdag(m)
        p = mmul(p, m)
      c = madd(c, p)
    let fom = mscale(mTAH(mscale(c, 1.0/6.0)), hcCloverSign*4.0/sqrt(3.0))
    let om = omega(hexagons[h])
    for a in 1..<nDim:
      for b in 0..<a:
        if om[a][b] != 0.0:
          result[pairIndex(a, b)] =
            madd(result[pairIndex(a, b)], mscale(fom, 0.375*om[a][b]))

proc eqRef(g: auto): tuple[e, q: float, f: seq[array[6, Mat]]] =
  ## brute-force avgE, Q and the full Fhat table (index 2*site + sub)
  let
    lo = g.hl.lo
    geom = g.hl.geom
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
  result.e = esum/float(g.hl.nSites)
  result.q = 0.5*qsum
  result.f = fs

# ---------------------------------------------------------------------------
# Abelian configurations
# ---------------------------------------------------------------------------

proc phaseIntegral(a, b: float): float =
  if abs(b) < 1e-12: cos(a)
  else: (sin(a+b) - sin(a))/b

proc linkPhase(x, n, eps, pv: array[4, float]): float =
  var a = 0.0
  var b = 0.0
  var en = 0.0
  for mu in 0..<4:
    a += pv[mu]*x[mu]
    b += pv[mu]*n[mu]
    en += eps[mu]*n[mu]
  en*phaseIntegral(a, b)

template setPhase(m: untyped, th: float) =
  ## m := exp(i th T), T = diag(1,-1,0); assumes m starts as the identity
  m[0, 0].re := cos(th)
  m[0, 0].im := sin(th)
  m[1, 1].re := cos(th)
  m[1, 1].im := -sin(th)

proc setAbelianHc(hg: auto, eps, pv: array[4, float]) =
  ## plane wave A_mu = eps_mu cos(p.x), exact line integrals
  hg.unit
  let lo = hg.lo
  for i in lo.sites:
    var y, yb: array[4, float]
    for mu in 0..<4:
      y[mu] = lo.coords[mu][i].float
      yb[mu] = y[mu] + 0.5
    for mu in 0..<4:
      var n: array[4, float]
      n[mu] = 1.0
      setPhase(hg.uA[mu]{i}, linkPhase(y, n, eps, pv))
      setPhase(hg.uB[mu]{i}, linkPhase(yb, n, eps, pv))
    for d in 0..<nDiag:
      var n: array[4, float]
      for mu in 0..<4:
        n[mu] = float((d shr mu) and 1) - 0.5
      setPhase(hg.uD[d]{i}, linkPhase(yb, n, eps, pv))

proc fluxPhase(x, n: array[4, float]; f1, f2: float; l0, l2: int): float =
  ## exact line integral of A_1 = f1 x0, A_3 = f2 x2 from x to x+n in the
  ## fundamental domain, with the transition-function corrections for paths
  ## crossing the x0 = L0 or x2 = L2 boundary (endpoint rule, see header)
  result = n[1]*f1*(x[0] + 0.5*n[0]) + n[3]*f2*(x[2] + 0.5*n[2])
  if x[0] + n[0] >= float(l0) - 1e-9:
    result -= f1*float(l0)*(x[1] + n[1])
  if x[2] + n[2] >= float(l2) - 1e-9:
    result -= f2*float(l2)*(x[3] + n[3])

proc setFluxHc(hg: auto, n1, n2: int) =
  ## constant field strength F_01 = 2 pi n1/(L0 L1), F_23 = 2 pi n2/(L2 L3)
  hg.unit
  let
    lo = hg.lo
    geom = hg.hl.geom
    f1 = 2.0*PI*float(n1)/float(geom[0]*geom[1])
    f2 = 2.0*PI*float(n2)/float(geom[2]*geom[3])
  for i in lo.sites:
    var y, yb: array[4, float]
    for mu in 0..<4:
      y[mu] = lo.coords[mu][i].float
      yb[mu] = y[mu] + 0.5
    for mu in 0..<4:
      var n: array[4, float]
      n[mu] = 1.0
      setPhase(hg.uA[mu]{i}, fluxPhase(y, n, f1, f2, geom[0], geom[2]))
      setPhase(hg.uB[mu]{i}, fluxPhase(yb, n, f1, f2, geom[0], geom[2]))
    for d in 0..<nDiag:
      var n: array[4, float]
      for mu in 0..<4:
        n[mu] = float((d shr mu) and 1) - 0.5
      setPhase(hg.uD[d]{i}, fluxPhase(yb, n, f1, f2, geom[0], geom[2]))

# ---------------------------------------------------------------------------

suite "hctopo":

  test "1. field implementation == brute-force FORMULATION 4.1-4.3":
    let hl = newHcLayout([4, 4, 4, 6])
    let lo = hl.lo
    var r = lo.newRNGField(RngMilc6, 97531'u64)
    var g = newHcGauge(hl)
    threads:
      g.random r
    var w = newHcTopoWork(g)
    let (e, q) = hcEQ(w, g)
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
    let hl1 = newHcLayoutX([4, 4, 4, 6], 1)
    var r1 = hl1.lo.newRNGField(RngMilc6, 97531'u64)
    var g1 = newHcGauge(hl1)
    threads:
      g1.random r1
    var w1 = newHcTopoWork(g1)
    let (e1, q1) = hcEQ(w1, g1)
    let (er1, qr1, fr1) = eqRef(g1)
    discard fr1
    ok(&"V=1 layout: avgE (rel {abs(e1-er1)/er1:.2e}), Q ({abs(q1-qr1):.2e})",
       abs(e1-er1) < 1e-12*abs(er1) and abs(q1-qr1) < 1e-12*(abs(qr1)+1.0))
    # and on [2,4,4,6]: L=2 makes the +-1 shift chains wrap
    let hl2 = newHcLayout([2, 4, 4, 6])
    var r2 = hl2.lo.newRNGField(RngMilc6, 13570'u64)
    var g2 = newHcGauge(hl2)
    threads:
      g2.random r2
    var w2 = newHcTopoWork(g2)
    let (e2, q2) = hcEQ(w2, g2)
    let (er2, qr2, fr2) = eqRef(g2)
    discard fr2
    ok(&"[2,4,4,6] (shift wraps): avgE (rel {abs(e2-er2)/er2:.2e}), " &
       &"Q ({abs(q2-qr2):.2e})",
       abs(e2-er2) < 1e-12*abs(er2) and abs(q2-qr2) < 1e-12*(abs(qr2)+1.0))

  test "2. weak-field Fhat and E vs exact continuum (O(p^2 a^2) scaling)":
    # A_mu = eps_mu cos(p.x) embedded via T: exact F_ab(x) = -(eps_b p_a -
    # eps_a p_b) sin(p.x); the clover must give Fhat_ab = +i F_ab T site by
    # site.  This validates 4/sqrt3, 3/8, ring orientation and hcCloverSign.
    const eps0 = 1e-4
    let ehat = [0.0, 1.0, -0.6, 0.3]
    proc weakCase(ns, k: int): tuple[fdev, erat: float] =
      var pv, eps: array[4, float]
      pv[0] = 2.0*PI*float(k)/float(ns)
      for mu in 0..<4: eps[mu] = eps0*ehat[mu]
      let hl = newHcLayout([ns, ns, ns, ns])
      let lo = hl.lo
      var g = newHcGauge(hl)
      setAbelianHc(g, eps, pv)
      var w = newHcTopoWork(g)
      let (eMeas, qMeas) = hcEQ(w, g)
      discard qMeas
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
      eex = eex/float(hl.nSites)
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
      let hl = newHcLayout([lsize, lsize, lsize, lsize])
      var g = newHcGauge(hl)
      setFluxHc(g, n1, n2)
      var w = newHcTopoWork(g)
      let (e, q) = hcEQ(w, g)
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
    let hl = newHcLayout([4, 4, 4, 6])
    let lo = hl.lo
    var r = lo.newRNGField(RngMilc6, 8642097'u64)
    var g = newHcGauge(hl)
    threads:
      g.warm(0.5, r)
    var w = newHcTopoWork(g)
    let (e0, q0) = hcEQ(w, g)
    var vA = lo.ColorMatrix(nc)
    var vB = lo.ColorMatrix(nc)
    threads:
      vA.randomSU r
      vB.randomSU r
    g.gaugeTransform(vA, vB)
    let (e1, q1) = hcEQ(w, g)
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
    let hl = newHcLayout([8, 8, 8, 8])
    var r = hl.lo.newRNGField(RngMilc6, 555777999'u64)
    var g = newHcGauge(hl)
    var p = newOneOf(g)
    var wt = newHcTopoWork(g)
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
      # U := exp(eta X) U on all 24 link fields
      threads:
        for mu in 0..<nDim:
          for e in g.uA[mu]:
            g.uA[mu][e] := exp(eta*p.uA[mu][e])*g.uA[mu][e]
          for e in g.uB[mu]:
            g.uB[mu][e] := exp(eta*p.uB[mu][e])*g.uB[mu][e]
        for d in 0..<nDiag:
          for e in g.uD[d]:
            g.uD[d][e] := exp(eta*p.uD[d][e])*g.uD[d][e]
      let qex = 2.0*float(n1*n2)
      let (e0, q0) = hcEQ(wt, g)
      discard e0
      var q1, q2: float
      var nstep = 0
      g.hcGaugeFlow(0.02, hcFlowCflow):
        inc nstep
        if nstep == 50:
          q1 = hcEQ(wt, g).q
        elif nstep == 100:
          q2 = hcEQ(wt, g).q
          break
      echo &"  {icfg:3d}  ({n1:2d},{n2:2d})  {qex:7.2f}  {q0:9.4f}  {q1:9.4f}  {q2:9.4f}"
      maxdev = max(maxdev, abs(q2 - qex))
      maxplat = max(maxplat, abs(q2 - q1))
    ok(&"flowed Q returns to the exact sector integer " &
       &"(max |Q - Q_exact| = {maxdev:.4f} < 0.05)", maxdev < 0.05)
    ok(&"Q plateaus in t (max |Q(2)-Q(1)| = {maxplat:.4f} < 0.02)",
       maxplat < 0.02)

  test "6. integer-Q clustering on flowed rough configs (secondary, loose)":
    # Generic rough (warm 0.85) configs on 8^4, flowed deep.  Their lumps sit
    # at the cutoff scale, so -- exactly as task C found on the cubic side --
    # Q keeps drifting through dislocations and the clustering is LOOSE; the
    # sharp normalisation statements are tests 3 and 5.  Here we check the
    # distribution favours integers over half-odd-integers (the factor-1/2
    # failure mode) and report the histogram.
    let hl = newHcLayout([8, 8, 8, 8])
    const ncfg = 12
    var r = hl.lo.newRNGField(RngMilc6, 1122334455'u64)
    var g = newHcGauge(hl)
    var wt = newHcTopoWork(g)
    var qend, qmid: seq[float]
    echo "  cfg   Q(t=3)    Q(t=7)    Q(t=10)   dist(Q(10),int)"
    for icfg in 0..<ncfg:
      threads:
        g.warm(0.85, r)
      var q1, q2, q3: float
      var nstep = 0
      g.hcGaugeFlow(0.02, hcFlowCflow):
        inc nstep
        if nstep == 150:
          q1 = hcEQ(wt, g).q
        elif nstep == 350:
          q2 = hcEQ(wt, g).q
        elif nstep == 500:
          q3 = hcEQ(wt, g).q
          break
      qmid.add q2
      qend.add q3
      echo &"  {icfg:3d}  {q1:8.4f}  {q2:8.4f}  {q3:8.4f}   {abs(q3-round(q3)):.4f}"
    var hist: array[-8..8, int]
    var dInt = 0.0
    var dHalfOdd = 0.0
    var nnz = 0
    var nstable = 0
    for i in 0..<ncfg:
      let k = int(round(qend[i]))
      if k >= -8 and k <= 8: inc hist[k]
      let di = abs(qend[i]-round(qend[i]))
      dInt += di
      dHalfOdd += 0.5 - di
      if abs(qend[i]) > 0.5: inc nnz
      if abs(qend[i]-qmid[i]) < 0.25: inc nstable
    dInt = dInt/ncfg.float
    dHalfOdd = dHalfOdd/ncfg.float
    var hs = ""
    for k in -8..8:
      if hist[k] > 0: hs &= &"  {k}:{hist[k]}"
    echo &"  histogram of round(Q(t=10)):{hs}"
    echo &"  mean |Q - nearest integer|      = {dInt:.4f}   (uniform: 0.25)"
    echo &"  mean |Q - nearest half-odd int| = {dHalfOdd:.4f}"
    echo &"  {nnz}/{ncfg} configs with |Q| > 0.5;  ",
         &"{nstable}/{ncfg} with |Q(10)-Q(7)| < 0.25"
    echo "  (looser than tests 3/5 by construction -- cutoff-scale lumps keep"
    echo "   dislocating, cf. task C's identical finding on the cubic side)"
    ok(&"nonzero sectors sampled ({nnz}/{ncfg})", nnz > 0)
    ok(&"Q clusters at integers (mean dist {dInt:.3f} < 0.2; uniform 0.25)",
       dInt < 0.2)
    ok(&"distribution favours integers over half-odd-integers " &
       &"({dInt:.3f} < {dHalfOdd:.3f})", dInt < dHalfOdd)

qexFinalize()
