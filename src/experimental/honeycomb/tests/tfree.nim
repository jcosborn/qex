## Tests for the free momentum-space Wilson-Dirac operator (task F).
##
##   export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
##   cd build_mac && make run src/experimental/honeycomb/tests/tfree.nim
##
## Checks every number in FORMULATION section 5.3 plus the internal
## consistency of the three independent representations of the operator
## (direct neighbour sum, closed form, 8x8 cell-momentum matrix).

import std/[math, complex, strformat, unittest]
import ../hcfree

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

proc pr(s: string) = echo "    ", s

# ---------------------------------------------------------------------------

suite "gamma matrices (DeGrand-Rossi)":

  test "Clifford algebra and gamma5":
    var worst = 0.0
    for mu in 0..3:
      for nu in 0..3:
        let ac = mmul(gammaMat[mu], gammaMat[nu])
        let ca = mmul(gammaMat[nu], gammaMat[mu])
        for i in 0..3:
          for j in 0..3:
            var want = complex64(0.0)
            if mu == nu and i == j: want = complex64(2.0)
            worst = max(worst, abs(ac[i][j] + ca[i][j] - want))
      # hermiticity
      worst = max(worst, maxAbsDiff(adj(gammaMat[mu]), gammaMat[mu]))
    pr &"max |{{g_mu,g_nu}} - 2 delta| and |g_mu^dag - g_mu| = {worst:.3e}"
    check worst < 1e-14

    let g5 = mmul(mmul(gammaMat[0], gammaMat[1]), mmul(gammaMat[2], gammaMat[3]))
    let d5 = maxAbsDiff(g5, gamma5Mat)
    pr &"|g0 g1 g2 g3 - gamma5| = {d5:.3e}"
    check d5 < 1e-14
    var wa = 0.0
    for mu in 0..3:
      let a = mmul(mmul(gamma5Mat, gammaMat[mu]), gamma5Mat)
      for i in 0..3:
        for j in 0..3: wa = max(wa, abs(a[i][j] + gammaMat[mu][i][j]))
    pr &"max |g5 g_mu g5 + g_mu| = {wa:.3e}"
    check wa < 1e-14

# ---------------------------------------------------------------------------

suite "direct sum vs closed form":

  test "M(p) and K(p) agree for random p":
    var st = 1234567891'u64
    proc rnd(): float =
      st = st*6364136223846793005'u64 + 1442695040888963407'u64
      float((st shr 11).int64)*(1.0/9007199254740992.0)
    for lat in [lCubic, lHoneycomb]:
      var wm = 0.0
      var wk = 0.0
      for it in 0..<20000:
        var p: Vec4
        for mu in 0..3: p[mu] = 25.0*(rnd() - 0.5)
        let r = 0.3 + rnd()
        wm = max(wm, abs(freeM(p, lat, r) - freeMClosed(p, lat, r)))
        let
          k1 = freeK(p, lat)
          k2 = freeKClosed(p, lat)
        for mu in 0..3: wk = max(wk, abs(k1[mu] - k2[mu]))
      pr &"{lat}: max |M_sum - M_closed| = {wm:.3e}, max |K_sum - K_closed| = {wk:.3e}"
      check wm < 1e-12
      check wk < 1e-12

  test "sum_i n_mu n_nu = 6 delta (2 delta for cubic)":
    for lat in [lCubic, lHoneycomb]:
      let want = if lat == lCubic: 2.0 else: 6.0
      var worst = 0.0
      for mu in 0..3:
        for nu in 0..3:
          var s = 0.0
          for n in dirsOf(lat): s += n[mu]*n[nu]
          worst = max(worst, abs(s - (if mu == nu: want else: 0.0)))
      pr &"{lat}: max deviation from {want} delta_munu = {worst:.3e}"
      check worst < 1e-14

  test "matrix form: freeD4 == freeD4Direct":
    var st = 987654321'u64
    proc rnd(): float =
      st = st*6364136223846793005'u64 + 1442695040888963407'u64
      float((st shr 11).int64)*(1.0/9007199254740992.0)
    for lat in [lCubic, lHoneycomb]:
      var worst = 0.0
      for it in 0..<2000:
        var p: Vec4
        for mu in 0..3: p[mu] = 25.0*(rnd() - 0.5)
        let r = 0.3 + rnd()
        let m = rnd() - 0.5
        worst = max(worst, maxAbsDiff(freeD4(p, lat, r, m),
                                      freeD4Direct(p, lat, r, m)))
      pr &"{lat}: max |D(p) - sum_i (g.n_i - r) e^(ip.n)| = {worst:.3e}"
      check worst < 1e-12

  test "eigenvalues M +- i|K| annihilate det(D - lambda)":
    var st = 55555'u64
    proc rnd(): float =
      st = st*6364136223846793005'u64 + 1442695040888963407'u64
      float((st shr 11).int64)*(1.0/9007199254740992.0)
    for lat in [lCubic, lHoneycomb]:
      var worst = 0.0
      for it in 0..<500:
        var p: Vec4
        for mu in 0..3: p[mu] = 25.0*(rnd() - 0.5)
        let d = freeD4(p, lat, 1.0, 0.1)
        for e in freeEigs(p, lat, 1.0, 0.1):
          var a = d
          for i in 0..3: a[i][i] = a[i][i] - e
          worst = max(worst, abs(det(a)))
      pr &"{lat}: max |det(D - lambda)| = {worst:.3e}"
      check worst < 1e-11

# ---------------------------------------------------------------------------

suite "small-p limits (fixes the 1/6 and the Wilson normalisation)":

  test "|K|/|p| -> 1 and M/(r p^2/2) -> 1":
    let dirn: array[3, Vec4] = [[1.0, 0.0, 0.0, 0.0],
                                [1.0, 1.0, 1.0, 1.0],
                                [0.7, -0.3, 1.3, 0.2]]
    for lat in [lCubic, lHoneycomb]:
      for d in dirn:
        var nrm = 0.0
        for mu in 0..3: nrm += d[mu]*d[mu]
        nrm = sqrt(nrm)
        for eps in [1e-2, 1e-3, 1e-4]:
          var p: Vec4
          for mu in 0..3: p[mu] = eps*d[mu]/nrm
          let
            rk = absK(p, lat)/eps
            rm = freeMClosed(p, lat, 1.0)/(0.5*eps*eps)
          if eps == 1e-4:
            pr &"{lat} dir {d}: |K|/|p| = {rk:.10f}, M/(r p^2/2) = {rm:.10f}"
          check abs(rk - 1.0) < 20.0*eps*eps
          check abs(rm - 1.0) < 20.0*eps*eps

# ---------------------------------------------------------------------------

proc scanMax(lat: HcLat; f: proc (p: Vec4): float; ng: int):
    tuple[val: float, p: Vec4] =
  ## coarse grid over a full period box, then a shrinking local search
  let per = if lat == lCubic: 2.0*PI else: 4.0*PI
  result.val = -1e300
  var p: Vec4
  for n0 in 0..<ng:
    p[0] = per*n0.float/ng.float
    for n1 in 0..<ng:
      p[1] = per*n1.float/ng.float
      for n2 in 0..<ng:
        p[2] = per*n2.float/ng.float
        for n3 in 0..<ng:
          p[3] = per*n3.float/ng.float
          let v = f(p)
          if v > result.val:
            result.val = v
            result.p = p
  var h = per/ng.float
  while h > 1e-13:
    var improved = true
    while improved:
      improved = false
      for mu in 0..3:
        for s in [-1.0, 1.0]:
          var q = result.p
          q[mu] += s*h
          let v = f(q)
          if v > result.val:
            result.val = v
            result.p = q
            improved = true
    h *= 0.5

suite "extreme eigenvalues (slide 14 / FORMULATION 5.3)":

  test "max Re lambda":
    for lat in [lCubic, lHoneycomb]:
      let want = if lat == lCubic: 8.0 else: 16.0/3.0
      let (v, p) = scanMax(lat, proc (p: Vec4): float = freeMClosed(p, lat, 1.0), 16)
      pr &"{lat}: max Re lambda = {v:.12f} (want {want:.12f}) at p = " &
         &"({p[0]:.5f},{p[1]:.5f},{p[2]:.5f},{p[3]:.5f})"
      check abs(v - want) < 1e-11

  test "max |Im lambda|":
    for lat in [lCubic, lHoneycomb]:
      let want = if lat == lCubic: 2.0 else: 1.4678898250138706
      let (v, p) = scanMax(lat, proc (p: Vec4): float = absK(p, lat), 16)
      pr &"{lat}: max |Im lambda| = {v:.12f} (want {want:.12f}) at p = " &
         &"({p[0]:.5f},{p[1]:.5f},{p[2]:.5f},{p[3]:.5f})"
      check abs(v - want) < 1e-9

  test "16-cell |Im lambda| at the two quoted maxima":
    let x = 2.0*arccos(0.5*(sqrt(3.0) - 1.0))
    let a = absK([x, 0.0, 0.0, 0.0], lHoneycomb)
    let b = absK([0.5*x, 0.5*x, 0.5*x, 0.5*x], lHoneycomb)
    let c = absK([2.0*PI/3.0, 0.0, 0.0, 0.0], lHoneycomb)
    let exact = pow(3.0, 0.25)*(1.0 + sqrt(3.0))/sqrt(6.0)
    pr &"x = 2 arccos((sqrt3-1)/2) = {x:.12f} (x/2 = {0.5*x:.6f})"
    pr &"|K| at (x,0,0,0)         = {a:.12f}"
    pr &"|K| at (x/2,x/2,x/2,x/2) = {b:.12f}"
    pr &"closed form 3^(1/4)(1+sqrt3)/sqrt6 = {exact:.12f}"
    pr &"|K| at (2pi/3,0,0,0)     = {c:.12f}   <- NOT a maximum"
    check abs(a - 1.4678898250138706) < 1e-13
    check abs(b - 1.4678898250138706) < 1e-13
    check abs(exact - 1.4678898250138706) < 1e-13
    check abs(c - 1.4433756729740645) < 1e-13

# ---------------------------------------------------------------------------

suite "16-cell: 8x8 cell-momentum form":

  test "gamma5 D gamma5 = D^dagger":
    var st = 24680'u64
    proc rnd(): float =
      st = st*6364136223846793005'u64 + 1442695040888963407'u64
      float((st shr 11).int64)*(1.0/9007199254740992.0)
    let g5 = gamma5x8()
    var worst = 0.0
    for it in 0..<500:
      var k: Vec4
      for mu in 0..3: k[mu] = 2.0*PI*rnd()
      let
        d = freeD8(k, 0.8, 0.2)
        l = mmul(mmul(g5, d), g5)
      worst = max(worst, maxAbsDiff(l, adj(d)))
    pr &"max |g5 D g5 - D^dag| = {worst:.3e}"
    check worst < 1e-14

  test "8x8 spectrum = union of the two true-momentum branches":
    var st = 13579'u64
    proc rnd(): float =
      st = st*6364136223846793005'u64 + 1442695040888963407'u64
      float((st shr 11).int64)*(1.0/9007199254740992.0)
    var worstTr = 0.0
    var worstDet = 0.0
    for it in 0..<300:
      var k: Vec4
      for mu in 0..3: k[mu] = 2.0*PI*rnd()
      let
        r = 0.5 + rnd()
        m = rnd() - 0.5
        d = freeD8(k, r, m)
      # predicted multiset: 4 eigenvalues per branch, each doubly degenerate
      var lam: seq[Complex64]
      for p in branchMomenta(k):
        for e in freeEigs(p, lHoneycomb, r, m):
          lam.add e
          lam.add e
      doAssert lam.len == 8
      # Newton's identities: matching tr(D^n), n = 1..8, fixes the multiset
      var a = d
      for n in 1..8:
        if n > 1: a = mmul(a, d)
        var s = complex64(0.0)
        for e in lam:
          var t = complex64(1.0)
          for j in 1..n: t = t*e
          s = s + t
        worstTr = max(worstTr, abs(tr(a) - s)/max(1.0, abs(s)))
      for e in lam:
        var b = d
        for i in 0..7: b[i][i] = b[i][i] - e
        worstDet = max(worstDet, abs(det(b)))
    pr &"max rel |tr D^n - sum lambda^n|, n=1..8 : {worstTr:.3e}"
    pr &"max |det(D8 - lambda)|                  : {worstDet:.3e}"
    check worstTr < 1e-11
    check worstDet < 1e-9

# ---------------------------------------------------------------------------

suite "Brillouin zone bookkeeping":

  test "momentum count = number of sites":
    for (ns, nt) in [(4, 4), (4, 6), (3, 5)]:
      var nc = 0
      for p in momenta(lCubic, ns, nt): inc nc
      var nh = 0
      for p in momenta(lHoneycomb, ns, nt): inc nh
      var nh2 = 0
      for p in momentaAlt(lHoneycomb, ns, nt): inc nh2
      pr &"ns={ns} nt={nt}: cubic {nc} (= {ns*ns*ns*nt}), " &
         &"16-cell {nh} / {nh2} (= {2*ns*ns*ns*nt})"
      check nc == ns*ns*ns*nt
      check nh == 2*ns*ns*ns*nt
      check nh2 == 2*ns*ns*ns*nt

  test "the two 16-cell momentum sets give the same spectral sums":
    for (ns, nt) in [(4, 4), (3, 5), (5, 4)]:
      var s1, s2: array[4, float]
      for p in momenta(lHoneycomb, ns, nt):
        let e = freeEigs(p, lHoneycomb, 1.0, 0.0)
        s1[0] += ln(abs2(e[0]))
        s1[1] += e[0].re
        s1[2] += e[0].re*e[0].re
        s1[3] += abs2(e[0])*abs2(e[0])
      for p in momentaAlt(lHoneycomb, ns, nt):
        let e = freeEigs(p, lHoneycomb, 1.0, 0.0)
        s2[0] += ln(abs2(e[0]))
        s2[1] += e[0].re
        s2[2] += e[0].re*e[0].re
        s2[3] += abs2(e[0])*abs2(e[0])
      var worst = 0.0
      for i in 0..3: worst = max(worst, abs(s1[i] - s2[i])/max(1.0, abs(s1[i])))
      pr &"ns={ns} nt={nt}: max rel difference of 4 spectral sums = {worst:.3e}"
      check worst < 1e-12

  test "massless spectrum has exactly one zero mode (p = 0)":
    for lat in [lCubic, lHoneycomb]:
      let ns = 6
      let nt = 6
      var nz = 0
      var minNonzero = 1e300
      for p in momenta(lat, ns, nt, antiperiodicTime = false):
        let a = abs2(freeEigs(p, lat)[0])
        if a < 1e-20: inc nz
        else: minNonzero = min(minNonzero, sqrt(a))
      pr &"{lat}: zero modes on {ns}^3x{nt} (periodic) = {nz}, " &
         &"smallest nonzero |lambda| = {minNonzero:.6f}"
      check nz == 1

# ---------------------------------------------------------------------------

suite "transfer-matrix modes (the free-energy kernel)":

  test "P(x) reproduces M^2+|K|^2":
    var st = 777'u64
    proc rnd(): float =
      st = st*6364136223846793005'u64 + 1442695040888963407'u64
      float((st shr 11).int64)*(1.0/9007199254740992.0)
    for lat in [lCubic, lHoneycomb]:
      var worst = 0.0
      for it in 0..<5000:
        var q: array[3, float]
        for j in 0..2: q[j] = 2.0*PI*rnd()
        let th = 2.0*PI*rnd()
        var p: Vec4
        for j in 0..2: p[j] = q[j]
        p[3] = if lat == lCubic: th else: 2.0*th
        let
          mm = freeMClosed(p, lat, 1.0)
          kk = freeKClosed(p, lat)
        var want = mm*mm
        for mu in 0..3: want += kk[mu]*kk[mu]
        worst = max(worst, abs(polyEval(polyB(q, lat), cos(th)) - want))
      pr &"{lat}: max |P(cos theta) - (M^2+|K|^2)| = {worst:.3e}"
      check worst < 1e-11

  test "Matsubara identity  (1/N')sum_n ln P - int ln P = (2/N') sum_i ln|1+v_i^N'|":
    var st = 31337'u64
    proc rnd(): float =
      st = st*6364136223846793005'u64 + 1442695040888963407'u64
      float((st shr 11).int64)*(1.0/9007199254740992.0)
    for lat in [lCubic, lHoneycomb]:
      var worst = 0.0
      for it in 0..<200:
        var q: array[3, float]
        for j in 0..2: q[j] = 0.05 + (2.0*PI - 0.1)*rnd()
        let b = polyB(q, lat)
        # reference integral: trapezoid in theta, spectrally accurate
        const nth = 4096
        var iref = 0.0
        for n in 0..<nth: iref += ln(polyEval(b, cos(2.0*PI*n.float/nth.float)))
        iref /= nth.float
        for nt in [3, 4, 6, 10]:
          let np = modeExponent(lat, nt)
          var s = 0.0
          for n in 0..<np: s += ln(polyEval(b, cos(2.0*PI*(n.float + 0.5)/np.float)))
          s /= np.float
          let rhs = (2.0/np.float)*pressureIntegrand(q, lat, nt)
          worst = max(worst, abs((s - iref) - rhs))
      pr &"{lat}: worst deviation = {worst:.3e}"
      check worst < 1e-12

  test "end-to-end: direct momentum sum reproduces the pressure formula":
    ## `p(T) = (1/(Nt Ns^3)) sum_p ln det D(p)` with `det = (M^2+|K|^2)^2`,
    ## against `p(T) - p(T') = (4/Nt)<h_Nt> - (4/Nt')<h_Nt'>` on the same
    ## spatial grid.  This fixes every factor of 2 in `O = Nt^3 <h>`.
    const ns = 6
    for lat in [lCubic, lHoneycomb]:
      proc fdirect(nt: int): float =
        for p in momentaAlt(lat, ns, nt):
          let
            mm = freeMClosed(p, lat, 1.0)
            kk = freeKClosed(p, lat)
          var a = mm*mm
          for mu in 0..3: a += kk[mu]*kk[mu]
          result += 2.0*ln(a)
        result /= nt.float*(ns.float^3)
      proc hav(nt: int): float =
        var q: array[3, float]
        for n0 in 0..<ns:
          q[0] = 2.0*PI*n0.float/ns.float
          for n1 in 0..<ns:
            q[1] = 2.0*PI*n1.float/ns.float
            for n2 in 0..<ns:
              q[2] = 2.0*PI*n2.float/ns.float
              result += pressureIntegrand(q, lat, nt)
        result /= ns.float^3
      var worst = 0.0
      for (a, b) in [(3, 4), (4, 6), (5, 8)]:
        let
          lhs = fdirect(a) - fdirect(b)
          rhs = 4.0/a.float*hav(a) - 4.0/b.float*hav(b)
        worst = max(worst, abs(lhs - rhs))
      pr &"{lat}: max |p(T)-p(T') direct - formula| = {worst:.3e}"
      check worst < 1e-12

  test "pressure ratio on a fixed ng=64 grid (pinning test)":
    ## cheap regression values; the converged numbers come from
    ## hcFreePressure with Richardson extrapolation over ng = 48..1024
    const ng = 64
    for (lat, nt, want) in [(lCubic, 4, 3.82826097), (lHoneycomb, 4, 1.07204731)]:
      var s = 0.0
      var q: array[3, float]
      for n0 in 0..<ng:
        q[0] = 2.0*PI*n0.float/ng.float
        for n1 in 0..<ng:
          q[1] = 2.0*PI*n1.float/ng.float
          for n2 in 0..<ng:
            q[2] = 2.0*PI*n2.float/ng.float
            s += pressureIntegrand(q, lat, nt)
      let o = nt.float^3*s/(ng.float^3)/contPressure
      pr &"{lat} Nt={nt}: O/O_cont(ng={ng}) = {o:.8f} (expect {want:.8f})"
      check abs(o - want) < 1e-7

  test "light mode: 2E -> |q| as q -> 0 (16-cell), E -> |q| (cubic)":
    for lat in [lCubic, lHoneycomb]:
      for eps in [1e-2, 1e-3]:
        let q: array[3, float] = [eps*0.6, eps*0.8, 0.0]
        var best = 0.0
        for v in modeV(q, lat):
          best = max(best, abs(v))
        let
          e = -ln(best)
          omega = if lat == lCubic: e else: 2.0*e
        if eps == 1e-3:
          pr &"{lat}: omega/|q| at |q|={eps} is {omega/eps:.10f}"
        check abs(omega/eps - 1.0) < 20.0*eps*eps
