## stout smearing (hcflow.stout, cubic.smearN).
##
## 1. smeared links exactly SU(3) on both lattices;
## 2. gauge covariance smear(U^g) = (smear U)^g;
## 3. rho = 0 is the identity;
## 4. smoothing: triangleSum / plaq strictly increase over 6 steps at
##    rho = 0.05 on a warm configuration (the sequences are printed);
## 5. THE NORMALISATION PIN: for a weak Abelian plane wave (exact line
##    integrals, the taction/tflow recipe) n stout steps multiply the mode
##    amplitude by exp(-t_eff p^2), t_eff = n rho kappa.  kappa is measured
##    from the per-step decay of the (quadratic) gauge action, three momenta,
##    the O(p^2) artifact removed by quadratic extrapolation:
##      cubic     kappa = 1    (exact; per-step factor even pinned against the
##                              closed form (1 - rho phat^2) of the linearised
##                              Euler/stout map)
##      honeycomb kappa = 1/3  (exact given hcflow's cflow = 6: one stout step
##                              = one Euler flow step of eps = rho * 2/cflow)
##    plus a rho-independence check of the extrapolated honeycomb kappa.
##    ==> rho equivalence for task D4: rho_hc = 3 rho_cubic at equal smearing.

import std/[math, strformat, unittest]
import qex except epsilon
import physics/qcdTypes
import gauge, gauge/gaugefix
import ../hcgeom
import ../hcgauge
import ../hcaction
import ../hcflow
import ../cubic
import helpers

const cubicStoutKappa = 1.0

qexInit()

let seed = 246810121'u64

proc cubicNorm2Diff(a, b: auto): float =
  for mu in 0..<a.len:
    result += norm2diff(a[mu], b[mu])

proc cubicNorm2(a: auto): float =
  for mu in 0..<a.len:
    result += norm2(a[mu])

proc cubicCheckSU(g: auto): tuple[avg, max: float] =
  var res: tuple[avg, max: float]
  threads:
    let d = checkSU(g)
    threadMaster: res = d
  res

proc hcStoutRate(ns, k: int, rho: float, nsteps: int):
    tuple[rate, drift, p2: float] =
  ## per-step heat-kernel decay rate of the plane-wave mode under stout, in
  ## units of rho: rate = ln(S_{n-1}/S_n)/(2 rho) -> kappa p^2 (1 + O(p^2))
  var pv, epsv: array[4, float]
  pv[0] = 2.0*PI*float(k)/float(ns)
  for mu in 0..<4: epsv[mu] = epsAmp*ehat[mu]
  let hl = newLayout([ns, ns, ns, ns])
  var g = newHcGauge(hl)
  setAbelianHc(g, epsv, pv)
  var w = newActionWork(g)
  var sc = newOneOf(g)
  var ss = @[action(w, 1.0, g)]
  for n in 0..<nsteps:
    stout(w, g, rho, sc)
    ss.add action(w, 1.0, g)
  let (r, d) = plateauRate(ss, rho)
  (r, d, pv[0]*pv[0])

proc cubicStoutRate(ns, k: int, rho: float, nsteps: int):
    tuple[rate, drift, p2, perstep: float] =
  ## same on the cubic lattice with QEX StoutSmear.
  ## `perstep` = the late-window per-step ln(S ratio)/2, to compare with the
  ## exact linearised Euler/stout map factor -ln(1 - rho phat^2).
  var pv, epsv: array[4, float]
  pv[0] = 2.0*PI*float(k)/float(ns)
  for mu in 0..<4: epsv[mu] = epsAmp*ehat[mu]
  let lo = newLayout(@[ns, ns, ns, ns])
  var g = lo.newGauge
  setAbelianCubic(g, lo, epsv, pv)
  var st = newStoutSmear(lo, rho)
  let gact = GaugeActionCoeffs(plaq: 1.0)
  var ss = @[wilsonAction(gact, g)]
  for n in 0..<nsteps:
    st.smear(g, g)
    ss.add wilsonAction(gact, g)
  let (r, d) = plateauRate(ss, rho)
  (r, d, pv[0]*pv[0], r*rho)

suite "hcstout (task D2 part 1)":

  test "1. smeared links exactly SU(3), both lattices":
    block:                      # honeycomb
      let hl = newLayout([4, 4, 4, 6])
      var r = hl.newRNGField(RngMilc6, seed)
      var g = newHcGauge(hl)
      threads:
        g.warm(0.5, r)
      var w = newActionWork(g)
      var sc = newOneOf(g)
      stout(w, g, 0.05, sc, 6)
      let d = g.checkSU
      echo &"  honeycomb [4,4,4,6], 6 steps rho=0.05: checkSU avg {d.avg:.2e} max {d.max:.2e}"
      ok(&"honeycomb smeared links SU(3) to 1e-10 (max {d.max:.2e})",
         d.max < 1e-10)
    block:                      # cubic
      let lo = newLayout(@[8, 8, 8, 8])
      var r = lo.newRNGField(RngMilc6, seed + 1)
      var g = lo.newGauge
      threads:
        g.random r
      var st = newStoutSmear(lo, 0.05)
      st.smearN(g, g, 6)
      let d = cubicCheckSU(g)
      echo &"  cubic 8^4 random, 6 steps rho=0.05: checkSU avg {d.avg:.2e} max {d.max:.2e}"
      ok(&"cubic smeared links SU(3) to 1e-10 (max {d.max:.2e})", d.max < 1e-10)

  test "2. gauge covariance: smear(U^g) = (smear U)^g":
    block:                      # honeycomb
      let hl = newLayout([4, 4, 4, 6])
      let lo = hl
      var r = lo.newRNGField(RngMilc6, seed + 2)
      var g = newHcGauge(hl)
      var g2 = newHcGauge(hl)
      var s1 = newHcGauge(hl)
      threads:
        g.random r
        g2 := g
      var vA = lo.ColorMatrix(nc)
      var vB = lo.ColorMatrix(nc)
      threads:
        vA.randomSU r
        vB.randomSU r
      var w = newActionWork(g)
      var sc = newOneOf(g)
      threads:
        s1 := g
      stout(w, s1, 0.05, sc)
      s1.gaugeTransform(vA, vB) # s1 = (smear U)^g
      g2.gaugeTransform(vA, vB)
      stout(w, g2, 0.05, sc)
      let rel = sqrt(norm2diff(g2, s1)/norm2(s1))
      echo &"  honeycomb |smear(U^g) - (smear U)^g| / |.| = {rel:.3e}"
      ok(&"honeycomb stout gauge covariant to 1e-12 ({rel:.3e})", rel < 1e-12)
    block:                      # cubic
      let lo = newLayout(@[4, 4, 4, 8])
      var r = lo.newRNGField(RngMilc6, seed + 3)
      var g = lo.newGauge
      var g2 = lo.newGauge
      var s1 = lo.newGauge
      var t = lo.ColorMatrix(nc)
      threads:
        g.random r
        t.randomSU r
      var st = newStoutSmear(lo, 0.05)
      st.smear(g, s1)
      gaugeTransform(s1, s1, t)     # s1 <- t s1 t(x+mu)^dag
      gaugeTransform(g2, g, t)
      st.smear(g2, g2)
      let rel = sqrt(cubicNorm2Diff(g2, s1)/cubicNorm2(s1))
      echo &"  cubic |smear(U^g) - (smear U)^g| / |.| = {rel:.3e}"
      ok(&"cubic stout gauge covariant to 1e-12 ({rel:.3e})", rel < 1e-12)

  test "3. rho = 0 is the identity":
    block:
      let hl = newLayout([4, 4, 4, 4])
      var r = hl.newRNGField(RngMilc6, seed + 4)
      var g = newHcGauge(hl)
      var g1 = newHcGauge(hl)
      threads:
        g.random r
      var w = newActionWork(g)
      var sc = newOneOf(g)
      threads:
        g1 := g
      stout(w, g1, 0.0, sc)
      let d = norm2diff(g1, g)
      echo &"  honeycomb rho=0: |smear(U) - U|^2 = {d:.3e}"
      ok("honeycomb rho=0 exact identity", d == 0.0)
    block:
      let lo = newLayout(@[4, 4, 4, 4])
      var r = lo.newRNGField(RngMilc6, seed + 5)
      var g = lo.newGauge
      var g1 = lo.newGauge
      threads:
        g.random r
      var st = newStoutSmear(lo, 0.0)
      st.smear(g, g1)
      let d = cubicNorm2Diff(g1, g)
      echo &"  cubic rho=0: |smear(U) - U|^2 = {d:.3e}"
      ok("cubic rho=0 exact identity", d == 0.0)

  test "cubic alpha changes the next smear":
    let lo = newLayout(@[4, 4, 4, 4])
    var r = lo.newRNGField(RngMilc6, seed + 10)
    var g = lo.newGauge
    var gs = lo.newGauge
    threads:
      g.random r
    var st = newStoutSmear(lo, 0.05)
    st.smear(g, gs)
    check cubicNorm2Diff(gs, g) > 0.0
    st.alpha = 0.0
    st.smear(g, gs)
    check cubicNorm2Diff(gs, g) == 0.0

  test "cubic repeated smearing copies or repeats the same step":
    let lo = newLayout(@[4, 4, 4, 4])
    var r = lo.newRNGField(RngMilc6, seed + 11)
    var g = lo.newGauge
    var gs = lo.newGauge
    var gt = lo.newGauge
    threads:
      g.random r
      gs.unit
    var st = newStoutSmear(lo, 0.05)
    st.smearN(g, gs, 0)
    check cubicNorm2Diff(gs, g) == 0.0
    st.smearN(gs, gs, 0)
    check cubicNorm2Diff(gs, g) == 0.0
    st.smearN(g, gs, 3)
    st.smear(g, gt)
    st.smear(gt, gt)
    st.smear(gt, gt)
    check cubicNorm2Diff(gs, gt) == 0.0
    st.smearN(g, gt, 0)
    st.smearN(gt, gt, 3)
    check cubicNorm2Diff(gs, gt) == 0.0

  test "4. smoothing: 6 steps at rho=0.05 strictly increase the loop sums":
    block:                      # honeycomb triangleSum
      let hl = newLayout([8, 8, 8, 8])
      var r = hl.newRNGField(RngMilc6, seed + 6)
      var g = newHcGauge(hl)
      threads:
        g.warm(0.35, r)
      var w = newActionWork(g)
      var sc = newOneOf(g)
      var ts = @[g.triangleSum]
      for n in 0..<6:
        stout(w, g, 0.05, sc)
        ts.add g.triangleSum
      var s = &"{ts[0]:.6f}"
      for i in 1..<ts.len: s &= &" -> {ts[i]:.6f}"
      echo "  honeycomb triangleSum (warm 8^4): ", s
      var mono = true
      for i in 1..<ts.len:
        if ts[i] <= ts[i-1]: mono = false
      ok("honeycomb triangleSum strictly increases over 6 steps", mono)
    block:                      # cubic plaq
      let lo = newLayout(@[8, 8, 8, 8])
      var r = lo.newRNGField(RngMilc6, seed + 7)
      var g = lo.newGauge
      threads:
        g.warm(0.35, r)
      var st = newStoutSmear(lo, 0.05)
      var ps = @[g.plaq.sum]
      for n in 0..<6:
        st.smear(g, g)
        ps.add g.plaq.sum
      var s = &"{ps[0]:.6f}"
      for i in 1..<ps.len: s &= &" -> {ps[i]:.6f}"
      echo "  cubic plaq (warm 8^4): ", s
      var mono = true
      for i in 1..<ps.len:
        if ps[i] <= ps[i-1]: mono = false
      ok("cubic plaq strictly increases over 6 steps", mono)

  test "5. heat-kernel normalisation: measure kappa on both lattices":
    const rho = 0.05
    # --- cubic: kappa = 1 (textbook: stout step = Euler Wilson-flow step) ---
    echo &"  cubic, rho = {rho}:"
    echo "  Ns k   p^2       rate        rate/p^2   perstep/exact  drift"
    var xs, ys: array[3, float]
    let cases = [(12, 1, 60), (8, 1, 60), (12, 2, 60)]
    var worstStep = 0.0
    for i in 0..<3:
      let (ns, k, nst) = cases[i]
      let (r, d, p2, perstep) = cubicStoutRate(ns, k, rho, nst)
      let p = sqrt(p2)
      let ph2 = 4.0*sin(0.5*p)^2
      # exact linearised Euler/stout map: amplitude factor (1 - rho phat^2)
      let stepExact = -ln(1.0 - rho*ph2)
      xs[i] = p2
      ys[i] = r/p2
      worstStep = max(worstStep, abs(perstep/stepExact - 1.0))
      echo &"  {ns:2d} {k:2d}  {p2:.5f}  {r:.8f}  {r/p2:.6f}   {perstep/stepExact:.8f}   {d:.2e}"
    ok(&"cubic per-step factor = -ln(1 - rho phat^2) to 1e-3 " &
       &"(worst dev {worstStep:.2e})", worstStep < 1e-3)
    let cc = fit3(xs, ys)
    echo &"  cubic kappa = rate/p^2 (p->0): {cc[0]:.6f}   (O(p^2) coeff {cc[1]:.4f})"
    ok(&"cubic kappa = 1 (measured {cc[0]:.6f}, |dev| {abs(cc[0]-1.0):.1e} < 3e-3)",
       abs(cc[0] - 1.0) < 3e-3)
    # --- honeycomb: kappa = 1/3 (exact given hcflow's cflow = 6) ---
    # 80 steps: t_eff = 1.33, the last-quarter window starts at t_eff = 1
    # where the optical transients (rate ~ cflow) are gone
    echo &"  honeycomb, rho = {rho}:"
    echo "  Ns k   p^2       rate        rate/p^2    drift"
    var xh, yh: array[3, float]
    let hcases = [(12, 1, 80), (8, 1, 80), (12, 2, 80)]
    for i in 0..<3:
      let (ns, k, nst) = hcases[i]
      let (r, d, p2) = hcStoutRate(ns, k, rho, nst)
      xh[i] = p2
      yh[i] = r/p2
      echo &"  {ns:2d} {k:2d}  {p2:.5f}  {r:.8f}  {r/p2:.7f}   {d:.2e}"
    let ch = fit3(xh, yh)
    let kappaHc = ch[0]
    echo &"  honeycomb kappa = rate/p^2 (p->0): {kappaHc:.7f}   (O(p^2) coeff {ch[1]:.4f})"
    # simple-rational scan (tflow recipe): smallest denominator within 1e-3
    var bestN, bestD = 0
    var bestErr = 1.0
    for den in 1..12:
      let num = int(round(kappaHc*den.float))
      if num < 1: continue
      let e = abs(kappaHc - num.float/den.float)
      if e < bestErr:
        bestErr = e
        bestN = num
        bestD = den
      if e < 1e-3: break
    echo &"  closest simple rational: {bestN}/{bestD} (|dev| = {bestErr:.1e})"
    ok(&"honeycomb kappa is the SIMPLE RATIONAL 1/3 (measured {kappaHc:.6f}, " &
       &"|dev| {abs(kappaHc - 1.0/3.0):.1e} < 1e-3)",
       bestN == 1 and bestD == 3 and abs(kappaHc - 1.0/3.0) < 1e-3)
    ok(&"shipped stoutKappa = {stoutKappa:.6f} agrees",
       abs(kappaHc - stoutKappa) < 1e-3)
    ok(&"shipped cubicStoutKappa = {cubicStoutKappa:.1f} agrees",
       abs(cc[0] - cubicStoutKappa) < 3e-3)
    # rho-independence of the extrapolated honeycomb kappa (Euler-step error
    # is O(rho p^2), i.e. it only feeds the artifact slope, not the intercept)
    block:
      let (rA, _, p2A) = hcStoutRate(12, 1, 0.5*rho, 160)
      let (rB, _, p2B) = hcStoutRate(8, 1, 0.5*rho, 160)
      let kHalf = rA/p2A - p2A*((rB/p2B - rA/p2A)/(p2B - p2A))
      echo &"  honeycomb kappa at rho = {0.5*rho} (2-pt extrap): {kHalf:.7f}"
      ok(&"kappa rho-independent (|kappa(rho/2) - 1/3| = " &
         &"{abs(kHalf - 1.0/3.0):.1e} < 1e-3)", abs(kHalf - 1.0/3.0) < 1e-3)
    echo ""
    echo "  ==> t_eff per step: cubic rho, honeycomb rho/3."
    echo "      Equal smearing radius sqrt(8 t_eff) needs rho_hc = 3 rho_cubic;"
    echo &"      the paper's 6 steps at rho = 0.05: t_eff = {6*0.05*cubicStoutKappa:.3f} a^2 (cubic)"
    echo &"      vs {6*0.05*stoutKappa:.3f} a^2 (honeycomb) in the plain MP staple-sum convention."

qexFinalize()
