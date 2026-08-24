## Task A test: `hcaction.nim` (triangle action, staple derivative, force).
##
##   export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
##   cd build_mac && make src/experimental/honeycomb/tests/taction.nim && \
##     OMP_NUM_THREADS=4 ./bin/taction
##
## The make-or-break test is the finite-difference force check: for random
## traceless anti-Hermitian momenta P on all 24 link fields,
##   d/ds S(exp(s P) U)|_0  must equal  sum_l redot(P_l, f_l)
## to ~1e-7 or better, on a warm (non-unit) configuration, for the full
## momentum and restricted to each link kind separately.  The same identity is
## first verified for QEX's own `gaugeForce` on a cubic lattice, so the force
## convention is pinned to QEX's, not assumed.

import math, strformat, unittest
import qex except epsilon
import physics/qcdTypes
import algorithms/numdiff
import ../hcgeom
import ../hclayout
import ../hcgauge
import ../hcaction

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

template ok(msg: string, cond: bool) =
  let c = cond
  if c: echo "PASS: ", msg
  else: echo "FAIL: ", msg
  check c

qexInit()

const nc = getDefaultNc()
let
  beta = 5.7
  seed = 135792468'u64

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

proc applyExpHc(gt, g0, p: auto, t: float) =
  ## gt := exp(t*p) * g0 on all 24 link fields (the HMC `mdt` update)
  threads:
    for mu in 0..<nDim:
      for e in gt.uA[mu]:
        gt.uA[mu][e] := exp(t*p.uA[mu][e])*g0.uA[mu][e]
      for e in gt.uB[mu]:
        gt.uB[mu][e] := exp(t*p.uB[mu][e])*g0.uB[mu][e]
    for d in 0..<nDiag:
      for e in gt.uD[d]:
        gt.uD[d][e] := exp(t*p.uD[d][e])*g0.uD[d][e]

proc setMom(p: auto, r: auto, sel: int) =
  ## random TAH momenta on a subset of the 24 link fields:
  ## 0 = all, 1 = uA only, 2 = uB only, 3 = uD only, 4 = uD[5] only
  threads:
    eachLink(p, u):
      u := 0
    case sel
    of 0:
      eachLink(p, u):
        u.randomTAH r
    of 1:
      for mu in 0..<nDim: p.uA[mu].randomTAH r
    of 2:
      for mu in 0..<nDim: p.uB[mu].randomTAH r
    of 3:
      for d in 0..<nDiag: p.uD[d].randomTAH r
    else:
      p.uD[5].randomTAH r

proc diffNorm2(a, b, scratch: auto): float =
  ## sum_l |a_l - b_l|^2 over all 24 link fields
  threads:
    for mu in 0..<nDim:
      scratch.uA[mu] := a.uA[mu] - b.uA[mu]
      scratch.uB[mu] := a.uB[mu] - b.uB[mu]
    for d in 0..<nDiag:
      scratch.uD[d] := a.uD[d] - b.uD[d]
  redot(scratch, scratch)

proc fdForce(geom: seq[int], seed: uint64, sel: int, label: string,
             V: static int = VLEN): tuple[rel, nde: float] =
  ## finite-difference check of hcForce on a warm configuration
  let hl = newHcLayoutX(geom, V)
  let lo = hl.lo
  var r = lo.newRNGField(RngMilc6, seed)
  var g0 = newHcGauge(hl)
  var gt = newOneOf(g0)
  var p = newOneOf(g0)
  var f = newOneOf(g0)
  threads:
    g0.warm(0.35, r)
  setMom(p, r, sel)
  var w = newHcActionWork(g0)
  hcForce(w, beta, g0, f)
  let pf = redot(p, f)
  proc sAt(t: float): float =
    applyExpHc(gt, g0, p, t)
    hcAction(w, beta, gt)
  var dsdt, e: float
  ndiff(dsdt, e, sAt, 0.0, 0.1, ordMax = 7)
  let rel = abs(dsdt - pf)/max(1e-300, max(abs(dsdt), abs(pf)))
  echo &"  {label:11s}: dS/dt = {dsdt: .13e}   sum redot(p,f) = {pf: .13e}"
  echo &"  {label:11s}: rel diff = {rel:.3e}  (ndiff err estimate {e:.1e})"
  (rel, e)

# --- weak Abelian plane wave, exact straight-line integrals -----------------

proc phaseIntegral(a, b: float): float =
  ## int_0^1 cos(a + t b) dt
  if abs(b) < 1e-12: cos(a)
  else: (sin(a+b) - sin(a))/b

proc linkPhase(x, n, eps, pv: array[4, float]): float =
  ## int_x^{x+n} A.dl  for  A_mu(x) = eps_mu cos(p.x)
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
  ## honeycomb links U = exp(i phi T) from exact line integrals of A
  hg.unit
  let lo = hg.lo
  for i in lo.sites:
    var y, yb: array[4, float]
    for mu in 0..<4:
      y[mu] = lo.coords[mu][i].float    # A(y) at integer coords
      yb[mu] = y[mu] + 0.5              # B(y) at half-integer coords
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

proc setAbelianCubic(g: auto, lo: auto, eps, pv: array[4, float]) =
  ## cubic links from the same continuum field
  g.unit
  for i in lo.sites:
    var y: array[4, float]
    for mu in 0..<4:
      y[mu] = lo.coords[mu][i].float
    for mu in 0..<4:
      var n: array[4, float]
      n[mu] = 1.0
      setPhase(g[mu]{i}, linkPhase(y, n, eps, pv))

proc wilsonAction(gact: GaugeActionCoeffs, g: auto): float =
  ## S = beta sum_{x,mu<nu} (1 - Re Tr P/N); gaugeAction1 omits the constant
  gact.gaugeAction1(g) + gact.plaq*6.0*float(g[0].l.physVol)

# ---------------------------------------------------------------------------

suite "hcaction":

  test "1. QEX gaugeForce satisfies the redot convention (cubic reference)":
    # d/ds S(exp(sP)U)|_0 = sum_mu redot(P_mu, F_mu) for QEX's own Wilson
    # action and force -- this is the convention hcForce must reproduce.
    let lo = newLayout(@[4, 4, 4, 4])
    var r = lo.newRNGField(RngMilc6, seed + 11)
    var g = lo.newGauge
    var gt = lo.newGauge
    var p = lo.newGauge
    var f = lo.newGauge
    threads:
      g.random r
      p.randomTAH r
    let gact = GaugeActionCoeffs(plaq: beta)
    gact.gaugeForce(g, f)
    var pf = 0.0
    threads:
      var s = 0.0
      for mu in 0..<4:
        s += redot(p[mu], f[mu])
      threadMaster: pf = s
    proc sAt(t: float): float =
      threads:
        for mu in 0..<4:
          for e in gt[mu]:
            gt[mu][e] := exp(t*p[mu][e])*g[mu][e]
      gact.gaugeAction1(gt)
    var dsdt, e: float
    ndiff(dsdt, e, sAt, 0.0, 0.1, ordMax = 7)
    let rel = abs(dsdt - pf)/max(abs(dsdt), abs(pf))
    echo &"  dS/dt = {dsdt: .13e}   sum redot(p,f) = {pf: .13e}"
    echo &"  rel diff = {rel:.3e}  (ndiff err estimate {e:.1e})"
    ok(&"QEX gaugeForce convention: rel diff {rel:.2e} < 1e-8", rel < 1e-8)

  test "2. hcAction: unit gauge, triangleSum consistency, sanity":
    let hl = newHcLayout([4, 4, 4, 6])
    let lo = hl.lo
    var r = lo.newRNGField(RngMilc6, seed + 22)
    var g = newHcGauge(hl)
    var w = newHcActionWork(g)
    let sTri = float(nTriPerSite*hl.nSites)
    # unit gauge: S = 0, force = 0
    let s0 = hcAction(w, beta, g)
    echo &"  hcAction(unit) = {s0:.3e}  (scale beta/2*nTri = {0.5*beta*sTri:.4g})"
    ok(&"unit gauge S = 0 ({abs(s0):.2e})", abs(s0) < 1e-8)
    var f = newOneOf(g)
    hcForce(w, beta, g, f)
    let f0n = redot(f, f)
    echo &"  |force(unit)|^2 = {f0n:.3e}"
    ok(&"unit gauge force = 0 ({f0n:.2e})", f0n < 1e-24)
    # warm configuration: independent cross-check against Task L's
    # triangleSum (verified in tgauge.nim against a brute-force reference)
    threads:
      g.warm(0.35, r)
    let sw = hcAction(w, beta, g)
    let swRef = 0.5*beta*sTri*(1.0 - g.triangleSum)
    let dw = abs(sw - swRef)/abs(swRef)
    echo &"  warm:   hcAction = {sw:.15g}"
    echo &"          (beta/2)*32*nSites*(1-triangleSum) = {swRef:.15g}"
    ok(&"warm action matches triangleSum path (rel {dw:.2e})", dw < 1e-10)
    # random configuration
    threads:
      g.random r
    let sr = hcAction(w, beta, g)
    let srRef = 0.5*beta*sTri*(1.0 - g.triangleSum)
    let dr = abs(sr - srRef)/abs(srRef)
    echo &"  random: hcAction = {sr:.15g}"
    echo &"          (beta/2)*32*nSites*(1-triangleSum) = {srRef:.15g}"
    ok(&"random action matches triangleSum path (rel {dr:.2e})", dr < 1e-10)
    ok(&"hcAction > 0 on random ({sr:.6g}) and warm ({sw:.6g})",
       sr > 0.0 and sw > 0.0)
    ok(&"warm action below random action ({sw:.4g} < {sr:.4g})", sw < sr)
    # convenience overloads agree with the work-based ones
    let sc = hcAction(beta, g)
    ok(&"hcAction convenience overload agrees ({abs(sc-sr):.2e})",
       abs(sc - sr) <= 1e-12*abs(sr))
    var f2 = newOneOf(g)
    var scr = newOneOf(g)
    hcForce(w, beta, g, f)
    hcForce(beta, g, f2)
    let dfc = diffNorm2(f, f2, scr)
    ok(&"hcForce convenience overload agrees ({dfc:.2e})", dfc < 1e-20)
    # hcActionDeriv -> contract equals hcForce (structure check)
    var f3 = newOneOf(g)
    hcActionDeriv(w, beta, g, f3)
    threads:
      for mu in 0..<nDim:
        for e in f3.uA[mu]:
          let s = g.uA[mu][e] * f3.uA[mu][e].adj
          f3.uA[mu][e].projectTAH s
        for e in f3.uB[mu]:
          let s = g.uB[mu][e] * f3.uB[mu][e].adj
          f3.uB[mu][e].projectTAH s
      for d in 0..<nDiag:
        for e in f3.uD[d]:
          let s = g.uD[d][e] * f3.uD[d][e].adj
          f3.uD[d][e].projectTAH s
    let dfd = diffNorm2(f, f3, scr)
    ok(&"hcActionDeriv + projectTAH(U D^dag) == hcForce ({dfd:.2e})",
       dfd < 1e-20)

  test "3. finite-difference force check, [4,4,4,4] cells":
    let g4 = @[4, 4, 4, 4]
    let (r0, _) = fdForce(g4, seed + 33, 0, "all links")
    let (r1, _) = fdForce(g4, seed + 34, 1, "uA only")
    let (r2, _) = fdForce(g4, seed + 35, 2, "uB only")
    let (r3, _) = fdForce(g4, seed + 36, 3, "uD only")
    let (r4, _) = fdForce(g4, seed + 37, 4, "uD[5] only")
    ok(&"FD force, all links  (rel {r0:.2e})", r0 < 1e-7)
    ok(&"FD force, uA only    (rel {r1:.2e})", r1 < 1e-7)
    ok(&"FD force, uB only    (rel {r2:.2e})", r2 < 1e-7)
    ok(&"FD force, uD only    (rel {r3:.2e})", r3 < 1e-7)
    ok(&"FD force, uD[5] only (rel {r4:.2e})", r4 < 1e-7)

  test "4. finite-difference force check, [2,4,4,6] cells":
    let g2 = @[2, 4, 4, 6]
    let (r0, _) = fdForce(g2, seed + 44, 0, "all links")
    let (r4, _) = fdForce(g2, seed + 45, 4, "uD[5] only")
    ok(&"FD force, all links, L=2 wraps (rel {r0:.2e})", r0 < 1e-7)
    ok(&"FD force, uD[5] only, L=2 wraps (rel {r4:.2e})", r4 < 1e-7)
    # and once on a non-vectorised layout, to catch SIMD-lane bugs
    let (rv, _) = fdForce(@[4, 4, 4, 4], seed + 46, 0, "all, V=1", 1)
    ok(&"FD force, all links, V = 1 layout (rel {rv:.2e})", rv < 1e-7)

  test "5. gauge invariance of hcAction and covariance of hcForce":
    let hl = newHcLayout([4, 4, 4, 6])
    let lo = hl.lo
    var r = lo.newRNGField(RngMilc6, seed + 55)
    var g = newHcGauge(hl)
    threads:
      g.random r
    var w = newHcActionWork(g)
    let s0 = hcAction(w, beta, g)
    var f0 = newOneOf(g)
    hcForce(w, beta, g, f0)
    var vA = lo.ColorMatrix(nc)
    var vB = lo.ColorMatrix(nc)
    threads:
      vA.randomSU r
      vB.randomSU r
    g.gaugeTransform(vA, vB)
    let s1 = hcAction(w, beta, g)
    let ds = abs(s1 - s0)/abs(s0)
    echo &"  hcAction before = {s0:.15g}"
    echo &"  hcAction after  = {s1:.15g}"
    ok(&"gauge invariant to 1e-12 (rel {ds:.2e})", ds < 1e-12)
    # force covariance: f -> V(start) f V(start)^dag
    var f1 = newOneOf(g)
    hcForce(w, beta, g, f1)
    var t = vA.newOneOf
    threads:
      for mu in 0..<nDim:
        t := vA * f0.uA[mu]
        f0.uA[mu] := t * vA.adj
        t := vB * f0.uB[mu]
        f0.uB[mu] := t * vB.adj
      for d in 0..<nDiag:
        t := vB * f0.uD[d]
        f0.uD[d] := t * vB.adj
    var scr = newOneOf(g)
    let d2 = diffNorm2(f1, f0, scr)
    let n2 = redot(f1, f1)
    echo &"  |f' - V f V^dag|^2 = {d2:.3e}   |f'|^2 = {n2:.6g}"
    ok(&"force is gauge covariant (rel {d2/n2:.2e})", d2/n2 < 1e-18)

  test "6. classical continuum limit / beta normalisation vs cubic Wilson":
    # Weak Abelian plane wave A_mu(x) = eps_mu cos(p.x) embedded via
    # T = diag(1,-1,0), links from exact straight-line integrals, same beta on
    # an Ns^4-cell honeycomb and an Ns^4 cubic lattice.  The ratio must -> 1
    # as eps -> 0 and p -> 0 with O((pa)^2) deviations.
    let gact = GaugeActionCoeffs(plaq: beta)
    proc measure(ns: int, k: int, pd: int, epsHat: array[4, float],
                 eps0: float): tuple[sh, sc, pred: float] =
      var pv: array[4, float]
      pv[pd] = 2.0*PI*float(k)/float(ns)
      var eps: array[4, float]
      for mu in 0..<4: eps[mu] = eps0*epsHat[mu]
      # honeycomb
      let hl = newHcLayout([ns, ns, ns, ns])
      var hg = newHcGauge(hl)
      setAbelianHc(hg, eps, pv)
      var w = newHcActionWork(hg)
      let sh = hcAction(w, beta, hg)
      # cubic
      let lc = newLayout(@[ns, ns, ns, ns])
      var gc = lc.newGauge
      setAbelianCubic(gc, lc, eps, pv)
      let sc = wilsonAction(gact, gc)
      # exact small-eps cubic value for this transverse single mode
      var e2 = 0.0
      for mu in 0..<4:
        if mu != pd: e2 += eps[mu]*eps[mu]
      doAssert abs(eps[pd]) < 1e-30   # transverse setup
      let pred = (2.0*beta/3.0)*sin(0.5*pv[pd])^2*e2*float(lc.physVol)
      (sh, sc, pred)

    echo "  Ns  dir k   eps0     S_16cell        S_cubic         ratio",
         "        |ratio-1|   |ratio-1|/p^2"
    var devs: array[5, float]
    var ratios: array[5, float]
    var p2s: array[5, float]
    let cases = [
      (12, 1, 0, [0.0, 1.0, -0.6, 0.3], 1e-3),   # A1: smallest p
      (12, 2, 0, [0.0, 1.0, -0.6, 0.3], 1e-3),   # A2: doubled p
      (12, 1, 0, [0.0, 1.0, -0.6, 0.3], 3e-2),   # A3: larger amplitude
      (12, 1, 3, [0.5, -1.0, 0.25, 0.0], 1e-3),  # B1: p along time
      (8,  1, 0, [0.0, 1.0, -0.6, 0.3], 1e-3)]   # C1: coarser p
    var cubicPredOK = true
    for i in 0..<cases.len:
      let (ns, k, pd, eh, e0) = cases[i]
      let (sh, sc, pred) = measure(ns, k, pd, eh, e0)
      let ratio = sh/sc
      let p = 2.0*PI*float(k)/float(ns)
      devs[i] = abs(ratio - 1.0)
      ratios[i] = ratio
      p2s[i] = p*p
      echo &"  {ns:3d} {pd:3d} {k:2d}  {e0:7.0e}  {sh:.8e}  {sc:.8e}  ",
           &"{ratio:.8f}  {devs[i]:.4e}  {devs[i]/(p*p):.5f}"
      if e0 <= 1e-3 and abs(sc/pred - 1.0) > 1e-4:
        cubicPredOK = false
        echo &"    cubic S = {sc:.8e} vs analytic {pred:.8e} ",
             &"(rel {abs(sc/pred-1.0):.2e})"
    ok("cubic action matches the analytic small-eps value to 1e-4",
       cubicPredOK)
    ok(&"ratio -> 1 at the smallest p (|ratio-1| = {devs[0]:.3e} < 0.25)",
       devs[0] < 0.25)
    let scale21 = devs[1]/devs[0]
    ok(&"momentum scaling k=1 -> k=2: dev ratio {scale21:.3f} ~ 4 (O(p^2))",
       scale21 > 2.8 and scale21 < 5.8)
    let scaleNs = devs[4]/devs[0]
    ok(&"volume scaling Ns=12 -> Ns=8: dev ratio {scaleNs:.3f} ~ 2.25 (O(p^2))",
       scaleNs > 1.6 and scaleNs < 3.2)
    let dAmp = abs(ratios[2] - ratios[0])
    ok(&"amplitude independence: |ratio(3e-2) - ratio(1e-3)| = {dAmp:.2e} " &
       "(O(eps^2))", dAmp < 5e-3)
    ok(&"time-direction case consistent (|ratio-1| = {devs[3]:.3e})",
       devs[3] < 0.25)

qexFinalize()
