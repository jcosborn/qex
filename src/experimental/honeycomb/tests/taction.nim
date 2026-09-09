## `hcaction.nim` (triangle action, staple derivative, force).
##
## The make-or-break test is the finite-difference force check: for random
## traceless anti-Hermitian momenta P on all 24 link fields,
##   d/ds S(exp(s P) U)|_0  must equal  sum_l redot(P_l, f_l)
## to ~1e-7 or better, on a warm (non-unit) configuration, for the full
## momentum and restricted to each link kind separately.  The same identity is
## first verified for QEX's own `gaugeForce` on a cubic lattice, so the force
## convention is pinned to QEX's, not assumed.

import std/[math, strformat, unittest]
import qex except epsilon
import physics/qcdTypes
import algorithms/numdiff
import ../hcgeom
import ../hcgauge
import ../hcaction
import helpers

qexInit()

let
  beta = 5.7
  seed = 135792468'u64

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

proc setMom(p: auto, r: auto, sel: int) =
  ## random TAH momenta on a subset of the 24 link fields:
  ## 0 = all, 1 = uA only, 2 = uB only, 3 = uD only, 4 = uD[5] only
  threads:
    for u in p.links: u := 0
    case sel
    of 0:
      p.randomTAH r
    of 1:
      for mu in 0..<nDim: p.uA[mu].randomTAH r
    of 2:
      for mu in 0..<nDim: p.uB[mu].randomTAH r
    of 3:
      for d in 0..<nDiag: p.uD[d].randomTAH r
    else:
      p.uD[5].randomTAH r

proc fdForce(geom: seq[int], seed: uint64, sel: int, label: string,
             V: static int = VLEN): tuple[rel, nde: float] =
  ## finite-difference check of force on a warm configuration
  let hl = newLayout(geom, V)
  let lo = hl
  var r = lo.newRNGField(RngMilc6, seed)
  var g0 = newHcGauge(hl)
  var gt = newOneOf(g0)
  var p = newOneOf(g0)
  var f = newOneOf(g0)
  threads:
    g0.warm(0.35, r)
  setMom(p, r, sel)
  var w = newActionWork(g0)
  force(w, beta, g0, f)
  let pf = redot(p, f)
  proc sAt(t: float): float =
    expMul(gt, g0, p, t)
    action(w, beta, gt)
  var dsdt, e: float
  ndiff(dsdt, e, sAt, 0.0, 0.1, ordMax = 7)
  let rel = abs(dsdt - pf)/max(abs(dsdt), abs(pf))
  echo &"  {label:11s}: dS/dt = {dsdt: .13e}   sum redot(p,f) = {pf: .13e}"
  echo &"  {label:11s}: rel diff = {rel:.3e}  (ndiff err estimate {e:.1e})"
  (rel, e)

# ---------------------------------------------------------------------------

suite "hcaction":

  test "1. QEX gaugeForce satisfies the redot convention (cubic reference)":
    # d/ds S(exp(sP)U)|_0 = sum_mu redot(P_mu, F_mu) for QEX's own Wilson
    # action and force -- this is the convention force must reproduce.
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
    let hl = newLayout([4, 4, 4, 6])
    let lo = hl
    var r = lo.newRNGField(RngMilc6, seed + 22)
    var g = newHcGauge(hl)
    var w = newActionWork(g)
    let sTri = float(nTriPerSite*2*hl.physVol)
    # unit gauge: S = 0, force = 0
    let s0 = action(w, beta, g)
    echo &"  action(unit) = {s0:.3e}  (scale beta/2*nTri = {0.5*beta*sTri:.4g})"
    ok(&"unit gauge S = 0 ({abs(s0):.2e})", abs(s0) < 1e-8)
    var f = newOneOf(g)
    force(w, beta, g, f)
    let f0n = redot(f, f)
    echo &"  |force(unit)|^2 = {f0n:.3e}"
    ok(&"unit gauge force = 0 ({f0n:.2e})", f0n < 1e-24)
    # warm configuration: independent cross-check against hcgauge's
    # triangleSum (verified in tgauge.nim against a brute-force reference)
    threads:
      g.warm(0.35, r)
    let sw = action(w, beta, g)
    let swRef = 0.5*beta*sTri*(1.0 - g.triangleSum)
    let dw = abs(sw - swRef)/abs(swRef)
    echo &"  warm:   hcAction = {sw:.15g}"
    echo &"          (beta/2)*32*nSites*(1-triangleSum) = {swRef:.15g}"
    ok(&"warm action matches triangleSum path (rel {dw:.2e})", dw < 1e-10)
    # random configuration
    threads:
      g.random r
    let sr = action(w, beta, g)
    let srRef = 0.5*beta*sTri*(1.0 - g.triangleSum)
    let dr = abs(sr - srRef)/abs(srRef)
    echo &"  random: hcAction = {sr:.15g}"
    echo &"          (beta/2)*32*nSites*(1-triangleSum) = {srRef:.15g}"
    ok(&"random action matches triangleSum path (rel {dr:.2e})", dr < 1e-10)
    ok(&"hcAction > 0 on random ({sr:.6g}) and warm ({sw:.6g})",
       sr > 0.0 and sw > 0.0)
    ok(&"warm action below random action ({sw:.4g} < {sr:.4g})", sw < sr)
    force(w, beta, g, f)
    # actionDeriv -> contract equals force (structure check)
    var f3 = newOneOf(g)
    actionDeriv(w, beta, g, f3)
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
    let dfd = norm2diff(f, f3)
    ok(&"actionDeriv + projectTAH(U D^dag) == force ({dfd:.2e})",
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

  test "5. gauge invariance of hcAction and covariance of force":
    let hl = newLayout([4, 4, 4, 6])
    let lo = hl
    var r = lo.newRNGField(RngMilc6, seed + 55)
    var g = newHcGauge(hl)
    threads:
      g.random r
    var w = newActionWork(g)
    let s0 = action(w, beta, g)
    var f0 = newOneOf(g)
    force(w, beta, g, f0)
    var vA = lo.ColorMatrix(nc)
    var vB = lo.ColorMatrix(nc)
    threads:
      vA.randomSU r
      vB.randomSU r
    g.gaugeTransform(vA, vB)
    let s1 = action(w, beta, g)
    let ds = abs(s1 - s0)/abs(s0)
    echo &"  hcAction before = {s0:.15g}"
    echo &"  hcAction after  = {s1:.15g}"
    ok(&"gauge invariant to 1e-12 (rel {ds:.2e})", ds < 1e-12)
    # force covariance: f -> V(start) f V(start)^dag
    var f1 = newOneOf(g)
    force(w, beta, g, f1)
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
    let d2 = norm2diff(f1, f0)
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
      let hl = newLayout([ns, ns, ns, ns])
      var hg = newHcGauge(hl)
      setAbelianHc(hg, eps, pv)
      var w = newActionWork(hg)
      let sh = action(w, beta, hg)
      # cubic
      let lc = newLayout(@[ns, ns, ns, ns])
      var gc = lc.newGauge
      setAbelianCubic(gc, lc, eps, pv)
      let sc = wilsonAction(gact, gc)
      # exact small-eps cubic value for this transverse single mode
      var e2 = 0.0
      for mu in 0..<4:
        if mu != pd: e2 += eps[mu]*eps[mu]
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
