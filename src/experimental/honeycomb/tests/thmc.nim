## Task M test: `hchmc.nim` (HMC), `hcheatbath.nim` (Cabibbo-Marinari
## heatbath + overrelaxation) and `hcio.nim` (configuration I/O).
##
##   export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
##   cd build_mac && make src/experimental/honeycomb/tests/thmc.nim && \
##     OMP_NUM_THREADS=4 ./bin/thmc
##
## The decisive checks:
##   * the update-schedule structure (every triangle has exactly 1 axis edge;
##     its 2 diagonal edges have different delta) is verified from
##     hcgeom.triPath, not assumed;
##   * the heatbath staple fields are pinned bit-tightly to hcActionDeriv
##     (beta = 2N), and the local action form -(beta/2N) ReTr(U Sigma^dag) is
##     pinned by single-link finite differences of hcAction;
##   * the Kennedy-Pendleton SU(2) sampler is tested against the exact
##     moments <a0> = I2(alpha)/I1(alpha);
##   * long heatbath+OR and HMC runs at the same beta must agree on
##     <triangleSum> within errors.
##
## beta = 8.0 is used throughout: from the heatbath scan (doc/plots/
## hb_scan.dat, 4^4 cells) <triangleSum>(8.0) = 0.5893, inside the requested
## 0.5..0.65 crossover window, safely above the steep crossover at beta ~ 7,
## with tau_int ~ 1 for the heatbath.

import std/[math, strformat, strutils, os, unittest]
import qex except epsilon
import physics/qcdTypes
import ../hcgeom
import ../hclayout
import ../hcgauge
import ../hcaction
import ../hchmc
import ../hcheatbath
import ../hcio
import ../hcanalysis

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

template ok(msg: string, cond: bool) =
  let c = cond
  if c: echo "PASS: ", msg
  else: echo "FAIL: ", msg
  check c

qexInit()

const nc = getDefaultNc()
let
  beta = 8.0            ## from the hb_scan crossover window (module docs)
  seed = 246813579'u64

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

type Mat = array[nc, array[nc, array[2, float]]]

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

proc setMat(f: auto, idx: int, m: Mat) =
  for a in 0..<nc:
    for b in 0..<nc:
      f{idx}[a, b].re := m[a][b][0]
      f{idx}[a, b].im := m[a][b][1]

proc reTrMulAdj(a, b: Mat): float =
  ## Re Tr(a b^dag)
  for i in 0..<nc:
    for j in 0..<nc:
      result += a[i][j][0]*b[i][j][0] + a[i][j][1]*b[i][j][1]

proc diffNorm2(a, b, scratch: auto): float =
  threads:
    for mu in 0..<nDim:
      scratch.uA[mu] := a.uA[mu] - b.uA[mu]
      scratch.uB[mu] := a.uB[mu] - b.uB[mu]
    for d in 0..<nDiag:
      scratch.uD[d] := a.uD[d] - b.uD[d]
  redot(scratch, scratch)

template triSumOf(w, beta, g: untyped): float =
  1.0 - hcAction(w, beta, g)/(0.5*beta*float(nTriPerSite*g.hl.nSites))

proc binErr(x: seq[float]): tuple[m, e, tau: float] =
  ## mean and autocorrelation-corrected error
  let tau = if x.len > 8: autocorrTime(x) else: 1.0
  let (m, e0) = jackknifeMean(x)
  (m, e0*sqrt(max(1.0, 2.0*tau)), tau)

proc besselIRatio21(a: float): float =
  ## I2(a)/I1(a) by the ascending series (fine for 0 < a <~ 100)
  let x = 0.5*a
  var t = x            # I1 k=0 term: x^1/(0! 1!)
  var s1 = t
  for k in 1..300:
    t *= x*x/(k.float*(k.float + 1.0))
    s1 += t
    if t < 1e-18*s1: break
  t = 0.5*x*x          # I2 k=0 term: x^2/(0! 2!)
  var s2 = t
  for k in 1..300:
    t *= x*x/(k.float*(k.float + 2.0))
    s2 += t
    if t < 1e-18*s2: break
  s2/s1

proc a0MeanQuad(alpha: float): float =
  ## <a0> for P(a0) ~ sqrt(1-a0^2) exp(alpha a0), by Simpson quadrature in
  ## theta (a0 = cos theta removes the endpoint sqrt singularities; the
  ## integrand is scaled by exp(-alpha), which cancels in the ratio)
  let n = 20000
  let h = PI/n.float
  var s0, s1 = 0.0
  for i in 0..n:
    let th = i.float*h
    let x = cos(th)
    let st = sin(th)
    let f = st*st*exp(alpha*(x - 1.0))
    let wgt = if i == 0 or i == n: 1.0 elif (i and 1) == 1: 4.0 else: 2.0
    s0 += wgt*f
    s1 += wgt*x*f
  s1/s0

# ---------------------------------------------------------------------------

suite "thmc":

  test "1. update-schedule structure from hcgeom.triPath":
    # (a) every triangle has exactly ONE axis edge -> axis staples contain no
    #     axis links -> all 8 axis fields can be updated simultaneously;
    # (b) the two diagonal edges of any triangle carry different delta
    #     indices -> each uD[delta] field can be updated in one pass.
    var nTri = 0
    var okA = true
    var okB = true
    for sub in 0..1:
      for t in apexTris:
        let path = triPath(Site(cell: [0, 0, 0, 0], sub: sub), t)
        var axes = 0
        var ds: seq[int]
        for l in path:
          if l.kind == lkD: ds.add l.idx
          else: inc axes
        if axes != 1 or ds.len != 2: okA = false
        if ds.len == 2 and ds[0] == ds[1]: okB = false
        inc nTri
    ok(&"scanned {nTri} triangles: every one has exactly 1 axis + 2 diagonal edges", okA and nTri == 64)
    ok("the 2 diagonal edges of a triangle always have different delta", okB)

  test "2. heatbath staple fields match hcActionDeriv(beta=2N)":
    let hl = newHcLayout([4, 4, 4, 6])
    var r = hl.lo.newRNGField(MRG32k3a, seed + 11)
    var g = newHcGauge(hl)
    threads:
      g.random r
    var w = newHcActionWork(g)
    var f = newOneOf(g)
    # hcActionDeriv prefactor is beta/(2N): beta = 2N gives the raw staple
    # sums Sigma_l = sum_{k=1}^{8} V_k
    hcActionDeriv(w, 2.0*float(nc), g, f)
    var hb = newHcHeatbath(g, beta)
    axisStaples(hb, g)
    threads:
      hb.refreshTrees
    var t = g.uA[0].newOneOf
    var worst = 0.0
    for mu in 0..<nDim:
      var dA, dB, nA, nB: float
      threads:
        t := hb.stA[mu] - f.uA[mu]
        threadMaster: dA = 0.0
      dA = t.norm2
      nA = f.uA[mu].norm2
      threads:
        t := hb.stB[mu] - f.uB[mu]
      dB = t.norm2
      nB = f.uB[mu].norm2
      worst = max(worst, max(dA/nA, dB/nB))
    for d in 0..<nDiag:
      dStaple(hb, g, d)
      var dd, nd: float
      threads:
        t := hb.stD - f.uD[d]
      dd = t.norm2
      nd = f.uD[d].norm2
      worst = max(worst, dd/nd)
    echo &"  worst relative |staple - hcActionDeriv|^2 over 24 fields: {worst:.3e}"
    ok(&"heatbath staples == hcActionDeriv(2N) staple sums (rel^2 {worst:.2e})",
       worst < 1e-24)

  test "3. local action is -(beta/2N) ReTr(U Sigma^dag): single-link pin":
    let hl = newHcLayout([4, 4, 4, 6])
    var r = hl.lo.newRNGField(MRG32k3a, seed + 22)
    var g = newHcGauge(hl)
    threads:
      g.random r
    var w = newHcActionWork(g)
    var f = newOneOf(g)
    hcActionDeriv(w, 2.0*float(nc), g, f)   # f = raw Sigma per link
    let cLoc = beta/(2.0*float(nc))
    let s0 = hcAction(w, beta, g)
    var worst = 0.0
    proc pinEcho(lbl: string, dsMeas, dsPred, rel: float) =
      echo &"  {lbl}: dS = {dsMeas: .10e}  pred = {dsPred: .10e}  rel {rel:.2e}"
    # replace U at one site of one field by the U from another site (stays
    # unitary), predict dS from Sigma, compare with the full action
    template pin(uf, sf, lbl: untyped) =
      block:
        let i = 3
        let j = 3 + hl.lo.V   # a different site
        let uOld = getMat(uf, i)
        let uNew = getMat(uf, j)
        let sig = getMat(sf, i)
        let dsPred = -cLoc*(reTrMulAdj(uNew, sig) - reTrMulAdj(uOld, sig))
        setMat(uf, i, uNew)
        let s1 = hcAction(w, beta, g)
        setMat(uf, i, uOld)   # restore
        let dsMeas = s1 - s0
        let rel = abs(dsMeas - dsPred)/max(1e-30, abs(dsMeas))
        pinEcho(lbl, dsMeas, dsPred, rel)
        worst = max(worst, rel)
    pin(g.uA[2], f.uA[2], "uA[2]")
    pin(g.uB[1], f.uB[1], "uB[1]")
    pin(g.uD[5], f.uD[5], "uD[5]")
    pin(g.uD[10], f.uD[10], "uD[10]")
    ok(&"single-link dS matches -(beta/2N) ReTr(dU Sigma^dag) (worst rel {worst:.2e})",
       worst < 1e-10)

  test "4. HMC reversibility":
    let hl = newHcLayout([4, 4, 4, 6])
    var r = hl.lo.newRNGField(MRG32k3a, seed + 33)
    var g = newHcGauge(hl)
    threads:
      g.warm(0.35, r)
    for algo in ["leapfrog", "2MN", "4MN5FV"]:
      var h = newHcHmc(g, beta, 1.0, 8, algo)
      let rc = h.revCheck(g, r)
      echo &"  {algo:9s}: dH_fwd {rc.dHf: .6e}  |dH_fwd+dH_bwd| {abs(rc.sumdH):.3e}  ",
           &"per-link |U_back-U_0| {rc.linkDiff:.3e}"
      ok(&"{algo}: |dH_fwd+dH_bwd| = {abs(rc.sumdH):.2e} < 1e-8",
         abs(rc.sumdH) < 1e-8)
      ok(&"{algo}: per-link |U_back-U_0| = {rc.linkDiff:.2e} < 1e-10",
         rc.linkDiff < 1e-10)

  test "5. dH ~ O(eps^2) for leapfrog at fixed tau":
    let hl = newHcLayout([4, 4, 4, 4])
    var r = hl.lo.newRNGField(MRG32k3a, seed + 44)
    var g = newHcGauge(hl)
    # semi-realistic start: short heatbath thermalisation
    var hb = newHcHeatbath(g, beta)
    threads:
      g.random r
    for n in 1..40:
      hb.update(g, r, 1)
    var h = newHcHmc(g, beta, 1.0, 8, "leapfrog")
    var gs = newOneOf(g)
    var ps = newOneOf(g)
    let p = h.p
    var r2 = hl.lo.newRNGField(MRG32k3a, seed + 45)
    threads:
      eachLink(p, u):
        u.randomTAH r2
      gs := g
      ps := p
    proc dhAt(ns: int): float =
      h.nsteps = ns
      h.sched = mdSchedule("leapfrog", ns)
      threads:
        g := gs
        h.p := ps
      let (_, _, h0) = h.hamiltonian(g)
      h.integrate(g)
      let (_, _, h1) = h.hamiltonian(g)
      threads:
        g := gs
      h1 - h0
    let d16 = dhAt(16)
    let d32 = dhAt(32)
    let d64 = dhAt(64)
    let r1 = d16/d32
    let r2r = d32/d64
    echo &"  dH(nsteps=16) = {d16: .6e}"
    echo &"  dH(nsteps=32) = {d32: .6e}   ratio {r1:.3f}"
    echo &"  dH(nsteps=64) = {d64: .6e}   ratio {r2r:.3f}"
    # leapfrog is symmetric: dH = c2 eps^2 + c4 eps^4 + ..., so the ratio
    # approaches 4 from above as eps shrinks
    ok(&"halving eps quarters dH: ratios {r1:.2f}, {r2r:.2f} in [3.3,5.3]",
       r1 > 3.3 and r1 < 5.3 and r2r > 3.3 and r2r < 5.3)

  test "6. Kennedy-Pendleton SU(2) sampler vs exact moments":
    # P(a0) ~ sqrt(1-a0^2) exp(alpha a0)  =>  <a0> = I2(alpha)/I1(alpha).
    # alpha = 0.5 exercises the plain-rejection branch, the rest the KP
    # branch (switch at alpha = 1).
    var R: MRG32k3a
    R.seed(seed, 4242)
    var worst = 0.0
    for alpha in [0.5, 2.0, 8.0, 20.0]:
      let exact = besselIRatio21(alpha)
      let quad = a0MeanQuad(alpha)
      doAssert abs(exact - quad) < 1e-8   # two independent exact evaluations
      const nSamp = 4_000_000
      var s, s2 = 0.0
      for k in 0..<nSamp:
        let x = sampleA0(R, alpha)
        s += x
        s2 += x*x
      let m = s/nSamp.float
      let sig = sqrt((s2/nSamp.float - m*m)/nSamp.float)
      let d = abs(m - exact)
      echo &"  alpha {alpha:5.1f}: <a0> = {m:.6f} +- {sig:.6f}   I2/I1 = {exact:.6f}   |diff| = {d:.2e} ({d/sig:.1f} sigma)"
      worst = max(worst, d)
      ok(&"alpha {alpha}: |<a0> - I2/I1| = {d:.1e} < max(1e-3, 5 sigma)",
         d < max(1e-3, 5.0*sig))
    # the full quaternion update with a rotated staple: a ~ exp(c ReTr(a W))
    # with W = kq * (random SU(2)); ReTr(aW)/(2 kq) must average to I2/I1.
    block:
      let kq = 3.7
      let c = 0.9
      let alpha = 2.0*c*kq
      let exact = besselIRatio21(alpha)
      # random unit quaternion for the staple direction
      var v: Quat
      block:
        var n2 = 0.0
        for i in 0..3:
          v[i] = gaussian(R)
          n2 += v[i]*v[i]
        let s = 1.0/sqrt(n2)
        for i in 0..3: v[i] *= s
      let w: Quat = [kq*v[0], kq*v[1], kq*v[2], kq*v[3]]
      const nSamp = 2_000_000
      var s, s2 = 0.0
      for k in 0..<nSamp:
        let a = su2HeatbathQuat(R, w, c)
        # ReTr(a W) = 2 (a0 w0 - a1 w1 - a2 w2 - a3 w3)
        let b0 = (a[0]*w[0] - a[1]*w[1] - a[2]*w[2] - a[3]*w[3])/kq
        s += b0
        s2 += b0*b0
      let m = s/nSamp.float
      let sig = sqrt((s2/nSamp.float - m*m)/nSamp.float)
      let d = abs(m - exact)
      echo &"  rotated staple, alpha {alpha:.2f}: <ReTr(aW)>/2kq = {m:.6f} +- {sig:.6f}  I2/I1 = {exact:.6f}"
      ok(&"full quaternion update reproduces I2/I1 ({d:.1e} < max(1e-3, 5 sigma))",
         d < max(1e-3, 5.0*sig))

  test "7. OR sweep: action invariant, links change":
    let hl = newHcLayout([4, 4, 4, 4])
    var r = hl.lo.newRNGField(MRG32k3a, seed + 66)
    var g = newHcGauge(hl)
    var hb = newHcHeatbath(g, beta)
    threads:
      g.random r
    for n in 1..30:
      hb.update(g, r, 1)
    var g0 = newOneOf(g)
    var scr = newOneOf(g)
    threads:
      g0 := g
    let s0 = hcAction(hb.w, beta, g)
    hb.orSweep(g)
    let s1 = hcAction(hb.w, beta, g)
    let rel = abs(s1 - s0)/abs(s0)
    let du = diffNorm2(g, g0, scr)/float(nLinks(hl))
    echo &"  S before = {s0:.10f}"
    echo &"  S after  = {s1:.10f}   rel change {rel:.3e}"
    echo &"  per-link |U'-U|^2 = {du:.4f}  (0 would mean no update)"
    ok(&"one OR sweep changes hcAction by {rel:.1e} < 1e-8 relative", rel < 1e-8)
    ok(&"links change substantially (per-link |U'-U|^2 = {du:.2f} > 0.1)",
       du > 0.1)
    let dsu = g.checkSU
    echo &"  checkSU after OR: avg {dsu.avg:.2e} max {dsu.max:.2e}"
    # the quaternion row updates drift off SU(3) by O(1e-14) per sweep;
    # production streams reunit periodically (hcPureGauge does)
    ok("links stay in SU(3) to 1e-10 without reunit", dsu.max < 1e-10)

  test "8. HMC: <exp(-dH)> = 1, acceptance; heatbath agrees on <triangleSum>":
    let hl = newHcLayout([4, 4, 4, 4])
    var r = hl.lo.newRNGField(MRG32k3a, seed + 77)
    var R: MRG32k3a
    R.seed(seed, 8888)
    var g = newHcGauge(hl)
    var hb = newHcHeatbath(g, beta)
    # thermalise once with the heatbath, keep a copy as common start
    threads:
      g.random r
    for n in 1..150:
      hb.update(g, r, 2)
    var gTherm = newOneOf(g)
    threads:
      gTherm := g
    # ---- HMC stream --------------------------------------------------
    var h = newHcHmc(g, beta, 1.0, 10, "2MN")
    var dhs, edhs, tsH: seq[float]
    var nacc = 0
    const nTraj = 300
    for n in 1..nTraj:
      let (dH, acc, accepted) = h.trajectory(g, r, R)
      if accepted: inc nacc
      dhs.add dH
      edhs.add acc
      tsH.add triSumOf(h.w, beta, g)
    let accRate = nacc.float/nTraj.float
    let (edhM, edhE, _) = binErr(edhs)
    let (tsHM, tsHE, tauH) = binErr(tsH)
    echo &"  HMC: 2MN tau 1 nsteps 10, {nTraj} trajectories, beta {beta}"
    echo &"  acceptance = {100.0*accRate:.2f} %   <dH> = {dhs.mean:.5f} +- {stderrMean(dhs):.5f}"
    echo &"  <exp(-dH)> = {edhM:.5f} +- {edhE:.5f}"
    echo &"  <triangleSum>_HMC = {tsHM:.6f} +- {tsHE:.6f}  (tau_int {tauH:.2f})"
    ok(&"acceptance {100.0*accRate:.1f} % > 70 %", accRate > 0.70)
    ok(&"<exp(-dH)> - 1 = {edhM-1.0:.2e} within 2 sigma ({edhE:.1e})",
       abs(edhM - 1.0) < 2.0*edhE)
    # ---- heatbath stream from the same thermalised start -------------
    threads:
      g := gTherm
    var tsB: seq[float]
    const nUpd = 800
    for n in 1..nUpd:
      hb.update(g, r, 2)
      tsB.add triSumOf(hb.w, beta, g)
    let (tsBM, tsBE, tauB) = binErr(tsB)
    echo &"  HB: 1 heatbath + 2 OR per update, {nUpd} updates"
    echo &"  <triangleSum>_HB  = {tsBM:.6f} +- {tsBE:.6f}  (tau_int {tauB:.2f})"
    let dTs = abs(tsHM - tsBM)
    let sig = sqrt(tsHE*tsHE + tsBE*tsBE)
    echo &"  |HMC - HB| = {dTs:.6f} = {dTs/sig:.2f} combined sigma"
    ok(&"HMC and heatbath <triangleSum> agree within 3 sigma ({dTs/sig:.2f})",
       dTs < 3.0*sig)

  test "9. I/O round trip is bit exact":
    let hl = newHcLayout([4, 4, 4, 6])
    var r = hl.lo.newRNGField(MRG32k3a, seed + 88)
    var g = newHcGauge(hl)
    threads:
      g.random r
    var w = newHcActionWork(g)
    let ts0 = triSumOf(w, beta, g)
    let s0 = hcAction(w, beta, g)
    let fn = getTempDir() / "thmc_io_test.lime"
    ok("saveHcGauge returns 0",
       0 == g.saveHcGauge(fn, beta = beta, traj = 137, info = "thmc test 9"))
    ok(&"stored cell geometry reads back as {hl.geom}", hcFileGeom(fn) == hl.geom)
    var g2 = newHcGauge(hl)
    let (st, meta) = g2.loadHcGauge(fn)
    ok("loadHcGauge returns 0", st == 0)
    ok(&"metadata round trip (beta {meta.beta}, traj {meta.traj}, geom {meta.geom})",
       meta.beta == beta and meta.traj == 137 and meta.geom == hl.geom and
       meta.version == 1)
    var scr = newOneOf(g)
    let d2 = diffNorm2(g, g2, scr)
    let ts1 = triSumOf(w, beta, g2)
    let s1 = hcAction(w, beta, g2)
    echo &"  sum|U-U'|^2 = {d2:.1e}   triangleSum {ts0:.16f} -> {ts1:.16f}"
    ok("configuration is bit exact (sum|U-U'|^2 == 0.0)", d2 == 0.0)
    ok("triangleSum and action identical",
       ts1 == ts0 and s1 == s0)
    removeFile fn

qexFinalize()
