## WP-G report: the free gauge current correlator (T1.5a/T1.5b) and the gradient-flow
## scale scan of slide 9 (T2.2).  Everything here is deterministic linear algebra
## except the flow scan, which uses the exact Gaussian heatbath.
##
## Job 1.  G_g(t) = (1/g^2) <J^t(t) J^t(0)>, J^t = Theta_tri/A_tri (V.12), averaged
##   over source triangles, from the double-CG pseudo-inverse (V.16)-(V.17).  Compared
##   to the continuum (V.14) via core/analytic's gaugeG/gaugeGPeriodic, then
##   Delta_eff by meas/fit.effMass and Delta_0 by meas/fit.plateauFit.  Run under all
##   three geometry conventions of `GeomConv`, which differ at O(abar^2) -- the paper's
##   section III says the flat one is "equally possible".
##
## Job 2.  E_s(s) sqrt(L) versus r/s on log-log, L x g^2 grid, on heatbath ensembles.
##
## Output: TSV under output/radial/gauge/.

import base
import std/[math, os, strformat, strutils]
import core/analytic
import ops/flow
import meas/[dataio, fit]

qexInit()
freezeTimers()

letParam:
  lev = 2            ## correlator: L = 1, 2, 4, ... up to this level
  nt = 120           ## time slices (paper: 120 for the gauge correlator)
  ttot = 16.0        ## temporal extent T; a_t = ttot/nt
  g2inv = 20.0       ## 1/g^2 (paper: 20).  G_g is g-independent; checked below.
  fitLo = 4.0        ## Delta_eff plateau fit window [fitLo, fitHi)
  fitHi = 8.0
  nsrc = 100000      ## cap on the number of source triangles averaged over
  r2req = 1e-26
  maxits = 400000
  flowLev = 4        ## flow scan: L = 1, 2, ... up to this level
  flowNt = 32
  flowAt = 0.2
  flowMax = 3.0      ## largest flow time
  nflow = 30         ## log-spaced flow times
  nconf = 24         ## heatbath configurations per point
  seed = 987654321
  doCorr = 1
  doCont = 1         ## O(a^2) scaling / continuum extrapolation over (L, L_t)
  contLev = 4
  contFlat = 0       ## also run the continuum fit in the exact-area and flat conventions
  doFlow = 1

installStandardParams()
echoParams()
processHelpParam()

const
  outDir = currentSourcePath().parentDir.parentDir.parentDir.parentDir /
           "output" / "radial" / "gauge"
  convTag: array[GeomConv, string] = ["sph", "exact", "flat"]

# Independent reference for L=1, nt=120, T=16 (a_t = 0.1333...), computed in Python by
# Fourier transforming in t and inverting the per-momentum (N_E + N_V) matrix with the
# gauge directions and the constant-temporal flat direction regularized away.
const
  refIdx = [0, 4, 8, 15, 23, 30, 38, 45]
  refG = [1.96741714e+00, 4.62879753e-01, 1.27395702e-01, 1.90361979e-02,
          3.18269428e-03, 7.89217296e-04, 1.74510635e-04, 4.82849711e-05]
  refD = [2.816897, 2.532089, 2.233189, 1.809259, 1.545205, 1.441194, 1.390088,
          1.371383]
  published = 1.33242        ## paper Sec. V B, Delta_0 at L=1, L_t=120

let at = ttot/float(nt)

proc levels(top: int): seq[int] =
  var l = 1
  while l <= top:
    result.add l
    l *= 2

# --- job 1: free gauge current correlator -----------------------------------

type CorrOut = object
  g: seq[float]              ## G_g(t), t = i*at, one period
  nused: int
  its: int
  maxr2: float

proc gaugeCorr(lv, ntv: int, atv, g2: float, conv: GeomConv, ns: int,
               escale = 1.0): CorrOut =
  ## <Theta_f(t) Theta_f(0)> = c_f(t)^T Mtilde^{-1} c_f(0), so one pseudo-inverse
  ## solve per source triangle gives the whole time dependence at once.
  ## The source is a row of the incidence matrix, hence exactly orthogonal to
  ## ker M, so the kernel-regularized operator can be used directly and the solve
  ## is stable down to the roundoff floor -- which matters, because G_g(T/2) is
  ## five orders of magnitude below G_g(0) and the plateau fit lives out there.
  let
    sph = newSphere(lv)
    lat = newLat(sph, ntv, atv)
  var bt = newBeta(lat, g2, conv)
  # `escale` rescales beta_l alone: at L = 1 every face and every edge is in one
  # icosahedral orbit, so beta_l/beta_tri is the ONLY shape parameter of the action
  # besides a_t, and this measures how far Delta_0 moves per unit of it.
  if escale != 1.0:
    for i in 0..<bt.edge.len: bt.edge[i] *= escale
  var reg = newRegOp(lat, bt)
  result.g = newSeq[float](ntv)
  var
    b = newGauge(lat)
    x = newGauge(lat)
  let stride = max(1, (sph.nf + ns - 1) div ns)
  var f = 0
  while f < sph.nf:
    triSource(lat, b, f, 0)
    let info = regSolve(lat, x, b, reg, r2req, maxits)
    result.its += info.iters
    result.maxr2 = max(result.maxr2, info.r2)
    let ia = 1.0/(bt.afac[f]*bt.afac[f])
    for t in 0..<ntv: result.g[t] += ia*plaqSpatial(lat, x, f, t)
    inc result.nused
    f += stride
  let sc = 1.0/(g2*float(result.nused))
  for t in 0..<ntv: result.g[t] *= sc

proc reportCorr(lv: int, conv: GeomConv, c: CorrOut) =
  let
    tag = convTag[conv]
    half = nt div 2
  var gh = newSeq[float](half + 1)
  for i in 0..half: gh[i] = c.g[i]
  let
    m = effMass(gh, at, ttot)
    i0 = int(round(fitLo/at))
    i1 = min(int(round(fitHi/at)), m.len)
    fit = plateauFit(m, i0, i1, at)
  echo ""
  echo &"--- T1.5a/T1.5b  L = {lv}  geometry = {tag}  nt = {nt}  T = {ttot}  " &
       &"at = {at:.6f}  1/g^2 = {g2inv}"
  echo &"    {c.nused} source triangles, {c.its} total CG iterations, " &
       &"max relative residual {c.maxr2:.2e}"
  echo "      t        G_g(t)        (V.14) periodic     ratio     Delta_eff"
  for i in 0..half:
    if i mod 5 != 0 and i > 8: continue
    let
      t = at*float(i)
      an = gaugeGPeriodic(t, ttot)
      d = if i < m.len: m[i] else: 0.0
    echo &"  {t:7.4f}  {c.g[i]:16.8e}  {an:16.8e}  {c.g[i]/an:9.5f}  {d:10.6f}"
  echo &"    plateau fit [{fitLo}, {fitHi}) -> Delta_0 = {fit.d0:.6f} +- {fit.ed0:.6f}" &
       &"   c = {fit.c:.4g}  Delta' = {fit.dp:.4f}  chi2/dof = {fit.chi2dof:.3e}" &
       &"  (dof {fit.dof}, {fit.iters} its, converged {fit.converged})"
  echo &"    published Delta_0 = {published}  (dev {100.0*(fit.d0/published - 1.0):+.3f}%)" &
       &"   continuum sqrt(2) = {sqrt(2.0):.6f}  (dev {100.0*(fit.d0/sqrt(2.0) - 1.0):+.3f}%)"
  if lv == 1 and nt == 120 and abs(ttot - 16.0) < 1e-12:
    echo "    against the independent Python oracle (L=1, nt=120, T=16):"
    echo "       t        G_g(t) here      oracle          rel.dev      " &
         "Delta_eff here   oracle"
    for k in 0..<refIdx.len:
      let
        i = refIdx[k]
        t = at*float(i)
        d = if i < m.len: m[i] else: 0.0
      echo &"  {t:7.4f}  {c.g[i]:15.8e}  {refG[k]:15.8e}  {c.g[i]/refG[k]-1.0:+11.3e}" &
           &"   {d:12.6f}  {refD[k]:11.6f}"
  var cols = @[newSeq[float](0), newSeq[float](0), newSeq[float](0), newSeq[float](0)]
  for i in 0..<nt:
    let t = at*float(i)
    cols[0].add t
    cols[1].add c.g[i]
    cols[2].add gaugeGPeriodic(t, ttot)
    cols[3].add (if i < m.len: m[i] else: NaN)
  writeTsv(outDir / &"gcorr_L{lv}_{tag}_nt{nt}.tsv",
           {"lattice": &"L{lv}", "geometry": tag, "nt": $nt, "at": &"{at:.17g}",
            "T": &"{ttot:.17g}", "g2inv": &"{g2inv:.17g}", "nsrc": $c.nused,
            "Delta0": &"{fit.d0:.17g}", "Delta0err": &"{fit.ed0:.17g}",
            "fitLo": &"{fitLo:.17g}", "fitHi": &"{fitHi:.17g}",
            "published_L1": &"{published:.17g}"},
           ["t", "G_g", "G_analytic", "Delta_eff"], cols)
  if lv == 1:
    echo "    fit-window scan (Delta_0 is stable in the window; it is not the " &
         "window that separates us from the published value)"
    echo "      window        Delta_0      +- err        chi2/dof     Delta'"
    for (a, b) in [(0.5, 8.0), (1.0, 8.0), (2.0, 8.0), (3.0, 8.0), (4.0, 8.0),
                   (5.0, 8.0), (6.0, 8.0), (4.0, 7.0), (4.0, 6.0), (2.0, 6.0)]:
      let
        j0 = int(round(a/at))
        j1 = min(int(round(b/at)), m.len)
        fw = plateauFit(m, j0, j1, at)
      echo &"    [{a:4.1f},{b:4.1f})   {fw.d0:11.6f}  {fw.ed0:11.6f}  " &
           &"{fw.chi2dof:11.3e}  {fw.dp:9.4f}"

proc delta0(lv, ntv: int, atv: float, conv: GeomConv,
            escale = 1.0): tuple[d, e, chi: float] =
  ## Delta_0 for one (L, L_t) point, same window and estimator as `reportCorr`.
  let
    c = gaugeCorr(lv, ntv, atv, 1.0/g2inv, conv, nsrc, escale)
    half = ntv div 2
  var gh = newSeq[float](half + 1)
  for i in 0..half: gh[i] = c.g[i]
  let
    tt = atv*float(ntv)
    m = effMass(gh, atv, tt)
    fit = plateauFit(m, int(round(fitLo/atv)), min(int(round(fitHi/atv)), m.len), atv)
  (fit.d0, fit.ed0, fit.chi2dof)

if doCorr != 0:
  let g2 = 1.0/g2inv
  echo ""
  echo "================ Job 1: free gauge current correlator (T1.5a, T1.5b) ========"
  for lv in levels(lev):
    for conv in [gcGeodesic, gcExactArea, gcFlat]:
      reportCorr(lv, conv, gaugeCorr(lv, nt, at, g2, conv, nsrc))
  # G_g must not depend on g^2: the 1/g^2 prefactor cancels the g^2 of Mtilde^{-1}.
  block:
    let
      a = gaugeCorr(1, 24, ttot/24.0, 1.0/g2inv, gcGeodesic, 4)
      b = gaugeCorr(1, 24, ttot/24.0, 3.7, gcGeodesic, 4)
    var e = 0.0
    for i in 0..<a.g.len: e = max(e, abs(a.g[i]/b.g[i] - 1.0))
    echo ""
    echo &"    g^2 independence of G_g (L=1, nt=24, g^2 = {1.0/g2inv:.4f} vs 3.7): " &
         &"max relative difference {e:.3e}"
  # How much O(a^2) freedom is there at L = 1?  Every face and every edge is in a
  # single icosahedral orbit there, so A_tri = 4pi/20 and A_e = 4pi/30 exactly under
  # any convention that tiles the sphere, and the whole convention dependence of the
  # action collapses to the single ratio beta_l/beta_tri.
  block:
    echo ""
    echo "    L=1 sensitivity of Delta_0 to beta_l/beta_tri (the only shape parameter"
    echo "    of the action at L=1 besides a_t); geodesic convention is escale = 1"
    echo "      escale    Delta_0     dDelta_0/dln(escale)"
    var prev = (0.0, 0.0)
    for s in [0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2]:
      let d = delta0(1, nt, at, gcGeodesic, s)
      var slope = NaN
      if prev[0] > 0.0: slope = (d.d - prev[1])/ln(s/prev[0])
      echo &"      {s:5.2f}   {d.d:10.6f}   {slope:12.5f}"
      prev = (s, d.d)
    let df = delta0(1, nt, at, gcFlat)
    echo &"      flat convention: Delta_0 = {df.d:.6f}, i.e. the O(abar^2) ambiguity" &
         &" at L=1 spans {100.0*(df.d/1.356697 - 1.0):.1f}%"
  # T = 12 variant (the paper's n_max table uses T = 12 at L_t = 120)
  block:
    let
      atv = 12.0/float(nt)
      c = gaugeCorr(1, nt, atv, 1.0/g2inv, gcGeodesic, nsrc)
      half = nt div 2
    var gh = newSeq[float](half + 1)
    for i in 0..half: gh[i] = c.g[i]
    let
      m = effMass(gh, atv, 12.0)
      fit = plateauFit(m, int(round(3.0/atv)), int(round(6.0/atv)), atv)
    echo ""
    echo &"    T = 12 variant (L=1, nt={nt}, at = {atv:.6f}), window [3, 6): " &
         &"Delta_0 = {fit.d0:.6f}"

# --- job 1c: continuum extrapolation (T1.5c, T1.5d) -------------------------

if doCont != 0:
  echo ""
  echo "================ Job 1c: O(a^2) scaling and Delta_0^cont (T1.5c, T1.5d) ===="
  const ntGrid = [48, 64, 96, 120]        ## paper's L_t grid at T = 16
  for conv in (if contFlat != 0: @[gcGeodesic, gcExactArea, gcFlat] else: @[gcGeodesic]):
    var as2, at2, y, e, lcol, ntcol: seq[float]
    echo ""
    echo &"  geometry = {convTag[conv]}, T = {ttot}, " &
         &"window [{fitLo}, {fitHi})"
    echo "    L   L_t      abar        a_t       Delta_0     dev(sqrt2)   chi2/dof"
    for lv in levels(contLev):
      let sph = newSphere(lv)
      for ntv in ntGrid:
        let
          atv = ttot/float(ntv)
          d = delta0(lv, ntv, atv, conv)
        echo &"  {lv:3d}{ntv:6d}  {sph.abar:10.6f} {atv:10.6f}  {d.d:11.6f}  " &
             &"{100.0*(d.d/sqrt(2.0)-1.0):+9.4f}%  {d.chi:10.3e}"
        as2.add sph.abar*sph.abar
        at2.add atv*atv
        y.add d.d
        e.add 1.0
        lcol.add float(lv)
        ntcol.add float(ntv)
    let pf = contFit2(as2, at2, y, e)
    # unit weights: rescale the parameter errors by sqrt(chi2/dof), the usual
    # unweighted-least-squares convention (there are no statistical errors here).
    let sc = sqrt(pf.chi2/float(pf.dof))
    echo &"    (V.7) Delta_0 = a + c_s abar^2 + c_t a_t^2 ->" &
         &"  a = {pf.a:.6f} +- {pf.ea*sc:.6f}   c_s = {pf.cs:+.5f}   c_t = {pf.ct:+.5f}"
    echo &"    published Delta_0^cont = 1.41409(18)   exact sqrt(2) = {sqrt(2.0):.6f}" &
         &"   dev {100.0*(pf.a/sqrt(2.0)-1.0):+.4f}%   rms residual {sc:.2e}" &
         &"   (dof {pf.dof})"
    writeTsv(outDir / &"scaling_{convTag[conv]}.tsv",
             {"T": &"{ttot:.17g}", "fitLo": &"{fitLo:.17g}", "fitHi": &"{fitHi:.17g}",
              "Delta0cont": &"{pf.a:.17g}", "Delta0conterr": &"{pf.ea*sc:.17g}",
              "cs": &"{pf.cs:.17g}", "ct": &"{pf.ct:.17g}",
              "published": "1.41409(18)"},
             ["L", "Lt", "abar2", "at2", "Delta0"], [lcol, ntcol, as2, at2, y])

# --- job 2: gradient-flow scale scan ----------------------------------------

proc flowScan(lv: int, g2: float, times: openArray[float], sed: int,
              hfac = 1.0): tuple[e, d: seq[float]] =
  ## <E_s(s)> over `nconf` exact-heatbath configurations.  The heatbath is at the
  ## quoted g^2; the flow itself is generated by the g-independent quadratic form
  ## (newBeta at g2 = 1), the standard Luscher normalization in which the flow time
  ## is a length^2.  Flowing with the 1/g^2 of (IV.26) instead only rescales s.
  let
    sph = newSphere(lv)
    lat = newLat(sph, flowNt, flowAt)
    bs = newBeta(lat, g2)
    bf = newBeta(lat, 1.0)
    flowH = 0.1/mDiagMax(lat, bf)   ## explicit RK needs h lambda_max = O(1)
  # `result` cannot be captured by the measurement closure (Nim memory safety),
  # so accumulate into a local and hand it over at the end.
  var
    acc = newSeq[float](times.len)
    acc2 = newSeq[float](times.len)
    r: Threefry4x64
    u = newGauge(lat)
  r.seedIndep(sed, 0)
  for n in 0..<nconf:
    let info = heatbath(lat, u, bs, r, 1e-20, maxits)
    if not info.converged:
      echo &"    WARNING heatbath did not converge, r2 = {info.r2:.3e}"
    var k = 0
    flowRun(lat, u, bf, times, hfac*flowH, RK4CK_2N,
      proc(t: float, v: Gauge) =
        let e = energyDensity(lat, v, bs)
        acc[k] += e
        acc2[k] += e*e
        inc k)
  let fn = float(nconf)
  var err = newSeq[float](times.len)
  for i in 0..<acc.len:
    acc[i] /= fn
    err[i] = sqrt(max(0.0, acc2[i]/fn - acc[i]*acc[i])/(fn - 1.0))
  (acc, err)

if doFlow != 0:
  echo ""
  echo "================ Job 2: gradient-flow scale scan (T2.2, slide 9) ==========="
  var times = newSeq[float](nflow)
  let
    smin = 0.5*flowAt*flowAt
    lo = ln smin
    hi = ln flowMax
  for i in 0..<nflow: times[i] = exp(lo + (hi - lo)*float(i)/float(nflow - 1))
  var
    colS, colE, colD, colY, colX, colL, colG: seq[float]
    sed = seed
    ref1: seq[float]
  for lv in levels(flowLev):
    let sph = newSphere(lv)
    for g2a in [0.5, 1.0, 1.5]:
      let
        g2 = g2a*float(lv)             ## g^2 R = (g^2 a) L, doc/03 Tier 2 grid
        (e, d) = flowScan(lv, g2, times, sed)
      if g2a == 0.5:
        # step-size convergence on the same ensemble at half the RK step, and the
        # reference curve the other two couplings are compared against
        let (h2, _) = flowScan(lv, g2, times, sed, 0.5)
        var de = 0.0
        for i in 0..<e.len: de = max(de, abs(e[i]/h2[i] - 1.0))
        echo &"    RK step convergence at L={lv}: max |E(h)/E(h/2) - 1| = {de:.3e}"
        ref1 = e
      else:
        var mp = 0.0
        for i in 0..<e.len: mp = max(mp, abs(e[i] - ref1[i])/d[i])
        echo &"    L={lv}, g^2 a = {g2a} vs 0.5, independent ensembles: " &
             &"max |dE|/sigma = {mp:.2f}   (the free theory is exactly " &
             &"g^2-independent, so this is a pure statistics check)"
      inc sed
      echo ""
      echo &"--- L = {lv}  abar = {sph.abar:.6f}  g^2 a = {g2a}  g^2 R = {g2:.3f}" &
           &"  nt = {flowNt}  at = {flowAt}  nconf = {nconf}"
      echo "        s        r/s       E_s(s)       err       E_s sqrt(L)  dlnE/dln(r/s)"
      for i in 0..<nflow:
        var sl = NaN
        if i > 0 and i < nflow-1:
          sl = ln(e[i+1]/e[i-1])/ln(times[i-1]/times[i+1])
        echo &"  {times[i]:9.5f}  {1.0/times[i]:9.4f}  {e[i]:12.6e} {d[i]:10.3e}  " &
             &"{e[i]*sqrt(float(lv)):12.6e}  {sl:12.5f}"
        colS.add times[i]
        colX.add 1.0/times[i]
        colE.add e[i]
        colD.add d[i]
        colY.add e[i]*sqrt(float(lv))
        colL.add float(lv)
        colG.add g2
  writeTsv(outDir / &"flowscan_nt{flowNt}.tsv",
           {"at": &"{flowAt:.17g}", "nt": $flowNt, "nconf": $nconf,
            "r": "1.0 (sphere radius R)",
            "flow": "generated by newBeta(g2=1), i.e. coupling-independent flow time",
            "Es": "spatial action density S_s/(4 pi T)"},
           ["s", "r_over_s", "E_s", "E_s_err", "E_s_sqrtL", "L", "g2R"],
           [colS, colX, colE, colD, colY, colL, colG])

echo ""
echo &"wrote {outDir}"

processSaveParams()
writeParamFile()
qexFinalize()
