## WP-K spectra figures (Tier 1, arXiv:2510.03085): Fig. 4 (free D_W spectra),
## Fig. 5 (the flat-limit formula (IV.8)), Fig. 6 (D_W vs D_ov and the generalized
## eigenvalues), plus the slide-8 free-legend table min|D_W - 1| under both kappa
## conventions.
##
## Everything is free-field, so the spectra are computed per antiperiodic Matsubara
## mode with the machinery shared from rfree.nim (spatialDense/modeBlock/ovFromDw,
## the construction twilson pins against the full dense operator to 5.8e-14).
## Eigenvalues per mode via zgeev (eigens/linalgFuncs.zgeigs).
##
## Fig. 6 conventions: raw D_lat of (IV.1), M = 1, and the generalized eigenvalues
## use the CORRECTED volume weight diag(volbar/volw) of doc/02 section 3.2 (i.e.
## eigenvalues of diag(Abar_y/A_y) D_lat; the weight the paper prints is inverted,
## WP-E's correction).
##
## Slide-8 legends (see doc/06 "THE COUPLING CONVENTION" and resolved open question
## 2): the published values 1.154 / 1.010 / 0.965 match the FLAT kappa
## (l*_1 + l*_2)/abar, at a_t = 0.2 for L = 2, 4 and a_t = 16/120 for L = 1; the
## tree-wide production convention is the exact kite area.  Both are emitted, exact
## as primary.
##
## Output: TSV under output/radial/free/.

import base
import std/[math, os, strformat, times]
import core/[analytic, lattice]
import meas/dataio
import rfree

qexInit()
freezeTimers()

letParam:
  doFig4 = 1
  doFig5 = 1
  doFig6 = 1
  doLegend = 1

installStandardParams()
echoParams()
processHelpParam()

const outDir = currentSourcePath().parentDir.parentDir.parentDir.parentDir /
               "output" / "radial" / "free"
createDir outDir

# --- Fig. 4: free D_W spectra, T = 4 -----------------------------------------

const fig4Cases = [(1, 24), (2, 24), (4, 24), (2, 16), (2, 48)]

if doFig4 != 0:
  echo ""
  echo "================ Fig. 4: free D_W spectra (T = 4, raw D_lat units) =========="
  for (lv, lt) in fig4Cases:
    let
      t0 = epochTime()
      sph = newSphere(lv)
      at = 4.0/float(lt)
      (l, dsp) = spatialDense(sph, at)
      nd = 2*sph.nv
    var
      re, im, idx: seq[float]
      x: seq[Complex64]
      reMin = 1e300
      reMax = -1e300
      d1 = 1e300
    for n in 0..<lt:
      modeBlockInto(l, dsp, PI*float(2*n + 1)/float(lt), x)
      let ev = eigvals(x, nd)
      for z in ev:
        re.add z.re
        im.add z.im
        idx.add float(n)
        reMin = min(reMin, z.re)
        reMax = max(reMax, z.re)
        d1 = min(d1, abs(z - complex64(1.0, 0.0)))
    writeTsv(outDir / &"fig4_L{lv}_Lt{lt}.tsv",
             {"source": "arXiv:2510.03085 Fig. 4", "lattice": &"L{lv}",
              "Lt": $lt, "T": "4", "at": &"{at:.17g}", "abar": &"{sph.abar:.17g}",
              "normalization": "raw D_lat of (IV.1), exact-kappa convention"},
             ["re", "im", "matsubara"], [re, im, idx])
    echo &"  L={lv} Lt={lt:2d}: {re.len} eigenvalues, Re in [{reMin:.6f}, " &
         &"{reMax:.6f}], min|lam-1| = {d1:.6f}   ({epochTime() - t0:.1f} s)"

# --- Fig. 5: flat-limit spectrum (IV.8) ---------------------------------------

if doFig5 != 0:
  echo ""
  echo "================ Fig. 5: flat spectrum (IV.8), kappa = 1/sqrt3 =============="
  for (lv, lt) in fig4Cases:
    let
      sph = newSphere(lv)
      at = 4.0/float(lt)
      kap = 1.0/sqrt(3.0)
      kapT = 0.5*sqrt(3.0)*sph.abar/at
      fl = flatSpectrum(kap, kapT, 16, 16, lt)
    var re, im, ktc: seq[float]
    for i, z in fl:
      re.add z.re
      im.add z.im
      ktc.add float((i div 2) mod lt)
    writeTsv(outDir / &"fig5_flat_L{lv}_Lt{lt}.tsv",
             {"source": "arXiv:2510.03085 Fig. 5", "matches": &"L{lv}_Lt{lt}",
              "kappa": &"{kap:.17g}", "kappaT": &"{kapT:.17g}",
              "grid": "16 x 16 x Lt, both branches"},
             ["re", "im", "kt"], [re, im, ktc])
    echo &"  L={lv} Lt={lt:2d}: kappa' = {kapT:.6f}, {fl.len} eigenvalues"

# --- Fig. 6: D_W vs D_ov and the generalized eigenvalues ----------------------

if doFig6 != 0:
  echo ""
  echo "================ Fig. 6: D_W vs D_ov, T = 4, L = 4, Lt = 24, M = 1 =========="
  let
    t0 = epochTime()
    sph = newSphere(4)
    lt = 24
    at = 4.0/float(lt)
    (l, dsp) = spatialDense(sph, at)
    nd = 2*sph.nv
  # corrected volume weight (doc/02 section 3.2): rows scaled by volbar/volw
  var wv = newSeq[float](nd)
  for i in 0..<nd: wv[i] = l.volbar/l.volw[i div 2]
  proc rowScale(a: seq[Complex64]): seq[Complex64] =
    result = a
    for j in 0..<nd:
      for i in 0..<nd: result[i + nd*j] = wv[i]*result[i + nd*j]
  var
    wRe, wIm, wIdx, oRe, oIm, oIdx: seq[float]
    gwRe, gwIm, gwIdx, goRe, goIm, goIdx: seq[float]
    x: seq[Complex64]
    circDev = 0.0
    d1w = 1e300
  for n in 0..<lt:
    let k = PI*float(2*n + 1)/float(lt)
    modeBlockInto(l, dsp, k, x)
    let dov = ovFromDw(x, nd)          # D_ov(k), M = 1
    var a = x
    for z in eigvals(a, nd):
      wRe.add z.re
      wIm.add z.im
      wIdx.add float(n)
      d1w = min(d1w, abs(z - complex64(1.0, 0.0)))
    a = dov
    for z in eigvals(a, nd):
      oRe.add z.re
      oIm.add z.im
      oIdx.add float(n)
      circDev = max(circDev, abs(abs(z - complex64(1.0, 0.0)) - 1.0))
    a = rowScale(x)
    for z in eigvals(a, nd):
      gwRe.add z.re
      gwIm.add z.im
      gwIdx.add float(n)
    a = rowScale(dov)
    for z in eigvals(a, nd):
      goRe.add z.re
      goIm.add z.im
      goIdx.add float(n)
  let hdr = {"source": "arXiv:2510.03085 Fig. 6", "lattice": "L4", "Lt": $lt,
             "T": "4", "at": &"{at:.17g}", "M": "1",
             "weight": "gen files: rows scaled by volbar/volw (doc/02 3.2 corrected)"}
  writeTsv(outDir / "fig6_dw.tsv", hdr, ["re", "im", "matsubara"],
           [wRe, wIm, wIdx])
  writeTsv(outDir / "fig6_dov.tsv", hdr, ["re", "im", "matsubara"],
           [oRe, oIm, oIdx])
  writeTsv(outDir / "fig6_dw_gen.tsv", hdr, ["re", "im", "matsubara"],
           [gwRe, gwIm, gwIdx])
  writeTsv(outDir / "fig6_dov_gen.tsv", hdr, ["re", "im", "matsubara"],
           [goRe, goIm, goIdx])
  echo &"  {wRe.len} eigenvalues per operator; min|D_W - 1| = {d1w:.6f}; " &
       &"D_ov Ginsparg-Wilson circle max ||lam-1|-1| = {circDev:.3e}   " &
       &"({epochTime() - t0:.1f} s)"

# --- slide-8 free legend: min|D_W - 1| under both kappa conventions -----------

if doLegend != 0:
  echo ""
  echo "================ slide-8 legend: min|D_W - 1|, raw D_lat, M = 1 ============="
  # (L, at, Lt, published legend or 0).  The a_t = 16/120 rows are the L=1
  # sub-question of doc/06 open question 2: the published 1.154 only reproduces
  # near the free-limit paper's temporal spacing, and min|D_W - 1| also depends
  # on Lt (the k-grid), so both Lt = 60 (the oracle's scan) and 120 are listed.
  const cases = [(1, 0.2, 60, 0.0), (2, 0.2, 60, 1.010), (4, 0.2, 60, 0.965),
                 (1, 16.0/120.0, 60, 1.154), (1, 16.0/120.0, 120, 1.154)]
  # pinned reference values from doc/06 (toverlap production row and the a_t scan)
  echo "  known pins: exact L=1 at=0.2 -> 1.1770 (toverlap); flat L=1 at=0.2 ->"
  echo "  1.1234 (WP-E); flat L=2 -> 1.0105; flat L=4 -> 0.9643; flat L=1"
  echo "  at=0.1333 -> 1.1582 (main's oracle).  Published legends: 1.154/1.010/0.965."
  var cl, cat, clt, cconv, cmin: seq[float]
  echo "    L      a_t   Lt   conv    min|D_W-1|   published"
  for (lv, at, lt, pub) in cases:
    let sph = newSphere(lv)
    for flat in [false, true]:
      let
        (l, dsp) = spatialDense(sph, at, flatKap = flat)
        nd = 2*sph.nv
      var
        x: seq[Complex64]
        d1 = 1e300
      for n in 0..<lt:
        modeBlockInto(l, dsp, PI*float(2*n + 1)/float(lt), x)
        for z in eigvals(x, nd):
          d1 = min(d1, abs(z - complex64(1.0, 0.0)))
      let tag = if flat: "flat " else: "exact"
      let ps = if pub > 0.0 and flat: &"{pub:.3f}" else: "-"
      echo &"  {lv:3d}  {at:7.4f}  {lt:3d}   {tag}  {d1:12.6f}   {ps}"
      cl.add float(lv)
      cat.add at
      clt.add float(lt)
      cconv.add (if flat: 1.0 else: 0.0)
      cmin.add d1
  writeTsv(outDir / "slide8_legend.tsv",
           {"source": "slides p.8 legend, free column",
            "quantity": "min |eig(D_W) - 1| over all Matsubara modes, raw D_lat",
            "conv": "0 = exact kite area (production), 1 = flat (l*1+l*2)/abar",
            "published": "L=1: 1.154 (at approx 0.1333), L=2: 1.010, L=4: 0.965"},
           ["L", "at", "Lt", "conv", "minD1"], [cl, cat, clt, cconv, cmin])

echo ""
echo &"wrote {outDir}"

processSaveParams()
writeParamFile()
qexFinalize()
