## the shift-invert spectrum pipeline of spectrum.nim (hcspec)
## on 4^4 lattices:
##
##   1. FREE FIELD, both lattices: unit gauge + antiperiodic time, massless
##      r=1 cSW=1 clover operator (the clover term vanishes at F = 0), the
##      shift-invert Arnoldi eigenvalues lambda = sigma + 1/mu must match the
##      exact free spectra: hcfree.freeD8 (honeycomb, 8x8 blocks at the cell
##      momenta, k3 -> k3 + pi/Nt) and the closed cubic form
##      m + sum(1-cos p) +- i|sin p|.
##   2. ROUGH CONFIG (honeycomb, warm + 2 stout steps + setBC): all pairs
##      converge; the converged set is conjugation symmetric
##      (gamma5-hermiticity); complex modes have chirality ~ 0 (biorthogonality
##      forces <v|g5|v> = 0 exactly for Im lambda != 0 -- the pair partners'
##      chiralities are trivially equal); direct residuals are small.
##   3. DETERMINISM: the same run repeated gives bit-identical eigenvalues and
##      chiralities (serial driver ops + deterministic start vectors).
##   4. ATIYAH-SINGER constant-flux config (tclover/ttopo construction),
##      both lattices: the hexagon/1x1-clover charge is Q ~= 2 n1 n2, the
##      operator has 2|n1 n2| near-zero REAL modes whose chiralities all have
##      one sign, complex modes have |chi| ~ 0, and
##      qDiracSign*(n+ - n-) = round(Q) -- this PINS the Q_Dirac sign
##      convention used by both drivers.
##   5. mode bookkeeping (measureModes, summarize, specLines) on hand-built
##      candidates without an eigensolve.

import std/[math, complex, sequtils, strformat, strutils, unittest]
import qex except epsilon
import physics/qcdTypes
import gauge
import ../hcspec
import ../hcflow
import ../cubic
import ../hcfree
import helpers

qexInit()
doAssert nRanks == 1, "tspectrum is a single-rank test"

# ---------------------------------------------------------------------------
# honeycomb pipeline pieces (4^4 cells)
# ---------------------------------------------------------------------------

let hl = newLayout([4, 4, 4, 4])
let lo = hl
var hg = newHcGauge(hl)
var cw = newHcWilson(hg, 1.0)
var proto = newHcFermion(hl)
type HF = typeof(proto)
var wt = newTopoWork(hg)

proc hcStartVec(count: ref uint64): proc (v: var HF) =
  result = proc (v: var HF) =
    inc count[]
    let salt = sm64(count[])*0x10000'u64
    for i in lo.sites:
      for sp in 0..3:
        for c in 0..2:
          let k = salt + uint64(i)*48 + uint64(sp)*12 + uint64(c)*4
          setC(v.a{i}[sp][c], u01(k), u01(k+1))
          setC(v.b{i}[sp][c], u01(k+2), u01(k+3))

proc hcEigs(sigma: float; nev, ncvv: int; innerR2: float):
    tuple[modes: seq[SpecMode], napply: int] =
  ## shift-invert eigensolve on the current contents of hg (via cw)
  var stats = SiStats()
  var cnt = new uint64
  var op = newShiftInvertOp[HF](
    applyM = proc (rr: var HF; x: HF) = cw.D(rr, x, -sigma, 1.0),
    applyMdag = proc (rr: var HF; x: HF) = cw.Ddag(rr, x, -sigma, 1.0),
    newVec = proc (): HF = newOneOf(proto),
    startVec = hcStartVec(cnt),
    r2req = innerR2, maxits = 4000, stats = stats)
  let (mus, vecs, _, napply) = arnoldi(op, nev, ncvv, 1e-9, 60, "LM", 0)
  result.napply = napply
  result.modes = measureModes(mus, vecs, sigma,
    proc(r: var HF; x: HF) = cw.D(r, x, 0.0, 1.0))

proc cubicEigs(ccw: auto; clo: Layout; sigma: float; nev, ncvv: int;
               innerR2: float; saltMul: uint64; residcut = -1.0): seq[SpecMode] =
  ## the same on a cubic CubicWilson operator
  var cproto = clo.DiracFermion()
  cproto := 0
  type DF = typeof(cproto)
  var stats = SiStats()
  var cnt = 0'u64
  var op = newShiftInvertOp[DF](
    applyM = proc (rr: var DF; x: DF) = ccw.D(rr, x, -sigma),
    applyMdag = proc (rr: var DF; x: DF) = ccw.Ddag(rr, x, -sigma),
    newVec = proc (): DF =
      result = newOneOf(cproto)
      result := 0,
    startVec = proc (v: var DF) =
      inc cnt
      let salt = sm64(cnt)*saltMul
      for i in clo.sites:
        for sp in 0..3:
          for c in 0..2:
            let k = salt + uint64(i)*24 + uint64(sp)*6 + uint64(c)*2
            setC(v{i}[sp][c], u01(k), u01(k+1)),
    r2req = innerR2, maxits = 4000, stats = stats)
  let (mus, vecs, _, _) = arnoldi(op, nev, ncvv, 1e-9, 60, "LM", 0)
  measureModes(mus, vecs, sigma, proc(r: var DF; x: DF) = ccw.D(r, x, 0.0),
               residcut)

proc lams(modes: seq[SpecMode]): seq[Complex64] = modes.mapIt(it.lam)

proc realModes(modes: seq[SpecMode]): tuple[minChi, maxRe, maxChiC: float] =
  ## min |chi| and max |Re lam| over the real modes (0 if none: the domain
  ## bounds), max |chi| over the complex modes; echoes the real modes
  var chis, res: seq[float]
  for m in modes:
    if abs(m.lam.im) < 1e-6:
      chis.add abs(m.chi)
      res.add abs(m.lam.re)
      echo &"  real mode: lam = {m.lam.re:.6f}  chi = {m.chi:.6f}"
    else:
      result.maxChiC = max(result.maxChiC, abs(m.chi))
  if chis.len > 0:
    result.minChi = min(chis)
    result.maxRe = max(res)

template checkModes(vecs: untyped) =
  ## measureModes / summarize on hand-built candidates with D = 2 I; the last
  ## candidate has lambda = 5 + 4i, hence residual 5
  type F = typeof(vecs[0])
  let sig = -0.25
  let mu = complex64(1.0/(2.0 - sig), 0.0)
  let bad = complex64(5.0, 4.0)
  let mus = @[mu, mu, mu, complex64(1.0, 0.0)/(bad - complex64(sig, 0.0))]
  let vs = @[vecs[0], vecs[1], vecs[2], vecs[2]]
  let n2 = vnorm2(vecs[2])
  var calls = 0
  let applyD = proc(r: var F; x: F) =
    inc calls
    copyP(r, x)
    scaleP(r, 2.0)
  let modes = measureModes(mus, vs, sig, applyD)
  check modes.len == 4
  for i in 0..2:
    check abs(modes[i].lam - complex64(2.0, 0.0)) < 1e-13
    check modes[i].resid < 1e-13
  check abs(modes[0].chi - 1.0) < 1e-13
  check abs(modes[1].chi + 1.0) < 1e-13
  check abs(modes[2].chi - 0.6) < 1e-13
  check abs(modes[3].lam - bad) < 1e-13
  check abs(modes[3].resid - 5.0) < 1e-13
  check vnorm2(vecs[2]) == n2
  let kept = measureModes(mus, vs, sig, applyD, 0.5)
  check kept.len == 3
  let s = summarize(kept, 1e-6, sig, 100.0)
  check s.nconv == 3
  check s.nreal == 3
  check s.nplus == 2
  check s.nminus == 1
  check s.qdirac == -1.0
  check abs(s.sumChiReal - 0.6) < 1e-13
  # lambda = 1 gives residual exactly 1, including the cutoff boundary
  let edge = measureModes(@[complex64(1.0, 0.0)], @[vecs[0]], 0.0, applyD, 1.0)
  check edge.len == 1
  check edge[0].resid == 1.0
  let n = calls
  check measureModes(newSeq[Complex64](), newSeq[F](), sig, applyD).len == 0
  check calls == n

var modesA: seq[SpecMode]    # test 2 result, reused by test 3

suite "spectrum pipeline (4^4)":

  test "1a. free honeycomb clover operator vs hcfree.freeD8 (anti time)":
    hg.unit
    threads:
      hg.setBC
    cw.gaugeRefresh
    # exact spectrum: freeD8 at the cell momenta with k3 -> k3 + pi/Nt
    var exact: seq[Complex64]
    for k in cellMomenta(4, 4, antiperiodicTime = true):
      let d8 = freeD8(k, 1.0, 0.0)
      var z = newZMat(8, 8)
      for a in 0..<8:
        for b in 0..<8:
          z[a, b] = d8[a][b]
      for e in zeig(z).w: exact.add e
    echo &"  {exact.len} exact eigenvalues (x3 color)"
    let sigma = -0.3
    let (modes, napply) = hcEigs(sigma, 12, 48, 1e-18)
    let s = summarize(modes, 1e-6, sigma, 100.0)
    let wn = worstNearest(modes.lams, exact)
    # is the closest-to-sigma exact shell captured?
    var lmin = exact[0]
    for e in exact:
      if abs(e - complex64(sigma, 0.0)) < abs(lmin - complex64(sigma, 0.0)):
        lmin = e
    let dcap = worstNearest([lmin], modes.lams)
    echo &"  {modes.len} converged, {napply} applies, worst direct resid {s.worstResid:.2e}"
    echo &"  worst |lambda - exact| = {wn:.3e}; closest-to-sigma exact shell missed by {dcap:.3e}"
    ok(&"free hc: 12 pairs converged", modes.len == 12)
    ok(&"free hc: direct residuals < 1e-7 (worst {s.worstResid:.1e})", s.worstResid < 1e-7)
    ok(&"free hc: eigenvalues match freeD8 to 1e-7 (worst {wn:.1e})", wn < 1e-7)
    ok(&"free hc: the closest-to-sigma shell is captured ({dcap:.1e})", dcap < 1e-7)
    ok("free hc: conj pairing < 1e-7", conjGap(modes.lams) < 1e-7)

  test "1b. free cubic clover operator vs closed form (anti time)":
    let clat = @[4, 4, 4, 4]
    let clo = clat.newLayout
    var cg = clo.newGauge
    for mu in 0..<4: cg[mu] := 1
    cg.setBC
    var ccw = newCubicWilson(cg, 1.0)
    var exact: seq[Complex64]
    block:
      var kco = @[0, 0, 0, 0]
      while true:
        var mm = 0.0
        var k2 = 0.0
        for mu in 0..<4:
          let off = if mu == 3: 0.5 else: 0.0
          let p = 2.0*PI*(float(kco[mu]) + off)/float(clat[mu])
          mm += 1.0 - cos(p)
          k2 += sin(p)*sin(p)
        exact.add complex64(mm, sqrt(k2))
        exact.add complex64(mm, -sqrt(k2))
        var mu = 0
        while mu < 4:
          inc kco[mu]
          if kco[mu] < clat[mu]: break
          kco[mu] = 0
          inc mu
        if mu == 4: break
    let sigma = -0.3
    let modes = cubicEigs(ccw, clo, sigma, 12, 48, 1e-18, 0x10000'u64)
    let s = summarize(modes, 1e-6, sigma, 100.0)
    let wn = worstNearest(modes.lams, exact)
    echo &"  {modes.len} converged, worst direct resid {s.worstResid:.2e}"
    echo &"  worst |lambda - exact| = {wn:.3e}"
    ok("free cubic: 12 pairs converged", modes.len == 12)
    ok(&"free cubic: direct residuals < 1e-7 (worst {s.worstResid:.1e})", s.worstResid < 1e-7)
    ok(&"free cubic: eigenvalues match closed form to 1e-7 (worst {wn:.1e})", wn < 1e-7)

  test "2. rough config (warm + 2 stout + setBC): pairing, chirality, resids":
    var r = lo.newRNGField(MRG32k3a, 13579'u64)
    threads:
      hg.warm(0.55, r)
      hg.reunit
    var sw = newActionWork(hg)
    var sc = newOneOf(hg)
    stout(sw, hg, 0.05, sc, 2)
    threads:
      hg.setBC
    cw.gaugeRefresh
    let (modes, _) = hcEigs(-0.45, 12, 48, 1e-16)
    modesA = modes
    let s = summarize(modes, 1e-6, -0.45, 100.0)
    let cp = conjGap(modes.lams)
    var maxChiC = 0.0        # complex modes: chi must vanish
    for m in modes:
      if abs(m.lam.im) >= 1e-6:
        maxChiC = max(maxChiC, abs(m.chi))
    echo &"  {modes.len} converged, worst direct resid {s.worstResid:.2e}, conj pairing {cp:.2e}"
    echo &"  {s.nreal} real modes; max |chi| over complex modes = {maxChiC:.2e}"
    ok("rough hc: all 12 pairs converged", modes.len == 12)
    ok(&"rough hc: direct residuals < 1e-6 (worst {s.worstResid:.1e})", s.worstResid < 1e-6)
    ok(&"rough hc: converged set conjugation symmetric ({cp:.1e})", cp < 1e-6)
    ok(&"rough hc: complex-pair chirality ~ 0, both partners (max {maxChiC:.1e})",
       maxChiC < 1e-4)

  test "3. determinism: identical rerun (same seeds, same thread count)":
    let (modes, _) = hcEigs(-0.45, 12, 48, 1e-16)    # hg unchanged: same operator
    var dv = 0.0
    var dc = 0.0
    for i, m in modes:
      dv = max(dv, abs(m.lam - modesA[i].lam))
      dc = max(dc, abs(m.chi - modesA[i].chi))
    echo &"  max |dlambda| = {dv:.3e}, max |dchi| = {dc:.3e}"
    ok("determinism: eigenvalues bit-identical", dv == 0.0)
    ok("determinism: chiralities bit-identical", dc == 0.0)

  test "4a. honeycomb flux config (n1 = n2 = 1): Q_Dirac vs hexagon-clover Q":
    setFluxHc(hg, 1, 1)
    let qF = EQ(wt, hg).q
    # Antiperiodic time, as in production.  Without it the color-charge-0
    # component (T = diag(1,-1,0)) is a FREE fermion whose p = 0 modes sit at
    # exactly lambda = 0 with arbitrary chirality mixing (observed); the index
    # zero modes of the charge +-1 components exist for either spin structure.
    threads:
      hg.setBC
    cw.gaugeRefresh
    echo &"  hexagon-clover Q = {qF:.6f} (exact 2 n1 n2 = 2, artifact 1/L^4)"
    let (modes, _) = hcEigs(-0.25, 10, 40, 1e-16)
    let s = summarize(modes, 1e-6, -0.25, 100.0)
    let rm = realModes(modes)
    echo &"  nreal = {s.nreal}, n+ = {s.nplus}, n- = {s.nminus}, sum chi = {s.sumChiReal:.4f}"
    echo &"  Q_Dirac = qDiracSign*(n+ - n-) = {s.qdirac:.1f}   vs   round(Q_flow-style) = {round(qF):.1f}"
    ok("flux hc: exactly 2 real near-zero modes in the window", s.nreal == 2)
    ok(&"flux hc: real modes near zero (max |Re| {rm.maxRe:.2e})", rm.maxRe < 0.05)
    ok(&"flux hc: real-mode chirality one sign, |chi| > 0.99 (min {rm.minChi:.4f})",
       rm.minChi > 0.99 and (s.nplus == 0 or s.nminus == 0))
    ok(&"flux hc: complex modes chi ~ 0 (max {rm.maxChiC:.1e})", rm.maxChiC < 1e-4)
    ok(&"flux hc: Q_Dirac == round(Q) == 2 (PINS qDiracSign = {qDiracSign})",
       abs(s.qdirac - round(qF)) < 0.5)
    ok(&"flux hc: sum of real-mode chiralities integer-ish " &
       &"(|{s.sumChiReal:.3f}| vs 2, within 0.3)", abs(abs(s.sumChiReal) - 2.0) < 0.3)

  test "4b. cubic flux config (n1 = n2 = 1): same counting":
    let clat = @[4, 4, 4, 4]
    let clo = clat.newLayout
    var cg = clo.newGauge
    setFluxCubic(cg, clo, 1, 1)
    let f = cg.fmunu 1
    let qF = f.topoQ
    echo &"  1x1-clover Q = {qF:.6f} (exact 2, artifact O(f^2))"
    cg.setBC                    # antiperiodic time (see 4a)
    var ccw = newCubicWilson(cg, 1.0)
    let sigma = -0.25
    let modes = cubicEigs(ccw, clo, sigma, 10, 40, 1e-16, 0x30000'u64, 1e-7)
    let s = summarize(modes, 1e-6, sigma, 100.0)
    let rm = realModes(modes)
    echo &"  nreal = {s.nreal}, n+ = {s.nplus}, n- = {s.nminus}, sum chi = {s.sumChiReal:.4f}, Q_Dirac = {s.qdirac:.1f}"
    ok("flux cubic: 2 near-zero real modes found", s.nreal == 2)
    ok(&"flux cubic: one-sign chirality, |chi| > 0.9 (min {rm.minChi:.4f})",
       rm.minChi > 0.9 and (s.nplus == 0 or s.nminus == 0))
    ok(&"flux cubic: Q_Dirac == round(Q) == 2 (same sign as honeycomb)",
       abs(s.qdirac - round(qF)) < 0.5)

  test "5. mode bookkeeping without an eigensolve":
    block:                      # honeycomb candidates
      var vs = @[newHcFermion(hl), newHcFermion(hl), newHcFermion(hl)]
      # gamma5 = diag(1, 1, -1, -1); mixed amplitudes 2 and 1 give chi = 3/5
      for i in hl.sites:
        setC(vs[0].a{i}[0][0], 3.0, 0.0)
        setC(vs[0].b{i}[0][0], 3.0, 0.0)
        setC(vs[1].a{i}[2][0], 2.0, 0.0)
        setC(vs[1].b{i}[2][0], 2.0, 0.0)
        setC(vs[2].a{i}[0][0], 2.0, 0.0)
        setC(vs[2].b{i}[0][0], 2.0, 0.0)
        setC(vs[2].a{i}[2][0], 1.0, 0.0)
        setC(vs[2].b{i}[2][0], 1.0, 0.0)
      checkModes(vs)
    block:                      # cubic candidates
      let clo = newLayout(@[4, 4, 4, 4])
      var vs = @[clo.DiracFermion(), clo.DiracFermion(), clo.DiracFermion()]
      for v in vs: v := 0
      for i in clo.sites:
        setC(vs[0]{i}[0][0], 3.0, 0.0)
        setC(vs[1]{i}[2][0], 2.0, 0.0)
        setC(vs[2]{i}[0][0], 2.0, 0.0)
        setC(vs[2]{i}[2][0], 1.0, 0.0)
      checkModes(vs)
    block:                      # report fields
      let modes = @[
        SpecMode(lam: complex64(0.5, 0.0), chi: -0.9, resid: 2e-5),
        SpecMode(lam: complex64(0.6, 0.8), chi: 0.0, resid: 3e-6),
        SpecMode(lam: complex64(0.6, -0.8), chi: 0.0, resid: 4e-6)]
      let st = SiStats(totIts: 96, maxUsedIts: 20, nHitMax: 1, worstDirect: 3e-5)
      let rows = specLines(modes, 1e-6, 0.0, 100.0, 7, 0.05, 0.75, 0.7, 2, 12, st, 1.25)
      check rows.len == 5
      for i in 0..2:
        let cols = rows[i].splitWhitespace
        check cols.len == 8
        check cols[0..3] == @["EIG", "7", "0.050", $i]
        check cols[4].parseFloat == modes[i].lam.re
        check cols[5].parseFloat == modes[i].lam.im
        check cols[6].parseFloat == modes[i].chi
        check cols[7].parseFloat == modes[i].resid
      let cfg = rows[3].splitWhitespace
      check cfg.len == 15
      check cfg[0..2] == @["CFG", "7", "0.050"]
      check cfg[3].parseFloat == 0.5
      check cfg[4].parseFloat == 0.0
      check cfg[5..^1] == @["3", "1", "0", "1", "1.0", "0.750000", "0.700000", "12", "96", "1.25"]
      let ext = rows[4].splitWhitespace
      check ext.len == 27
      check ext[0..2] == @["CFGX", "7", "0.050"]
      let keys = ["reach", "recut", "realabove", "conj", "worstresid", "imgapC", "imgapR", "sumchi", "nbad", "cgmax", "cghitmax", "sidirect"]
      let vals = [1.0, 100.0, 0.0, 0.0, 2e-5, 0.8, 0.0, -0.9, 2.0, 20.0, 1.0, 3e-5]
      for i, key in keys:
        check ext[3 + 2*i] == key
        check ext[4 + 2*i].parseFloat == vals[i]
      ok("specLines: EIG/CFG/CFGX rows carry the mode and summary fields", rows.len == 5)

qexFinalize()
