## Task D4 test: the shift-invert spectrum pipeline of hcSpectrum.nim /
## cubicSpectrum.nim, on 4^4 lattices (fast, < 2 min):
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
##
##   export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
##   cd build_mac && make src/experimental/honeycomb/tests/tspectrum.nim
##   OMP_NUM_THREADS=4 ./bin/tspectrum

import std/[math, complex, sequtils, strformat]
import qex except epsilon
import physics/qcdTypes
import gauge
import ../hcSpectrum
import ../hcfree

qexInit()
doAssert nRanks == 1, "tspectrum is a single-rank test"

var nFail = 0
proc ok(msg: string; cond: bool) =
  if cond:
    echo "PASS: ", msg
  else:
    echo "FAIL: ", msg
    inc nFail

template section(msg: string) =
  echo ""
  echo "=== ", msg, " ==="

proc worstNearest(modes: openArray[SpecMode]; exact: openArray[Complex64]): float =
  for m in modes:
    var d = 1e300
    for e in exact: d = min(d, abs(m.lam - e))
    result = max(result, d)

# ---------------------------------------------------------------------------
# constant-flux Abelian backgrounds (tclover/ttopo/refCubicMeas recipes)
# ---------------------------------------------------------------------------

template setPhase(m: untyped, th: float) =
  ## m := exp(i th T), T = diag(1,-1,0); assumes m starts as the identity
  m[0, 0].re := cos(th)
  m[0, 0].im := sin(th)
  m[1, 1].re := cos(th)
  m[1, 1].im := -sin(th)

proc fluxPhase(x, n: array[4, float]; f1, f2: float; l0, l2: int): float =
  ## exact line integral of A_1 = f1 x0, A_3 = f2 x2 from x to x+n, with the
  ## transition-function corrections at the x0 = L0 / x2 = L2 seams
  result = n[1]*f1*(x[0] + 0.5*n[0]) + n[3]*f2*(x[2] + 0.5*n[2])
  if x[0] + n[0] >= float(l0) - 1e-9:
    result -= f1*float(l0)*(x[1] + n[1])
  if x[2] + n[2] >= float(l2) - 1e-9:
    result -= f2*float(l2)*(x[3] + n[3])

proc setFluxHc(hg: auto, n1, n2: int) =
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

proc setFluxCubic(g: auto, lo: auto, n1, n2: int) =
  let lat = lo.physGeom
  let
    f1 = 2.0*PI*float(n1)/float(lat[0]*lat[1])
    f2 = 2.0*PI*float(n2)/float(lat[2]*lat[3])
  g.unit
  for i in lo.sites:
    let
      x0 = lo.coords[0][i]
      x1 = lo.coords[1][i]
      x2 = lo.coords[2][i]
      x3 = lo.coords[3][i]
    setPhase(g[1]{i}, f1*x0.float)
    if x0 == lat[0]-1:
      setPhase(g[0]{i}, -f1*float(lat[0]*x1))
    setPhase(g[3]{i}, f2*x2.float)
    if x2 == lat[2]-1:
      setPhase(g[2]{i}, -f2*float(lat[2]*x3))

# ---------------------------------------------------------------------------
# honeycomb pipeline pieces (4^4 cells)
# ---------------------------------------------------------------------------

let hl = newHcLayout([4, 4, 4, 4])
let lo = hl.lo
var hg = newHcGauge(hl)
var cw = newHcCloverWilson(hg, 1.0)
var proto = newHcFermion(hl)
type HF = typeof(proto)
var wt = newHcTopoWork(hg)

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

# ---------------------------------------------------------------------------
# 1a. free honeycomb
# ---------------------------------------------------------------------------

section "1a. free honeycomb clover operator vs hcfree.freeD8 (4^4, anti time)"
block:
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
  let s = summarize(modes, 1e-6)
  let wn = worstNearest(modes, exact)
  # is the closest-to-sigma exact shell captured?
  var dmin = 1e300
  var lmin: Complex64
  for e in exact:
    if abs(e - complex64(sigma, 0.0)) < dmin:
      dmin = abs(e - complex64(sigma, 0.0))
      lmin = e
  var dcap = 1e300
  for m in modes: dcap = min(dcap, abs(m.lam - lmin))
  echo &"  {modes.len} converged, {napply} applies, worst direct resid {s.worstResid:.2e}"
  echo &"  worst |lambda - exact| = {wn:.3e}; closest-to-sigma exact shell missed by {dcap:.3e}"
  ok(&"free hc: 12 pairs converged", modes.len == 12)
  ok(&"free hc: direct residuals < 1e-7 (worst {s.worstResid:.1e})", s.worstResid < 1e-7)
  ok(&"free hc: eigenvalues match freeD8 to 1e-7 (worst {wn:.1e})", wn < 1e-7)
  ok(&"free hc: the closest-to-sigma shell is captured ({dcap:.1e})", dcap < 1e-7)
  ok("free hc: conj pairing < 1e-7", worstConjPairing(modes.mapIt(it.lam)) < 1e-7)

# ---------------------------------------------------------------------------
# 1b. free cubic
# ---------------------------------------------------------------------------

section "1b. free cubic clover operator vs closed form (4^4, anti time)"
block:
  let clat = @[4, 4, 4, 4]
  let clo = clat.newLayout
  var cg = clo.newGauge
  for mu in 0..<4: cg[mu] := 1
  cg.setBC
  var ccw = newCubicCloverWilson(cg, 1.0)
  var cproto = clo.DiracFermion()
  cproto := 0
  type DF = typeof(cproto)
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
  var stats = SiStats()
  var cnt = 0'u64
  let sigma = -0.3
  var op = newShiftInvertOp[DF](
    applyM = proc (rr: var DF; x: DF) = ccw.D(rr, x, 0.3),
    applyMdag = proc (rr: var DF; x: DF) = ccw.Ddag(rr, x, 0.3),
    newVec = proc (): DF =
      result = newOneOf(cproto)
      result := 0,
    startVec = proc (v: var DF) =
      inc cnt
      let salt = sm64(cnt)*0x10000'u64
      for i in clo.sites:
        for sp in 0..3:
          for c in 0..2:
            let k = salt + uint64(i)*24 + uint64(sp)*6 + uint64(c)*2
            setC(v{i}[sp][c], u01(k), u01(k+1)),
    r2req = 1e-18, maxits = 4000, stats = stats)
  let (mus, vecs, _, napply) = arnoldi(op, 12, 48, 1e-9, 60, "LM", 0)
  let modes = measureModes(mus, vecs, sigma,
    proc(r: var DF; x: DF) = ccw.D(r, x, 0.0))
  let s = summarize(modes, 1e-6)
  let wn = worstNearest(modes, exact)
  echo &"  {modes.len} converged, {napply} applies, worst direct resid {s.worstResid:.2e}"
  echo &"  worst |lambda - exact| = {wn:.3e}"
  ok("free cubic: 12 pairs converged", modes.len == 12)
  ok(&"free cubic: direct residuals < 1e-7 (worst {s.worstResid:.1e})", s.worstResid < 1e-7)
  ok(&"free cubic: eigenvalues match closed form to 1e-7 (worst {wn:.1e})", wn < 1e-7)

# ---------------------------------------------------------------------------
# 2 + 3. rough honeycomb config: conjugation symmetry, chirality of complex
#        modes, determinism
# ---------------------------------------------------------------------------

section "2. rough config (warm + 2 stout + setBC): pairing, chirality, resids"
var modesA: seq[SpecMode]
block:
  var r = lo.newRNGField(MRG32k3a, 13579'u64)
  threads:
    hg.warm(0.55, r)
    hg.reunit
  var st = newHcStout(hg, 0.05)
  st.smearN(hg, hg, 2)
  threads:
    hg.setBC
  cw.gaugeRefresh
  let (modes, _) = hcEigs(-0.45, 12, 48, 1e-16)
  modesA = modes
  let s = summarize(modes, 1e-6)
  let cp = worstConjPairing(modes.mapIt(it.lam))
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

section "3. determinism: identical rerun (same seeds, same thread count)"
block:
  let (modes, _) = hcEigs(-0.45, 12, 48, 1e-16)    # hg unchanged: same operator
  var dv = 0.0
  var dc = 0.0
  for i, m in modes:
    dv = max(dv, abs(m.lam - modesA[i].lam))
    dc = max(dc, abs(m.chi - modesA[i].chi))
  echo &"  max |dlambda| = {dv:.3e}, max |dchi| = {dc:.3e}"
  ok("determinism: eigenvalues bit-identical", dv == 0.0)
  ok("determinism: chiralities bit-identical", dc == 0.0)

# ---------------------------------------------------------------------------
# 4. Atiyah-Singer constant-flux configs: Q_Dirac counting and the sign pin
# ---------------------------------------------------------------------------

section "4a. honeycomb flux config (n1 = n2 = 1): Q_Dirac vs hexagon-clover Q"
block:
  setFluxHc(hg, 1, 1)
  let (eF, qF) = hcEQ(wt, hg)
  discard eF
  # Antiperiodic time, as in production.  Without it the color-charge-0
  # component (T = diag(1,-1,0)) is a FREE fermion whose p = 0 modes sit at
  # exactly lambda = 0 with arbitrary chirality mixing (observed); the index
  # zero modes of the charge +-1 components exist for either spin structure.
  threads:
    hg.setBC
  cw.gaugeRefresh
  echo &"  hexagon-clover Q = {qF:.6f} (exact 2 n1 n2 = 2, artifact 1/L^4)"
  let (modes, _) = hcEigs(-0.25, 10, 40, 1e-16)
  let s = summarize(modes, 1e-6)
  var maxChiC = 0.0
  var minAbsChiR = 1e300
  var maxReR = 0.0
  for m in modes:
    if abs(m.lam.im) < 1e-6:
      minAbsChiR = min(minAbsChiR, abs(m.chi))
      maxReR = max(maxReR, abs(m.lam.re))
      echo &"  real mode: lam = {m.lam.re:.6f}  chi = {m.chi:.6f}"
    else:
      maxChiC = max(maxChiC, abs(m.chi))
  echo &"  nreal = {s.nreal}, n+ = {s.nplus}, n- = {s.nminus}, sum chi = {s.sumChiReal:.4f}"
  echo &"  Q_Dirac = qDiracSign*(n+ - n-) = {s.qdirac:.1f}   vs   round(Q_flow-style) = {round(qF):.1f}"
  ok("flux hc: exactly 2 real near-zero modes in the window", s.nreal == 2)
  ok(&"flux hc: real modes near zero (max |Re| {maxReR:.2e})", maxReR < 0.05)
  ok(&"flux hc: real-mode chirality one sign, |chi| > 0.99 (min {minAbsChiR:.4f})",
     minAbsChiR > 0.99 and (s.nplus == 0 or s.nminus == 0))
  ok(&"flux hc: complex modes chi ~ 0 (max {maxChiC:.1e})", maxChiC < 1e-4)
  ok(&"flux hc: Q_Dirac == round(Q) == 2 (PINS qDiracSign = {qDiracSign})",
     abs(s.qdirac - round(qF)) < 0.5)
  ok(&"flux hc: sum of real-mode chiralities integer-ish " &
     &"(|{s.sumChiReal:.3f}| vs 2, within 0.3)", abs(abs(s.sumChiReal) - 2.0) < 0.3)

section "4b. cubic flux config (n1 = n2 = 1): same counting (informational-strict)"
block:
  let clat = @[4, 4, 4, 4]
  let clo = clat.newLayout
  var cg = clo.newGauge
  setFluxCubic(cg, clo, 1, 1)
  let f = cg.fmunu 1
  let qF = f.topoQ
  echo &"  1x1-clover Q = {qF:.6f} (exact 2, artifact O(f^2))"
  cg.setBC                    # antiperiodic time (see 4a)
  var ccw = newCubicCloverWilson(cg, 1.0)
  var cproto = clo.DiracFermion()
  cproto := 0
  type DF = typeof(cproto)
  var stats = SiStats()
  var cnt = 0'u64
  let sigma = -0.25
  var op = newShiftInvertOp[DF](
    applyM = proc (rr: var DF; x: DF) = ccw.D(rr, x, -sigma),
    applyMdag = proc (rr: var DF; x: DF) = ccw.Ddag(rr, x, -sigma),
    newVec = proc (): DF =
      result = newOneOf(cproto)
      result := 0,
    startVec = proc (v: var DF) =
      inc cnt
      let salt = sm64(cnt)*0x30000'u64
      for i in clo.sites:
        for sp in 0..3:
          for c in 0..2:
            let k = salt + uint64(i)*24 + uint64(sp)*6 + uint64(c)*2
            setC(v{i}[sp][c], u01(k), u01(k+1)),
    r2req = 1e-16, maxits = 4000, stats = stats)
  let (mus, vecs, _, _) = arnoldi(op, 10, 40, 1e-9, 60, "LM", 0)
  let modes = measureModes(mus, vecs, sigma,
    proc(r: var DF; x: DF) = ccw.D(r, x, 0.0), 1e-7)
  let s = summarize(modes, 1e-6)
  var minAbsChiR = 1e300
  for m in modes:
    if abs(m.lam.im) < 1e-6:
      minAbsChiR = min(minAbsChiR, abs(m.chi))
      echo &"  real mode: lam = {m.lam.re:.6f}  chi = {m.chi:.6f}"
  echo &"  nreal = {s.nreal}, n+ = {s.nplus}, n- = {s.nminus}, sum chi = {s.sumChiReal:.4f}, Q_Dirac = {s.qdirac:.1f}"
  ok("flux cubic: 2 near-zero real modes found", s.nreal == 2)
  ok(&"flux cubic: one-sign chirality, |chi| > 0.9 (min {minAbsChiR:.4f})",
     minAbsChiR > 0.9 and (s.nplus == 0 or s.nminus == 0))
  ok(&"flux cubic: Q_Dirac == round(Q) == 2 (same sign as honeycomb)",
     abs(s.qdirac - round(qF)) < 0.5)

section "summary"
if nFail > 0:
  echo "FAILED checks: ", nFail
else:
  echo "all checks passed"
qexFinalize()
quit(if nFail > 0: 1 else: 0)
