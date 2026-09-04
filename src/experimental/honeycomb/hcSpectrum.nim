## hcSpectrum -- low-lying Wilson-Dirac spectrum, chirality, Q_Dirac and
## Q_flow on the 16-cell honeycomb (task **D4**; slides 16/18/19 = paper
## Figs. 5-8).
##
## The module has two halves:
##
## 1. Importable machinery, shared with `cubicSpectrum.nim` and
##    `tests/tspectrum.nim`:
##    * the QEX-Field/HcFermion vector-space adapter for `hcarnoldi.arnoldi`
##      (the six mixin ops of tests/tarnoldi.nim, generalised to both vector
##      types; axpy/scale/copy/zero are threaded, dot/norm2 are serial so a
##      run is bit-reproducible at fixed OMP_NUM_THREADS);
##    * `newShiftInvertOp` -- shift-invert Arnoldi operator
##      r = (D - sigma)^{-1} x with real sigma OUTSIDE the spectrum to the
##      LEFT, computed by CG on the normal equations (CGNR):
##      (M^dag M) y = M^dag x with M = D(m = -sigma) (the shift is a mass
##      shift, no extra work).  Eigenvalues come back as mu = 1/(lambda-sigma),
##      i.e. lambda = sigma + 1/mu; "LM" on the inverted operator = closest
##      to sigma = the physical low modes (which sit at Re lambda ~ +0.03..0.1;
##      the brief's picture).  Direct residuals |D v - lambda v|/|v| are
##      recomputed with the EXACT operator afterwards (the inexact inner solve
##      makes arnoldi's own residuals optimistic).
##    * mode bookkeeping: chirality, real-mode classification, Q_Dirac.
##
## 2. `when isMainModule`: the honeycomb driver.  Generates its own heatbath
##    ensemble in place (task M machinery; deterministic given -seed), and per
##    configuration measures
##      a. gradient flow on the PERIODIC, UNsmeared links: Q(t) and t^2E(t) on
##         a fixed grid to -flowtmax; Q_flow = Q(t = -t0use) by linear
##         interpolation (t0use = the matched-point t0/a^2, nominal 1.917);
##         the per-config and ensemble-average t^2E = 0.3 crossings are
##         recorded as a scale check (CAVEAT: at L = 8, L/sqrt(t0) ~ 5.8 is
##         below the >~9 bound of tasks C/W, so the measured t0 is
##         finite-volume biased -- the matching leans on task R's calibration);
##      b. 6 stout steps (rho = 0.05, the paper's parameters; NOTE from task
##         D2 that in the plain MP convention this is a ~sqrt(3)x smaller
##         smearing radius on the honeycomb: t_eff per step = rho/3 vs rho.
##         -rho2/-rho2every measures a subset at a second rho for robustness),
##         antiperiodic time BC, then the massless r=1 cSW=1 clover operator
##         (task D2) and the shift-invert eigensolve;
##      c. per mode: lambda, chirality chi = Re<v|g5|v>/<v|v>, direct residual;
##         real modes (|Im lambda| < -epsreal), Q_Dirac = qDiracSign*(n+ - n-).
##
## Output: parseable lines (also copied to -outfile, flushed per config):
##   EIG  icfg rho k re(lam) im(lam) chi resid
##   FLOWQ icfg t0cfg qflow q(1.5) q(2.4) qraw0
##   CFG  icfg rho lam0re lam0im nconv nreal nplus nminus qdirac qflow
##        napply cgits secs
##   ENST2E t <t2E>   /  ENST0 t0ens   (end of run)
##
##   export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
##   cd build_mac && make src/experimental/honeycomb/hcSpectrum.nim
##   OMP_NUM_THREADS=4 ./bin/hcSpectrum -geom:8,8,8,8 -beta:7.22 -ncfg:100 \
##     -nev:32 -sigma:-0.25 -outfile:hc722.log

import std/[math, complex, os, strformat, strutils, monotimes, times]
import qex except epsilon
import physics/qcdTypes
import hcgeom, hclayout, hcgauge, hcaction, hcflow, hctopo, hcio
import hcwilson, hcstout, hcclover
import hcanalysis
import hcarnoldi

export hcarnoldi, hcwilson, hcclover, hcstout, hcflow, hctopo, hcio

const qDiracSign* = -1.0
  ## Q_Dirac = qDiracSign * (n+ - n-), with n+/n- the number of real
  ## low-lying eigenvalues with positive/negative chirality.  The sign is
  ## PINNED by tests/tspectrum.nim test 4 on the exact constant-flux
  ## Atiyah-Singer background (hexagon-clover Q = +2 n1 n2, task W): the two
  ## near-zero real modes come out with chirality -sign(n1 n2) in our
  ## DeGrand-Rossi/hctopo conventions, so n+ - n- = -Q there.

# ---------------------------------------------------------------------------
# scalar/lane helpers (tests/tarnoldi.nim conventions)
# ---------------------------------------------------------------------------

template toF*(x: untyped): float =
  block:
    var v: float
    v := x
    v

template setC*(dest: untyped; a, b: float) =
  dest.re := a
  dest.im := b

proc sm64*(x: uint64): uint64 =
  ## splitmix64: deterministic start-vector noise
  var z = x + 0x9e3779b97f4a7c15'u64
  z = (z xor (z shr 30)) * 0xbf58476d1ce4e5b9'u64
  z = (z xor (z shr 27)) * 0x94d049bb133111eb'u64
  z xor (z shr 31)

proc u01*(x: uint64): float =   # in [-0.5, 0.5)
  float(sm64(x) shr 11) * (1.0/9007199254740992.0) - 0.5

# ---------------------------------------------------------------------------
# vector-space adapter: QEX Field AND HcFermion
# All driver-level vector ops are SERIAL (single-threaded), like
# tests/tarnoldi.nim: on this box a `threads:` fork-join costs ~0.3 ms under
# load, far more than a serial memory pass over an 8^4 fermion (~0.1 ms), so
# threading them makes the CG inner loop 2-6x SLOWER (measured).  Only the
# operator applications thread internally.  Serial ops also make runs
# bit-reproducible independent of OMP_NUM_THREADS for everything except the
# matvec's own reduction order.
# ---------------------------------------------------------------------------

proc axpyP*[F](y: HcFermion[F]; a: Complex64; x: HcFermion[F]) =
  let z = newComplex(a.re, a.im)
  y.a += z*x.a
  y.b += z*x.b

proc axpyP*(y: auto; a: Complex64; x: auto) =
  let z = newComplex(a.re, a.im)
  y += z*x

proc scaleP*[F](v: HcFermion[F]; s: float) =
  v.a := s*v.a
  v.b := s*v.b

proc scaleP*(v: auto; s: float) =
  v := s*v

proc xpbyP*[F](p: HcFermion[F]; x: HcFermion[F]; b: float) =
  ## p = x + b*p (CG search-direction update, single pass)
  p.a := x.a + b*p.a
  p.b := x.b + b*p.b

proc xpbyP*(p: auto; x: auto; b: float) =
  p := x + b*p

proc copyP*(dst: auto; src: auto) =
  dst := src

proc zeroP*(v: auto) =
  v := 0

template vcopy*(dst, src: untyped) = copyP(dst, src)
template vzero*(v: untyped) = zeroP(v)
template vscale*(v: untyped; s: float) = scaleP(v, s)
template vaxpy*(y: untyped; a: Complex64; x: untyped) = axpyP(y, a, x)
template vdot*(x, y: untyped): Complex64 =
  block:
    let d = dot(x, y)          # conjugates the FIRST argument (QEX convention)
    complex64(toF d.re, toF d.im)
template vnorm2*(x: untyped): float = x.norm2

# ---------------------------------------------------------------------------
# shift-invert operator
# ---------------------------------------------------------------------------

type SiStats* = ref object
  ## accumulated inner-solver cost/quality (reset per configuration)
  nSolve*: int          ## outer applies
  totIts*: int          ## total CGNR iterations
  maxUsedIts*: int      ## largest single-solve iteration count
  nHitMax*: int         ## solves that hit maxits (inner tolerance not met)
  worstDirect*: float   ## worst directly-checked |M y - x|/|x|
  nDirect*: int         ## how many solves were directly checked

proc reset*(s: SiStats) =
  s.nSolve = 0
  s.totIts = 0
  s.maxUsedIts = 0
  s.nHitMax = 0
  s.worstDirect = 0.0
  s.nDirect = 0

proc newShiftInvertOp*[F](
    applyM: proc (r: var F; x: F);      # r = (D - sigma) x  [= D at m = -sigma]
    applyMdag: proc (r: var F; x: F);   # adjoint of applyM
    newVec: proc (): F;                 # fresh zeroed vector
    startVec: proc (v: var F);          # deterministic start vector
    r2req: float;                       # CGNR stop: |r|^2 <= r2req |M^dag x|^2
    maxits: int;
    stats: SiStats;
    nDirectCheck = 3                    # directly verify the first N solves
  ): ArnoldiOp[F] =
  ## ArnoldiOp computing r = (D - sigma)^{-1} x by CG on the normal equations.
  ## All vector work uses the adapter above (deterministic at fixed threads).
  var rhs = newVec()
  var rr = newVec()
  var pp = newVec()
  var mp = newVec()
  var ap = newVec()
  result = ArnoldiOp[F](
    newVec: newVec,
    start: startVec,
    apply: proc (y: var F; x: F) =
      applyMdag(rhs, x)
      zeroP(y)
      copyP(rr, rhs)
      copyP(pp, rr)
      let rhs2 = vnorm2(rhs)
      let stop = r2req*rhs2
      var rho = rhs2
      var its = 0
      while rho > stop and its < maxits:
        applyM(mp, pp)
        applyMdag(ap, mp)
        let pap = vnorm2(mp)             # <p, M^dag M p> = |M p|^2
        if pap <= 0.0: break
        let alpha = rho/pap
        axpyP(y, complex64(alpha, 0.0), pp)
        axpyP(rr, complex64(-alpha, 0.0), ap)
        let rhoNew = vnorm2(rr)
        let beta = rhoNew/rho
        rho = rhoNew
        xpbyP(pp, rr, beta)
        inc its
      stats.totIts += its
      inc stats.nSolve
      if its > stats.maxUsedIts: stats.maxUsedIts = its
      if its >= maxits: inc stats.nHitMax
      if stats.nDirect < nDirectCheck:   # true-residual spot check
        inc stats.nDirect
        applyM(mp, y)
        axpyP(mp, complex64(-1.0, 0.0), x)
        let d = sqrt(vnorm2(mp)/max(vnorm2(x), 1e-300))
        if d > stats.worstDirect: stats.worstDirect = d
  )

# ---------------------------------------------------------------------------
# mode bookkeeping
# ---------------------------------------------------------------------------

type SpecMode* = object
  lam*: Complex64      ## eigenvalue of D (= sigma + 1/mu)
  chi*: float          ## Re <v|gamma5|v> / <v|v>
  resid*: float        ## direct |D v - lam v| / |v| (exact operator)

proc lamFromMu*(mu: Complex64; sigma: float): Complex64 =
  complex64(sigma, 0.0) + complex64(1.0, 0.0)/mu

proc measureModes*[F](mus: openArray[Complex64]; vecs: openArray[F];
                      sigma: float; applyD: proc(r: var F; x: F);
                      residcut = -1.0): seq[SpecMode] =
  ## vecs contains nonzero eigenvectors paired with mus. applyD evaluates D.
  ## A negative residcut retains all candidates; otherwise keep direct residuals <= residcut.
  if mus.len != vecs.len:
    raise newException(ValueError, "eigenvalues and eigenvectors must have equal lengths")
  if mus.len == 0: return
  var tmp = newOneOf(vecs[0])
  for i, mu in mus:
    var m = SpecMode(lam: lamFromMu(mu, sigma))
    let v = vecs[i]
    let n2 = vnorm2(v)
    applyD(tmp, v)
    axpyP(tmp, -m.lam, v)
    m.resid = sqrt(vnorm2(tmp)/n2)
    when F is HcFermion:
      applyGamma5(tmp, v)
    else:
      for e in tmp:
        tmp[e] := gamma5 * v[e]
    m.chi = vdot(v, tmp).re/n2
    if residcut >= 0.0 and m.resid > residcut: continue
    result.add m

type SpecSummary* = object
  lam0*: Complex64     ## converged mode with the smallest Re lambda
  nconv*, nreal*, nplus*, nminus*: int
  qdirac*: float       ## qDiracSign*(nplus - nminus)
  sumChiReal*: float   ## sum of chiralities over the real modes
  minAbsImComplex*: float  ## smallest |Im lam| among complex-classified modes
  maxAbsImReal*: float     ## largest  |Im lam| among real-classified modes
  worstResid*: float

proc summarize*(modes: openArray[SpecMode]; epsReal: float): SpecSummary =
  result.lam0 = complex64(1e300, 0.0)
  result.minAbsImComplex = 1e300
  for m in modes:
    inc result.nconv
    if m.resid > result.worstResid: result.worstResid = m.resid
    if m.lam.re < result.lam0.re: result.lam0 = m.lam
    if abs(m.lam.im) < epsReal:
      inc result.nreal
      result.sumChiReal += m.chi
      if m.chi > 0.0: inc result.nplus else: inc result.nminus
      if abs(m.lam.im) > result.maxAbsImReal:
        result.maxAbsImReal = abs(m.lam.im)
    else:
      if abs(m.lam.im) < result.minAbsImComplex:
        result.minAbsImComplex = abs(m.lam.im)
  result.qdirac = qDiracSign*float(result.nplus - result.nminus)

proc worstConjPairing*(vals: openArray[Complex64]): float =
  ## worst distance from conj(v) to the nearest member of vals
  ## (gamma5-hermiticity check: the converged set must be conj symmetric)
  for v in vals:
    var d = 1e300
    for u in vals:
      d = min(d, abs(conjugate(v) - u))
    result = max(result, d)

proc specLines*(modes: openArray[SpecMode]; epsreal: float; icfg: int;
                rho, qflow: float; nbad, napply: int; stats: SiStats;
                secs: float): seq[string] =
  ## secs measures computation through measureModes, before formatting or writing.
  let s = summarize(modes, epsreal)
  var vals: seq[Complex64]
  for k, m in modes:
    vals.add m.lam
    result.add &"EIG {icfg} {rho:.3f} {k} {m.lam.re:.10g} {m.lam.im:.10g} {m.chi:.8f} {m.resid:.3e}"
  let cp = worstConjPairing(vals)
  result.add &"CFG {icfg} {rho:.3f} {s.lam0.re:.10g} {s.lam0.im:.10g} " &
             &"{s.nconv} {s.nreal} {s.nplus} {s.nminus} {s.qdirac:.1f} {qflow:.6f} " &
             &"{napply} {stats.totIts} {secs:.2f}"
  result.add &"CFGX {icfg} {rho:.3f} conj {cp:.3e} worstresid {s.worstResid:.3e} " &
             &"imgapC {s.minAbsImComplex:.3e} imgapR {s.maxAbsImReal:.3e} " &
             &"sumchi {s.sumChiReal:.4f} nbad {nbad} cgmax {stats.maxUsedIts} " &
             &"cghitmax {stats.nHitMax} sidirect {stats.worstDirect:.3e}"

# ===========================================================================
# the honeycomb driver
# ===========================================================================

when isMainModule:
  import hcheatbath

  qexInit()
  tic()

  letParam:
    geom = @[8, 8, 8, 8]     ## cell geometry (VLEN-compatible: 4^4, 8^4, 12^4)
    beta = 7.22              ## triangle-action beta (task R: t0/a^2 ~ 1.9)
    seed: uint64 = 20260821'u64
    nwarm = 300              ## heatbath thermalisation updates
    stride = 10              ## heatbath updates (1 HB + norsweeps OR) between configs
    norsweeps = 3
    ncfg = 100               ## measured configurations
    firstcfg = 1             ## skip measurements for icfg < firstcfg (restart aid)
    nstout = 6               ## stout steps
    rho = 0.05               ## stout rho (paper value)
    rho2 = 0.15              ## second rho for the robustness subset (0 = off)
    rho2every = 0            ## measure rho2 on every Nth config (0 = off)
    mass = 0.0               ## valence quark mass
    rw = 1.0                 ## Wilson r
    cSW = 1.0                ## clover coefficient
    sigma = -0.25            ## shift-invert point (real, left of the spectrum)
    nev = 32                 ## converged pairs requested
    ncv = 0                  ## Krylov size (0 = 3*nev)
    tol = 1e-7               ## outer (Arnoldi) tolerance
    maxRestarts = 60
    residcut = 2e-5          ## exclude modes with direct residual above this
    innerR2 = 1e-12          ## CGNR |r|^2/|rhs|^2 stop
    innerMaxIts = 4000
    epsreal = 1e-6           ## |Im lambda| below this = real mode
    doflow: bool = 1
    doeigs: bool = 1
    floweps = 0.05           ## RK3 step (continuum-normalised flow time)
    flowmeas = 2             ## measure (E,Q) every this many RK steps
    flowtmax = 2.6           ## flow to this t (fixed grid, no early stop)
    t0use = 1.917            ## Q_flow = Q(t0use); nominal matched t0/a^2
    verb = 0                 ## arnoldi verbosity
    outfile = ""             ## copy of all output lines
    bench = 0                ## time D_c/axpy/dot at this geometry, then exit
    showTimers: bool = 0

  installHelpParam()
  echoParams()
  echo "rank ", myRank, "/", nRanks
  threads: echo "thread ", threadNum, "/", numThreads

  let ncvUse = if ncv > 0: ncv else: 3*nev
  var fh: File
  let haveOut = outfile.len > 0 and myRank == 0
  if haveOut: fh = open(outfile, fmWrite)
  proc emit(line: string) =
    echo line
    if haveOut:
      fh.writeLine line
      fh.flushFile

  emit &"# hcSpectrum geom {geom} beta {beta} seed {seed} nwarm {nwarm} stride {stride}"
  emit &"# nstout {nstout} rho {rho} rho2 {rho2} rho2every {rho2every} cSW {cSW} mass {mass} rw {rw}"
  emit &"# sigma {sigma} nev {nev} ncv {ncvUse} tol {tol} innerR2 {innerR2} epsreal {epsreal}"
  emit &"# floweps {floweps} flowmeas {flowmeas} flowtmax {flowtmax} t0use {t0use}"
  emit &"# qDiracSign {qDiracSign}"

  let hl = newHcLayout(geom)
  var r = hl.lo.newRNGField(MRG32k3a, seed)
  var g = newHcGauge(hl)       # the Markov-chain configuration (periodic)
  var gf = newOneOf(g)         # flow copy
  var gs = newOneOf(g)         # smeared + BC copy (the operator lives on this)
  threads:
    gf := g
    gs := g
  var hb = newHcHeatbath(g, beta)
  var st1 = newHcStout(gs, rho)
  var st2 = newHcStout(gs, if rho2 > 0: rho2 else: rho)
  var wt = newHcTopoWork(gf)
  var cw = newHcCloverWilson(gs, cSW)
  let lo = hl.lo

  # fermion workspace
  var proto = newHcFermion(hl)
  type HF = typeof(proto)

  # deterministic start vectors
  var startCount = 0'u64
  proc startVec(v: var HF) =
    inc startCount
    let salt = sm64(startCount)*0x10000'u64
    for i in lo.sites:
      for sp in 0..3:
        for c in 0..2:
          let k = salt + uint64(i)*48 + uint64(sp)*12 + uint64(c)*4
          setC(v.a{i}[sp][c], u01(k), u01(k+1))
          setC(v.b{i}[sp][c], u01(k+2), u01(k+3))

  var stats = SiStats()
  let msolve = mass - sigma    # M = D - sigma = D at mass m - sigma
  var op = newShiftInvertOp[HF](
    applyM = proc (rr: var HF; x: HF) = cw.D(rr, x, msolve, rw),
    applyMdag = proc (rr: var HF; x: HF) = cw.Ddag(rr, x, msolve, rw),
    newVec = proc (): HF = newOneOf(proto),
    startVec = startVec,
    r2req = innerR2, maxits = innerMaxIts, stats = stats)

  # flow bookkeeping (fixed grid so the ensemble average is trivial)
  let nFlowSteps = int(round(flowtmax/floweps))
  var flowT: seq[float]
  var sumT2E: seq[float]
  var nFlowCfg = 0
  var lastQflow = 0.0          # set by measureFlow, consumed by the CFG line

  proc measureFlow(icfg: int) =
    threads:
      gf := g
    var ts = @[0.0]
    var t2Es = @[0.0]
    var qs: seq[float]
    block:
      let (e0, q0) = hcEQ(wt, gf)
      discard e0
      qs.add q0
    var nstepF = 0
    gf.hcGaugeFlow(nFlowSteps, floweps, hcFlowCflow):
      inc nstepF
      if nstepF mod flowmeas == 0:
        let (e, q) = hcEQ(wt, gf)
        ts.add wflowT
        t2Es.add wflowT*wflowT*e
        qs.add q
    # ensemble accumulation on the fixed grid
    if flowT.len == 0:
      flowT = ts
      sumT2E = newSeq[float](ts.len)
    inc nFlowCfg
    for j in 0..<min(t2Es.len, sumT2E.len):
      sumT2E[j] += t2Es[j]
    let t0cfg = findT0(ts, t2Es, 0.3, 1)
    let qflow = interpAt(ts, qs, t0use)
    let q15 = interpAt(ts, qs, 1.5)
    let q24 = interpAt(ts, qs, min(2.4, ts[^1]))
    lastQflow = qflow            # consumed by the CFG summary line
    emit &"FLOWQ {icfg} {t0cfg:.6f} {qflow:.6f} {q15:.6f} {q24:.6f} {qs[0]:.6f}"

  proc measureSpec(icfg: int; rhoTag: float; st: typeof(st1)) =
    let tw0 = getMonoTime()
    threads:
      gs := g
    st.smearN(gs, gs, nstout)
    threads:
      gs.setBC
    cw.gaugeRefresh
    stats.reset
    let (mus, vecs, _, napply) =
      arnoldi(op, nev, ncvUse, tol, maxRestarts, "LM", verb)
    let modes = measureModes(mus, vecs, sigma,
      proc(r: var HF; x: HF) = cw.D(r, x, mass, rw), residcut)
    let secs = (getMonoTime() - tw0).inMicroseconds.float*1e-6
    for line in specLines(modes, epsreal, icfg, rhoTag, lastQflow,
                          mus.len - modes.len, napply, stats, secs):
      emit line

  if bench > 0:
    var x1 = newOneOf(proto)
    var y1 = newOneOf(proto)
    startVec(x1)
    template timeIt(nm: string; n: int; body: untyped) =
      block:
        body                       # warm up
        let t0b = getMonoTime()
        for it in 1..n: body
        let dt = (getMonoTime() - t0b).inMicroseconds.float*1e-3/float(n)
        echo "BENCH ", nm, ": ", dt.formatFloat(ffDecimal, 3), " ms"
    timeIt("D_c apply", bench):
      cw.D(y1, x1, msolve, rw)
    timeIt("Ddag_c apply", bench):
      cw.Ddag(y1, x1, msolve, rw)
    timeIt("axpy", bench):
      axpyP(y1, complex64(0.5, 0.0), x1)
    timeIt("xpby", bench):
      xpbyP(y1, x1, 0.5)
    timeIt("norm2", bench):
      discard vnorm2(y1)
    timeIt("dot", bench):
      discard vdot(x1, y1)
    timeIt("CG iteration (2 matvec + 3 vec + 2 red)", bench):
      cw.D(y1, x1, msolve, rw)
      cw.Ddag(y1, x1, msolve, rw)
      axpyP(y1, complex64(1e-8, 0.0), x1)
      axpyP(y1, complex64(-1e-8, 0.0), x1)
      xpbyP(y1, x1, 1e-8)
      discard vnorm2(y1)
      discard vnorm2(x1)
    qexFinalize()
    quit 0

  # ------------------------------------------------------------------
  # the run
  # ------------------------------------------------------------------
  toc("setup")
  echo "# thermalising: ", nwarm, " heatbath updates"
  for u in 1..nwarm:
    hb.update(g, r, norsweeps)
  threads:
    g.reunit
  toc("warmup")
  emit &"# thermalised: triangleSum = {g.triangleSum:.8f}"

  let trun0 = getMonoTime()
  for icfg in 1..ncfg:
    for u in 1..stride:
      hb.update(g, r, norsweeps)
    threads:
      g.reunit
    if icfg < firstcfg: continue
    emit &"TS {icfg} {g.triangleSum:.8f}"
    lastQflow = 0.0
    if doflow: measureFlow(icfg)
    if doeigs:
      measureSpec(icfg, rho, st1)
      if rho2every > 0 and rho2 > 0 and icfg mod rho2every == 0:
        measureSpec(icfg, rho2, st2)
    GC_fullCollect()
    let el = (getMonoTime() - trun0).inMicroseconds.float*1e-6
    emit &"# done cfg {icfg}/{ncfg}, elapsed {el:.1f} s ({el/float(icfg-firstcfg+1):.1f} s/cfg)"

  if nFlowCfg > 0:
    for j in 0..<flowT.len:
      emit &"ENST2E {flowT[j]:.4f} {sumT2E[j]/float(nFlowCfg):.8f}"
    var meanT2E = newSeq[float](flowT.len)
    for j in 0..<flowT.len: meanT2E[j] = sumT2E[j]/float(nFlowCfg)
    let t0ens = findT0(flowT, meanT2E, 0.3, 1)
    emit &"ENST0 {t0ens:.6f}   # ensemble <t^2E> = 0.3 crossing; " &
         "finite-volume biased at L/sqrt(t0) < 9"
  toc("run")
  if haveOut: fh.close
  if showTimers: echoTimers()
  qexFinalize()
