## Low-mode machinery shared by spectrum.nim and its tests: the vector-space
## adapter for hcarnoldi (QEX Field or HcFermion), the shift-invert operator,
## and the per-mode bookkeeping (chirality, real modes, Q_Dirac).
##
## Shift-invert: r = (D - sigma)^{-1} x with real sigma left of the spectrum,
## by CG on the normal equations (M^dag M) y = M^dag x, M = D(m - sigma).
## Eigenvalues return as mu = 1/(lambda - sigma); "LM" on the inverted
## operator = nearest to sigma.  Direct residuals |D v - lambda v|/|v| are
## recomputed with the exact operator.  Vector ops are serial: a threads
## fork/join costs more than a pass over an 8^4 fermion, and serial reductions
## keep runs bit-reproducible at any thread count.
##
##   chi = Re <v|gamma5|v> / <v|v>,   real mode: |Im lambda| < epsReal,
##   Q_Dirac = qDiracSign (n+ - n-) over real modes with Re lambda < recut,
##   reach = sigma + max_i |lambda_i - sigma|: the real modes below recut are
##   complete iff recut <= reach.
## qDiracSign = -1 is pinned by tests/tspectrum.nim on the constant-flux
## background (Q = +2 n1 n2, zero modes have chirality -sign(n1 n2)).

import std/[math, complex, strformat]
import qex except epsilon
import physics/qcdTypes
import hcwilson
import hcarnoldi

export hcarnoldi, hcwilson

const qDiracSign* = -1.0

# ---------------------------------------------------------------------------
# scalar/lane helpers
# ---------------------------------------------------------------------------

proc sm64*(x: uint64): uint64 =
  ## splitmix64 for deterministic start vectors
  var z = x + 0x9e3779b97f4a7c15'u64
  z = (z xor (z shr 30)) * 0xbf58476d1ce4e5b9'u64
  z = (z xor (z shr 27)) * 0x94d049bb133111eb'u64
  z xor (z shr 31)

proc u01*(x: uint64): float =
  ## in [-0.5, 0.5)
  float(sm64(x) shr 11) * (1.0/9007199254740992.0) - 0.5

# ---------------------------------------------------------------------------
# vector-space adapter: QEX Field and HcFermion
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
  ## p = x + b p
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
    let d = dot(x, y)          # conjugates the first argument
    complex64(toF d.re, toF d.im)
template vnorm2*(x: untyped): float = x.norm2

# ---------------------------------------------------------------------------
# shift-invert operator
# ---------------------------------------------------------------------------

type SiStats* = ref object
  ## inner-solver cost and quality, reset per configuration
  nSolve*: int          ## outer applies
  totIts*: int          ## total CGNR iterations
  maxUsedIts*: int
  nHitMax*: int         ## solves that hit maxits
  worstDirect*: float   ## worst directly checked |M y - x|/|x|
  nDirect*: int

proc reset*(s: SiStats) =
  s.nSolve = 0
  s.totIts = 0
  s.maxUsedIts = 0
  s.nHitMax = 0
  s.worstDirect = 0.0
  s.nDirect = 0

proc newShiftInvertOp*[F](
    applyM: proc (r: var F; x: F);      # r = (D - sigma) x
    applyMdag: proc (r: var F; x: F);
    newVec: proc (): F;                 # fresh zeroed vector
    startVec: proc (v: var F);          # deterministic start vector
    r2req: float;                       # stop: |r|^2 <= r2req |M^dag x|^2
    maxits: int;
    stats: SiStats;
    nDirectCheck = 3                    # directly verify the first N solves
  ): ArnoldiOp[F] =
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
        let pap = vnorm2(mp)
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
      if stats.nDirect < nDirectCheck:
        inc stats.nDirect
        applyM(mp, y)
        axpyP(mp, complex64(-1.0, 0.0), x)
        let d = sqrt(vnorm2(mp)/vnorm2(x))
        if d > stats.worstDirect: stats.worstDirect = d
  )

# ---------------------------------------------------------------------------
# mode bookkeeping
# ---------------------------------------------------------------------------

type SpecMode* = object
  lam*: Complex64      ## eigenvalue of D (= sigma + 1/mu)
  chi*: float          ## Re <v|gamma5|v> / <v|v>
  resid*: float        ## direct |D v - lam v| / |v|

proc lamFromMu*(mu: Complex64; sigma: float): Complex64 =
  complex64(sigma, 0.0) + complex64(1.0, 0.0)/mu

proc measureModes*[F](mus: openArray[Complex64]; vecs: openArray[F];
                      sigma: float; applyD: proc(r: var F; x: F);
                      residcut = -1.0): seq[SpecMode] =
  ## eigenvalue, chirality and direct residual of each (mu, v) pair; a
  ## negative residcut keeps all, otherwise modes with resid <= residcut
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
  lam0*: Complex64     ## mode with the smallest Re lambda
  reach*: float        ## sigma + max |lambda - sigma|
  nconv*, nreal*, nplus*, nminus*: int
  nrealAbove*: int     ## real modes with Re lambda >= recut (not counted)
  qdirac*: float       ## qDiracSign (nplus - nminus)
  sumChiReal*: float
  minAbsImComplex*: float  ## smallest |Im lam| among complex modes
  maxAbsImReal*: float     ## largest |Im lam| among real modes
  worstResid*: float

proc summarize*(modes: openArray[SpecMode]; epsReal, sigma, recut: float):
    SpecSummary =
  ## modes must be non-empty; recut is the real-mode window (see module docs)
  result.lam0 = modes[0].lam
  result.reach = sigma
  var haveComplex = false
  for m in modes:
    inc result.nconv
    if m.resid > result.worstResid: result.worstResid = m.resid
    if m.lam.re < result.lam0.re: result.lam0 = m.lam
    result.reach = max(result.reach, sigma + abs(m.lam - complex64(sigma, 0.0)))
    if abs(m.lam.im) < epsReal:
      if m.lam.re >= recut:
        inc result.nrealAbove
        continue
      inc result.nreal
      result.sumChiReal += m.chi
      if m.chi > 0.0: inc result.nplus else: inc result.nminus
      if abs(m.lam.im) > result.maxAbsImReal:
        result.maxAbsImReal = abs(m.lam.im)
    else:
      if not haveComplex or abs(m.lam.im) < result.minAbsImComplex:
        result.minAbsImComplex = abs(m.lam.im)
        haveComplex = true
  result.qdirac = qDiracSign*float(result.nplus - result.nminus)

proc conjGap*(vals: openArray[Complex64]): float =
  ## worst distance from conj(v) to the nearest member of vals
  ## (gamma5-hermiticity: the converged set must be conjugation symmetric)
  for v in vals:
    var d = abs(conjugate(v) - vals[0])
    for u in vals:
      d = min(d, abs(conjugate(v) - u))
    result = max(result, d)

proc specLines*(modes: openArray[SpecMode]; epsreal, sigma, recut: float;
                icfg: int; rho, qflow, qflowT0: float; nbad, napply: int;
                stats: SiStats; secs: float): seq[string] =
  ## EIG lines per mode, then CFG (summary) and CFGX (diagnostics)
  let s = summarize(modes, epsreal, sigma, recut)
  var vals: seq[Complex64]
  for k, m in modes:
    vals.add m.lam
    result.add &"EIG {icfg} {rho:.3f} {k} {m.lam.re:.10g} {m.lam.im:.10g} {m.chi:.8f} {m.resid:.3e}"
  let cp = conjGap(vals)
  result.add &"CFG {icfg} {rho:.3f} {s.lam0.re:.10g} {s.lam0.im:.10g} " &
             &"{s.nconv} {s.nreal} {s.nplus} {s.nminus} {s.qdirac:.1f} {qflow:.6f} " &
             &"{qflowT0:.6f} {napply} {stats.totIts} {secs:.2f}"
  result.add &"CFGX {icfg} {rho:.3f} reach {s.reach:.6f} recut {recut:.6f} " &
             &"realabove {s.nrealAbove} conj {cp:.3e} worstresid {s.worstResid:.3e} " &
             &"imgapC {s.minAbsImComplex:.3e} imgapR {s.maxAbsImReal:.3e} " &
             &"sumchi {s.sumChiReal:.4f} nbad {nbad} cgmax {stats.maxUsedIts} " &
             &"cghitmax {stats.nHitMax} sidirect {stats.worstDirect:.3e}"
