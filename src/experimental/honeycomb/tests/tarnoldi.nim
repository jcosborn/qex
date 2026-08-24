## Task D3 test: the Krylov-Schur Arnoldi eigensolver `hcarnoldi.nim` driven
## through QEX fields, verified on operators whose spectrum is known exactly:
##
##   0. QEX field-op conventions (dot conjugation, norm2, complex axpy) are
##      *verified*, not assumed, against a manual site loop;
##   1. a diagonal operator on a lattice Complex field (1024 sites, known
##      values): "LM"/"SM", values to 1e-10, direct residuals < 1e-8;
##   2. the free cubic Wilson-Dirac operator (unit gauge, antiperiodic time,
##      m = 0.3) on 4^4 and 6^4.  The exact spectrum is extracted FROM THE
##      OPERATOR ITSELF by plane-wave projection (no normalization assumed;
##      a completeness check certifies the projection is exact), then "LM"
##      and "SM" Arnoldi runs are matched against it, plus a gamma5
##      (conjugation-symmetry) check;
##   3. shift-invert: Arnoldi "LM" on D^{-1} (CG on the normal equations,
##      y = (D^dag D)^{-1} D^dag x), 1/lambda vs the smallest exact
##      eigenvalues, with full apply/iteration cost reporting.
##
## Build/run (LAPACK comes from `-framework Accelerate` by default on macOS,
## so NO extra flags are needed; elsewhere add e.g. `:hcLapackLib=-llapack`):
##
##   export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
##   cd build_mac
##   make src/experimental/honeycomb/tests/tarnoldi.nim
##   OMP_NUM_THREADS=4 ./bin/tarnoldi
##
## All field operations in this test (and in the hcarnoldi driver) run
## serially — no `threads:` blocks — which is correct QEX usage (ops fall
## back to a single thread) and deterministic.  The lattices here are tiny;
## a production driver would thread the operator closure.

import math, strformat
import std/complex
import qex except epsilon
import physics/wilsonD
import solvers/cg
import ../hcarnoldi

qexInit()
doAssert nRanks == 1, "tarnoldi is a single-rank test"

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

# ---------------------------------------------------------------------------
# deterministic splitmix64 helpers
# ---------------------------------------------------------------------------

proc sm64(x: uint64): uint64 =
  var z = x + 0x9e3779b97f4a7c15'u64
  z = (z xor (z shr 30)) * 0xbf58476d1ce4e5b9'u64
  z = (z xor (z shr 27)) * 0x94d049bb133111eb'u64
  z xor (z shr 31)

proc u01(x: uint64): float =   # in [-0.5, 0.5)
  float(sm64(x) shr 11) * (1.0/9007199254740992.0) - 0.5

# ---------------------------------------------------------------------------
# QEX single-site element access helpers (lane proxies are not raw floats)
# ---------------------------------------------------------------------------

template toF(x: untyped): float =
  block:
    var v: float
    v := x
    v

template setC(dest: untyped; a, b: float) =
  dest.re := a
  dest.im := b

template getC(src: untyped): Complex64 =
  block:
    var a, b: float
    a := src.re
    b := src.im
    complex64(a, b)

# ---------------------------------------------------------------------------
# the mixin vector-space adapter that hcarnoldi.arnoldi needs for QEX fields
# (any QEX Field type works; reuse this block in hcSpectrum later)
# ---------------------------------------------------------------------------

template vcopy(dst, src: untyped) = dst := src
template vzero(v: untyped) = v := 0
template vscale(v: untyped; s: float) = v := s*v
template vaxpy(y: untyped; a: Complex64; x: untyped) =
  y += newComplex(a.re, a.im)*x
template vdot(x, y: untyped): Complex64 =
  block:
    let d = dot(x, y)          # verified below: = sum conj(x)*y
    complex64(toF d.re, toF d.im)
template vnorm2(x: untyped): float = x.norm2

# ---------------------------------------------------------------------------
# generic spectrum bookkeeping
# ---------------------------------------------------------------------------

proc worstNearest(vals, exact: openArray[Complex64]): float =
  ## worst over vals of the distance to the nearest exact eigenvalue
  for v in vals:
    var d = 1e300
    for e in exact: d = min(d, abs(v - e))
    result = max(result, d)

proc worstConjPair(vals: openArray[Complex64]): float =
  ## worst over vals of the distance from conj(v) to the nearest val
  for v in vals:
    var d = 1e300
    for u in vals: d = min(d, abs(conjugate(v) - u))
    result = max(result, d)

proc shellCoverage(exact, conv: openArray[Complex64]; which: string;
                   nev: int): float =
  ## distinct values (cluster tolerance 1e-6) among the first nev exact
  ## eigenvalues in `which` order; worst distance from any of them to the
  ## nearest converged value.  This is the completeness check.
  let ord = eigOrder(@exact, which)
  var reps: seq[Complex64]
  for i in 0..<min(nev, exact.len):
    let v = exact[ord[i]]
    var found = false
    for r in reps:
      if abs(v - r) < 1e-6:
        found = true
        break
    if not found: reps.add v
  for r in reps:
    var d = 1e300
    for c in conv: d = min(d, abs(c - r))
    result = max(result, d)

proc convOnly(vals: seq[Complex64]; resids: seq[float];
              tol: float): seq[Complex64] =
  for i in 0..<vals.len:
    if resids[i] <= tol: result.add vals[i]

# ---------------------------------------------------------------------------
# 0 + 1: conventions and the diagonal operator on a lattice Complex field
# ---------------------------------------------------------------------------

proc diagTests =
  section "test 0: QEX field-op conventions (verified, not assumed)"
  var lo = newLayout(@[4, 4, 8, 8])
  let nsite = lo.nSites
  type CF = type(lo.Complex())
  var x = lo.Complex()
  var y = lo.Complex()
  var xs = newSeq[Complex64](nsite)
  var ys = newSeq[Complex64](nsite)
  for i in lo.sites:
    xs[i] = complex64(u01(uint64(i)*4+0), u01(uint64(i)*4+1))
    ys[i] = complex64(u01(uint64(i)*4+2), u01(uint64(i)*4+3))
    setC(x{i}, xs[i].re, xs[i].im)
    setC(y{i}, ys[i].re, ys[i].im)
  # dot convention: conjugate-linear in the FIRST argument
  var dman = complex64(0.0, 0.0)
  var n2man = 0.0
  for i in 0..<nsite:
    dman += conjugate(xs[i])*ys[i]
    n2man += abs2(xs[i])
  let dq = vdot(x, y)
  let n2q = vnorm2(x)
  echo &"  dot(x,y)   QEX ({dq.re:.15g},{dq.im:.15g})  manual ({dman.re:.15g},{dman.im:.15g})"
  echo &"  norm2(x)   QEX {n2q:.15g}  manual {n2man:.15g}"
  ok(&"dot(x,y) = sum conj(x) y  (diff {abs(dq-dman):.2e})", abs(dq - dman) < 1e-12*n2man)
  ok(&"norm2 = sum |x|^2  (diff {abs(n2q-n2man):.2e})", abs(n2q - n2man) < 1e-12*n2man)
  # vaxpy with a complex coefficient
  let a = complex64(0.3, -0.7)
  vaxpy(y, a, x)
  var worst = 0.0
  for i in lo.sites:
    let want = ys[i] + a*xs[i]
    worst = max(worst, abs(getC(y{i}) - want))
  ok(&"vaxpy(y, a, x) matches manual complex axpy (worst {worst:.2e})", worst < 1e-14)
  # vscale
  vscale(y, 0.25)
  var worst2 = 0.0
  for i in lo.sites:
    let want = 0.25*(ys[i] + a*xs[i])
    worst2 = max(worst2, abs(getC(y{i}) - want))
  ok(&"vscale(y, s) matches (worst {worst2:.2e})", worst2 < 1e-14)

  section "test 1: diagonal operator, 1024 sites, LM and SM"
  # Well-separated values: 8 small (0.02..0.16), bulk in [1,5], 8 large
  # (10,12,..,24); angles spread by the golden ratio over the sector
  # (-1,1) rad.  The sector matters for "SM": restarted Arnoldi (like
  # ARPACK) can only reach smallest-magnitude eigenvalues if the origin is
  # not ENCLOSED by the spectrum — with full-circle angles the SM run
  # stagnates by design, not by bug (interior targets need shift-invert,
  # cf. test 3).  A Wilson-Dirac spectrum at m > 0 lives in Re > 0, which
  # this sector mimics.
  proc cvalue(k, n: int): Complex64 =
    let ang = 2.0*((0.6180339887498949*float(k+1)) mod 1.0) - 1.0
    var r: float
    if k < 8: r = 0.02*float(k+1)
    elif k >= n-8: r = 10.0 + 2.0*float(k-(n-8))
    else: r = 1.0 + 4.0*float(k-8)/float(n-17)
    complex64(r*cos(ang), r*sin(ang))
  var cf = lo.Complex()
  var exact = newSeq[Complex64](nsite)
  for i in lo.sites:
    exact[i] = cvalue(i, nsite)
    setC(cf{i}, exact[i].re, exact[i].im)
  var startCount = 0'u64
  var op = ArnoldiOp[CF](
    apply: proc(r: var CF; z: CF) =
      r := cf*z,
    newVec: proc(): CF =
      result = lo.Complex()
      result := 0,
    start: proc(v: var CF) =
      inc startCount
      let salt = startCount*0x100000'u64
      for i in lo.sites:
        setC(v{i}, u01(salt + uint64(i)*2), u01(salt + uint64(i)*2 + 1)))

  for which in ["LM", "SM"]:
    let (vals, vecs, resids, nap) = arnoldi(op, 8, 24, 1e-12, 300, which, verb=1)
    let dord = eigOrder(exact, which)
    let vord = eigOrder(vals, which)
    var worstV = 0.0
    var worstR = 0.0
    for i in 0..<8:
      worstV = max(worstV, abs(vals[vord[i]] - exact[dord[i]]))
      worstR = max(worstR, resids[i])
    echo &"  {which}: worst |ritz - exact| = {worstV:.3e}, worst direct resid = {worstR:.3e}, {nap} applies"
    ok(&"diag {which}: 8 values to 1e-10 (worst {worstV:.1e})", worstV < 1e-10)
    ok(&"diag {which}: direct residuals < 1e-8 (worst {worstR:.1e})", worstR < 1e-8)
    doAssert vecs.len == 8

# ---------------------------------------------------------------------------
# 2: free cubic Wilson-Dirac
# ---------------------------------------------------------------------------

proc makeWilson(lat: seq[int]; V: static int): auto =
  # V = SIMD inner length: 6^4 is not layoutable with the default VLEN = 4
  # (odd outer dimensions), so the 6^4 test uses V = 1.  The parameterless
  # constructors (lo.DiracFermion(), newWilson(g)) are VLEN-only, so build
  # the prototype fermion field for general V explicitly (the element type
  # must be Simd-wrapped even for V = 1, or the shift machinery breaks).
  var lo = newLayout(lat, V)
  var g: seq[type(lo.ColorMatrix(3))]
  for mu in 0..<4:
    g.add lo.ColorMatrix(3)
    g[mu] := 1
  g.setBC                     # antiperiodic in time (direction 3)
  type C = type(lo.newDComplexV)
  type DF = Spin[VectorArray[4, Color[VectorArray[3, C]]]]
  var proto = lo.newField(DF)
  proto := 0
  var s = newWilson(@g, proto)
  (lo: lo, g: g, s: s, proto: proto)

proc wilsonExact(w: auto; mass: float): seq[Complex64] =
  ## Extract the exact free spectrum from the operator itself: for every
  ## lattice momentum (half-integer in time, matching setBC) build the four
  ## spin plane waves psi_s(x) = e_s (x) e_{color 0} e^{i p.x}, apply D once,
  ## and project D(p)_{s's} = <psi_{s'}|D psi_s>/V.  A completeness check
  ## certifies that D psi_s lies in the plane-wave block exactly, i.e. that
  ## the momentum set (and hence the BC handling) is right.  Eigenvalues of
  ## the 4x4 blocks (via LAPACK) form the exact spectrum (color adds a
  ## trivial multiplicity 3).
  let lo = w.lo
  let s = w.s
  let geom = lo.physGeom
  let vol = float(lo.nSites)
  type F = type(w.proto)
  var phi: array[4, F]
  var dphi: array[4, F]
  for sp in 0..3:
    phi[sp] = newOneOf(w.proto)
    dphi[sp] = newOneOf(w.proto)
  var worstComp = 0.0          # completeness defect
  var worstNorm = 0.0          # plane-wave normalization defect
  var worstForm = 0.0          # informational: deviation from the closed form
  var kco = newSeq[int](4)
  var nmom = 0
  while true:
    var p: array[4, float]
    for mu in 0..<4:
      let off = if mu == 3: 0.5 else: 0.0   # antiperiodic time
      p[mu] = 2.0*PI*(float(kco[mu]) + off)/float(geom[mu])
    for sp in 0..3:
      phi[sp] := 0
      for i in lo.sites:
        var th = 0.0
        for mu in 0..<4: th += p[mu]*float(lo.coords[mu][i])
        setC(phi[sp]{i}[sp][0], cos(th), sin(th))
      s.D(dphi[sp], phi[sp], mass)
    worstNorm = max(worstNorm, abs(vnorm2(phi[0]) - vol)/vol)
    var m4 = newZMat(4, 4)
    for sp in 0..3:
      for sq in 0..3:
        let z = vdot(phi[sq], dphi[sp])
        m4[sq, sp] = complex64(z.re/vol, z.im/vol)
    for sp in 0..3:            # completeness: |D psi|^2 == V sum_q |m4[q,sp]|^2
      var pr = 0.0
      for sq in 0..3: pr += abs2(m4[sq, sp])
      let nn = vnorm2(dphi[sp])
      worstComp = max(worstComp, abs(nn - vol*pr)/max(nn, 1e-30))
    let ev = zeig(m4).w
    for e in ev: result.add e
    # informational cross-check against the standard Wilson closed form
    var mm = mass
    var k2 = 0.0
    for mu in 0..<4:
      mm += 1.0 - cos(p[mu])
      k2 += sin(p[mu])*sin(p[mu])
    let cplus = complex64(mm, sqrt(k2))
    for e in ev:
      worstForm = max(worstForm, min(abs(e - cplus), abs(e - conjugate(cplus))))
    inc nmom
    # next momentum
    var mu = 0
    while mu < 4:
      inc kco[mu]
      if kco[mu] < geom[mu]: break
      kco[mu] = 0
      inc mu
    if mu == 4: break
  echo &"  {nmom} momenta, {result.len} exact eigenvalues (x3 color)"
  echo &"  plane-wave norm defect {worstNorm:.2e}, projection completeness defect {worstComp:.2e}"
  echo &"  closed-form m+sum(1-cos p) +- i|sin p| deviation (informational): {worstForm:.2e}"
  ok(&"plane waves normalized (defect {worstNorm:.1e})", worstNorm < 1e-12)
  ok(&"plane-wave projection is complete (defect {worstComp:.1e})", worstComp < 1e-10)
  ok("exact spectrum is conjugation symmetric (gamma5)",
     worstConjPair(result) < 1e-10)

proc newDiracOp(w: auto; mass: float; seed: uint64): auto =
  let lo = w.lo
  let s = w.s
  let proto = w.proto
  type F = type(w.proto)
  var startCount = seed
  ArnoldiOp[F](
    apply: proc(r: var F; z: F) =
      s.D(r, z, mass),
    newVec: proc(): F =
      result = newOneOf(proto)
      result := 0,
    start: proc(v: var F) =
      inc startCount
      let salt = startCount*0x1000000'u64
      for i in lo.sites:
        for sp in 0..3:
          for c in 0..2:
            let k = salt + uint64(i)*24 + uint64(sp)*6 + uint64(c)*2
            setC(v{i}[sp][c], u01(k), u01(k+1)))

proc wilsonRun(lo: auto; op: auto; exact: seq[Complex64]; name, which: string;
               nev, ncv, maxRestarts: int): seq[Complex64] =
  let (vals, _, resids, nap) = arnoldi(op, nev, ncv, 1e-10, maxRestarts,
                                       which, verb=1)
  let conv = convOnly(vals, resids, 1e-8)
  var worstR = 0.0
  for r in resids: worstR = max(worstR, r)
  let wn = worstNearest(conv, exact)
  let cp = worstConjPair(conv)
  let sc = shellCoverage(exact, conv, which, nev)
  echo &"  {name} {which}: {conv.len}/{nev} converged, {nap} applies " &
       &"({float(nap)/float(max(1,conv.len)):.0f} per converged pair), worst resid {worstR:.2e}"
  echo &"  worst |ritz - exact| = {wn:.3e}, conj-pairing = {cp:.3e}, shell coverage = {sc:.3e}"
  ok(&"{name} {which}: all {nev} pairs converged (resid<1e-8)", conv.len == nev)
  ok(&"{name} {which}: every converged value matches an exact one to 1e-8 (worst {wn:.1e})",
     wn < 1e-8)
  ok(&"{name} {which}: converged spectrum conjugation symmetric (worst {cp:.1e})",
     cp < 1e-8)
  ok(&"{name} {which}: first-{nev} exact shells all found (worst {sc:.1e})",
     sc < 1e-8)
  conv

proc wilsonTests(lat: seq[int]; V: static int; doSM: bool; smRestarts: int) =
  let name = &"wilson {lat[0]}^4"
  section &"test 2: free Wilson-Dirac, {name}, unit gauge, antiperiodic time, m = 0.3"
  let mass = 0.3
  let w = makeWilson(lat, V)
  let exact = wilsonExact(w, mass)
  var op = newDiracOp(w, mass, 12345'u64)
  discard wilsonRun(w.lo, op, exact, name, "LM", 16, 48, 500)
  if doSM:
    discard wilsonRun(w.lo, op, exact, name, "SM", 16, 48, smRestarts)

# ---------------------------------------------------------------------------
# 3: shift-invert — Arnoldi LM on D^{-1}
# ---------------------------------------------------------------------------

proc shiftInvertTest(lat: seq[int]) =
  section &"test 3: shift-invert (LM on D^-1), wilson {lat[0]}^4, m = 0.3"
  let mass = 0.3
  let w = makeWilson(lat, 4)
  let lo = w.lo
  let s = w.s
  let exact = wilsonExact(w, mass)
  type F = type(w.proto)
  var t3 = newOneOf(w.proto)   # D b, inside the normal operator
  var rhs = newOneOf(w.proto)  # D^dag x
  var chk = newOneOf(w.proto)
  proc normalOp(a: F; b: F) =  # a = D^dag D b, called inside cg's threads
    threadBarrier()
    s.D(t3, b, mass)
    threadBarrier()
    s.Ddag(a, t3, mass)
    threadBarrier()
  var cgop = (apply: normalOp, precon: cpNone)
  var totalCGits = 0
  var nSolve = 0
  var worstInv = 0.0           # direct check |D y - x|/|x|
  var startCount = 777'u64
  var op = ArnoldiOp[F](
    apply: proc(r: var F; x: F) =
      # r = D^{-1} x = (D^dag D)^{-1} D^dag x   (exact for any D)
      s.Ddag(rhs, x, mass)
      r := 0
      var sp = initSolverParams()
      sp.r2req = 1e-24
      sp.maxits = 10000
      sp.verbosity = 0
      sp.subset.layoutSubset(lo, "all")
      cgSolve(r, rhs, cgop, sp)
      totalCGits += sp.finalIterations
      inc nSolve
      if nSolve <= 3:          # verify the inversion directly a few times
        s.D(chk, r, mass)
        chk -= x
        worstInv = max(worstInv, sqrt(chk.norm2/x.norm2)),
    newVec: proc(): F =
      result = newOneOf(t3)
      result := 0,
    start: proc(v: var F) =
      inc startCount
      let salt = startCount*0x1000000'u64
      for i in lo.sites:
        for sp in 0..3:
          for c in 0..2:
            let k = salt + uint64(i)*24 + uint64(sp)*6 + uint64(c)*2
            setC(v{i}[sp][c], u01(k), u01(k+1)))

  let nev = 16
  let (vals, _, resids, nap) = arnoldi(op, nev, 48, 1e-10, 100, "LM", verb=1)
  echo &"  first solves: |D y - x|/|x| = {worstInv:.2e} (CG on normal equations, r2req 1e-24)"
  ok(&"D^-1 apply verified directly ({worstInv:.1e})", worstInv < 1e-9)
  let conv = convOnly(vals, resids, 1e-8)
  var inv: seq[Complex64]
  for v in conv: inv.add complex64(1.0, 0.0)/v
  let wn = worstNearest(inv, exact)
  let sc = shellCoverage(exact, inv, "SM", nev)
  var worstR = 0.0
  for r in resids: worstR = max(worstR, r)
  echo &"  {conv.len}/{nev} converged, worst resid {worstR:.2e}"
  echo &"  worst |1/ritz - exact| = {wn:.3e}, smallest-shell coverage = {sc:.3e}"
  echo &"  cost: {nap} D^-1 applies, {totalCGits} CG iterations total " &
       &"({float(totalCGits)/float(max(1,nSolve)):.1f}/solve, 2 D-applies per CG it + 1 per solve)"
  echo &"  => about {2*totalCGits + nSolve} Wilson D applies, " &
       &"{float(2*totalCGits+nSolve)/float(max(1,conv.len)):.0f} per converged pair"
  ok(&"shift-invert: all {nev} pairs converged", conv.len == nev)
  ok(&"shift-invert: 1/lambda matches exact smallest to 1e-8 (worst {wn:.1e})",
     wn < 1e-8)
  ok(&"shift-invert: smallest exact shells covered (worst {sc:.1e})", sc < 1e-8)

# ---------------------------------------------------------------------------

diagTests()
wilsonTests(@[4, 4, 4, 4], 4, doSM = true, smRestarts = 4000)
wilsonTests(@[6, 6, 6, 6], 1, doSM = true, smRestarts = 6000)
shiftInvertTest(@[4, 4, 4, 4])

section "summary"
if nFail > 0:
  echo "FAILED checks: ", nFail
else:
  echo "all checks passed"
qexFinalize()
quit(if nFail > 0: 1 else: 0)
