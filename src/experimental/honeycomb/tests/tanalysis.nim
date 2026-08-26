## Standalone test of `hcanalysis.nim`.  Prints PASS/FAIL lines, exits non-zero
## on any failure.  No QEX dependency.
##
##   make src/experimental/honeycomb/tests/tanalysis.nim && ./bin/tanalysis

import std/[math, strformat, random]
import ../hcanalysis

var nfail = 0
var ntest = 0

proc check(name: string, ok: bool, info = "") =
  inc ntest
  if ok:
    echo "PASS  ", name, (if info.len > 0: "   " & info else: "")
  else:
    inc nfail
    echo "FAIL  ", name, (if info.len > 0: "   " & info else: "")

proc close(a, b, tol: float): bool = abs(a-b) <= tol

# --------------------------------------------------------------------------
echo "== findT0 / findW0 / findCrossing on analytic curves =="

# 1. Exactly linear curve: linear interpolation must be exact.
block:
  var t, y: seq[float]
  for i in 0..200:
    let tt = 0.01*i.float
    t.add tt
    y.add 0.25*tt            # crosses 0.3 at t = 1.2 exactly
  let t0 = findT0(t, y)
  check("findT0 linear curve, exact crossing",
        close(t0, 1.2, 1e-12), &"t0 = {t0:.15g} (exact 1.2)")

# 2. t^2 E = A t^3  (a smooth analytic curve).  Crossing of 0.3 at
#    t = (0.3/A)^(1/3).  Linear interpolation is O(dt^2); cubic is O(dt^4).
block:
  const A = 0.11
  let exact = pow(0.3/A, 1.0/3.0)
  for dt in [0.05, 0.025]:
    var t, y: seq[float]
    var i = 0
    while 0.5 + dt*i.float <= 3.0:
      let tt = 0.5 + dt*i.float
      t.add tt
      y.add A*tt*tt*tt
      inc i
    let l = findT0(t, y, 0.3, 1)
    let c = findT0(t, y, 0.3, 3)
    check(&"findT0 cubic curve dt={dt} linear",
          abs(l-exact) < 3.0*dt*dt, &"err = {abs(l-exact):.3e}")
    check(&"findT0 cubic curve dt={dt} cubic ",
          abs(c-exact) < 1e-9, &"err = {abs(c-exact):.3e}")

# 3. Sanity: no crossing -> noCrossing sentinel (not NaN: the build uses
#    -ffast-math, under which NaN comparisons are unreliable).
block:
  let t = @[0.1, 0.2, 0.3]
  let y = @[0.01, 0.02, 0.03]
  let r = findT0(t, y)
  check("findT0 returns noCrossing when target is never reached",
        r == noCrossing, &"got {r}")

# 4. derivT2E + findW0 on t^2E = A t^3: W(t) = t d/dt(A t^3) = 3A t^3.
#    So W = 0.3 at t = (0.1/A)^(1/3).
block:
  const A = 0.11
  let exact = pow(0.1/A, 1.0/3.0)
  var t, y: seq[float]
  var i = 0
  while 0.3 + 0.01*i.float <= 3.0:
    let tt = 0.3 + 0.01*i.float
    t.add tt
    y.add A*tt*tt*tt
    inc i
  let w = derivT2E(t, y)
  let w0sq = findW0(t, w, 0.3)
  check("findW0 (returns w0^2) on t^2E = A t^3",
        close(w0sq, exact, 1e-3), &"w0^2 = {w0sq:.8f} exact {exact:.8f}")

# --------------------------------------------------------------------------
echo ""
echo "== jackknife =="

# 5. jackknife of the mean must equal the ordinary standard error, exactly.
block:
  var r = initRand(20260821)
  var x: seq[float]
  for i in 0..<137: x.add r.gauss(mu = 2.5, sigma = 0.75)
  let (m, e) = jackknifeMean(x)
  let m0 = x.mean
  let e0 = stderrMean(x)
  check("jackknife mean == sample mean",
        close(m, m0, 1e-14), &"{m:.12f} vs {m0:.12f}")
  check("jackknife error == sqrt(var/n) for f=mean",
        abs(e-e0) < 1e-14*max(1.0, e0), &"{e:.12e} vs {e0:.12e}")

# 6. Non-linear estimator: jackknife of <x^2> - <x>^2 against the delta-method
#    error on a large sample (loose tolerance, this is a consistency check).
block:
  var r = initRand(11)
  var x: seq[float]
  for i in 0..<4000: x.add r.gauss(mu = 0.0, sigma = 1.0)
  proc v(s: openArray[float]): float =
    var m = 0.0
    var m2 = 0.0
    for a in s:
      m += a
      m2 += a*a
    m = m/s.len.float
    m2 = m2/s.len.float
    m2 - m*m
  let (m, e) = jackknife(x, v)
  # for a Gaussian, var(sigma^2_hat) = 2 sigma^4/n
  let expect = sqrt(2.0/4000.0)
  check("jackknife error of the variance ~ sqrt(2/n)",
        abs(e-expect) < 0.25*expect, &"jk {e:.5f} expected ~{expect:.5f} (val {m:.4f})")

# 7. Binned jackknife on an artificially block-correlated series: with bin
#    equal to the block length, the error must grow to the block-level error.
block:
  var r = initRand(777)
  var x: seq[float]
  const nb = 200
  const bl = 8
  for i in 0..<nb:
    let c = r.gauss(mu = 1.0, sigma = 1.0)
    for j in 0..<bl: x.add c + 0.001*r.gauss()
  let e1 = jackknifeMean(x, 1).err
  let e8 = jackknifeMean(x, bl).err
  check("binned jackknife inflates error on correlated data by ~sqrt(bin)",
        e8/e1 > 2.0 and e8/e1 < 4.0, &"e1 = {e1:.5f} e8 = {e8:.5f} ratio {e8/e1:.3f}")

# A trailing partial block is not exchangeable with the full blocks.  For the
# mean, unequal-delete pseudovalues are exactly the deleted block means, making
# this small example independently calculable.
block:
  let x = @[1.0, 2.0, 4.0, 8.0, 16.0]  # block means 1.5, 6, 16 for bin=2
  let (m, e) = jackknifeMean(x, 2)
  let
    jm = x.mean
    expectVar = ((2.0/3.0)*(1.5-jm)^2 +
                 (2.0/3.0)*(6.0-jm)^2 +
                 0.25*(16.0-jm)^2)/3.0
    expect = sqrt(expectVar)
  check("unequal-delete jackknife weights a partial block",
        close(m, jm, 1e-14) and close(e, expect, 1e-14),
        &"mean {m:.12f}, err {e:.12f}, expected {expect:.12f}")

# --------------------------------------------------------------------------
echo ""
echo "== autocorrTime =="

# 8. Uncorrelated white noise -> tau ~ 0.5.
block:
  var r = initRand(5)
  var x: seq[float]
  for i in 0..<20000: x.add r.gauss()
  let tau = autocorrTime(x)
  check("autocorrTime of white noise ~ 0.5",
        abs(tau-0.5) < 0.15, &"tau = {tau:.4f}")

# 9. AR(1) with x_{n+1} = a x_n + noise -> tau_int = 1/2 (1+a)/(1-a).
block:
  const a = 0.8
  let exact = 0.5*(1.0+a)/(1.0-a)
  var r = initRand(6)
  var x: seq[float]
  var v = 0.0
  for i in 0..<400000:
    v = a*v + r.gauss()
    if i >= 1000: x.add v
  let tau = autocorrTime(x)
  check("autocorrTime of AR(1) a=0.8",
        abs(tau-exact) < 0.15*exact, &"tau = {tau:.4f} exact {exact:.4f}")

# --------------------------------------------------------------------------
echo ""
echo "== fitPoly =="

# 10. Noiseless synthetic data: exact recovery of the coefficients.
block:
  const c0 = 6.65
  const c1 = -1.37
  var x, y, dy: seq[float]
  for i in 0..4:
    let xx = 0.10 + 0.06*i.float
    x.add xx
    y.add c0 + c1*xx
    dy.add 0.05
  let (co, er, cd) = fitPoly(x, y, dy, [0, 1])
  check("fitPoly [0,1] recovers exact coefficients",
        close(co[0], c0, 1e-10) and close(co[1], c1, 1e-10),
        &"c0 = {co[0]:.12f} c1 = {co[1]:.12f} chi2/dof = {cd:.3e}")
  check("fitPoly chi^2 = 0 on exact data", cd < 1e-18, &"chi2/dof = {cd:.3e}")
  check("fitPoly error on the intercept is positive and sane",
        er[0] > 0.0 and er[0] < 1.0, &"dc0 = {er[0]:.5f}")

# 11. Quadratic (O(a^4)) form, exact recovery.
block:
  const c = [6.78, -0.05, -0.55]
  var x, y, dy: seq[float]
  for i in 0..6:
    let xx = 0.05 + 0.3*i.float
    x.add xx
    y.add c[0] + c[1]*xx + c[2]*xx*xx
    dy.add 0.02
  let (co, _, cd) = fitPoly(x, y, dy, [0, 1, 2])
  var ok = true
  for k in 0..2:
    if not close(co[k], c[k], 1e-8): ok = false
  check("fitPoly [0,1,2] recovers exact coefficients", ok,
        &"c = {co[0]:.10f} {co[1]:.10f} {co[2]:.10f} chi2/dof = {cd:.2e}")

# 12. Gaussian noise: chi^2/dof ~ 1 and the intercept is within ~1-2 sigma of
#     the truth; also check the error scales as 1/sqrt(nrepeat).
block:
  const c0 = 6.65
  const c1 = -1.37
  const sig = 0.04
  var r = initRand(4242)
  var nIn = 0
  const nrep = 400
  var chis: seq[float]
  for rep in 0..<nrep:
    var x, y, dy: seq[float]
    for i in 0..5:
      let xx = 0.10 + 0.05*i.float
      x.add xx
      y.add c0 + c1*xx + r.gauss(mu = 0.0, sigma = sig)
      dy.add sig
    let (co, er, cd) = fitPoly(x, y, dy, [0, 1])
    chis.add cd
    if abs(co[0]-c0) < er[0]: inc nIn
  let frac = nIn.float/nrep.float
  check("fitPoly chi^2/dof ~ 1 with Gaussian noise",
        abs(chis.mean - 1.0) < 0.15, &"<chi2/dof> = {chis.mean:.4f}")
  check("fitPoly 1-sigma coverage of the intercept ~ 68%",
        abs(frac-0.6827) < 0.06, &"coverage = {frac*100.0:.1f}%")

# 13. Weighting actually matters: a badly-measured outlier must be ignored.
block:
  var x = @[0.1, 0.2, 0.3, 0.4]
  var y = @[6.5, 6.3, 6.1, 99.0]
  var dy = @[0.02, 0.02, 0.02, 1.0e6]
  let (co, _, _) = fitPoly(x, y, dy, [0, 1])
  check("fitPoly downweights a huge-error point",
        close(co[0], 6.7, 0.02) and close(co[1], -2.0, 0.2),
        &"c0 = {co[0]:.5f} c1 = {co[1]:.5f}")

# 14. evalPoly agrees with the fit at the data points.
block:
  let co = @[6.65, -1.37]
  check("evalPoly", close(evalPoly(co, [0, 1], 0.2), 6.65-1.37*0.2, 1e-14))

# --------------------------------------------------------------------------
echo ""
echo "== topoQ normalisation constant =="
# Verdict (see hcanalysis.nim and doc/RESULTS_CUBIC.md): QEX's topoQ is right,
# so the fix factor is 1.  The *numerical* proof lives in
# `refCubicMeas -abeliantest`, which needs a QEX layout and so cannot run here.
check("qexTopoQNormFix == 1 (QEX topoQ needs no correction)",
      qexTopoQNormFix == 1.0)
check("fixTopoQ is the identity", close(fixTopoQ(0.37), 0.37, 1e-15))

# --------------------------------------------------------------------------
echo ""
echo &"{ntest-nfail}/{ntest} tests passed"
if nfail > 0:
  echo "FAILED"
  quit(1)
echo "ALL PASS"
