## Test of `hcanalysis.nim`: flow-scale finders, jackknife, autocorrelation
## time and polynomial fits on analytic and synthetic data.

import std/[math, strformat, random, unittest]
import ../hcanalysis
import helpers

proc close(a, b, tol: float): bool = abs(a-b) <= tol

suite "hcanalysis":

  test "1. findT0 / findW0 / findCrossing on analytic curves":
    # exactly linear curve: linear interpolation is exact
    block:
      var t, y: seq[float]
      for i in 0..200:
        let tt = 0.01*i.float
        t.add tt
        y.add 0.25*tt            # crosses 0.3 at t = 1.2 exactly
      let t0 = findT0(t, y)
      ok(&"findT0 linear curve, exact crossing: t0 = {t0:.15g} (exact 1.2)",
         close(t0, 1.2, 1e-12))
    # t^2 E = A t^3: crossing of 0.3 at t = (0.3/A)^(1/3); linear
    # interpolation is O(dt^2), cubic O(dt^4)
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
        ok(&"findT0 cubic curve dt={dt} linear: err = {abs(l-exact):.3e}",
           abs(l-exact) < 3.0*dt*dt)
        ok(&"findT0 cubic curve dt={dt} cubic: err = {abs(c-exact):.3e}",
           abs(c-exact) < 1e-9)
    # no crossing -> noCrossing sentinel (not NaN: the build uses -ffast-math,
    # under which NaN comparisons are unreliable)
    block:
      let t = @[0.1, 0.2, 0.3]
      let y = @[0.01, 0.02, 0.03]
      let r = findT0(t, y)
      ok(&"findT0 returns noCrossing when the target is never reached (got {r})",
         r == noCrossing)
    # derivT2E + findW0 on t^2E = A t^3: W(t) = t d/dt(A t^3) = 3A t^3, so
    # W = 0.3 at t = (0.1/A)^(1/3)
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
      ok(&"findW0 (returns w0^2) on t^2E = A t^3: w0^2 = {w0sq:.8f} exact {exact:.8f}",
         close(w0sq, exact, 1e-3))
    block:
      var t, y: seq[float]
      for i in 0..20:
        let tt = 0.1*i.float
        t.add tt
        y.add 0.15*tt*tt
      # W(t) = t d/dt(0.15 t^2) = 0.3 t^2, so w0^2 = 1
      let w0sq = findW0(t, derivT2E(t, y))
      ok(&"findW0 recovers the quadratic crossing: w0^2 = {w0sq:.12f}",
         close(w0sq, 1.0, 1e-12))
    block:
      let t = @[0.0, 0.1, 0.4, 1.1, 1.6, 2.0]
      var y: seq[float]
      for tt in t: y.add 0.15*tt*tt
      let w = derivT2E(t, y)
      var exact = true
      for i in 1..<t.len-1:
        exact = exact and close(w[i], 0.3*t[i]*t[i], 1e-12)
      ok("derivT2E is exact on a nonuniform quadratic grid at interior points", exact)
    block:
      let x = @[1.0, 2.0, 4.0]
      let y = @[3.0, 5.0, 9.0]
      ok("interpAt within a bracket", close(interpAt(x, y, 3.0), 7.0, 1e-14))
      ok("interpAt at and beyond endpoints",
         interpAt(x, y, 1.0) == 3.0 and interpAt(x, y, 0.0) == 3.0 and
         interpAt(x, y, 4.0) == 9.0 and interpAt(x, y, 5.0) == 9.0)
      ok("interpAt empty inputs", interpAt([], [], 1.0) == 0.0)
      ok("interpAt singleton inputs",
         interpAt([2.0], [5.0], 1.0) == 5.0 and interpAt([2.0], [5.0], 3.0) == 5.0)

  test "2. jackknife":
    # jackknife of the mean equals the ordinary standard error exactly
    block:
      var r = initRand(20260821)
      var x: seq[float]
      for i in 0..<137: x.add r.gauss(mu = 2.5, sigma = 0.75)
      let (m, e) = jackknifeMean(x)
      let m0 = x.mean
      let e0 = stderrMean(x)
      ok(&"jackknife mean == sample mean: {m:.12f} vs {m0:.12f}",
         close(m, m0, 1e-14))
      ok(&"jackknife error == sqrt(var/n) for f=mean: {e:.12e} vs {e0:.12e}",
         abs(e-e0) < 1e-14*max(1.0, e0))
    # non-linear estimator: jackknife of <x^2> - <x>^2 against the
    # delta-method error on a large sample (loose, a consistency check)
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
      ok(&"jackknife error of the variance ~ sqrt(2/n): jk {e:.5f} expected ~{expect:.5f} (val {m:.4f})",
         abs(e-expect) < 0.25*expect)
    # binned jackknife on a block-correlated series: with bin equal to the
    # block length the error grows to the block-level error
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
      ok(&"binned jackknife inflates the error on correlated data by ~sqrt(bin): e1 = {e1:.5f} e8 = {e8:.5f} ratio {e8/e1:.3f}",
         e8/e1 > 2.0 and e8/e1 < 4.0)
    # a trailing partial block is not exchangeable with the full blocks; for
    # the mean the unequal-delete pseudovalues are the deleted block means
    block:
      let x = @[1.0, 2.0, 4.0, 8.0, 16.0]  # block means 1.5, 6, 16 for bin=2
      let (m, e) = jackknifeMean(x, 2)
      let
        jm = x.mean
        expectVar = ((2.0/3.0)*(1.5-jm)^2 +
                     (2.0/3.0)*(6.0-jm)^2 +
                     0.25*(16.0-jm)^2)/3.0
        expect = sqrt(expectVar)
      ok(&"unequal-delete jackknife weights a partial block: mean {m:.12f}, err {e:.12f}, expected {expect:.12f}",
         close(m, jm, 1e-14) and close(e, expect, 1e-14))

  test "3. autocorrTime":
    block:
      var short = true
      for n in 0..3:
        short = short and autocorrTimeW(newSeq[float](n)) == (0.5, 0.0, 0)
      ok("autocorrTimeW for 0..3 samples", short)
      let x = @[2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
      let (tau, err, w) = autocorrTimeW(x)
      ok("autocorrTimeW constant series",
         tau == 0.5 and w == 3 and close(err, 0.5*sqrt(14.0/8.0), 1e-14))
    block:
      let x = @[1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0]
      ok("autocorr alternating series", autocorr(x) == @[1.0, -1.0, 1.0, -1.0, 1.0])
      let (tau, err, w) = autocorrTimeW(x)
      ok("autocorrTimeW clamps negative contributions before closing the window",
         tau == 0.5 and w == 3 and close(err, 0.5*sqrt(14.0/8.0), 1e-14))
      let (tm, em, wm) = autocorrTimeW(x, 100.0)
      ok("autocorrTimeW stops at half the series length when the window stays open",
         tm == 1.5 and wm == 4 and close(em, 2.25, 1e-14))
    # uncorrelated white noise -> tau ~ 0.5
    block:
      var r = initRand(5)
      var x: seq[float]
      for i in 0..<20000: x.add r.gauss()
      let tau = autocorrTime(x)
      ok(&"autocorrTime of white noise ~ 0.5: tau = {tau:.4f}", abs(tau-0.5) < 0.15)
    # AR(1) with x_{n+1} = a x_n + noise -> tau_int = 1/2 (1+a)/(1-a)
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
      ok(&"autocorrTime of AR(1) a=0.8: tau = {tau:.4f} exact {exact:.4f}",
         abs(tau-exact) < 0.15*exact)

  test "4. fitPoly":
    # noiseless synthetic data: exact recovery of the coefficients
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
      ok(&"fitPoly [0,1] recovers exact coefficients: c0 = {co[0]:.12f} c1 = {co[1]:.12f} chi2/dof = {cd:.3e}",
         close(co[0], c0, 1e-10) and close(co[1], c1, 1e-10))
      ok(&"fitPoly chi^2 = 0 on exact data: chi2/dof = {cd:.3e}", cd < 1e-18)
      ok(&"fitPoly error on the intercept is positive and sane: dc0 = {er[0]:.5f}",
         er[0] > 0.0 and er[0] < 1.0)
    # quadratic (O(a^4)) form, exact recovery
    block:
      const c = [6.78, -0.05, -0.55]
      var x, y, dy: seq[float]
      for i in 0..6:
        let xx = 0.05 + 0.3*i.float
        x.add xx
        y.add c[0] + c[1]*xx + c[2]*xx*xx
        dy.add 0.02
      let (co, _, cd) = fitPoly(x, y, dy, [0, 1, 2])
      var exact = true
      for k in 0..2:
        if not close(co[k], c[k], 1e-8): exact = false
      ok(&"fitPoly [0,1,2] recovers exact coefficients: c = {co[0]:.10f} {co[1]:.10f} {co[2]:.10f} chi2/dof = {cd:.2e}",
         exact)
    # Gaussian noise: chi^2/dof ~ 1 and 68% coverage of the intercept
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
      ok(&"fitPoly chi^2/dof ~ 1 with Gaussian noise: <chi2/dof> = {chis.mean:.4f}",
         abs(chis.mean - 1.0) < 0.15)
      ok(&"fitPoly 1-sigma coverage of the intercept ~ 68%: coverage = {frac*100.0:.1f}%",
         abs(frac-0.6827) < 0.06)
    # weighting matters: a badly measured outlier is ignored
    block:
      var x = @[0.1, 0.2, 0.3, 0.4]
      var y = @[6.5, 6.3, 6.1, 99.0]
      var dy = @[0.02, 0.02, 0.02, 1.0e6]
      let (co, _, _) = fitPoly(x, y, dy, [0, 1])
      ok(&"fitPoly downweights a huge-error point: c0 = {co[0]:.5f} c1 = {co[1]:.5f}",
         close(co[0], 6.7, 0.02) and close(co[1], -2.0, 0.2))
    block:
      let co = @[6.65, -1.37]
      ok("evalPoly", close(evalPoly(co, [0, 1], 0.2), 6.65-1.37*0.2, 1e-14))
