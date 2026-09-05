#RUNCMD $RUN1

import std/[math, os, random, tables, unittest]
import ../meas/[fit, dataio]
import utils/resample

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

template note(s: varargs[string, `$`]) =
  ## Achieved accuracy, echoed so the numbers quoted in doc/06-status.md are checkable.
  var m = "  measured:"
  for x in s: m.add " " & x
  echo m

# --- (V.4)-(V.5) --------------------------------------------------------------

suite "effMass":
  test "exact cosh recovers Delta":
    # c(t) = cosh(D (T/2 - t))  =>  f(t) = D (T/2 - t)  =>  Delta_eff = D exactly
    for cs in [(1.0, 16.0, 168), (0.5, 12.0, 120), (2.0, 8.0, 64), (sqrt 2.0, 16.0, 96)]:
      let
        (d, tt, nt) = cs
        at = tt/float(nt)
      var c = newSeq[float](nt+1)
      for i in 0..nt: c[i] = cosh(d*(0.5*tt - at*float(i)))
      let m = effMass(c, at, tt)
      check m.len == nt
      var worst = 0.0
      for i in 0..<(nt div 2): worst = max(worst, abs(m[i] - d))   # t + at <= T/2
      note "effMass D=", d, " T=", tt, " Lt=", nt, " worst|dD|=", worst
      check worst < 1e-12

  test "reference point is c(T/2), not c(t+at)":
    # on a pure exponential the log-ratio mass is exactly d everywhere, while (V.5)
    # runs away near T/2 because arccosh(y) is not ln(2y) there.  Pins the convention.
    const
      d = 1.3
      at = 0.1
      tt = 8.0
      nt = 80
    var c = newSeq[float](nt+1)
    for i in 0..nt: c[i] = exp(-d*at*float(i))
    let m = effMass(c, at, tt)
    proc want(i: int): float =
      -(arccosh(exp(d*(0.5*tt - at*float(i+1)))) -
        arccosh(exp(d*(0.5*tt - at*float(i)))))/at
    var worst = 0.0
    for i in 0..<(nt div 2): worst = max(worst, abs(m[i] - want(i)))
    note "effMass vs closed-form (V.5) worst=", worst, " m[0]=", m[0], " m[38]=", m[38]
    check worst < 1e-10
    check abs(m[0] - d) < 1e-3         # far from T/2 the two conventions agree
    check abs(m[38] - d) > 0.5         # two steps from T/2 they do not

# --- (V.6) --------------------------------------------------------------------

suite "plateauFit":
  test "constant data has an unidentified exponential gap":
    let f = plateauFit([2.0, 2, 2, 2, 2, 2, 2, 2, 2, 2], 0, 10,
                       e = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
    check f.status == fitSingular

  test "insufficient data and iteration exhaustion are explicit":
    check plateauFit([1.0, 1.1, 1.2], 0, 3).status == fitShort
    var m = newSeq[float](10)
    for i in 0..<m.len: m[i] = 1.0 + 0.8*exp(-0.9*float(i))
    let f = plateauFit(m, 0, m.len, maxit = 0)
    check f.status == fitLimit
    check f.iters == 0

  test "noiseless recovery of (Delta_0, c, Delta')":
    # (d0, c, dp, at, t0, t1); the first row is the paper's window 4 <= t < 8, T=16, Lt=168
    for cs in [(1.0, 0.5, 1.3, 16.0/168.0, 42, 84),
               (1.41, -0.7, 0.9, 0.1, 0, 41),
               (0.953918, 2.0, 2.5, 0.2, 5, 40)]:
      let (d0, c, dp, at, t0, t1) = cs
      var m = newSeq[float](t1)
      for i in 0..<t1: m[i] = d0 + c*exp(-dp*at*float(i))
      let f = plateauFit(m, t0, t1, at)
      note "plateauFit t=[", at*float(t0), ",", at*float(t1), ") |dd0|=", abs(f.d0 - d0),
        " |dc/c|=", abs(f.c/c - 1.0), " |ddp|=", abs(f.dp - dp), " iters=", f.iters
      check f.status == fitOk
      check f.dof == t1 - t0 - 3
      check abs(f.d0 - d0) < 1e-9
      check abs(f.c - c) < 1e-9*abs(c)
      check abs(f.dp - dp) < 1e-9

  test "quoted error is the right size under known Gaussian noise":
    const
      d0 = 1.0
      c = 0.8
      dp = 0.9
      at = 0.1
      t1 = 41
      sig = 2e-3
      nrep = 400
    var
      r = initRand 20260821
      sp, sp2 = 0.0
      nbad = 0
      m = newSeq[float](t1)
      e = newSeq[float](t1)
    for i in 0..<t1: e[i] = sig
    for k in 0..<nrep:
      for i in 0..<t1: m[i] = d0 + c*exp(-dp*at*float(i)) + r.gauss(0.0, sig)
      let f = plateauFit(m, 0, t1, at, e)
      if f.status != fitOk: inc nbad
      let p = (f.d0 - d0)/f.ed0
      sp += p
      sp2 += p*p
    let
      mp = sp/float(nrep)
      rp = sqrt(sp2/float(nrep))
    note "plateauFit pull: mean=", mp, " rms=", rp, " over ", nrep, " noise draws, nonconverged=", nbad
    check nbad == 0
    check abs(mp) < 3.0/sqrt(float nrep)     # unbiased
    check abs(rp - 1.0) < 0.12               # error scale correct

  test "chi2/dof of a clean fit is O(1)":
    var
      r = initRand 777
      m = newSeq[float](41)
      e = newSeq[float](41)
    for i in 0..<41:
      m[i] = 1.0 + 0.8*exp(-0.9*0.1*float(i)) + r.gauss(0.0, 2e-3)
      e[i] = 2e-3
    let f = plateauFit(m, 0, 41, 0.1, e)
    note "plateauFit chi2/dof=", f.chi2dof
    check f.chi2dof > 0.4
    check f.chi2dof < 2.0

suite "local dimensions":
  test "logarithm availability comes from its input domain":
    let
      a = effLocal([4.0, 2.0, 1.0], 0.5)
      b = effLocal([-4.0, -2.0, -1.0], 0.5)
      c = effLocal([4.0, 0.0, -1.0, 2.0], 0.5)
      g = effLocal([-4.0, -2.0, -1.0], 0.5, positive = true)
    for i in 0..1:
      check a[i].ok and b[i].ok
      check abs(a[i].v - 2.0*ln(2.0)) < 1e-14
      check a[i].v == b[i].v
      check not g[i].ok
    for x in c: check not x.ok

  test "plateau selection uses the longest available run and keeps the first tie":
    var m = newSeq[Estimate](11)
    for i in 0..<m.len:
      let d = if i < 5: 1.2 else: 2.3
      m[i] = Estimate(v: d + 0.7*exp(-0.8*0.2*float(i)), ok: i != 5)
    let a = deltaFit(m, 0, m.len, 0.2)
    check a.ok
    check abs(a.v - 1.2) < 1e-9
    m[1].ok = false
    let b = deltaFit(m, 0, m.len, 0.2)
    check b.ok
    check abs(b.v - 2.3) < 1e-9
    check not deltaFit(m, 0, 5, 0.2).ok

suite "configuration jackknife":
  test "unequal groups use the shared pseudo-value convention":
    let
      xs = @[0.0, 0.0, 0.0, 0.0, 100.0]
      cs = @[@[0.0], @[0.0], @[0.0], @[0.0], @[100.0]]
      ids = @[10, 20, 30, 40, 50]
    proc mean(xs: Ensemble[seq[float]]): float =
      for i in 0..<xs.len: result += xs[i]/float(xs.len)
    for bs in [1, 2]:
      let j = jkFrom(cs, proc(c: seq[float]): float = c[0], ids, bs)
      let got = jkStat(j)
      let want = xs.jackknife(bs, mean)
      check got.v.ok and got.e.ok
      check abs(got.v.v - want.mean) < 1e-13
      check abs(got.e.v - want.stdev) < 1e-13
      for i in 0..<j.reps.len:
        check j.reps[i].ok
        check abs(j.reps[i].v - want.jksamples[i]) < 1e-13
    let st = jkStat(jkFrom(cs, proc(c: seq[float]): float = c[0], ids, 2))
    check abs(st.e.v - 80.0/3.0) < 1e-13

  test "one group evaluates only the full sample and has no error":
    for cs in [@[@[3.0]], @[@[2.0], @[4.0]]]:
      var calls = 0
      var ids: seq[int]
      for i in 0..<cs.len: ids.add i
      let j = jkFrom(cs, proc(c: seq[float]): float =
        inc calls
        c[0], ids, 3)
      let st = jkStat(j)
      check calls == 1
      check j.reps.len == 0
      check st.v.ok and st.v.v == 3.0
      check not st.e.ok
    let st = jack([3.0], bs = 1)
    check st.mean == 3.0
    check not st.hasErr

  test "failed replicas keep their positions and invalidate uncertainty":
    let cs = @[@[1.0], @[2.0], @[3.0], @[4.0]]
    let ids = @[1, 2, 3, 4]
    let a = jkFrom(cs, proc(c: seq[float]): Estimate =
      Estimate(v: c[0], ok: c[0] != 3.0), ids)
    let b = jkFrom(cs, proc(c: seq[float]): float = c[0], ids)
    let st = jkStat(a)
    check a.reps.len == 4
    check not a.reps[0].ok
    for i in 1..3: check a.reps[i].ok
    check st.v.ok and st.v.v == 2.5
    check not st.e.ok
    let av = jkStat(jkMean([a, b]))
    check av.v.ok and av.v.v == 2.5
    check not av.e.ok
    let rat = ratioVE(a, b)
    check rat.v.ok and rat.v.v == 1.0
    check not rat.e.ok
    let bad = jkFrom(cs, proc(c: seq[float]): Estimate =
      Estimate(v: c[0], ok: c[0] != 2.5), ids)
    check not jkStat(bad).v.ok
    check not jkStat(bad).e.ok

  test "ratios pair groups only when sample identities and blocks match":
    let cs = @[@[1.0], @[2.0], @[4.0], @[8.0], @[16.0]]
    let ids = @[1, 2, 3, 4, 5]
    let a = jkFrom(cs, proc(c: seq[float]): float = 2.0*c[0], ids, 2)
    let b = jkFrom(cs, proc(c: seq[float]): float = c[0], ids, 2)
    let rat = ratioVE(a, b)
    check rat.v.ok and rat.e.ok
    check rat.v.v == 2.0
    check rat.e.v == 0.0
    let dis = jkFrom(cs, proc(c: seq[float]): float = c[0], @[6, 7, 8, 9, 10], 2)
    let dr = ratioVE(a, dis)
    let sa = jkStat(a)
    let sb = jkStat(dis)
    let want = 2.0*sqrt((sa.e.v/sa.v.v)^2 + (sb.e.v/sb.v.v)^2)
    check dr.v.ok and dr.e.ok
    check abs(dr.e.v - want) < 1e-13
    for other in [jkFrom(cs, proc(c: seq[float]): float = c[0], ids, 1),
                  jkFrom(cs, proc(c: seq[float]): float = c[0], @[2, 1, 3, 4, 5], 2),
                  jkFrom(cs, proc(c: seq[float]): float = c[0], @[5, 6, 7, 8, 9], 2)]:
      let r = ratioVE(a, other)
      check r.v.ok and r.v.v == 2.0
      check not r.e.ok
      expect ValueError: discard jkMean([a, other])

# --- (V.7), one dimension -----------------------------------------------------

suite "contFit":
  test "exact line to 1e-14":
    let
      x = [0.0, 0.010, 0.040, 0.090, 0.160, 0.250]
      a = 1.4142135623730951
      b = -0.375
    var y, e: array[6, float]
    for i in 0..5:
      y[i] = a + b*x[i]
      e[i] = 0.001*(1.0 + 0.3*float(i))
    let f = contFit(x, y, e)
    note "contFit |da|=", abs(f.a - a), " |db|=", abs(f.b - b), " chi2=", f.chi2
    check abs(f.a - a) < 1e-14
    check abs(f.b - b) < 1e-14
    check f.chi2 < 1e-20
    check f.dof == 4

  test "errors scale with the input errors":
    let x = [0.019, 0.076, 0.30, 1.23]
    var y, e1, e4: array[4, float]
    for i in 0..3:
      y[i] = 1.0 - 0.2*x[i] + 0.001*float(i*i - 2)
      e1[i] = 0.002
      e4[i] = 4.0*e1[i]
    let
      f1 = contFit(x, y, e1)
      f4 = contFit(x, y, e4)
    check abs(f4.a - f1.a) < 1e-14
    check abs(f4.ea/f1.ea - 4.0) < 1e-12
    check abs(f4.eb/f1.eb - 4.0) < 1e-12
    check abs(f4.cab/f1.cab - 16.0) < 1e-12
    check abs(f4.chi2/f1.chi2 - 1.0/16.0) < 1e-12

# --- (V.7), the actual two-dimensional ansatz ---------------------------------

suite "contFit2":
  const
    asv = [1.10715, 0.55, 0.276, 0.138]
    atv = [0.20, 0.15, 0.10, 0.08, 0.0666667]
    a = 0.999998
    cs = -0.043
    ct = 0.117

  proc grid(): tuple[as2, at2, y: seq[float]] =
    for s in asv:
      for t in atv:
        result.as2.add s*s
        result.at2.add t*t
        result.y.add a + cs*s*s + ct*t*t

  test "exact plane on a 4x5 grid to 1e-12":
    let g = grid()
    var e = newSeq[float](g.y.len)
    for i in 0..<e.len: e[i] = 1e-4*(1.0 + 0.05*float(i))
    let f = contFit2(g.as2, g.at2, g.y, e)
    note "contFit2 |da|=", abs(f.a - a), " |dcs|=", abs(f.cs - cs), " |dct|=", abs(f.ct - ct), " chi2=", f.chi2
    check f.dof == 17
    check abs(f.a - a) < 1e-12
    check abs(f.cs - cs) < 1e-12
    check abs(f.ct - ct) < 1e-12
    check f.chi2 < 1e-16
    # covariance is symmetric and consistent with the quoted errors
    check abs(f.cov[0][1] - f.cov[1][0]) < 1e-30
    check abs(f.cov[0][2] - f.cov[2][0]) < 1e-30
    check abs(f.cov[1][2] - f.cov[2][1]) < 1e-30
    check abs(f.ea - sqrt(f.cov[0][0])) < 1e-30
    check abs(f.ecs - sqrt(f.cov[1][1])) < 1e-30
    check abs(f.ect - sqrt(f.cov[2][2])) < 1e-30

  test "errors scale with the input errors":
    let g = grid()
    var e1, e7 = newSeq[float](g.y.len)
    for i in 0..<e1.len:
      e1[i] = 1e-4*(1.0 + 0.05*float(i))
      e7[i] = 7.0*e1[i]
    let
      f1 = contFit2(g.as2, g.at2, g.y, e1)
      f7 = contFit2(g.as2, g.at2, g.y, e7)
    check abs(f7.a - f1.a) < 1e-12
    check abs(f7.ea/f1.ea - 7.0) < 1e-12
    check abs(f7.ecs/f1.ecs - 7.0) < 1e-12
    check abs(f7.ect/f1.ect - 7.0) < 1e-12
    for i in 0..2:
      for j in 0..2:
        check abs(f7.cov[i][j] - 49.0*f1.cov[i][j]) <= 1e-12*abs(49.0*f1.cov[i][j])

  test "a rank deficient design is rejected, not silently fitted":
    # every lattice at the same a_t: the plane is not determined
    let
      as2 = [1.0, 0.25, 0.0625, 0.015625]
      at2 = [0.04, 0.04, 0.04, 0.04]
      y = [0.95, 0.99, 0.999, 0.9999]
      e = [1e-3, 1e-3, 1e-3, 1e-3]
    expect ValueError:
      discard contFit2(as2, at2, y, e)

  test "noise gives a sensible chi2/dof and pull":
    let g = grid()
    var
      r = initRand 4242
      e = newSeq[float](g.y.len)
      yn = newSeq[float](g.y.len)
      sp2 = 0.0
    for i in 0..<e.len: e[i] = 2e-5
    const nrep = 400
    for k in 0..<nrep:
      for i in 0..<yn.len: yn[i] = g.y[i] + r.gauss(0.0, e[i])
      let f = contFit2(g.as2, g.at2, yn, e)
      let p = (f.a - a)/f.ea
      sp2 += p*p
    let rp = sqrt(sp2/float(nrep))
    note "contFit2 rms pull=", rp, " over ", nrep, " noise draws"
    check abs(rp - 1.0) < 0.12

# --- (V.9) --------------------------------------------------------------------

suite "nmaxFit":
  # (V.3): G(t; nmax) = sum_{n=0..nmax} (n+1) exp(-(n+1) t) / (4 pi)
  proc fermion(ts: seq[float], nmax: int): seq[float] =
    result = newSeq[float](ts.len)
    for i, t in ts:
      var s = 0.0
      for n in 0..nmax: s += float(n+1)*exp(-float(n+1)*t)
      result[i] = s/(4.0*PI)

  # (V.14): G_g(t; nmax) = sum_{n=1..nmax} sqrt(n(n+1)) (2n+1) exp(-sqrt(n(n+1)) t) / (8 pi)
  proc gauge(ts: seq[float], nmax: int): seq[float] =
    result = newSeq[float](ts.len)
    for i, t in ts:
      var s = 0.0
      for n in 1..nmax:
        let l = sqrt(float(n)*float(n+1))
        s += l*float(2*n+1)*exp(-l*t)
      result[i] = s/(8.0*PI)

  test "fermion tower: exact integer nmax and normalization":
    let
      at = 12.0/168.0
      nrm = 3.75
    var ts = newSeq[float](84)
    for i in 0..<84: ts[i] = at*float(i+1)      # t = at .. 6, t=0 excluded (divergent)
    for ntrue in [6, 10, 19]:
      var
        r = initRand(1000 + ntrue)
        g = fermion(ts, ntrue)
      for i in 0..<g.len: g[i] *= nrm*(1.0 + 1e-4*r.gauss(0.0, 1.0))
      let f = nmaxFit(g, proc(nm: int): seq[float] = fermion(ts, nm), 1..40)
      note "nmaxFit fermion true=", ntrue, " fit=", f.nmax, " |dC/C|=", abs(f.c/nrm - 1.0),
        " res/dof=", f.resDof
      check f.nmax == ntrue
      check abs(f.c/nrm - 1.0) < 1e-4
      check f.dof == 82
      check f.resDof < 1e-7                     # ~ (1e-4)^2

  test "gauge tower: exact integer nmax and normalization":
    let
      at = 12.0/120.0
      nrm = 0.125
    var ts = newSeq[float](60)
    for i in 0..<60: ts[i] = at*float(i+1)
    for ntrue in [3, 8, 18]:
      var
        r = initRand(2000 + ntrue)
        g = gauge(ts, ntrue)
      for i in 0..<g.len: g[i] *= nrm*(1.0 + 1e-4*r.gauss(0.0, 1.0))
      let f = nmaxFit(g, proc(nm: int): seq[float] = gauge(ts, nm), 1..40)
      note "nmaxFit gauge true=", ntrue, " fit=", f.nmax, " |dC/C|=", abs(f.c/nrm - 1.0),
        " res/dof=", f.resDof
      check f.nmax == ntrue
      check abs(f.c/nrm - 1.0) < 1e-4
      check f.dof == 58
      check f.resDof < 1e-7

  test "explicit errors switch off the relative weighting":
    let at = 0.1
    var ts = newSeq[float](40)
    for i in 0..<40: ts[i] = at*float(i+1)
    let g = fermion(ts, 8)
    var e = newSeq[float](40)
    for i in 0..<40: e[i] = 1e-6*g[i]
    let f = nmaxFit(g, proc(nm: int): seq[float] = fermion(ts, nm), 1..30, e)
    check f.nmax == 8
    check abs(f.c - 1.0) < 1e-10

# --- jackknife / autocorrelation ---------------------------------------------

proc naiveErr(x: seq[float]): float =
  var m = 0.0
  for a in x: m += a
  m /= float(x.len)
  var s2 = 0.0
  for a in x: s2 += (a - m)*(a - m)
  sqrt(s2/(float(x.len - 1)*float(x.len)))

suite "jack":
  test "AR(1) alpha=0.8 recovers 2 tau = 9":
    const
      alpha = 0.8
      n = 32768
    var
      r = initRand 31415
      v = 0.0
      x = newSeq[float](n)
    for k in 0..<4000: v = alpha*v + r.gauss(0.0, 1.0)     # burn in
    for i in 0..<n:
      v = alpha*v + r.gauss(0.0, 1.0)
      x[i] = v
    let
      st = jack x
      nv = naiveErr x
    note "jack AR(1) a=0.8: tau2=", st.tau2, " tau2p=", st.tau2p, " b=", st.blockSize,
      " nblock=", st.nblock, " err/naive=", st.err/nv, " neff=", st.neff
    check abs(st.tau2 - 9.0) < 1.2
    check abs(st.tau2p - 9.0) < 1.5
    check st.n == n
    check st.blockSize == max(1, int(ceil(st.tau2)))
    check st.nblock == (n + st.blockSize - 1) div st.blockSize
    check not st.provisional
    check abs(st.neff - float(n)/st.tau2) < 1e-9
    # the delete-b jackknife on AR(1) has E[var] = var_naive * S(b)/b with
    # S(b) = b + 2 sum_{k=1}^{b-1} (b-k) alpha^k -- b = 2 tau_int does NOT reach sqrt(2 tau)
    var sb = float(st.blockSize)
    for k in 1..<st.blockSize: sb += 2.0*float(st.blockSize - k)*pow(alpha, float k)
    let pred = sqrt(sb/float(st.blockSize))
    note "jack AR(1) err/naive=", st.err/nv, " blocked prediction=", pred,
      " sqrt(2 tau)=", sqrt 9.0
    check st.err/nv > 2.0
    check abs(st.err/nv/pred - 1.0) < 0.12

  test "uncorrelated series gives 2 tau ~ 1 and the naive error":
    const n = 16384
    var
      r = initRand 2718
      x = newSeq[float](n)
    for i in 0..<n: x[i] = r.gauss(0.0, 1.0)
    let
      st = jack x
      nv = naiveErr x
    note "jack iid: tau2=", st.tau2, " tau2p=", st.tau2p, " b=", st.blockSize,
      " err/naive=", st.err/nv
    check abs(st.tau2 - 1.0) < 0.1
    check abs(st.tau2p - 1.0) < 0.1
    check abs(st.err/nv - 1.0) < 0.05
    # block size 1 is the delete-1 jackknife, identically the naive error
    let s1 = jack(x, bs = 1)
    note "jack delete-1 err/naive - 1 = ", s1.err/nv - 1.0
    check s1.blockSize == 1
    check s1.nblock == n
    check abs(s1.err/nv - 1.0) < 1e-10
    check abs(s1.bias) < 1e-12

  test "stride and provisional are reported":
    var x = newSeq[float](6)
    for i in 0..<6: x[i] = float(i)
    let st = jack(x, stride = 5, bs = 3)
    check st.stride == 5
    check st.blockSize == 3
    check st.nblock == 2
    check st.provisional
    check abs(st.mean - 2.5) < 1e-13

# --- TSV ----------------------------------------------------------------------

suite "dataio":
  test "round trip is bit identical and keeps metadata":
    let
      path = getTempDir()/"radial_tfit_roundtrip.tsv"
      meta = @[("lattice", "L4"), ("nt", "168"), ("at", $(16.0/168.0)),
               ("note", "free limit, T=16"), ("eq", "a=b=c")]
      names = @["t", "g", "e"]
    var cols = newSeq[seq[float]](3)
    for i in 0..<32:
      cols[0].add float(i)
      cols[1].add 1.0/(16.0*PI*sinh(0.5*(0.1 + 0.1*float(i)))^2)
      cols[2].add (if i mod 3 == 0: -1.0/3.0 else: 1e-300*float(i+1))
    cols[1].add 1e300
    cols[0].add 0.0
    cols[2].add PI
    writeTsv(path, meta, names, cols)
    let got = readTsv path
    let header = readTsvMeta path
    note "tsv round trip rows=", cols[0].len, " cols=", cols.len, " meta=", meta.len
    check got.names == names
    check got.meta.len == meta.len
    check header.len == meta.len
    for (k, v) in meta:
      check got.meta[k] == v
      check header[k] == v
    check not header.hasKey("columns")
    check got.cols.len == 3
    for j in 0..2:
      check got.cols[j].len == cols[j].len
      for i in 0..<cols[j].len:
        check got.cols[j][i] == cols[j][i]         # bit identical, that is what %.17g is for
    removeFile path

  test "comments and blank lines are tolerated":
    let path = getTempDir()/"radial_tfit_comments.tsv"
    writeFile(path, """
# a free-form comment with no key
# key=value

# columns=x y
1	2

3	4
""")
    let got = readTsv path
    check got.names == @["x", "y"]
    check got.meta.len == 1
    check got.meta["key"] == "value"
    check got.cols == @[@[1.0, 3.0], @[2.0, 4.0]]
    removeFile path
