## Fitting and statistics for the radial (S^2 x R) analysis.
##
## Conventions are doc/02-formulation.md section 7.3, i.e. the paper's (V.4)-(V.9):
##   f(t)           = arccosh[ G(t) / G(T/2) ]                        (V.4)
##   Delta_eff(t)   = -[ f(t+a_t) - f(t) ] / a_t                      (V.5)
##   Delta_eff(t)  ~= Delta_0 + c exp(-Delta' t)                      (V.6)
##   Delta_0(as,at) = Delta_0^cont + c_s abar_s^2 + c_t a_t^2         (V.7)
##   g(t)          ~= C * G(t; nmax), minimized over the integer nmax (V.9)
##
## Error convention: where per-point errors are supplied the parameter covariance is
## (J^T W J)^{-1} with W = diag(1/e^2), unscaled.  `plateauFit` called without errors
## uses unit weights and scales the covariance by chi2/dof, the usual unweighted rule.

import std/math
import utils/resample

type
  PlateauFit* = object
    d0*, c*, dp*: float               ## Delta_0, c, Delta' of (V.6)
    ed0*, ec*, edp*: float            ## 1 sigma parameter errors
    cov*: array[3, array[3, float]]   ## full covariance, order (d0, c, dp)
    chi2*: float
    dof*: int                         ## npoint - 3
    iters*: int
    converged*: bool

  LineFit* = object
    a*, b*: float                     ## intercept, slope
    ea*, eb*: float
    cab*: float                       ## cov(a, b)
    chi2*: float
    dof*: int                         ## npoint - 2

  PlaneFit* = object
    a*, cs*, ct*: float               ## Delta_0^cont, c_s, c_t of (V.7)
    ea*, ecs*, ect*: float
    cov*: array[3, array[3, float]]   ## full covariance, order (a, cs, ct)
    chi2*: float
    dof*: int                         ## npoint - 3

  NmaxFit* = object
    nmax*: int                        ## best truncation
    c*: float                         ## fitted overall normalization
    res*: float                       ## weighted residual sum at `nmax`
    dof*: int                         ## npoint - 2, the two fitted quantities being C and nmax
    resDof*: float                    ## res/dof

  SeriesStat* = object
    mean*: float
    err*: float                       ## blocked jackknife error of the mean
    bias*: float                      ## Quenouille bias estimate
    tau2*: float                      ## 2 tau_int, Wolff automatic window
    tau2p*: float                     ## 2 tau_int, positive-sequence truncation
    neff*: float                      ## n / max(1, tau2)
    n*: int
    stride*: int                      ## measurement separation of the samples; recorded only
    blockSize*: int
    nblock*: int
    provisional*: bool                ## nblock < 4: `err` is an order of magnitude, nothing more

func chi2dof*(f: PlateauFit|LineFit|PlaneFit): float = f.chi2/float(f.dof)

proc inv3(a: var array[3, array[3, float]]): bool =
  ## In-place 3x3 inverse, Gauss-Jordan with partial pivoting.  False if singular.
  var b: array[3, array[3, float]]
  for i in 0..2: b[i][i] = 1.0
  for k in 0..2:
    var p = k
    for i in k+1..2:
      if abs(a[i][k]) > abs(a[p][k]): p = i
    if a[p][k] == 0.0: return false
    if p != k:
      swap a[p], a[k]
      swap b[p], b[k]
    let d = 1.0/a[k][k]
    for j in 0..2:
      a[k][j] *= d
      b[k][j] *= d
    for i in 0..2:
      if i != k:
        let f = a[i][k]
        if f != 0.0:
          for j in 0..2:
            a[i][j] -= f*a[k][j]
            b[i][j] -= f*b[k][j]
  a = b
  true

proc sym3(c: var array[3, array[3, float]]) =
  ## A covariance must be symmetric; neither the Gauss-Jordan inverse nor the
  ## M C M^T change of variables preserves that bit for bit.
  for i in 0..2:
    for j in i+1..2:
      let s = 0.5*(c[i][j] + c[j][i])
      c[i][j] = s
      c[j][i] = s

proc effMass*(c: openArray[float], at, T: float): seq[float] =
  ## (V.4)-(V.5).  `c` is sampled at t = i*at; the result is Delta_eff at the same
  ## t for i = 0 .. c.len-2.  This is NOT the log-ratio effective mass: the reference
  ## point is the single value c(T/2), at index round(T/(2 at)).
  ## Only t + at <= T/2 is meaningful; beyond that f(t) folds back.
  let
    ih = int(round(0.5*T/at))
    n = c.len
  var f = newSeq[float](n)
  for i in 0..<n: f[i] = arccosh(c[i]/c[ih])   # exact 1.0 at i = ih, so f[ih] = 0
  result = newSeq[float](n-1)
  for i in 0..<n-1: result[i] = -(f[i+1] - f[i])/at

proc plateauFit*(m: openArray[float], t0, t1: int, at = 1.0,
                 e: openArray[float] = [], maxit = 200): PlateauFit =
  ## Fit Delta_eff(t) = Delta_0 + c exp(-Delta' t), (V.6), over the index window
  ## [t0, t1) with t = at*i.  The paper's fermion window is 4 <= t < 8 = T/2.
  ## `e` is indexed like `m` (e[i] pairs with m[i]); empty means unit weights.
  ##
  ## The exponential is fitted about the window origin tref = at*t0 and c is mapped
  ## back afterwards.  Without that shift the 4 <= t < 8 window is badly conditioned:
  ## exp(-Delta' t) is ~1e-2 across the whole window and c is ~1e2 times the signal.
  ## Started from a variable-projection scan over Delta' (exact linear solve for
  ## (Delta_0, c) at each Delta'), geometric over [0.01/span, 10/at], then
  ## Levenberg-Marquardt on all three.  Only the scan assumes Delta' > 0, which (V.6)
  ## does too -- it is an excited-state gap; LM afterwards is unconstrained.
  let
    n = t1 - t0
    tref = at*float(t0)
  var
    u = newSeq[float](n)
    y = newSeq[float](n)
    w = newSeq[float](n)
  for i in 0..<n:
    u[i] = at*float(i)
    y[i] = m[t0+i]
    w[i] = if e.len == 0: 1.0 else: 1.0/(e[t0+i]*e[t0+i])

  proc chi2At(d0, cp, dp: float): float =
    for i in 0..<n:
      let r = y[i] - (d0 + cp*exp(-dp*u[i]))
      result += w[i]*r*r

  proc linAt(dp: float): tuple[d0, cp: float] =
    ## (Delta_0, c) by weighted least squares at fixed Delta'.
    var s, sx, sxx, sy, sxy = 0.0
    for i in 0..<n:
      let x = exp(-dp*u[i])
      s += w[i]
      sx += w[i]*x
      sxx += w[i]*x*x
      sy += w[i]*y[i]
      sxy += w[i]*x*y[i]
    let d = s*sxx - sx*sx
    if d == 0.0: (sy/s, 0.0)
    else: ((sxx*sy - sx*sxy)/d, (s*sxy - sx*sy)/d)

  proc normal(p: array[3, float], jtj: var array[3, array[3, float]],
              jtr: var array[3, float]) =
    for a in 0..2:
      jtr[a] = 0.0
      for b in 0..2: jtj[a][b] = 0.0
    for i in 0..<n:
      let
        ex = exp(-p[2]*u[i])
        g = [1.0, ex, -p[1]*u[i]*ex]     # dF/d(d0, cp, dp)
        r = y[i] - (p[0] + p[1]*ex)
      for a in 0..2:
        jtr[a] += w[i]*r*g[a]
        for b in a..2: jtj[a][b] += w[i]*g[a]*g[b]
    for a in 1..2:
      for b in 0..<a: jtj[a][b] = jtj[b][a]

  # Start from the constant model, a finite upper bound on chi2 that any candidate must beat.
  let span = max(u[n-1], at)
  var sw, swy = 0.0
  for i in 0..<n:
    sw += w[i]
    swy += w[i]*y[i]
  var p = [swy/sw, 0.0, 1.0/span]
  var chi2 = chi2At(p[0], p[1], p[2])

  const ngrid = 400
  let
    dpLo = 0.01/span
    dpHi = 10.0/at
    fac = pow(dpHi/dpLo, 1.0/float(ngrid-1))
  var dp = dpLo
  for k in 0..<ngrid:
    let (d0, cp) = linAt(dp)
    let c2 = chi2At(d0, cp, dp)
    if c2 < chi2:
      chi2 = c2
      p = [d0, cp, dp]
    dp *= fac

  var
    jtj: array[3, array[3, float]]
    jtr: array[3, float]
    lam = 1e-3
    it = 0
  while it < maxit:
    inc it
    normal(p, jtj, jtr)
    var aa = jtj
    for a in 0..2: aa[a][a] *= 1.0 + lam
    if inv3 aa:
      var d: array[3, float]
      for a in 0..2:
        var s = 0.0
        for b in 0..2: s += aa[a][b]*jtr[b]
        d[a] = s
      let c2 = chi2At(p[0]+d[0], p[1]+d[1], p[2]+d[2])
      if c2 < chi2:
        var rel = 0.0
        for a in 0..2: rel = max(rel, abs(d[a])/(abs(p[a]) + 1e-12))
        for a in 0..2: p[a] += d[a]
        chi2 = c2
        lam = max(0.1*lam, 1e-14)
        if rel < 1e-13:
          result.converged = true
          break
        continue
    lam *= 10.0
    if lam > 1e14:                     # step rejected at the roundoff floor
      result.converged = true
      break

  normal(p, jtj, jtr)
  result.chi2 = chi2
  result.dof = n - 3
  result.iters = it
  var cv = jtj
  if inv3 cv:
    if e.len == 0 and result.dof > 0:
      let f = chi2/float(result.dof)
      for a in 0..2:
        for b in 0..2: cv[a][b] *= f
    # cp = c exp(-dp tref), so (d0, cp, dp) -> (d0, c, dp) with c = cp exp(+dp tref)
    let
      ee = exp(p[2]*tref)
      cc = p[1]*ee
    var mm, tm: array[3, array[3, float]]
    mm[0][0] = 1.0
    mm[1][1] = ee
    mm[1][2] = tref*cc
    mm[2][2] = 1.0
    for a in 0..2:
      for b in 0..2:
        var s = 0.0
        for k in 0..2: s += mm[a][k]*cv[k][b]
        tm[a][b] = s
    for a in 0..2:
      for b in 0..2:
        var s = 0.0
        for k in 0..2: s += tm[a][k]*mm[b][k]
        result.cov[a][b] = s
    sym3 result.cov
  result.d0 = p[0]
  result.c = p[1]*exp(p[2]*tref)
  result.dp = p[2]
  result.ed0 = sqrt result.cov[0][0]
  result.ec = sqrt result.cov[1][1]
  result.edp = sqrt result.cov[2][2]

proc contFit*(x, y, e: openArray[float]): LineFit =
  ## Inverse-variance weighted straight line y = a + b x.  For the O(a^2) continuum
  ## extrapolation pass x = abar_s^2 (or a_t^2); the intercept is the continuum limit.
  ## Centered form (Numerical Recipes `fit`): no cancellation when the x are far from 0.
  var s, sx, sy = 0.0
  for i in 0..<x.len:
    let w = 1.0/(e[i]*e[i])
    s += w
    sx += w*x[i]
    sy += w*y[i]
  let xb = sx/s
  var st2, sty = 0.0
  for i in 0..<x.len:
    let
      w = 1.0/(e[i]*e[i])
      t = x[i] - xb
    st2 += w*t*t
    sty += w*t*y[i]
  result.b = sty/st2
  result.a = sy/s - xb*result.b
  result.ea = sqrt(1.0/s + xb*xb/st2)
  result.eb = sqrt(1.0/st2)
  result.cab = -xb/st2
  for i in 0..<x.len:
    let r = y[i] - result.a - result.b*x[i]
    result.chi2 += r*r/(e[i]*e[i])
  result.dof = x.len - 2

proc contFit2*(as2, at2, y, e: openArray[float]): PlaneFit =
  ## (V.7): weighted three-parameter linear fit y = a + c_s*as2 + c_t*at2 over a
  ## two-dimensional grid of lattices.  `as2` and `at2` are ALREADY SQUARED, i.e.
  ## abar_s^2 and a_t^2; `a` is Delta_0^cont.
  ## The two columns are centered on their weighted means before the normal equations
  ## are formed and the result is mapped back, which keeps the 3x3 system well
  ## conditioned when abar_s^2 and a_t^2 are both far from zero.
  let n = y.len
  var w = newSeq[float](n)
  var s, s1, s2 = 0.0
  for i in 0..<n:
    w[i] = 1.0/(e[i]*e[i])
    s += w[i]
    s1 += w[i]*as2[i]
    s2 += w[i]*at2[i]
  let
    m1 = s1/s
    m2 = s2/s
  var
    aa: array[3, array[3, float]]
    rhs: array[3, float]
  for i in 0..<n:
    let g = [1.0, as2[i] - m1, at2[i] - m2]
    for a in 0..2:
      rhs[a] += w[i]*y[i]*g[a]
      for b in 0..2: aa[a][b] += w[i]*g[a]*g[b]
  var cv = aa
  if not inv3 cv:                       # would otherwise return finite garbage
    raise newException(ValueError,
      "contFit2: the (abar_s^2, a_t^2) points do not span a plane")
  var bt: array[3, float]
  for a in 0..2:
    for b in 0..2: bt[a] += cv[a][b]*rhs[b]
  for i in 0..<n:
    let r = y[i] - (bt[0] + bt[1]*(as2[i] - m1) + bt[2]*(at2[i] - m2))
    result.chi2 += w[i]*r*r
  result.dof = n - 3
  result.a = bt[0] - m1*bt[1] - m2*bt[2]
  result.cs = bt[1]
  result.ct = bt[2]
  # (a', cs, ct) -> (a, cs, ct)
  var mm, tm: array[3, array[3, float]]
  mm[0][0] = 1.0
  mm[0][1] = -m1
  mm[0][2] = -m2
  mm[1][1] = 1.0
  mm[2][2] = 1.0
  for a in 0..2:
    for b in 0..2:
      var q = 0.0
      for k in 0..2: q += mm[a][k]*cv[k][b]
      tm[a][b] = q
  for a in 0..2:
    for b in 0..2:
      var q = 0.0
      for k in 0..2: q += tm[a][k]*mm[b][k]
      result.cov[a][b] = q
  sym3 result.cov
  result.ea = sqrt result.cov[0][0]
  result.ecs = sqrt result.cov[1][1]
  result.ect = sqrt result.cov[2][2]

proc nmaxFit*(g: openArray[float], model: proc(nmax: int): seq[float],
              nmaxRange: Slice[int], e: openArray[float] = []): NmaxFit =
  ## (V.9): minimize over the integer truncation `nmax` and a floating overall
  ## normalization C,
  ##   res(nmax) = sum_t w_t (g_t - C m_t)^2,   C = sum w g m / sum w m^2,
  ## with m = model(nmax).  DOF = g.len - 2 (C and nmax).
  ## An empty `e` selects the RELATIVE residual w_t = 1/g_t^2.  That is the scale-free
  ## choice for the deterministic free-limit correlators, which carry no statistical
  ## error and span several orders of magnitude across the fit range.
  let n = g.len
  var w = newSeq[float](n)
  for i in 0..<n:
    w[i] = if e.len == 0: 1.0/(g[i]*g[i]) else: 1.0/(e[i]*e[i])
  result.dof = n - 2
  for nm in nmaxRange.a..nmaxRange.b:
    let m = model(nm)
    var a, b = 0.0
    for i in 0..<n:
      a += w[i]*g[i]*m[i]
      b += w[i]*m[i]*m[i]
    let c = a/b
    var r = 0.0
    for i in 0..<n:
      let d = g[i] - c*m[i]
      r += w[i]*d*d
    if nm == nmaxRange.a or r < result.res:
      result.nmax = nm
      result.c = c
      result.res = r
  result.resDof = result.res/float(result.dof)

proc seriesMean(xs: Ensemble[seq[float]]): float =
  var m = 0.0
  for i in 0..<xs.len: m += (xs[i] - m)/float(i+1)
  m

proc jack*(s: openArray[float], stride = 1, bs = 0): SeriesStat =
  ## Blocked jackknife mean and error with the block size taken from the measured
  ## integrated autocorrelation time: blockSize = ceil(max(1, 2 tau_int)).  `bs > 0`
  ## overrides it.  Fewer than 4 blocks sets `provisional`.
  ##
  ## `tau2`/`tau2p` are in units of samples; multiply by `stride` for trajectory units.
  ## `stride` itself is only recorded so a result can be interpreted on its own.
  ##
  ## Note that a block size of exactly 2 tau_int recovers only part of the asymptotic
  ## variance -- for AR(1) the blocked error is sqrt(S(b)/b) times the naive one with
  ## S(b) = b + 2 sum_{k<b} (b-k) rho(k), which for rho = 0.8^k is 2.27 rather than
  ## sqrt(2 tau_int) = 3.  Use `bs` explicitly if a larger block is wanted.
  var xs = newSeq[float](s.len)
  for i in 0..<s.len: xs[i] = s[i]
  let
    n = xs.len
    tau2 = xs.jackknife(n, intAutocorr).mean       # one block: a single full-series pass
    tau2p = xs.jackknife(n, intAutocorrPositive).mean
    b = if bs > 0: bs else: max(1, int(ceil(max(1.0, tau2))))
    st = xs.jackknife(b, seriesMean)
    nb = (n + b - 1) div b
  SeriesStat(mean: st.mean, err: st.stdev, bias: st.bias,
             tau2: tau2, tau2p: tau2p, neff: float(n)/max(1.0, tau2),
             n: n, stride: stride, blockSize: b, nblock: nb, provisional: nb < 4)
