## Lattice-independent analysis utilities for the 16-cell honeycomb project.
##
## Nothing here knows about QEX fields or about the lattice geometry, so the
## same code serves the cubic reference pipeline (`refCubicGen`/`refCubicMeas`)
## and the honeycomb measurements.  Pure Nim, only `math`/`strutils`.
##
## Contents
##   * flow-scale finders   `findT0`, `findW0`, `findCrossing`, `interpAt`
##   * statistics           `mean`, `variance`, `stddev`, `stderrMean`,
##                          `autocorrTime`, `autocorr`
##   * resampling           `jackknife` (delete-1, optionally binned)
##   * fitting              `fitPoly` (weighted linear least squares in
##                          arbitrary integer powers of x)
##
## Test: `src/experimental/honeycomb/tests/tanalysis.nim`.

import std/[math, strutils]

export math.sqrt

const
  noCrossing* = -1.0
    ## Sentinel returned by `findCrossing`/`findT0`/`findW0` when the target is
    ## never reached.  A sentinel rather than `NaN` because QEX is built with
    ## `-Ofast -ffast-math`, under which `x != x` is optimised away.

# ---------------------------------------------------------------- basic stats

proc mean*(x: openArray[float]): float =
  ## Arithmetic mean.  0 for an empty input.
  if x.len == 0: return 0.0
  var s = 0.0
  for v in x: s += v
  s / x.len.float

proc variance*(x: openArray[float]): float =
  ## Unbiased (n-1) sample variance.
  if x.len < 2: return 0.0
  let m = x.mean
  var s = 0.0
  for v in x: s += (v-m)*(v-m)
  s / (x.len-1).float

proc stddev*(x: openArray[float]): float =
  ## Unbiased sample standard deviation.
  sqrt x.variance

proc stderrMean*(x: openArray[float]): float =
  ## Naive standard error of the mean, sqrt(var/n).  Ignores autocorrelation;
  ## multiply by sqrt(2 tau_int) to correct (see `autocorrTime`).
  if x.len < 2: return 0.0
  sqrt(x.variance / x.len.float)

proc cov(x: openArray[float], m: float, t: int): float =
  var s = 0.0
  for i in 0..<(x.len-t): s += (x[i]-m)*(x[i+t]-m)
  s / float(x.len-t)

proc autocorr*(x: openArray[float], tmax = -1): seq[float] =
  ## Normalised autocorrelation function `rho(t) = Gamma(t)/Gamma(0)`,
  ## `t = 0 .. tmax` (default `min(n-1, n div 2)`).
  let n = x.len
  if n < 2: return @[1.0]
  let tm = if tmax >= 0: min(tmax, n-1) else: min(n-1, n div 2)
  let m = x.mean
  var g = newSeq[float](tm+1)
  for t in 0..tm:
    g[t] = cov(x, m, t)
  let g0 = g[0]
  if g0 == 0.0:
    for t in 0..tm: g[t] = if t == 0: 1.0 else: 0.0
    return g
  for t in 0..tm: g[t] = g[t]/g0
  g

proc autocorrTimeW*(x: openArray[float], c = 5.0):
    tuple[tau, dtau: float, window: int] =
  ## Integrated autocorrelation time with the Madras-Sokal automatic window:
  ## `tau_int(W) = 1/2 + sum_{t=1..W} rho(t)`, W the smallest value with
  ## `W >= c*tau_int(W)`.  `dtau` is the Madras-Sokal error estimate
  ## `tau*sqrt(2(2W+1)/n)`.
  let n = x.len
  if n < 4: return (0.5, 0.0, 0)
  let m = x.mean
  let g0 = cov(x, m, 0)
  var tau = 0.5
  var w = n div 2
  for t in 1..w:
    if g0 != 0.0: tau += cov(x, m, t)/g0
    if tau < 0.5: tau = 0.5      # guard against noise driving tau negative
    if float(t) >= c*tau:
      w = t
      break
  let dtau = tau*sqrt(2.0*float(2*w+1)/float(n))
  (tau, dtau, w)

proc autocorrTime*(x: openArray[float], c = 5.0): float =
  ## Integrated autocorrelation time (see `autocorrTimeW`).  A completely
  ## uncorrelated series gives 0.5; the number of effectively independent
  ## samples is `n/(2 tau)`.
  x.autocorrTimeW(c).tau

proc binned*(x: openArray[float], b: int): seq[float] =
  ## Block averages of `x` in bins of `b` consecutive entries.  A trailing
  ## partial bin is dropped.
  if b <= 1: return @x
  let nb = x.len div b
  result = newSeq[float](nb)
  for i in 0..<nb:
    var s = 0.0
    for j in 0..<b: s += x[i*b+j]
    result[i] = s/b.float

# ------------------------------------------------------------------ jackknife

proc jackknife*[T](x: openArray[T], f: proc(s: openArray[T]): float,
                   bin = 1): tuple[mean, err: float] =
  ## Delete-`bin` jackknife of the estimator `f`.
  ##
  ## `mean` is the full-sample estimate `f(x)`.  Equal-sized blocks use the
  ## standard `sqrt((nb-1)/nb * sum_i (theta_i - theta_bar)^2)` error.  A
  ## trailing partial block is handled with the unequal-delete (`delete-m_j`)
  ## pseudovalue variance, which reduces exactly to that formula when all
  ## blocks have equal size.  For `f = mean` and `bin = 1` this reproduces the
  ## ordinary standard error of the mean exactly.
  ##
  ## `bin > 1` deletes whole blocks of `bin` consecutive entries, which is the
  ## standard way to absorb autocorrelation.  A trailing partial block is kept
  ## as its own, correctly weighted shorter block.
  let n = x.len
  if n == 0: return (0.0, 0.0)
  result.mean = f(x)
  if n == 1: return
  let b = max(1, bin)
  let nb = (n + b - 1) div b
  if nb < 2: return
  var pseudo = newSeq[float](nb)
  var hs = newSeq[float](nb)
  var sample = newSeq[T](0)
  for i in 0..<nb:
    let lo = i*b
    let hi = min(n, lo+b)
    sample.setLen(0)
    for j in 0..<n:
      if j < lo or j >= hi: sample.add x[j]
    let h = n.float/(hi-lo).float
    let deleted = f(sample)
    hs[i] = h
    pseudo[i] = h*result.mean - (h-1.0)*deleted
  var jackmean = 0.0
  for i in 0..<nb: jackmean += pseudo[i]/hs[i]
  var v = 0.0
  for i in 0..<nb:
    let d = pseudo[i] - jackmean
    v += d*d/(hs[i]-1.0)
  result.err = sqrt(max(0.0, v/nb.float))

proc jackknifeMean*(x: openArray[float], bin = 1): tuple[mean, err: float] =
  ## Convenience wrapper: jackknife of the plain mean.
  jackknife(x, proc(s: openArray[float]): float = s.mean, bin)

# -------------------------------------------------------------- flow  scales

proc interpAt*(x, y: openArray[float], x0: float): float =
  ## Linear interpolation of y(x), clamped at the ends; 0 for empty inputs.
  ## x must be increasing, and x and y must have equal lengths.
  if x.len != y.len:
    raise newException(ValueError, "interpAt requires x and y of equal length")
  let n = x.len
  if n == 0: return 0.0
  if x0 <= x[0]: return y[0]
  for i in 1..<n:
    if x[i] >= x0:
      let f = (x0-x[i-1])/(x[i]-x[i-1])
      return y[i-1] + f*(y[i]-y[i-1])
  y[n-1]

proc lagrange(xs, ys: openArray[float], x: float): float =
  ## Lagrange interpolation through all supplied points.
  var s = 0.0
  for i in 0..<xs.len:
    var p = ys[i]
    for j in 0..<xs.len:
      if j != i: p *= (x-xs[j])/(xs[i]-xs[j])
    s += p
  s

proc findCrossing*(x, y: openArray[float], target: float, order = 1): float =
  ## Flow-time (or generic abscissa) at which the series `y(x)` first crosses
  ## `target` from below.  `x` must be increasing; x and y must have equal lengths.
  ##
  ## `order = 1` linear interpolation between the bracketing pair (the usual
  ## convention for t0).  `order = 3` fits a cubic through the two points on
  ## either side and solves it by bisection inside the bracket; this reduces the
  ## discretisation error of the *interpolation* from O(dt^2) to O(dt^4).
  ##
  ## Returns `noCrossing` (= -1) if the target is never reached.
  if x.len != y.len:
    raise newException(ValueError, "findCrossing requires x and y of equal length")
  let n = x.len
  if n < 2: return noCrossing
  var k = -1
  for i in 1..<n:
    if y[i-1] < target and y[i] >= target:
      k = i
      break
  if k < 0: return noCrossing
  if order <= 1 or n < 4:
    let f = (target - y[k-1])/(y[k] - y[k-1])
    return x[k-1] + f*(x[k]-x[k-1])
  # cubic through k-2,k-1,k,k+1 (shifted to stay inside the data)
  var i0 = k-2
  if i0 < 0: i0 = 0
  if i0+4 > n: i0 = n-4
  var xs, ys: array[4, float]
  for j in 0..3:
    xs[j] = x[i0+j]
    ys[j] = y[i0+j]
  var a = x[k-1]
  var b = x[k]
  # p(a)-target < 0 <= p(b)-target by construction of the bracket
  for _ in 0..<60:
    let m = 0.5*(a+b)
    if lagrange(xs, ys, m) < target: a = m else: b = m
  0.5*(a+b)

proc findT0*(t: openArray[float], t2E: openArray[float], target = 0.3,
             order = 1): float =
  ## Lüscher's `t0`: the flow time where `t^2 <E>(t) = target` (0.3 for SU(3)).
  ## Linear interpolation of the flow-time series by default.
  ## Returns `noCrossing` (= -1) if the series never reaches `target`.
  findCrossing(t, t2E, target, order)

proc findW0*(t: openArray[float], tdt2E: openArray[float], target = 0.3,
             order = 1): float =
  ## The BMW `w0` scale, defined by `W(t) = t d/dt [t^2 <E>] = target`.
  ##
  ## **Returns `w0^2`, i.e. the flow time at the crossing** (so that it is the
  ## direct analogue of `findT0`); take `sqrt` of the result to get `w0`.
  ## Returns `noCrossing` (= -1) if the target is never reached.
  findCrossing(t, tdt2E, target, order)

proc derivT2E*(t, t2E: openArray[float]): seq[float] =
  ## `W(t) = t d/dt [t^2 <E>]` by centred differences on a (not necessarily
  ## uniform) grid; one-sided at the ends.
  ## t must be increasing, and t and t2E must have equal lengths.
  if t.len != t2E.len:
    raise newException(ValueError, "derivT2E requires t and t2E of equal length")
  let n = t.len
  result = newSeq[float](n)
  if n < 2: return
  for i in 0..<n:
    var d: float
    if i == 0:
      d = (t2E[1]-t2E[0])/(t[1]-t[0])
    elif i == n-1:
      d = (t2E[n-1]-t2E[n-2])/(t[n-1]-t[n-2])
    else:
      let h1 = t[i]-t[i-1]
      let h2 = t[i+1]-t[i]
      d = (h1*h1*t2E[i+1] + (h2*h2-h1*h1)*t2E[i] - h2*h2*t2E[i-1]) /
          (h1*h2*(h1+h2))
    result[i] = t[i]*d

# -------------------------------------------------------------------- fitting

proc invertSym(a: seq[seq[float]]): tuple[inv: seq[seq[float]], ok: bool] =
  ## Gauss-Jordan inverse with partial pivoting.  `a` is small (n <= ~8).
  let n = a.len
  var m = newSeq[seq[float]](n)
  var inv = newSeq[seq[float]](n)
  for i in 0..<n:
    m[i] = a[i]
    inv[i] = newSeq[float](n)
    inv[i][i] = 1.0
  for c in 0..<n:
    var p = c
    for r in c+1..<n:
      if abs(m[r][c]) > abs(m[p][c]): p = r
    if abs(m[p][c]) < 1e-300: return (inv, false)
    if p != c:
      swap(m[p], m[c])
      swap(inv[p], inv[c])
    let d = 1.0/m[c][c]
    for j in 0..<n:
      m[c][j] *= d
      inv[c][j] *= d
    for r in 0..<n:
      if r == c: continue
      let f = m[r][c]
      if f == 0.0: continue
      for j in 0..<n:
        m[r][j] -= f*m[c][j]
        inv[r][j] -= f*inv[c][j]
  (inv, true)

proc fitPolyCov*(x, y, dy: openArray[float], powers: openArray[int]):
    tuple[coef, err: seq[float], cov: seq[seq[float]], chisq: float,
          dof: int, chisqDof: float] =
  ## Weighted linear least squares of `y +- dy` to `sum_k c_k x^{p_k}`.
  ## Returns the coefficients, their errors, the full covariance matrix, and
  ## chi^2 / dof.  `dy` entries must be > 0; x, y and dy must have equal lengths.
  if x.len != y.len or x.len != dy.len:
    raise newException(ValueError, "fitPolyCov requires x, y and dy of equal length")
  let n = x.len
  let np = powers.len
  result.coef = newSeq[float](np)
  result.err = newSeq[float](np)
  result.cov = newSeq[seq[float]](np)
  for i in 0..<np: result.cov[i] = newSeq[float](np)
  if n == 0 or np == 0: return
  var a = newSeq[seq[float]](n)
  for i in 0..<n:
    a[i] = newSeq[float](np)
    for k in 0..<np:
      a[i][k] = pow(x[i], powers[k].float)
  var mm = newSeq[seq[float]](np)
  for k in 0..<np: mm[k] = newSeq[float](np)
  var b = newSeq[float](np)
  for i in 0..<n:
    let w = 1.0/(dy[i]*dy[i])
    for k in 0..<np:
      b[k] += w*a[i][k]*y[i]
      for l in 0..<np:
        mm[k][l] += w*a[i][k]*a[i][l]
  let (inv, ok) = invertSym(mm)
  if not ok:
    # singular normal-equation matrix: leave coefficients zero and flag with a
    # negative chi^2 (NaN is unusable under -ffast-math).
    result.chisq = -1.0
    result.chisqDof = -1.0
    return
  result.cov = inv
  for k in 0..<np:
    var s = 0.0
    for l in 0..<np: s += inv[k][l]*b[l]
    result.coef[k] = s
    result.err[k] = sqrt(max(0.0, inv[k][k]))
  var chi = 0.0
  for i in 0..<n:
    var f = 0.0
    for k in 0..<np: f += result.coef[k]*a[i][k]
    let r = (y[i]-f)/dy[i]
    chi += r*r
  result.chisq = chi
  result.dof = n - np
  result.chisqDof = if result.dof > 0: chi/result.dof.float else: -1.0

proc fitPoly*(x, y, dy: openArray[float], powers: openArray[int]):
    tuple[coef, err: seq[float], chisqDof: float] =
  ## Weighted linear least squares in the given (integer) powers of `x`.
  ## `powers = [0, 1]` is the O(a^2) continuum extrapolation of slide 10 when
  ## `x = a^2/t0`; `powers = [0, 1, 2]` is the O(a^4) form.
  let r = fitPolyCov(x, y, dy, powers)
  (r.coef, r.err, r.chisqDof)

proc evalPoly*(coef: openArray[float], powers: openArray[int], x: float): float =
  ## Evaluate `sum_k coef[k] x^{powers[k]}`; coef and powers must have equal lengths.
  if coef.len != powers.len:
    raise newException(ValueError, "evalPoly requires coef and powers of equal length")
  var s = 0.0
  for k in 0..<coef.len:
    s += coef[k]*pow(x, powers[k].float)
  s

# ------------------------------------------------------------------- file I/O

proc readColumns*(fn: string): seq[seq[float]] =
  ## Read a whitespace-separated numeric table, skipping blank lines and lines
  ## whose first non-blank character is '#'.  Returns rows.
  result = @[]
  for line in lines(fn):
    let s = line.strip
    if s.len == 0 or s[0] == '#': continue
    var row: seq[float] = @[]
    var bad = false
    for tok in s.splitWhitespace:
      try: row.add parseFloat(tok)
      except ValueError: bad = true
    if not bad and row.len > 0: result.add row

when isMainModule:
  # Tiny self-demo so that `make run .../hcanalysis.nim` does something useful.
  var t, e: seq[float]
  var i = 0
  while i <= 200:
    let tt = 0.01*i.float
    t.add tt
    e.add 0.3*(tt/1.234)          # t^2E linear in t, crossing 0.3 at t=1.234
    inc i
  echo "findT0 (linear) = ", findT0(t, e)
  echo "mean/stddev of {1..10}: ",
       mean([1.0,2,3,4,5,6,7,8,9,10]), " ", stddev([1.0,2,3,4,5,6,7,8,9,10])
