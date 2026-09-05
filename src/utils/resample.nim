## Resampling and statistics with autocorrelations
##
## THE JACKKNIFE, THE BOOTSTRAP, AND OTHER RESAMPLING PLANS
## - BRADLEY EFRON, 1980

import math

type
  JackknifeStat*[Value] = object
    mean*: Value  ## the expectation value
    jksamples*: seq[Value]  ## expectation values of the jackknife samples
    bias*: Value  ## Quenouille's bias estimate
      # for "grouped jackknife" with n samples, and block size k
      # E = m + m1/n + O(n^{-2})
      # Ej = m + m1/(n-k) + O((n-k)^{-2})
      # m = (n E - (n-k) Ej) / k
      # bias = E - m = (n-k)/k (Ej-E)
      # In the unequal last-block case (n % k != 0),
      # let k_i be the delete size for group i and define
      #   h_i = n/k_i
      #   zeta_i = (n/k_i) * E - ((n-k_i)/k_i) * E_{(i)}
      # Then E_jack = sum_i zeta_i/h_i, bias = E - E_jack, and
      #   Var(E) ≈ (1/g) * sum_i (zeta_i-E_jack)^2/(h_i-1).
      # For equal blocks h_i=g, reducing to the usual grouped-jackknife
      # formulas.
    stdev*: Value  ## the jackknife estimate of the standard deviation of the expectation value
    hasStdev*: bool  ## whether stdev is defined

  EnsembleKind = enum
    EKoriginal, EKjackknife
  Ensemble*[D] = object
    case kind: EnsembleKind
    of EKoriginal:
      discard
    of EKjackknife:
      skiplow:int
      blocksize:int
    data:ptr D

func wrapEnsemble[D](ensemble:ptr D):Ensemble[D] =
  Ensemble[D](kind:EKoriginal, data:ensemble)

func jackknifeSample[D](ensemble:ptr D, blocksize:int, index:int):Ensemble[D] =
  Ensemble[D](kind:EKjackknife, skiplow:blocksize*index, blocksize:blocksize, data:ensemble)

func `[]`*[D](sample:Ensemble[D], i:int):auto =
  case sample.kind:
  of EKoriginal:
    sample.data[][i]
  of EKjackknife:
    if i<sample.skiplow:
      sample.data[][i]
    else:
      sample.data[][i+sample.blocksize]

func len*[D](sample:Ensemble[D]):int =
  let n = len(sample.data[])
  case sample.kind:
  of EKoriginal:
    n
  of EKjackknife:
    if sample.skiplow+sample.blocksize>n:
      sample.skiplow
    else:
      n - sample.blocksize

proc intAc(xs: Ensemble[seq[float]], mean: float, wolff: bool): float =
  ## rho(t) = sum_i dx_i*dx_(i+t)/(n*C(0)); return 2*tau_int.
  ## Ulli Wolff's automatic window (CPC 156 (2004) 143), S = 1.5.
  let n = xs.len
  if n < 2: return 1.0
  var c0 = 0.0
  for i in 0..<n:
    let d = xs[i]-mean
    c0 += d*d
  c0 /= float(n)
  if c0 == 0.0: return 1.0
  var tau = 0.5
  let wmax = if wolff: n div 2 else: n-1
  for lag in 1..wmax:
    var c = 0.0
    for i in 0..<(n-lag):
      c += (xs[i]-mean)*(xs[i+lag]-mean)
    let rho = c/(float(n)*c0)
    if not wolff and rho <= 0.0: break
    tau += rho
    if wolff:
      if tau <= 0.5: break
      let
        w = float(lag)
        tw = 1.5/ln((2.0*tau+1.0)/(2.0*tau-1.0))
      if exp(-w/tw) < tw/sqrt(w*float(n)): break
    elif float(lag) >= 6.0*tau:
      break
  2.0*tau

proc intAc(xs: Ensemble[seq[float]], wolff: bool): float =
  var mean = 0.0
  for i in 0..<xs.len:
    mean += (xs[i]-mean)/float(i+1)
  intAc(xs, mean, wolff)

proc intAutocorr*(xs: Ensemble[seq[float]]): float =
  intAc(xs, true)

proc intAutocorr*(xs: Ensemble[seq[float]], mean: float): float =
  intAc(xs, mean, true)

proc intAutocorrPositive*(xs: Ensemble[seq[float]]): float =
  intAc(xs, false)

proc jackknife*[V](mean: V, reps: seq[V], n, blocksize: int): JackknifeStat[V] =
  ## Reduce saved delete-block estimates; reps[i] omits group i from n samples.
  if blocksize <= 0 or n < 0:
    raise newException(ValueError, "jackknife: require positive blocksize and nonnegative sample count")
  let g = reps.len
  if g != (n + blocksize - 1) div blocksize:
    raise newException(ValueError, "jackknife: replica count does not match the groups")
  var
    zeta = newSeq[V](g)
    hs = newSeq[float](g)
  for i in 0..<g:
    let ki = min(blocksize, n - i*blocksize)
    hs[i] = float(n)/float(ki)
    zeta[i] = hs[i]*mean - (hs[i]-1.0)*reps[i]
  var jackmean = V(0)
  for i in 0..<g:
    jackmean += zeta[i]/hs[i]
  var variance = V(0)
  if g > 1:
    for i in 0..<g:
      let dz = zeta[i] - jackmean
      variance += dz*dz/(hs[i]-1.0)
    variance /= float(g)
  JackknifeStat[V](mean:mean, jksamples:reps,
    bias:mean - jackmean,
    stdev:if g > 1: sqrt(variance) else: V(0),
    hasStdev:g > 1)

proc jackknife*[D,V,A](ensemble:D, blocksize:int, estimator:proc(x:Ensemble[D], arg:A):V, arg:A): auto =
  ## Perform grouped jackknife with blocksize
  ## The expectation value: estimator(ensemble, arg)
  let
    m = estimator(wrapEnsemble(addr ensemble), arg)
    n = ensemble.len
    g = (n+blocksize-1) div blocksize
  var jk = newseq[V](g)
  for i in 0..<g:
    jk[i] = estimator(jackknifeSample(addr ensemble, blocksize, i), arg)
  jackknife(m, jk, n, blocksize)

proc jackknife*[D,V](ensemble:D, blocksize:int, estimator:proc(x:Ensemble[D]):V): auto =
  ## Perform grouped jackknife with blocksize
  ## The expectation value: estimator(ensemble)
  let
    m = estimator(wrapEnsemble(addr ensemble))
    n = ensemble.len
    g = (n+blocksize-1) div blocksize
  var jk = newseq[V](g)
  for i in 0..<g:
    jk[i] = estimator(jackknifeSample(addr ensemble, blocksize, i))
  jackknife(m, jk, n, blocksize)

type WeightedValue = tuple[x, w: float]

proc weightedMean(xs: Ensemble[seq[WeightedValue]]): float =
  var s, w = 0.0
  for i in 0..<xs.len:
    let x = xs[i]
    s += x.w*x.x
    w += x.w
  if w == 0.0: NaN else: s/w

proc weightedRms(xs: Ensemble[seq[WeightedValue]]): float =
  var s, w = 0.0
  for i in 0..<xs.len:
    let x = xs[i]
    s += x.w*x.x*x.x
    w += x.w
  if w == 0.0: NaN else: sqrt(s/w)

proc weightedJackknife*(values, weights: seq[float]; blocksize: int; isRms = false): JackknifeStat[float] =
  ## Weighted block jackknife; zero weights keep their block positions.
  if values.len != weights.len:
    raise newException(ValueError, "weightedJackknife values/weights length mismatch")
  var xs = newSeq[WeightedValue](values.len)
  for i in 0..<xs.len:
    xs[i] = (values[i], weights[i])
  result = xs.jackknife(blocksize, if isRms: weightedRms else: weightedMean)
  var
    lo = 0
    blocks = 0
  while lo < weights.len:
    var w = 0.0
    for i in lo..<min(lo + blocksize, weights.len):
      w += weights[i]
    if w != 0.0: inc blocks
    lo += blocksize
  if blocks < 2:
    result.stdev = 0.0
    result.hasStdev = false
