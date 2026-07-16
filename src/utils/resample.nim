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
      #   zeta_i = (n/k_i) * E - ((n-k_i)/k_i) * E_{(i)}
      # Then bias = E - mean(zeta_i) and the variance estimate is
      #   Var(E) ≈ (1/(g (g-1))) * sum_i (zeta_i - mean(zeta))^2
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

proc intAutocorr*(xs: Ensemble[seq[float]]): float =
  ## rho(t) = sum_i dx_i*dx_(i+t)/(n*C(0)).
  ## tau = 1/2 + sum_t rho(t); stop at rho <= 0 or t >= 6*tau; return 2*tau.
  let n = xs.len
  if n < 2: return 1.0
  var m = 0.0
  for i in 0..<n:
    m += (xs[i]-m)/float(i+1)
  var c0 = 0.0
  for i in 0..<n:
    let d = xs[i]-m
    c0 += d*d
  c0 /= float(n)
  if c0 == 0.0: return 1.0
  var tau = 0.5
  for lag in 1..<n:
    var c = 0.0
    for i in 0..<(n-lag):
      c += (xs[i]-m)*(xs[i+lag]-m)
    let rho = c/(float(n)*c0)
    if rho <= 0.0: break
    tau += rho
    if float(lag) >= 6.0*tau: break
  2.0*tau

proc jackknife*[D,V,A](ensemble:D, blocksize:int, estimator:proc(x:Ensemble[D], arg:A):V, arg:A): auto =
  ## Perform grouped jackknife with blocksize
  ## The expectation value: estimator(ensemble, arg)
  let
    m = estimator(wrapEnsemble(addr ensemble), arg)
    n = ensemble.len
    g = (n+blocksize-1) div blocksize
  var jk = newseq[V](g)
  # Pseudo-value accumulation for unequal delete sizes k_i
  var zmean = V(0)  # mean of pseudo-values ζ_i
  var zM2   = V(0)  # sum of squared deviations for ζ_i
  for i in 0..<g:
    let j = jackknifeSample(addr ensemble, blocksize, i)
    let e = estimator(j, arg)
    jk[i] = e
    let skiplow = i*blocksize
    let ki = if skiplow + blocksize > n: n - skiplow else: blocksize
    let kif = float(ki)
    let zeta = (float(n)/kif)*m - (float(n-ki)/kif)*e
    let dz = zeta - zmean
    let dzn = dz / float(i+1)
    zmean += dzn
    zM2 += dz * dzn * float(i)
  JackknifeStat[V](mean:m, jksamples:jk,
    bias:m - zmean,
    stdev:sqrt(zM2 / (float(g) * float(max(1, g-1)))),
    hasStdev:g > 1)

proc jackknife*[D,V](ensemble:D, blocksize:int, estimator:proc(x:Ensemble[D]):V): auto =
  ## Perform grouped jackknife with blocksize
  ## The expectation value: estimator(ensemble)
  let
    m = estimator(wrapEnsemble(addr ensemble))
    n = ensemble.len
    g = (n+blocksize-1) div blocksize
  var jk = newseq[V](g)
  # Pseudo-value accumulation for unequal delete sizes k_i
  var zmean = V(0)  # mean of pseudo-values ζ_i
  var zM2   = V(0)  # sum of squared deviations for ζ_i
  for i in 0..<g:
    let j = jackknifeSample(addr ensemble, blocksize, i)
    let e = estimator(j)
    jk[i] = e
    let skiplow = i*blocksize
    let ki = if skiplow + blocksize > n: n - skiplow else: blocksize
    let kif = float(ki)
    let zeta = (float(n)/kif)*m - (float(n-ki)/kif)*e
    let dz = zeta - zmean
    let dzn = dz / float(i+1)
    zmean += dzn
    zM2 += dz * dzn * float(i)
  JackknifeStat[V](mean:m, jksamples:jk,
    bias:m - zmean,
    stdev:sqrt(zM2 / (float(g) * float(max(1, g-1)))),
    hasStdev:g > 1)

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
