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
    stdev:sqrt(zM2 / (float(g) * float(max(1, g-1)))))

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
    stdev:sqrt(zM2 / (float(g) * float(max(1, g-1)))))

when isMainModule:
  import std/stats
  import qex
  import rng
  import utils/test

  qexInit()

  let mytest = newQEXTest("jackknife")

  proc meanEst[D](xs:Ensemble[D]):float =
    var m:typeof(xs[0]) = 0
    let n = xs.len
    for i in 0..<n:
      let x = xs[i]
      m += (x-m)/float(i+1)
    m

  proc meanErrEst[D](xs:Ensemble[D], bs:int):float =
    let jkstat = jackknife(xs, bs, meanEst)
    jkstat.stdev

  let nconf = 1024
  var xs = newseq[float](nconf)
  var r: MRG32k3a
  r.seed(7654321, 1)
  for j in 0..<xs.len:
    xs[j] = r.gaussian
  let mean0 = mean(xs)
  let stdev0 = sqrt(varianceS(xs)/float(nconf))
  let sampletest = mytest.newTest("samplesize=" & $nconf)

  proc testbs(bs:int) =
    let testbs = sampletest.newTest("blocksize=" & $bs)
    block:
      let jkstat = jackknife(xs, bs, meanEst)
      testbs.assertAlmostEqual(mean0, jkstat.mean)
      testbs.assertAlmostEqual(stdev0, jkstat.stdev, absTol=if bs==1: 1e-13 else: 2e-3)
    block:
      let test2 = testbs.newTest("nested")
      let jkstat = jackknife(xs, bs, meanErrEst, bs)
      test2.assertAlmostEqual(stdev0, jkstat.mean, absTol=if bs==1: 1e-13 else: 2e-3)
      test2.assertAlmostEqual(sqrt(2.0)/float(nconf), jkstat.stdev, absTol=2e-3)

  testbs(1)
  testbs(3)
  testbs(8)

  # -------------------------------------------------------------
  # A small AR(1) test
  # -------------------------------------------------------------
  block:
    let arTest = mytest.newTest("AR(1) sequence test")

    # AR(1) parameters
    let alpha = 0.8    # correlation coefficient
    let N = 32768      # length of sequence
    let noiseVar = 1.0 # variance of the driving Gaussian
    # Analytical results for an AR(1) of the form
    # x_{n+1} = alpha * x_n + eps_n,    eps_n ~ Normal(0, noiseVar)
    # mean = 0
    # variance = noiseVar / (1 - alpha^2)
    # lag-1 autocorr = alpha
    # integrated autocorr time (assuming alpha > 0) = (1 + alpha) / (1 - alpha)

    let anaMean = 0.0
    let anaVar = noiseVar / (1 - alpha^2)
    let anaLag1 = alpha
    let anaIntAc = (1 + alpha) / (1 - alpha)

    # Generate the AR(1) sequence
    var r: MRG32k3a
    r.seed(987654321, 1)
    var x = newSeq[float](N)
    for i in 1..<N:
      x[i] = alpha * x[i-1] + sqrt(noiseVar) * r.gaussian

    # ---- Define local estimators that take an Ensemble ----

    proc varEst[D](xs: Ensemble[D]): float =
      ## Unbiased sample variance
      let m = meanEst(xs)
      var s = 0.0
      for i in 0..<xs.len:
        let diff = xs[i] - m
        s += diff * diff
      s / float(xs.len - 1)

    proc autocovariance[D](xs: Ensemble[D], lag: int): float =
      ## Compute sample autocovariance at given lag
      let n = xs.len
      if lag >= n:
        return 0.0
      let m = meanEst(xs)
      var c = 0.0
      for i in 0..<(n - lag):
        c += (xs[i] - m) * (xs[i + lag] - m)
      c / float(n - lag)

    proc intAutocorr[D](xs: Ensemble[D], maxLag: int): float =
      ## Very rough estimator of integrated autocorrelation length:
      ## sum rho(k) for k=0..maxLag, where rho(k) = C(k)/C(0).
      let c0 = autocovariance(xs, 0)
      if c0 == 0.0:
        return 1.0
      var sumRho = 1.0  # start at lag=0
      for k in 1..maxLag:
        let ck = autocovariance(xs, k)
        let rho = ck / c0
        # a naive stopping criterion:
        if rho < 0.0:
          break
        sumRho += 2.0 * rho
      sumRho

    # ---- Now do jackknife for each statistic ----
    # Choose a blocksize; in practice, it depends on correlation.
    let blockSize = 64

    let jkMean = jackknife(x, blockSize, meanEst)
    let jkVar  = jackknife(x, blockSize, varEst)
    let jkLag1 = jackknife(x, blockSize, proc(xs:Ensemble[seq[float]]):float=autocovariance(xs, 1)/autocovariance(xs, 0))
    let jkIntAc = jackknife(x, blockSize, proc(xs:Ensemble[seq[float]]):float=intAutocorr(xs, 200))

    # Compare with analytics
    arTest.assertAlmostEqual(anaMean, jkMean.mean, absTol=0.05)
    arTest.assertAlmostEqual(anaVar,  jkVar.mean,  absTol=0.1)
    arTest.assertAlmostEqual(anaLag1, jkLag1.mean, absTol=0.01)
    arTest.assertAlmostEqual(anaIntAc, jkIntAc.mean, absTol=1.0)

    # Print everything
    echo "AR(1) parameter alpha     = ", alpha
    echo "Sequence length           = ", N
    echo "Block size used           = ", blockSize
    echo "Analytical mean           = ", anaMean
    echo "Jackknife mean estimate   = ", jkMean.mean, " ± ", jkMean.stdev
    echo "Analytical variance       = ", anaVar
    echo "Jackknife variance est    = ", jkVar.mean,  " ± ", jkVar.stdev
    echo "Analytical lag-1 corr     = ", anaLag1
    echo "Jackknife lag-1 corr est  = ", jkLag1.mean, " ± ", jkLag1.stdev
    echo "Analytical int ac time    = ", anaIntAc
    echo "Jackknife int ac time est = ", jkIntAc.mean, " ± ", jkIntAc.stdev
    echo "-------------------------------------------------------"

  mytest.qexFinalize
