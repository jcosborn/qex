#RUNCMD $RUN1

import std/stats
import qex
import rng
import utils/resample
import utils/test

qexInit()

let mytest = newQEXTest("jackknife")

proc meanEst[D](xs: Ensemble[D]): float =
  var m: typeof(xs[0]) = 0
  let n = xs.len
  for i in 0..<n:
    let x = xs[i]
    m += (x-m)/float(i+1)
  m

proc meanErrEst[D](xs: Ensemble[D], bs: int): float =
  let jkstat = jackknife(xs, bs, meanEst)
  jkstat.stdev

let nconf = 1024
var xs = newSeq[float](nconf)
var r: MRG32k3a
r.seed(7654321, 1)
for j in 0..<xs.len:
  xs[j] = r.gaussian
let mean0 = mean(xs)
let stdev0 = sqrt(varianceS(xs)/float(nconf))
let sampletest = mytest.newTest("samplesize=" & $nconf)

proc testbs(bs: int) =
  let testbs = sampletest.newTest("blocksize=" & $bs)
  block:
    let jkstat = jackknife(xs, bs, meanEst)
    testbs.assertAlmostEqual(mean0, jkstat.mean)
    testbs.assertAlmostEqual(stdev0, jkstat.stdev, absTol=if bs == 1: 1e-13 else: 2e-3)
  block:
    let test2 = testbs.newTest("nested")
    let jkstat = jackknife(xs, bs, meanErrEst, bs)
    test2.assertAlmostEqual(stdev0, jkstat.mean, absTol=if bs == 1: 1e-13 else: 2e-3)
    test2.assertAlmostEqual(sqrt(2.0)/float(nconf), jkstat.stdev, absTol=2e-3)

testbs(1)
testbs(3)
testbs(8)

block:
  let test = mytest.newTest("autocorrelation positive window")
  let
    constant = @[2.0, 2.0, 2.0, 2.0]
    alternating = @[1.0, -1.0, 1.0, -1.0]
    onePositive = @[1.0, 1.0, -1.0, -1.0]
    onePositiveSmall = @[1.0e-8, 1.0e-8, -1.0e-8, -1.0e-8]
  test.assertAlmostEqual(1.0, constant.jackknife(1, intAutocorr).mean)
  test.assertAlmostEqual(1.0, alternating.jackknife(1, intAutocorr).mean)
  test.assertAlmostEqual(1.5, onePositive.jackknife(1, intAutocorr).mean)
  test.assertAlmostEqual(1.5, onePositiveSmall.jackknife(1, intAutocorr).mean)

block:
  let test = mytest.newTest("autocorrelation window past 200 lags")
  var x = newSeq[float](2048)
  for i in 0..<x.len:
    x[i] = sin(2.0*PI*float(i)/float(x.len))
  let tau = x.jackknife(x.len, intAutocorr).mean
  test.assertAlmostEqual(1.0, if tau > 500.0: 1.0 else: 0.0)

block:
  let test = mytest.newTest("weighted jackknife preserves original blocks")
  let
    changed = @[true, true, false, false, true, false, true, false]
    x = @[0.0, 0.0, 0.0, 0.0, 10.0, 10.0, 20.0, 20.0]
  for wantChanged in [false, true]:
    var w = newSeq[float](x.len)
    for i in 0..<x.len:
      if changed[i] == wantChanged:
        w[i] = 1.0
    let
      m = x.weightedJackknife(w, 2)
      r = x.weightedJackknife(w, 2, isRms = true)
    test.assertAlmostEqual(7.5, m.mean)
    test.assertAlmostEqual(7.386290792181598, m.stdev)
    test.assertAlmostEqual(11.180339887498949, r.mean)
    test.assertAlmostEqual(6.1708850248690785, r.stdev)

block:
  let test = mytest.newTest("weighted jackknife undefined one-block error")
  let x = @[1.0, 2.0, 3.0, 4.0]
  let s = x.weightedJackknife(@[1.0, 1.0, 0.0, 0.0], 2)
  test.assertAlmostEqual(1.5, s.mean)
  test.assertAlmostEqual(1.0, if classify(s.stdev) == fcNan: 1.0 else: 0.0)

block:
  let arTest = mytest.newTest("AR(1) sequence test")
  let
    alpha = 0.8
    n = 32768
    noiseVar = 1.0
    anaMean = 0.0
    anaVar = noiseVar/(1-alpha^2)
    anaLag1 = alpha
    anaIntAc = (1+alpha)/(1-alpha)

  var r: MRG32k3a
  r.seed(987654321, 1)
  var x = newSeq[float](n)
  for i in 1..<n:
    x[i] = alpha*x[i-1] + sqrt(noiseVar)*r.gaussian

  proc varEst[D](xs: Ensemble[D]): float =
    let m = meanEst(xs)
    var s = 0.0
    for i in 0..<xs.len:
      let d = xs[i]-m
      s += d*d
    s/float(xs.len-1)

  proc autocovariance[D](xs: Ensemble[D], lag: int): float =
    let n = xs.len
    if lag >= n:
      return 0.0
    let m = meanEst(xs)
    var c = 0.0
    for i in 0..<(n-lag):
      c += (xs[i]-m)*(xs[i+lag]-m)
    c/float(n-lag)

  let blockSize = 64
  let jkMean = jackknife(x, blockSize, meanEst)
  let jkVar = jackknife(x, blockSize, varEst)
  let jkLag1 = jackknife(x, blockSize, proc(xs: Ensemble[seq[float]]): float = autocovariance(xs, 1)/autocovariance(xs, 0))
  let jkIntAc = jackknife(x, blockSize, intAutocorr)

  arTest.assertAlmostEqual(anaMean, jkMean.mean, absTol=0.05)
  arTest.assertAlmostEqual(anaVar, jkVar.mean, absTol=0.1)
  arTest.assertAlmostEqual(anaLag1, jkLag1.mean, absTol=0.01)
  arTest.assertAlmostEqual(anaIntAc, jkIntAc.mean, absTol=1.0)

  echo "AR(1) parameter alpha     = ", alpha
  echo "Sequence length           = ", n
  echo "Block size used           = ", blockSize
  echo "Analytical mean           = ", anaMean
  echo "Jackknife mean estimate   = ", jkMean.mean, " ± ", jkMean.stdev
  echo "Analytical variance       = ", anaVar
  echo "Jackknife variance est    = ", jkVar.mean, " ± ", jkVar.stdev
  echo "Analytical lag-1 corr     = ", anaLag1
  echo "Jackknife lag-1 corr est  = ", jkLag1.mean, " ± ", jkLag1.stdev
  echo "Analytical int ac time    = ", anaIntAc
  echo "Jackknife int ac time est = ", jkIntAc.mean, " ± ", jkIntAc.stdev

mytest.qexFinalize
