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
  var jkm = V(0)
  var jk2 = V(0)
  for i in 0..<g:
    let j = jackknifeSample(addr ensemble, blocksize, i)
    let e = estimator(j, arg)
    let delta = e - jkm
    let deltan = delta / float(i+1)
    jkm += deltan
    jk2 += delta * deltan * float(i)
    jk[i] = e
  JackknifeStat[V](mean:m, jksamples:jk,
    bias:float(n-blocksize)/float(blocksize)*(jkm-m),
    stdev:sqrt(float(n-blocksize)*jk2/float(n)))

proc jackknife*[D,V](ensemble:D, blocksize:int, estimator:proc(x:Ensemble[D]):V): auto =
  ## Perform grouped jackknife with blocksize
  ## The expectation value: estimator(ensemble)
  let
    m = estimator(wrapEnsemble(addr ensemble))
    n = ensemble.len
    g = (n+blocksize-1) div blocksize
  var jk = newseq[V](g)
  var jkm = V(0)
  var jk2 = V(0)
  for i in 0..<g:
    let j = jackknifeSample(addr ensemble, blocksize, i)
    let e = estimator(j)
    let delta = e - jkm
    let deltan = delta / float(i+1)
    jkm += deltan
    jk2 += delta * deltan * float(i)
    jk[i] = e
  JackknifeStat[V](mean:m, jksamples:jk,
    bias:float(n-blocksize)/float(blocksize)*(jkm-m),
    stdev:sqrt(float(n-blocksize)*jk2/float(n)))

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

  mytest.qexFinalize
