## HMC statistics shared across gauge theories.

import math
import utils/resample
import integrator

proc mean*(xs: Ensemble[seq[float]]): float =
  var m = 0.0
  for i in 0..<xs.len:
    let x = xs[i]
    m += (x - m)/float(i+1)
  m

proc rms*(xs: Ensemble[seq[float]]): float =
  var m = 0.0
  for i in 0..<xs.len:
    let x2 = xs[i]*xs[i]
    m += (x2 - m)/float(i+1)
  sqrt(m)

proc extrema*(xs: openArray[float]): tuple[lo, hi: float] =
  result = (xs[0], xs[0])
  for i in 1..<xs.len:
    result.lo = min(result.lo, xs[i])
    result.hi = max(result.hi, xs[i])

proc echoMdStats*(xs: seq[MdForceStats]; blockSize: int; tag: string = ""; selected: seq[bool] = @[]) =
  var
    count = 0
    rmsMax = 0.0
    fminMin = 1e300
    fmaxMax = 0.0
    rmsVals = newSeq[float](xs.len)
    fminVals = newSeq[float](xs.len)
    fmaxVals = newSeq[float](xs.len)
    weights = newSeq[float](xs.len)
  for i, x in xs:
    rmsVals[i] = x.rmsMean
    fminVals[i] = x.fminMean
    fmaxVals[i] = x.fmaxMean
    if selected.len == 0 or selected[i]:
      weights[i] = float(x.count)
      if x.count > 0:
        count += x.count
        rmsMax = max(rmsMax, x.rmsMax)
        fminMin = min(fminMin, x.fminMin)
        fmaxMax = max(fmaxMax, x.fmaxMax)
  if count == 0: return
  proc je(vals: seq[float]): string =
    let x = vals.weightedJackknife(weights, blockSize)
    $x.mean & " ± " & (if classify(x.stdev) == fcNan: "n/a" else: $x.stdev)
  let prefix = if tag.len == 0: "" else: tag & " "
  echo prefix, "MD forces: ", count
  echo prefix, "fRMS mean: ", je(rmsVals)
  echo prefix, "fRMS max: ", rmsMax
  echo prefix, "fMin mean: ", je(fminVals)
  echo prefix, "fMin min: ", fminMin
  echo prefix, "fMax mean: ", je(fmaxVals)
  echo prefix, "fMax max: ", fmaxMax
