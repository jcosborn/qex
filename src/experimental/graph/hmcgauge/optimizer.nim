from ../core/base import raiseValueError
from math import cos, PI, pow, sqrt

type
  AdamW* = object
    stepScale, beta1, beta2, epsilon, weightDecay: float
    firstMoment, secondMoment: seq[float]
  AdamStepStat* = object
    firstMoment*, secondMoment*: float

proc requireOptimizerShape(opt: AdamW,
                           param: openArray[float],
                           grad: openArray[float]) =
  if param.len != grad.len:
    raiseValueError(
      "optimizer parameter/gradient count mismatch: " &
      $param.len & " vs " & $grad.len)
  if opt.firstMoment.len != param.len:
    raiseValueError(
      "optimizer first-moment count mismatch: " &
      $opt.firstMoment.len & " vs " & $param.len)
  if opt.secondMoment.len != param.len:
    raiseValueError(
      "optimizer second-moment count mismatch: " &
      $opt.secondMoment.len & " vs " & $param.len)

proc requirePositive(label: string, value: float) =
  if value <= 0.0:
    raiseValueError(label & " must be > 0, got " & $value)

proc requireNonNegative(label: string, value: float) =
  if value < 0.0:
    raiseValueError(label & " must be >= 0, got " & $value)

proc requireMomentumCoeff(label: string, value: float) =
  if value < 0.0 or value >= 1.0:
    raiseValueError(label & " must satisfy 0 <= value < 1, got " & $value)

proc initAdamW*(param: openArray[float],
                stepScale = 0.001,
                beta1 = 0.9,
                beta2 = 0.999,
                epsilon = 1e-8,
                weightDecay = 0.01): AdamW =
  requirePositive("optimizer stepScale", stepScale)
  requireMomentumCoeff("optimizer beta1", beta1)
  requireMomentumCoeff("optimizer beta2", beta2)
  requirePositive("optimizer epsilon", epsilon)
  requireNonNegative("optimizer weightDecay", weightDecay)
  result = AdamW(
    stepScale: stepScale,
    beta1: beta1,
    beta2: beta2,
    epsilon: epsilon,
    weightDecay: weightDecay)
  result.firstMoment = newSeq[float](param.len)
  result.secondMoment = newSeq[float](param.len)

proc optimize*(opt: var AdamW,
               param: var seq[float],
               grad: openArray[float],
               t: int,
               lr: float): seq[AdamStepStat] =
  opt.requireOptimizerShape(param, grad)
  if t <= 0:
    raiseValueError("optimizer step must be >= 1, got " & $t)
  requireNonNegative("optimizer learning rate", lr)
  result = newSeq[AdamStepStat](grad.len)
  let
    stepScale = opt.stepScale
    b1 = opt.beta1
    b2 = opt.beta2
    sb1 = 1.0 - b1
    sb2 = 1.0 - b2
    sb1t = 1.0 - pow(b1, t.float)
    sb2t = 1.0 - pow(b2, t.float)
    weightDecay = opt.weightDecay
    epsilon = opt.epsilon
  for i in 0..<grad.len:
    opt.firstMoment[i] = b1 * opt.firstMoment[i] + sb1 * grad[i]
    opt.secondMoment[i] = b2 * opt.secondMoment[i] + sb2 * (grad[i] * grad[i])
    let firstMoment = opt.firstMoment[i] / sb1t
    let secondMoment = sqrt(opt.secondMoment[i] / sb2t)
    result[i] = AdamStepStat(
      firstMoment: firstMoment,
      secondMoment: secondMoment)
    let update = stepScale * firstMoment / (secondMoment + epsilon)
    param[i] = param[i] - lr * (update + weightDecay * param[i])

func warmUpCosDecay*(t,
                     twarm,
                     tmax: int,
                     lrmax: float,
                     lrmin = 0.0): float =
  if tmax <= 0:
    return lrmin

  let clampedT =
    if t < 0:
      0
    elif t > tmax:
      tmax
    else:
      t
  let warmSteps =
    if twarm < 0:
      0
    elif twarm > tmax:
      tmax
    else:
      twarm

  if warmSteps > 0 and clampedT <= warmSteps:
    return lrmin + (lrmax - lrmin) * clampedT.float / warmSteps.float

  if warmSteps >= tmax:
    return lrmax

  if warmSteps == 0:
    if clampedT <= 1 or tmax == 1:
      return lrmax
    return lrmin + 0.5 * (lrmax - lrmin) *
      (1.0 + cos(PI * float(clampedT - 1) / float(tmax - 1)))

  lrmin + 0.5 * (lrmax - lrmin) *
    (1.0 + cos(PI * float(clampedT - warmSteps) / float(tmax - warmSteps)))
