import qex
import ../core
import ../scalar
import config, trajectory
import optimizer

type
  TrainingState* = object
    learned: seq[LearnedParameter]
    optimizer: AdamW

proc formatFloatValues*(label: string, values: openArray[float]): string =
  result = label
  for value in values:
    result &= " " & $value

proc parameterValues*(training: TrainingState): seq[float] =
  result = newSeq[float](training.learned.len)
  for i in 0..<result.len:
    result[i] = training.learned[i].node.getfloat

proc gradientValues(training: TrainingState): seq[float] =
  result = newSeq[float](training.learned.len)
  for i in 0..<result.len:
    result[i] = training.learned[i].gradientExpr.eval.getfloat

proc applyParameterValues(training: var TrainingState,
                          values: openArray[float]) =
  if values.len != training.learned.len:
    raiseValueError(
      "training parameter count mismatch: expected " &
      $training.learned.len & ", got " & $values.len)
  for i in 0..<values.len:
    training.learned[i].node.update values[i]

proc logOptimizerStats(training: TrainingState,
                       stats: openArray[AdamStepStat]) =
  ## Log per-parameter optimizer state without adding logging to the optimizer.
  let label = training.optimizer.optimizerName
  for i in 0..<stats.len:
    echo label, ": ", i, " ", stats[i].firstMoment, " ", stats[i].secondMoment

proc initTrainingState*(graph: TrajectoryGraph,
                        weightDecay: float): TrainingState =
  result.learned = graph.learnedParameters
  result.optimizer = initAdamW(result.parameterValues, weightDecay = weightDecay)

proc trainStep*(training: var TrainingState,
                graph: TrajectoryGraph,
                config: RunConfig,
                traj: int) =
  tic()
  let trainingStep = traj - config.trajsThermo
  echo "tloss: ", graph.lossValue
  let gradients = training.gradientValues
  echo formatFloatValues("grad:", gradients)
  var parameters = training.parameterValues
  let learningRate = warmUpCosDecay(
    trainingStep,
    config.trajsTrainlrWarm,
    config.trajsTrain,
    config.lrmax,
    config.lrmin)
  echo "lr: ", learningRate
  let optimizerStats = training.optimizer.optimize(
    parameters,
    gradients,
    trainingStep,
    learningRate)
  training.logOptimizerStats(optimizerStats)
  training.applyParameterValues(parameters)
  echo formatFloatValues("param:", parameters)
  toc("training")
