import qex
import ../core
import ../scalar
import config, trajectory
import optimizer

type
  TrainingState* = object
    learned*: seq[LearnedParameter]
    optimizer: AdamW

proc formatNamedValues(training: TrainingState,
                       label: string,
                       values: openArray[float]): string =
  if values.len != training.learned.len:
    raiseValueError(
      "training value count mismatch: expected " &
      $training.learned.len & ", got " & $values.len)
  result = label
  for i in 0..<values.len:
    result &= " " & training.learned[i].name & "=" & $values[i]

proc parameterValues*(training: TrainingState): seq[float] =
  result = newSeq[float](training.learned.len)
  for i in 0..<result.len:
    result[i] = training.learned[i].node.sval

proc gradientValues(training: TrainingState): seq[float] =
  result = newSeq[float](training.learned.len)
  for i in 0..<result.len:
    result[i] = training.learned[i].gradientExpr.eval.sval

proc applyParameterValues(training: var TrainingState,
                          values: openArray[float]) =
  if values.len != training.learned.len:
    raiseValueError(
      "training parameter count mismatch: expected " &
      $training.learned.len & ", got " & $values.len)
  for i in 0..<values.len:
    training.learned[i].node.update values[i]

proc initTrainingState*(graph: TrajectoryGraph,
                        weightDecay: float): TrainingState =
  result.learned = graph.learnedParameters
  result.optimizer = initAdamW(result.parameterValues, weightDecay = weightDecay)

proc formatParameterValues*(training: TrainingState): string =
  training.formatNamedValues("param:", training.parameterValues)

proc trainStep*(training: var TrainingState,
                config: RunConfig,
                traj: int) =
  tic()
  let trainingStep = traj - config.trajsThermo
  let gradients = training.gradientValues
  echo training.formatNamedValues("grad:", gradients)
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
  let optimizerLabel =
    if config.weightDecay == 0.0:
      "Adam"
    else:
      "AdamW"
  for i in 0..<optimizerStats.len:
    echo optimizerLabel, ": ", training.learned[i].name, " ",
      optimizerStats[i].firstMoment, " ", optimizerStats[i].secondMoment
  training.applyParameterValues(parameters)
  echo training.formatNamedValues("param:", parameters)
  toc("training")
