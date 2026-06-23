import integrator
from ../core/base import raiseValueError

type
  TrajectoryPhase* = enum
    tpThermo, tpTrain, tpInfer
  RunConfig* = object
    gaugefile*, savefile*: string
    dt*, lrmax*, lrmin*, weightDecay*: float
    trajsThermo*, trajsTrain*, trajsTrainlrWarm*, trajsInfer*: int
    savefreq*, gsteps*: int
    alwaysAccept*: bool
    integratorCoeffs*: IntegratorCoeffs

proc validateRunConfig*(config: RunConfig) =
  if config.dt <= 0.0:
    raiseValueError("dt must be > 0, got " & $config.dt)
  if config.lrmin < 0.0:
    raiseValueError("lrmin must be >= 0, got " & $config.lrmin)
  if config.lrmax < 0.0:
    raiseValueError("lrmax must be >= 0, got " & $config.lrmax)
  if config.trajsThermo < 0:
    raiseValueError("trajsThermo must be >= 0, got " & $config.trajsThermo)
  if config.trajsTrain < 0:
    raiseValueError("trajsTrain must be >= 0, got " & $config.trajsTrain)
  if config.trajsTrainlrWarm < 0:
    raiseValueError(
      "trajsTrainlrWarm must be >= 0, got " & $config.trajsTrainlrWarm)
  if config.trajsTrainlrWarm > config.trajsTrain:
    raiseValueError(
      "trajsTrainlrWarm must be <= trajsTrain, got " &
      $config.trajsTrainlrWarm & " > " & $config.trajsTrain)
  if config.trajsInfer < 0:
    raiseValueError("trajsInfer must be >= 0, got " & $config.trajsInfer)
  if config.savefreq < 0:
    raiseValueError("savefreq must be >= 0, got " & $config.savefreq)
  if config.gsteps < 1:
    raiseValueError("gsteps must be >= 1, got " & $config.gsteps)
  if config.lrmin > config.lrmax:
    raiseValueError(
      "lrmin must be <= lrmax, got " & $config.lrmin & " > " & $config.lrmax)
  if config.weightDecay < 0.0:
    raiseValueError("weightDecay must be >= 0, got " & $config.weightDecay)

proc totalTrajs*(config: RunConfig): int =
  config.trajsThermo + config.trajsTrain + config.trajsInfer

proc trajectoryPhase*(config: RunConfig,
                      traj: int): TrajectoryPhase =
  let total = config.totalTrajs
  if traj < 1 or traj > total:
    raiseValueError(
      "trajectory index must satisfy 1 <= traj <= " & $total &
      ", got " & $traj)
  if traj <= config.trajsThermo:
    return tpThermo
  if traj <= config.trajsThermo + config.trajsTrain:
    return tpTrain
  tpInfer
