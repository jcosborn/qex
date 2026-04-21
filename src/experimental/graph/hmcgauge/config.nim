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
    integratorCoeffs*: IntegratorCoeffs
    alwaysAccept*: bool

proc requireAtLeast(label: string,
                    value: int,
                    minimum: int) =
  if value < minimum:
    raiseValueError(label & " must be >= " & $minimum & ", got " & $value)

proc requirePositive(label: string,
                     value: float) =
  if value <= 0.0:
    raiseValueError(label & " must be > 0, got " & $value)

proc requireNonNegative(label: string,
                        value: float) =
  if value < 0.0:
    raiseValueError(label & " must be >= 0, got " & $value)

proc validateRunConfig*(config: RunConfig) =
  requirePositive("dt", config.dt)
  requireNonNegative("lrmin", config.lrmin)
  requireNonNegative("lrmax", config.lrmax)
  requireAtLeast("trajsThermo", config.trajsThermo, 0)
  requireAtLeast("trajsTrain", config.trajsTrain, 0)
  requireAtLeast("trajsTrainlrWarm", config.trajsTrainlrWarm, 0)
  requireAtLeast("trajsInfer", config.trajsInfer, 0)
  requireAtLeast("savefreq", config.savefreq, 0)
  requireAtLeast("gsteps", config.gsteps, 1)
  if config.lrmin > config.lrmax:
    raiseValueError(
      "lrmin must be <= lrmax, got " & $config.lrmin & " > " & $config.lrmax)
  if config.weightDecay < 0.0:
    raiseValueError("weightDecay must be >= 0, got " & $config.weightDecay)

proc totalTrajs*(config: RunConfig): int =
  config.trajsThermo + config.trajsTrain + config.trajsInfer

proc trajectoryPhase*(config: RunConfig,
                      traj: int): TrajectoryPhase =
  if traj <= config.trajsThermo:
    return tpThermo
  if traj <= config.trajsThermo + config.trajsTrain:
    return tpTrain
  tpInfer

proc phaseLabel*(config: RunConfig,
                 phase: TrajectoryPhase,
                 traj: int): string =
  case phase
  of tpThermo:
    "Thermalization step: " & $traj
  of tpTrain:
    "Training step: " & $(traj - config.trajsThermo)
  of tpInfer:
    "Inference step: " & $(traj - (config.trajsThermo + config.trajsTrain))
