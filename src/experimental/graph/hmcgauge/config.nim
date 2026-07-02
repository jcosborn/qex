import integrator
from ../core/base import raiseValueError

type
  TrajectoryPhase* = enum
    tpThermo, tpTrain, tpInfer
  RunConfig* = object
    gaugefile*, savefile*: string
    dt*, lrmax*, lrmin*, weightDecay*: float
    # Trajectory schedule shared by every gauge HMC: `trajsThermo` thermalization
    # trajectories then `trajs` production trajectories. The first `trajsForceAcc`
    # trajectories overall (counted from the first trajectory; default trajsThermo
    # div 8, hence within thermalization) are force-accepted to escape the hot start;
    # the first `trajsTrain` production trajectories do training (the rest inference).
    trajsThermo*, trajs*, trajsForceAcc*: int
    trajsTrain*, trajsTrainlrWarm*: int
    savefreq*, gsteps*: int
    revCheckFreq*: int   # run a reversibility check every this many trajectories (0=off)
    integratorCoeffs*: IntegratorCoeffs

proc validateRunConfig*(config: RunConfig) =
  if config.dt <= 0.0:
    raiseValueError("dt must be > 0, got " & $config.dt)
  if config.lrmin < 0.0:
    raiseValueError("lrmin must be >= 0, got " & $config.lrmin)
  if config.lrmax < 0.0:
    raiseValueError("lrmax must be >= 0, got " & $config.lrmax)
  if config.lrmin > config.lrmax:
    raiseValueError(
      "lrmin must be <= lrmax, got " & $config.lrmin & " > " & $config.lrmax)
  if config.weightDecay < 0.0:
    raiseValueError("weightDecay must be >= 0, got " & $config.weightDecay)
  if config.trajsThermo < 0:
    raiseValueError("trajsThermo must be >= 0, got " & $config.trajsThermo)
  if config.trajs < 0:
    raiseValueError("trajs must be >= 0, got " & $config.trajs)
  if config.trajsForceAcc < 0:
    raiseValueError("trajsForceAcc must be >= 0, got " & $config.trajsForceAcc)
  if config.trajsTrain < 0:
    raiseValueError("trajsTrain must be >= 0, got " & $config.trajsTrain)
  if config.trajsTrain > config.trajs:
    raiseValueError(
      "trajsTrain must be <= trajs, got " &
      $config.trajsTrain & " > " & $config.trajs)
  if config.trajsTrainlrWarm < 0:
    raiseValueError(
      "trajsTrainlrWarm must be >= 0, got " & $config.trajsTrainlrWarm)
  if config.trajsTrainlrWarm > config.trajsTrain:
    raiseValueError(
      "trajsTrainlrWarm must be <= trajsTrain, got " &
      $config.trajsTrainlrWarm & " > " & $config.trajsTrain)
  if config.savefreq < 0:
    raiseValueError("savefreq must be >= 0, got " & $config.savefreq)
  if config.gsteps < 1:
    raiseValueError("gsteps must be >= 1, got " & $config.gsteps)
  if config.revCheckFreq < 0:
    raiseValueError("revCheckFreq must be >= 0, got " & $config.revCheckFreq)

proc totalTrajs*(config: RunConfig): int =
  ## Total trajectories run: thermalization followed by production.
  config.trajsThermo + config.trajs

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
