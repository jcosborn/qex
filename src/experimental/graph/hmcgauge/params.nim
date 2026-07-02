import qex
from os import fileExists
import config, rng
import integrator

type
  GaugeParams* = object
    ## Common gauge HMC inputs.
    lat*: seq[int]
    beta*: float
    seed*: uint
    rng*: RngKind
    gaugefile*, savefile*: string
    savefreq*: int
    dt*: float
    gsteps*: int
    intalg*: string
    intcoeffs*: seq[float]
    trajsThermo*: int      # thermalization trajectories (discarded)
    trajs*: int            # production trajectories (measured)
    trajsForceAcc*: int    # leading trajectories force-accepted
    revCheckFreq*: int     # reversibility-check period in trajectories (0 = off)
    jkBlockSize*: int      # jackknife block size for the analysis
    showTimers*: bool      # echo the tic/toc timer report at the end
    showRunStats*: bool    # echo per-node graph run statistics at the end
    # GC tuning, applied as process globals in toRunConfig (see there for why).
    eagerGc*: bool         # keep the field allocator's forced-GC threshold (default off lifts it)
    gcStats*: bool         # print GC statistics each trajectory (sets VerboseGCStats)
  TrainParams* = object
    ## Inputs used only when training. The first `trajsTrain` production
    ## trajectories do training; the rest are inference. Grows as more parameters
    ## become trainable.
    trajsTrain*, trajsTrainlrWarm*: int
    lrmax*, lrmin*, weightDecay*: float
  StoutParams* = object
    ## U = f(V) parameters.
    rho*: float
    nsmear*: int

proc readGaugeInputs*(d: GaugeParams): GaugeParams =
  ## Parse the common gauge HMC parameters, using `d` for defaults. The lattice
  ## comes from the gauge file when one is loaded. `trajsForceAcc` and `jkBlockSize`
  ## honor a positive `d` value but otherwise fall back to a count-derived default
  ## (trajsThermo div 8, force-accepted from the first trajectory; and max(1, trajs
  ## div 64)); both remain overridable on the command line.
  letParam:
    gaugefile = d.gaugefile
    savefile = d.savefile
    savefreq = d.savefreq
    lat =
      if fileExists(gaugefile):
        getFileLattice gaugefile
      else:
        if gaugefile.len > 0:
          qexWarn "Nonexistent gauge file: ", gaugefile
        d.lat
    beta = d.beta
    dt = d.dt
    gsteps = d.gsteps
    seed = d.seed
    rng = $d.rng
    gintalg = d.intalg
    gintcoeffs = d.intcoeffs
    trajsThermo = d.trajsThermo
    trajs = d.trajs
    trajsForceAcc = if d.trajsForceAcc > 0: d.trajsForceAcc else: trajsThermo div 8
    revCheckFreq = d.revCheckFreq
    jkBlockSize = if d.jkBlockSize > 0: d.jkBlockSize else: max(1, trajs div 64)
    showTimers = d.showTimers
    showRunStats = d.showRunStats
    eagerGc = d.eagerGc
    gcStats = d.gcStats
  GaugeParams(
    lat: lat,
    beta: beta,
    seed: seed,
    rng: parseRngKind(rng),
    gaugefile: gaugefile,
    savefile: savefile,
    savefreq: savefreq,
    dt: dt,
    gsteps: gsteps,
    intalg: gintalg,
    intcoeffs: gintcoeffs,
    trajsThermo: trajsThermo,
    trajs: trajs,
    trajsForceAcc: trajsForceAcc,
    revCheckFreq: revCheckFreq,
    jkBlockSize: jkBlockSize,
    showTimers: showTimers,
    showRunStats: showRunStats,
    eagerGc: eagerGc,
    gcStats: gcStats)

proc readTrainParams*(d: TrainParams): TrainParams =
  ## Parse the training-only command-line parameters, using `d` for defaults.
  letParam:
    trajsTrain = d.trajsTrain
    trajsTrainlrWarm = d.trajsTrainlrWarm
    lrmax = d.lrmax
    lrmin = d.lrmin
    weightDecay = d.weightDecay
  TrainParams(
    trajsTrain: trajsTrain,
    trajsTrainlrWarm: trajsTrainlrWarm,
    lrmax: lrmax,
    lrmin: lrmin,
    weightDecay: weightDecay)

proc readStoutParams*(d: StoutParams): StoutParams =
  ## Parse the stout field-transformation parameters, using `d` for defaults.
  letParam:
    rho = d.rho
    nsmear = d.nsmear
  StoutParams(rho: rho, nsmear: nsmear)

proc integratorCoeffs*(gp: GaugeParams): IntegratorCoeffs =
  parseIntegratorCoeffs(parseIntegratorKind(gp.intalg), gp.intcoeffs)

proc toRunConfig*(gp: GaugeParams): RunConfig =
  ## Build the run config and apply GC flags.
  ## Long-lived graph fields can stay above rmGcThreshold, making every allocation
  ## trigger a full GC. Unless eagerGc is set, disable that threshold.
  if not gp.eagerGc: setRawMemGcThreshold(int.high)
  VerboseGCStats = gp.gcStats
  RunConfig(
    gaugefile: gp.gaugefile,
    savefile: gp.savefile,
    savefreq: gp.savefreq,
    dt: gp.dt,
    gsteps: gp.gsteps,
    integratorCoeffs: gp.integratorCoeffs,
    trajsThermo: gp.trajsThermo,
    trajs: gp.trajs,
    trajsForceAcc: gp.trajsForceAcc,
    revCheckFreq: gp.revCheckFreq)

proc toRunConfig*(gp: GaugeParams, tp: TrainParams): RunConfig =
  ## Run configuration including the training schedule.
  result = gp.toRunConfig
  result.trajsTrain = tp.trajsTrain
  result.trajsTrainlrWarm = tp.trajsTrainlrWarm
  result.lrmax = tp.lrmax
  result.lrmin = tp.lrmin
  result.weightDecay = tp.weightDecay
