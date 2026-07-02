import qex
import core
import scalar
import gauge
import gauge/shared as graphGauge
from hmcgauge/config import
  RunConfig, tpThermo, tpTrain, tpInfer,
  totalTrajs, trajectoryPhase, validateRunConfig
from hmcgauge/trajectory import
  TrajectoryGraph, MdForceStats, buildTrajectoryGraph, runHmc
from hmcgauge/rng import rkPhilox4x64, withRng
from hmcgauge/training import
  TrainingState, initTrainingState, formatParameterValues, trainStep
from hmcgauge/gauge_io import loadOrInitGauge, maybeSaveGauge
from hmcgauge/params import
  GaugeParams, TrainParams, readGaugeInputs, readTrainParams, toRunConfig

proc echoRuntimeBanner() =
  echo "rank ", myRank, "/", nRanks
  threads:
    echo "thread ", threadNum, "/", numThreads

proc runHmcGauge*() =
  qexInit()
  defer:
    qexFinalize()

  tic()
  let
    gp = readGaugeInputs(GaugeParams(
      lat: @[8, 8, 8, 16], beta: 5.4, seed: 1234567891, rng: rkPhilox4x64,
      savefile: "config", dt: 0.025, gsteps: 4, intalg: "2MN", trajs: 50))
    tp = readTrainParams(TrainParams(
      lrmax: 1.0, lrmin: 0.0001, trajsTrain: 50, trajsTrainlrWarm: 10))
    runConfig = toRunConfig(gp, tp)
  runConfig.validateRunConfig

  echoRuntimeBanner()

  installStandardParams()
  echoParams()
  processHelpParam()

  let grt = initGraphRuntime()
  withRng(gp.rng, R):
    let
      lo = gp.lat.newLayout
      gc = actWilson(scalar.toGvalue(grt, gp.beta))

    var
      randomField = lo.newRNGField(R, gp.seed)
      acceptRandom: R
      g = lo.newgauge
      p = lo.newgauge
    acceptRandom.seed(gp.seed, 987654321)

    g.loadOrInitGauge(runConfig.gaugefile)

    g.echoPlaq

    let action = proc(x: Ggauge): Gscalar = gaugeAction(gc, x)
    let graph =
      if runConfig.trajsTrain == 0:
        let force = proc(x: Ggauge): Ggauge = gaugeForce(gc, x)
        buildTrajectoryGraph(grt, g, p, action, runConfig, force = force)
      else:
        buildTrajectoryGraph(grt, g, p, action, runConfig)
    var trainer = initTrainingState(graph, runConfig.weightDecay)

    echo trainer.formatParameterValues

    toc("prep")

    # Pre-commit (training) and post-commit (measurement) hooks for the shared driver.
    proc onProposal(traj: int; dH, acc: float) =
      let loss = graph.lossExpr.eval.sval
      case runConfig.trajectoryPhase(traj)
      of tpThermo: echo "bloss: ", loss
      of tpTrain:
        echo "tloss: ", loss
        trainer.trainStep(runConfig, traj - runConfig.trajsThermo)
      of tpInfer: echo "iloss: ", loss
    proc measureTraj(traj: int; dH, acc: float; accepted: bool; forceStats: MdForceStats) =
      let currentGauge = graph.initialState.gauge.gaugeSnapshot
      currentGauge.echoPlaq
      currentGauge.maybeSaveGauge(runConfig, traj)
    runHmc(graph, runConfig, randomField, acceptRandom, measureTraj, onProposal)

  toc()

  if gp.showTimers: echoTimers()
  if gp.showRunStats: grt.echoRunStats
  processSaveParams()
  writeParamFile()

when isMainModule:
  runHmcGauge()
