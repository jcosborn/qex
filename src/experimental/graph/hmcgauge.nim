import qex
import core
import gauge
from hmcgauge/config import
  RunConfig, TrajectoryPhase, tpThermo, tpTrain, tpInfer,
  totalTrajs, trajectoryPhase, phaseLabel, validateRunConfig
from hmcgauge/trajectory import
  buildTrajectoryGraph, resampleMomentum, echoTrajectoryHamiltonians,
  shouldAccept, logAcceptance, lossValue, commitAcceptedTrajectory,
  currentGauge, finalGaugeSnapshot
from hmcgauge/training import
  initTrainingState, formatFloatValues, parameterValues, trainStep
from hmcgauge/gauge_io import loadOrInitGauge, maybeSaveGauge
from hmcgauge/params import
  installHmcGaugeParams, readHmcGaugeInputs, persistHmcGaugeParams

proc echoRuntimeBanner() =
  echo "rank ", myRank, "/", nRanks
  threads:
    echo "thread ", threadNum, "/", numThreads

proc runTrajectoryPhase(graph: auto,
                        trainer: var auto,
                        runConfig: RunConfig,
                        phase: TrajectoryPhase,
                        traj: int) =
  case phase
  of tpThermo:
    echo "bloss: ", graph.lossValue
  of tpTrain:
    trainer.trainStep(graph, runConfig, traj)
  of tpInfer:
    echo "iloss: ", graph.lossValue

proc finishTrajectory(graph: auto,
                      runConfig: RunConfig,
                      traj: int,
                      acceptedGauge: Ggauge) =
  if acceptedGauge != nil:
    graph.commitAcceptedTrajectory(acceptedGauge)

  graph.currentGauge.echoPlaq
  graph.currentGauge.maybeSaveGauge(runConfig, traj)

proc runTrajectory(random: var RngMilc6,
                   graph: auto,
                   trainer: var auto,
                   runConfig: RunConfig,
                   traj: int) =
  tic("traj")
  let phase = runConfig.trajectoryPhase(traj)
  echo runConfig.phaseLabel(phase, traj)
  graph.echoTrajectoryHamiltonians
  let acceptDraw = random.uniform
  let accepted = graph.shouldAccept(acceptDraw, runConfig.alwaysAccept)
  graph.logAcceptance(
    accepted,
    acceptDraw,
    runConfig.alwaysAccept)
  let acceptedGauge: Ggauge =
    if accepted:
      graph.finalGaugeSnapshot
    else:
      nil
  qexGC "traj done"

  toc("forward end")

  runTrajectoryPhase(graph, trainer, runConfig, phase, traj)
  finishTrajectory(graph, runConfig, traj, acceptedGauge)
  qexLog "traj ",traj," secs: ",getElapsedTime()
  toc("traj end")

proc runHmcGauge*() =
  qexInit()
  defer:
    qexFinalize()

  tic()
  let inputs = readHmcGaugeInputs()
  inputs.config.validateRunConfig

  echoRuntimeBanner()

  installHmcGaugeParams()

  let
    runConfig = inputs.config
    lo = inputs.lat.newLayout
    grt = initGraphRuntime()
    gc = actWilson(grt, inputs.beta)

  var r = lo.newRNGField(RngMilc6, inputs.seed)
  var R:RngMilc6  # global RNG
  R.seed(inputs.seed, 987654321)

  var
    g = lo.newgauge
    p = lo.newgauge

  g.loadOrInitGauge(runConfig.gaugefile)

  g.echoPlaq

  let graph = buildTrajectoryGraph(grt, g, p, gc, runConfig)
  var trainer = initTrainingState(graph, runConfig.weightDecay)

  echo formatFloatValues("param:", trainer.parameterValues)

  toc("prep")

  for traj in 1..runConfig.totalTrajs:
    graph.resampleMomentum(r)
    runTrajectory(R, graph, trainer, runConfig, traj)

  toc()

  persistHmcGaugeParams()

when isMainModule:
  runHmcGauge()
