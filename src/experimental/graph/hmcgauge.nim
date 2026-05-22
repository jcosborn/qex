import qex
import core
import scalar
import gauge
import gauge/shared as graphGauge
from hmcgauge/config import
  RunConfig, tpThermo, tpTrain, tpInfer,
  totalTrajs, trajectoryPhase, phaseLabel, validateRunConfig
from hmcgauge/trajectory import
  TrajectoryGraph, buildTrajectoryGraph, resampleMomentum,
  commitAcceptedTrajectory
from hmcgauge/training import
  TrainingState, initTrainingState, formatParameterValues, trainStep
from hmcgauge/gauge_io import loadOrInitGauge, maybeSaveGauge
from hmcgauge/params import readHmcGaugeInputs

proc echoRuntimeBanner() =
  echo "rank ", myRank, "/", nRanks
  threads:
    echo "thread ", threadNum, "/", numThreads

proc runTrajectory(random: var RngMilc6,
                   graph: TrajectoryGraph,
                   trainer: var TrainingState,
                   runConfig: RunConfig,
                   traj: int) =
  tic("traj")
  let phase = runConfig.trajectoryPhase(traj)
  echo runConfig.phaseLabel(phase, traj)
  echo "Begin H: ", graph.initialState.hamiltonian.eval.sval,
    "  Sg: ", graph.initialState.gaugeAction.eval.sval,
    "  T: ", graph.initialState.kinetic.eval.sval
  echo "End H: ", graph.finalState.hamiltonian.eval.sval,
    "  Sg: ", graph.finalState.gaugeAction.eval.sval,
    "  T: ", graph.finalState.kinetic.eval.sval
  let acceptDraw = random.uniform
  let deltaHamiltonian = graph.deltaHamiltonian.eval.sval
  let acceptanceProbability = graph.acceptanceExpr.eval.sval
  let accepted = acceptDraw <= acceptanceProbability or runConfig.alwaysAccept
  if accepted:
    echo "ACCEPT:  dH: ", deltaHamiltonian,
      "  exp(-dH): ", acceptanceProbability,
      "  r: ", acceptDraw,
      (if runConfig.alwaysAccept: " (ignored)" else: "")
  else:
    echo "REJECT:  dH: ", deltaHamiltonian,
      "  exp(-dH): ", acceptanceProbability,
      "  r: ", acceptDraw
  var acceptedGauge: graphGauge.Gauge
  if accepted:
    discard graph.finalState.gauge.eval
    acceptedGauge = graph.finalState.gauge.gaugeSnapshot
  qexGC "traj done"

  toc("forward end")
  let loss = graph.lossExpr.eval.sval

  case phase
  of tpThermo:
    echo "bloss: ", loss
  of tpTrain:
    echo "tloss: ", loss
    trainer.trainStep(runConfig, traj)
  of tpInfer:
    echo "iloss: ", loss

  if accepted:
    graph.commitAcceptedTrajectory(acceptedGauge)

  let currentGauge = graph.initialState.gauge.gaugeSnapshot
  currentGauge.echoPlaq
  currentGauge.maybeSaveGauge(runConfig, traj)
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

  installStandardParams()
  echoParams()
  processHelpParam()

  let
    runConfig = inputs.config
    lo = inputs.lat.newLayout
    grt = initGraphRuntime()
    gc = actWilson(scalar.toGvalue(grt, inputs.beta))

  var momentumRandom = lo.newRNGField(RngMilc6, inputs.seed)
  var acceptRandom: RngMilc6
  acceptRandom.seed(inputs.seed, 987654321)

  var
    g = lo.newgauge
    p = lo.newgauge

  g.loadOrInitGauge(runConfig.gaugefile)

  g.echoPlaq

  let graph = buildTrajectoryGraph(grt, g, p, gc, runConfig)
  var trainer = initTrainingState(graph, runConfig.weightDecay)

  echo trainer.formatParameterValues

  toc("prep")

  for traj in 1..runConfig.totalTrajs:
    graph.resampleMomentum(momentumRandom)
    runTrajectory(acceptRandom, graph, trainer, runConfig, traj)

  toc()

  processSaveParams()
  writeParamFile()

when isMainModule:
  runHmcGauge()
