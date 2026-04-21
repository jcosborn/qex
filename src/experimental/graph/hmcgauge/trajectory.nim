import qex
import ../[core, scalar, gauge]
import ../gauge/shared as graphGauge
import config, integrator

type
  TrajectoryState = object
    gauge: Gvalue
    momentum: Gvalue
    gaugeAction: Gvalue
    kinetic: Gvalue
    hamiltonian: Gvalue
  LearnedParameter* = object
    node*: Gvalue
    gradientExpr*: Gvalue
  TrajectoryGraph* = object
    initialState: TrajectoryState
    finalState: TrajectoryState
    deltaHamiltonian: Gvalue
    acceptanceExpr: Gvalue
    lossExpr: Gvalue
    learnedParameters*: seq[LearnedParameter]

proc initLearnedParameters(lossExpr: Gvalue,
                           parameterNodes: openArray[Gvalue]): seq[LearnedParameter] =
  result = newSeq[LearnedParameter](parameterNodes.len)
  for i in 0..<result.len:
    result[i] = LearnedParameter(
      node: parameterNodes[i],
      gradientExpr: lossExpr.grad parameterNodes[i])

proc lossValue*(graph: TrajectoryGraph): float =
  graph.lossExpr.eval.getfloat

proc acceptanceProbability*(graph: TrajectoryGraph): float =
  graph.acceptanceExpr.eval.getfloat

proc currentGauge*(graph: TrajectoryGraph): graphGauge.Gauge =
  graph.initialState.gauge.getgauge

proc currentMomentum*(graph: TrajectoryGraph): graphGauge.Gauge =
  graph.initialState.momentum.getgauge

template resampleMomentum*(graphValue, randomFieldValue: untyped) =
  block:
    threads:
      graphValue.initialState.momentum.getgauge.randomTAH randomFieldValue
    graphValue.initialState.momentum.updated

proc buildTrajectoryState(gc: Gvalue,
                          gauge: Gvalue,
                          momentum: Gvalue): TrajectoryState =
  let kineticScale = scalarLeafLike(momentum, 0.5)
  result.gauge = gauge
  result.momentum = momentum
  result.gaugeAction = gc.gaugeAction(gauge)
  result.kinetic = kineticScale * momentum.norm2
  result.hamiltonian = result.gaugeAction + result.kinetic

proc logTrajectoryState(label: string,
                        state: TrajectoryState) =
  echo label, " H: ", state.hamiltonian.eval,
    "  Sg: ", state.gaugeAction.eval,
    "  T: ", state.kinetic.eval

proc buildTrajectoryGraph*(grt: GraphRuntime,
                           g, p: auto,
                           gc: Gvalue,
                           config: RunConfig): TrajectoryGraph =
  let gdt = toGvalue(grt, config.dt)
  var parameterNodes = @[Gvalue(gdt)]
  result.initialState = buildTrajectoryState(gc, toGvalue(grt, g), toGvalue(grt, p))
  let tau = float(config.gsteps) * gdt
  let (g1, p1, learnedCoeffs) = integrateGauge(
    gc,
    result.initialState.gauge,
    result.initialState.momentum,
    gdt,
    config.gsteps,
    config.integratorCoeffs)
  result.finalState = buildTrajectoryState(gc, g1, p1)
  result.deltaHamiltonian =
    result.finalState.hamiltonian - result.initialState.hamiltonian
  let deltaZero = scalarLeafLike(result.deltaHamiltonian, 0.0)
  let acceptOne = scalarLeafLike(result.deltaHamiltonian, 1.0)
  result.acceptanceExpr =
    cond(result.deltaHamiltonian < deltaZero, acceptOne, exp(-result.deltaHamiltonian))
  result.lossExpr = -result.acceptanceExpr * (tau * tau)
  parameterNodes.add learnedCoeffs
  result.learnedParameters = initLearnedParameters(result.lossExpr, parameterNodes)

proc echoTrajectoryHamiltonians*(graph: TrajectoryGraph) =
  logTrajectoryState("Begin", graph.initialState)
  logTrajectoryState("End", graph.finalState)

proc logAcceptance*(graph: TrajectoryGraph,
                    accepted: bool,
                    draw: float,
                    alwaysAccept: bool) =
  if accepted:
    echo "ACCEPT:  dH: ", graph.deltaHamiltonian.eval,
      "  exp(-dH): ", graph.acceptanceExpr.eval,
      "  r: ", draw,
      (if alwaysAccept: " (ignored)" else: "")
  else:
    echo "REJECT:  dH: ", graph.deltaHamiltonian.eval,
      "  exp(-dH): ", graph.acceptanceExpr.eval,
      "  r: ", draw

proc shouldAccept*(graph: TrajectoryGraph,
                   draw: float,
                   alwaysAccept: bool): bool =
  draw <= graph.acceptanceProbability or alwaysAccept

proc commitAcceptedTrajectory*(graph: TrajectoryGraph) =
  graph.initialState.gauge.valCopy graph.finalState.gauge
  graph.initialState.gauge.getgauge.reunitGauge
  graph.initialState.gauge.updated
