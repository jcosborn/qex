import qex
import ../[core, scalar, gauge]
import ../gauge/shared as graphGauge
import config, integrator

type
  TrajectoryState* = object
    gauge*: Ggauge
    momentum*: Ggauge
    gaugeAction*: Gscalar
    kinetic*: Gscalar
    hamiltonian*: Gscalar
  LearnedParameter* = object
    name*: string
    node*: Gscalar
    gradientExpr*: Gscalar
  TrajectoryGraph* = object
    initialState*: TrajectoryState
    finalState*: TrajectoryState
    deltaHamiltonian*: Gscalar
    acceptanceExpr*: Gscalar
    lossExpr*: Gscalar
    learnedParameters*: seq[LearnedParameter]

proc initLearnedParameters(lossExpr: Gscalar,
                           parameters: openArray[LearnedCoeff]): seq[LearnedParameter] =
  result = newSeq[LearnedParameter](parameters.len)
  for i in 0..<result.len:
    result[i] = LearnedParameter(
      name: parameters[i].name,
      node: parameters[i].node,
      gradientExpr: lossExpr.grad parameters[i].node)

proc resampleMomentum*(graph: TrajectoryGraph, randomField: var auto) =
  let randomFieldPtr = addr randomField
  mutateGauge(graph.initialState.momentum, momentumStorage):
    threads:
      momentumStorage.randomTAH randomFieldPtr[]

proc buildTrajectoryState(gc: Gactcoeff,
                          gauge: Ggauge,
                          momentum: Ggauge): TrajectoryState =
  let kineticScale = scalarLeafLike(momentum, 0.5)
  result.gauge = gauge
  result.momentum = momentum
  result.gaugeAction = gc.gaugeAction(gauge)
  result.kinetic = kineticScale * momentum.norm2
  result.hamiltonian = result.gaugeAction + result.kinetic

proc buildTrajectoryGraph*(grt: GraphRuntime,
                           g, p: graphGauge.Gauge,
                           gc: Gactcoeff,
                           config: RunConfig): TrajectoryGraph =
  let gdt = toGvalue(grt, config.dt)
  var parameters = @[LearnedCoeff(name: "dt", node: gdt)]
  result.initialState = buildTrajectoryState(gc, toGvalue(grt, g), toGvalue(grt, p))
  let tau = float(config.gsteps) * gdt
  let integrated = integrateGauge(
    gc,
    result.initialState.gauge,
    result.initialState.momentum,
    gdt,
    config.gsteps,
    config.integratorCoeffs)
  result.finalState = buildTrajectoryState(
    gc,
    integrated.gauge,
    integrated.momentum)
  result.deltaHamiltonian =
    result.finalState.hamiltonian - result.initialState.hamiltonian
  let deltaZero = scalarLeafLike(result.deltaHamiltonian, 0.0)
  let acceptOne = scalarLeafLike(result.deltaHamiltonian, 1.0)
  result.acceptanceExpr =
    cond(result.deltaHamiltonian < deltaZero, acceptOne, exp(-result.deltaHamiltonian))
  result.lossExpr = -result.acceptanceExpr * (tau * tau)
  parameters.add integrated.learnedCoeffs
  result.learnedParameters = initLearnedParameters(result.lossExpr, parameters)

proc commitAcceptedTrajectory*(graph: TrajectoryGraph,
                               finalGauge: graphGauge.Gauge) =
  graph.initialState.gauge.update finalGauge
  graph.initialState.gauge.unsafeGaugeStorage.reunitGauge
