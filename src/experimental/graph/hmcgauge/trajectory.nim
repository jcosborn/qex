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
  TrajectoryGraph* = object
    initialState*: TrajectoryState
    finalState*: TrajectoryState
    deltaHamiltonian*: Gscalar
    acceptanceExpr*: Gscalar
    lossExpr*: Gscalar
    learnedParameters*: seq[LearnedParameter]

proc resampleMomentum*(graph: TrajectoryGraph, randomField: var auto) =
  let randomFieldPtr = addr randomField
  mutateGauge(graph.initialState.momentum, momentumStorage):
    threads:
      momentumStorage.randomTAH randomFieldPtr[]

proc buildTrajectoryState(gc: Gactcoeff,
                          gauge: Ggauge,
                          momentum: Ggauge): TrajectoryState =
  result.gauge = gauge
  result.momentum = momentum
  result.gaugeAction = gc.gaugeAction(gauge)
  result.kinetic = toGvalue(momentum.runtime, 0.5) * momentum.norm2
  result.hamiltonian = result.gaugeAction + result.kinetic

proc buildTrajectoryGraph*(grt: GraphRuntime,
                           g, p: graphGauge.Gauge,
                           gc: Gactcoeff,
                           config: RunConfig): TrajectoryGraph =
  let gdt = toGvalue(grt, config.dt)
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
  let deltaZero = toGvalue(result.deltaHamiltonian.runtime, 0.0)
  let acceptOne = toGvalue(result.deltaHamiltonian.runtime, 1.0)
  result.acceptanceExpr =
    cond(result.deltaHamiltonian < deltaZero, acceptOne, exp(-result.deltaHamiltonian))
  result.lossExpr = -result.acceptanceExpr * (tau * tau)
  # Pair each learned parameter (dt + integrator coefficients) with its gradient
  # expression. gradientExpr needs lossExpr, so it is filled once lossExpr exists.
  result.learnedParameters = @[LearnedParameter(name: "dt", node: gdt)]
  for c in integrated.learnedCoeffs:
    result.learnedParameters.add c
  for p in mitems(result.learnedParameters):
    p.gradientExpr = result.lossExpr.grad p.node

proc commitAcceptedTrajectory*(graph: TrajectoryGraph,
                               finalGauge: graphGauge.Gauge) =
  finalGauge.reunitGauge
  graph.initialState.gauge.update finalGauge
