import qex
import ../[core, scalar, gauge]
import ../gauge/shared as graphGauge
import config, integrator
export GaugeAction, GaugeForceFn, MdForceStats

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
    lossExpr*: Gscalar
    learnedParameters*: seq[LearnedParameter]
    mdForces*: seq[Ggauge]

proc resampleMomentum*(graph: TrajectoryGraph, randomField: var auto) =
  let randomFieldPtr = addr randomField
  mutateGauge(graph.initialState.momentum, momentumStorage):
    threads:
      momentumStorage.randomTAH randomFieldPtr[]

proc buildTrajectoryState(action: GaugeAction,
                          gauge: Ggauge,
                          momentum: Ggauge): TrajectoryState =
  result.gauge = gauge
  result.momentum = momentum
  result.gaugeAction = action(gauge)
  result.kinetic = toGvalue(momentum.runtime, 0.5) * momentum.norm2
  result.hamiltonian = result.gaugeAction + result.kinetic

proc buildTrajectoryGraph*(grt: GraphRuntime,
                           g, p: graphGauge.Gauge,
                           action: GaugeAction,
                           config: RunConfig,
                           buildTraining = true,
                           force: GaugeForceFn = nil): TrajectoryGraph =
  ## Build HMC from S(g); the same S defines the force and both Hamiltonians.
  ## If not buildTraining, omit the loss and parameter-gradient graph.
  let gdt = toGvalue(grt, config.dt)
  result.initialState = buildTrajectoryState(action, toGvalue(grt, g), toGvalue(grt, p))
  let integrated = integrateGauge(
    action,
    result.initialState.gauge,
    result.initialState.momentum,
    gdt,
    config.gsteps,
    config.integratorCoeffs,
    force)
  result.mdForces = integrated.forces
  result.finalState = buildTrajectoryState(
    action,
    integrated.gauge,
    integrated.momentum)
  result.deltaHamiltonian =
    result.finalState.hamiltonian - result.initialState.hamiltonian
  if buildTraining:
    # Acceptance probability min(1, exp(-dH)) — the training reward. (The sampler's
    # accept test in runHmc uses the raw, uncapped exp(-dH) directly, which is also
    # the ⟨exp(-dH)⟩ diagnostic, so this capped expression is training-only.)
    let deltaZero = toGvalue(result.deltaHamiltonian.runtime, 0.0)
    let acceptOne = toGvalue(result.deltaHamiltonian.runtime, 1.0)
    let acceptanceExpr =
      cond(result.deltaHamiltonian < deltaZero, acceptOne, exp(-result.deltaHamiltonian))
    let tau = float(config.gsteps) * gdt
    result.lossExpr = -acceptanceExpr * (tau * tau)
    # Pair each learned parameter (dt + integrator coefficients) with its gradient
    # expression. gradientExpr needs lossExpr, so it is filled once lossExpr exists.
    result.learnedParameters = @[LearnedParameter(name: "dt", node: gdt)]
    for c in integrated.learnedCoeffs:
      result.learnedParameters.add c
    for lp in mitems(result.learnedParameters):
      lp.gradientExpr = result.lossExpr.grad lp.node

proc commitAcceptedTrajectory*(graph: TrajectoryGraph,
                               finalGauge: graphGauge.Gauge) =
  finalGauge.reunitGauge
  graph.initialState.gauge.update finalGauge

proc reversibilityCheck*(graph: TrajectoryGraph) =
  ## Integrate (g_f,-p_f); report g_r-g_0, p_r+p_0, and Hamiltonian drift.
  ## Restore the initial leaves after the check.
  let
    h0 = graph.initialState.hamiltonian.eval.sval
    s0 = graph.initialState.gaugeAction.sval
    t0 = graph.initialState.kinetic.sval
  discard graph.finalState.gauge.eval
  discard graph.finalState.momentum.eval
  let
    g0 = graph.initialState.gauge.gaugeSnapshot
    p0 = graph.initialState.momentum.gaugeSnapshot
    g1 = graph.finalState.gauge.gaugeSnapshot
  var p1 = graph.finalState.momentum.gaugeSnapshot   # negated below to start the reverse leg
  threads:
    for mu in 0..<p1.len:
      p1[mu] := -1*p1[mu]
  graph.initialState.gauge.update g1
  graph.initialState.momentum.update p1
  let
    h1 = graph.finalState.hamiltonian.eval.sval
    dH = h1 - h0
    dS = graph.finalState.gaugeAction.sval - s0
    dT = graph.finalState.kinetic.sval - t0
  var
    gd = graph.finalState.gauge.gaugeSnapshot      # reverse-final gauge
    pd = graph.finalState.momentum.gaugeSnapshot   # reverse-final momentum
  graph.initialState.gauge.update g0          # restore the forward initial state
  graph.initialState.momentum.update p0
  # per-link round-trip differences: gd = g_rev - g0, pd = p_rev + p0 (p_rev ≈ -p0).
  # Reduce per-site SIMD norms with simdSum/simdMax then threadRankSum/threadRankMax,
  # the same pattern as gauge/checkSU — no per-lane loop, and cross-thread/rank correct.
  var dg2sum, dp2sum, dg2max, dp2max = 0.0
  threads:
    var sg, sp, mg, mp = 0.0
    for mu in 0..<gd.len:
      gd[mu] -= g0[mu]
      pd[mu] += p0[mu]
    threadBarrier()
    for mu in 0..<gd.len:
      for s in gd[mu]:
        let gn = gd[mu][s].norm2
        let pn = pd[mu][s].norm2
        sg += gn.simdSum
        sp += pn.simdSum
        let gm = gn.simdMax
        if mg < gm: mg = gm
        let pm = pn.simdMax
        if mp < pm: mp = pm
    sg.threadRankSum; sp.threadRankSum
    mg.threadRankMax; mp.threadRankMax
    threadSingle:
      dg2sum = sg; dp2sum = sp; dg2max = mg; dp2max = mp
  let
    nl = float(gd.len * gd[0].l.physVol)
    dgRMS = sqrt(dg2sum / nl)
    dpRMS = sqrt(dp2sum / nl)
    dgMax = sqrt(dg2max)
    dpMax = sqrt(dp2max)
  qexLog "Reversibility: dH: ", dH, "  dS: ", dS, "  dT: ", dT,
    "  dgRMS: ", dgRMS, "  dgMax: ", dgMax, "  dpRMS: ", dpRMS, "  dpMax: ", dpMax
  if abs(dH) > 1e-8 * (abs(h0) + 1.0):
    qexWarn "broken reversibility (|dH|/|H0| > 1e-8): dH: ", dH,
      "  dgRMS: ", dgRMS, "  dpRMS: ", dpRMS

proc runHmc*[R: RNG](graph: TrajectoryGraph;
                     runConfig: RunConfig;
                     randomField: var Field[1, R];
                     randomSerial: var R;
                     measure: proc(traj: int; dH, acc: float; accepted: bool;
                                   forceStats: MdForceStats);
                     onProposal: proc(traj: int; dH, acc: float) = nil) =
  ## Resample, integrate, accept/reject, then call the hooks.
  ## onProposal runs before measure; save the accepted gauge before either hook.
  for traj in 1 .. runConfig.totalTrajs:
    tic("traj")
    echo "Begin traj: ", traj
    graph.resampleMomentum(randomField)
    var h0, s0, t0, h1, s1, t1, dH: float
    block:
      tic("Hamiltonian eval")
      h0 = graph.initialState.hamiltonian.eval.sval
      s0 = graph.initialState.gaugeAction.sval
      t0 = graph.initialState.kinetic.sval
      h1 = graph.finalState.hamiltonian.eval.sval
      s1 = graph.finalState.gaugeAction.sval
      t1 = graph.finalState.kinetic.sval
      dH = h1 - h0
      toc("Hamiltonian eval end")
    echo "Begin H: ", h0, "  S: ", s0, "  T: ", t0
    echo "End H: ", h1, "  S: ", s1, "  T: ", t1
    if runConfig.revCheckFreq > 0 and traj mod runConfig.revCheckFreq == 0:
      graph.reversibilityCheck
    let
      forceStats = graph.mdForces.mdForceStats
      acc = exp(-dH)
      accr = randomSerial.uniform
      forced = traj <= runConfig.trajsForceAcc
      accepted = forced or accr <= acc
    var acceptedGauge: graphGauge.Gauge
    if accepted:
      discard graph.finalState.gauge.eval
      acceptedGauge = graph.finalState.gauge.gaugeSnapshot
    if accepted:
      echo (if forced: "ACCEPT(FORCE)" else: "ACCEPT"),
        ":  dH: ", dH, "  exp(-dH): ", acc, "  r: ", accr
    else:
      echo "REJECT:  dH: ", dH, "  exp(-dH): ", acc, "  r: ", accr
    echo "MD forces: n=", forceStats.count,
      "  fRMS mean/max: ", forceStats.rmsMean, " / ", forceStats.rmsMax,
      "  fMin mean/min: ", forceStats.fminMean, " / ", forceStats.fminMin,
      "  fMax mean/max: ", forceStats.fmaxMean, " / ", forceStats.fmaxMax
    if onProposal != nil: onProposal(traj, dH, acc)
    if accepted: graph.commitAcceptedTrajectory(acceptedGauge)
    measure(traj, dH, acc, accepted, forceStats)
    qexGC "traj done"
    qexLog "traj ", traj, " secs: ", getElapsedTime()
    toc("traj end")
