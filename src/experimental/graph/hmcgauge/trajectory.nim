import qex
import ../[core, scalar, multi, gauge, plan]
import ../gauge/types as graphGauge
import config, integrator
export GaugeAction, GaugeForceFn, MdForceStats, IntegrationEventKind, IntegrationEvent

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
    trace*: seq[IntegrationEvent]
      ## Integration events when built with trace = true; empty otherwise.
  Proposal* = object
    ## Gauge values are borrowed through onProposal; save an owned gaugeSnapshot.
    ## loss is defined when lossExpr exists; gradients exist only during training.
    gauge*, view*: Ggauge
    dH*, acc*, loss*: float
    gradients*: seq[float]

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
                           force: GaugeForceFn = nil,
                           parameters: openArray[LearnedParameter] = [],
                           trace = false): TrajectoryGraph =
  ## Build HMC from S(g); the same S defines the force and both Hamiltonians.
  ## If not buildTraining, omit the loss and parameter-gradient graph.
  ## parameters adds exposed action/flow parameters to dt and integrator training.
  ## trace keeps every integration event as graph nodes for tests and debugging.
  let gdt = toGvalue(grt, config.dt)
  result.initialState = buildTrajectoryState(action, toGvalue(grt, g), toGvalue(grt, p))
  let integrated = integrateGauge(
    action,
    result.initialState.gauge,
    result.initialState.momentum,
    gdt,
    config.gsteps,
    config.integratorCoeffs,
    force,
    trace)
  result.mdForces = integrated.forces
  result.trace = integrated.trace
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
    for p in parameters:
      discard sharedGraphRuntime([Gvalue(result.lossExpr), Gvalue(p.node)], "trajectory parameter")
      for q in result.learnedParameters:
        if p.node.nodeKey == q.node.nodeKey or p.name == q.name:
          raiseValueError("trajectory learned parameters require distinct nodes and names: " & p.name)
      result.learnedParameters.add p
    for lp in mitems(result.learnedParameters):
      lp.gradientExpr = result.lossExpr.grad lp.node

proc commitAcceptedTrajectory*(graph: TrajectoryGraph,
                               finalGauge: graphGauge.Gauge) =
  finalGauge.reunitGauge
  graph.initialState.gauge.update finalGauge

proc reversibilityCheck(graph: TrajectoryGraph;
                        finalGauge, finalMomentum: Ggauge;
                        h0, s0, t0: float;
                        planned = true) =
  ## Integrate (g_f,-p_f); report g_r-g_0, p_r+p_0, and Hamiltonian drift.
  ## Restore the initial leaves after the check.
  let
    g0 = graph.initialState.gauge.gaugeSnapshot
    p0 = graph.initialState.momentum.gaugeSnapshot
    g1 = finalGauge.gaugeSnapshot
  var p1 = finalMomentum.gaugeSnapshot   # negated below to start the reverse leg
  threads:
    for mu in 0..<p1.len:
      p1[mu] := -1*p1[mu]
  var
    dH, dS, dT: float
    gd, pd: graphGauge.Gauge
    reverse: GraphPlan
  try:
    if planned:
      reverse = plan(graph.finalState.gauge, graph.finalState.momentum,
        graph.finalState.hamiltonian, graph.finalState.gaugeAction, graph.finalState.kinetic)
    graph.initialState.gauge.update g1
    graph.initialState.momentum.update p1
    if reverse != nil:
      let values = reverse.eval
      gd = Ggauge(values[0]).gaugeSnapshot
      pd = Ggauge(values[1]).gaugeSnapshot
      dH = Gscalar(values[2]).sval - h0
      dS = Gscalar(values[3]).sval - s0
      dT = Gscalar(values[4]).sval - t0
    else:
      # Explicit retained reference path for standalone checks.
      dH = graph.finalState.hamiltonian.eval.sval - h0
      dS = graph.finalState.gaugeAction.sval - s0
      dT = graph.finalState.kinetic.sval - t0
      gd = graph.finalState.gauge.gaugeSnapshot
      pd = graph.finalState.momentum.gaugeSnapshot
  finally:
    graph.initialState.gauge.update g0
    graph.initialState.momentum.update p0
    if reverse != nil: reverse.clear
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

proc reversibilityCheck*(graph: TrajectoryGraph) =
  ## Standalone retained reference; runHmc supplies its planned forward values.
  let
    h0 = graph.initialState.hamiltonian.eval.sval
    s0 = graph.initialState.gaugeAction.sval
    t0 = graph.initialState.kinetic.sval
  discard graph.finalState.gauge.eval
  discard graph.finalState.momentum.eval
  graph.reversibilityCheck(graph.finalState.gauge, graph.finalState.momentum,
    h0, s0, t0, planned = false)

proc runHmc*[R: RNG](graph: TrajectoryGraph;
                     runConfig: RunConfig;
                     randomField: var Field[1, R];
                     randomSerial: var R;
                     measure: proc(traj: int; dH, acc: float; accepted: bool;
                                   forceStats: MdForceStats);
                     onProposal: proc(traj: int; proposal: Proposal) = nil;
                     proposalView: Ggauge = nil) =
  ## One shared proposal evaluation; callbacks consume its values before commit.
  if runConfig.trajsTrain > 0 and graph.lossExpr == nil:
    raiseValueError("training trajectories require buildTraining=true")
  var roots = @[Gvalue(graph.initialState.hamiltonian),
    Gvalue(graph.initialState.gaugeAction), Gvalue(graph.initialState.kinetic)]
  for force in graph.mdForces:
    roots.add force.forceStats
  roots.add graph.finalState.gauge
  roots.add graph.finalState.hamiltonian
  roots.add graph.finalState.gaugeAction
  roots.add graph.finalState.kinetic
  let reverse = runConfig.revCheckFreq > 0
  if reverse: roots.add graph.finalState.momentum
  if proposalView != nil: roots.add proposalView
  if graph.lossExpr != nil: roots.add graph.lossExpr
  let common = roots.len
  var
    active: GraphPlan
    training = false
  defer:
    if active != nil: active.clear
  for traj in 1 .. runConfig.totalTrajs:
    let train = runConfig.trajectoryPhase(traj) == tpTrain
    if active == nil or train != training:
      if active != nil: active.clear
      roots.setLen(common)
      if train:
        for lp in graph.learnedParameters: roots.add lp.gradientExpr
      active = plan(roots)
      training = train
    tic("traj")
    echo "Begin traj: ", traj
    graph.resampleMomentum(randomField)
    let
      values = block:
        tic("proposal eval")
        let evaluated = active.eval
        toc("proposal eval end")
        evaluated
      h0 = Gscalar(values[0]).sval
      s0 = Gscalar(values[1]).sval
      t0 = Gscalar(values[2]).sval
    var
      i = 3
      forces = newSeq[tuple[rms, fmin, fmax: float]](graph.mdForces.len)
    for f in mitems(forces):
      let diag = Gmulti(values[i])
      f = (rms: Gscalar(diag.storedSlot(0)).sval,
        fmin: Gscalar(diag.storedSlot(1)).sval,
        fmax: Gscalar(diag.storedSlot(2)).sval)
      inc i
    var proposal = Proposal(gauge: Ggauge(values[i]))
    let
      h1 = Gscalar(values[i+1]).sval
      s1 = Gscalar(values[i+2]).sval
      t1 = Gscalar(values[i+3]).sval
    i += 4
    var momentum: Ggauge
    if reverse:
      momentum = Ggauge(values[i])
      inc i
    proposal.view = proposal.gauge
    if proposalView != nil:
      proposal.view = Ggauge(values[i])
      inc i
    if graph.lossExpr != nil:
      proposal.loss = Gscalar(values[i]).sval
      inc i
    if train:
      proposal.gradients = newSeq[float](graph.learnedParameters.len)
      for gradient in mitems(proposal.gradients):
        gradient = Gscalar(values[i]).sval
        inc i
    proposal.dH = h1 - h0
    proposal.acc = exp(-proposal.dH)
    let
      forceStats = forces.mdForceStats
      dH = proposal.dH
      acc = proposal.acc
    echo "Begin H: ", h0, "  S: ", s0, "  T: ", t0
    echo "End H: ", h1, "  S: ", s1, "  T: ", t1
    let
      accr = randomSerial.uniform
      forced = traj <= runConfig.trajsForceAcc
      accepted = forced or accr <= acc
    var acceptedGauge: graphGauge.Gauge
    if accepted:
      acceptedGauge = proposal.gauge.gaugeSnapshot
    if reverse and traj mod runConfig.revCheckFreq == 0:
      graph.reversibilityCheck(proposal.gauge, momentum, h0, s0, t0)
    if accepted:
      echo (if forced: "ACCEPT(FORCE)" else: "ACCEPT"),
        ":  dH: ", dH, "  exp(-dH): ", acc, "  r: ", accr
    else:
      echo "REJECT:  dH: ", dH, "  exp(-dH): ", acc, "  r: ", accr
    echo "MD forces: n=", forceStats.count,
      "  fRMS mean/max: ", forceStats.rmsMean, " / ", forceStats.rmsMax,
      "  fMin mean/min: ", forceStats.fminMean, " / ", forceStats.fminMin,
      "  fMax mean/max: ", forceStats.fmaxMean, " / ", forceStats.fmaxMax
    if onProposal != nil: onProposal(traj, proposal)
    if accepted: graph.commitAcceptedTrajectory(acceptedGauge)
    measure(traj, dH, acc, accepted, forceStats)
    qexGC "traj done"
    qexLog "traj ", traj, " secs: ", getElapsedTime()
    toc("traj end")
