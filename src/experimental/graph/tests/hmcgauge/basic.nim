suite "hmcgauge":
  proc act(gc: Gactcoeff): GaugeAction =
    (proc(g: Ggauge): Gscalar = gaugeAction(gc, g))
  proc validRunConfig(): RunConfig =
    RunConfig(
      dt: 0.025, lrmax: 1.0, lrmin: 0.0001, weightDecay: 0.0,
      trajsThermo: 0, trajs: 50, trajsForceAcc: 0, trajsTrain: 50, trajsTrainlrWarm: 10,
      savefreq: 0, gsteps: 4,
      integratorCoeffs: parseIntegratorCoeffs(ik2MN, []))

  proc validIntegratorInputs(): tuple[gc: Gactcoeff, g0: Ggauge, p0: Ggauge, dt: Gscalar] =
    (gc: actWilson(scalar.toGvalue(grt, 6.0)), g0: grt.toGvalue(g), p0: grt.toGvalue(p), dt: grt.toGvalue(0.025))

  proc integrateTest(inputs: tuple[gc: Gactcoeff, g0: Ggauge, p0: Ggauge, dt: Gscalar]; coeffs: IntegratorCoeffs;
                     steps = 1): IntegrationResult =
    integrateGauge(act(inputs.gc), inputs.g0, inputs.p0, inputs.dt, steps, coeffs)

  template expectIntegrateError(inputs, coeffs, steps: untyped) =
    expect(GraphValueError):
      discard integrateTest(inputs, coeffs, steps)

  template expectInvalidConfig(body: untyped) =
    block:
      var bad {.inject.} = validRunConfig()
      body
      expect(GraphValueError):
        bad.validateRunConfig

  proc learnedNames(graph: TrajectoryGraph): seq[string] =
    for learned in graph.learnedParameters:
      result.add learned.name

  test "RNG names select the supported HMC generators":
    check default(GaugeParams).rng == rkPhilox4x64
    check parseRngKind("Philox4x64") == rkPhilox4x64
    check parseRngKind("PHILOX4X64") == rkPhilox4x64
    check parseRngKind("Threefry4x64") == rkThreefry4x64
    check parseRngKind("THREEFRY4X64") == rkThreefry4x64
    check parseRngKind("MRG32K3A") == rkMrg32k3a
    for kind in RngKind:
      check parseRngKind($kind) == kind
    expect(ValueError):
      discard parseRngKind("unknown")
    expect(ValueError):
      discard parseRngKind("RngMilc6")

  test "RNG dispatch selects the native type":
    for kind in RngKind:
      var selected = rkPhilox4x64
      withRng(kind, R):
        when R is Philox4x64: selected = rkPhilox4x64
        elif R is Threefry4x64: selected = rkThreefry4x64
        elif R is MRG32k3a: selected = rkMrg32k3a
      check selected == kind

  test "trajectoryPhase rejects indexes outside the configured run":
    let config = validRunConfig()

    expect(GraphValueError):
      discard config.trajectoryPhase(-1)
    expect(GraphValueError):
      discard config.trajectoryPhase(0)
    expect(GraphValueError):
      discard config.trajectoryPhase(config.totalTrajs + 1)

  test "trajectoryPhase maps one-based phase boundaries":
    var config = validRunConfig()
    config.trajsThermo = 2
    config.trajsTrain = 3
    config.trajs = 7        # 3 training + 4 inference

    check config.trajectoryPhase(1) == tpThermo
    check config.trajectoryPhase(2) == tpThermo
    check config.trajectoryPhase(3) == tpTrain
    check config.trajectoryPhase(5) == tpTrain
    check config.trajectoryPhase(6) == tpInfer
    check config.trajectoryPhase(9) == tpInfer

  test "optimizer rejects nonpositive step before bias correction":
    var optimizer = initAdamW([1.0], weightDecay = 0.0)
    var parameters = @[1.0]
    let gradients = @[0.5]

    expect(GraphValueError):
      discard optimizer.optimize(parameters, gradients, 0, 0.1)
    expect(GraphValueError):
      discard optimizer.optimize(parameters, gradients, -1, 0.1)

  test "optimizer rejects invalid hyperparameters and inputs":
    expect(GraphValueError):
      discard initAdamW([1.0], stepScale = 0.0)
    expect(GraphValueError):
      discard initAdamW([1.0], beta1 = -0.1)
    expect(GraphValueError):
      discard initAdamW([1.0], beta1 = 1.0)
    expect(GraphValueError):
      discard initAdamW([1.0], beta2 = -0.1)
    expect(GraphValueError):
      discard initAdamW([1.0], beta2 = 1.0)

    var optimizer = initAdamW([1.0], weightDecay = 0.0)
    var parameters = @[1.0]
    expect(GraphValueError):
      discard optimizer.optimize(parameters, @[0.5, 0.25], 1, 0.1)
    expect(GraphValueError):
      discard optimizer.optimize(parameters, @[0.5], 1, -0.1)

  test "warmUpCosDecay handles validated schedule boundaries":
    proc close(a: float, b: float): bool =
      abs(a - b) < 1e-12

    check close(warmUpCosDecay(5, 1, 0, 1.0, 0.2), 0.2)
    check close(warmUpCosDecay(5, 5, 10, 1.0, 0.2), 1.0)
    check close(warmUpCosDecay(0, 0, 10, 1.0, 0.2), 1.0)
    check close(warmUpCosDecay(1, 0, 10, 1.0, 0.2), 1.0)
    check warmUpCosDecay(2, 0, 10, 1.0, 0.2) < 1.0
    check close(warmUpCosDecay(10, 0, 10, 1.0, 0.2), 0.2)

  test "integrateGauge rejects nonpositive step count before gauge ops":
    let inputs = validIntegratorInputs()

    expectIntegrateError(inputs, parseIntegratorCoeffs(ik2MN, []), 0)
    expectIntegrateError(inputs, parseIntegratorCoeffs(ik2MN, []), -1)

  test "integrateGauge rejects step count before building a spec":
    let inputs = validIntegratorInputs()
    var failed = false

    try:
      discard integrateGauge(act(inputs.gc), inputs.g0, inputs.p0, inputs.dt, 0, IntegratorCoeffs())
    except GraphValueError as e:
      failed = true
      check e.msg.contains("integrator step count")

    check failed

  test "integrator variants keep kind and learned coefficient counts":
    let inputs = validIntegratorInputs()

    let twoMN = parseIntegratorCoeffs(ik2MN, [])
    check twoMN.kind == ik2MN
    check twoMN.lambda != 0.0
    let twoMNResult: IntegrationResult = integrateTest(inputs, twoMN)
    check twoMNResult.gauge != nil
    check twoMNResult.momentum != nil
    check twoMNResult.gauge.runtime == inputs.g0.runtime
    check twoMNResult.momentum.runtime == inputs.p0.runtime
    check twoMNResult.learnedCoeffs.len == 1
    check twoMNResult.forces.len == 2

    let fourMN3F1GP = parseIntegratorCoeffs(ik4MN3F1GP, [])
    check fourMN3F1GP.kind == ik4MN3F1GP
    check fourMN3F1GP.lambda != 0.0
    check fourMN3F1GP.theta != 0.0
    check fourMN3F1GP.chi != 0.0
    let fourMN3F1GPResult = integrateTest(inputs, fourMN3F1GP)
    check fourMN3F1GPResult.learnedCoeffs.len == 3
    check fourMN3F1GPResult.forces.len == 4

    let fourMN5F2GP = parseIntegratorCoeffs(ik4MN5F2GP, [])
    check fourMN5F2GP.kind == ik4MN5F2GP
    check fourMN5F2GP.rho != 0.0
    check fourMN5F2GP.theta != 0.0
    check fourMN5F2GP.vtheta != 0.0
    check fourMN5F2GP.lambda != 0.0
    check fourMN5F2GP.xi != 0.0
    let fourMN5F2GPResult = integrateTest(inputs, fourMN5F2GP)
    check fourMN5F2GPResult.learnedCoeffs.len == 5
    check fourMN5F2GPResult.forces.len == 7

  test "integrator variants accept expanded multi-step schedules":
    let inputs = validIntegratorInputs()

    proc checkExpanded(coeffs: IntegratorCoeffs; learnedLen, forcesPerStep: int) =
      let oneStep = integrateTest(inputs, coeffs)
      let twoStep = integrateTest(inputs, coeffs, 2)
      let chained = integrateGauge(act(inputs.gc), oneStep.gauge, oneStep.momentum, inputs.dt, 1, coeffs)
      check oneStep.learnedCoeffs.len == learnedLen
      check twoStep.learnedCoeffs.len == learnedLen
      check oneStep.forces.len == forcesPerStep
      check twoStep.forces.len == 2 * forcesPerStep
      norm2(twoStep.gauge - chained.gauge) :< 1e-16
      norm2(twoStep.momentum - chained.momentum) :< 1e-16

    block:
      let coeffs = parseIntegratorCoeffs(ik2MN, [])
      checkExpanded(coeffs, 1, 2)

    block:
      let coeffs = parseIntegratorCoeffs(ik4MN3F1GP, [])
      checkExpanded(coeffs, 3, 4)

    block:
      let coeffs = parseIntegratorCoeffs(ik4MN5F2GP, [])
      checkExpanded(coeffs, 5, 7)

  test "integrator direct force matches action differentiation":
    let
      inputs = validIntegratorInputs()
      coeffs = parseIntegratorCoeffs(ik2MN, [])
      automatic = integrateTest(inputs, coeffs)
    proc force(x: Ggauge): Ggauge = gaugeForce(inputs.gc, x)
    let direct = integrateGauge(act(inputs.gc), inputs.g0, inputs.p0, inputs.dt, 1, coeffs, force)
    norm2(direct.gauge - automatic.gauge) :< 1e-16
    norm2(direct.momentum - automatic.momentum) :< 1e-16
    check direct.forces.len == automatic.forces.len
    for i in 0..<direct.forces.len:
      norm2(direct.forces[i] - automatic.forces[i]) :< 1e-16

  test "integrator reuses evaluated nodes for MD force statistics":
    let inputs = validIntegratorInputs()
    let result = integrateTest(inputs, parseIntegratorCoeffs(ik2MN, []), 2)
    discard result.momentum.eval
    var runs = newSeq[int](result.forces.len)
    for i, force in result.forces:
      runs[i] = force.runCount
      check runs[i] > 0
    let stats = result.forces.mdForceStats

    check stats.count == 4
    check stats.rmsMean >= 0.0
    check stats.rmsMax >= stats.rmsMean
    check stats.fminMin > 0.0
    check stats.fminMean >= stats.fminMin and stats.fminMean <= stats.rmsMean
    check stats.fmaxMean >= stats.rmsMean
    check stats.fmaxMax >= stats.fmaxMean
    for i, force in result.forces:
      check force.runCount == runs[i]

  test "force extrema include exact zero magnitudes":
    let force = grt.toGvalue(zeroGaugeLike(g))
    let stats = force.forceRmsMinMax(float(g.len * lo.physVol))

    check stats.rms == 0.0
    check stats.fmin == 0.0
    check stats.fmax == 0.0

  test "4MN3F1GP rejects partial coefficient tuples":
    let inputs = validIntegratorInputs()

    expectIntegrateError(inputs, parseIntegratorCoeffs(ik4MN3F1GP, [0.0]), 1)

  test "4MN3F1GP accepts explicit full finite tuple":
    let inputs = validIntegratorInputs()
    let result = integrateTest(inputs, parseIntegratorCoeffs(ik4MN3F1GP, [0.0, 0.25, 0.0]))

    check result.learnedCoeffs.len == 3
    check result.learnedCoeffs[0].name == "lambda"
    check result.learnedCoeffs[0].node.sval == 0.0
    check result.learnedCoeffs[1].name == "theta"
    check result.learnedCoeffs[1].node.sval == 0.25
    check result.learnedCoeffs[2].name == "chi"
    check result.learnedCoeffs[2].node.sval == 0.0

  test "4MN5F2GP rejects partial coefficient tuples":
    let inputs = validIntegratorInputs()

    expectIntegrateError(inputs, parseIntegratorCoeffs(ik4MN5F2GP, [0.1, 0.2, 0.1, 0.0]), 1)

  test "4MN5F2GP accepts explicit full finite tuple":
    let inputs = validIntegratorInputs()
    let result = integrateTest(inputs, parseIntegratorCoeffs(ik4MN5F2GP, [0.1, 0.2, 0.1, 0.0, 0.0]))

    check result.learnedCoeffs.len == 5
    check result.learnedCoeffs[0].name == "rho"
    check result.learnedCoeffs[0].node.sval == 0.1
    check result.learnedCoeffs[1].name == "theta"
    check result.learnedCoeffs[1].node.sval == 0.2
    check result.learnedCoeffs[2].name == "vtheta"
    check result.learnedCoeffs[2].node.sval == 0.1
    check result.learnedCoeffs[3].name == "lambda"
    check result.learnedCoeffs[3].node.sval == 0.0
    check result.learnedCoeffs[4].name == "xi"
    check result.learnedCoeffs[4].node.sval == 0.0

  test "integrator learned coefficients inherit dt runtime":
    let grt = initGraphRuntime()
    let gc = actWilson(scalar.toGvalue(grt, 6.0))
    let g0 = gauge.toGvalue(grt, g)
    let p0 = gauge.toGvalue(grt, p)
    let dt = scalar.toGvalue(grt, 0.025)
    let result = integrateGauge(act(gc), g0, p0, dt, 1, parseIntegratorCoeffs(ik4MN3F1GP, []))

    check result.learnedCoeffs.len == 3
    for coeff in result.learnedCoeffs:
      check coeff.node.runtime == grt

  test "trajectory learned parameters keep explicit names":
    block:
      var config = validRunConfig()
      config.integratorCoeffs = parseIntegratorCoeffs(ik2MN, [])
      let graph = buildTrajectoryGraph(grt, g, p, act(actWilson(scalar.toGvalue(grt, 6.0))), config)
      check graph.learnedNames == @["dt", "lambda"]
      check graph.mdForces.len == 2 * config.gsteps

    block:
      var config = validRunConfig()
      config.integratorCoeffs = parseIntegratorCoeffs(ik4MN3F1GP, [])
      let graph = buildTrajectoryGraph(grt, g, p, act(actWilson(scalar.toGvalue(grt, 6.0))), config)
      check graph.learnedNames == @["dt", "lambda", "theta", "chi"]
      check graph.mdForces.len == 4 * config.gsteps

    block:
      var config = validRunConfig()
      config.integratorCoeffs = parseIntegratorCoeffs(ik4MN5F2GP, [])
      let graph = buildTrajectoryGraph(grt, g, p, act(actWilson(scalar.toGvalue(grt, 6.0))), config)
      check graph.learnedNames == @["dt", "rho", "theta", "vtheta", "lambda", "xi"]
      check graph.mdForces.len == 7 * config.gsteps

  test "trajectory resamples graph-owned momentum":
    let graph = buildTrajectoryGraph(grt, g, p, act(actWilson(scalar.toGvalue(grt, 6.0))), validRunConfig())
    var before = 0.0
    for mu in graph.initialState.momentum.gaugeSnapshot:
      before += mu.norm2

    graph.resampleMomentum(r)
    var after = 0.0
    for mu in graph.initialState.momentum.gaugeSnapshot:
      after += mu.norm2

    check graph.initialState.gauge.gaugeSnapshot.len == g.len
    check graph.initialState.momentum.gaugeSnapshot.len == p.len
    discard graph.finalState.gauge.eval
    check graph.finalState.gauge.gaugeSnapshot.len == g.len
    check after != before
    discard graph.lossExpr.eval.sval

  test "accepted trajectory commit uses a pre-training final gauge snapshot":
    let graph = buildTrajectoryGraph(grt, g, p, act(actWilson(scalar.toGvalue(grt, 6.0))), validRunConfig())
    discard graph.finalState.gauge.eval
    let acceptedGauge = graph.finalState.gauge.gaugeSnapshot
    let expectedGauge = grt.toGvalue(acceptedGauge)
    mutateGauge(expectedGauge, storage):
      storage.reunitGauge

    for learned in graph.learnedParameters:
      learned.node.update learned.node.sval + 0.001

    graph.commitAcceptedTrajectory(acceptedGauge)

    let currentGauge = graph.initialState.gauge.gaugeSnapshot
    check currentGauge.len == acceptedGauge.len
    norm2(grt.toGvalue(currentGauge) - expectedGauge) :< 1e-26
    discard graph.lossExpr.eval.sval

  test "accepted trajectory commit marks freshness":
    let graph = buildTrajectoryGraph(grt, g, p, act(actWilson(scalar.toGvalue(grt, 6.0))), validRunConfig())
    discard graph.finalState.gauge.eval
    let acceptedGauge = graph.finalState.gauge.gaugeSnapshot
    let epochBeforeCommit = grt.graphEpochCounter

    graph.commitAcceptedTrajectory(acceptedGauge)

    check grt.graphEpochCounter == epochBeforeCommit + 1
    check graph.initialState.gauge.gaugeSnapshot.len == acceptedGauge.len
    discard graph.lossExpr.eval.sval

  test "training step updates existing learned parameters":
    let config = validRunConfig()
    let graph = buildTrajectoryGraph(grt, g, p, act(actWilson(scalar.toGvalue(grt, 6.0))), config)
    var trainer = initTrainingState(graph, config.weightDecay)
    let before = trainer.parameterValues

    trainer.trainStep(config, 1)

    let after = trainer.parameterValues
    check after.len == before.len
    check after.len == graph.learnedParameters.len
    check after != before

  test "training step rejects indexes outside training phase":
    let config = validRunConfig()
    let graph = buildTrajectoryGraph(grt, g, p, act(actWilson(scalar.toGvalue(grt, 6.0))), config)
    var trainer = initTrainingState(graph, config.weightDecay)

    expect(GraphValueError):
      trainer.trainStep(config, 0)
    expect(GraphValueError):
      trainer.trainStep(config, config.trajsTrain + 1)

  test "nonempty missing gauge file fails instead of silently uniting":
    var localGauge = lo.newgauge

    expect(GraphValueError):
      localGauge.loadOrInitGauge("/definitely/missing/qex-graph-test.lime")

  test "run config validation rejects invalid ranges":
    expectInvalidConfig:
      bad.dt = 0.0
    expectInvalidConfig:
      bad.lrmin = -0.1
    expectInvalidConfig:
      bad.lrmax = -0.1
    expectInvalidConfig:
      bad.trajsThermo = -1
    expectInvalidConfig:
      bad.trajsTrain = -1
    expectInvalidConfig:
      bad.trajsTrainlrWarm = -1
    expectInvalidConfig:
      bad.trajsTrainlrWarm = bad.trajsTrain + 1
    expectInvalidConfig:
      bad.trajs = -1
    expectInvalidConfig:
      bad.trajsForceAcc = -1
    expectInvalidConfig:
      bad.trajsTrain = bad.trajs + 1
    expectInvalidConfig:
      bad.savefreq = -1
    expectInvalidConfig:
      bad.gsteps = 0
    expectInvalidConfig:
      bad.lrmin = 2.0
      bad.lrmax = 1.0
    expectInvalidConfig:
      bad.weightDecay = -0.1
