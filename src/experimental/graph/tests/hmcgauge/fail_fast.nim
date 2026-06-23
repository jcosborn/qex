suite "hmcgauge fail-fast":
  proc validRunConfig(): RunConfig =
    RunConfig(
      dt: 0.025,
      lrmax: 1.0,
      lrmin: 0.0001,
      weightDecay: 0.0,
      trajsThermo: 0,
      trajsTrain: 50,
      trajsTrainlrWarm: 10,
      trajsInfer: 0,
      savefreq: 0,
      gsteps: 4,
      alwaysAccept: false,
      integratorCoeffs: parseIntegratorCoeffs(ik2MN, []))

  proc validIntegratorInputs(): tuple[gc: Gactcoeff, g0: Ggauge, p0: Ggauge, dt: Gscalar] =
    (
      gc: actWilson(scalar.toGvalue(grt, 6.0)),
      g0: grt.toGvalue(g),
      p0: grt.toGvalue(p),
      dt: grt.toGvalue(0.025))

  proc integrateTest(
      inputs: tuple[gc: Gactcoeff, g0: Ggauge, p0: Ggauge, dt: Gscalar],
      coeffs: IntegratorCoeffs,
      steps = 1): IntegrationResult =
    integrateGauge(
      inputs.gc,
      inputs.g0,
      inputs.p0,
      inputs.dt,
      steps,
      coeffs)

  template expectIntegrateError(inputs: untyped,
                                coeffs: untyped,
                                steps: untyped) =
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
    config.trajsInfer = 4

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
      discard integrateGauge(
        inputs.gc,
        inputs.g0,
        inputs.p0,
        inputs.dt,
        0,
        IntegratorCoeffs())
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

    let fourMN3F1GP = parseIntegratorCoeffs(ik4MN3F1GP, [])
    check fourMN3F1GP.kind == ik4MN3F1GP
    check fourMN3F1GP.lambda != 0.0
    check fourMN3F1GP.theta != 0.0
    check fourMN3F1GP.chi != 0.0
    let fourMN3F1GPResult = integrateTest(inputs, fourMN3F1GP)
    check fourMN3F1GPResult.learnedCoeffs.len == 3

    let fourMN5F2GP = parseIntegratorCoeffs(ik4MN5F2GP, [])
    check fourMN5F2GP.kind == ik4MN5F2GP
    check fourMN5F2GP.rho != 0.0
    check fourMN5F2GP.theta != 0.0
    check fourMN5F2GP.vtheta != 0.0
    check fourMN5F2GP.lambda != 0.0
    check fourMN5F2GP.xi != 0.0
    let fourMN5F2GPResult = integrateTest(inputs, fourMN5F2GP)
    check fourMN5F2GPResult.learnedCoeffs.len == 5

  test "integrator variants accept expanded multi-step schedules":
    let inputs = validIntegratorInputs()

    proc checkExpanded(coeffs: IntegratorCoeffs, learnedLen: int) =
      let oneStep = integrateTest(inputs, coeffs)
      let twoStep = integrateTest(inputs, coeffs, 2)
      let chained = integrateGauge(
        inputs.gc,
        oneStep.gauge,
        oneStep.momentum,
        inputs.dt,
        1,
        coeffs)
      check oneStep.learnedCoeffs.len == learnedLen
      check twoStep.learnedCoeffs.len == learnedLen
      norm2(twoStep.gauge - chained.gauge) :< 1e-16
      norm2(twoStep.momentum - chained.momentum) :< 1e-16

    block:
      let coeffs = parseIntegratorCoeffs(ik2MN, [])
      checkExpanded(coeffs, 1)

    block:
      let coeffs = parseIntegratorCoeffs(ik4MN3F1GP, [])
      checkExpanded(coeffs, 3)

    block:
      let coeffs = parseIntegratorCoeffs(ik4MN5F2GP, [])
      checkExpanded(coeffs, 5)

  test "4MN3F1GP rejects partial coefficient tuples":
    let inputs = validIntegratorInputs()

    expectIntegrateError(inputs, parseIntegratorCoeffs(ik4MN3F1GP, [0.0]), 1)

  test "4MN3F1GP accepts explicit full finite tuple":
    let inputs = validIntegratorInputs()
    let result = integrateTest(
      inputs,
      parseIntegratorCoeffs(ik4MN3F1GP, [0.0, 0.25, 0.0]))

    check result.learnedCoeffs.len == 3
    check result.learnedCoeffs[0].name == "lambda"
    check result.learnedCoeffs[0].node.sval == 0.0
    check result.learnedCoeffs[1].name == "theta"
    check result.learnedCoeffs[1].node.sval == 0.25
    check result.learnedCoeffs[2].name == "chi"
    check result.learnedCoeffs[2].node.sval == 0.0

  test "4MN5F2GP rejects partial coefficient tuples":
    let inputs = validIntegratorInputs()

    expectIntegrateError(
      inputs,
      parseIntegratorCoeffs(ik4MN5F2GP, [0.1, 0.2, 0.1, 0.0]),
      1)

  test "4MN5F2GP accepts explicit full finite tuple":
    let inputs = validIntegratorInputs()
    let result = integrateTest(
      inputs,
      parseIntegratorCoeffs(ik4MN5F2GP, [0.1, 0.2, 0.1, 0.0, 0.0]))

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
    let result = integrateGauge(
      gc,
      g0,
      p0,
      dt,
      1,
      parseIntegratorCoeffs(ik4MN3F1GP, []))

    check result.learnedCoeffs.len == 3
    for coeff in result.learnedCoeffs:
      check coeff.node.runtime == grt

  test "trajectory learned parameters keep explicit names":
    block:
      var config = validRunConfig()
      config.integratorCoeffs = parseIntegratorCoeffs(ik2MN, [])
      let graph = buildTrajectoryGraph(
        grt,
        g,
        p,
        actWilson(scalar.toGvalue(grt, 6.0)),
        config)
      check graph.learnedNames == @["dt", "lambda"]

    block:
      var config = validRunConfig()
      config.integratorCoeffs = parseIntegratorCoeffs(ik4MN3F1GP, [])
      let graph = buildTrajectoryGraph(
        grt,
        g,
        p,
        actWilson(scalar.toGvalue(grt, 6.0)),
        config)
      check graph.learnedNames == @["dt", "lambda", "theta", "chi"]

    block:
      var config = validRunConfig()
      config.integratorCoeffs = parseIntegratorCoeffs(ik4MN5F2GP, [])
      let graph = buildTrajectoryGraph(
        grt,
        g,
        p,
        actWilson(scalar.toGvalue(grt, 6.0)),
        config)
      check graph.learnedNames == @[
        "dt", "rho", "theta", "vtheta", "lambda", "xi"]

  test "trajectory resamples graph-owned momentum":
    let graph = buildTrajectoryGraph(
      grt,
      g,
      p,
      actWilson(scalar.toGvalue(grt, 6.0)),
      validRunConfig())
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
    let graph = buildTrajectoryGraph(
      grt,
      g,
      p,
      actWilson(scalar.toGvalue(grt, 6.0)),
      validRunConfig())
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
    let graph = buildTrajectoryGraph(
      grt,
      g,
      p,
      actWilson(scalar.toGvalue(grt, 6.0)),
      validRunConfig())
    discard graph.finalState.gauge.eval
    let acceptedGauge = graph.finalState.gauge.gaugeSnapshot
    let epochBeforeCommit = grt.graphEpochCounter

    graph.commitAcceptedTrajectory(acceptedGauge)

    check grt.graphEpochCounter == epochBeforeCommit + 1
    check graph.initialState.gauge.gaugeSnapshot.len == acceptedGauge.len
    discard graph.lossExpr.eval.sval

  test "training step updates existing learned parameters":
    let config = validRunConfig()
    let graph = buildTrajectoryGraph(
      grt,
      g,
      p,
      actWilson(scalar.toGvalue(grt, 6.0)),
      config)
    var trainer = initTrainingState(graph, config.weightDecay)
    let before = trainer.parameterValues

    trainer.trainStep(config, 1)

    let after = trainer.parameterValues
    check after.len == before.len
    check after.len == graph.learnedParameters.len
    check after != before

  test "training step rejects indexes outside training phase":
    let config = validRunConfig()
    let graph = buildTrajectoryGraph(
      grt,
      g,
      p,
      actWilson(scalar.toGvalue(grt, 6.0)),
      config)
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
      bad.trajsInfer = -1
    expectInvalidConfig:
      bad.savefreq = -1
    expectInvalidConfig:
      bad.gsteps = 0
    expectInvalidConfig:
      bad.lrmin = 2.0
      bad.lrmax = 1.0
    expectInvalidConfig:
      bad.weightDecay = -0.1
