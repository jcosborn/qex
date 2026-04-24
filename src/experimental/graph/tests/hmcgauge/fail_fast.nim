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
      integratorCoeffs: parseIntegratorCoeffs(ik2MN, []),
      alwaysAccept: false)

  proc validIntegratorInputs(): tuple[gc: Gactcoeff, g0: Ggauge, p0: Ggauge, dt: Gscalar] =
    (
      gc: grt.actWilson(6.0),
      g0: grt.toGvalue(g),
      p0: grt.toGvalue(p),
      dt: grt.toGvalue(0.025))

  test "optimizer rejects nonpositive step before bias correction":
    var optimizer = initAdamW([1.0], weightDecay = 0.0)
    var parameters = @[1.0]
    let gradients = @[0.5]

    expect(GraphValueError):
      discard optimizer.optimize(parameters, gradients, 0, 0.1)
    expect(GraphValueError):
      discard optimizer.optimize(parameters, gradients, -1, 0.1)

  test "integrateGauge rejects nonpositive step count before gauge ops":
    let inputs = validIntegratorInputs()

    expect(GraphValueError):
      discard integrateGauge(
        inputs.gc,
        inputs.g0,
        inputs.p0,
        inputs.dt,
        0,
        parseIntegratorCoeffs(ik2MN, []))
    expect(GraphValueError):
      discard integrateGauge(
        inputs.gc,
        inputs.g0,
        inputs.p0,
        inputs.dt,
        -1,
        parseIntegratorCoeffs(ik2MN, []))

  test "4MN3F1GP rejects partial coefficient tuples":
    let inputs = validIntegratorInputs()

    expect(GraphValueError):
      discard integrateGauge(
        inputs.gc,
        inputs.g0,
        inputs.p0,
        inputs.dt,
        1,
        parseIntegratorCoeffs(ik4MN3F1GP, [0.0]))

  test "4MN3F1GP accepts explicit full finite tuple":
    let inputs = validIntegratorInputs()
    let (g1, p1, learnedCoeffs) = integrateGauge(
      inputs.gc,
      inputs.g0,
      inputs.p0,
      inputs.dt,
      1,
      parseIntegratorCoeffs(ik4MN3F1GP, [0.0, 0.25, 0.0]))

    check g1 != nil
    check p1 != nil
    check learnedCoeffs.len == 3

  test "4MN5F2GP rejects partial coefficient tuples":
    let inputs = validIntegratorInputs()

    expect(GraphValueError):
      discard integrateGauge(
        inputs.gc,
        inputs.g0,
        inputs.p0,
        inputs.dt,
        1,
        parseIntegratorCoeffs(ik4MN5F2GP, [0.1, 0.2, 0.1, 0.0]))

  test "4MN5F2GP accepts explicit full finite tuple":
    let inputs = validIntegratorInputs()
    let (g1, p1, learnedCoeffs) = integrateGauge(
      inputs.gc,
      inputs.g0,
      inputs.p0,
      inputs.dt,
      1,
      parseIntegratorCoeffs(ik4MN5F2GP, [0.1, 0.2, 0.1, 0.0, 0.0]))

    check g1 != nil
    check p1 != nil
    check learnedCoeffs.len == 5

  test "integrator learned coefficients inherit dt runtime":
    let grt = initGraphRuntime()
    let gc = actWilson(scalar.toGvalue(grt, 6.0))
    let g0 = gauge.toGvalue(grt, g)
    let p0 = gauge.toGvalue(grt, p)
    let dt = scalar.toGvalue(grt, 0.025)
    let (_, _, learnedCoeffs) = integrateGauge(
      gc,
      g0,
      p0,
      dt,
      1,
      parseIntegratorCoeffs(ik4MN3F1GP, []))

    check learnedCoeffs.len == 3
    for coeff in learnedCoeffs:
      check coeff.runtime == grt

  test "trajectory resamples graph-owned momentum":
    let graph = buildTrajectoryGraph(grt, g, p, grt.actWilson(6.0), validRunConfig())
    var before = 0.0
    for mu in graph.currentMomentum:
      before += mu.norm2

    graph.resampleMomentum(r)
    var after = 0.0
    for mu in graph.currentMomentum:
      after += mu.norm2

    check graph.currentGauge.len == g.len
    check after != before
    discard graph.lossValue

  test "accepted trajectory commit uses a pre-training final gauge snapshot":
    let graph = buildTrajectoryGraph(grt, g, p, grt.actWilson(6.0), validRunConfig())
    let acceptedGauge = graph.finalGaugeSnapshot
    let expectedGauge = acceptedGauge.gaugeNodeLike
    expectedGauge.valCopy acceptedGauge
    expectedGauge.getgauge.reunitGauge
    expectedGauge.updated

    for learned in graph.learnedParameters:
      learned.node.update learned.node.getfloat + 0.001

    graph.commitAcceptedTrajectory(acceptedGauge)

    check graph.currentGauge.len == acceptedGauge.getgauge.len
    norm2(grt.toGvalue(graph.currentGauge) - expectedGauge) :< 1e-26
    discard graph.lossValue

  test "run config validation rejects invalid ranges":
    var bad = validRunConfig()
    bad.trajsTrain = -1
    expect(GraphValueError):
      bad.validateRunConfig

    bad = validRunConfig()
    bad.savefreq = -1
    expect(GraphValueError):
      bad.validateRunConfig

    bad = validRunConfig()
    bad.gsteps = 0
    expect(GraphValueError):
      bad.validateRunConfig

    bad = validRunConfig()
    bad.lrmin = 2.0
    bad.lrmax = 1.0
    expect(GraphValueError):
      bad.validateRunConfig
