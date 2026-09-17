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
                     steps = 1; trace = false): IntegrationResult =
    integrateGauge(act(inputs.gc), inputs.g0, inputs.p0, inputs.dt, steps, coeffs, trace = trace)

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

  proc forceValues(x: Gmulti): tuple[rms, fmin, fmax: float] =
    (rms: Gscalar(x.storedSlot(0)).sval,
      fmin: Gscalar(x.storedSlot(1)).sval,
      fmax: Gscalar(x.storedSlot(2)).sval)

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
    expectIntegrateError(inputs, parseIntegratorCoeffs(ik2MNp, []), 0)

  test "both 2MN orderings share the minimal-norm lambda and have explicit names":
    check parseIntegratorKind("2MNp") == ik2MNp
    check parseIntegratorCoeffs(ik2MNp, []).lambda == 0.1931833275037836
    check parseIntegratorCoeffs(ik2MN, []).lambda == 0.1931833275037836
    check parseIntegratorCoeffs(ik2MNp, [0.21]).lambda == 0.21
    expect(GraphValueError):
      discard parseIntegratorCoeffs(ik2MNp, [0.1,0.2])

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
    check twoMNResult.trace.len == 0
    check integrateTest(inputs, twoMN, trace = true).trace.len == 7

    let twoMNp = integrateTest(inputs,parseIntegratorCoeffs(ik2MNp, []),trace = true)
    check twoMNp.learnedCoeffs.len == 1
    check twoMNp.forces.len == 3
    check twoMNp.trace.len == 8
    check twoMNp.trace[0].kind == ieForce
    check twoMNp.trace[0].gauge.nodeKey == inputs.g0.nodeKey
    check twoMNp.trace[1].kind == ieKick
    check twoMNp.trace[2].kind == ieDrift
    check twoMNp.trace[2].momentum.nodeKey == twoMNp.trace[1].momentum.nodeKey

    let fourMN3F1GP = parseIntegratorCoeffs(ik4MN3F1GP, [])
    check fourMN3F1GP.kind == ik4MN3F1GP
    check fourMN3F1GP.lambda != 0.0
    check fourMN3F1GP.theta != 0.0
    check fourMN3F1GP.chi != 0.0
    let fourMN3F1GPResult = integrateTest(inputs, fourMN3F1GP)
    check fourMN3F1GPResult.learnedCoeffs.len == 3
    check fourMN3F1GPResult.forces.len == 4
    check fourMN3F1GPResult.trace.len == 0
    check integrateTest(inputs, fourMN3F1GP, trace = true).trace.len == 12

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
    check integrateTest(inputs, fourMN5F2GP, trace = true).trace.len == 20

  test "integrator variants accept expanded multi-step schedules":
    let inputs = validIntegratorInputs()

    proc checkExpanded(coeffs: IntegratorCoeffs; learnedLen, forces1, forces2: int) =
      let oneStep = integrateTest(inputs, coeffs)
      let twoStep = integrateTest(inputs, coeffs, 2)
      let chained = integrateGauge(act(inputs.gc), oneStep.gauge, oneStep.momentum, inputs.dt, 1, coeffs)
      check oneStep.learnedCoeffs.len == learnedLen
      check twoStep.learnedCoeffs.len == learnedLen
      check oneStep.forces.len == forces1
      check twoStep.forces.len == forces2
      norm2(twoStep.gauge - chained.gauge) :< 1e-16
      norm2(twoStep.momentum - chained.momentum) :< 1e-16

    block:
      let coeffs = parseIntegratorCoeffs(ik2MN, [])
      checkExpanded(coeffs, 1, 2, 4)

    block:
      let coeffs = parseIntegratorCoeffs(ik2MNp, [])
      checkExpanded(coeffs, 1, 3, 5)

    block:
      let coeffs = parseIntegratorCoeffs(ik4MN3F1GP, [])
      checkExpanded(coeffs, 3, 4, 8)

    block:
      let coeffs = parseIntegratorCoeffs(ik4MN5F2GP, [])
      checkExpanded(coeffs, 5, 7, 14)

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

    let mp = parseIntegratorCoeffs(ik2MNp, [])
    let ma = integrateTest(inputs,mp)
    let md = integrateGauge(act(inputs.gc),inputs.g0,inputs.p0,inputs.dt,1,mp,force)
    norm2(md.gauge-ma.gauge) :< 1e-16
    norm2(md.momentum-ma.momentum) :< 1e-16

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

  test "copied MD force triples aggregate means and extrema":
    let stats = [(rms: 2.0, fmin: 1.0, fmax: 3.0),
      (rms: 4.0, fmin: 0.0, fmax: 7.0)].mdForceStats
    check stats.count == 2
    check stats.rmsMean == 3.0
    check stats.rmsMax == 4.0
    check stats.fminMean == 0.5
    check stats.fminMin == 0.0
    check stats.fmaxMean == 5.0
    check stats.fmaxMax == 7.0
    check [(rms: 0.0, fmin: 0.0, fmax: 0.0)].mdForceStats == MdForceStats(count: 1)

  test "MD force statistics reject empty force and scalar lists":
    expect(GraphValueError):
      discard newSeq[Ggauge]().mdForceStats
    expect(GraphValueError):
      discard newSeq[tuple[rms, fmin, fmax: float]]().mdForceStats

  test "planned force diagnostics match direct reductions and isolate source caches":
    let
      lo = lat.newLayout
      g = lo.newgauge
      p = lo.newgauge
      dof = float(g.len * lo.physVol)
    var rng = lo.newRNGField(Philox4x64, 921260914u64)
    threads:
      g.random rng
      p.randomTAH rng
    for kind in IntegratorKind:
      let
        grt = initGraphRuntime()
        gc = actWilson(grt.toGvalue(6.0))
        dt = grt.toGvalue(0.025)
        traj = integrateGauge(act(gc), grt.toGvalue(g), grt.toGvalue(p), dt,
          2, parseIntegratorCoeffs(kind, []))
      var
        stats: seq[Gmulti]
        roots: seq[Gvalue]
      for force in traj.forces:
        let diag = force.forceStats
        stats.add diag
        roots.add diag
      roots.add traj.gauge
      let planned = plan(roots)
      var
        runs = newSeq[int](stats.len)
        forceRuns = newSeq[int](stats.len)
        previous = newSeq[tuple[rms, fmin, fmax: float]](stats.len)
      for pass in 0..1:
        discard planned.eval
        var values = newSeq[tuple[rms, fmin, fmax: float]](stats.len)
        for i, diag in stats:
          check diag.runCount == runs[i]
          check traj.forces[i].runCount == forceRuns[i]
          if pass == 0:
            check not diag.valueReady
            check not traj.forces[i].valueReady
            check not traj.forces[i].hasStorage
          else:
            check diag.forceValues == previous[i]
          let value = Gmulti(planned[i])
          check value.len == 3
          check value.bufferBytes == 0
          for j in 0..<3:
            check value.storedSlot(j) of Gscalar
          values[i] = value.forceValues
        if pass > 0:
          check values != previous
        let forwards = planned.stats.forwards
        discard planned.eval
        check planned.stats.forwards == forwards

        for i, diag in stats:
          let want = traj.forces[i].forceRmsMinMax(dof)
          discard diag.eval
          let direct = diag.forceValues
          for (got, expected) in [
              (values[i].rms, want.rms), (values[i].fmin, want.fmin), (values[i].fmax, want.fmax),
              (direct.rms, want.rms), (direct.fmin, want.fmin), (direct.fmax, want.fmax)]:
            check abs(got - expected) < 1e-12 * (1.0 + abs(expected))
          runs[i] = diag.runCount
          forceRuns[i] = traj.forces[i].runCount
          previous[i] = direct
        let
          got = values.mdForceStats
          want = traj.forces.mdForceStats
        check got.count == want.count
        for (a, b) in [
            (got.rmsMean, want.rmsMean), (got.rmsMax, want.rmsMax),
            (got.fminMean, want.fminMean), (got.fminMin, want.fminMin),
            (got.fmaxMean, want.fmaxMean), (got.fmaxMax, want.fmaxMax)]:
          check abs(a - b) < 1e-12 * (1.0 + abs(b))
        if pass == 0:
          dt.update 0.0375
      planned.clear

  test "force extrema include exact zero magnitudes":
    let force = grt.toGvalue(zeroGaugeLike(g))
    let stats = force.forceRmsMinMax(float(g.len * lo.physVol))

    check stats.rms == 0.0
    check stats.fmin == 0.0
    check stats.fmax == 0.0
    let diag = force.forceStats
    discard diag.eval
    check diag.forceValues == stats
    expect(GraphError):
      discard grad(Gscalar(diag[0]), force)

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

  test "reversibility failure restores both initial leaves":
    let
      grt = initGraphRuntime()
      lo = lat.newLayout
      g = lo.newgauge
      p = lo.newgauge
      gc = actWilson(grt.toGvalue(6.0))
    var rng = lo.newRNGField(Philox4x64, 921260912u64)
    threads:
      g.random rng
      p.randomTAH rng
    var config = validRunConfig()
    config.gsteps = 1
    var
      armed, failed = false
      gi, pi: Ggauge
      ge, pe: graphGaugeShared.Gauge
    proc forward(v: Gvalue) =
      if armed:
        ge = gi.gaugeSnapshot
        pe = pi.gaugeSnapshot
        raiseError("requested reverse force failure")
      v.valCopy(v.inputs[0])
    proc force(x: Ggauge): Ggauge =
      let f = gaugeForce(gc, x)
      graphNode(f.gaugeNodeLike, [Gvalue(f)],
        Gfunc(forward: forward, bufferMode: bmFull, name: "reverse force failure"),
        "reverse force failure")
    let graph = buildTrajectoryGraph(grt, g, p, act(gc), config,
      buildTraining = false, force = force)
    gi = graph.initialState.gauge
    pi = graph.initialState.momentum
    let
      g0 = gi.gaugeSnapshot
      p0 = pi.gaugeSnapshot
    discard graph.finalState.hamiltonian.eval
    let
      gf = graph.finalState.gauge.gaugeSnapshot
      pf = graph.finalState.momentum.gaugeSnapshot
    # Arming leaves the forward caches current; reverse input updates trigger failure.
    armed = true
    try:
      graph.reversibilityCheck
    except GraphError as e:
      failed = true
      check e.msg == "requested reverse force failure"
    check failed
    norm2(grt.toGvalue(ge) - grt.toGvalue(gf)) :< 1e-26
    norm2(grt.toGvalue(pe) + grt.toGvalue(pf)) :< 1e-26
    norm2(gi - grt.toGvalue(g0)) :< 1e-26
    norm2(pi - grt.toGvalue(p0)) :< 1e-26

  test "planned proposals match retained values for integrators and loss branches":
    let
      lo = lat.newLayout
      g = lo.newgauge
      p = lo.newgauge
    var rng = lo.newRNGField(Philox4x64, 921260915u64)
    threads:
      g.random rng
      p.randomTAH rng
    for kind in IntegratorKind:
      for bias in [-2.0, 2.0]:
        var config = validRunConfig()
        config.gsteps = 1
        config.trajs = 1
        config.trajsTrain = 1
        config.trajsTrainlrWarm = 0
        config.trajsForceAcc = 1
        config.integratorCoeffs = parseIntegratorCoeffs(kind, [])
        let
          rt = initGraphRuntime()
          rr = initGraphRuntime()
          gc = actWilson(rt.toGvalue(6.0))
          rc = actWilson(rr.toGvalue(6.0))
        # An initial-energy offset selects each acceptance-loss branch while its
        # constant derivative leaves the integrator force unchanged.
        proc action(x: Ggauge): Gscalar =
          result = gaugeAction(gc, x)
          if x.gfunc == nil: result = result + bias
        proc referenceAction(x: Ggauge): Gscalar =
          result = gaugeAction(rc, x)
          if x.gfunc == nil: result = result + bias
        let
          graph = buildTrajectoryGraph(rt, g, p, action, config)
          reference = buildTrajectoryGraph(rr, g, p, referenceAction, config)
          view = adj(graph.finalState.gauge)
        var
          random = lo.newRNGField(Philox4x64, 921260916u64)
          referenceRandom = lo.newRNGField(Philox4x64, 921260916u64)
          serial: Philox4x64
        serial.seed(921260917u64, 0)
        reference.resampleMomentum(referenceRandom)
        let
          dh = reference.deltaHamiltonian.eval.sval
          loss = reference.lossExpr.eval.sval
          stats = reference.mdForces.mdForceStats
          finalGauge = reference.finalState.gauge.gaugeSnapshot
        check (dh < 0.0) == (bias > 0.0)
        var gradients = newSeq[float](reference.learnedParameters.len)
        for i, learned in reference.learnedParameters:
          gradients[i] = learned.gradientExpr.eval.sval
        var
          roots = @[Gvalue(graph.initialState.hamiltonian),
            Gvalue(graph.finalState.hamiltonian), Gvalue(view), Gvalue(graph.lossExpr)]
          proposalCalled, measureCalled = false
        for learned in graph.learnedParameters: roots.add learned.gradientExpr
        let originals = graphValues(roots)
        proc onProposal(traj: int; proposal: Proposal) =
          proposalCalled = true
          check traj == 1
          check abs(proposal.dH - dh) < 1e-10
          check abs(proposal.acc - exp(-dh)) < 1e-10 * (1.0 + exp(-dh))
          check abs(proposal.loss - loss) < 1e-11 * (1.0 + abs(loss))
          check proposal.gradients.len == gradients.len
          for i, expected in gradients:
            check abs(proposal.gradients[i] - expected) < 1e-9 * (1.0 + abs(expected))
          norm2(proposal.gauge - rt.toGvalue(finalGauge)) :< 1e-20
          norm2(proposal.view - adj(rt.toGvalue(finalGauge))) :< 1e-20
          for original in originals: check original.runCount == 0
          for learned in graph.learnedParameters:
            learned.node.update learned.node.sval + 0.001
        proc measure(traj: int; dH, acc: float; accepted: bool; forceStats: MdForceStats) =
          measureCalled = true
          check proposalCalled
          check accepted
          check forceStats.count == stats.count
          for (got, expected) in [
              (forceStats.rmsMean, stats.rmsMean), (forceStats.rmsMax, stats.rmsMax),
              (forceStats.fminMean, stats.fminMean), (forceStats.fminMin, stats.fminMin),
              (forceStats.fmaxMean, stats.fmaxMean), (forceStats.fmaxMax, stats.fmaxMax)]:
            check abs(got - expected) < 1e-11 * (1.0 + abs(expected))
          finalGauge.reunitGauge
          norm2(graph.initialState.gauge - rt.toGvalue(finalGauge)) :< 1e-20
          for original in originals: check original.runCount == 0
        runHmc(graph, config, random, serial, measure, onProposal, proposalView = view)
        check proposalCalled and measureCalled

  test "driver phases evaluate loss each time and gradients only for training":
    let
      rt = initGraphRuntime()
      lo = lat.newLayout
      g = lo.newgauge
      p = lo.newgauge
      gc = actWilson(rt.toGvalue(6.0))
    var random = lo.newRNGField(Philox4x64, 921260918u64)
    threads:
      g.random random
      p.randomTAH random
    var config = validRunConfig()
    config.gsteps = 1
    config.trajsThermo = 1
    config.trajs = 2
    config.trajsTrain = 1
    config.trajsTrainlrWarm = 0
    config.trajsForceAcc = 3
    config.revCheckFreq = 2
    var
      graph = buildTrajectoryGraph(rt, g, p, act(gc), config)
      proposalRuns, lossRuns, gradientRuns, proposals, measures: int
      serial: Philox4x64
      expected: graphGaugeShared.Gauge
      copied: seq[float]
    serial.seed(921260919u64, 0)
    proc proposalForward(v: Gvalue) =
      inc proposalRuns
      v.valCopy(v.inputs[0])
    proc lossForward(v: Gvalue) =
      inc lossRuns
      v.valCopy(v.inputs[0])
    proc gradientForward(v: Gvalue) =
      inc gradientRuns
      v.valCopy(v.inputs[0])
    graph.finalState.gauge = graphNode(graph.finalState.gauge.gaugeNodeLike,
      [Gvalue(graph.finalState.gauge)],
      Gfunc(forward: proposalForward, bufferMode: bmFull, name: "proposal count"), "proposal count")
    graph.lossExpr = graphNode(graph.lossExpr.scalarNodeLike, [Gvalue(graph.lossExpr)],
      Gfunc(forward: lossForward, bufferMode: bmFull, name: "loss count"), "loss count")
    for learned in mitems(graph.learnedParameters):
      learned.gradientExpr = graphNode(learned.gradientExpr.scalarNodeLike,
        [Gvalue(learned.gradientExpr)],
        Gfunc(forward: gradientForward, bufferMode: bmFull, name: "gradient count"), "gradient count")
    var trainer = initTrainingState(graph, config.weightDecay)
    proc onProposal(traj: int; proposal: Proposal) =
      inc proposals
      check proposalRuns == traj + traj div 2
      check lossRuns == traj
      check proposal.view == proposal.gauge
      let
        tau = float(config.gsteps) * graph.learnedParameters[0].node.sval
        loss = -min(1.0, proposal.acc) * tau * tau
      check abs(proposal.loss - loss) < 1e-12
      expected = proposal.gauge.gaugeSnapshot
      expected.reunitGauge
      if config.trajectoryPhase(traj) == tpTrain:
        check proposal.gradients.len == graph.learnedParameters.len
        copied = proposal.gradients
        trainer.trainStep(config, 1, proposal.gradients)
        check proposal.gradients == copied
      else:
        check proposal.gradients.len == 0
      check gradientRuns == (if traj < 2: 0 else: graph.learnedParameters.len)
      check graph.finalState.gauge.runCount == 0
      check graph.finalState.hamiltonian.runCount == 0
      check graph.lossExpr.runCount == 0
      for learned in graph.learnedParameters: check learned.gradientExpr.runCount == 0
    proc measure(traj: int; dH, acc: float; accepted: bool; forceStats: MdForceStats) =
      inc measures
      check proposals == traj
      check accepted
      norm2(graph.initialState.gauge - rt.toGvalue(expected)) :< 1e-26
    runHmc(graph, config, random, serial, measure, onProposal)
    check proposals == 3
    check measures == 3
    check gradientRuns == graph.learnedParameters.len
    check copied.len == graph.learnedParameters.len

  test "driver accept reject and forced acceptance preserve the chosen snapshot":
    let
      lo = lat.newLayout
      g = lo.newgauge
      p = lo.newgauge
    var rng = lo.newRNGField(Philox4x64, 921260920u64)
    threads:
      g.random rng
      p.randomTAH rng
    for mode in 0..2:
      let rt = initGraphRuntime()
      var config = validRunConfig()
      config.gsteps = 1
      config.trajs = 1
      config.trajsTrain = 0
      config.trajsTrainlrWarm = 0
      config.trajsForceAcc = (if mode == 2: 1 else: 0)
      var graph = buildTrajectoryGraph(rt, g, p,
        act(actWilson(rt.toGvalue(6.0))), config, buildTraining = false)
      graph.finalState.hamiltonian = graph.finalState.hamiltonian +
        (if mode == 0: -100.0 else: 100.0)
      var
        random = lo.newRNGField(Philox4x64, 921260921u64)
        serial: Philox4x64
        proposed: graphGaugeShared.Gauge
        proposals, measures: int
      serial.seed(921260922u64, 0)
      proc onProposal(traj: int; proposal: Proposal) =
        inc proposals
        check proposal.gradients.len == 0
        check graph.lossExpr == nil
        check graph.learnedParameters.len == 0
        check (proposal.dH < 0.0) == (mode == 0)
        proposed = proposal.gauge.gaugeSnapshot
        proposed.reunitGauge
        # The callback can write a borrowed result; commit uses its earlier copy.
        proposal.gauge.update g
      proc measure(traj: int; dH, acc: float; accepted: bool; forceStats: MdForceStats) =
        inc measures
        check proposals == 1
        check accepted == (mode != 1)
        let expected = if accepted: proposed else: g
        norm2(graph.initialState.gauge - rt.toGvalue(expected)) :< 1e-26
        check graph.finalState.hamiltonian.runCount == 0
        check graph.finalState.gauge.runCount == 0
      runHmc(graph, config, random, serial, measure, onProposal)
      check measures == 1

  test "seeded trajectories match retained acceptance state and RNG continuation":
    let
      lo = lat.newLayout
      g = lo.newgauge
      p = lo.newgauge
      rt = initGraphRuntime()
      rr = initGraphRuntime()
      shifts = [-100.0, 0.0, 100.0]
      shift = rt.toGvalue(shifts[0])
      referenceShift = rr.toGvalue(shifts[0])
    var rng = lo.newRNGField(Philox4x64, 921260925u64)
    threads:
      g.random rng
      p.randomTAH rng
    var config = validRunConfig()
    config.gsteps = 1
    config.trajs = shifts.len
    config.trajsTrain = 0
    config.trajsTrainlrWarm = 0
    config.revCheckFreq = 2
    var
      graph = buildTrajectoryGraph(rt, g, p,
        act(actWilson(rt.toGvalue(6.0))), config, buildTraining = false)
      reference = buildTrajectoryGraph(rr, g, p,
        act(actWilson(rr.toGvalue(6.0))), config, buildTraining = false)
      random = lo.newRNGField(Philox4x64, 921260926u64)
      referenceRandom = lo.newRNGField(Philox4x64, 921260926u64)
      serial, referenceSerial: Philox4x64
      delta, acc: float
      accepted: bool
      stats: MdForceStats
      proposals, measures: int
    serial.seed(921260927u64, 0)
    referenceSerial.seed(921260927u64, 0)
    # Finite offsets select acceptance and rejection. The middle trajectory has
    # no offset for its reverse check.
    graph.finalState.hamiltonian = graph.finalState.hamiltonian + shift
    reference.finalState.hamiltonian = reference.finalState.hamiltonian + referenceShift
    proc onProposal(traj: int; proposal: Proposal) =
      inc proposals
      reference.resampleMomentum(referenceRandom)
      let h0 = reference.initialState.hamiltonian.eval.sval
      delta = reference.finalState.hamiltonian.eval.sval - h0
      acc = exp(-delta)
      stats = reference.mdForces.mdForceStats
      let finalGauge = reference.finalState.gauge.gaugeSnapshot
      accepted = referenceSerial.uniform <= acc
      check abs(proposal.dH - delta) < 1e-10
      check abs(proposal.acc - acc) < 1e-11 * (1.0 + acc)
      check proposal.gradients.len == 0
      check norm2(proposal.gauge - rt.toGvalue(finalGauge)).eval.sval < 1e-20
      if traj == 1: check accepted
      if traj == shifts.len: check not accepted
      if traj mod config.revCheckFreq == 0:
        reference.reversibilityCheck
      if accepted:
        reference.commitAcceptedTrajectory(finalGauge)
    proc measure(traj: int; dH, probability: float; chosen: bool; forceStats: MdForceStats) =
      inc measures
      check proposals == traj
      check chosen == accepted
      check abs(dH - delta) < 1e-10
      check abs(probability - acc) < 1e-11 * (1.0 + acc)
      check forceStats.count == stats.count
      for (got, expected) in [
          (forceStats.rmsMean, stats.rmsMean), (forceStats.rmsMax, stats.rmsMax),
          (forceStats.fminMean, stats.fminMean), (forceStats.fminMin, stats.fminMin),
          (forceStats.fmaxMean, stats.fmaxMean), (forceStats.fmaxMax, stats.fmaxMax)]:
        check abs(got - expected) < 1e-11 * (1.0 + abs(expected))
      check norm2(graph.initialState.gauge -
        rt.toGvalue(reference.initialState.gauge.gaugeSnapshot)).eval.sval < 1e-20
      check norm2(graph.initialState.momentum -
        rt.toGvalue(reference.initialState.momentum.gaugeSnapshot)).eval.sval < 1e-20
      if traj < shifts.len:
        shift.update shifts[traj]
        referenceShift.update shifts[traj]
    runHmc(graph, config, random, serial, measure, onProposal)
    check proposals == shifts.len
    check measures == shifts.len
    for draw in 0..<3:
      check serial.uniform == referenceSerial.uniform

  test "planned reverse failure restores resampled leaves and preserves forward storage":
    let
      rt = initGraphRuntime()
      lo = lat.newLayout
      g = lo.newgauge
      p = lo.newgauge
      gc = actWilson(rt.toGvalue(6.0))
    var random = lo.newRNGField(Philox4x64, 921260923u64)
    threads:
      g.random random
      p.randomTAH random
    var config = validRunConfig()
    config.gsteps = 1
    config.trajs = 1
    config.trajsTrain = 0
    config.trajsTrainlrWarm = 0
    config.trajsForceAcc = 1
    config.revCheckFreq = 1
    var
      calls: int
      failed, proposalCalled, measureCalled = false
      gi, pi: Ggauge
      g0, p0, forwardStorage, forwardSnapshot: graphGaugeShared.Gauge
      serial: Philox4x64
    serial.seed(921260924u64, 0)
    proc forceForward(v: Gvalue) =
      if calls == 0:
        g0 = gi.gaugeSnapshot
        p0 = pi.gaugeSnapshot
      inc calls
      if calls > 2: raiseError("requested planned reverse force failure")
      v.valCopy(v.inputs[0])
    proc force(x: Ggauge): Ggauge =
      let f = gaugeForce(gc, x)
      graphNode(f.gaugeNodeLike, [Gvalue(f)],
        Gfunc(forward: forceForward, bufferMode: bmFull, name: "planned reverse failure"),
        "planned reverse failure")
    proc viewForward(v: Gvalue) =
      v.valCopy(v.inputs[0])
      forwardStorage = Ggauge(v).gval
      forwardSnapshot = Ggauge(v).gaugeSnapshot
    let graph = buildTrajectoryGraph(rt, g, p, act(gc), config,
      buildTraining = false, force = force)
    gi = graph.initialState.gauge
    pi = graph.initialState.momentum
    let view = graphNode(graph.finalState.gauge.gaugeNodeLike,
      [Gvalue(graph.finalState.gauge)],
      Gfunc(forward: viewForward, bufferMode: bmFull, name: "save forward storage"),
      "save forward storage")
    proc onProposal(traj: int; proposal: Proposal) =
      proposalCalled = true
    proc measure(traj: int; dH, acc: float; accepted: bool; forceStats: MdForceStats) =
      measureCalled = true
    try:
      runHmc(graph, config, random, serial, measure, onProposal, proposalView = view)
    except GraphError as e:
      failed = true
      check e.msg == "requested planned reverse force failure"
    check failed
    check calls == 3
    check not proposalCalled
    check not measureCalled
    check forwardSnapshot.len == g.len
    norm2(gi - rt.toGvalue(g0)) :< 1e-26
    norm2(pi - rt.toGvalue(p0)) :< 1e-26
    norm2(rt.toGvalue(forwardStorage) - rt.toGvalue(forwardSnapshot)) :< 1e-26
    check graph.finalState.hamiltonian.runCount == 0
    check graph.finalState.gauge.runCount == 0

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
    let gradients = @[0.5, -0.25]
    trainer.trainStep(config, 1, gradients)

    let after = trainer.parameterValues
    check after.len == before.len
    check after.len == graph.learnedParameters.len
    check after != before
    for learned in graph.learnedParameters:
      check learned.gradientExpr.runCount == 0
    check graph.finalState.hamiltonian.runCount == 0

  test "training step rejects indexes outside training phase":
    let config = validRunConfig()
    let graph = buildTrajectoryGraph(grt, g, p, act(actWilson(scalar.toGvalue(grt, 6.0))), config)
    var trainer = initTrainingState(graph, config.weightDecay)

    expect(GraphValueError):
      trainer.trainStep(config, 0, @[0.5, -0.25])
    expect(GraphValueError):
      trainer.trainStep(config, config.trajsTrain + 1, @[0.5, -0.25])
    expect(GraphValueError):
      trainer.trainStep(config, 1, @[0.5])

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
