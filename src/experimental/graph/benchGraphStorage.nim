## Graph storage measurements in fresh processes; see
## graph_validation.md for the root and counter contract and the
## required paired HMC fingerprint comparison.
import qex
import base/[alignedMem, omp]
import std/[math, tables]
when defined(macosx) or defined(linux):
  import std/posix
import ./[core, scalar, gauge, functional, multi, plan]
import ./gauge/types
import ./hmcgauge/[config, integrator, trajectory]

letParam:
  # chain, branch, gradient, or hmc
  `case` = "chain"
  # direct or planned
  mode = "planned"
  # chain only: number of scalar-times-gauge nodes
  steps = 24
  # calls per phase after the single warm call
  reps = (if `case` == "hmc": 3 else: 20)
  # chain only: allow the scale operator to reuse its gauge input
  inplace = false
  # chain, branch and gradient only; hmc starts with beta = 6
  alpha = 1.01
  # deterministic input increment in changed phases
  delta = 0.0001
let nth = block:
  letParam:
    # overrides OMP_NUM_THREADS
    threads = 1
  threads
let planned = case mode
  of "direct": false
  of "planned": true
  else: raiseValueError("mode must be direct or planned")

doAssert reps > 0, "reps must be positive"
doAssert nth > 0, "threads must be positive"
doAssert `case` in ["chain", "branch", "gradient", "hmc"], "case must be chain, branch, gradient, or hmc"
if `case` == "chain":
  doAssert steps > 0, "steps must be positive for case=chain"
ompSetNumThreads(nth.cint)
qexInit()
echo "graph-storage timing_run_kind=",
  (if defined(graphPlanMemory) or defined(nimAllocStats): "instrumented" else: "ordinary"),
  " graph_plan_memory=", defined(graphPlanMemory),
  " nim_allocation_counts=", defined(nimAllocStats),
  " raw_peak_scope=process-lifetime rss_peak_scope=process-lifetime",
  " teardown_measured=false"

proc forwards(rt: GraphRuntime): int =
  for s in rt.runStatsByNode.values:
    result += s.count

proc rssPeak(): int =
  when defined(macosx) or defined(linux):
    var ru: Rusage
    discard getrusage(RUSAGE_SELF, addr ru)
    when defined(macosx):
      result = int(ru.ru_maxrss)
    else:
      result = 1024 * int(ru.ru_maxrss)
  else:
    result = -1

proc chain() =
  letParam:
    # global lattice; the default is derived from an 8^4 local lattice
    lat = latticeFromLocalLattice(@[8,8,8,8],nRanks)
  let
    lo = lat.newLayout
    rt = initGraphRuntime()
    # A single external gauge allocation is preserved through every execution.
    x = Ggauge(runtime: rt, gval: lo.newGauge).assignStableNodeId
    a = rt.toGvalue(alpha)
  threads:
    for f in x.gval:
      f := 1.0
  x.updated
  let initialNorm = norm2(x.gval)
  GC_fullCollect()
  let
    baseUsed = getRawMemUsed()
    baseAllocated = getRawMemAllocated()
    basePeak = getRawMemMaxUsed()
    baseManaged = getOccupiedMem()

  proc scaleF(v: Gvalue) =
    let
      a = Gscalar(v.inputs[0])
      x = Ggauge(v.inputs[1])
      z = Ggauge(v)
    threads:
      for mu in 0..<z.gval.len:
        z.gval[mu] := a.sval * x.gval[mu]

  let fn = Gfunc(forward: scaleF, name: "storage scalar times gauge",
    bufferMode: bmFull, inplace: (if inplace: @[1] else: @[]))
  when defined(nimAllocStats):
    let initialAllocs = getAllocStats()
  tic("construct graph")
  var y = x
  for _ in 0..<steps:
    y = graphNode(y.gaugeNodeLike, @[Gvalue(a), Gvalue(y)], fn, fn.name)
  var p: GraphPlan
  if planned:
    p = plan(y)
  let constructionSeconds = getElapsedTime()
  toc("construction")
  let constructionBytes = getRawMemAllocated()-baseAllocated

  echo "graph-storage benchmark=fixed-chain mode=", mode,
    " inplace=", inplace, " rank=", myRank, " ranks=", nRanks,
    " threads=", nth, " lattice=", lat, " rank_geometry=", lo.rankGeom,
    " steps=", steps, " reps=", reps,
    " alpha=", alpha, " delta=", delta,
    " evaluation_roots=1 root_carrier=gauge roots=chain.final.gauge",
    " reference=analytic-norm"
  echo "graph-storage input_raw_bytes=", x.bufferBytes,
    " baseline_raw_current=", baseUsed, " baseline_raw_peak=", basePeak,
    " baseline_raw_allocated=", baseAllocated,
    " raw_gc_threshold=", getRawMemGcThreshold(), " input_norm2=", initialNorm
  echo "graph-storage construction_seconds=", constructionSeconds,
    " construction_raw_bytes=", constructionBytes,
    " current_above_inputs=", getRawMemUsed()-baseUsed,
    " process_raw_peak_above_baseline=", max(0, getRawMemMaxUsed()-baseUsed),
    " cumulative_raw_allocated=", getRawMemAllocated(),
    " managed_occupied=", getOccupiedMem(),
    " managed_occupied_above_inputs=", getOccupiedMem()-baseManaged,
    " managed_heap=", getTotalMem(), " process_rss_peak_bytes=", rssPeak(),
    " forwards=", forwards(rt)
  when defined(nimAllocStats):
    echo "graph-storage construction_nim_allocator_delta=", getAllocStats()-initialAllocs

  proc run(): Ggauge =
    if planned:
      discard p.eval()
      Ggauge(p[0])
    else:
      y.eval

  var inputScale = 1.0

  proc phase(label: string, count: int, updateAlpha: bool, updateGauge = false, offset = 0) =
    let
      priorForwards = forwards(rt)
      priorAllocated = getRawMemAllocated()
    when defined(nimAllocStats):
      let priorAllocs = getAllocStats()
    var
      secs = 0.0
      evalAllocated = 0
      checksum = 0.0
      analytic = 0.0
      maxRelativeError = 0.0
    for i in 0..<count:
      if updateAlpha:
        a.update(alpha + delta*float(offset+i+1))
      if updateGauge:
        inputScale = 1.0 + delta*float(offset+i+1)
        threads:
          for f in x.gval:
            f := inputScale
        x.updated
      let raw = getRawMemAllocated()
      tic("evaluate graph")
      let output = run()
      secs += getElapsedTime()
      toc("evaluation")
      evalAllocated += getRawMemAllocated()-raw
      # ||a^steps U||^2 = a^(2 steps) ||U||^2; verification is outside the timer.
      let value = norm2(output.gval)
      let want = initialNorm * inputScale*inputScale * pow(a.sval, float(2*steps))
      let err = abs(value-want)/max(1.0, abs(want))
      doAssert err < 1e-10, "graph chain disagrees with its analytic norm"
      maxRelativeError = max(maxRelativeError, err)
      checksum += value
      analytic += want
    let
      executed = forwards(rt)-priorForwards
      allocated = getRawMemAllocated()-priorAllocated
      current = getRawMemUsed()-baseUsed
      peak = max(0, getRawMemMaxUsed()-baseUsed)
      managed = getOccupiedMem()
    when defined(nimAllocStats):
      let phaseAllocs = getAllocStats()-priorAllocs
    if updateAlpha or updateGauge or label == "warm":
      doAssert executed == count*steps, "each changed chain node must run once"
    if label != "warm":
      doAssert evalAllocated == 0, "a warmed chain allocated raw storage"
      if not updateAlpha and not updateGauge:
        doAssert executed == 0, "an unchanged chain repeated forwards"
    let inputNorm = norm2(x.gval)
    let wantInput = initialNorm * inputScale*inputScale
    doAssert abs(inputNorm-wantInput) < 1e-11*max(1.0, wantInput),
      "chain execution modified the preserved input"
    # Record both allocator occupancy and occupancy after collection. Timing above
    # includes evaluation only; it excludes checksum reductions and this collection.
    GC_fullCollect()
    echo "graph-storage phase=", label, " calls=", count,
      " eval_seconds=", secs, " seconds_per_call=", secs/float(count),
      " eval_raw_allocated=", evalAllocated, " phase_raw_allocated=", allocated,
      " allocated_above_inputs=", getRawMemAllocated()-baseAllocated,
      " cumulative_raw_allocated=", getRawMemAllocated(),
      " current_above_inputs=", current,
      " current_after_gc_above_inputs=", getRawMemUsed()-baseUsed,
      " process_raw_peak_above_baseline=", peak,
      " managed_occupied_before_gc=", managed,
      " managed_occupied_after_gc=", getOccupiedMem(),
      " managed_heap=", getTotalMem(), " process_rss_peak_bytes=", rssPeak(),
      " forwards=", executed, " cumulative_forwards=", forwards(rt),
      " checksum=", checksum, " analytic_checksum=", analytic,
      " max_relative_error=", maxRelativeError
    when defined(nimAllocStats):
      echo "graph-storage phase=", label, " phase_nim_allocator_delta=", phaseAllocs
    if p != nil:
      echo "graph-storage arena phase=", label,
        " runs=", p.stats.runs, " forwards=", p.stats.forwards,
        " reuses=", p.stats.reuses, " buffers=", p.stats.buffers,
        " bytes=", p.stats.arenaBytes, " peak_live_bytes=", p.stats.peakLiveBytes,
        " workspaces=", p.stats.workspaces, " workspace_bytes=", p.stats.workspaceBytes

  phase("warm", 1, false)
  phase("unchanged", reps, false)
  phase("updated", reps, true)
  phase("all-inputs-changed", reps, true, true, reps)
  phase("partly-changed", reps, false, true, 2*reps)
  rt.echoRunStats()

proc fill(x: Ggauge, value: float) =
  threads:
    for f in x.gval:
      f := value
  x.updated

proc fillHmc(x: Ggauge, scale: float, momentum: bool) =
  let g = x.gval
  const nc = g[0][0].nrows
  threads:
    for mu in 0..<g.len:
      g[mu] := (if momentum: 0.0 else: 1.0)
      for e in g[mu]:
        let t = scale * math.sin(float(3*e+2*mu+1))
        if momentum:
          g[mu][e][0,0].im := t
          when nc > 1:
            g[mu][e][1,1].im := -t
        else:
          g[mu][e][0,0].re := math.cos(t)
          g[mu][e][0,0].im := math.sin(t)
          when nc > 1:
            g[mu][e][1,1].re := math.cos(t)
            g[mu][e][1,1].im := -math.sin(t)
  x.updated

proc fingerprint(v: Gvalue): seq[float] =
  if v of Gscalar:
    return @[Gscalar(v).sval]
  doAssert v of Ggauge, "benchmark fingerprints require scalar or gauge roots"
  let g = Ggauge(v).gval
  const nc = g[0][0].nrows
  var re, im: float
  threads:
    var sr, si = 0.0
    for mu in 0..<g.len:
      for e in g[mu]:
        for row in 0..<nc:
          for col in 0..<nc:
            let w = float(1+mu+3*e+7*row+11*col)
            sr += w * g[mu][e][row,col].re.simdSum
            si += w * g[mu][e][row,col].im.simdSum
    sr.threadRankSum
    si.threadRankSum
    threadMaster:
      re = sr
      im = si
  @[norm2(g), re, im]

proc difference(a, b: openArray[float]): float =
  doAssert a.len == b.len, "benchmark reference arity differs"
  for i in 0..<a.len:
    result = max(result, abs(a[i]-b[i])/max(1.0, abs(b[i])))

proc matched() =
  letParam:
    # global lattice; the default is derived from a 4^4 local lattice
    lat = latticeFromLocalLattice(@[4,4,4,4],nRanks)
  let
    lo = lat.newLayout
    rt = initGraphRuntime()
  var
    x = Ggauge(runtime: rt, gval: lo.newGauge).assignStableNodeId
    y = Ggauge(runtime: rt, gval: lo.newGauge).assignStableNodeId
    xv = 1.0
    yv = 0.75
    av = (if `case` == "hmc": 6.0 else: alpha)
    bv = 0.625
  let
    a = rt.toGvalue(av)
    b = if `case` == "hmc": nil else: rt.toGvalue(bv)
    sel = if `case` == "branch": rt.toGvalue(1) else: nil
  if `case` == "hmc":
    x.fillHmc(0.05, false)
    y.fillHmc(0.35, true)
  else:
    x.fill(xv)
    y.fill(yv)
  var
    wantX = fingerprint(x)
    wantY = fingerprint(y)
  let unit = wantX
  GC_fullCollect()
  let
    baseUsed = getRawMemUsed()
    baseAllocated = getRawMemAllocated()
    basePeak = getRawMemMaxUsed()
    baseManaged = getOccupiedMem()
  when defined(nimAllocStats):
    let initialAllocs = getAllocStats()
  var
    roots: seq[Gvalue]
    names: seq[string]
    traj: TrajectoryGraph
    initialParams: seq[float]
  proc root(name: string, v: Gvalue) =
    names.add name
    roots.add v

  tic("construct graph")
  case `case`
  of "branch":
    let
      v = Ggauge(x.newOneOf)
      fn = cond(sel, lambda(v, a*v*v+y), lambda(v, b*v+y*y))
      z = Ggauge(apply(fn, x))
    root("application", z)
    root("independent", b*y)
    root("application.norm2", z.norm2)
  of "gradient":
    let
      z = a*x+b*y
      loss = z.norm2
    root("value", z)
    root("loss", loss)
    root("dLoss/dAlpha", loss.grad(a))
    root("dLoss/dBeta", loss.grad(b))
    root("dLoss/dX", loss.grad(x))
    root("dLoss/dY", loss.grad(y))
  of "hmc":
    let
      coeff = actWilson(a)
      cfg = RunConfig(dt: 0.05, gsteps: 1, trajs: 1, trajsTrain: 1,
        integratorCoeffs: parseIntegratorCoeffs(ik2MN, []))
    proc action(v: Ggauge): Gscalar = gaugeAction(coeff, v)
    traj = buildTrajectoryGraph(rt, x.gval, y.gval, action, cfg, buildTraining = true)
    # The current constructor takes numerical gauges and owns copies of its leaves.
    x = traj.initialState.gauge
    y = traj.initialState.momentum
    root("initial.H", traj.initialState.hamiltonian)
    root("initial.S", traj.initialState.gaugeAction)
    root("initial.T", traj.initialState.kinetic)
    for i, force in traj.mdForces:
      root("force[" & $i & "].norm2", force.norm2)
    root("final.gauge", traj.finalState.gauge)
    root("final.momentum", traj.finalState.momentum)
    root("final.H", traj.finalState.hamiltonian)
    root("final.S", traj.finalState.gaugeAction)
    root("final.T", traj.finalState.kinetic)
    root("dH", traj.deltaHamiltonian)
    root("loss", traj.lossExpr)
    for lp in traj.learnedParameters:
      initialParams.add lp.node.sval
      root("dLoss/d" & lp.name, lp.gradientExpr)
  else:
    raiseValueError("unknown matched-work case")
  # The carrier itself is part of the matched work in both modes.
  let bundle = multiValues("storage " & `case` & " roots", roots)
  let p = if planned: plan(bundle) else: nil
  let
    constructionSeconds = getElapsedTime()
    constructionBytes = getRawMemAllocated()-baseAllocated
    constructionUsed = getRawMemUsed()
    constructionManaged = getOccupiedMem()
  toc("construction")
  when defined(nimAllocStats):
    let constructionAllocs = getAllocStats()-initialAllocs
  GC_fullCollect()
  echo "graph-storage benchmark=matched-work case=", `case`, " mode=", mode,
    " inplace=stock-operator-contracts",
    " rank=", myRank, " ranks=", nRanks, " threads=", nth,
    " lattice=", lat, " rank_geometry=", lo.rankGeom,
    " reps=", reps, " delta=", delta,
    " evaluation_roots=1 root_carrier=multiValues root_slots=", names.len,
    " roots=", names,
    " ignored_options=", (if `case` == "hmc": "steps,inplace,alpha" else: "steps,inplace")
  echo "graph-storage fingerprint_layout=scalar:value,gauge:norm2+weighted_re+weighted_im",
    " checksum_weight=call_index_plus_one",
    " reference=", (if `case` == "hmc": "energy-loss-identities" else: "analytic"),
    " fingerprint_comparison=", (if `case` == "hmc": "external-required" else: "analytic-in-process")
  if `case` == "hmc":
    echo "graph-storage action=Wilson beta=6 integrator=2MN dt=0.05 gsteps=1",
      " workload=proposal-graph diagnostics=force_norm2 final_momentum=included",
      " training=loss+dt+lambda-gradients",
      " all_updates=gauge,momentum,beta,dt,lambda partly_updates=momentum",
      " construction_includes_input_copies=true"
  else:
    echo "graph-storage alpha=", av, " beta=", bv, " x=", xv, " y=", yv,
      " all_updates=x,y,alpha,beta", (if `case` == "branch": ",selector" else: ""),
      " partly_updates=x"
  echo "graph-storage input_raw_bytes=", x.bufferBytes+y.bufferBytes,
    " baseline_raw_current=", baseUsed, " baseline_raw_peak=", basePeak,
    " baseline_raw_allocated=", baseAllocated, " baseline_managed_occupied=", baseManaged,
    " raw_gc_threshold=", getRawMemGcThreshold()
  echo "graph-storage construction_seconds=", constructionSeconds,
    " construction_raw_bytes=", constructionBytes,
    " current_above_inputs=", constructionUsed-baseUsed,
    " current_after_gc_above_inputs=", getRawMemUsed()-baseUsed,
    " process_raw_peak_above_baseline=", max(0, getRawMemMaxUsed()-baseUsed),
    " cumulative_raw_allocated=", getRawMemAllocated(),
    " managed_occupied_before_gc=", constructionManaged,
    " managed_occupied_after_gc=", getOccupiedMem(), " managed_heap=", getTotalMem(),
    " process_rss_peak_bytes=", rssPeak(), " forwards=", forwards(rt)
  when defined(nimAllocStats):
    echo "graph-storage construction_nim_allocator_delta=", constructionAllocs

  proc change(label: string, i: int) =
    if label notin ["all-inputs-changed", "partly-changed"]:
      return
    let k = float((if label == "partly-changed": reps else: 0)+i+1)
    if `case` == "hmc":
      y.fillHmc(0.35+delta*k, true)
      if label == "all-inputs-changed":
        x.fillHmc(0.05+delta*k, false)
        a.update(6.0+delta*k)
        for j, lp in traj.learnedParameters:
          lp.node.update(initialParams[j]+delta*k*float(j+1))
    else:
      xv = 1.0+delta*k
      x.fill(xv)
      if label == "all-inputs-changed":
        yv = 0.75+2.0*delta*k
        av = alpha+delta*k
        bv = 0.625+3.0*delta*k
        y.fill(yv)
        a.update(av)
        b.update(bv)
        if `case` == "branch":
          sel.update(1-i mod 2)
    wantX = fingerprint(x)
    wantY = fingerprint(y)

  proc scaled(value: float): seq[float] =
    @[unit[0]*value*value, unit[1]*value, unit[2]*value]

  proc reference(): seq[seq[float]] =
    let d = unit[0]
    if `case` == "branch":
      let z = if sel.ival != 0: av*xv*xv+yv else: bv*xv+yv*yv
      return @[scaled(z), scaled(bv*yv), @[d*z*z]]
    let z = av*xv+bv*yv
    @[scaled(z), @[d*z*z], @[2.0*d*xv*z], @[2.0*d*yv*z],
      scaled(2.0*av*z), scaled(2.0*bv*z)]

  proc run(): Gmulti =
    if p != nil:
      discard p.eval()
      Gmulti(p[0])
    else:
      bundle.eval

  proc phase(label: string, count: int) =
    let
      priorForwards = forwards(rt)
      priorAllocated = getRawMemAllocated()
    when defined(nimAllocStats):
      let priorAllocs = getAllocStats()
    var
      secs = 0.0
      evalAllocated = 0
      maxRelativeError = 0.0
      checksums = newSeq[seq[float]](names.len)
      analytic = newSeq[seq[float]](names.len)
      last: seq[seq[float]]
    for i in 0..<count:
      change(label, i)
      let raw = getRawMemAllocated()
      tic("evaluate graph")
      let output = run()
      secs += getElapsedTime()
      toc("evaluation")
      evalAllocated += getRawMemAllocated()-raw
      # All numerical reductions, references and input checks are outside timers.
      last = newSeq[seq[float]](names.len)
      let want = if `case` == "hmc": newSeq[seq[float]]() else: reference()
      for k in 0..<names.len:
        last[k] = fingerprint(output.storedSlot(k))
        if checksums[k].len == 0:
          checksums[k] = newSeq[float](last[k].len)
          analytic[k] = newSeq[float](last[k].len)
        for j, value in last[k]:
          doAssert classify(value) notin {fcNan, fcInf, fcNegInf},
            "benchmark root fingerprint is not finite: " & names[k]
          checksums[k][j] += float(i+1)*value
          if `case` != "hmc":
            analytic[k][j] += float(i+1)*want[k][j]
        if `case` != "hmc":
          let err = difference(last[k], want[k])
          doAssert err < 1e-10, "benchmark root disagrees with its analytic reference: " & names[k]
          maxRelativeError = max(maxRelativeError, err)
      if `case` == "hmc":
        let
          f = 3+traj.mdForces.len
          h0 = last[0][0]
          h1 = last[f+2][0]
          dh = h1-h0
          dt = traj.learnedParameters[0].node.sval
          loss = -min(1.0, math.exp(-dh))*dt*dt
          err = difference([h0, h1, last[f+5][0], last[f+6][0]],
            [last[1][0]+last[2][0], last[f+3][0]+last[f+4][0], dh, loss])
        doAssert err < 1e-10, "HMC energy/loss identities differ"
        maxRelativeError = max(maxRelativeError, err)
      doAssert difference(fingerprint(x), wantX) < 1e-11 and
        difference(fingerprint(y), wantY) < 1e-11, "execution modified a preserved input"
    let
      executed = forwards(rt)-priorForwards
      allocated = getRawMemAllocated()-priorAllocated
      current = getRawMemUsed()-baseUsed
      peak = max(0, getRawMemMaxUsed()-baseUsed)
      managed = getOccupiedMem()
    when defined(nimAllocStats):
      let phaseAllocs = getAllocStats()-priorAllocs
    if label == "unchanged":
      doAssert executed == 0, "an unchanged root bundle repeated forwards"
      doAssert evalAllocated == 0, "an unchanged root bundle allocated raw storage"
    if label != "unchanged":
      doAssert executed > 0, "changed root work did not execute"
    # A direct run can cache untouched subgraphs after partial input changes.
    # Report actual forwards; equality across modes is not an incremental contract.
    GC_fullCollect()
    echo "graph-storage case=", `case`, " phase=", label, " calls=", count,
      " eval_seconds=", secs, " seconds_per_call=", secs/float(count),
      " eval_raw_allocated=", evalAllocated, " phase_raw_allocated=", allocated,
      " allocated_above_inputs=", getRawMemAllocated()-baseAllocated,
      " cumulative_raw_allocated=", getRawMemAllocated(),
      " current_above_inputs=", current,
      " current_after_gc_above_inputs=", getRawMemUsed()-baseUsed,
      " process_raw_peak_above_baseline=", peak,
      " managed_occupied_before_gc=", managed,
      " managed_occupied_after_gc=", getOccupiedMem(), " managed_heap=", getTotalMem(),
      " process_rss_peak_bytes=", rssPeak(),
      " forwards=", executed, " cumulative_forwards=", forwards(rt),
      " error_reference=", (if `case` == "hmc": "energy-loss-identities" else: "analytic"),
      " max_relative_error=", maxRelativeError
    when defined(nimAllocStats):
      echo "graph-storage phase=", label, " phase_nim_allocator_delta=", phaseAllocs
    for k, name in names:
      echo "graph-storage case=", `case`, " phase=", label,
        " root=", k, " name=", name, " fingerprint_checksum=", checksums[k],
        " fingerprint_last=", last[k]
      if `case` != "hmc":
        echo "graph-storage case=", `case`, " phase=", label,
          " root=", k, " analytic_checksum=", analytic[k]
    if p != nil:
      echo "graph-storage arena case=", `case`, " phase=", label,
        " runs=", p.stats.runs, " forwards=", p.stats.forwards,
        " reuses=", p.stats.reuses, " buffers=", p.stats.buffers,
        " bytes=", p.stats.arenaBytes, " peak_live_bytes=", p.stats.peakLiveBytes,
        " workspaces=", p.stats.workspaces, " workspace_bytes=", p.stats.workspaceBytes

  phase("warm", 1)
  phase("unchanged", reps)
  phase("all-inputs-changed", reps)
  phase("partly-changed", reps)
  rt.echoRunStats()

if `case` == "chain":
  chain()
else:
  matched()
qexFinalize()
