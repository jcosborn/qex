import ../[core, scalar, multi, gauge]
import layout, physics/qcdTypes   # threads / simd reductions for forceRmsMinMax
from math import sqrt

type
  IntegratorKind* = enum
    ik2MN, ik4MN3F1GP, ik4MN5F2GP, ik2MNp
  IntegratorCoeffs* = object
    kind*: IntegratorKind
    rho*: float
    theta*: float
    vtheta*: float
    lambda*: float
    chi*: float
    xi*: float
  LearnedParameter* = object
    name*: string
    node*: Gscalar
    gradientExpr*: Gscalar
  IntegrationEventKind* = enum
    ieForce, ieKick, ieDrift, ieShift
  IntegrationEvent* = object
    ## Graph-node references, not stored values. A force event holds the gauge
    ## and the momentum before its kick; the kick holds the gauge of the state
    ## and the force node with the momentum after it. A shift event holds the
    ## gauge displaced along a force for a force-gradient kick, that force and
    ## the displacement; the force event that follows is at the shifted gauge.
    kind*: IntegrationEventKind
    step*: int # Repeated loop label: 0 for the leading events, i inside the
               # loop including the seam, n for the final events.
    gauge*, momentum*: Ggauge
    force*: Ggauge
    coefficient*: Gscalar # Signed applied coefficient; nil for force events.
  IntegrationResult* = object
    gauge*: Ggauge
    momentum*: Ggauge
    learnedCoeffs*: seq[LearnedParameter]
    forces*: seq[Ggauge] # Force nodes in integration order for every kind; the
                         # one home of force statistics.
    trace*: seq[IntegrationEvent] # Every event when traced; empty otherwise.
  MdForceStats* = object
    count*: int
    rmsMean*, rmsMax*: float
    fminMean*, fminMin*: float
    fmaxMean*, fmaxMax*: float
  GaugeAction* = proc(g: Ggauge): Gscalar
    ## S(g); the integrator differentiates it for the MD force.
  GaugeForceFn* = proc(g: Ggauge): Ggauge

proc gradForce*(action: GaugeAction, g: Ggauge, force: GaugeForceFn = nil): Ggauge =
  ## F = projectTAH(grad(S(g),g)*g†).
  ## `force`, when given, computes the same F directly.
  if force.isNil:
    contractProjTAH(grad(action(g), g), g)
  else:
    force(g)

proc gradForce(action: GaugeAction; g: Ggauge; forces: var seq[Ggauge]; force: GaugeForceFn): Ggauge =
  result = gradForce(action, g, force)
  forces.add result

proc forceRmsMinMaxValue(force: Ggauge; dof: float): tuple[rms, fmin, fmax: float] =
  let fs = force.gval
  var s2 = 0.0
  var n2 = 1e300
  var m2 = 0.0
  threads:
    var ls = 0.0
    var ln = 1e300
    var lm = 0.0
    for mu in 0..<fs.len:
      for x in fs[mu]:
        let n = fs[mu][x].norm2
        ls += n.simdSum
        let nn = n.simdMin
        if nn < ln: ln = nn
        let mm = n.simdMax
        if lm < mm: lm = mm
    ls.threadRankSum
    ln = -ln
    ln.threadRankMax  # min(x) = -max(-x)
    lm.threadRankMax
    threadSingle:
      s2 = ls
      n2 = -ln
      m2 = lm
  (rms: sqrt(s2/dof), fmin: sqrt(n2), fmax: sqrt(m2))

proc forceRmsMinMax*(force: Ggauge; dof: float): tuple[rms, fmin, fmax: float] =
  ## RMS, min, and max over links of the MD force magnitude |F| at the force node's
  ## current value. `dof` is the total link degrees of freedom (Σ_μ vol).
  discard force.eval
  forceRmsMinMaxValue(force, dof)

proc forceStatsForward(v: Gvalue) =
  let
    force = Ggauge(v.inputs[0])
    dof = float(force.gval.len * force.gval[0].l.physVol)
    f = forceRmsMinMaxValue(force, dof)
    z = Gmulti(v)
  Gscalar(z.storedSlot(0)).sval = f.rms
  Gscalar(z.storedSlot(1)).sval = f.fmin
  Gscalar(z.storedSlot(2)).sval = f.fmax

let forceStatsFunc = Gfunc(
  forward: forceStatsForward, bufferMode: bmFull, name: "forceStats")

proc forceStats*(force: Ggauge): Gmulti =
  ## Diagnostic sink with RMS/min/max scalar slots; differentiation is undefined.
  ## Schedule each root near its force producer to avoid retaining force history.
  let s = force.scalarNodeLike
  newMultiOutputNode([Gvalue(s), Gvalue(s), Gvalue(s)], [Gvalue(force)],
    forceStatsFunc, "forceStats")

proc mdForceStats*(forces: openArray[tuple[rms, fmin, fmax: float]]): MdForceStats =
  ## Aggregate copied RMS/min/max triples from integrator forces.
  if forces.len == 0:
    raiseValueError("MD force statistics require at least one force")
  result.count = forces.len
  result.rmsMax = forces[0].rms
  result.fminMin = forces[0].fmin
  result.fmaxMax = forces[0].fmax
  for f in forces:
    result.rmsMean += f.rms
    result.fminMean += f.fmin
    result.fmaxMean += f.fmax
    result.rmsMax = max(result.rmsMax, f.rms)
    result.fminMin = min(result.fminMin, f.fmin)
    result.fmaxMax = max(result.fmaxMax, f.fmax)
  result.rmsMean /= float(result.count)
  result.fminMean /= float(result.count)
  result.fmaxMean /= float(result.count)

proc mdForceStats*(forces: openArray[Ggauge]): MdForceStats =
  ## Statistics from retained integrator force nodes.
  if forces.len == 0:
    raiseValueError("MD force statistics require at least one force")
  discard forces[^1].eval
  let dof = float(forces[0].gval.len * forces[0].gval[0].l.physVol)
  var values = newSeq[tuple[rms, fmin, fmax: float]](forces.len)
  for i, force in forces:
    values[i] = forceRmsMinMaxValue(force, dof)
  values.mdForceStats

proc requireCoeffCountOrDefault(label: string,
                                values: openArray[float],
                                expected: int) =
  if values.len notin {0, expected}:
    raiseValueError(
      label & " expects either 0 or " & $expected &
      " coefficient values, got " & $values.len)

proc parseIntegratorKind*(name: string): IntegratorKind =
  case name
  of "2MN":
    ik2MN
  of "2MNp":
    ik2MNp
  of "4MN3F1GP":
    ik4MN3F1GP
  of "4MN5F2GP":
    ik4MN5F2GP
  else:
    raiseValueError("unknown intalg: " & name & "; expected 2MN, 2MNp, 4MN3F1GP, or 4MN5F2GP")

proc parseIntegratorCoeffs*(kind: IntegratorKind,
                            values: openArray[float]): IntegratorCoeffs =
  case kind
  of ik2MN, ik2MNp:
    requireCoeffCountOrDefault(if kind == ik2MN: "2MN" else: "2MNp", values, 1)
    result = IntegratorCoeffs(kind: kind)
    result.lambda = if values.len == 0: 0.1931833275037836 else: values[0]
  of ik4MN3F1GP:
    # Force-gradient family: defaults are all-or-default by design; partial
    # positional completion is unsupported. Keep the derived formulas explicit.
    requireCoeffCountOrDefault("4MN3F1GP", values, 3)
    result = IntegratorCoeffs(kind: ik4MN3F1GP)
    if values.len == 0:
      result.lambda = 0.2470939580390842
      result.theta = 0.5 - 1.0 / sqrt(24.0 * result.lambda)
      let numer = 1.0 - sqrt(6.0 * result.lambda) * (1.0 - result.lambda)
      let scale = 20.0 / (1.0 - 2.0 * result.lambda)
      result.chi = (numer / 12.0) * scale
    else:
      result.lambda = values[0]
      result.theta = values[1]
      result.chi = values[2]
  of ik4MN5F2GP:
    requireCoeffCountOrDefault("4MN5F2GP", values, 5)
    result = IntegratorCoeffs(kind: ik4MN5F2GP)
    if values.len == 0:
      result.rho = 0.06419108866816235
      result.theta = 0.1919807940455741
      result.vtheta = 0.1518179640276466
      result.lambda = 0.2158369476787619
      # Keep the scale factors explicit so the learned coefficient formula is auditable.
      result.xi = 0.0009628905212024874 * (2.0 / result.lambda * 20.0)
    else:
      result.rho = values[0]
      result.theta = values[1]
      result.vtheta = values[2]
      result.lambda = values[3]
      result.xi = values[4]

type Stepper = object
  ## One MD state with the shared kick, drift and shift bookkeeping. Kicks and
  ## shifts apply the negated coefficient: p -= c F and g_s = exp(-c F) g.
  action: GaugeAction
  force: GaugeForceFn
  g, p: Ggauge
  forces: seq[Ggauge]
  trace: seq[IntegrationEvent]
  traced: bool

proc event(s: var Stepper; kind: IntegrationEventKind; step: int; force: Ggauge = nil;
           coefficient: Gscalar = nil; gauge: Ggauge = nil) =
  if s.traced:
    s.trace.add IntegrationEvent(kind: kind, step: step,
      gauge: (if gauge == nil: s.g else: gauge), momentum: s.p,
      force: force, coefficient: coefficient)

proc kick(s: var Stepper; c: Gscalar; step: int) =
  let f = gradForce(s.action, s.g, s.forces, s.force)
  s.event(ieForce, step, f)
  let mc = -c
  s.p = axpy(mc, f, s.p)
  s.event(ieKick, step, f, mc)

proc drift(s: var Stepper; c: Gscalar; step: int) =
  s.g = axexpmuly(c, s.p, s.g)
  s.event(ieDrift, step, coefficient = c)

proc gradKick(s: var Stepper; c, shift: Gscalar; step: int) =
  ## Kick with the force at the gauge displaced against its own force by shift.
  let fg = gradForce(s.action, s.g, s.forces, s.force)
  s.event(ieForce, step, fg)
  let ms = -shift
  let gs = axexpmuly(ms, fg, s.g)
  s.event(ieShift, step, fg, ms, gs)
  let f = gradForce(s.action, gs, s.forces, s.force)
  s.event(ieForce, step, f, gauge = gs)
  let mc = -c
  s.p = axpy(mc, f, s.p)
  s.event(ieKick, step, f, mc)

proc integrate2MN(action: GaugeAction,
                  g0, p0: Ggauge,
                  dt: Gscalar,
                  n: int,
                  coeffs: IntegratorCoeffs,
                  force: GaugeForceFn,
                  traced: bool): IntegrationResult =
  ## Second-order minimal norm, O(a) [I(h) O(b) I(h) O(2a)]^(n-1) I(h) O(b) I(h) O(a)
  ## with a = lambda dt, b = (1-2 lambda) dt, h = dt/2: O is the kick and I the
  ## drift for ik2MNp (momentum first), the reverse for ik2MN.
  let lambda = toGvalue(dt.runtime, coeffs.lambda)
  let first = lambda * dt
  let between = 2.0 * first
  let middle = dt - between
  let half = 0.5 * dt
  let mf = coeffs.kind == ik2MNp
  var s = Stepper(action: action, force: force, g: g0, p: p0, traced: traced)
  proc outer(c: Gscalar; step: int) =
    if mf: s.kick(c, step) else: s.drift(c, step)
  proc inner(step: int) =
    if mf: s.drift(half, step) else: s.kick(half, step)
  outer(first, 0)
  for i in 0..<n:
    inner(i)
    outer(middle, i)
    inner(i)
    if i+1 < n: outer(between, i)
  outer(first, n)
  IntegrationResult(gauge: s.g, momentum: s.p,
    learnedCoeffs: @[LearnedParameter(name: "lambda", node: lambda)],
    forces: s.forces, trace: s.trace)

proc integrate4MN3F1GP(action: GaugeAction,
                       g0, p0: Ggauge,
                       dt: Gscalar,
                       n: int,
                       coeffs: IntegratorCoeffs,
                       force: GaugeForceFn,
                       traced: bool): IntegrationResult =
  let lambda = toGvalue(dt.runtime, coeffs.lambda)
  let theta = toGvalue(dt.runtime, coeffs.theta)
  let chi = toGvalue(dt.runtime, coeffs.chi)
  let a0 = theta * dt
  let a02 = 2.0 * a0
  let a1 = 0.5 * dt - a0
  let b0 = lambda * dt
  let b1 = dt - 2.0 * b0
  let c1 = 0.1 * chi * (dt * dt)
  var s = Stepper(action: action, force: force, g: g0, p: p0, traced: traced)
  for i in 0..<n:
    s.drift(if i == 0: a0 else: a02, i)
    s.kick(b0, i)
    s.drift(a1, i)
    s.gradKick(b1, c1, i)
    s.drift(a1, i)
    s.kick(b0, i)
  s.drift(a0, n)
  IntegrationResult(gauge: s.g, momentum: s.p,
    learnedCoeffs: @[
      LearnedParameter(name: "lambda", node: lambda),
      LearnedParameter(name: "theta", node: theta),
      LearnedParameter(name: "chi", node: chi)],
    forces: s.forces, trace: s.trace)

proc integrate4MN5F2GP(action: GaugeAction,
                       g0, p0: Ggauge,
                       dt: Gscalar,
                       n: int,
                       coeffs: IntegratorCoeffs,
                       force: GaugeForceFn,
                       traced: bool): IntegrationResult =
  let rho = toGvalue(dt.runtime, coeffs.rho)
  let theta = toGvalue(dt.runtime, coeffs.theta)
  let vtheta = toGvalue(dt.runtime, coeffs.vtheta)
  let lambda = toGvalue(dt.runtime, coeffs.lambda)
  let xi = toGvalue(dt.runtime, coeffs.xi)
  let a0 = rho * dt
  let a02 = 2.0 * a0
  let a1 = theta * dt
  let a2 = (0.5 - (theta + rho)) * dt
  let b1 = lambda * dt
  let b0 = vtheta * dt
  let b2 = (1.0 - 2.0 * (lambda + vtheta)) * dt
  let c1 = 0.05 * xi * (dt * dt)
  var s = Stepper(action: action, force: force, g: g0, p: p0, traced: traced)
  for i in 0..<n:
    s.drift(if i == 0: a0 else: a02, i)
    s.kick(b0, i)
    s.drift(a1, i)
    s.gradKick(b1, c1, i)
    s.drift(a2, i)
    s.kick(b2, i)
    s.drift(a2, i)
    s.gradKick(b1, c1, i)
    s.drift(a1, i)
    s.kick(b0, i)
  s.drift(a0, n)
  IntegrationResult(gauge: s.g, momentum: s.p,
    learnedCoeffs: @[
      LearnedParameter(name: "rho", node: rho),
      LearnedParameter(name: "theta", node: theta),
      LearnedParameter(name: "vtheta", node: vtheta),
      LearnedParameter(name: "lambda", node: lambda),
      LearnedParameter(name: "xi", node: xi)],
    forces: s.forces, trace: s.trace)

proc integrateGauge*(action: GaugeAction,
                     g0: Ggauge,
                     p0: Ggauge,
                     dt: Gscalar,
                     n: int,
                     coeffs: IntegratorCoeffs,
                     force: GaugeForceFn = nil,
                     trace = false): IntegrationResult =
  ## trace records every force, kick, drift and shift as graph-node events.
  if n <= 0:
    raiseValueError("integrator step count must be >= 1, got " & $n)
  case coeffs.kind
  of ik2MN, ik2MNp:
    integrate2MN(action, g0, p0, dt, n, coeffs, force, trace)
  of ik4MN3F1GP:
    integrate4MN3F1GP(action, g0, p0, dt, n, coeffs, force, trace)
  of ik4MN5F2GP:
    integrate4MN5F2GP(action, g0, p0, dt, n, coeffs, force, trace)
