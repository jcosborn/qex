import ../[core, scalar, gauge]
import layout, physics/qcdTypes   # threads / simd reductions for forceRmsMinMax
from math import sqrt

type
  IntegratorKind* = enum
    ik2MN, ik4MN3F1GP, ik4MN5F2GP
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
  IntegrationResult* = object
    gauge*: Ggauge
    momentum*: Ggauge
    learnedCoeffs*: seq[LearnedParameter]
    forces*: seq[Ggauge]
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

proc mdForceStats*(forces: openArray[Ggauge]): MdForceStats =
  ## Statistics from retained integrator force nodes.
  if forces.len == 0:
    raiseValueError("MD force statistics require at least one force")
  discard forces[^1].eval
  let dof = float(forces[0].gval.len * forces[0].gval[0].l.physVol)
  result.count = forces.len
  result.fminMin = 1e300
  for force in forces:
    let f = forceRmsMinMaxValue(force, dof)
    result.rmsMean += f.rms
    result.fminMean += f.fmin
    result.fmaxMean += f.fmax
    result.rmsMax = max(result.rmsMax, f.rms)
    result.fminMin = min(result.fminMin, f.fmin)
    result.fmaxMax = max(result.fmaxMax, f.fmax)
  result.rmsMean /= float(result.count)
  result.fminMean /= float(result.count)
  result.fmaxMean /= float(result.count)

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
  of "4MN3F1GP":
    ik4MN3F1GP
  of "4MN5F2GP":
    ik4MN5F2GP
  else:
    raiseValueError("unknown intalg: " & name)

proc parseIntegratorCoeffs*(kind: IntegratorKind,
                            values: openArray[float]): IntegratorCoeffs =
  case kind
  of ik2MN:
    requireCoeffCountOrDefault("2MN", values, 1)
    result = IntegratorCoeffs(kind: ik2MN)
    result.lambda =
      if values.len == 0:
        0.1931833275037836
      else:
        values[0]
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

proc integrate2MN(action: GaugeAction,
                  g0: Ggauge,
                  p0: Ggauge,
                  dt: Gscalar,
                  n: int,
                  coeffs: IntegratorCoeffs,
                  force: GaugeForceFn): IntegrationResult =
  let lambda = toGvalue(dt.runtime, coeffs.lambda)
  let h = 0.5 * dt
  let t05 = lambda * dt
  let t0 = 2.0 * t05
  let t1 = dt - t0
  let mh = -h
  var g = g0
  var p = p0
  var forces: seq[Ggauge]
  for i in 0..<n:
    g = axexpmuly(if i == 0: t05 else: t0, p, g)
    p = axpy(mh, gradForce(action, g, forces, force), p)
    g = axexpmuly(t1, p, g)
    p = axpy(mh, gradForce(action, g, forces, force), p)
  g = axexpmuly(t05, p, g)
  IntegrationResult(
    gauge: g,
    momentum: p,
    learnedCoeffs: @[LearnedParameter(name: "lambda", node: lambda)],
    forces: forces)

proc integrate4MN3F1GP(action: GaugeAction,
                       g0: Ggauge,
                       p0: Ggauge,
                       dt: Gscalar,
                       n: int,
                       coeffs: IntegratorCoeffs,
                       force: GaugeForceFn): IntegrationResult =
  let lambda = toGvalue(dt.runtime, coeffs.lambda)
  let theta = toGvalue(dt.runtime, coeffs.theta)
  let chi = toGvalue(dt.runtime, coeffs.chi)
  let a0 = theta * dt
  let a02 = 2.0 * a0
  let a1 = 0.5 * dt - a0
  let b0 = lambda * dt
  let b1 = dt - 2.0 * b0
  let c1 = 0.1 * chi * (dt * dt)
  let mb0 = -b0
  let mb1 = -b1
  var g = g0
  var p = p0
  var forces: seq[Ggauge]
  for i in 0..<n:
    g = axexpmuly(if i == 0: a0 else: a02, p, g)
    p = axpy(mb0, gradForce(action, g, forces, force), p)
    g = axexpmuly(a1, p, g)
    let fg = gradForce(action, g, forces, force)
    p = axpy(mb1, gradForce(action, axexpmuly(-c1, fg, g), forces, force), p)
    g = axexpmuly(a1, p, g)
    p = axpy(mb0, gradForce(action, g, forces, force), p)
  g = axexpmuly(a0, p, g)
  IntegrationResult(
    gauge: g,
    momentum: p,
    learnedCoeffs: @[
      LearnedParameter(name: "lambda", node: lambda),
      LearnedParameter(name: "theta", node: theta),
      LearnedParameter(name: "chi", node: chi)],
    forces: forces)

proc integrate4MN5F2GP(action: GaugeAction,
                       g0: Ggauge,
                       p0: Ggauge,
                       dt: Gscalar,
                       n: int,
                       coeffs: IntegratorCoeffs,
                       force: GaugeForceFn): IntegrationResult =
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
  let mb0 = -b0
  let mb1 = -b1
  let mb2 = -b2
  var g = g0
  var p = p0
  var forces: seq[Ggauge]
  for i in 0..<n:
    g = axexpmuly(if i == 0: a0 else: a02, p, g)
    p = axpy(mb0, gradForce(action, g, forces, force), p)
    g = axexpmuly(a1, p, g)
    block:
      let fg = gradForce(action, g, forces, force)
      p = axpy(mb1, gradForce(action, axexpmuly(-c1, fg, g), forces, force), p)
    g = axexpmuly(a2, p, g)
    p = axpy(mb2, gradForce(action, g, forces, force), p)
    g = axexpmuly(a2, p, g)
    block:
      let fg = gradForce(action, g, forces, force)
      p = axpy(mb1, gradForce(action, axexpmuly(-c1, fg, g), forces, force), p)
    g = axexpmuly(a1, p, g)
    p = axpy(mb0, gradForce(action, g, forces, force), p)
  g = axexpmuly(a0, p, g)
  IntegrationResult(
    gauge: g,
    momentum: p,
    learnedCoeffs: @[
      LearnedParameter(name: "rho", node: rho),
      LearnedParameter(name: "theta", node: theta),
      LearnedParameter(name: "vtheta", node: vtheta),
      LearnedParameter(name: "lambda", node: lambda),
      LearnedParameter(name: "xi", node: xi)],
    forces: forces)

proc integrateGauge*(action: GaugeAction,
                     g0: Ggauge,
                     p0: Ggauge,
                     dt: Gscalar,
                     n: int,
                     coeffs: IntegratorCoeffs,
                     force: GaugeForceFn = nil): IntegrationResult =
  if n <= 0:
    raiseValueError("integrator step count must be >= 1, got " & $n)
  case coeffs.kind
  of ik2MN:
    integrate2MN(action, g0, p0, dt, n, coeffs, force)
  of ik4MN3F1GP:
    integrate4MN3F1GP(action, g0, p0, dt, n, coeffs, force)
  of ik4MN5F2GP:
    integrate4MN5F2GP(action, g0, p0, dt, n, coeffs, force)
