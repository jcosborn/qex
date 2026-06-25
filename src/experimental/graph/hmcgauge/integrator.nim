import ../[core, scalar, gauge]
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

proc integrate2MN(gc: Gactcoeff,
                  g0: Ggauge,
                  p0: Ggauge,
                  dt: Gscalar,
                  n: int,
                  coeffs: IntegratorCoeffs): IntegrationResult =
  let lambda = toGvalue(dt.runtime, coeffs.lambda)
  let h = 0.5 * dt
  let t05 = lambda * dt
  let t0 = 2.0 * t05
  let t1 = dt - t0
  var g = g0
  var p = p0
  for i in 0..<n:
    g = axexpmuly(if i == 0: t05 else: t0, p, g)
    p = p - h * gaugeForce(gc, g)
    g = axexpmuly(t1, p, g)
    p = p - h * gaugeForce(gc, g)
  g = axexpmuly(t05, p, g)
  IntegrationResult(
    gauge: g,
    momentum: p,
    learnedCoeffs: @[LearnedParameter(name: "lambda", node: lambda)])

proc integrate4MN3F1GP(gc: Gactcoeff,
                       g0: Ggauge,
                       p0: Ggauge,
                       dt: Gscalar,
                       n: int,
                       coeffs: IntegratorCoeffs): IntegrationResult =
  let lambda = toGvalue(dt.runtime, coeffs.lambda)
  let theta = toGvalue(dt.runtime, coeffs.theta)
  let chi = toGvalue(dt.runtime, coeffs.chi)
  let a0 = theta * dt
  let a02 = 2.0 * a0
  let a1 = 0.5 * dt - a0
  let b0 = lambda * dt
  let b1 = dt - 2.0 * b0
  let c1 = 0.1 * chi * (dt * dt)
  var g = g0
  var p = p0
  for i in 0..<n:
    g = axexpmuly(if i == 0: a0 else: a02, p, g)
    p = p - b0 * gaugeForce(gc, g)
    g = axexpmuly(a1, p, g)
    let fg = gaugeForce(gc, g)
    p = p - b1 * gaugeForce(gc, axexpmuly(-c1, fg, g))
    g = axexpmuly(a1, p, g)
    p = p - b0 * gaugeForce(gc, g)
  g = axexpmuly(a0, p, g)
  IntegrationResult(
    gauge: g,
    momentum: p,
    learnedCoeffs: @[
      LearnedParameter(name: "lambda", node: lambda),
      LearnedParameter(name: "theta", node: theta),
      LearnedParameter(name: "chi", node: chi)])

proc integrate4MN5F2GP(gc: Gactcoeff,
                       g0: Ggauge,
                       p0: Ggauge,
                       dt: Gscalar,
                       n: int,
                       coeffs: IntegratorCoeffs): IntegrationResult =
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
  var g = g0
  var p = p0
  for i in 0..<n:
    g = axexpmuly(if i == 0: a0 else: a02, p, g)
    p = p - b0 * gaugeForce(gc, g)
    g = axexpmuly(a1, p, g)
    block:
      let fg = gaugeForce(gc, g)
      p = p - b1 * gaugeForce(gc, axexpmuly(-c1, fg, g))
    g = axexpmuly(a2, p, g)
    p = p - b2 * gaugeForce(gc, g)
    g = axexpmuly(a2, p, g)
    block:
      let fg = gaugeForce(gc, g)
      p = p - b1 * gaugeForce(gc, axexpmuly(-c1, fg, g))
    g = axexpmuly(a1, p, g)
    p = p - b0 * gaugeForce(gc, g)
  g = axexpmuly(a0, p, g)
  IntegrationResult(
    gauge: g,
    momentum: p,
    learnedCoeffs: @[
      LearnedParameter(name: "rho", node: rho),
      LearnedParameter(name: "theta", node: theta),
      LearnedParameter(name: "vtheta", node: vtheta),
      LearnedParameter(name: "lambda", node: lambda),
      LearnedParameter(name: "xi", node: xi)])

proc integrateGauge*(gc: Gactcoeff,
                     g0: Ggauge,
                     p0: Ggauge,
                     dt: Gscalar,
                     n: int,
                     coeffs: IntegratorCoeffs): IntegrationResult =
  if n <= 0:
    raiseValueError("integrator step count must be >= 1, got " & $n)
  case coeffs.kind
  of ik2MN:
    integrate2MN(gc, g0, p0, dt, n, coeffs)
  of ik4MN3F1GP:
    integrate4MN3F1GP(gc, g0, p0, dt, n, coeffs)
  of ik4MN5F2GP:
    integrate4MN5F2GP(gc, g0, p0, dt, n, coeffs)
