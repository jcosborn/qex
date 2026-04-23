import ../[core, scalar, gauge]
from math import sqrt, classify, fcNan, fcInf, fcNegInf

type
  IntegratorKind* = enum
    ik2MN, ik4MN3F1GP, ik4MN5F2GP
  Integrator2MNCoeffs* = object
    lambda*: float
  Integrator4MN3F1GPCoeffs* = object
    lambda*: float
    theta*: float
    chi*: float
  Integrator4MN5F2GPCoeffs* = object
    rho*: float
    theta*: float
    vtheta*: float
    lambda*: float
    xi*: float
  IntegratorCoeffs* = object
    kind*: IntegratorKind
    twoMN*: Integrator2MNCoeffs
    fourMN3F1GP*: Integrator4MN3F1GPCoeffs
    fourMN5F2GP*: Integrator4MN5F2GPCoeffs
  IntegratorStepKind = enum
    iskDrift, iskKick, iskForceGradKick
  IntegratorStep = object
    kind: IntegratorStepKind
    coeff: Gscalar
    gradCoeff: Gscalar
  IntegratorRunSpec = object
    leadDrift, repeatLeadDrift, tailDrift: Gscalar
    steps: seq[IntegratorStep]
    learnedCoeffs: seq[Gscalar]

proc getCoeffValue(coeffs: openArray[float],
                   i: int,
                   defaultValue: float): float =
  if coeffs.len <= i:
    return defaultValue
  coeffs[i]

proc requireCoeffCount(label: string,
                       values: openArray[float],
                       maximum: int) =
  if values.len > maximum:
    raiseValueError(
      label & " accepts at most " & $maximum &
      " coefficient values, got " & $values.len)

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

proc driftStep(coeff: Gscalar): IntegratorStep =
  IntegratorStep(kind: iskDrift, coeff: coeff)

proc kickStep(coeff: Gscalar): IntegratorStep =
  IntegratorStep(kind: iskKick, coeff: coeff)

proc forceGradKickStep(coeff: Gscalar,
                       gradCoeff: Gscalar): IntegratorStep =
  IntegratorStep(kind: iskForceGradKick, coeff: coeff, gradCoeff: gradCoeff)

proc applyIntegratorStep(gc: Gactcoeff,
                         g: var Ggauge,
                         p: var Ggauge,
                         step: IntegratorStep) =
  case step.kind
  of iskDrift:
    g = axexpmuly(step.coeff, p, g)
  of iskKick:
    p = p - step.coeff * gaugeForce(gc, g)
  of iskForceGradKick:
    let fg = gaugeForce(gc, g)
    p = p - step.coeff * gaugeForce(gc, axexpmuly(-step.gradCoeff, fg, g))

proc runIntegrator(gc: Gactcoeff,
                   g0: Ggauge,
                   p0: Ggauge,
                   leadDrift: Gscalar,
                   repeatLeadDrift: Gscalar,
                   tailDrift: Gscalar,
                   n: int,
                   steps: openArray[IntegratorStep]): (Ggauge, Ggauge) =
  var g = g0
  var p = p0
  g = axexpmuly(leadDrift, p, g)
  for i in 0..<n:
    if i > 0:
      g = axexpmuly(repeatLeadDrift, p, g)
    for step in steps:
      applyIntegratorStep(gc, g, p, step)
  g = axexpmuly(tailDrift, p, g)
  (g, p)

proc initIntegratorRunSpec(leadDrift, repeatLeadDrift, tailDrift: Gscalar,
                           steps: openArray[IntegratorStep],
                           learnedCoeffs: openArray[Gscalar]): IntegratorRunSpec =
  result.leadDrift = leadDrift
  result.repeatLeadDrift = repeatLeadDrift
  result.tailDrift = tailDrift
  result.steps = @steps
  result.learnedCoeffs = @learnedCoeffs

proc requireIntegratorStepCount(n: int): int =
  if n <= 0:
    raiseValueError("integrator step count must be >= 1, got " & $n)
  n

proc requireFiniteCoeff(label: string,
                        value: float): float =
  if classify(value) in {fcNan, fcInf, fcNegInf}:
    raiseValueError(
      label & " must be finite after completion, got " & $value)
  value

proc default4MN3F1GPTheta(lambdaValue: float): float =
  ## Default completion for the force-gradient integrator family.
  ## The tuple is all-or-default by design; partial positional completion is unsupported.
  0.5 - 1.0 / sqrt(24.0 * lambdaValue)

proc default4MN3F1GPChi(lambdaValue: float): float =
  let numer = 1.0 - sqrt(6.0 * lambdaValue) * (1.0 - lambdaValue)
  let scale = 20.0 / (1.0 - 2.0 * lambdaValue)
  (numer / 12.0) * scale

proc complete2MNCoeffs(coeffs: openArray[float]): Integrator2MNCoeffs =
  requireCoeffCount("2MN", coeffs, 1)
  result.lambda = coeffs.getCoeffValue(0, 0.1931833275037836)
  result.lambda = requireFiniteCoeff("2MN.lambda", result.lambda)

proc complete4MN3F1GPCoeffs(
    coeffs: openArray[float]): Integrator4MN3F1GPCoeffs =
  requireCoeffCountOrDefault("4MN3F1GP", coeffs, 3)
  if coeffs.len == 0:
    result.lambda = 0.2470939580390842
    result.theta = default4MN3F1GPTheta(result.lambda)
    result.chi = default4MN3F1GPChi(result.lambda)
  else:
    result.lambda = coeffs[0]
    result.theta = coeffs[1]
    result.chi = coeffs[2]
  result.lambda = requireFiniteCoeff("4MN3F1GP.lambda", result.lambda)
  result.theta = requireFiniteCoeff("4MN3F1GP.theta", result.theta)
  result.chi = requireFiniteCoeff("4MN3F1GP.chi", result.chi)

proc default4MN5F2GPXi(lambdaValue: float): float =
  ## Keep the scale factors explicit so learned coefficient formulas are auditable.
  0.0009628905212024874 * (2.0 / lambdaValue * 20.0)

proc complete4MN5F2GPCoeffs(
    coeffs: openArray[float]): Integrator4MN5F2GPCoeffs =
  requireCoeffCountOrDefault("4MN5F2GP", coeffs, 5)
  if coeffs.len == 0:
    result.rho = 0.06419108866816235
    result.theta = 0.1919807940455741
    result.vtheta = 0.1518179640276466
    result.lambda = 0.2158369476787619
    result.xi = default4MN5F2GPXi(result.lambda)
  else:
    result.rho = coeffs[0]
    result.theta = coeffs[1]
    result.vtheta = coeffs[2]
    result.lambda = coeffs[3]
    result.xi = coeffs[4]
  result.rho = requireFiniteCoeff("4MN5F2GP.rho", result.rho)
  result.theta = requireFiniteCoeff("4MN5F2GP.theta", result.theta)
  result.vtheta = requireFiniteCoeff("4MN5F2GP.vtheta", result.vtheta)
  result.lambda = requireFiniteCoeff("4MN5F2GP.lambda", result.lambda)
  result.xi = requireFiniteCoeff("4MN5F2GP.xi", result.xi)

proc parseIntegratorCoeffs*(kind: IntegratorKind,
                            values: openArray[float]): IntegratorCoeffs =
  result.kind = kind
  case kind
  of ik2MN:
    result.twoMN = complete2MNCoeffs(values)
  of ik4MN3F1GP:
    result.fourMN3F1GP = complete4MN3F1GPCoeffs(values)
  of ik4MN5F2GP:
    result.fourMN5F2GP = complete4MN5F2GPCoeffs(values)

proc build2MNSpec(dt: Gscalar,
                  coeffs: Integrator2MNCoeffs): IntegratorRunSpec =
  let lambda = scalarLeafLike(dt, coeffs.lambda)
  let h = 0.5 * dt
  let t05 = lambda * dt
  let t0 = 2.0 * t05
  let t1 = dt - t0
  initIntegratorRunSpec(t05, t0, t05,
    [kickStep(h), driftStep(t1), kickStep(h)],
    [lambda])

proc build4MN3F1GPSpec(dt: Gscalar,
                       coeffs: Integrator4MN3F1GPCoeffs): IntegratorRunSpec =
  let lambda = scalarLeafLike(dt, coeffs.lambda)
  let theta = scalarLeafLike(dt, coeffs.theta)
  let chi = scalarLeafLike(dt, coeffs.chi)
  let a0 = theta * dt
  let a02 = 2.0 * a0
  let a1 = 0.5 * dt - a0
  let b0 = lambda * dt
  let b1 = dt - 2.0 * b0
  let c1 = 0.1 * chi * (dt * dt)
  initIntegratorRunSpec(a0, a02, a0,
    [kickStep(b0), driftStep(a1), forceGradKickStep(b1, c1),
     driftStep(a1), kickStep(b0)],
    [lambda, theta, chi])

proc build4MN5F2GPSpec(dt: Gscalar,
                       coeffs: Integrator4MN5F2GPCoeffs): IntegratorRunSpec =
  let rho = scalarLeafLike(dt, coeffs.rho)
  let theta = scalarLeafLike(dt, coeffs.theta)
  let vtheta = scalarLeafLike(dt, coeffs.vtheta)
  let lambda = scalarLeafLike(dt, coeffs.lambda)
  let xi = scalarLeafLike(dt, coeffs.xi)
  let a0 = rho * dt
  let a02 = 2.0 * a0
  let a1 = theta * dt
  let a2 = (0.5 - (theta + rho)) * dt
  let b1 = lambda * dt
  let b0 = vtheta * dt
  let b2 = (1.0 - 2.0 * (lambda + vtheta)) * dt
  let c1 = 0.05 * xi * (dt * dt)
  initIntegratorRunSpec(a0, a02, a0,
    [kickStep(b0), driftStep(a1), forceGradKickStep(b1, c1),
     driftStep(a2), kickStep(b2), driftStep(a2),
     forceGradKickStep(b1, c1), driftStep(a1), kickStep(b0)],
    [rho, theta, vtheta, lambda, xi])

proc buildIntegratorRunSpec(dt: Gscalar,
                            coeffs: IntegratorCoeffs): IntegratorRunSpec =
  case coeffs.kind
  of ik2MN:
    build2MNSpec(dt, coeffs.twoMN)
  of ik4MN3F1GP:
    build4MN3F1GPSpec(dt, coeffs.fourMN3F1GP)
  of ik4MN5F2GP:
    build4MN5F2GPSpec(dt, coeffs.fourMN5F2GP)

proc integrateGauge*(gc: Gactcoeff,
                     g0: Ggauge,
                     p0: Ggauge,
                     dt: Gscalar,
                     n: int,
                     coeffs: IntegratorCoeffs): (Ggauge, Ggauge, seq[Gscalar]) =
  let spec = buildIntegratorRunSpec(dt, coeffs)
  let (g, p) = runIntegrator(
    gc,
    g0,
    p0,
    spec.leadDrift,
    spec.repeatLeadDrift,
    spec.tailDrift,
    n.requireIntegratorStepCount,
    spec.steps)
  (g, p, spec.learnedCoeffs)
