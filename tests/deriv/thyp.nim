import qex, physics/qcdTypes, algorithms/numdiff, gauge/hypsmear, gauge/gaugeAction, gauge/smearutil
import maths/matrixFunctions
import ./utils

qexInit()

var env = initGaugeTestEnv()
let lo = env.lo
var g = env.g
var r = env.r

type
  GaugeField = typeof(g)

proc makeChain(): GaugeField =
  var c = lo.newGauge
  c.randomTAH r
  c

proc makeDir(): GaugeField =
  var d = lo.newGauge
  d.randomTAH r
  d

proc makeProbe(): GaugeField =
  var p = lo.newGauge
  p.gaussian r
  p

var coef = HypCoefs(alpha1:0.4, alpha2:0.5, alpha3:0.5)
var info = PerfInfo()

let gc = GaugeActionCoeffs(plaq:6.0)

var shapeHs = lo.newHypSmear(coef)
var shapeSg = lo.newGauge
shapeHs.smear(g, shapeSg)

type
  PairField = type(newOneOf(shapeHs.state.l1x))

let adjCoefCombos = [
  HypCoefs(alpha1:0.0, alpha2:0.0, alpha3:0.0),  # 0: All off
  HypCoefs(alpha1:0.4, alpha2:0.0, alpha3:0.0),  # 1: Only L1
  HypCoefs(alpha1:0.0, alpha2:0.5, alpha3:0.0),  # 2: Only L2
  HypCoefs(alpha1:0.0, alpha2:0.0, alpha3:0.5),  # 3: Only L3
  HypCoefs(alpha1:0.4, alpha2:0.5, alpha3:0.0),  # 4: L1+L2
  HypCoefs(alpha1:0.4, alpha2:0.0, alpha3:0.5),  # 5: L1+L3
  HypCoefs(alpha1:0.0, alpha2:0.5, alpha3:0.5),  # 6: L2+L3
  HypCoefs(alpha1:0.4, alpha2:0.5, alpha3:0.5),  # 7: All on
]

let adjCoefNames = [
  "All off (trivial)",
  "Only L1",
  "Only L2",
  "Only L3",
  "L1+L2",
  "L1+L3",
  "L2+L3",
  "All on",
]

echo "S(g): ", gc.gaugeAction1(g)

# Production hypsmear using closure-based derivative

proc hypAction(g:auto):auto =
  var sg = lo.newGauge
  coef.smear(g, sg)
  gc.gaugeAction1(sg)

proc hypForce(g:auto, f:auto) =
  var sg = lo.newGauge
  var ds = lo.newGauge
  var fn = coef.smearGetForce(g, sg, info)
  gc.gaugeActionDeriv(sg, ds)
  fn(f, ds)
  contractProjectTAH(g, f)

# Preallocated hypsmear environment derivative

proc hypActionPrealloc(g:auto):auto =
  var sg = lo.newGauge
  var hs = lo.newHypSmear(coef)
  hs.smear(g, sg)
  gc.gaugeAction1(sg)

proc hypForcePrealloc(g:auto, f:auto) =
  var sg = lo.newGauge
  var ds = lo.newGauge
  var hs = lo.newHypSmear(coef)
  hs.smear(g, sg)
  gc.gaugeActionDeriv(sg, ds)
  hs.smearVJP(f, ds)
  contractProjectTAH(g, f)

proc testSmearForceConsistency(): int =
  let cfg = initTestConfig(
    name = "smearGetForce/prealloc consistency",
    samples = 1,
    atol = 1e-14,
    rtol = 1e-12,
    verbose = true
  )
  cfg.reportBanner
  var sgClosure = lo.newGauge
  var sgPrealloc = lo.newGauge
  var chain = lo.newGauge
  var fClosure = lo.newGauge
  var fPrealloc = lo.newGauge
  var diff = lo.newGauge
  var hs = lo.newHypSmear(coef)
  var fn = coef.smearGetForce(g, sgClosure, info)

  hs.smear(g, sgPrealloc)
  gc.gaugeActionDeriv(sgClosure, chain)
  fn(fClosure, chain)
  gc.gaugeActionDeriv(sgPrealloc, chain)
  hs.smearVJP(fPrealloc, chain)
  contractProjectTAH(g, fClosure)
  contractProjectTAH(g, fPrealloc)

  threads:
    for mu in 0..<4:
      diff[mu] := sgPrealloc[mu] - sgClosure[mu]
  let smearRel = sqrt(norm2Like(diff) / max(norm2Like(sgClosure), 1e-30))

  threads:
    for mu in 0..<4:
      diff[mu] := fPrealloc[mu] - fClosure[mu]
  let forceRel = sqrt(norm2Like(diff) / max(norm2Like(fClosure), 1e-30))

  let smearOk = reportCompare(cfg, 0, smearRel, 0.0, "smear output rel diff")
  let forceOk = reportCompare(cfg, 0, forceRel, 0.0, "force rel diff")
  if not (smearOk and forceOk):
    result = 1

proc restoreHypState(hs: var auto; sg: GaugeField) =
  hs.smear(g, sg)

proc primeHypReverse(hs: var auto; sg, deriv: GaugeField; chain: GaugeField) =
  restoreHypState(hs, sg)
  hs.smearVJP(deriv, chain)

template forPairsTest(mu, nu: untyped; body: untyped) =
  for mu in 0..<4:
    for nu in 0..<4:
      if nu != mu:
        body

template forThirdTest(mu, nu, a, b: untyped; body: untyped) =
  for a in 0..<4:
    if a != mu and a != nu:
      let b = 1 + 2 + 3 - mu - nu - a
      discard b
      body

proc projectPairField(xfield: PairField): PairField =
  var projected = newOneOf(xfield)
  threads:
    forPairsTest(mu, nu):
      projectU(projected[mu,nu], xfield[mu,nu])
  projected

proc currentHypL1(hs: var auto): PairField =
  when compiles(hs.state.l1):
    result = hs.state.l1
  else:
    result = projectPairField(hs.state.l1x)

proc currentHypL2(hs: var auto): PairField =
  when compiles(hs.state.l2):
    result = hs.state.l2
  else:
    result = projectPairField(hs.state.l2x)

proc makePairFieldGaussianLike(shape: PairField): PairField =
  var field = newOneOf(shape)
  threads:
    forPairsTest(mu, nu):
      field[mu,nu].gaussian r
  field

proc perturbGaugeFieldAdd(x: GaugeField; dx: GaugeField; eps: float; res: var GaugeField) =
  let resPtr = cast[ptr GaugeField](unsafeAddr(res))
  threads:
    for mu in 0..<4:
      resPtr[][mu] := x[mu]
      resPtr[][mu] += eps * dx[mu]

proc perturbPairFieldAdd(x: PairField; dx: PairField; eps: float; res: var PairField) =
  let resPtr = cast[ptr PairField](unsafeAddr(res))
  threads:
    forPairsTest(mu, nu):
      resPtr[][mu,nu] := x[mu,nu]
      resPtr[][mu,nu] += eps * dx[mu,nu]

proc zeroGaugeField(): GaugeField =
  var z = lo.newGauge
  threads:
    for mu in 0..<4:
      z[mu] := 0
  z

proc zeroPairField(shape: PairField): PairField =
  var z = newOneOf(shape)
  threads:
    forPairsTest(mu, nu):
      z[mu,nu] := 0
  z

proc zeroInPlace(x: var GaugeField) =
  let xPtr = cast[ptr GaugeField](unsafeAddr(x))
  threads:
    for mu in 0..<4:
      xPtr[][mu] := 0

proc zeroInPlace(x: var PairField) =
  let xPtr = cast[ptr PairField](unsafeAddr(x))
  threads:
    forPairsTest(mu, nu):
      xPtr[][mu,nu] := 0

proc initShiftedGaugeTangents(dx: GaugeField): auto =
  var s: array[4,array[4,Shifter[type(dx[0]),type(dx[0][0])]]]
  forPairsTest(mu, nu):
    s[mu][nu] = newShifter(dx[mu], nu, 1)
  threads:
    forPairsTest(mu, nu):
      discard s[mu][nu] ^*! dx[mu]
  s

proc initShiftedStage2Tangents(dl1: PairField): auto =
  var s: array[4, array[4, array[4, array[2, Shifter[type(dl1[0,1]),type(dl1[0,1][0])]]]]]
  forPairsTest(mu, nu):
    forThirdTest(mu, nu, a, b):
      s[mu][nu][a][0] = newShifter(dl1[a,b], mu, 1)
      s[mu][nu][a][1] = newShifter(dl1[mu,b], a, 1)
  threads:
    forPairsTest(mu, nu):
      forThirdTest(mu, nu, a, b):
        discard s[mu][nu][a][0] ^*! dl1[a,b]
        discard s[mu][nu][a][1] ^*! dl1[mu,b]
  s

proc initShiftedOutputTangents(dl2: PairField): auto =
  var s: array[4,array[4,Shifter[type(dl2[0,1]),type(dl2[0,1][0])]]]
  forPairsTest(mu, nu):
    s[mu][nu] = newShifter(dl2[mu,nu], nu, 1)
  threads:
    forPairsTest(mu, nu):
      discard s[mu][nu] ^*! dl2[mu,nu]
  s

proc runMixedAdjointIdentity(
    cfg: TestConfig,
    makeBase: auto,
    makeDir: auto,
    makeYbar: auto,
    jvp: auto,
    vjp: auto
  ): TestResult =
  type
    X = typeof(makeBase(0))
    DX = typeof(makeDir(0, default(X)))
    Y = typeof(makeYbar(0, default(X)))
  reportBanner(cfg)
  result.name = cfg.name
  for sample in 0..<cfg.samples:
    let base = makeBase(sample)
    let dx = makeDir(sample, base)
    let ybar = makeYbar(sample, base)
    var dy = allocLike(ybar)
    var dxbar = allocLike(dx)
    jvp(base, dx, dy)
    vjp(base, ybar, dxbar)
    let lhs = dotLike(ybar, dy)
    let rhs = dotLike(dxbar, dx)
    let err = abs(lhs - rhs)
    let refVal = max(1.0, max(abs(lhs), abs(rhs)))
    let rel = err / refVal
    let ok = err <= max(cfg.atol, cfg.rtol * refVal)
    reportSimple(cfg, sample, lhs, rhs, err, rel, ok)
    if ok: inc result.passed else: inc result.failed

template contractedForceJvp(hs, x, forcebar, chain, dg: untyped;
                            dchain: untyped = @[]): untyped =
  ## Compose the clean uncontracted `smearHVP` with the outer
  ## `contractProjectTAH(x, ·)` derivative.
  block:
    var sg = lo.newGauge
    hs.smear(x, sg)

    var deriv = lo.newGauge
    hs.smearVJP(deriv, chain)

    # `smearVJP` reuses reverse scratch, so rebuild the cached forward state
    # before evaluating the uncontracted HVP.
    hs.smear(x, sg)

    var dderiv = lo.newGauge
    hs.smearHVP(dderiv, chain, dg, dchain)

    var derivbar = lo.newGauge
    var gbarContract = lo.newGauge
    contractProjectTAHVJP(x, deriv, forcebar, derivbar, gbarContract)

    var dgEff = lo.newGauge
    threads:
      for mu in 0..<4:
        dgEff[mu] := dg[mu] * x[mu]

    dotLike(derivbar, dderiv) + dotLike(gbarContract, dgEff)

proc runSmearHVPCase(name: string; samples: int): int =
  var hs = lo.newHypSmear(coef)
  var sg = lo.newGauge
  hs.smear(g, sg)
  let chain = makeChain()
  let cfg = initTestConfig(
    name = name,
    samples = samples,
    eps = 1e-2,
    atol = 1e-6,
    rtol = 1e-6,
    verbose = true
  )
  let res = runGaugeHVPFD(
    env,
    grad = proc(x: GaugeField; res: var GaugeField) =
      var sg2 = lo.newGauge
      hs.smear(x, sg2)
      hs.smearVJP(res, chain),
    hvp = proc(x: GaugeField; dx: GaugeField; res: var GaugeField) =
      hs.smear(x, sg)
      hs.smearHVP(res, chain, dx),
    cfg = cfg
  )
  res.failed

proc testSmearDeriv2(): int =
  result = runSmearHVPCase(
    name = "smearHVP (uncontracted HYP second derivative)",
    samples = 5
  )

## Level isolation tests for debugging smearHVP
## Test with different alpha coefficients to isolate which level has the bug
proc testSmearDeriv2Level(levelName: string, testCoef: HypCoefs): int =
  var hsTest = lo.newHypSmear(testCoef)
  var sgTest = lo.newGauge
  hsTest.smear(g, sgTest)
  let chainTest = makeChain()
  let cfg = initTestConfig(
    name = "smearHVP: " & levelName,
    samples = 3,
    eps = 1e-2,
    scale = 2.0,
    atol = 1e-6,
    rtol = 1e-6,
    fdFactor = 3.0,
    verbose = true
  )
  let res = runGaugeHVPFD(
    env,
    grad = proc(x: GaugeField; res: var GaugeField) =
      var sg2 = lo.newGauge
      hsTest.smear(x, sg2)
      hsTest.smearVJP(res, chainTest),
    hvp = proc(x: GaugeField; dx: GaugeField; res: var GaugeField) =
      hsTest.smear(x, sgTest)
      hsTest.smearHVP(res, chainTest, dx),
    cfg = cfg
  )
  result = res.failed

proc testGradForceFullChain(): int =
  var hs = lo.newHypSmear(coef)
  var sg = lo.newGauge
  hs.smear(g, sg)

  var p = makeDir()

  let cfg = initTestConfig(
    name = "grad <p, force(g)> (full chain)",
    samples = 5,
    eps = 1e-2,
    scale = 2.0,
    atol = 1e-5,
    rtol = 1e-5,
    fdFactor = 3.0,
    verbose = true
  )

  let res = runFDScalarJVPDefault(
    cfg = cfg,
    makeBase = proc(sample: int): GaugeField =
      g,
    makeDir = proc(sample: int; x: GaugeField): GaugeField =
      var dg = lo.newGauge
      dg.randomTAH r
      dg,
    f = proc(x: GaugeField): float =
      var hsx = lo.newHypSmear(coef)
      var sgx = lo.newGauge
      hsx.smear(x, sgx)
      var chain = lo.newGauge
      gc.gaugeActionDeriv(sgx, chain)
      var force = lo.newGauge
      hsx.smearVJP(force, chain)
      contractProjectTAH(x, force)
      var s = 0.0
      for mu in 0..<4:
        s += redot(p[mu], force[mu])
      s,
    jvp = proc(x: GaugeField; dx: GaugeField): float =
      hs.smear(x, sg)
      var chain0 = lo.newGauge
      gc.gaugeActionDeriv(sg, chain0)
      var dsg = lo.newGauge
      hs.smearJVP(dsg, dx, sg)
      var dchain = lo.newGauge
      threads:
        for mu in 0..<4:
          dchain[mu] := 0
      gc.gaugeDerivDeriv2(sg, dsg, dchain)
      threads:
        for mu in 0..<4:
          dchain[mu] := -1.0 * dchain[mu]
      contractedForceJvp(hs, x, p, chain0, dx, dchain)
  )
  result = res.failed

## Test: smearJVP alone (compare against FD of smear output)

proc testSmearTangent(): int =
  var hs = lo.newHypSmear(coef)
  var sg = lo.newGauge
  let mkProbe = makeProbe
  let cfg = initTestConfig(
    name = "smearJVP (forward tangent of HYP smearing)",
    samples = 3,
    eps = 1e-2,
    scale = 2.0,
    atol = 1e-5,
    rtol = 1e-5,
    fdFactor = 3.0,
    verbose = true
  )
  let res = runGaugeJVPFD(
    env = env,
    f = proc(x: GaugeField; res: var GaugeField) =
      hs.smear(x, res),
    jvp = proc(x: GaugeField; dx: GaugeField; res: var GaugeField) =
      hs.smear(x, sg)
      hs.smearJVP(res, dx, sg),
    cfg = cfg,
    makeProbe = proc(sample: int; x: GaugeField): GaugeField =
      mkProbe()
  )
  result = res.failed

## Test: gaugeDerivDeriv2 alone (verify it computes d/dε[gaugeActionDeriv(g+ε·h)])

proc testGaugeDerivDeriv2(): int =
  var probe = makeProbe()
  let mkProbe = makeProbe
  let cfg = initTestConfig(
    name = "gaugeDerivDeriv2 (second derivative of gauge action)",
    samples = 3,
    eps = 1e-2,
    scale = 2.0,
    atol = 1e-5,
    rtol = 1e-4,
    fdFactor = 3.0,
    verbose = true
  )
  let res = runGaugeJVPFD(
    env = env,
    f = proc(x: GaugeField; res: var GaugeField) =
      gc.gaugeActionDeriv(x, res),
    jvp = proc(x: GaugeField; dx: GaugeField; res: var GaugeField) =
      var dgEff = lo.newGauge
      let resPtr = cast[ptr GaugeField](unsafeAddr(res))
      threads:
        for mu in 0..<4:
          dgEff[mu] := dx[mu] * x[mu]
          resPtr[][mu] := 0
      gc.gaugeDerivDeriv2(x, dgEff, res)
      threads:
        for mu in 0..<4:
          resPtr[][mu] := -resPtr[][mu],
    cfg = cfg,
    makeProbe = proc(sample: int; x: GaugeField): GaugeField =
      mkProbe()
  )
  result = res.failed

## Test: Fixed-chain contracted force directional derivative.
## This keeps a scalar sanity check for the outer contraction while the HYP API
## itself stays uncontracted.

proc testSmearDeriv2FixedChain(): int =
  var hs = lo.newHypSmear(coef)
  let chain = makeChain()

  let p = makeDir()

  let cfg = initTestConfig(
    name = "smearHVP fixed chain (sanity check)",
    samples = 3,
    eps = 1e-2,
    scale = 2.0,
    atol = 1e-5,
    rtol = 1e-5,
    fdFactor = 3.0,
    verbose = true
  )
  let res = runFDScalarJVPDefault(
    cfg = cfg,
    makeBase = proc(sample: int): GaugeField =
      g,
    makeDir = proc(sample: int; x: GaugeField): GaugeField =
      var dg = lo.newGauge
      dg.randomTAH r
      dg,
    f = proc(x: GaugeField): float =
      var hsx = lo.newHypSmear(coef)
      var sgx = lo.newGauge
      hsx.smear(x, sgx)
      var force = lo.newGauge
      hsx.smearVJP(force, chain)
      contractProjectTAH(x, force)
      var s = 0.0
      for mu in 0..<4:
        s += redot(p[mu], force[mu])
      s,
    jvp = proc(x: GaugeField; dx: GaugeField): float =
      contractedForceJvp(hs, x, p, chain, dx)
  )
  result = res.failed

## Test: Adjoint identity for smearJVP/smearVJP
## Verifies: <probe, smearJVP(dg)> = <smearVJP(probe), dg*g>

proc testSmearTangentAdjoint(): int =
  let cfg = initTestConfig(
    name = "smearJVP/smearVJP adjoint identity",
    samples = 3,
    atol = 1e-12,
    rtol = 1e-10,
    verbose = true
  )
  var hs = lo.newHypSmear(coef)
  var sg = lo.newGauge
  let res = runAdjointIdentity(
    cfg = cfg,
    makeBase = proc(sample: int): GaugeField =
      g,
    makeDir = proc(sample: int; x: GaugeField): GaugeField =
      var dg = lo.newGauge
      dg.randomTAH r
      dg,
    jvp = proc(x: GaugeField; dx: GaugeField; res: var GaugeField) =
      var dg = lo.newGauge
      for mu in 0..<4:
        dg[mu] := dx[mu] * x[mu].adj
      hs.smear(x, sg)
      hs.smearJVP(res, dg, sg),
    vjp = proc(x: GaugeField; ybar: GaugeField; res: var GaugeField) =
      hs.smear(x, sg)
      hs.smearVJP(res, ybar)
  )
  result = res.failed

## Test: Adjoint identity for smearVJP w.r.t. chain
## Verifies: <derivbar, smearVJP(g, chain)> = <smearVJPChain(derivbar), chain>

proc testSmearDerivAdjChain(): int =
  let cfg = initTestConfig(
    name = "smearVJPChain adjoint identity",
    samples = 3,
    atol = 1e-12,
    rtol = 1e-10,
    verbose = true
  )
  var hs = lo.newHypSmear(coef)
  var sg = lo.newGauge
  hs.smear(g, sg)
  let res = runGaugeAdjointIdentity(
    env,
    jvp = proc(x: GaugeField; dx: GaugeField; res: var GaugeField) =
      hs.smear(x, sg)  # ensure internal state is current
      hs.smearVJP(res, dx),
    vjp = proc(x: GaugeField; ybar: GaugeField; res: var GaugeField) =
      hs.smear(x, sg)
      hs.smearVJPChain(res, ybar),
    cfg = cfg,
    makeDir = proc(sample: int; x: GaugeField): GaugeField =
      var chain = lo.newGauge
      chain.gaussian r
      chain,
  )
  result = res.failed

proc testHypPairStage1JVP(): int =
  var hs = lo.newHypSmear(coef)
  let cfg = initTestConfig(
    name = "L1 stage JVP",
    samples = 3,
    eps = 1e-3,
    scale = 2.0,
    atol = 1e-5,
    rtol = 1e-4,
    fdFactor = 3.0,
    verbose = true
  )
  let res = runFDJVP(
    cfg = cfg,
    makeBase = proc(sample: int): GaugeField =
      g,
    makeDir = proc(sample: int; x: GaugeField): GaugeField =
      var dx = lo.newGauge
      dx.gaussian r
      dx,
    makeProbe = proc(sample: int; x: GaugeField): PairField =
      var hsProbe = lo.newHypSmear(coef)
      var sgProbe = lo.newGauge
      hsProbe.smear(x, sgProbe)
      makePairFieldGaussianLike(currentHypL1(hsProbe)),
    perturb = proc(x: GaugeField; dx: GaugeField; eps: float; res: var GaugeField) =
      perturbGaugeFieldAdd(x, dx, eps, res),
    f = proc(x: GaugeField; res: var PairField) =
      var sg = lo.newGauge
      hs.smear(x, sg)
      var pre = newOneOf(hs.state.l1x)
      hypPairStageL1(hs, pre, res),
    jvp = proc(x: GaugeField; dx: GaugeField; res: var PairField) =
      var sg = lo.newGauge
      hs.smear(x, sg)
      var pre0 = newOneOf(hs.state.l1x)
      var proj0 = newOneOf(hs.state.l1x)
      hypPairStageL1(hs, pre0, proj0)
      let ds1 = initShiftedGaugeTangents(dx)
      var dpre = newOneOf(pre0)
      hypPairStageJVPL1(hs, dx, ds1, pre0, proj0, dx, dpre, res)
  )
  result = res.failed

proc testHypPairStage2JVP(): int =
  var hs = lo.newHypSmear(coef)
  var sg = lo.newGauge
  hs.smear(g, sg)
  let l1 = currentHypL1(hs)
  let l2 = currentHypL2(hs)
  let dgEffZero = zeroGaugeField()
  let cfg = initTestConfig(
    name = "L2 stage JVP (pair input)",
    samples = 3,
    eps = 1e-3,
    scale = 2.0,
    atol = 1e-5,
    rtol = 1e-4,
    fdFactor = 3.0,
    verbose = true
  )
  let res = runFDJVP(
    cfg = cfg,
    makeBase = proc(sample: int): PairField =
      l1,
    makeDir = proc(sample: int; x: PairField): PairField =
      makePairFieldGaussianLike(x),
    makeProbe = proc(sample: int; x: PairField): PairField =
      makePairFieldGaussianLike(l2),
    perturb = proc(x: PairField; dx: PairField; eps: float; res: var PairField) =
      perturbPairFieldAdd(x, dx, eps, res),
    f = proc(x: PairField; res: var PairField) =
      var pre = newOneOf(hs.state.l2x)
      hypPairStageL2(hs, x, pre, res),
    jvp = proc(x: PairField; dx: PairField; res: var PairField) =
      var pre0 = newOneOf(hs.state.l2x)
      var proj0 = newOneOf(hs.state.l2x)
      hypPairStageL2(hs, x, pre0, proj0)
      let dsl1 = initShiftedStage2Tangents(dx)
      var dpre = newOneOf(pre0)
      hypPairStageJVPL2(hs, x, dx, dsl1, pre0, proj0, dgEffZero, dpre, res)
  )
  result = res.failed

proc testHypOutputStageJVP(): int =
  var hs = lo.newHypSmear(coef)
  var sg = lo.newGauge
  hs.smear(g, sg)
  let l2 = currentHypL2(hs)
  let dgEffZero = zeroGaugeField()
  let cfg = initTestConfig(
    name = "Output stage JVP (pair input)",
    samples = 3,
    eps = 1e-3,
    scale = 2.0,
    atol = 1e-5,
    rtol = 1e-4,
    fdFactor = 3.0,
    verbose = true
  )
  let res = runFDJVP(
    cfg = cfg,
    makeBase = proc(sample: int): PairField =
      l2,
    makeDir = proc(sample: int; x: PairField): PairField =
      makePairFieldGaussianLike(x),
    makeProbe = proc(sample: int; x: PairField): GaugeField =
      makeProbe(),
    perturb = proc(x: PairField; dx: PairField; eps: float; res: var PairField) =
      perturbPairFieldAdd(x, dx, eps, res),
    f = proc(x: PairField; res: var GaugeField) =
      var proj = lo.newGauge
      hypOutputStage(hs, x, res, proj),
    jvp = proc(x: PairField; dx: PairField; res: var GaugeField) =
      let dsl2 = initShiftedOutputTangents(dx)
      hypOutputStageJVP(hs, x, dx, dsl2, dgEffZero, res)
  )
  result = res.failed

proc testHypPairStage1ChainAdjoint(): int =
  var hs = lo.newHypSmear(coef)
  let cfg = initTestConfig(
    name = "L1 stage staple/chain adjoint",
    samples = 3,
    atol = 1e-12,
    rtol = 1e-10,
    verbose = true
  )
  let res = runMixedAdjointIdentity(
    cfg = cfg,
    makeBase = proc(sample: int): GaugeField =
      g,
    makeDir = proc(sample: int; x: GaugeField): PairField =
      var hsProbe = lo.newHypSmear(coef)
      var sgProbe = lo.newGauge
      hsProbe.smear(x, sgProbe)
      makePairFieldGaussianLike(currentHypL1(hsProbe)),
    makeYbar = proc(sample: int; x: GaugeField): GaugeField =
      makeProbe(),
    jvp = proc(x: GaugeField; dx: PairField; res: var GaugeField) =
      var sg = lo.newGauge
      hs.smear(x, sg)
      zeroInPlace(res)
      hypPairStageStapleVJPL1(hs, dx, res),
    vjp = proc(x: GaugeField; ybar: GaugeField; res: var PairField) =
      var sg = lo.newGauge
      hs.smear(x, sg)
      zeroInPlace(res)
      hypPairStageChainVJPL1(hs, ybar, res)
  )
  result = res.failed

proc testHypPairStage2ChainAdjoint(): int =
  var hs = lo.newHypSmear(coef)
  var sg = lo.newGauge
  hs.smear(g, sg)
  let l1 = currentHypL1(hs)
  let l2 = currentHypL2(hs)
  let cfg = initTestConfig(
    name = "L2 stage staple/chain adjoint",
    samples = 3,
    atol = 1e-12,
    rtol = 1e-10,
    verbose = true
  )
  let res = runMixedAdjointIdentity(
    cfg = cfg,
    makeBase = proc(sample: int): int =
      0,
    makeDir = proc(sample: int; x: int): PairField =
      makePairFieldGaussianLike(l2),
    makeYbar = proc(sample: int; x: int): PairField =
      makePairFieldGaussianLike(l1),
    jvp = proc(x: int; dx: PairField; res: var PairField) =
      zeroInPlace(res)
      hypPairStageStapleVJPL2(hs, l1, dx, res),
    vjp = proc(x: int; ybar: PairField; res: var PairField) =
      zeroInPlace(res)
      hypPairStageChainVJPL2(hs, l1, ybar, res)
  )
  result = res.failed

proc testHypOutputStageChainAdjoint(): int =
  var hs = lo.newHypSmear(coef)
  var sg = lo.newGauge
  hs.smear(g, sg)
  let l2 = currentHypL2(hs)
  let cfg = initTestConfig(
    name = "Output stage staple/chain adjoint",
    samples = 3,
    atol = 1e-12,
    rtol = 1e-10,
    verbose = true
  )
  let res = runMixedAdjointIdentity(
    cfg = cfg,
    makeBase = proc(sample: int): int =
      0,
    makeDir = proc(sample: int; x: int): GaugeField =
      makeProbe(),
    makeYbar = proc(sample: int; x: int): PairField =
      makePairFieldGaussianLike(l2),
    jvp = proc(x: int; dx: GaugeField; res: var PairField) =
      zeroInPlace(res)
      hypOutputStageStapleVJP(hs, l2, dx, res),
    vjp = proc(x: int; ybar: PairField; res: var GaugeField) =
      zeroInPlace(res)
      hypOutputStageChainVJP(hs, l2, ybar, res)
  )
  result = res.failed

## Test: Pure reverse-mode gradient for <p, force(g)>
## This test verifies the full reverse-mode chain:
##   T(g) = <p, contractProjectTAH(g, smearVJP(g, gaugeActionDeriv(smear(g))))>
## Using only backward propagation (adjoint operations).
##
## The pure reverse-mode chain is:
##   Forward: sg = smear(g), chain = gaugeActionDeriv(sg), deriv = smearVJP(chain)
##   Backward: forcebar = p
##             (derivbar, gbar_contract) = contractProjectTAHVJP(g, deriv, forcebar)
##             chainbar = smearVJPChain(derivbar)
##             sgbar = gaugeActionDeriv_adj(chainbar)  [self-adjoint structure]
##             gbar = smearVJP(sgbar)  [by smearJVP/smearVJP duality]
##   Gradient: <gbar, dg*g>

proc testGradForceReverseModeWithCoef(testCoef: HypCoefs): TestResult =
  ## Test pure reverse mode with a specific coefficient set.
  var summary = TestResult(name: "pure reverse-mode FD")
  let eps = 1e-2

  var hs = lo.newHypSmear(testCoef)
  var sg = lo.newGauge
  hs.smear(g, sg)

  var p = makeDir()

  let cfgFD = initTestConfig(
    name = "pure reverse-mode FD",
    samples = 1,
    eps = eps,
    scale = 2.0,
    atol = 1e-5,
    rtol = 1e-5,
    fdFactor = 3.0,
    verbose = false
  )
  let cfgPure = initTestConfig(
    name = "Pure vs FoR",
    samples = 1,
    atol = 1e-8,
    rtol = 1e-8,
    verbose = false
  )

  for n in 0..<3:
    var dg = makeDir()

    var fd, fde: float
    var sAnaFoR = 0.0
    discard runFDScalarJVPDefault(
      cfg = cfgFD,
      makeBase = proc(sample: int): GaugeField =
        g,
      makeDir = proc(sample: int; x: GaugeField): GaugeField =
        dg,
      f = proc(x: GaugeField): float =
        var hsx = lo.newHypSmear(testCoef)
        var sgx = lo.newGauge
        hsx.smear(x, sgx)
        var chain = lo.newGauge
        gc.gaugeActionDeriv(sgx, chain)
        var force = lo.newGauge
        hsx.smearVJP(force, chain)
        contractProjectTAH(x, force)
        var s = 0.0
        for mu in 0..<4:
          s += redot(p[mu], force[mu])
        s,
      jvp = proc(x: GaugeField; dx: GaugeField): float =
        hs.smear(x, sg)
        var chain0 = lo.newGauge
        gc.gaugeActionDeriv(sg, chain0)

        var dsg = lo.newGauge
        hs.smearJVP(dsg, dx, sg)

        var dchain = lo.newGauge
        threads:
          for mu in 0..<4:
            dchain[mu] := 0
        gc.gaugeDerivDeriv2(sg, dsg, dchain)
        threads:
          for mu in 0..<4:
            dchain[mu] := -1.0 * dchain[mu]
        sAnaFoR = contractedForceJvp(hs, x, p, chain0, dx, dchain)
        sAnaFoR,
      accept = captureFD(addr fd, addr fde)
    )

    # ========== PURE REVERSE MODE ==========
    # Forward pass:
    restoreHypState(hs, sg)
    var chain = lo.newGauge
    gc.gaugeActionDeriv(sg, chain)
    var deriv = lo.newGauge
    hs.smearVJP(deriv, chain)
    # Note: deriv is pre-contraction; force = contractProjectTAH(g, deriv)

    # Backward pass:
    # Step 1: contractProjectTAHVJP gives both the adjoint w.r.t. the
    # pre-contraction force (`derivbar`) and the direct gauge contribution
    # from differentiating the outer contraction itself (`gbar_contract`).
    var derivbar = lo.newGauge
    var gbar_contract = lo.newGauge
    contractProjectTAHVJP(g, deriv, p, derivbar, gbar_contract)

    # Step 2: chainbar = smearVJPChain(derivbar)
    restoreHypState(hs, sg)
    var chainbar = lo.newGauge
    hs.smearVJPChain(chainbar, derivbar)

    # Step 2b: gbar_internal = smearHVPVJP(derivbar, chain)
    # This captures the contribution from internal HYP state dependence.
    restoreHypState(hs, sg)
    var gbar_internal = lo.newGauge
    hs.smearHVPVJP(gbar_internal, derivbar, chain)

    # Step 3: sgbar = adjoint of (sg → gaugeActionDeriv(sg)) applied to chainbar
    # For plaquette action, the Hessian is symmetric, so:
    #   <chainbar, gaugeDerivDeriv2(sg, dsg)> = <gaugeDerivDeriv2(sg, chainbar), dsg>
    # Thus sgbar = gaugeDerivDeriv2(sg, chainbar)
    var sgbar = lo.newGauge
    threads:
      for mu in 0..<4:
        sgbar[mu] := 0
    gc.gaugeDerivDeriv2(sg, chainbar, sgbar)
    # Negate due to force sign convention
    threads:
      for mu in 0..<4:
        sgbar[mu] := -1.0 * sgbar[mu]

    # Step 4: gbar_chain = smearJVP_adj(sgbar) = smearVJP(sgbar)
    # From testSmearTangentAdjoint: <probe, smearJVP(dg, sg)> = <smearVJP(probe), dg*g>
    # So smearJVP† = smearVJP in the appropriate sense.
    restoreHypState(hs, sg)
    var gbar_chain = lo.newGauge
    hs.smearVJP(gbar_chain, sgbar)

    # Step 6: gbar = gbar_chain + gbar_internal + gbar_contract
    var gbar = lo.newGauge
    threads:
      for mu in 0..<4:
        gbar[mu] := gbar_chain[mu] + gbar_internal[mu] + gbar_contract[mu]

    # Compute gradient: <gbar, dg*g>
    var dgEff = lo.newGauge
    threads:
      for mu in 0..<4:
        dgEff[mu] := dg[mu] * g[mu]

    var sAnaPure = 0.0
    for mu in 0..<4:
      sAnaPure += redot(gbar[mu], dgEff[mu])

    # Compare FoR with FD (always should pass)
    let errFoR = abs(sAnaFoR - fd)
    let refVal = max(1.0, abs(fd))
    let okFoR = acceptTol(cfgFD, errFoR, fde, refVal)
    let errPure = abs(sAnaPure - sAnaFoR)
    let refPure = max(1.0, max(abs(sAnaPure), abs(sAnaFoR)))
    let relPure = errPure / refPure
    let okPure = errPure <= max(cfgPure.atol, cfgPure.rtol * refPure)
    echo "    sample ", n,
         " | FoR vs FD: ", (if okFoR: "PASS" else: "FAIL"),
         " (err=", errFoR, ", fdErr=", fde, ")",
         " | Pure vs FoR: ", (if okPure: "PASS" else: "FAIL"),
         " (err=", errPure, ", rel=", relPure, ")"
    if okFoR: inc summary.passed else: inc summary.failed

    if okPure: inc summary.passed else: inc summary.failed

  result = summary


proc testGradForceReverseMode(): int =
  echo "Testing all 8 coefficient combinations to isolate bugs...\n"

  var labels: seq[string] = @[]
  for idx in 0..<adjCoefCombos.len:
    let testCoef = adjCoefCombos[idx]
    let l1 = if testCoef.alpha1 > 0: "+" else: "0"
    let l2 = if testCoef.alpha2 > 0: "+" else: "0"
    let l3 = if testCoef.alpha3 > 0: "+" else: "0"
    labels.add("(L1=" & l1 & " L2=" & l2 & " L3=" & l3 & ") " & adjCoefNames[idx])

  let summary = runCoefSweep(
    name = "pure reverse-mode gradient for <p, force(g)>",
    coefs = adjCoefCombos,
    labels = labels,
    runOne = proc(testCoef: HypCoefs): TestResult =
      testGradForceReverseModeWithCoef(testCoef)
  )
  result = summary.failed

## Minimal test for dgEff conversion: dgEff = dg * g, gbar = dgEffbar * g†
proc testDgEffConversion(): int =
  let cfg = initTestConfig(
    name = "dgEff ↔ gbar conversion",
    samples = 3,
    atol = 1e-12,
    rtol = 1e-10,
    verbose = true
  )
  let res = runAdjointIdentity(
    cfg = cfg,
    makeBase = proc(sample: int): GaugeField =
      g,
    makeDir = proc(sample: int; x: GaugeField): GaugeField =
      var dg = lo.newGauge
      dg.gaussian r
      dg,
    jvp = proc(x: GaugeField; dx: GaugeField; res: var GaugeField) =
      for mu in 0..<4:
        res[mu] := dx[mu] * x[mu],
    vjp = proc(x: GaugeField; ybar: GaugeField; res: var GaugeField) =
      for mu in 0..<4:
        res[mu] := ybar[mu] * x[mu].adj
  )
  result = res.failed

## Test that the direct contribution adjoint is correct
## Forward: dderiv = ma * dflx, Adjoint: dflxbar = ma * derivbar
proc testDirectContribution(): int =
  let cfg = initTestConfig(
    name = "direct contribution (scaling)",
    samples = 3,
    atol = 1e-12,
    rtol = 1e-10,
    verbose = true
  )
  let ma = 0.6
  let res = runAdjointIdentity(
    cfg = cfg,
    makeBase = proc(sample: int): int =
      0,
    makeDir = proc(sample: int; x: int): GaugeField =
      var dflx = lo.newGauge
      dflx.gaussian r
      dflx,
    jvp = proc(x: int; dx: GaugeField; res: var GaugeField) =
      for mu in 0..<4:
        res[mu] := ma * dx[mu],
    vjp = proc(x: int; ybar: GaugeField; res: var GaugeField) =
      for mu in 0..<4:
        res[mu] := ma * ybar[mu]
  )
  result = res.failed

proc testSmearDeriv2AdjInternal(coef: HypCoefs): int =
  ## Test internal adjoint for smearHVPVJP with a specific coefficient set.
  let cfg = initTestConfig(
    name = "smearHVPVJP adjoint identity (internal)",
    samples = 3,
    atol = 1e-10,
    rtol = 1e-8,
    verbose = true
  )

  var hs = lo.newHypSmear(coef)
  var sg = lo.newGauge
  var deriv0 = lo.newGauge

  let res = runAdjointIdentity(
    cfg = cfg,
    makeBase = proc(sample: int): GaugeField =
      var chain = lo.newGauge
      chain.randomTAH r
      chain,
    makeDir = proc(sample: int; x: GaugeField): GaugeField =
      var dg = lo.newGauge
      dg.randomTAH r
      dg,
    jvp = proc(x: GaugeField; dx: GaugeField; res: var GaugeField) =
      var dg = lo.newGauge
      for mu in 0..<4:
        dg[mu] := dx[mu] * g[mu].adj
      primeHypReverse(hs, sg, deriv0, x)
      hs.smearHVP(res, x, dg),
    vjp = proc(x: GaugeField; ybar: GaugeField; res: var GaugeField) =
      primeHypReverse(hs, sg, deriv0, x)
      hs.smearHVPVJP(res, ybar, x)
  )
  result = res.failed

proc testLevel3Direct(): int =
  let cfg = initTestConfig(
    name = "Level 3 direct contribution",
    samples = 3,
    atol = 1e-12,
    rtol = 1e-10,
    verbose = true
  )

  let ma3 = 1 - coef.alpha3

  var hs = lo.newHypSmear(coef)
  var sg = lo.newGauge
  hs.smear(g, sg)
  let flx = hs.state.flx

  let res = runAdjointIdentity(
    cfg = cfg,
    makeBase = proc(sample: int): GaugeField =
      var chain = lo.newGauge
      chain.gaussian r
      chain,
    makeDir = proc(sample: int; x: GaugeField): GaugeField =
      var dgEff = lo.newGauge
      dgEff.gaussian r
      dgEff,
    jvp = proc(x: GaugeField; dx: GaugeField; res: var GaugeField) =
      var dflx = lo.newGauge
      threads:
        for mu in 0..<4:
          dflx[mu] := ma3 * dx[mu]
      var dfc = lo.newGauge
      threads:
        for mu in 0..<4:
          for i in dfc[mu]:
            dfc[mu][i].projectUHVP(flx[mu][i], x[mu][i], dflx[mu][i])
      let resPtr = cast[ptr GaugeField](unsafeAddr(res))
      threads:
        for mu in 0..<4:
          resPtr[][mu] := ma3 * dfc[mu],
    vjp = proc(x: GaugeField; ybar: GaugeField; res: var GaugeField) =
      var dfcbar = lo.newGauge
      threads:
        for mu in 0..<4:
          dfcbar[mu] := ma3 * ybar[mu]
      var dflxbar = lo.newGauge
      threads:
        for mu in 0..<4:
          for i in dflxbar[mu]:
            dflxbar[mu][i].projectUHVPVJP_dx(flx[mu][i], x[mu][i], dfcbar[mu][i])
      let resPtr = cast[ptr GaugeField](unsafeAddr(res))
      threads:
        for mu in 0..<4:
          resPtr[][mu] := ma3 * dflxbar[mu]
  )
  result = res.failed

proc testLevel2Projection(): int =
  let cfg = initTestConfig(
    name = "Level 2 projections only",
    samples = 3,
    atol = 1e-12,
    rtol = 1e-10,
    verbose = true
  )

  var hs = lo.newHypSmear(coef)
  var sg = lo.newGauge
  hs.smear(g, sg)

  type Fld = type(g[0])
  let l2x = hs.state.l2x
  var l2: typeof(hs.state.l2x)
  when compiles(hs.state.l2):
    l2 = hs.state.l2
  else:
    l2 = newFieldArray2(lo, Fld, [4,4], mu!=nu)
    threads:
      for mu in 0..<4:
        for nu in 0..<4:
          if nu != mu:
            for e in l2[mu,nu]:
              l2[mu,nu][e].projectU l2x[mu,nu][e]
  var chain0L2 = newOneOf(l2x)
  threads:
    for mu in 0..<4:
      for nu in 0..<4:
        if nu != mu:
          chain0L2[mu,nu].gaussian r

  type FA = typeof(l2x)
  type Y = Pair[FA, FA]

  let res = runAdjointIdentity(
    cfg = cfg,
    makeBase = proc(sample: int): int =
      0,
    makeDir = proc(sample: int; x: int): FA =
      var dl2x = newOneOf(l2x)
      threads:
        for mu in 0..<4:
          for nu in 0..<4:
            if nu != mu:
              dl2x[mu,nu].gaussian r
      dl2x,
    jvp = proc(x: int; dx: FA; res: var Y) =
      var dl2 = newOneOf(dx)
      var dfl2 = newOneOf(dx)
      threads:
        for mu in 0..<4:
          for nu in 0..<4:
            if nu != mu:
              for i in dl2[mu,nu]:
                dl2[mu,nu][i].projectUJVP(l2[mu,nu][i], l2x[mu,nu][i], dx[mu,nu][i])
      threads:
        for mu in 0..<4:
          for nu in 0..<4:
            if nu != mu:
              for i in dfl2[mu,nu]:
                dfl2[mu,nu][i].projectUHVPu(l2[mu,nu][i], l2x[mu,nu][i], chain0L2[mu,nu][i], dx[mu,nu][i])
      res.a := dfl2
      res.b := dl2,
    vjp = proc(x: int; ybar: Y; res: var FA) =
      var dl2xbar = newOneOf(l2x)
      threads:
        for mu in 0..<4:
          for nu in 0..<4:
            if nu != mu:
              var tmp {.noinit.}: type(dl2xbar[mu,nu][0])
              for i in dl2xbar[mu,nu]:
                tmp.projectUVJP(l2[mu,nu][i], l2x[mu,nu][i], ybar.b[mu,nu][i])
                dl2xbar[mu,nu][i] := tmp
      threads:
        for mu in 0..<4:
          for nu in 0..<4:
            if nu != mu:
              var tmp {.noinit.}: type(dl2xbar[mu,nu][0])
              for i in dl2xbar[mu,nu]:
                tmp.projectUHVPVJP_dx(l2[mu,nu][i], l2x[mu,nu][i], chain0L2[mu,nu][i], ybar.a[mu,nu][i])
                dl2xbar[mu,nu][i] += tmp
      res := dl2xbar,
  )
  result = res.failed

var fail = 0

# Section A: baseline FD (smear action/force)
fail += runGaugeActionForceFD(
  env,
  hypAction, hypForce,
  actionForceCfgDefault("hypAction/hypForce")
).failed
fail += runGaugeActionForceFD(
  env,
  hypActionPrealloc, hypForcePrealloc,
  actionForceCfgDefault("hypActionPrealloc/hypForcePrealloc")
).failed
fail += testSmearForceConsistency()

# Section B: JVP/HVP FD
fail += testSmearTangent()
fail += testGaugeDerivDeriv2()
fail += testSmearDeriv2()
fail += testSmearDeriv2FixedChain()

# Section C: smearHVP level sweep (physics coefficients)
let levelCoefs = [
  HypCoefs(alpha1: 0.0, alpha2: 0.0, alpha3: 0.3),
  HypCoefs(alpha1: 0.75, alpha2: 0.0, alpha3: 0.3),
  HypCoefs(alpha1: 0.0, alpha2: 0.6, alpha3: 0.3),
  HypCoefs(alpha1: 0.0, alpha2: 0.6, alpha3: 0.0),
  HypCoefs(alpha1: 0.0, alpha2: 0.0, alpha3: 0.0),
  HypCoefs(alpha1: 0.75, alpha2: 0.0, alpha3: 0.0),
  HypCoefs(alpha1: 0.75, alpha2: 0.6, alpha3: 0.0),
  HypCoefs(alpha1: 0.75, alpha2: 0.6, alpha3: 0.3),
]
let levelLabels = [
  "Level 3 only (α₁=α₂=0, α₃=0.3)",
  "L1+L3 only (α₂=0)",
  "L2+L3 only (α₁=0)",
  "L2 only (α₁=α₃=0)",
  "Level 0 (all α=0, identity)",
  "Level 1 only (α₂=α₃=0)",
  "Levels 1+2 (α₃=0)",
  "Full HYP (α₁=0.75, α₂=0.6, α₃=0.3)",
]
for idx in 0..<levelCoefs.len:
  fail += testSmearDeriv2Level(levelLabels[idx], levelCoefs[idx])

# Full chain gradient
fail += testGradForceFullChain()

# Section C: adjoint identities + reverse mode
fail += testSmearTangentAdjoint()
fail += testSmearDerivAdjChain()
fail += testGradForceReverseMode()

# Section D: individual stage operators
fail += testHypPairStage1JVP()
fail += testHypPairStage2JVP()
fail += testHypOutputStageJVP()
fail += testHypPairStage1ChainAdjoint()
fail += testHypPairStage2ChainAdjoint()
fail += testHypOutputStageChainAdjoint()

# Section E: term-by-term / internal adjoints
fail += testDgEffConversion()
fail += testDirectContribution()
block:
  let res = runCoefSweep(
    name = "smearHVPVJP coefficient sweep",
    coefs = adjCoefCombos,
    labels = adjCoefNames,
    runOne = proc(c: HypCoefs): TestResult =
      let failed = testSmearDeriv2AdjInternal(c)
      result.name = "smearHVPVJP combo"
      if failed == 0:
        result.passed = 1
      else:
        result.failed = failed
  )
  fail += res.failed

# Stage-specific adjoint checks
fail += testLevel3Direct()
fail += testLevel2Projection()

if fail == 0:
  qexFinalize()
else:
  qexAbort(fail)
