import math
import random
import maths/[matrixFunctions, projUderiv, types, matrixConcept, complexNumbers]
import ./utils

var gRng = initRand(123456)

type
  C64 = ComplexType[float64]
  M = MatrixArray[3,3,C64]
  Dir = object
    dx: M
    dc: M

proc randUnitMat(rng: var Rand; s = 1.0): M =
  var x = randMat[M](rng, s)
  let nx = sqrt(x.norm2)
  if nx != 0.0:
    x *= (1.0 / nx)
  x

proc projectUHVPRef(res: var M; x: M; c: M; dx: M; dc: M; eps = 1e-20) =
  discard eps
  var z, y: M
  projectUrsqrt(z, x, eps)
  inverse(y, z)

  var u: M
  u := x * z

  var r0: M
  r0 := c * z

  var b: M
  b := u.adj * r0

  var s: M
  var t2: M
  t2 := b + b.adj
  sylsolve(s, y, t2)

  var da, tmp: M
  tmp := x.adj * dx
  da := tmp + tmp.adj

  var dy: M
  sylsolve(dy, y, da)

  var dz: M
  tmp := dy * z
  dz := z * tmp
  dz := -dz

  var dr0: M
  dr0 := dc * z
  tmp := c * dz
  dr0 += tmp

  var du: M
  du := dx * z
  tmp := x * dz
  du += tmp

  var db: M
  db := du.adj * r0
  tmp := u.adj * dr0
  db += tmp

  var q: M
  q := dy * s

  var rhs: M
  rhs := db + db.adj
  rhs -= q
  rhs -= q.adj

  var ds: M
  sylsolve(ds, y, rhs)

  res := dr0
  tmp := dx * s
  res -= tmp
  tmp := x * ds
  res -= tmp

proc runProjHvpCase(
    caseName: string,
    cfg: TestConfig,
    x, c, dx, dc, probe: M,
    f: auto,
    jvp: auto
  ): TestResult =
  echo "  ", caseName
  runFDJVPDefault(
    cfg = cfg,
    makeBase = proc(sample: int): Pair[M, M] =
      makePair(x, c),
    makeDir = proc(sample: int; x: Pair[M, M]): Pair[M, M] =
      makePair(dx, dc),
    makeProbe = proc(sample: int; x: Pair[M, M]): M =
      probe,
    f = f,
    jvp = jvp
  )

proc runProjHvpCaseCustom(
    caseName: string,
    cfg: TestConfig,
    x, c, dx, dc, probe: M,
    perturb: auto,
    f: auto,
    jvp: auto
  ): TestResult =
  echo "  ", caseName
  runFDJVP(
    cfg = cfg,
    makeBase = proc(sample: int): Pair[M, M] =
      makePair(x, c),
    makeDir = proc(sample: int; x: Pair[M, M]): Pair[M, M] =
      makePair(dx, dc),
    makeProbe = proc(sample: int; x: Pair[M, M]): M =
      probe,
    perturb = perturb,
    f = f,
    jvp = jvp
  )

proc testProjTangent(rng: var Rand): int =
  let cfg = initTestConfig(
    name = "projectU JVP (tangent)",
    samples = 10,
    eps = 0.1,
    atol = 1e-8,
    rtol = 1e-6,
    verbose = true
  )
  let res = runFDJVPDefault(
    cfg = cfg,
    makeBase = proc(sample: int): M =
      let sx = 1.0 + 0.1 * float(sample)
      randUnitMat(gRng, sx),
    makeDir = proc(sample: int; x: M): M =
      let sdx = 0.7 + 0.05 * float(sample)
      randUnitMat(gRng, sdx),
    makeProbe = proc(sample: int; x: M): M =
      randMat[M](gRng, 1.0),
    f = proc(x: M; res: var M) =
      projectU(res, x),
    jvp = proc(x: M; dx: M; res: var M) =
      var u: M
      projectU(u, x)
      projectUJVP(res, u, x, dx)
  )
  result = res.failed

proc testProjDeriv2FD(rng: var Rand): int =
  echo "\n=== projectUHVP FD checks ==="
  var fail = 0
  let eps = 0.1

  var x0 = randMat[M](gRng)
  var dx0 = randMat[M](gRng)
  var c0 = randMat[M](gRng)
  var dc0 = randMat[M](gRng)
  var a0 = randMat[M](gRng)
  var zeroMat: M
  zeroMat := 0

  let cfgHvp = initTestConfig(
    name = "projectUHVP FD",
    samples = 1,
    eps = eps,
    atol = 1e-6,
    rtol = 1e-5,
    verbose = false
  )
  fail += runProjHvpCase(
    "Test 0: Reference HVP",
    cfgHvp,
    x0, c0, dx0, dc0, a0,
    f = proc(x: Pair[M, M]; res: var M) =
      projectUVJP(res, x.a, x.b),
    jvp = proc(x: Pair[M, M]; dx: Pair[M, M]; res: var M) =
      projectUHVPRef(res, x.a, x.b, dx.a, dx.b)
  ).failed

  fail += runProjHvpCase(
    "Test 1: Random matrices with dc",
    cfgHvp,
    x0, c0, dx0, dc0, a0,
    f = proc(x: Pair[M, M]; res: var M) =
      projectUVJP(res, x.a, x.b),
    jvp = proc(x: Pair[M, M]; dx: Pair[M, M]; res: var M) =
      var u: M
      projectU(u, x.a)
      projectUHVPu(res, u, x.a, x.b, dx.a, dx.b)
  ).failed

  fail += runProjHvpCase(
    "Test 2: Random matrices, dc=0 (fixed chain)",
    cfgHvp,
    x0, c0, dx0, zeroMat, a0,
    f = proc(x: Pair[M, M]; res: var M) =
      projectUVJP(res, x.a, x.b),
    jvp = proc(x: Pair[M, M]; dx: Pair[M, M]; res: var M) =
      var u: M
      projectU(u, x.a)
      projectUHVPu(res, u, x.a, x.b, dx.a, dx.b)
  ).failed

  var xg = randMat[M](gRng)
  var g: M
  projectU(g, xg)
  var chain = randMat[M](gRng)
  var pDir = randTAH[M](gRng)
  var a = randMat[M](gRng)
  fail += runProjHvpCaseCustom(
    "Test 3: Unitary input with multiplicative tangent",
    cfgHvp,
    g, chain, pDir, zeroMat, a,
    perturb = proc(x: Pair[M, M]; dx: Pair[M, M]; eps: float; res: var Pair[M, M]) =
      res.a := exp(eps * dx.a) * x.a
      res.b := x.b,
    f = proc(x: Pair[M, M]; res: var M) =
      projectUVJP(res, x.a, x.b),
    jvp = proc(x: Pair[M, M]; dx: Pair[M, M]; res: var M) =
      var u: M
      projectU(u, x.a)
      let dxMult = dx.a * x.a
      projectUHVPu(res, u, x.a, x.b, dxMult, zeroMat)
  ).failed

  var chainTAH = randTAH[M](gRng)
  var pTAH = randTAH[M](gRng)
  fail += runProjHvpCaseCustom(
    "Test 4: Unitary input, TAH chain, multiplicative tangent",
    cfgHvp,
    g, chainTAH, pTAH, zeroMat, pTAH,
    perturb = proc(x: Pair[M, M]; dx: Pair[M, M]; eps: float; res: var Pair[M, M]) =
      res.a := exp(eps * dx.a) * x.a
      res.b := x.b,
    f = proc(x: Pair[M, M]; res: var M) =
      projectUVJP(res, x.a, x.b),
    jvp = proc(x: Pair[M, M]; dx: Pair[M, M]; res: var M) =
      var u: M
      projectU(u, x.a)
      let dxMult = dx.a * x.a
      projectUHVPu(res, u, x.a, x.b, dxMult, zeroMat)
  ).failed

  let cfgDxDc = initTestConfig(
    name = "projectUHVP with dx and dc",
    samples = 1,
    eps = 0.01,
    atol = 0.1,
    rtol = 0.01,
    verbose = true
  )
  block:
    let x = randMat[M](gRng)
    let c = randMat[M](gRng)
    let dx = randMat[M](gRng)
    let dc = randMat[M](gRng)
    let probe = randMat[M](gRng)
    fail += runProjHvpCase(
      "Test 5: dx and dc",
      cfgDxDc,
      x, c, dx, dc, probe,
      f = proc(x: Pair[M, M]; res: var M) =
        projectUVJP(res, x.a, x.b),
      jvp = proc(x: Pair[M, M]; dx: Pair[M, M]; res: var M) =
        projectUHVP(res, x.a, x.b, dx.a, dx.b)
    ).failed

  let cfgDcOnly = initTestConfig(
    name = "projectUHVP dc-only",
    samples = 1,
    eps = 0.01,
    atol = 0.1,
    rtol = 0.01,
    verbose = true
  )
  block:
    let x = randMat[M](gRng)
    let c = randMat[M](gRng)
    var dxZero: M
    dxZero := 0
    let dc = randMat[M](gRng)
    let probe = randMat[M](gRng)
    fail += runProjHvpCase(
      "Test 6: dc-only",
      cfgDcOnly,
      x, c, dxZero, dc, probe,
      f = proc(x: Pair[M, M]; res: var M) =
        projectUVJP(res, x.a, x.b),
      jvp = proc(x: Pair[M, M]; dx: Pair[M, M]; res: var M) =
        projectUHVP(res, x.a, x.b, dx.a, dx.b)
    ).failed

  result = fail

proc testProjDeriv2TemplateVsDirect(rng: var Rand): int =
  echo "\n=== projectUHVP template vs direct ==="
  var fail = 0
  let eps = 0.1

  var xg = randMat[M](gRng)
  var g: M
  projectU(g, xg)
  var chainTAH = randTAH[M](gRng)
  var pTAH = randTAH[M](gRng)

  var dxMult = pTAH * g

  var result1: M
  projectUHVP(result1, g, chainTAH, dxMult)

  var u_g: M
  projectU(u_g, g)
  var dcZero: M
  dcZero := 0
  var result2: M
  projectUHVPu(result2, u_g, g, chainTAH, dxMult, dcZero)

  let norm1 = result1.norm2
  let norm2 = result2.norm2
  let diff = (result1 - result2).norm2

  echo "  |result1 (template)|^2 = ", norm1
  echo "  |result2 (direct)|^2  = ", norm2
  echo "  |result1 - result2|^2 = ", diff

  let cfg = initTestConfig(
    name = "projectUHVP template vs direct (FD)",
    samples = 1,
    eps = eps,
    atol = 1e-6,
    rtol = 1e-5,
    verbose = false
  )
  var fd5, fde5: float
  block:
    let res = runFDJVP(
      cfg = cfg,
      makeBase = proc(sample: int): Pair[M, M] =
        makePair(g, chainTAH),
      makeDir = proc(sample: int; x: Pair[M, M]): Pair[M, M] =
        var z: M
        z := 0
        makePair(pTAH, z),
      makeProbe = proc(sample: int; x: Pair[M, M]): M =
        pTAH,
      perturb = proc(x: Pair[M, M]; dx: Pair[M, M]; eps: float; res: var Pair[M, M]) =
        res.a := exp(eps * dx.a) * x.a
        res.b := x.b,
      f = proc(x: Pair[M, M]; res: var M) =
        projectUVJP(res, x.a, x.b),
      jvp = proc(x: Pair[M, M]; dx: Pair[M, M]; res: var M) =
        var u: M
        projectU(u, x.a)
        var dc0zero: M
        dc0zero := 0
        let dxMult = dx.a * x.a
        projectUHVPu(res, u, x.a, x.b, dxMult, dc0zero),
      accept = captureFD(addr fd5, addr fde5)
    )
    fail += res.failed

  let ana1 = redot(pTAH, result1)
  let ana2 = redot(pTAH, result2)

  var cfgCmp = cfg
  cfgCmp.verbose = true
  cfgCmp.atol = max(cfg.atol, cfg.fdFactor * fde5)
  cfgCmp.rtol = 0.0
  if not reportCompare(cfgCmp, 0, fd5, ana1, "Template vs FD"): inc fail
  if not reportCompare(cfgCmp, 0, fd5, ana2, "Direct vs FD"): inc fail

  if diff < 1e-20:
    echo "  Template and Direct are IDENTICAL"
  else:
    echo "  Template and Direct DIFFER"
    inc fail

  result = fail

proc testProjDeriv2Linearity(rng: var Rand): int =
  let cfg = initTestConfig(
    name = "projectUHVP linearity",
    samples = 3,
    atol = 1e-10,
    rtol = 1e-10,
    verbose = true
  )
  var x0 = randMat[M](gRng)
  var c0 = randMat[M](gRng)
  let res = runLinearity(
    cfg = cfg,
    makeVec = proc(sample: int): Dir =
      Dir(dx: randMat[M](gRng), dc: randMat[M](gRng)),
    makeProbe = proc(sample: int): M =
      randMat[M](gRng),
    combine = proc(v1, v2: Dir; a, b: float; res: var Dir) =
      res.dx = a * v1.dx + b * v2.dx
      res.dc = a * v1.dc + b * v2.dc,
    apply = proc(v: Dir; res: var M) =
      projectUHVP(res, x0, c0, v.dx, v.dc)
  )
  result = res.failed

proc testProjAdjDx(rng: var Rand): int =
  let cfg = initTestConfig(
    name = "projectUHVP adjoint w.r.t dx",
    samples = 5,
    atol = 1e-10,
    rtol = 1e-10,
    verbose = true
  )
  let res = runAdjointIdentity(
    cfg = cfg,
    makeBase = proc(sample: int): Pair[M, M] =
      makePair(randMat[M](gRng), randMat[M](gRng)),
    makeDir = proc(sample: int; x: Pair[M, M]): M =
      randMat[M](gRng),
    jvp = proc(x: Pair[M, M]; dx: M; res: var M) =
      var u: M
      projectU(u, x.a)
      var dcZero: M
      dcZero := 0
      projectUHVPu(res, u, x.a, x.b, dx, dcZero),
    vjp = proc(x: Pair[M, M]; ybar: M; res: var M) =
      var u: M
      projectU(u, x.a)
      projectUHVPVJP_dx(res, u, x.a, x.b, ybar)
  )
  result = res.failed

proc testProjAdjDc(rng: var Rand): int =
  let cfg = initTestConfig(
    name = "projectUHVP adjoint w.r.t dc",
    samples = 5,
    atol = 1e-10,
    rtol = 1e-10,
    verbose = true
  )
  let res = runAdjointIdentity(
    cfg = cfg,
    makeBase = proc(sample: int): Pair[M, M] =
      makePair(randMat[M](gRng), randMat[M](gRng)),
    makeDir = proc(sample: int; x: Pair[M, M]): M =
      randMat[M](gRng),
    jvp = proc(x: Pair[M, M]; dc: M; res: var M) =
      var u: M
      projectU(u, x.a)
      var dxZero: M
      dxZero := 0
      projectUHVPu(res, u, x.a, x.b, dxZero, dc),
    vjp = proc(x: Pair[M, M]; ybar: M; res: var M) =
      var u: M
      projectU(u, x.a)
      projectUHVPVJP_dc(res, u, x.a, x.b, ybar)
  )
  result = res.failed

proc testProjAdjChain(rng: var Rand): int =
  let cfg = initTestConfig(
    name = "projectUVJP chain adjoint",
    samples = 5,
    atol = 1e-10,
    rtol = 1e-10,
    verbose = true
  )
  let res = runAdjointIdentity(
    cfg = cfg,
    makeBase = proc(sample: int): M =
      randMat[M](gRng),
    makeDir = proc(sample: int; x: M): M =
      randMat[M](gRng),
    jvp = proc(x: M; dx: M; res: var M) =
      projectUVJP(res, x, dx),
    vjp = proc(x: M; ybar: M; res: var M) =
      projectUVJPChain(res, x, ybar)
  )
  result = res.failed

proc testProjJVPAdjoint(rng: var Rand): int =
  let cfg = initTestConfig(
    name = "projectUJVP adjoint identity",
    samples = 5,
    atol = 1e-12,
    rtol = 1e-10,
    verbose = true
  )
  let res = runAdjointIdentity(
    cfg = cfg,
    makeBase = proc(sample: int): M =
      randUnitMat(gRng),
    makeDir = proc(sample: int; x: M): M =
      randUnitMat(gRng),
    jvp = proc(x: M; dx: M; res: var M) =
      var u: M
      projectU(u, x)
      projectUJVP(res, u, x, dx),
    vjp = proc(x: M; ybar: M; res: var M) =
      var u: M
      projectU(u, x)
      projectUVJP(res, u, x, ybar)
  )
  result = res.failed

proc testAdjointTerms(rng: var Rand): int =
  echo "\n=== projectUHVP adjoint term checks ==="
  var fail = 0

  # Term 1: r = -dx S, S Hermitian
  block:
    let cfg = initTestConfig(
      name = "Term -dx S",
      samples = 1,
      atol = 1e-12,
      rtol = 1e-12,
      verbose = false
    )
    let res = runAdjointIdentity(
      cfg = cfg,
      makeBase = proc(sample: int): M =
        hermitianize(randMat[M](gRng)),
      makeDir = proc(sample: int; x: M): M =
        randMat[M](gRng),
      jvp = proc(x: M; dx: M; res: var M) =
        res := -dx * x,
      vjp = proc(x: M; ybar: M; res: var M) =
        res := -ybar * x
    )
    reportSummary(res)
    fail += res.failed

  # Term 2: r = Z dx† R0, Z Hermitian
  block:
    let cfg = initTestConfig(
      name = "Term Z dx† R0",
      samples = 1,
      atol = 1e-12,
      rtol = 1e-12,
      verbose = false
    )
    let res = runAdjointIdentity(
      cfg = cfg,
      makeBase = proc(sample: int): Pair[M, M] =
        makePair(hermitianize(randMat[M](gRng)), randMat[M](gRng)),
      makeDir = proc(sample: int; x: Pair[M, M]): M =
        randMat[M](gRng),
      jvp = proc(x: Pair[M, M]; dx: M; res: var M) =
        res := x.a * dx.adj * x.b,
      vjp = proc(x: Pair[M, M]; ybar: M; res: var M) =
        res := x.b * ybar.adj * x.a
    )
    reportSummary(res)
    fail += res.failed

  # Term 3: r = X† dx
  block:
    let cfg = initTestConfig(
      name = "Term X† dx",
      samples = 1,
      atol = 1e-12,
      rtol = 1e-12,
      verbose = false
    )
    let res = runAdjointIdentity(
      cfg = cfg,
      makeBase = proc(sample: int): M =
        randMat[M](gRng),
      makeDir = proc(sample: int; x: M): M =
        randMat[M](gRng),
      jvp = proc(x: M; dx: M; res: var M) =
        res := x.adj * dx,
      vjp = proc(x: M; ybar: M; res: var M) =
        res := x * ybar
    )
    reportSummary(res)
    fail += res.failed

  # Term 4: r = dx† X
  block:
    let cfg = initTestConfig(
      name = "Term dx† X",
      samples = 1,
      atol = 1e-12,
      rtol = 1e-12,
      verbose = false
    )
    let res = runAdjointIdentity(
      cfg = cfg,
      makeBase = proc(sample: int): M =
        randMat[M](gRng),
      makeDir = proc(sample: int; x: M): M =
        randMat[M](gRng),
      jvp = proc(x: M; dx: M; res: var M) =
        res := dx.adj * x,
      vjp = proc(x: M; ybar: M; res: var M) =
        res := x * ybar.adj
    )
    reportSummary(res)
    fail += res.failed

  # Term 5: r = dZ X† R0, dZ Hermitian
  block:
    let cfg = initTestConfig(
      name = "Term dZ X† R0",
      samples = 1,
      atol = 1e-12,
      rtol = 1e-12,
      verbose = false
    )
    let res = runAdjointIdentity(
      cfg = cfg,
      makeBase = proc(sample: int): Pair[M, M] =
        makePair(randMat[M](gRng), randMat[M](gRng)),
      makeDir = proc(sample: int; x: Pair[M, M]): M =
        hermitianize(randMat[M](gRng)),
      jvp = proc(x: Pair[M, M]; dx: M; res: var M) =
        res := dx * x.a.adj * x.b,
      vjp = proc(x: Pair[M, M]; ybar: M; res: var M) =
        res := ybar * x.b.adj * x.a
    )
    reportSummary(res)
    fail += res.failed

  # Term 6: Sylvester solve adjoint
  block:
    let cfg = initTestConfig(
      name = "Sylvester adjoint",
      samples = 1,
      atol = 1e-12,
      rtol = 1e-12,
      verbose = false
    )
    let res = runAdjointIdentity(
      cfg = cfg,
      makeBase = proc(sample: int): M =
        let x = randMat[M](gRng)
        x.adj * x,
      makeDir = proc(sample: int; x: M): M =
        randMat[M](gRng),
      jvp = proc(x: M; dx: M; res: var M) =
        sylsolve(res, x, dx),
      vjp = proc(x: M; ybar: M; res: var M) =
        sylsolve(res, x, ybar)
    )
    reportSummary(res)
    fail += res.failed

  result = fail

var fail = 0

# Section A: JVP/HVP FD
fail += testProjTangent(gRng)
fail += testProjDeriv2FD(gRng)
fail += testProjDeriv2TemplateVsDirect(gRng)
fail += testProjDeriv2Linearity(gRng)

# Section B: adjoint identities
fail += testProjAdjDx(gRng)
fail += testProjAdjDc(gRng)
fail += testProjAdjChain(gRng)
fail += testProjJVPAdjoint(gRng)

# Section C: term-by-term adjoint checks
fail += testAdjointTerms(gRng)

if fail == 0:
  echo "\nAll tests passed"
else:
  echo "\n", fail, " tests failed"
  quit(fail)
