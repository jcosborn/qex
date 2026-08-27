import math
import os
import random
import sequtils
import macros
import algorithms/numdiff
import maths/matrixFunctions
import qex
import gauge/gaugeUtils
import gauge/stoutsmear

const
  defaultSamples = 5
  defaultEps = 1e-2
  defaultScale = 2.0
  defaultAtol = 1e-6
  defaultRtol = 1e-5
  defaultFdFactor = 3.0
  defaultVerbose = true
  actionForceDefaultEps = 0.125
  actionForceDefaultAtol = 5e-8
  actionForceDefaultRtol = 1e-7
  actionForceDefaultFdFactor = 4.0

type
  TestConfig* = object
    name*: string
    samples*: int
    eps*: float
    scale*: float
    atol*: float
    rtol*: float
    fdFactor*: float
    verbose*: bool

  TestResult* = object
    name*: string
    passed*: int
    failed*: int

  GaugeTestEnv*[L, G, R] = object
    lo*: L
    g*: G
    r*: R
    lat*: seq[int]
    gaugefile*: string
    seed*: uint64
    hasFile*: bool

  Pair*[A, B] = object
    a*: A
    b*: B


  AcceptFD* = proc(
    cfg: TestConfig;
    ana, fd, err, fdErr, refVal: float;
    sample: int
  ): bool

proc initTestConfig*(
    name: string;
    samples = defaultSamples;
    eps = defaultEps;
    scale = defaultScale;
    atol = defaultAtol;
    rtol = defaultRtol;
    fdFactor = defaultFdFactor;
    verbose = defaultVerbose
  ): TestConfig =
  result.name = name
  result.samples = samples
  result.eps = eps
  result.scale = scale
  result.atol = atol
  result.rtol = rtol
  result.fdFactor = fdFactor
  result.verbose = verbose

proc add(acc: var TestResult, other: TestResult) =
  if acc.name.len == 0:
    acc.name = other.name
  acc.passed += other.passed
  acc.failed += other.failed

proc reportSummary*(r: TestResult) =
  if r.name.len > 0:
    if r.failed == 0:
      echo r.name, ": all ", r.passed, " tests passed"
    else:
      echo r.name, ": ", r.failed, " failed (", r.passed, " passed)"

proc reportBanner*(name: string) =
  if name.len > 0:
    echo "### Testing ", name

proc reportBanner*(cfg: TestConfig) =
  reportBanner(cfg.name)

proc makePair*[A, B](a: A; b: B): Pair[A, B] =
  Pair[A, B](a: a, b: b)

proc seedFromDir*[Y, DX](dx: DX): Y =
  when Y is Pair:
    when DX is Pair:
      result = makePair(
        seedFromDir[typeof(default(Y).a), typeof(default(DX).a)](dx.a),
        seedFromDir[typeof(default(Y).b), typeof(default(DX).b)](dx.b)
      )
    else:
      result = makePair(
        seedFromDir[typeof(default(Y).a), DX](dx),
        seedFromDir[typeof(default(Y).b), DX](dx)
      )
  else:
    when Y is DX:
      result = allocLike(dx)
      copyLike(result, dx)
    elif DX is Pair:
      result = seedFromDir[Y, typeof(default(DX).a)](dx.a)
    else:
      {.error: "seedFromDir requires Y to be DX or derived from Pair DX".}

proc allocLike*[T](x: T): T =
  when T is Pair:
    Pair[typeof(allocLike(x.a)), typeof(allocLike(x.b))](
      a: allocLike(x.a),
      b: allocLike(x.b)
    )
  elif T is seq:
    var y: T
    y.setLen(x.len)
    for i in 0..<x.len:
      y[i] = allocLike(x[i])
    y
  elif T is array:
    var y: T
    for i in 0..<y.len:
      y[i] = allocLike(x[i])
    y
  elif compiles(newOneOf(x)):
    newOneOf(x)
  else:
    var y: T
    y


proc zeroLike*[T](x: T): T =
  when T is Pair:
    Pair[typeof(zeroLike(x.a)), typeof(zeroLike(x.b))](
      a: zeroLike(x.a),
      b: zeroLike(x.b)
    )
  elif T is seq:
    var y = allocLike(x)
    for i in 0..<x.len:
      y[i] = zeroLike(x[i])
    y
  elif T is array:
    var y: T
    for i in 0..<y.len:
      y[i] = zeroLike(x[i])
    y
  elif T is FieldArray:
    var y = allocLike(x)
    for i in 0..<y.arr.len:
      if y.arr[i] != nil:
        y.arr[i] := 0
    y
  else:
    var y = allocLike(x)
    when compiles(y := 0):
      y := 0
    elif compiles(assignM(y, 0)):
      assignM(y, 0)
    elif compiles(y = 0):
      y = 0
    else:
      discard
    y


proc copyLike*[T](dst: var T; src: T) =
  when T is Pair:
    copyLike(dst.a, src.a)
    copyLike(dst.b, src.b)
  elif compiles(dst := src):
    dst := src
  elif compiles(assignM(dst, src)):
    assignM(dst, src)
  else:
    dst = src


proc dotLike*[T](a, b: T): float =
  when T is Pair:
    dotLike(a.a, b.a) + dotLike(a.b, b.b)
  elif T is FieldArray:
    dotFieldArray(a, b)
  elif compiles(redot(a, b)):
    redot(a, b)
  elif T is seq or T is array:
    var s = 0.0
    for i in 0..<a.len:
      s += dotLike(a[i], b[i])
    s
  else:
    a * b

proc dotFieldArray*[N:static[int],V:static[int],T](a, b: FieldArray[N,V,T]): float =
  var s = 0.0
  when N == 1:
    for i in 0..<a.shape[0]:
      let fa = a[i]
      if fa != nil:
        s += redot(fa, b[i])
  elif N == 2:
    for i in 0..<a.shape[0]:
      for j in 0..<a.shape[1]:
        let fa = a[i, j]
        if fa != nil:
          s += redot(fa, b[i, j])
  else:
    for i in 0..<a.arr.len:
      let fa = a.arr[i]
      if fa != nil:
        s += redot(fa, b.arr[i])
  s

proc dotLike*[N:static[int],V:static[int],T](a, b: FieldArray[N,V,T]): float =
  dotFieldArray(a, b)

proc norm2Like*[T](x: T): float =
  dotLike(x, x)

proc perturbLike*[T](x: T; dx: T; eps: float; res: var T; useExp: static[bool] = true) =
  when T is Pair:
    perturbLike(x.a, dx.a, eps, res.a, useExp)
    perturbLike(x.b, dx.b, eps, res.b, useExp)
  else:
    when useExp and compiles(addNoiseExp(eps, dx, x, res)):
      addNoiseExp(eps, dx, x, res)
    elif useExp and compiles(addNoiseExpField(eps, dx, x, res)):
      addNoiseExpField(eps, dx, x, res)
    elif compiles(perturbAdd(x, dx, eps, res)):
      perturbAdd(x, dx, eps, res)
    elif compiles(copyLike(res, x)):
      copyLike(res, x)
      when compiles(res += eps * dx):
        res += eps * dx
      else:
        res = res + eps * dx
    else:
      res = x + eps * dx


proc acceptTol*(cfg: TestConfig; err, fdErr, refVal: float): bool =
  let tol = max(max(cfg.atol, cfg.rtol * refVal), cfg.fdFactor * fdErr)
  err <= tol

proc reportFD*(cfg: TestConfig; sample: int; ana, fd, err, fdErr: float; ok: bool) =
  if ok:
    if cfg.verbose:
      echo "Test ", sample, " Passed: Ana=", ana, " FD=", fd,
           " err=", err, " ±", fdErr
  else:
    echo "Test ", sample, " FAILED: Ana=", ana, " FD=", fd,
         " err=", err, " ±", fdErr

proc reportSimple*(cfg: TestConfig; sample: int; lhs, rhs, err, rel: float; ok: bool) =
  if ok:
    if cfg.verbose:
      echo "Test ", sample, " Passed: LHS=", lhs, " RHS=", rhs,
           " err=", err, " rel=", rel
  else:
    echo "Test ", sample, " FAILED: LHS=", lhs, " RHS=", rhs,
         " err=", err, " rel=", rel

proc reportCompare*(
    cfg: TestConfig;
    sample: int;
    lhs, rhs: float;
    label = ""
  ): bool =
  if label.len > 0:
    echo "- ", label
  let err = abs(lhs - rhs)
  let refVal = max(1.0, max(abs(lhs), abs(rhs)))
  let rel = err / refVal
  let ok = err <= max(cfg.atol, cfg.rtol * refVal)
  reportSimple(cfg, sample, lhs, rhs, err, rel, ok)
  ok

proc captureFD*(
    fdStore: ptr float;
    fdeStore: ptr float
  ): AcceptFD =
  proc accept(
      cfg: TestConfig;
      ana, fd, err, fdErr, refVal: float;
      sample: int
    ): bool =
    if fdStore != nil:
      fdStore[] = fd
    if fdeStore != nil:
      fdeStore[] = fdErr
    acceptTol(cfg, err, fdErr, refVal)
  accept

proc randMat*[MatT](rng: var Rand; scale = 1.0): MatT =
  var m: MatT
  for i in 0..<m.nrows:
    for j in 0..<m.ncols:
      m[i,j].re = scale * (rng.rand(1.0) - 0.5)
      m[i,j].im = scale * (rng.rand(1.0) - 0.5)
  m

proc randTAH*[MatT](rng: var Rand; scale = 1.0): MatT =
  var m = randMat[MatT](rng, scale)
  m.projectTAH
  m

proc hermitianize*[MatT](m: MatT): MatT =
  0.5 * (m + m.adj)

proc perturbAdd*[MatT](x: MatT; dx: MatT; eps: float; res: var MatT) =
  res := x
  res += eps * dx

# Runs finite-difference vs JVP for general X->Y maps (custom perturb/dot).
# Check: d/d eps <p, f(x + eps*dx)>|0  ==  <p, J dx>.
proc runFDJVP*(
    cfg: TestConfig,
    makeBase: auto,
    makeDir: auto,
    makeProbe: auto,
    perturb: auto,
    f: auto,
    jvp: auto,
    accept: AcceptFD = nil,
    ordMax: static int = 6
  ): TestResult =
  reportBanner(cfg)
  result.name = cfg.name
  for sample in 0..<cfg.samples:
    let base = makeBase(sample)
    let dx = makeDir(sample, base)
    let probe = makeProbe(sample, base)
    var xPert = allocLike(base)
    var y = allocLike(probe)
    var dy = allocLike(probe)
    proc F(eps: float): float =
      perturb(base, dx, eps, xPert)
      f(xPert, y)
      dotLike(probe, y)
    var fd, fde: float
    ndiff(fd, fde, F, 0.0, cfg.eps, scale = cfg.scale, ordMax = ordMax)
    jvp(base, dx, dy)
    let ana = dotLike(probe, dy)
    let err = abs(ana - fd)
    let refVal = max(1.0, abs(fd))
    let ok = if accept != nil:
      accept(cfg, ana, fd, err, fde, refVal, sample)
    else:
      acceptTol(cfg, err, fde, refVal)
    reportFD(cfg, sample, ana, fd, err, fde, ok)
    if ok: inc result.passed else: inc result.failed

# Convenience JVP FD check for common cases (default cfg + exp-perturb).
# Check: d/d eps <p, f(x + eps*dx)>|0  ==  <p, J dx>.
proc runFDJVPDefault*(
    cfg: TestConfig,
    makeBase: auto,
    makeDir: auto,
    makeProbe: auto,
    f: auto,
    jvp: auto,
    accept: AcceptFD = nil,
    ordMax: static int = 6
  ): TestResult =
  type X = typeof(makeBase(0))
  type DX = typeof(makeDir(0, default(X)))
  type Y = typeof(makeProbe(0, default(X)))
  proc perturb(x: X; dx: DX; eps: float; res: var X) =
    perturbLike(x, dx, eps, res, true)
  runFDJVP(
    cfg = cfg,
    makeBase = makeBase,
    makeDir = makeDir,
    makeProbe = makeProbe,
    perturb = perturb,
    f = f,
    jvp = jvp,
    accept = accept,
    ordMax = ordMax
  )

# Runs finite-difference vs VJP for general X->Y maps (custom perturb/dot).
# Check: d/d eps <ybar, f(x + eps*dx)>|0  ==  <J^T ybar, dx>.
# Runs finite-difference vs JVP for scalar-valued f (custom perturb).
# Check: d/d eps f(x + eps*dx)|0  ==  J dx (directional derivative).
proc runFDScalarJVP*(
    cfg: TestConfig,
    makeBase: auto,
    makeDir: auto,
    perturb: auto,
    f: auto,
    jvp: auto,
    accept: AcceptFD = nil,
    ordMax: static int = 6
  ): TestResult =
  reportBanner(cfg)
  result.name = cfg.name
  for sample in 0..<cfg.samples:
    let base = makeBase(sample)
    let dx = makeDir(sample, base)
    var xPert = allocLike(base)
    proc F(eps: float): float =
      perturb(base, dx, eps, xPert)
      var val: float
      when compiles(f(xPert, val)):
        f(xPert, val)
        val
      else:
        f(xPert)
    var fd, fde: float
    ndiff(fd, fde, F, 0.0, cfg.eps, scale = cfg.scale, ordMax = ordMax)
    var ana: float
    when compiles(jvp(base, dx, ana)):
      jvp(base, dx, ana)
    else:
      ana = jvp(base, dx)
    let err = abs(ana - fd)
    let refVal = max(1.0, abs(fd))
    let ok = if accept != nil:
      accept(cfg, ana, fd, err, fde, refVal, sample)
    else:
      acceptTol(cfg, err, fde, refVal)
    reportFD(cfg, sample, ana, fd, err, fde, ok)
    if ok: inc result.passed else: inc result.failed

# Runs finite-difference vs JVP for scalar f with default exp-perturb.
# Check: d/d eps f(x + eps*dx)|0  ==  J dx (directional derivative).
proc runFDScalarJVPDefault*(
    cfg: TestConfig,
    makeBase: auto,
    makeDir: auto,
    f: auto,
    jvp: auto,
    accept: AcceptFD = nil,
    ordMax: static int = 6
  ): TestResult =
  type X = typeof(makeBase(0))
  type DX = typeof(makeDir(0, default(X)))
  proc perturb(x: X; dx: DX; eps: float; res: var X) =
    perturbLike(x, dx, eps, res, true)
  runFDScalarJVP(
    cfg = cfg,
    makeBase = makeBase,
    makeDir = makeDir,
    perturb = perturb,
    f = f,
    jvp = jvp,
    accept = accept,
    ordMax = ordMax
  )

macro jvpResType(jvp: typed): untyped =
  let t = jvp.getTypeInst
  let procTy = if t.kind == nnkProcTy: t else: t.getTypeInst
  let params = procTy[0]
  if params.len < 4:
    error("jvp must have signature (x, dx, res: var Y)", jvp)
  var resType = params[3][1]
  if resType.kind == nnkVarTy:
    resType = resType[0]
  resType

# Checks adjoint identity <ybar, J dx> = <J^T ybar, dx>, with ybar generated internally.
# Check: <ybar, J dx>  ==  <J^T ybar, dx>.
proc runAdjointIdentity*(
    cfg: TestConfig,
    makeBase: auto,
    makeDir: auto,
    jvp: auto,
    vjp: auto
  ): TestResult =
  type X = typeof(makeBase(0))
  type DX = typeof(makeDir(0, default(X)))
  type Y = jvpResType(jvp)
  proc makeYbar(sample: int; base: X): Y =
    let s1 = sample + cfg.samples + 1
    let s2 = sample + cfg.samples + 2
    when compiles(result = makeDir(s1, base)):
      result = makeDir(s1, base)
    elif compiles(result = makePair(makeDir(s1, base), makeDir(s2, base))):
      result = makePair(makeDir(s1, base), makeDir(s2, base))
    else:
      let dx2 = makeDir(s1, base)
      result = seedFromDir[Y, DX](dx2)
      jvp(base, dx2, result)
  reportBanner(cfg)
  result.name = cfg.name
  for sample in 0..<cfg.samples:
    let base = makeBase(sample)
    let dx = makeDir(sample, base)
    var dy = seedFromDir[Y, DX](dx)
    let ybar = makeYbar(sample, base)
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

# Checks linearity: A(a v1 + b v2) == a A(v1) + b A(v2).
# Check: A(a v1 + b v2)  ==  a A(v1) + b A(v2).
proc runLinearity*(
    cfg: TestConfig,
    makeVec: auto,
    makeProbe: auto,
    combine: auto,
    apply: auto,
    a: float = 0.7,
    b: float = -1.3
  ): TestResult =
  reportBanner(cfg)
  result.name = cfg.name
  for sample in 0..<cfg.samples:
    let v1 = makeVec(sample)
    let v2 = makeVec(sample)
    let probe = makeProbe(sample)
    var v12 = allocLike(v1)
    combine(v1, v2, a, b, v12)
    var y1 = allocLike(probe)
    var y2 = allocLike(probe)
    var y12 = allocLike(probe)
    apply(v1, y1)
    apply(v2, y2)
    apply(v12, y12)
    let lhs = dotLike(probe, y12)
    let rhs = a * dotLike(probe, y1) + b * dotLike(probe, y2)
    let err = abs(lhs - rhs)
    let refVal = max(1.0, max(abs(lhs), abs(rhs)))
    let rel = err / refVal
    let ok = err <= max(cfg.atol, cfg.rtol * refVal)
    reportSimple(cfg, sample, lhs, rhs, err, rel, ok)
    if ok: inc result.passed else: inc result.failed

# Runs a sweep over coefficient cases and aggregates pass/fail.
# Check: passed = sum_i passed_i, failed = sum_i failed_i.
proc runCoefSweep*[T](
    name: string;
    coefs: openArray[T];
    labels: openArray[string];
    runOne: proc(coef: T): TestResult
  ): TestResult =
  reportBanner(name)
  result.name = name
  var passedCombos: seq[int] = @[]
  var failedCombos: seq[int] = @[]
  for idx in 0..<coefs.len:
    let res = runOne(coefs[idx])
    result.add(res)
    let status = if res.failed == 0: "PASS" else: "FAIL"
    if idx < labels.len:
      echo "Combo ", idx, " ", labels[idx], ": ", status
    else:
      echo "Combo ", idx, ": ", status
    if res.failed == 0:
      passedCombos.add(idx)
    else:
      failedCombos.add(idx)
  echo ""
  echo "Summary: Passed=", passedCombos, " Failed=", failedCombos
  reportSummary(result)

proc initGaugeTestEnv*(
    defaultLat: seq[int] = @[8,8,8,8],
    defaultSeed: uint64 = 4321'u64,
    randomize = true,
    preSmear = true,
    preSmearSteps = 10,
    preSmearAlpha = 0.1
  ): auto =
  letParam:
    gaugefile = ""
    lat =
      if fileExists(gaugefile):
        getFileLattice gaugefile
      else:
        if gaugefile.len > 0:
          qexWarn "Nonexistent gauge file: ", gaugefile
        defaultLat
    seed = defaultSeed
  echoParams()
  let lo = lat.newLayout
  var g = lo.newGauge
  let hasFile = gaugefile.len > 0 and fileExists(gaugefile)
  if hasFile:
    let status = g.loadGauge(gaugefile)
    if status != 0:
      qexWarn "Failed to load gauge file: ", gaugefile
  elif randomize:
    g.random
    if preSmear and preSmearSteps > 0:
      var ss = lo.newStoutSmear(preSmearAlpha)
      for _ in 0..<preSmearSteps:
        ss.smear(g, g)
  g.echoPlaq
  var r = newRNGField(MRG32k3a, lo, seed)
  var env: GaugeTestEnv[type(lo), type(g), type(r)]
  env.lo = lo
  env.g = g
  env.r = r
  env.lat = lat
  env.gaugefile = gaugefile
  env.seed = seed
  env.hasFile = hasFile
  result = env

proc newRandomTAHDir*[L, R](lo: L; r: var R): auto =
  var p = lo.newGauge
  p.randomTAH r
  p


proc addNoiseExp*[G](x: float; p: G; g: G; ng: G) =
  qexGC()
  threads:
    for mu in 0..<g.len:
      for e in g[mu]:
        let t = x * p[mu][e]
        ng[mu][e] := exp(t) * g[mu][e]
  qexGC()

proc addNoiseExpField*[F](x: float; p: F; g: F; ng: F) =
  qexGC()
  threads:
    for i in g:
      ng[i] := exp(x * p[i]) * g[i]
  qexGC()

proc scaleForceField*[G](a: var G; scale: float) =
  ## Force fields live in the algebra but share the same container as gauge links.
  if scale == 1.0 or a.len == 0:
    return
  let p = cast[ptr cArray[type(a[0])]](unsafeAddr(a[0]))
  let len = a.len
  threads:
    for mu in 0..<len:
      for e in p[mu]:
        p[mu][e] *= scale

proc acceptActionForceDefault*(
    cfg: TestConfig;
    ana, fd, err, fdErr, refVal: float;
    sample: int
  ): bool =
  discard ana
  discard fd
  discard sample
  acceptTol(cfg, err, fdErr, refVal)

proc actionForceCfgDefault*(
    name: string,
    samples = defaultSamples,
    eps = actionForceDefaultEps,
    scale = defaultScale
  ): TestConfig =
  initTestConfig(
    name = name,
    samples = samples,
    eps = eps,
    scale = scale,
    atol = actionForceDefaultAtol,
    rtol = actionForceDefaultRtol,
    fdFactor = actionForceDefaultFdFactor,
    verbose = true
  )

# FD check for gauge action vs gauge force (VJP with scalar output).
# Check: d/d eps S(g + eps*dg)|0  ==  <dg, force(g)> scaled by ybar.
proc runGaugeActionForceFD*[L, G, R](
    env: GaugeTestEnv[L, G, R],
    action: proc(g: G): float,
    deriv: proc(g: G; f: G),
    cfg: TestConfig,
    accept: AcceptFD = acceptActionForceDefault
  ): TestResult =
  type GaugeField = G
  let lo = env.lo
  let gptr = addr env.g
  let rptr = addr env.r
  var rng = initRand(env.seed.int)
  reportBanner(cfg)
  result.name = cfg.name
  for sample in 0..<cfg.samples:
    let base = gptr[]
    let dx = newRandomTAHDir(lo, rptr[])
    let ybar = 1.0 + 0.1 * (2.0 * rng.rand(1.0) - 1.0)
    var xPert = allocLike(base)
    var y: float
    var dxbar = allocLike(dx)
    proc F(eps: float): float =
      addNoiseExp(eps, dx, base, xPert)
      y = action(xPert)
      dotLike(ybar, y)
    var fd, fde: float
    ndiff(fd, fde, F, 0.0, cfg.eps, scale = cfg.scale, ordMax = 6)
    deriv(base, dxbar)
    scaleForceField(dxbar, ybar)
    let ana = dotLike(dxbar, dx)
    let err = abs(ana - fd)
    let refVal = max(1.0, abs(fd))
    let ok = if accept != nil:
      accept(cfg, ana, fd, err, fde, refVal, sample)
    else:
      acceptTol(cfg, err, fde, refVal)
    reportFD(cfg, sample, ana, fd, err, fde, ok)
    if ok: inc result.passed else: inc result.failed

# FD check for gauge HVP using a gradient and hvp callback.
# Check: d/d eps <p, grad(g + eps*dg)>|0  ==  <p, H dg>.
proc runGaugeHVPFD*[L, G, R](
    env: var GaugeTestEnv[L, G, R],
    grad: proc(x: G; res: var G),
    hvp: proc(x: G; dx: G; res: var G),
    cfg: TestConfig,
    makeProbe: proc(sample: int; x: G; dx: G): G = nil
  ): TestResult =
  type GaugeField = G
  let lo = env.lo
  let gptr = addr env.g
  let rptr = addr env.r
  var sharedDir: GaugeField
  result = runFDJVP(
    cfg = cfg,
    makeBase = proc(sample: int): GaugeField = gptr[],
    makeDir = proc(sample: int; x: GaugeField): GaugeField =
      sharedDir = newRandomTAHDir(lo, rptr[])
      sharedDir,
    makeProbe = proc(sample: int; x: GaugeField): GaugeField =
      if makeProbe != nil:
        makeProbe(sample, x, sharedDir)
      else:
        sharedDir,
    perturb = proc(x: GaugeField; dx: GaugeField; eps: float; res: var GaugeField) =
      addNoiseExp(eps, dx, x, res),
    f = grad,
    jvp = hvp,
    ordMax = 6
  )

# FD check for gauge JVP with field-valued outputs.
# Check: d/d eps <p, f(g + eps*dg)>|0  ==  <p, J dg>.
proc runGaugeJVPFD*[L, G, R](
    env: var GaugeTestEnv[L, G, R],
    f: proc(x: G; res: var G),
    jvp: proc(x: G; dx: G; res: var G),
    cfg: TestConfig,
    makeProbe: proc(sample: int; x: G): G = nil,
    makeDir: proc(sample: int; x: G): G = nil,
    perturb: proc(x: G; dx: G; eps: float; res: var G) = nil
  ): TestResult =
  type GaugeField = G
  let lo = env.lo
  let gptr = addr env.g
  let rptr = addr env.r
  var sharedDir: GaugeField
  proc dirProc(sample: int; x: GaugeField): GaugeField =
    if makeDir != nil:
      sharedDir = makeDir(sample, x)
    else:
      sharedDir = newRandomTAHDir(lo, rptr[])
    sharedDir
  proc probeProc(sample: int; x: GaugeField): GaugeField =
    if makeProbe != nil:
      makeProbe(sample, x)
    else:
      sharedDir
  proc perturbProc(x: GaugeField; dx: GaugeField; eps: float; res: var GaugeField) =
    if perturb != nil:
      perturb(x, dx, eps, res)
    else:
      addNoiseExp(eps, dx, x, res)
  result = runFDJVP(
    cfg = cfg,
    makeBase = proc(sample: int): GaugeField = gptr[],
    makeDir = dirProc,
    makeProbe = probeProc,
    perturb = perturbProc,
    f = f,
    jvp = jvp,
    ordMax = 6
  )

# Adjoint identity check specialized to gauge fields.
# Check: <ybar, J dg>  ==  <J^T ybar, dg>.
proc runGaugeAdjointIdentity*[L, G, R](
    env: var GaugeTestEnv[L, G, R],
    jvp: proc(x: G; dx: G; res: var G),
    vjp: proc(x: G; ybar: G; res: var G),
    cfg: TestConfig,
    makeDir: proc(sample: int; x: G): G = nil
  ): TestResult =
  type GaugeField = G
  let lo = env.lo
  let gptr = addr env.g
  let rptr = addr env.r
  proc dirProc(sample: int; x: GaugeField): GaugeField =
    if makeDir != nil:
      makeDir(sample, x)
    else:
      newRandomTAHDir(lo, rptr[])
  result = runAdjointIdentity(
    cfg = cfg,
    makeBase = proc(sample: int): GaugeField = gptr[],
    makeDir = dirProc,
    jvp = jvp,
    vjp = vjp
  )
