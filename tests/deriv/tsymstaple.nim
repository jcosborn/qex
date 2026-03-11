import qex, physics/qcdTypes, gauge/smearutil, algorithms/numdiff
import maths/matrixFunctions
import ./utils

qexInit()

var env = initGaugeTestEnv(preSmear = false)
let lo = env.lo
var g = env.g
var r = env.r

proc setupSymStaple(mu, nu: int; rng: var auto): auto =
  var g0 = lo.newGauge
  g0.random
  let g1 = g0[nu]
  let g2 = g0[mu]
  var c = g1.newOneOf
  c.randomTAH rng
  let s1 = newShifter(g1, mu, 1)
  let s2 = newShifter(g2, nu, 1)
  let s  = newShifter(c, nu, 1)
  let sm1 = newShifter(g1, nu, -1)
  let sm2 = newShifter(g2, mu, -1)
  var tm1 = g1.newOneOf
  var tm2 = g1.newOneOf
  (g1: g1, g2: g2, c: c, s1: s1, s2: s2, s: s, sm1: sm1, sm2: sm2, tm1: tm1, tm2: tm2)

type H3[T] = Pair[T, Pair[T, T]]

template makeH3(a, b, c: untyped): untyped =
  makePair(a, makePair(b, c))

template h2(x: untyped): untyped = x.b.a
template h3(x: untyped): untyped = x.b.b

proc scopedCfg(name: string; base: TestConfig): TestConfig =
  initTestConfig(
    name = name,
    samples = base.samples,
    atol = base.atol,
    rtol = base.rtol,
    verbose = base.verbose
  )

proc testSymStapleJvpFD(): int =
  echo "\n=== symStapleJVP isolated test (multiplicative) ==="
  let mu = 0
  let nu = 1
  let alp = 0.5
  var g1 = lo.newGauge
  var g2 = lo.newGauge
  g1.random
  g2.random
  type Fld = type(g1[0])
  proc makeBase(sample: int): Pair[type(g1), type(g2)] =
    makePair(g1, g2)
  proc makeDir(sample: int; x: Pair[type(g1), type(g2)]): Pair[type(g1), type(g2)] =
    var d1 = lo.newGauge
    var d2 = lo.newGauge
    threads:
      for d in 0..<4:
        d1[d] := 0
        d2[d] := 0
    d1[mu].randomTAH env.r
    d2[nu].randomTAH env.r
    makePair(d1, d2)
  proc makeProbe(sample: int; x: Pair[type(g1), type(g2)]): Fld =
    var p: Fld
    p.new(lo)
    p.gaussian env.r
    p
  proc f(x: Pair[type(g1), type(g2)]; res: var Fld) =
    var s1 = newShifter(x.a[mu], nu, 1)
    var s2 = newShifter(x.b[nu], mu, 1)
    var sm = newShifter(x.a[mu], mu, -1)
    var tm: Fld
    tm.new(lo)
    let resPtr = cast[ptr Fld](unsafeAddr(res))
    threads:
      discard s1 ^*! x.a[mu]
      discard s2 ^*! x.b[nu]
      resPtr[] := 0
      symStaple(resPtr[], alp, x.a[mu], x.b[nu], s1, s2, tm, sm)
  proc jvp(x: Pair[type(g1), type(g2)]; dx: Pair[type(g1), type(g2)]; res: var Fld) =
    var s1 = newShifter(x.a[mu], nu, 1)
    var s2 = newShifter(x.b[nu], mu, 1)
    var sm = newShifter(x.a[mu], mu, -1)
    var tm: Fld
    tm.new(lo)
    var dg1mu: Fld
    var dg2nu: Fld
    dg1mu.new(lo)
    dg2nu.new(lo)
    threads:
      dg1mu := dx.a[mu] * x.a[mu]
      dg2nu := dx.b[nu] * x.b[nu]
    var ds1 = newShifter(dg1mu, nu, 1)
    var ds2 = newShifter(dg2nu, mu, 1)
    let resPtr = cast[ptr Fld](unsafeAddr(res))
    threads:
      discard s1 ^*! x.a[mu]
      discard s2 ^*! x.b[nu]
      discard ds1 ^*! dg1mu
      discard ds2 ^*! dg2nu
      resPtr[] := 0
      symStapleJVP(resPtr[], alp, x.a[mu], x.b[nu], s1, s2, dg1mu, dg2nu, ds1, ds2, tm, sm)
  let cfg = initTestConfig(
    name = "symStapleJVP isolated test (multiplicative)",
    samples = 1,
    eps = 0.01,
    scale = 2.0,
    atol = 0.1,
    rtol = 0.01,
    fdFactor = 3.0,
    verbose = true
  )
  let res = runFDJVPDefault(
    cfg = cfg,
    makeBase = makeBase,
    makeDir = makeDir,
    makeProbe = makeProbe,
    f = f,
    jvp = jvp
  )
  result = res.failed

proc testSymStapleHvpFD(): int =
  echo "\n=== symStapleHVP isolated test (multiplicative) ==="
  let mu = 0
  let nu = 1
  var g1 = lo.newGauge
  var g2 = lo.newGauge
  g1.random
  g2.random
  var c = lo.newGauge
  c.randomTAH env.r
  type Fld = type(g1[0])
  type X = Pair[type(g1), Pair[type(g2), type(c)]]
  type Y = Pair[Fld, Fld]
  proc makeBase(sample: int): X =
    makePair(g1, makePair(g2, c))
  proc makeDir(sample: int; x: X): X =
    var d1 = lo.newGauge
    var d2 = lo.newGauge
    var dc = lo.newGauge
    for d in 0..<4:
      d1[d] := 0
      d2[d] := 0
      dc[d] := 0
    d1[mu].randomTAH env.r
    d2[nu].randomTAH env.r
    dc[mu].randomTAH env.r
    makePair(d1, makePair(d2, dc))
  proc makeProbe(sample: int; x: X): Y =
    var p1, p2: Fld
    p1.new(lo)
    p2.new(lo)
    p1.gaussian env.r
    p2.gaussian env.r
    makePair(p1, p2)
  proc perturb(x: X; dx: X; eps: float; res: var X) =
    addNoiseExp(eps, dx.a, x.a, res.a)
    addNoiseExp(eps, dx.b.a, x.b.a, res.b.a)
    for d in 0..<4:
      res.b.b[d] := x.b.b[d]
    res.b.b[mu] := x.b.b[mu] + eps * dx.b.b[mu]
  proc grad(x: X; res: var Y) =
    var s1 = newShifter(x.a[mu], nu, 1)
    var s2 = newShifter(x.b.a[nu], mu, 1)
    var sc = newShifter(x.b.b[mu], mu, 1)
    var sm1 = newShifter(x.a[mu], mu, -1)
    var sm2 = newShifter(x.a[mu], nu, -1)
    var tm1, tm2: Fld
    tm1.new(lo)
    tm2.new(lo)
    let resPtr = cast[ptr Y](unsafeAddr(res))
    threads:
      discard s1 ^*! x.a[mu]
      discard s2 ^*! x.b.a[nu]
      discard sc ^*! x.b.b[mu]
      resPtr[].a := 0
      resPtr[].b := 0
      symStapleVJP(resPtr[].a, resPtr[].b, x.a[mu], x.b.a[nu], s1, s2, x.b.b[mu], sc, tm1, tm2, sm1, sm2)
  proc hvp(x: X; dx: X; res: var Y) =
    var dg1mu: Fld
    var dg2nu: Fld
    dg1mu.new(lo)
    dg2nu.new(lo)
    threads:
      dg1mu := dx.a[mu] * x.a[mu]
      dg2nu := dx.b.a[nu] * x.b.a[nu]
    var s1 = newShifter(x.a[mu], nu, 1)
    var s2 = newShifter(x.b.a[nu], mu, 1)
    var sc = newShifter(x.b.b[mu], mu, 1)
    var sm1 = newShifter(x.a[mu], mu, -1)
    var sm2 = newShifter(x.a[mu], nu, -1)
    var tm1, tm2: Fld
    tm1.new(lo)
    tm2.new(lo)
    var ds1 = newShifter(dg1mu, nu, 1)
    var ds2 = newShifter(dg2nu, mu, 1)
    var dsc = newShifter(dx.b.b[mu], mu, 1)
    var acc1, acc2: Fld
    acc1.new(lo)
    acc2.new(lo)
    let resPtr = cast[ptr Y](unsafeAddr(res))
    threads:
      discard s1 ^*! x.a[mu]
      discard s2 ^*! x.b.a[nu]
      discard sc ^*! x.b.b[mu]
      discard ds1 ^*! dg1mu
      discard ds2 ^*! dg2nu
      discard dsc ^*! dx.b.b[mu]
      resPtr[].a := 0
      resPtr[].b := 0
      symStapleHVP(resPtr[].a, resPtr[].b, acc1, acc2, x.a[mu], x.b.a[nu], s1, s2, x.b.b[mu], sc,
                   dg1mu, dg2nu, dx.b.b[mu], ds1, ds2, dsc, tm1, tm2, sm1, sm2)
  let cfg = initTestConfig(
    name = "symStapleHVP isolated test (multiplicative)",
    samples = 1,
    eps = 0.01,
    scale = 2.0,
    atol = 1e-2,
    rtol = 1e-3,
    fdFactor = 3.0,
    verbose = true
  )
  let res = runFDJVP(
    cfg = cfg,
    makeBase = makeBase,
    makeDir = makeDir,
    makeProbe = makeProbe,
    perturb = perturb,
    f = grad,
    jvp = hvp
  )
  result = res.failed

proc testSymStapleDeriv2Full(): int =
  var rng = lo.newRNGField(MRG32k3a, 1234)
  # Pick a fixed direction (μ,ν) to define g1,g2 and the corresponding shifters.
  let mu = 1
  let nu = 2
  let setup = setupSymStaple(mu, nu, rng)
  let g1 = setup.g1
  let g2 = setup.g2
  let c = setup.c
  let s1 = setup.s1
  let s2 = setup.s2
  let s = setup.s
  let sm1 = setup.sm1
  let sm2 = setup.sm2
  let tm1 = setup.tm1
  let tm2 = setup.tm2
  type
    Fld = type(g1)
    X = Pair[Fld, Pair[Fld, Fld]]
    Y = Pair[Fld, Fld]
  let cfg = initTestConfig(
    name = "tsymstaple2(link)",
    samples = 5,
    eps = 1e-3,
    scale = 2.0,
    atol = 1e-8,
    rtol = 1e-6,
    fdFactor = 5.0,
    verbose = true
  )
  let res = runFDJVP(
    cfg = cfg,
    makeBase = proc(sample: int): X =
      makePair(g1, makePair(g2, c)),
    makeDir = proc(sample: int; x: X): X =
      var dg1 = g1.newOneOf
      var dg2 = g2.newOneOf
      var dc = c.newOneOf
      dg1.randomTAH rng
      dg2.randomTAH rng
      dc.randomTAH rng
      makePair(dg1, makePair(dg2, dc)),
    makeProbe = proc(sample: int; x: X): Y =
      var probe1 = g1.newOneOf
      var probe2 = g2.newOneOf
      probe1.randomTAH rng
      probe2.randomTAH rng
      makePair(probe1, probe2),
    perturb = proc(x: X; dx: X; eps: float; res: var X) =
      addNoiseExpField(eps, dx.a, x.a, res.a)
      addNoiseExpField(eps, dx.b.a, x.b.a, res.b.a)
      res.b.b := x.b.b
      res.b.b += eps * dx.b.b,
    f = proc(x: X; res: var Y) =
      discard s1 ^*! x.a
      discard s2 ^*! x.b.a
      discard s  ^*! x.b.b
      res.a := 0
      res.b := 0
      symStapleVJP(res.a, res.b, x.a, x.b.a, s1, s2, x.b.b, s, tm1, tm2, sm1, sm2),
    jvp = proc(x: X; dx: X; res: var Y) =
      var dg1Eff = g1.newOneOf
      var dg2Eff = g2.newOneOf
      dg1Eff := dx.a * x.a
      dg2Eff := dx.b.a * x.b.a
      var ds1 = newShifter(dg1Eff, mu, 1)
      var ds2 = newShifter(dg2Eff, nu, 1)
      var dsc = newShifter(dx.b.b, nu, 1)
      discard s1 ^*! x.a
      discard s2 ^*! x.b.a
      discard s  ^*! x.b.b
      discard ds1 ^*! dg1Eff
      discard ds2 ^*! dg2Eff
      discard dsc ^*! dx.b.b
      var acc1 = g1.newOneOf
      var acc2 = g2.newOneOf
      res.a := 0
      res.b := 0
      symStapleHVP(res.a, res.b, acc1, acc2, x.a, x.b.a, s1, s2, x.b.b, s,
                   dg1Eff, dg2Eff, dx.b.b, ds1, ds2, dsc, tm1, tm2, sm1, sm2)
  )
  result = res.failed

proc testIsolatedG1(): int =
  ## Test with only dg1 non-zero (dg2=0, dc=0)
  ## This isolates the contribution from varying g1 and its shift s1.
  let mu = 1
  let nu = 2
  var rng = lo.newRNGField(MRG32k3a, 5555)
  let rptr = addr rng
  let setup = setupSymStaple(mu, nu, rng)
  type
    Fld = type(setup.g1)
    X = H3[Fld]
    Y = Pair[Fld, Fld]
  let g1 = setup.g1
  let g2 = setup.g2
  let c = setup.c
  let s1 = setup.s1
  let s2 = setup.s2
  let s = setup.s
  let sm1 = setup.sm1
  let sm2 = setup.sm2
  let tm1 = setup.tm1
  let tm2 = setup.tm2

  var zero2 = g2.newOneOf
  var zeroc = c.newOneOf
  zero2 := 0
  zeroc := 0
  var ds2_zero = newShifter(zero2, nu, 1)
  var dsc_zero = newShifter(zeroc, nu, 1)

  let cfg = initTestConfig(
    name = "testIsolatedG1",
    samples = 3,
    atol = 1e-8,
    rtol = 0.1,
    verbose = true
  )
  let res = runFDJVPDefault(
    cfg = cfg,
    makeBase = proc(sample: int): X =
      makeH3(g1, g2, c),
    makeDir = proc(sample: int; x: X): X =
      var dg1 = g1.newOneOf
      dg1.randomTAH rptr[]
      makeH3(dg1, zero2, zeroc),
    makeProbe = proc(sample: int; x: X): Y =
      var probe1 = g1.newOneOf
      var probe2 = g2.newOneOf
      probe1.randomTAH rptr[]
      probe2.randomTAH rptr[]
      makePair(probe1, probe2),
    f = proc(x: X; res: var Y) =
      discard s1 ^*! x.a
      discard s2 ^*! h2(x)
      discard s  ^*! h3(x)
      res.a := 0
      res.b := 0
      symStapleVJP(res.a, res.b, x.a, h2(x), s1, s2, h3(x), s, tm1, tm2, sm1, sm2),
    jvp = proc(x: X; dx: X; res: var Y) =
      var dg1Eff = x.a.newOneOf
      dg1Eff := dx.a * x.a
      var ds1 = newShifter(dg1Eff, mu, 1)
      discard s1 ^*! x.a
      discard s2 ^*! h2(x)
      discard s  ^*! h3(x)
      discard ds1 ^*! dg1Eff
      discard ds2_zero ^*! zero2
      discard dsc_zero ^*! zeroc
      res.a := 0
      res.b := 0
      var acc1 = res.a.newOneOf
      var acc2 = res.b.newOneOf
      symStapleHVP(res.a, res.b, acc1, acc2, x.a, h2(x), s1, s2, h3(x), s,
                   dg1Eff, zero2, zeroc, ds1, ds2_zero, dsc_zero, tm1, tm2, sm1, sm2)
  )
  result = res.failed

proc testIsolatedG2(): int =
  ## Test with only dg2 non-zero (dg1=0, dc=0)
  ## This isolates the contribution from varying g2 and its shift s2.
  let mu = 1
  let nu = 2
  var rng = lo.newRNGField(MRG32k3a, 6666)
  let rptr = addr rng
  let setup = setupSymStaple(mu, nu, rng)
  type
    Fld = type(setup.g1)
    X = H3[Fld]
    Y = Pair[Fld, Fld]
  let g1 = setup.g1
  let g2 = setup.g2
  let c = setup.c
  let s1 = setup.s1
  let s2 = setup.s2
  let s = setup.s
  let sm1 = setup.sm1
  let sm2 = setup.sm2
  let tm1 = setup.tm1
  let tm2 = setup.tm2

  var zero1 = g1.newOneOf
  var zeroc = c.newOneOf
  zero1 := 0
  zeroc := 0
  var ds1_zero = newShifter(zero1, mu, 1)
  var dsc_zero = newShifter(zeroc, nu, 1)

  let cfg = initTestConfig(
    name = "testIsolatedG2",
    samples = 3,
    atol = 1e-8,
    rtol = 0.1,
    verbose = true
  )
  let res = runFDJVPDefault(
    cfg = cfg,
    makeBase = proc(sample: int): X =
      makeH3(g1, g2, c),
    makeDir = proc(sample: int; x: X): X =
      var dg2 = g2.newOneOf
      dg2.randomTAH rptr[]
      makeH3(zero1, dg2, zeroc),
    makeProbe = proc(sample: int; x: X): Y =
      var probe1 = g1.newOneOf
      var probe2 = g2.newOneOf
      probe1.randomTAH rptr[]
      probe2.randomTAH rptr[]
      makePair(probe1, probe2),
    f = proc(x: X; res: var Y) =
      discard s1 ^*! x.a
      discard s2 ^*! h2(x)
      discard s  ^*! h3(x)
      res.a := 0
      res.b := 0
      symStapleVJP(res.a, res.b, x.a, h2(x), s1, s2, h3(x), s, tm1, tm2, sm1, sm2),
    jvp = proc(x: X; dx: X; res: var Y) =
      var dg2Eff = h2(x).newOneOf
      dg2Eff := h2(dx) * h2(x)
      var ds2 = newShifter(dg2Eff, nu, 1)
      discard s1 ^*! x.a
      discard s2 ^*! h2(x)
      discard s  ^*! h3(x)
      discard ds1_zero ^*! zero1
      discard ds2 ^*! dg2Eff
      discard dsc_zero ^*! zeroc
      res.a := 0
      res.b := 0
      var acc1 = res.a.newOneOf
      var acc2 = res.b.newOneOf
      symStapleHVP(res.a, res.b, acc1, acc2, x.a, h2(x), s1, s2, h3(x), s,
                   zero1, dg2Eff, zeroc, ds1_zero, ds2, dsc_zero, tm1, tm2, sm1, sm2)
  )
  result = res.failed

proc testIsolatedC(): int =
  ## Test with only dc non-zero (dg1=0, dg2=0)
  ## This isolates the contribution from varying c and its shift s.
  let mu = 1
  let nu = 2
  var rng = lo.newRNGField(MRG32k3a, 7777)
  let rptr = addr rng
  let setup = setupSymStaple(mu, nu, rng)
  type
    Fld = type(setup.g1)
    X = H3[Fld]
    Y = Pair[Fld, Fld]
  let g1 = setup.g1
  let g2 = setup.g2
  let c = setup.c
  let s1 = setup.s1
  let s2 = setup.s2
  let s = setup.s
  let sm1 = setup.sm1
  let sm2 = setup.sm2
  let tm1 = setup.tm1
  let tm2 = setup.tm2

  var zero1 = g1.newOneOf
  var zero2 = g2.newOneOf
  zero1 := 0
  zero2 := 0
  var ds1_zero = newShifter(zero1, mu, 1)
  var ds2_zero = newShifter(zero2, nu, 1)

  let cfg = initTestConfig(
    name = "testIsolatedC",
    samples = 3,
    atol = 1e-8,
    rtol = 0.1,
    verbose = true
  )
  let res = runFDJVP(
    cfg = cfg,
    makeBase = proc(sample: int): X =
      makeH3(g1, g2, c),
    makeDir = proc(sample: int; x: X): X =
      var dc = c.newOneOf
      dc.randomTAH rptr[]
      makeH3(zero1, zero2, dc),
    makeProbe = proc(sample: int; x: X): Y =
      var probe1 = g1.newOneOf
      var probe2 = g2.newOneOf
      probe1.randomTAH rptr[]
      probe2.randomTAH rptr[]
      makePair(probe1, probe2),
    perturb = proc(x: X; dx: X; eps: float; res: var X) =
      addNoiseExpField(eps, dx.a, x.a, res.a)
      addNoiseExpField(eps, h2(dx), h2(x), res.b.a)
      res.b.b := h3(x)
      res.b.b += eps * h3(dx),
    f = proc(x: X; res: var Y) =
      discard s1 ^*! x.a
      discard s2 ^*! h2(x)
      discard s  ^*! h3(x)
      res.a := 0
      res.b := 0
      symStapleVJP(res.a, res.b, x.a, h2(x), s1, s2, h3(x), s, tm1, tm2, sm1, sm2),
    jvp = proc(x: X; dx: X; res: var Y) =
      var dsc = newShifter(h3(dx), nu, 1)
      discard s1 ^*! x.a
      discard s2 ^*! h2(x)
      discard s  ^*! h3(x)
      discard ds1_zero ^*! zero1
      discard ds2_zero ^*! zero2
      discard dsc ^*! h3(dx)
      res.a := 0
      res.b := 0
      var acc1 = res.a.newOneOf
      var acc2 = res.b.newOneOf
      symStapleHVP(res.a, res.b, acc1, acc2, x.a, h2(x), s1, s2, h3(x), s,
                   zero1, zero2, h3(dx), ds1_zero, ds2_zero, dsc, tm1, tm2, sm1, sm2)
  )
  result = res.failed

proc testSymStapleDeriv(): int =
  ## Verify symStapleVJP is correct by comparing with FD of symStaple.
  ## This tests the first derivative, not the second.
  let mu = 1
  let nu = 2
  var rng = lo.newRNGField(MRG32k3a, 9999)
  let rptr = addr rng
  let setup = setupSymStaple(mu, nu, rng)
  type
    Fld = type(setup.g1)
    X = Pair[Fld, Fld]
  let g1 = setup.g1
  let g2 = setup.g2
  let c = setup.c
  let s1 = setup.s1
  let s2 = setup.s2
  let s = setup.s
  let sm1 = setup.sm1
  let sm2 = setup.sm2
  let tm1 = setup.tm1
  let sm = sm1

  let alp = 1.0
  let cfg = initTestConfig(
    name = "testSymStapleDeriv(g1)",
    samples = 3,
    atol = 1e-8,
    rtol = 0.1,
    verbose = true
  )

  let res = runFDScalarJVPDefault(
    cfg = cfg,
    makeBase = proc(sample: int): X =
      makePair(g1, g2),
    makeDir = proc(sample: int; x: X): X =
      var dg1 = g1.newOneOf
      dg1.randomTAH rptr[]
      var zero2 = g2.newOneOf
      zero2 := 0
      makePair(dg1, zero2),
    f = proc(x: X): float =
      var staple = g1.newOneOf
      staple := 0
      discard s1 ^*! x.a
      discard s2 ^*! x.b
      symStaple(staple, alp, x.a, x.b, s1, s2, tm1, sm)
      redot(c, staple),
    jvp = proc(x: X; dx: X): float =
      var dg1Eff = x.a.newOneOf
      dg1Eff := dx.a * x.a
      var f1 = g1.newOneOf
      var f2 = g2.newOneOf
      f1 := 0
      f2 := 0
      discard s1 ^*! x.a
      discard s2 ^*! x.b
      discard s  ^*! c
      symStapleVJP(f1, f2, x.a, x.b, s1, s2, c, s, tm1, tm1, sm1, sm2)
      redot(dg1Eff, f1)
  )
  result = res.failed

proc testSymStapleDerivAdjChain(): int =
  let cfg = initTestConfig(
    name = "symStapleVJPChain adjoint identity",
    samples = 3,
    atol = 1e-12,
    rtol = 1e-10,
    verbose = true
  )
  let mu = 0
  let nu = 1

  type Fld = type(g[0])

  var s1_nu_mu = newShifter(g[nu], mu, 1)
  var s1_mu_nu = newShifter(g[mu], nu, 1)
  var shiftedF2bar = newShifter(g[nu], nu, 1)
  var shiftedF1bar = newShifter(g[mu], mu, 1)
  var sm1_nu = newShifter(g[nu], nu, -1)
  var sm1_mu = newShifter(g[mu], mu, -1)

  var tm1, tm2: Fld
  tm1.new(lo)
  tm2.new(lo)

  let res = runAdjointIdentity(
    cfg = cfg,
    makeBase = proc(sample: int): int =
      0,
    makeDir = proc(sample: int; x: int): Fld =
      var c: Fld
      c.new(lo)
      threads:
        c.gaussian r
      c,
    jvp = proc(x: int; dx: Fld; res: var Pair[Fld, Fld]) =
      var sc = newShifter(dx, nu, 1)
      let resPtr = cast[ptr Pair[Fld, Fld]](unsafeAddr(res))
      threads:
        discard s1_nu_mu ^*! g[nu]
        discard s1_mu_nu ^*! g[mu]
        discard sc ^*! dx
        resPtr[].a := 0
        resPtr[].b := 0
        symStapleVJP(resPtr[].a, resPtr[].b, g[nu], g[mu], s1_nu_mu, s1_mu_nu, dx, sc, tm1, tm2, sm1_nu, sm1_mu),
    vjp = proc(x: int; ybar: Pair[Fld, Fld]; res: var Fld) =
      let resPtr = cast[ptr Fld](unsafeAddr(res))
      threads:
        discard s1_nu_mu ^*! g[nu]
        discard s1_mu_nu ^*! g[mu]
        resPtr[] := 0
        symStapleVJPChain(resPtr[], g[nu], g[mu], s1_nu_mu, s1_mu_nu,
                          ybar.a, ybar.b, shiftedF2bar, shiftedF1bar, tm1, tm2, sm1_nu),
  )
  result = res.failed

## Test each term of symStapleVJP individually
proc testSymStapleDerivTermByTerm(): int =
  echo "\n### Testing symStapleVJP terms individually ###"
  var fail = 0
  let cfgBase = initTestConfig(
    name = "symStapleVJP term",
    samples = 1,
    atol = 1e-10,
    rtol = 1e-10,
    verbose = true
  )

  let mu = 0
  let nu = 1

  type Fld = type(g[0])

  # Random inputs
  var c: Fld
  c.new(lo)
  threads:
    c.gaussian r
  proc makeBase(sample: int): int = 0
  proc makeDir(sample: int; x: int): Fld = c

  echo "\n  Testing Term 1: f1[x] += c[x] s1g[x] s2g[x]†"
  block:
    let res = runAdjointIdentity(
      cfg = cfgBase,
      makeBase = makeBase,
      makeDir = makeDir,
      jvp = proc(x: int; dx: Fld; res: var Fld) =
        var s1_nu_mu = newShifter(g[nu], mu, 1)
        var s1_mu_nu = newShifter(g[mu], nu, 1)
        var s1g, s2g: Fld
        s1g.new(lo)
        s2g.new(lo)
        let resPtr = cast[ptr Fld](unsafeAddr(res))
        threads:
          discard s1_nu_mu ^*! g[nu]
          discard s1_mu_nu ^*! g[mu]
          s1g := s1_nu_mu.field
          s2g := s1_mu_nu.field
          resPtr[] := dx * s1g * s2g.adj,
      vjp = proc(x: int; ybar: Fld; res: var Fld) =
        var s1_nu_mu = newShifter(g[nu], mu, 1)
        var s1_mu_nu = newShifter(g[mu], nu, 1)
        var s1g, s2g: Fld
        s1g.new(lo)
        s2g.new(lo)
        let resPtr = cast[ptr Fld](unsafeAddr(res))
        threads:
          discard s1_nu_mu ^*! g[nu]
          discard s1_mu_nu ^*! g[mu]
          s1g := s1_nu_mu.field
          s2g := s1_mu_nu.field
          resPtr[] := ybar * s2g * s1g.adj,
    )
    fail += res.failed

  echo "\n  Testing Term 2: f1[x] += g2[x] s1g[x] c[x-ν]†"
  block:
    let res = runAdjointIdentity(
      cfg = cfgBase,
      makeBase = makeBase,
      makeDir = makeDir,
      jvp = proc(x: int; dx: Fld; res: var Fld) =
        var s1_nu_mu = newShifter(g[nu], mu, 1)
        var sc = newShifter(dx, nu, 1)
        var s1g: Fld
        s1g.new(lo)
        let resPtr = cast[ptr Fld](unsafeAddr(res))
        threads:
          discard s1_nu_mu ^*! g[nu]
          discard sc ^*! dx
          s1g := s1_nu_mu.field
          resPtr[] := g[mu] * s1g * sc.field.adj,
      vjp = proc(x: int; ybar: Fld; res: var Fld) =
        var s1_nu_mu = newShifter(g[nu], mu, 1)
        var sm1_nu = newShifter(g[nu], nu, -1)
        var tm1, s1g: Fld
        tm1.new(lo)
        s1g.new(lo)
        let resPtr = cast[ptr Fld](unsafeAddr(res))
        threads:
          discard s1_nu_mu ^*! g[nu]
          s1g := s1_nu_mu.field
          tm1 := ybar.adj * g[mu] * s1g
          threadBarrier()
          discard sm1_nu ^*! tm1
          threadBarrier()
          resPtr[] := sm1_nu.field,
    )
    fail += res.failed

  echo "\n  Testing Term 3: f2[x] += g1[x] c[x-ν] s1g[x]†"
  block:
    let res = runAdjointIdentity(
      cfg = cfgBase,
      makeBase = makeBase,
      makeDir = makeDir,
      jvp = proc(x: int; dx: Fld; res: var Fld) =
        var s1_nu_mu = newShifter(g[nu], mu, 1)
        var sc = newShifter(dx, nu, 1)
        var s1g: Fld
        s1g.new(lo)
        let resPtr = cast[ptr Fld](unsafeAddr(res))
        threads:
          discard s1_nu_mu ^*! g[nu]
          discard sc ^*! dx
          s1g := s1_nu_mu.field
          resPtr[] := g[nu] * sc.field * s1g.adj,
      vjp = proc(x: int; ybar: Fld; res: var Fld) =
        var s1_nu_mu = newShifter(g[nu], mu, 1)
        var sm1_nu = newShifter(g[nu], nu, -1)
        var tm1, s1g: Fld
        tm1.new(lo)
        s1g.new(lo)
        let resPtr = cast[ptr Fld](unsafeAddr(res))
        threads:
          discard s1_nu_mu ^*! g[nu]
          s1g := s1_nu_mu.field
          tm1 := g[nu].adj * ybar * s1g
          threadBarrier()
          discard sm1_nu ^*! tm1
          threadBarrier()
          resPtr[] := sm1_nu.field,
    )
    fail += res.failed

  echo "\n  Testing Term 4: f1[x] += g2[x+μ]† g1[x+μ] c[x+μ-ν]"
  block:
    let res = runAdjointIdentity(
      cfg = cfgBase,
      makeBase = makeBase,
      makeDir = makeDir,
      jvp = proc(x: int; dx: Fld; res: var Fld) =
        var sc = newShifter(dx, nu, 1)
        var sm1_mu = newShifter(g[mu], mu, -1)
        var tm2: Fld
        tm2.new(lo)
        let resPtr = cast[ptr Fld](unsafeAddr(res))
        threads:
          discard sc ^*! dx
          tm2 := g[mu].adj * g[nu] * sc.field
          threadBarrier()
          discard sm1_mu ^*! tm2
          threadBarrier()
          resPtr[] := sm1_mu.field,
      vjp = proc(x: int; ybar: Fld; res: var Fld) =
        var s1_nu_mu = newShifter(g[nu], mu, 1)
        var sm1_nu = newShifter(g[nu], nu, -1)
        var tm1: Fld
        tm1.new(lo)
        let resPtr = cast[ptr Fld](unsafeAddr(res))
        threads:
          discard s1_nu_mu ^*! ybar
          threadBarrier()
          tm1 := g[nu].adj * g[mu] * s1_nu_mu.field
          threadBarrier()
          discard sm1_nu ^*! tm1
          threadBarrier()
          resPtr[] := sm1_nu.field,
    )
    fail += res.failed

  echo "\n  Testing Term 5: f1[x] += c[x+μ]† g1[x+μ] g2[x+μ-ν]"
  block:
    let res = runAdjointIdentity(
      cfg = cfgBase,
      makeBase = makeBase,
      makeDir = makeDir,
      jvp = proc(x: int; dx: Fld; res: var Fld) =
        var s1_mu_nu = newShifter(g[mu], nu, 1)
        var sm1_mu = newShifter(g[mu], mu, -1)
        var s2g, tm2: Fld
        s2g.new(lo)
        tm2.new(lo)
        let resPtr = cast[ptr Fld](unsafeAddr(res))
        threads:
          discard s1_mu_nu ^*! g[mu]
          s2g := s1_mu_nu.field
          tm2 := dx.adj * g[nu] * s2g
          threadBarrier()
          discard sm1_mu ^*! tm2
          threadBarrier()
          resPtr[] := sm1_mu.field,
      vjp = proc(x: int; ybar: Fld; res: var Fld) =
        var s1_nu_mu = newShifter(g[nu], mu, 1)
        var s1_mu_nu = newShifter(g[mu], nu, 1)
        var s2g: Fld
        s2g.new(lo)
        let resPtr = cast[ptr Fld](unsafeAddr(res))
        threads:
          discard s1_nu_mu ^*! ybar
          discard s1_mu_nu ^*! g[mu]
          s2g := s1_mu_nu.field
          resPtr[] := g[nu] * s2g * s1_nu_mu.field.adj,
    )
    fail += res.failed

  echo "\n  Testing Term 6: f2[x] += g1[x+ν]† c[x+ν] g1[x+ν-μ]"
  block:
    let res = runAdjointIdentity(
      cfg = cfgBase,
      makeBase = makeBase,
      makeDir = makeDir,
      jvp = proc(x: int; dx: Fld; res: var Fld) =
        var s1_nu_mu = newShifter(g[nu], mu, 1)
        var sm1_nu = newShifter(g[nu], nu, -1)
        var s1g, tm1: Fld
        s1g.new(lo)
        tm1.new(lo)
        let resPtr = cast[ptr Fld](unsafeAddr(res))
        threads:
          discard s1_nu_mu ^*! g[nu]
          s1g := s1_nu_mu.field
          tm1 := g[nu].adj * dx * s1g
          threadBarrier()
          discard sm1_nu ^*! tm1
          threadBarrier()
          resPtr[] := sm1_nu.field,
      vjp = proc(x: int; ybar: Fld; res: var Fld) =
        var s1_nu_mu = newShifter(g[nu], mu, 1)
        var s1_mu_nu = newShifter(g[mu], nu, 1)
        var s1g: Fld
        s1g.new(lo)
        let resPtr = cast[ptr Fld](unsafeAddr(res))
        threads:
          discard s1_mu_nu ^*! ybar
          discard s1_nu_mu ^*! g[nu]
          s1g := s1_nu_mu.field
          resPtr[] := g[nu] * s1_mu_nu.field * s1g.adj,
    )
    fail += res.failed

  echo "\n  Summary: ", fail, " terms failed"
  result = fail

## First, verify the basic adjoint formula for f = M c†
proc testBasicAdjoint(): int =
  type Fld = type(g[0])

  var M, c: Fld
  M.new(lo)
  c.new(lo)

  threads:
    M.gaussian r
    c.gaussian r
  let cfg = initTestConfig(
    name = "basic adjoint f = M c†",
    samples = 1,
    atol = 1e-12,
    rtol = 1e-10,
    verbose = true
  )
  let res = runAdjointIdentity(
    cfg = cfg,
    makeBase = proc(sample: int): int = 0,
    makeDir = proc(sample: int; x: int): Fld = c,
    jvp = proc(x: int; dx: Fld; res: var Fld) =
      let resPtr = cast[ptr Fld](unsafeAddr(res))
      threads:
        resPtr[] := M * dx.adj,
    vjp = proc(x: int; ybar: Fld; res: var Fld) =
      let resPtr = cast[ptr Fld](unsafeAddr(res))
      threads:
        resPtr[] := ybar.adj * M,
  )
  result = res.failed

## Test multiple (mu, nu) pairs

proc testSymStapleTangentAdjAllDirs(): int =
  echo "\n### Testing symStapleJVPVJP for all direction pairs ###"
  let cfgBase = initTestConfig(
    name = "symStapleJVPVJP adjoint identity (all dirs)",
    samples = 3,
    atol = 1e-12,
    rtol = 1e-10,
    verbose = false
  )
  var fail = 0
  let alp = 0.5

  type Fld = type(g[0])
  var tm1, tm2: Fld
  tm1.new(lo)
  tm2.new(lo)

  for mu in 0..<4:
    for nu in 0..<4:
      if nu == mu: continue

      var s1 = newShifter(g[nu], mu, 1)
      var s2 = newShifter(g[mu], nu, 1)
      var shiftedDsbar = newShifter(g[nu], nu, 1)
      var sm1 = newShifter(g[0], nu, -1)
      var sm2 = newShifter(g[0], mu, -1)

      let cfg = scopedCfg(
        "symStapleJVPVJP adjoint (mu=" & $mu & ", nu=" & $nu & ")",
        cfgBase
      )

      let res = runAdjointIdentity(
        cfg = cfg,
        makeBase = proc(sample: int): int =
          0,
        makeDir = proc(sample: int; x: int): Pair[Fld, Fld] =
          var dg1, dg2: Fld
          dg1.new(lo)
          dg2.new(lo)
          dg1.gaussian r
          dg2.gaussian r
          makePair(dg1, dg2),
        jvp = proc(x: int; dx: Pair[Fld, Fld]; res: var Fld) =
          var ds1 = newShifter(dx.a, mu, 1)
          var ds2 = newShifter(dx.b, nu, 1)
          discard s1 ^*! g[nu]
          discard s2 ^*! g[mu]
          discard ds1 ^*! dx.a
          discard ds2 ^*! dx.b
          res := 0
          symStapleJVP(res, alp, g[nu], g[mu], s1, s2, dx.a, dx.b, ds1, ds2, tm1, sm1),
        vjp = proc(x: int; ybar: Fld; res: var Pair[Fld, Fld]) =
          discard s1 ^*! g[nu]
          discard s2 ^*! g[mu]
          res.a := 0
          res.b := 0
          symStapleJVPVJP(res.a, res.b, alp, g[nu], g[mu], s1, s2,
                          ybar, shiftedDsbar, tm1, tm2, sm1, sm2),
      )
      if res.failed == 0:
        echo "  (mu=", mu, ", nu=", nu, ") Passed"
      else:
        echo "  (mu=", mu, ", nu=", nu, ") FAILED"
        inc fail

  result = fail

proc testSymStapleTangentAdjTerms(): int =
  echo "\n### Testing symStapleJVPVJP term adjoints (A–F) ###"
  let mu = 0
  let nu = 1
  let alp = 0.5

  type Fld = type(g[0])

  proc randFld(): Fld =
    var x: Fld
    x.new(lo)
    x.gaussian r
    x

  type TermDesc = object
    name: string
    jvp: proc(dx: Fld; res: var Fld)
    vjp: proc(ybar: Fld; res: var Fld)

  let terms = @[
    TermDesc(
      name: "Term A",
      jvp: proc(dx: Fld; res: var Fld) =
        var s1 = newShifter(g[nu], mu, 1)
        var sm = newShifter(g[0], nu, -1)
        var tm: Fld
        tm.new(lo)
        discard s1 ^*! g[nu]
        tm := dx.adj * g[mu] * s1.field
        discard sm ^*! tm
        res := alp * sm.field,
      vjp: proc(ybar: Fld; res: var Fld) =
        var sf = newShifter(g[0], nu, 1)
        var s1 = newShifter(g[nu], mu, 1)
        discard sf ^*! ybar
        discard s1 ^*! g[nu]
        res := alp * (g[mu] * s1.field * sf.field.adj)
    ),
    TermDesc(
      name: "Term B",
      jvp: proc(dx: Fld; res: var Fld) =
        var s1 = newShifter(g[nu], mu, 1)
        var s2 = newShifter(g[mu], nu, 1)
        discard s1 ^*! g[nu]
        discard s2 ^*! g[mu]
        res := alp * (dx * s2.field * s1.field.adj),
      vjp: proc(ybar: Fld; res: var Fld) =
        var s1 = newShifter(g[nu], mu, 1)
        var s2 = newShifter(g[mu], nu, 1)
        discard s1 ^*! g[nu]
        discard s2 ^*! g[mu]
        res := alp * (ybar * s1.field * s2.field.adj)
    ),
    TermDesc(
      name: "Term C",
      jvp: proc(dx: Fld; res: var Fld) =
        var s1 = newShifter(g[nu], mu, 1)
        var sm = newShifter(g[0], nu, -1)
        var tm: Fld
        tm.new(lo)
        discard s1 ^*! g[nu]
        tm := g[nu].adj * dx * s1.field
        discard sm ^*! tm
        res := alp * sm.field,
      vjp: proc(ybar: Fld; res: var Fld) =
        var sf = newShifter(g[0], nu, 1)
        var s1 = newShifter(g[nu], mu, 1)
        discard sf ^*! ybar
        discard s1 ^*! g[nu]
        res := alp * (g[nu] * sf.field * s1.field.adj)
    ),
    TermDesc(
      name: "Term D",
      jvp: proc(dx: Fld; res: var Fld) =
        var s1 = newShifter(g[nu], mu, 1)
        var sm1 = newShifter(g[0], nu, -1)
        var ds1 = newShifter(dx, mu, 1)
        var tm1: Fld
        tm1.new(lo)
        discard s1 ^*! g[nu]
        discard ds1 ^*! dx
        tm1 := g[nu].adj * g[mu] * ds1.field
        discard sm1 ^*! tm1
        res := alp * sm1.field,
      vjp: proc(ybar: Fld; res: var Fld) =
        var s2 = newShifter(g[mu], nu, 1)
        var sm2 = newShifter(g[0], mu, -1)
        var tm2: Fld
        tm2.new(lo)
        discard s2 ^*! ybar
        tm2 := g[mu].adj * g[nu] * s2.field
        discard sm2 ^*! tm2
        res := alp * sm2.field
    ),
    TermDesc(
      name: "Term E",
      jvp: proc(dx: Fld; res: var Fld) =
        var s1 = newShifter(g[nu], mu, 1)
        var ds2 = newShifter(dx, nu, 1)
        discard s1 ^*! g[nu]
        discard ds2 ^*! dx
        res := alp * (g[nu] * ds2.field * s1.field.adj),
      vjp: proc(ybar: Fld; res: var Fld) =
        var s1 = newShifter(g[nu], mu, 1)
        var sm1 = newShifter(g[0], nu, -1)
        var tm: Fld
        tm.new(lo)
        discard s1 ^*! g[nu]
        tm := g[nu].adj * ybar * s1.field
        discard sm1 ^*! tm
        res := alp * sm1.field
    ),
    TermDesc(
      name: "Term F",
      jvp: proc(dx: Fld; res: var Fld) =
        var s2 = newShifter(g[mu], nu, 1)
        var ds1 = newShifter(dx, mu, 1)
        discard s2 ^*! g[mu]
        discard ds1 ^*! dx
        res := alp * (g[nu] * s2.field * ds1.field.adj),
      vjp: proc(ybar: Fld; res: var Fld) =
        var s2 = newShifter(g[mu], nu, 1)
        var sm2 = newShifter(g[0], mu, -1)
        var tm: Fld
        tm.new(lo)
        discard s2 ^*! g[mu]
        tm := ybar.adj * g[nu] * s2.field
        discard sm2 ^*! tm
        res := alp * sm2.field
    )
  ]

  var fail = 0
  for term in terms:
    let cfg = initTestConfig(
      name = term.name,
      samples = 3,
      atol = 1e-12,
      rtol = 1e-10,
      verbose = true
    )
    let t = term
    let res = runAdjointIdentity(
      cfg = cfg,
      makeBase = proc(sample: int): int =
        0,
      makeDir = proc(sample: int; x: int): Fld =
        randFld(),
      jvp = proc(x: int; dx: Fld; res: var Fld) =
        t.jvp(dx, res),
      vjp = proc(x: int; ybar: Fld; res: var Fld) =
        t.vjp(ybar, res),
    )
    fail += res.failed

  ## Combined A+B+C manually
  block:
    echo "\n### Testing Combined A+B+C manually ###"
    var s1 = newShifter(g[nu], mu, 1)
    var s2 = newShifter(g[mu], nu, 1)
    var sm = newShifter(g[0], nu, -1)
    var sf = newShifter(g[0], nu, 1)
    var tm: Fld
    tm.new(lo)
    let cfg = initTestConfig(
      name = "Combined A+B+C",
      samples = 3,
      atol = 1e-12,
      rtol = 1e-10,
      verbose = true
    )
    let res = runAdjointIdentity(
      cfg = cfg,
      makeBase = proc(sample: int): int =
        0,
      makeDir = proc(sample: int; x: int): Pair[Fld, Fld] =
        makePair(randFld(), randFld()),
      jvp = proc(x: int; dx: Pair[Fld, Fld]; res: var Fld) =
        discard s1 ^*! g[nu]
        discard s2 ^*! g[mu]
        res := 0
        tm := dx.a.adj * g[mu] * s1.field
        discard sm ^*! tm
        res += alp * sm.field
        res += alp * (dx.a * s2.field * s1.field.adj)
        tm := g[nu].adj * dx.b * s1.field
        discard sm ^*! tm
        res += alp * sm.field,
      vjp = proc(x: int; ybar: Fld; res: var Pair[Fld, Fld]) =
        discard sf ^*! ybar
        discard s1 ^*! g[nu]
        discard s2 ^*! g[mu]
        res.a := 0
        res.b := 0
        res.a += alp * (g[mu] * s1.field * sf.field.adj)
        res.a += alp * (ybar * s1.field * s2.field.adj)
        res.b += alp * (g[nu] * sf.field * s1.field.adj),
    )
    fail += res.failed

  ## Using symStapleJVP and symStapleJVPVJP
  block:
    echo "\n### Testing with symStapleJVP and symStapleJVPVJP ###"
    var s1 = newShifter(g[nu], mu, 1)
    var s2 = newShifter(g[mu], nu, 1)
    var shiftedDsbar = newShifter(g[nu], nu, 1)
    var sm1 = newShifter(g[0], nu, -1)
    var sm2 = newShifter(g[0], mu, -1)
    var tm1, tm2: Fld
    tm1.new(lo)
    tm2.new(lo)
    let cfg = initTestConfig(
      name = "symStapleJVP/JVPVJP adjoint",
      samples = 3,
      atol = 1e-12,
      rtol = 1e-10,
      verbose = true
    )
    let res = runAdjointIdentity(
      cfg = cfg,
      makeBase = proc(sample: int): int =
        0,
      makeDir = proc(sample: int; x: int): Pair[Fld, Fld] =
        makePair(randFld(), randFld()),
      jvp = proc(x: int; dx: Pair[Fld, Fld]; res: var Fld) =
        var ds1 = newShifter(dx.a, mu, 1)
        var ds2 = newShifter(dx.b, nu, 1)
        discard s1 ^*! g[nu]
        discard s2 ^*! g[mu]
        discard ds1 ^*! dx.a
        discard ds2 ^*! dx.b
        res := 0
        symStapleJVP(res, alp, g[nu], g[mu], s1, s2, dx.a, dx.b, ds1, ds2, tm1, sm1),
      vjp = proc(x: int; ybar: Fld; res: var Pair[Fld, Fld]) =
        discard s1 ^*! g[nu]
        discard s2 ^*! g[mu]
        res.a := 0
        res.b := 0
        symStapleJVPVJP(res.a, res.b, alp, g[nu], g[mu], s1, s2,
                        ybar, shiftedDsbar, tm1, tm2, sm1, sm2),
    )
    fail += res.failed

  result = fail

## Test all (mu, nu) pairs

proc testSymStapleDeriv2AdjAllDirs(): int =
  echo "\n### Testing symStapleHVPVJP for all direction pairs ###"
  let cfgBase = initTestConfig(
    name = "symStapleHVPVJP adjoint identity (all dirs)",
    samples = 3,
    atol = 1e-12,
    rtol = 1e-10,
    verbose = false
  )
  var fail = 0

  type Fld = type(g[0])
  var tm1, tm2, acc1, acc2: Fld
  tm1.new(lo)
  tm2.new(lo)
  acc1.new(lo)
  acc2.new(lo)

  for mu in 0..<4:
    for nu in 0..<4:
      if nu == mu: continue

      var s1 = newShifter(g[nu], mu, 1)
      var s2 = newShifter(g[mu], nu, 1)
      var shiftedF2bar = newShifter(g[nu], nu, 1)
      var shiftedF1bar = newShifter(g[mu], mu, 1)
      var sm1 = newShifter(g[0], nu, -1)
      var sm2 = newShifter(g[0], mu, -1)

      var c: Fld
      c.new(lo)
      c.gaussian r
      var sc = newShifter(c, nu, 1)

      let cfg = scopedCfg(
        "symStapleHVPVJP adjoint (mu=" & $mu & ", nu=" & $nu & ")",
        cfgBase
      )

      let res = runAdjointIdentity(
        cfg = cfg,
        makeBase = proc(sample: int): int =
          0,
        makeDir = proc(sample: int; x: int): H3[Fld] =
          var dg1, dg2, dc: Fld
          dg1.new(lo)
          dg2.new(lo)
          dc.new(lo)
          dg1.gaussian r
          dg2.gaussian r
          dc.gaussian r
          makeH3(dg1, dg2, dc),
        jvp = proc(x: int; dx: H3[Fld]; res: var Pair[Fld, Fld]) =
          var ds1 = newShifter(dx.a, mu, 1)
          var ds2 = newShifter(h2(dx), nu, 1)
          var dsc = newShifter(h3(dx), nu, 1)
          discard s1 ^*! g[nu]
          discard s2 ^*! g[mu]
          discard sc ^*! c
          discard ds1 ^*! dx.a
          discard ds2 ^*! h2(dx)
          discard dsc ^*! h3(dx)
          res.a := 0
          res.b := 0
          acc1 := 0
          acc2 := 0
          symStapleHVP(res.a, res.b, acc1, acc2, g[nu], g[mu], s1, s2,
                       c, sc, dx.a, h2(dx), h3(dx), ds1, ds2, dsc, tm1, tm2, sm1, sm2),
        vjp = proc(x: int; ybar: Pair[Fld, Fld]; res: var H3[Fld]) =
          discard s1 ^*! g[nu]
          discard s2 ^*! g[mu]
          discard sc ^*! c
          res.a := 0
          h2(res) := 0
          h3(res) := 0
          symStapleHVPVJP(res.a, h2(res), h3(res), g[nu], g[mu], s1, s2,
                          c, sc, ybar.a, ybar.b, shiftedF2bar, shiftedF1bar, tm1, tm2, sm1, sm2),
      )
      if res.failed == 0:
        echo "  (mu=", mu, ", nu=", nu, ") Passed"
      else:
        echo "  (mu=", mu, ", nu=", nu, ") FAILED"
        inc fail

  result = fail

proc testSymStapleDeriv2AdjTerms(): int =
  echo "\n### Testing symStapleHVPVJP term adjoints (1–10) ###"
  let mu = 0
  let nu = 1
  type Fld = type(g[0])

  proc randFld(): Fld =
    var x: Fld
    x.new(lo)
    x.gaussian r
    x

  type TermDesc = object
    name: string
    jvp: proc(dx: Fld; res: var Fld)
    vjp: proc(ybar: Fld; res: var Fld)

  let terms = @[
    block:
      var s1 = newShifter(g[nu], mu, 1)
      var sc = newShifter(g[0], nu, 1)
      var c = randFld()
      TermDesc(
        name: "Term 1",
        jvp: proc(dx: Fld; res: var Fld) =
          discard s1 ^*! g[nu]
          discard sc ^*! c
          res := dx * s1.field * sc.field.adj,
        vjp: proc(ybar: Fld; res: var Fld) =
          discard s1 ^*! g[nu]
          discard sc ^*! c
          res := ybar * sc.field * s1.field.adj
      ),
    block:
      var sc = newShifter(g[0], nu, 1)
      var sm2 = newShifter(g[0], mu, -1)
      var c = randFld()
      var tm: Fld
      tm.new(lo)
      TermDesc(
        name: "Term 2",
        jvp: proc(dx: Fld; res: var Fld) =
          var ds1 = newShifter(dx, mu, 1)
          discard sc ^*! c
          discard ds1 ^*! dx
          res := g[mu] * ds1.field * sc.field.adj,
        vjp: proc(ybar: Fld; res: var Fld) =
          discard sc ^*! c
          tm := g[mu].adj * ybar * sc.field
          discard sm2 ^*! tm
          res := sm2.field
      ),
    block:
      var s1 = newShifter(g[nu], mu, 1)
      var sm1 = newShifter(g[0], nu, -1)
      var tm: Fld
      tm.new(lo)
      TermDesc(
        name: "Term 3",
        jvp: proc(dx: Fld; res: var Fld) =
          var dsc = newShifter(dx, nu, 1)
          discard s1 ^*! g[nu]
          discard dsc ^*! dx
          res := g[mu] * s1.field * dsc.field.adj,
        vjp: proc(ybar: Fld; res: var Fld) =
          discard s1 ^*! g[nu]
          tm := ybar.adj * g[mu] * s1.field
          discard sm1 ^*! tm
          res := sm1.field
      ),
    block:
      var s1 = newShifter(g[nu], mu, 1)
      var s2 = newShifter(g[mu], nu, 1)
      TermDesc(
        name: "Term 4",
        jvp: proc(dx: Fld; res: var Fld) =
          discard s1 ^*! g[nu]
          discard s2 ^*! g[mu]
          res := dx * s1.field * s2.field.adj,
        vjp: proc(ybar: Fld; res: var Fld) =
          discard s1 ^*! g[nu]
          discard s2 ^*! g[mu]
          res := ybar * s2.field * s1.field.adj
      ),
    block:
      var s2 = newShifter(g[mu], nu, 1)
      var sm2 = newShifter(g[0], mu, -1)
      var c = randFld()
      var tm: Fld
      tm.new(lo)
      TermDesc(
        name: "Term 5",
        jvp: proc(dx: Fld; res: var Fld) =
          var ds1 = newShifter(dx, mu, 1)
          discard s2 ^*! g[mu]
          discard ds1 ^*! dx
          res := c * ds1.field * s2.field.adj,
        vjp: proc(ybar: Fld; res: var Fld) =
          discard s2 ^*! g[mu]
          tm := c.adj * ybar * s2.field
          discard sm2 ^*! tm
          res := sm2.field
      ),
    block:
      var s1 = newShifter(g[nu], mu, 1)
      var sm1 = newShifter(g[0], nu, -1)
      var c = randFld()
      var tm: Fld
      tm.new(lo)
      TermDesc(
        name: "Term 6",
        jvp: proc(dx: Fld; res: var Fld) =
          var ds2 = newShifter(dx, nu, 1)
          discard s1 ^*! g[nu]
          discard ds2 ^*! dx
          res := c * s1.field * ds2.field.adj,
        vjp: proc(ybar: Fld; res: var Fld) =
          discard s1 ^*! g[nu]
          tm := ybar.adj * c * s1.field
          discard sm1 ^*! tm
          res := sm1.field
      ),
    block:
      var s1 = newShifter(g[nu], mu, 1)
      var sc = newShifter(g[0], nu, 1)
      var c = randFld()
      TermDesc(
        name: "Term 8",
        jvp: proc(dx: Fld; res: var Fld) =
          discard s1 ^*! g[nu]
          discard sc ^*! c
          res := dx * sc.field * s1.field.adj,
        vjp: proc(ybar: Fld; res: var Fld) =
          discard s1 ^*! g[nu]
          discard sc ^*! c
          res := ybar * s1.field * sc.field.adj
      ),
    block:
      var s1 = newShifter(g[nu], mu, 1)
      var sm1 = newShifter(g[0], nu, -1)
      var tm: Fld
      tm.new(lo)
      TermDesc(
        name: "Term 9",
        jvp: proc(dx: Fld; res: var Fld) =
          var dsc = newShifter(dx, nu, 1)
          discard s1 ^*! g[nu]
          discard dsc ^*! dx
          res := g[nu] * dsc.field * s1.field.adj,
        vjp: proc(ybar: Fld; res: var Fld) =
          discard s1 ^*! g[nu]
          tm := g[nu].adj * ybar * s1.field
          discard sm1 ^*! tm
          res := sm1.field
      ),
    block:
      var sc = newShifter(g[0], nu, 1)
      var sm2 = newShifter(g[0], mu, -1)
      var c = randFld()
      var tm: Fld
      tm.new(lo)
      TermDesc(
        name: "Term 10",
        jvp: proc(dx: Fld; res: var Fld) =
          var ds1 = newShifter(dx, mu, 1)
          discard sc ^*! c
          discard ds1 ^*! dx
          res := g[nu] * sc.field * ds1.field.adj,
        vjp: proc(ybar: Fld; res: var Fld) =
          discard sc ^*! c
          tm := ybar.adj * g[nu] * sc.field
          discard sm2 ^*! tm
          res := sm2.field
      )
  ]

  var fail = 0
  for term in terms:
    let cfg = initTestConfig(
      name = term.name,
      samples = 3,
      atol = 1e-12,
      rtol = 1e-10,
      verbose = true
    )
    let t = term
    let res = runAdjointIdentity(
      cfg = cfg,
      makeBase = proc(sample: int): int =
        0,
      makeDir = proc(sample: int; x: int): Fld =
        randFld(),
      jvp = proc(x: int; dx: Fld; res: var Fld) =
        t.jvp(dx, res),
      vjp = proc(x: int; ybar: Fld; res: var Fld) =
        t.vjp(ybar, res),
    )
    fail += res.failed

  result = fail

var fail = 0
# Section A: baseline FD
fail += testSymStapleJvpFD()
fail += testSymStapleHvpFD()
# Section B: JVP/HVP FD
fail += testSymStapleDeriv2Full()
fail += testIsolatedG1()
fail += testIsolatedG2()
fail += testIsolatedC()
fail += testSymStapleDeriv()
# Section C: adjoint identities
fail += testSymStapleDerivAdjChain()
fail += testSymStapleTangentAdjAllDirs()
fail += testSymStapleDeriv2AdjAllDirs()
# Section D: term-by-term checks
fail += testSymStapleDerivTermByTerm()
fail += testSymStapleTangentAdjTerms()
fail += testSymStapleDeriv2AdjTerms()
# Section E: diagnostics
fail += testBasicAdjoint()

if fail == 0:
  qexFinalize()
else:
  qexAbort(fail)
