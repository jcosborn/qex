import qex, physics/qcdTypes, gauge/gaugeAction, gauge/gaugeUtils
import ./utils

qexInit()

var env = initGaugeTestEnv()
let lo = env.lo

type
  GaugeField = typeof(env.g)
  ContractProjectTAHInput = Pair[GaugeField, GaugeField]
  ContractProjectTAHDir = Pair[GaugeField, GaugeField]
  ContractProjectTAHScalarBase = Pair[ContractProjectTAHInput, GaugeField]

let gc = GaugeActionCoeffs(plaq:6.0)
echo "S(g): ", gc.gaugeAction1(env.g)

proc randomLinearField(): GaugeField =
  var x = lo.newGauge
  x.gaussian env.r
  x

proc addLinearPerturb(dst: var GaugeField; x, dx: GaugeField; eps: float) =
  let d = cast[ptr cArray[type(dst[0])]](unsafeAddr(dst[0]))
  let xx = cast[ptr cArray[type(x[0])]](unsafeAddr(x[0]))
  let dd = cast[ptr cArray[type(dx[0])]](unsafeAddr(dx[0]))
  let nd = dst.len
  threads:
    for mu in 0..<nd:
      d[mu] := xx[mu]
      d[mu] += eps * dd[mu]

proc perturbContractProjectTAH(
    x: ContractProjectTAHInput;
    dx: ContractProjectTAHDir;
    eps: float;
    res: var ContractProjectTAHInput
  ) =
  addLinearPerturb(res.a, x.a, dx.a, eps)
  addLinearPerturb(res.b, x.b, dx.b, eps)

proc perturbContractProjectTAHScalar(
    x: ContractProjectTAHScalarBase;
    dx: ContractProjectTAHDir;
    eps: float;
    res: var ContractProjectTAHScalarBase
  ) =
  perturbContractProjectTAH(x.a, dx, eps, res.a)
  copyLike(res.b, x.b)

proc contractProjectTAHForward(res: var GaugeField; x: ContractProjectTAHInput) =
  copyLike(res, x.b)
  contractProjectTAH(x.a, res)

proc contractProjectTAHScalarForward(x: ContractProjectTAHScalarBase): float =
  var y = lo.newGauge
  contractProjectTAHForward(y, x.a)
  dotLike(x.b, y)

proc contractProjectTAHJVP(
    x: ContractProjectTAHInput;
    dx: ContractProjectTAHDir;
    res: var GaugeField
  ) =
  mixin adj
  let g = cast[ptr cArray[type(x.a[0])]](unsafeAddr(x.a[0]))
  let d = cast[ptr cArray[type(x.b[0])]](unsafeAddr(x.b[0]))
  let dg = cast[ptr cArray[type(dx.a[0])]](unsafeAddr(dx.a[0]))
  let dd = cast[ptr cArray[type(dx.b[0])]](unsafeAddr(dx.b[0]))
  let r = cast[ptr cArray[type(res[0])]](unsafeAddr(res[0]))
  let nd = res.len
  threads:
    for mu in 0..<nd:
      for e in r[mu]:
        var s {.noInit.}: type(r[mu][e])
        s := dg[mu][e] * d[mu][e].adj
        s += g[mu][e] * dd[mu][e].adj
        r[mu][e].projectTAH s

proc contractProjectTAHAdjoint(
    x: ContractProjectTAHInput;
    ybar: GaugeField;
    res: var ContractProjectTAHDir
  ) =
  contractProjectTAHVJP(x.a, x.b, ybar, res.b, res.a)

proc contractProjectTAHScalarJVP(
    x: ContractProjectTAHScalarBase;
    dx: ContractProjectTAHDir;
    res: var float
  ) =
  var dbar = lo.newGauge
  var gbar = lo.newGauge
  contractProjectTAHVJP(x.a.a, x.a.b, x.b, dbar, gbar)
  res = dotLike(gbar, dx.a) + dotLike(dbar, dx.b)

proc randomContractProjectTAHInput(): ContractProjectTAHInput =
  makePair(env.g, randomLinearField())

proc randomContractProjectTAHDir(): ContractProjectTAHDir =
  makePair(randomLinearField(), randomLinearField())

proc randomContractProjectTAHScalarBase(): ContractProjectTAHScalarBase =
  makePair(randomContractProjectTAHInput(), randomLinearField())

proc testContractProjectTAHVJPFD(): int =
  let cfg = initTestConfig(
    name = "contractProjectTAHVJP FD",
    samples = 5,
    eps = 1e-2,
    scale = 2.0,
    atol = 1e-7,
    rtol = 1e-6,
    fdFactor = 3.0,
    verbose = true
  )
  let res = runFDScalarJVP(
    cfg = cfg,
    makeBase = proc(sample: int): ContractProjectTAHScalarBase =
      randomContractProjectTAHScalarBase(),
    makeDir = proc(sample: int; x: ContractProjectTAHScalarBase): ContractProjectTAHDir =
      randomContractProjectTAHDir(),
    perturb = perturbContractProjectTAHScalar,
    f = contractProjectTAHScalarForward,
    jvp = contractProjectTAHScalarJVP
  )
  result = res.failed

proc testContractProjectTAHVJPAdjoint(): int =
  let cfg = initTestConfig(
    name = "contractProjectTAHVJP adjoint identity",
    samples = 5,
    atol = 1e-12,
    rtol = 1e-10,
    verbose = true
  )
  let res = runAdjointIdentity(
    cfg = cfg,
    makeBase = proc(sample: int): ContractProjectTAHInput =
      randomContractProjectTAHInput(),
    makeDir = proc(sample: int; x: ContractProjectTAHInput): ContractProjectTAHDir =
      randomContractProjectTAHDir(),
    jvp = contractProjectTAHJVP,
    vjp = contractProjectTAHAdjoint
  )
  result = res.failed

var fail = 0
# Section A: baseline FD
fail += runGaugeActionForceFD(
  env,
  action = proc(g: GaugeField): float =
    gc.gaugeAction1(g),
  deriv = proc(g: GaugeField; f: GaugeField) =
    gc.gaugeForce(g, f),
  cfg = actionForceCfgDefault("gaugeAction1/gaugeForce")
).failed
fail += testContractProjectTAHVJPFD()
fail += testContractProjectTAHVJPAdjoint()

if fail==0:
  qexFinalize()
else:
  qexAbort(fail)
