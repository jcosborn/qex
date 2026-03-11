import qex, physics/qcdTypes, gauge/stoutsmear
import ./utils

qexInit()

var env = initGaugeTestEnv()
let lo = env.lo
var g = env.g

type
  GaugeField = typeof(g)

var ss = lo.newStoutSmear(0.1)

let gc = GaugeActionCoeffs(plaq:6.0)
echo "S(g): ",gc.gaugeAction1(g)

var fail = 0

proc smearTest0[G](ss:var StoutSmear, gf:G, fl:G) =
  const nc = gf[0][0].nrows.float
  let
    alpha = -ss.alpha*nc  # negative from gaugeForce, and nc compensate force normalization
    f = ss.f
    ds = ss.ds
  ss.gf = gf
  gaugeActionDeriv(GaugeActionCoeffs(plaq:1.0), gf, ds)
  threads:
    for mu in 0..<f.len:
      for e in f[mu]:
        let s = gf[mu][e]*ds[mu][e].adj
        var t{.noinit.}: evalType(f[mu][e])
        t.projectTAH s
        f[mu][e] := t
        fl[mu][e] := exp(alpha*t)

proc smearTest0Deriv[G](ss:StoutSmear, deriv:G, chain:G) =
  const nc = chain[0][0].nrows.float
  let
    alpha = -ss.alpha*nc  # negative from gaugeForce, and nc compensate force normalization
    f = ss.f

  threads:
    for mu in 0..<f.len:
      for e in deriv[mu]:
        deriv[mu][e] := alpha*expDeriv(alpha*f[mu][e], chain[mu][e])
  gaugeForceDeriv(ss.gf, deriv, deriv, ss.ds, ss.cg)

proc smearedActionTest0(g:auto):auto =
  var sg = lo.newGauge
  ss.smearTest0(g, sg)
  gc.gaugeAction1(sg)
proc smearedForceTest0(g:auto, f:auto) =
  var sg = lo.newGauge
  var ds = lo.newGauge
  ss.smearTest0(g, sg)
  gc.gaugeActionDeriv(sg, ds)
  ss.smearTest0Deriv(f, ds)
  contractProjectTAH(g, f)


proc smearedAction(g:auto):auto =
  var sg = lo.newGauge
  ss.smear(g, sg)
  gc.gaugeAction1(sg)
proc smearedForce(g:auto, f:auto) =
  var sg = lo.newGauge
  var ds = lo.newGauge
  ss.smear(g, sg)
  gc.gaugeActionDeriv(sg, ds)
  ss.smearDeriv(f, ds)
  contractProjectTAH(g, f)


var s2 = lo.newStoutSmear(0.09)

proc smeared2Action(g:auto):auto =
  var sg = lo.newGauge
  var s2g = lo.newGauge
  ss.smear(g, sg)
  s2.smear(sg, s2g)
  gc.gaugeAction1(s2g)
proc smeared2Force(g:auto, f:auto) =
  var sg = lo.newGauge
  var ds = lo.newGauge
  var s2g = lo.newGauge
  var f2 = lo.newGauge
  ss.smear(g, sg)
  s2.smear(sg, s2g)
  gc.gaugeActionDeriv(s2g, ds)
  s2.smearDeriv(f2, ds)
  ss.smearDeriv(f, f2)
  contractProjectTAH(g, f)


var s3 = lo.newStoutSmear(0.12)

proc smeared3Action(g:auto):auto =
  var sg = lo.newGauge
  var s2g = lo.newGauge
  var s3g = lo.newGauge
  ss.smear(g, sg)
  s2.smear(sg, s2g)
  s3.smear(s2g, s3g)
  gc.gaugeAction1(s3g)
proc smeared3Force(g:auto, f:auto) =
  var sg = lo.newGauge
  var ds = lo.newGauge
  var s2g = lo.newGauge
  var s3g = lo.newGauge
  var f2 = lo.newGauge
  var f3 = lo.newGauge
  ss.smear(g, sg)
  s2.smear(sg, s2g)
  s3.smear(s2g, s3g)
  gc.gaugeActionDeriv(s3g, ds)
  s3.smearDeriv(f3, ds)
  s2.smearDeriv(f2, f3)
  ss.smearDeriv(f, f2)
  contractProjectTAH(g, f)

type
  ActionProc = proc(g: GaugeField): float
  ForceProc = proc(g: GaugeField; f: GaugeField)

let pipelines: seq[(string, ActionProc, ForceProc)] = @[
  ("smearedActionTest0/smearedForceTest0", smearedActionTest0, smearedForceTest0),
  ("smearedAction/smearedForce", smearedAction, smearedForce),
  ("smeared2Action/smeared2Force", smeared2Action, smeared2Force),
  ("smeared3Action/smeared3Force", smeared3Action, smeared3Force),
]

# Section A: baseline FD
for (name, action, force) in pipelines:
  fail += runGaugeActionForceFD(
    env,
    action, force,
    actionForceCfgDefault(name)
  ).failed

# echoTimers()

if fail==0:
  qexFinalize()
else:
  qexAbort(fail)
