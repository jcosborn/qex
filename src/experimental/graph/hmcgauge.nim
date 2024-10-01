import qex
import core, scalar, gauge
from os import fileExists
from strformat import `&`
from math import cos, PI

proc newOneOf(x: float): float = 0.0

type AdamW[Param] = object
  alpha, beta1, beta2, eps, lambda: float
  m: Param
  v: Param

func newAdamW[Param](param: Param, alpha = 0.001, beta1 = 0.9, beta2 = 0.999, eps = 1e-8, lambda = 0.01): AdamW[Param] =
  result = AdamW[Param](alpha: alpha, beta1: beta1, beta2: beta2, eps: eps, lambda: lambda)
  result.m = newOneOf param
  result.v = newOneOf param

func newAdam[Param](param: Param, alpha = 0.001, beta1 = 0.9, beta2 = 0.999, eps = 1e-8): AdamW[Param] =
  newAdamW[Param](param, alpha, beta1, beta2, eps, lambda = 0.0)

proc optimize[Param](opt: var AdamW[Param], param: var Param, grad: Param, t: int, lr: float)  =
  ## arXiv:1711.05101, standard Adam if lambda == 0, effectively scale grad by alpha/stdev(grad) for descent.
  ## Decay term is equivalent to an additional term of lambda/(2 scale) param^2 ~ lambda/(2 alpha) stdev(grad) param^2, added to the objective.
  ## Normalized weight decay suggests lambda = lambda_norm sqrt(b/BT), for batch size b, training number B, total epoch T.
  let
    a = opt.alpha
    b1 = opt.beta1
    b2 = opt.beta2
    sb1 = 1.0 - b1
    sb2 = 1.0 - b2
    sb1t = 1.0 - b1^t
    sb2t = 1.0 - b2^t
    dr = opt.lambda
    eps = opt.eps
  for i in 0..<grad.len:
    opt.m[i] = b1 * opt.m[i] + sb1 * grad[i]
    opt.v[i] = b2 * opt.v[i] + sb2 * (grad[i] * grad[i])
    let m = opt.m[i] / sb1t
    let v = sqrt(opt.v[i] / sb2t)
    if dr == 0:
      echo "Adam: ",i," ",m," ",v
      param[i] = param[i] - lr * (a * m / (v + eps))
    else:
      echo "AdamW: ",i," ",m," ",v
      param[i] = param[i] - lr * (a * m / (v + eps) + dr * param[i])

func warmUpCosDecay(t, twarm, tmax: int, lrmax: float, lrmin = 0.0): float =
  if t <= twarm:
    return lrmin + (lrmax - lrmin) * t.float / twarm.float
  else:
    return lrmin + 0.5 * (lrmax - lrmin) * (1.0 + cos(PI * float(t - twarm) / float(tmax - twarm)))

proc reunit(g:auto) =
  tic()
  threads:
    let d = g.checkSU
    threadBarrier()
    echo "unitary deviation avg: ",d.avg," max: ",d.max
    g.projectSU
    threadBarrier()
    let dd = g.checkSU
    echo "new unitary deviation avg: ",dd.avg," max: ",dd.max
  toc("reunit")

proc get(coeffs: openarray[float], i: int, d: float): Gvalue =
  toGvalue(if coeffs.len <= i or coeffs[i] == 0.0: d else: coeffs[i])

proc int2MN(gc, g0, p0, dt: Gvalue, n: int, coeffs: openarray[float]): (Gvalue, Gvalue, seq[Gvalue]) =
  let lambda = coeffs.get(0, 0.1931833275037836)
  var g = g0
  var p = p0
  let h = 0.5*dt
  let t05 = lambda*dt
  let t0 = 2.0*t05
  let t1 = dt-t0
  g = axexpmuly(t05, p, g)
  for i in 0..<n:
    if i>0:
      g = axexpmuly(t0, p, g)
    p = p - h * gaugeForce(gc, g)
    g = axexpmuly(t1, p, g)
    p = p - h * gaugeForce(gc, g)
  g = axexpmuly(t05, p, g)
  (g, p, @[lambda])

proc int4MN3F1GP(gc, g0, p0, dt: Gvalue, n: int, coeffs: openarray[float]): (Gvalue, Gvalue, seq[Gvalue]) =
  let lambda = coeffs.get(0, 0.2470939580390842)
  let theta = coeffs.get(1, 0.5 - 1.0 / sqrt(24.0 * lambda.getfloat))
  # scale the force gradient coeff to about the same order as the other
  let chi = coeffs.get(2, (1.0 - sqrt(6.0 * lambda.getfloat) * (1.0 - lambda.getfloat)) / 12.0 * (2.0 / (1.0 - 2.0*lambda.getfloat) * 10.0))
  var g = g0
  var p = p0
  let a0 = theta*dt
  let a02 = 2.0*a0
  let a1 = 0.5*dt - a0
  let b0 = lambda*dt
  let b1 = dt - 2.0*b0
  let c1 = 0.1*chi*(dt*dt)
  g = axexpmuly(a0, p, g)
  for i in 0..<n:
    if i>0:
      g = axexpmuly(a02, p, g)
    p = p - b0 * gaugeForce(gc, g)
    g = axexpmuly(a1, p, g)
    p = p - b1 * gaugeForce(gc, axexpmuly(-c1, gaugeForce(gc, g), g))
    g = axexpmuly(a1, p, g)
    p = p - b0 * gaugeForce(gc, g)
  g = axexpmuly(a0, p, g)
  (g, p, @[lambda, theta, chi])

proc int4MN5F2GP(gc, g0, p0, dt: Gvalue, n: int, coeffs: openarray[float]): (Gvalue, Gvalue, seq[Gvalue]) =
  let rho = coeffs.get(0, 0.06419108866816235)
  let theta = coeffs.get(1, 0.1919807940455741)
  let vtheta = coeffs.get(2, 0.1518179640276466)
  let lambda = coeffs.get(3, 0.2158369476787619)
  # scale the force gradient coeff to about the same order as the other
  let xi = coeffs.get(4, 0.0009628905212024874 * (2.0 / lambda.getfloat * 20.0))
  var g = g0
  var p = p0
  let a0 = rho*dt
  let a02 = 2.0*a0
  let a1 = theta*dt
  let a2 = (0.5-(theta+rho))*dt
  let b1 = lambda*dt
  let b0 = vtheta*dt
  let b2 = (1.0-2.0*(lambda+vtheta))*dt
  let c1 = 0.05*xi*(dt*dt)
  g = axexpmuly(a0, p, g)
  for i in 0..<n:
    if i>0:
      g = axexpmuly(a02, p, g)
    p = p - b0 * gaugeForce(gc, g)
    g = axexpmuly(a1, p, g)
    p = p - b1 * gaugeForce(gc, axexpmuly(-c1, gaugeForce(gc, g), g))
    g = axexpmuly(a2, p, g)
    p = p - b2 * gaugeForce(gc, g)
    g = axexpmuly(a2, p, g)
    p = p - b1 * gaugeForce(gc, axexpmuly(-c1, gaugeForce(gc, g), g))
    g = axexpmuly(a1, p, g)
    p = p - b0 * gaugeForce(gc, g)
  g = axexpmuly(a0, p, g)
  (g, p, @[rho, theta, vtheta, lambda, xi])

qexInit()

tic()

letParam:
  gaugefile = ""
  savefile = "config"
  savefreq = 0
  lat =
    if fileExists(gaugefile):
      getFileLattice gaugefile
    else:
      if gaugefile.len > 0:
        qexWarn "Nonexistent gauge file: ", gaugefile
      @[8,8,8,16]
  beta = 5.4
  dt = 0.025
  trajsThermo = 0
  trajsTrain = 50
  trajsTrainlrWarm = 10
  trajsInfer = 0
  lrmax = 1.0
  lrmin = 0.0001
  weightDecay = 0.0
  seed:uint = 1234567891
  gintalg = "2MN"
  lambda = @[0.0]
  gsteps = 4
  alwaysAccept:bool = 0

echo "rank ", myRank, "/", nRanks
threads: echo "thread ", threadNum, "/", numThreads

installStandardParams()
echoParams()
processHelpParam()

let
  lo = lat.newLayout
  vol = lo.physVol
  gc = actWilson(beta)

var r = lo.newRNGField(RngMilc6, seed)
var R:RngMilc6  # global RNG
R.seed(seed, 987654321)

var
  g = lo.newgauge
  p = lo.newgauge

if fileExists(gaugefile):
  tic("load")
  if 0 != g.loadGauge gaugefile:
    qexError "failed to load gauge file: ", gaugefile
  qexLog "loaded gauge from file: ", gaugefile," secs: ",getElapsedTime()
  toc("read")
  g.reunit
  toc("reunit")
else:
  #g.random r
  g.unit

g.echoPlaq

let gdt = toGvalue dt
var params = @[gdt]
let
  gg = toGvalue g
  gp = toGvalue p
  ga0 = gc.gaugeAction gg
  t0 = 0.5 * gp.norm2
  h0 = ga0 + t0
  tau = float(gsteps) * gdt
  (g1, p1, coeffs) = case gintalg
    of "2MN":
      int2MN(gc, gg, gp, gdt, gsteps, lambda)
    of "4MN3F1GP":
      int4MN3F1GP(gc, gg, gp, gdt, gsteps, lambda)
    of "4MN5F2GP":
      int4MN5F2GP(gc, gg, gp, gdt, gsteps, lambda)
    else:
      raise newException(ValueError, "unknown intalg: " & gintalg)
  ga1 = gc.gaugeAction g1
  t1 = 0.5 * p1.norm2
  h1 = ga1 + t1
  dH = h1 - h0
  acc = cond(dH<0.0, 1.0, exp(-dH))
  loss = -acc * (tau * tau)

params.add coeffs
var grads = newseq[Gvalue]()
for x in params:
  grads.add loss.grad x

var param = newseq[float]()
for x in params:
  param.add x.getfloat
var grad = param
var opt = newAdamW(param, lambda = weightDecay)

block:
  var ps = "param:"
  for i in 0..<params.len:
    ps &= " " & $param[i]
  echo ps

toc("prep")

for traj in 1..(trajsThermo+trajsTrain+trajsInfer):
  tic("traj")
  if traj <= trajsThermo:
    echo "Thermolization step: ",traj
  elif traj <= trajsThermo + trajsTrain:
    echo "Training step: ",traj-trajsThermo
  else:
    echo "Inference step: ",traj-(trajsThermo+trajsTrain)

  threads:
    p.randomTAH r
  gp.updated

  echo "Begin H: ",h0.eval,"  Sg: ",ga0.eval,"  T: ",t0.eval
  echo "End H: ",h1.eval,"  Sg: ",ga1.eval,"  T: ",t1.eval
  let accr = R.uniform
  if accr <= acc.eval.getfloat or alwaysAccept:
    echo "ACCEPT:  dH: ",dH.eval,"  exp(-dH): ",acc.eval,"  r: ",accr,(if alwaysAccept:" (ignored)" else:"")
  else:  # reject
    echo "REJECT:  dH: ",dH.eval,"  exp(-dH): ",acc.eval,"  r: ",accr
  qexGC "traj done"

  toc("forward end")

  if traj <= trajsThermo:
    echo "bloss: ",loss.eval
  elif traj <= trajsThermo + trajsTrain:
    # training
    tic()
    let t = traj - trajsThermo
    echo "tloss: ",loss.eval
    var gs = "grad:"
    for i in 0..<grads.len:
      grad[i] = grads[i].eval.getfloat
      gs &= " " & $grad[i]
    echo gs
    let lr = warmUpCosDecay(t, trajsTrainlrWarm, trajsTrain, lrmax, lrmin)
    echo "lr: ",lr
    opt.optimize(param, grad, t, lr)
    var ps = "param:"
    for i in 0..<params.len:
      ps &= " " & $param[i]
      params[i].update param[i]
    echo ps
    toc("training")
  else:
    echo "iloss: ",loss.eval

  if accr <= acc.eval.getfloat or alwaysAccept:
    gg.valCopy g1
    gg.getgauge.reunit
    gg.updated

  g.echoPlaq

  if savefreq > 0 and traj mod savefreq == 0:
    tic("save")
    let fn = savefile & &".{traj:05}.lime"
    if 0 != g.saveGauge(fn):
      qexError "Failed to save gauge to file: ",fn
    qexLog "saved gauge to file: ",fn," secs: ",getElapsedTime()
    toc("done")

  qexLog "traj ",traj," secs: ",getElapsedTime()
  toc("traj end")

toc()

processSaveParams()
writeParamFile()
qexFinalize()
