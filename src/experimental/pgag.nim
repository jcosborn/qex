import qex, gauge, physics/[qcdTypes,stagSolve]
import times, macros, algorithm
import hmc/metropolis
import observables/sources
import quda/qudaWrapper
import hmc/agradOps

qexinit()

let
  lat = intSeqParam("lat", @[8,8,8,8])
  beta = floatParam("beta", 5.6)
  tau0 = floatParam("tau", 0.04)
  lam0 = floatParam("lam", 0.19)
  sig0 = floatParam("sig", 0.3)
  rho0 = floatParam("rho", 0.2)
  theta0 = floatParam("theta", 0.2)
  ups0 = floatParam("ups", 0.2)
  upit = intParam("upit", 1)
  gsteps0 = intParam("gsteps", 32)
  nwarm = intParam("nwarm", 10)
  trajs = intParam("trajs", 10)
  seed0 = defaultComm.broadcast(int(1000*epochTime()))
  seed = uint64 intParam("seed", seed0)
  infn = stringParam("infn", "")
  outfn = stringParam("outfn", "")
var
  tau = tau0
  lam = lam0
  sig = sig0
  rho = rho0
  ups = ups0
  theta = theta0
  gsteps = gsteps0

macro echoparam(x: typed): untyped =
  let n = x.repr
  result = quote do:
    echo `n`, ": ", `x`

echoparam(beta)
echoparam(tau)
echoparam(lam)
echoparam(sig)
echoparam(rho)
echoparam(theta)
echoparam(ups)
echoparam(upit)
echoparam(gsteps)
echoparam(nwarm)
echoparam(trajs)
echoparam(seed)
echoparam(infn)
echoparam(outfn)

let
  gc = GaugeActionCoeffs(plaq: beta, adjplaq: 0)
  lo = lat.newLayout
  vol = lo.physVol

var r = lo.newRNGField(RngMilc6, seed)
var R:RngMilc6  # global RNG
R.seed(seed, 987654321)

var g = lo.newgauge
type
  Gauge = typeof(g)
  GaugeV = AgVar[Gauge]
template newGaugeV(c: AgTape, x: Gauge): auto = newGaugeFV(c, x)
case infn
of "":
  g.random r
  #g.unit
of "hot":
  g.random r
of "cold":
  g.unit
else:
  echo "Loading gauge field from file: ", infn
  let err = g.loadGauge infn

echo "plaq: ", 6.0*g.plaq
echo "gaugeAction2: ", g.gaugeAction2 gc
echo "actionA: ", gc.actionA g

var
  p = lo.newgauge
  f = lo.newgauge
  g0 = lo.newgauge

proc norm2*(x: seq): float =
  for i in 0..<x.len:
    result += x[i].norm2
    #echo result
  #result /= (x.len * x[0].l.physVol).float

proc norm2v*(x: seq): float =
  for i in 0..<x.len:
    result += x[i].norm2
    #echo result
  result /= (x.len * x[0].l.physVol).float

proc norm2subtract*(x: Field, y: float): float =
  var s: evalType(norm2(toDouble(x[0])))
  for i in x:
    s += x[i].toDouble.norm2 - y
  result = s.simdReduce
  x.l.threadRankSum(result)

proc norm2subtract*(x: seq, y: float): float =
  for i in 0..<x.len:
    result += x[i].norm2subtract(y)

proc forceX(c: GaugeActionCoeffs, g: auto, f: auto) =
  when defined(qudaDir):
    #var ff = f.newOneOf
    #gc.forceA(g, ff)
    gc.qudaGaugeForce(g, f)
    #echo "QEX:  ", ff[0].norm2
    #echo "QUDA: ", f[0].norm2
    #echo ff[0][0]
    #echo f[0][0]
    #echo (f[0]-ff[0]).norm2
  else:
    gc.forceA(g, f)

var gutime = 0.0
proc mdt(t: float) =
  tic()
  threads:
    for mu in 0..<g.len:
      for s in g[mu]:
        g[mu][s] := exp(t*p[mu][s])*g[mu][s]
  gutime += getElapsedTime()
  toc("mdt")
var gf2s = 0.0
var ngf2 = 0
var gftime = 0.0
proc mdv(t: float) =
  tic()
  gc.forceX(g, f)
  let gf2 = f.norm2v
  gf2s += gf2
  inc ngf2
  threads:
    for mu in 0..<f.len:
      p[mu] -= t*f[mu]
  gftime += getElapsedTime()
  toc("mdv")

proc mdvx(t: float) =
  tic()
  gc.forceX(g, f)
  let gf2 = f.norm2v
  gf2s += gf2
  inc ngf2
  threads:
    for mu in 0..<f.len:
      p[mu] -= t*f[mu]
  gftime += getElapsedTime()
  toc("mdv")

var
  gauges = newSeq[Gauge](0)
  moms = newSeq[Gauge](0)
  gaugevs = newSeq[GaugeV](0)
  momvs = newSeq[GaugeV](0)
  params = newSeq[FloatV](0)
  ptemps = newSeq[FloatV](0)
var tape = newAgTape()
proc pushParam(x: float): FloatV =
  params.add tape.newFloatV(x)
  result = params[^1]

gauges.add g0
gaugevs.add tape.newGaugeV(g0)
var p0 = lo.newGauge
moms.add p0
momvs.add tape.newGaugeV(p0)
#params.add tape.newFloatV(tau)
#let vtau = params[0]
let vtau = pushParam(tau)

proc pushTemp =
  ptemps.add tape.newFloatV
proc pushGauge =
  gauges.add lo.newGauge
  gaugevs.add tape.newGaugeV(gauges[^1])
proc pushMom =
  moms.add lo.newGauge
  momvs.add tape.newGaugeV(moms[^1])

proc addT(veps: FloatV) =
  pushGauge()
  exp(gaugevs[^1], veps, momvs[^1])
  pushGauge()
  mul(gaugevs[^1], gaugevs[^2], gaugevs[^3])

proc addGf(g: GaugeV) =
  pushMom()
  gc.gderiv(momvs[^1], g)
  pushMom()
  mulna(momvs[^1], g, momvs[^2])
  pushMom()
  projtah(momvs[^1], momvs[^2])

proc addGx(veps: FloatV, p,g: GaugeV) =
  addGf(g)
  pushMom()
  xpay(momvs[^1], p, veps, momvs[^2])

proc addG(veps: FloatV) =
  addGx(veps, momvs[^1], gaugevs[^1])

proc addGF(va, vb: FloatV) =
  let p = momvs[^1]
  addGf(gaugevs[^1])
  pushMom()
  exp(momvs[^1], vb, momvs[^2])
  pushMom()
  mul(momvs[^1], momvs[^2], gaugevs[^1])
  addGx(va, p, momvs[^1])

proc setupMDx =
  pushTemp()
  let vtau2 = ptemps[^1]
  mul(vtau2, 0.5, vtau)
  addT(vtau2)
  addG(vtau)
  addT(vtau2)

proc setupMDa =
  pushTemp()
  let vtau2 = ptemps[^1]
  mul(vtau2, 0.5, vtau)
  pushTemp()
  let vtau4 = ptemps[^1]
  mul(vtau4, 0.25, vtau)
  addT(vtau4)
  addG(vtau2)
  addT(vtau2)
  addG(vtau2)
  addT(vtau4)

proc `-`(x: FloatV, y: FloatV): FloatV =
  pushTemp()
  result = ptemps[^1]
  sub(result, x, y)
proc `*`(x: float, y: FloatV): FloatV =
  pushTemp()
  result = ptemps[^1]
  mul(result, x, y)
proc `*`(x: FloatV, y: FloatV): FloatV =
  pushTemp()
  result = ptemps[^1]
  mul(result, x, y)

proc setupMD2 =
  let lam = pushParam(lam0)
  let lamt = lam * vtau
  let sig = vtau - 2 * lamt
  let tau2 = 0.5 * vtau
  addT(lamt)
  addG(tau2)
  addT(sig)
  addG(tau2)
  addT(lamt)

proc setupMD2p =
  let lam = pushParam(lam0)
  let lamt = lam * vtau
  let sig = vtau - 2 * lamt
  let tau2 = 0.5 * vtau
  addG(lamt)
  addT(tau2)
  addG(sig)
  addT(tau2)
  addG(lamt)

proc setupMD2g =
  let lam = pushParam(lam0)
  let lamt = lam * vtau
  let sig = vtau - 2 * lamt
  let tau2 = 0.5 * vtau
  let rho = pushParam(rho0)
  let rhot = rho * vtau * vtau * vtau
  addG(lamt)
  addT(tau2)
  addGF(sig, rhot)
  addT(tau2)
  addG(lamt)

proc setupMD3 =
  let lam = pushParam(lam0)
  let lamt = lam * vtau
  let sig = pushParam(sig0)
  let sigt = sig * vtau
  let alp = 0.5*vtau - lamt
  let bet = vtau - 2*sigt
  addT(lamt)
  addG(sigt)
  addT(alp)
  addG(bet)
  addT(alp)
  addG(sigt)
  addT(lamt)

proc setupMD3g =
  let lam = pushParam(lam0)
  let lamt = lam * vtau
  let sig = pushParam(sig0)
  let sigt = sig * vtau
  let alp = 0.5*vtau - lamt
  let bet = vtau - 2*sigt
  let rho = pushParam(rho0)
  let rhot = rho * vtau * vtau * vtau
  addT(lamt)
  addG(sigt)
  addT(alp)
  addGF(bet, rhot)
  addT(alp)
  addG(sigt)
  addT(lamt)

proc setupMD3p =
  let lam = pushParam(lam0)
  let lamt = lam * vtau
  let sig = pushParam(sig0)
  let sigt = sig * vtau
  let alp = 0.5*vtau - lamt
  let bet = vtau - 2*sigt
  addG(lamt)
  addT(sigt)
  addG(alp)
  addT(bet)
  addG(alp)
  addT(sigt)
  addG(lamt)

# 2x OMF2: 0.1 0.25 0.3
proc setupMD4 =
  let lam = pushParam(lam0)
  let lamt = lam * vtau
  let sig = pushParam(sig0)
  let sigt = sig * vtau
  let rho = pushParam(rho0)
  let rhot = rho * vtau
  let alp = 0.5*vtau - sigt
  let bet = vtau - 2*lamt - 2*rhot
  addT(lamt)
  addG(sigt)
  addT(rhot)
  addG(alp)
  addT(bet)
  addG(alp)
  addT(rhot)
  addG(sigt)
  addT(lamt)

proc setupMD4g =
  let lam = pushParam(lam0)
  let lamt = lam * vtau
  let sig = pushParam(sig0)
  let sigt = sig * vtau
  let rho = pushParam(rho0)
  let rhot = rho * vtau
  let alp = 0.5*vtau - sigt
  let bet = vtau - 2*lamt - 2*rhot
  let theta = pushParam(theta0)
  let thetat = theta * vtau * vtau * vtau
  addG(lamt)
  addT(sigt)
  addG(rhot)
  addT(alp)
  addGF(bet, thetat)
  addT(alp)
  addG(rhot)
  addT(sigt)
  addG(lamt)

proc setupMD5 =
  let lam = pushParam(lam0)
  let lamt = lam * vtau
  let sig = pushParam(sig0)
  let sigt = sig * vtau
  let rho = pushParam(rho0)
  let rhot = rho * vtau
  let theta = pushParam(theta0)
  let thetat = theta * vtau
  let alp = 0.5*vtau - lamt - rhot
  let bet = vtau - 2*sigt - 2*thetat
  addT(lamt)
  addG(sigt)
  addT(rhot)
  addG(thetat)
  addT(alp)
  addG(bet)
  addT(alp)
  addG(thetat)
  addT(rhot)
  addG(sigt)
  addT(lamt)

# 6x LF: 0.08 0.17 0.17 0.17 0.17
# 3x OMF2: 0.06 0.17 0.2 0.17 0.13
proc setupMD6 =
  let lam = pushParam(lam0)
  let lamt = lam * vtau
  let sig = pushParam(sig0)
  let sigt = sig * vtau
  let rho = pushParam(rho0)
  let rhot = rho * vtau
  let theta = pushParam(theta0)
  let thetat = theta * vtau
  let ups = pushParam(ups0)
  let upst = ups * vtau
  let alp = 0.5*vtau - sigt - thetat
  let bet = vtau - 2*lamt - 2*rhot - 2*upst
  addT(lamt)
  addG(sigt)
  addT(rhot)
  addG(thetat)
  addT(upst)
  addG(alp)
  addT(bet)
  addG(alp)
  addT(upst)
  addG(thetat)
  addT(rhot)
  addG(sigt)
  addT(lamt)

proc update =
  tape.run
  for mu in 0..<g.len:
    g[mu] := gauges[^1][mu]
    p[mu] := moms[^1][mu]

var
  p2xv = tape.newFloatV
  p2v = tape.newFloatV
  gav = tape.newFloatV
  hv = tape.newFloatV
proc addAction(p: GaugeV, g: GaugeV) =
  norm2subtract(p2xv, p, 8.0)
  mul(p2v, 0.5, p2xv)
  gaction(gc, gav, g)
  add(hv, p2v, gav)

proc setupAction =
  tape.addTrack
  addAction(momvs[0], gaugevs[0])
  tape.addTrack
  addAction(momvs[^1], gaugevs[^1])
  tape.setTrack 0
  #echo tape

type
  Met = ref object of MetropolisRoot
    state: int

proc init(m: var Met) =
  m.new
  var r = MetropolisRoot m
  init(r)
  m.verbosity = 1

proc start*(m: var Met) =
  tic()
  m.state = 0
  threads:
    for i in 0..<g.len:
      g0[i] := g[i]
    p.randomTAH r
    for i in 0..<p.len:
      p0[i] := p[i]
    threadBarrier()
    #psi.gaussian r
    #threadBarrier()
    #g.rephase
    #threadBarrier()
    #stag.D(phi, psi, mass)
    #threadBarrier()
    #phi.odd := 0
    #g.rephase
  toc("init p, phi")

proc getH*(m: Met): float =
  tic()
  var p2 = 0.0
  var f2 = 0.0
  #threads:
  #  g.rephase
  #stag.solve(psi, phi, mass, spa)
  #toc("fa solve")
  #echo "psi e: ", psi.even.norm2
  #echo "psi o: ", psi.odd.norm2
  threads:
    #g.rephase
    var p2t = 0.0
    for i in 0..<p.len:
      #p2t += p[i].norm2
      p2t += p[i].norm2subtract(8.0)
    threadMaster: p2 = p2t
    #var psi2 = psi.norm2
    #var psi2 = psi.norm2subtract(3.0)
    #threadMaster: f2 = psi2
  let
    ga0 = gc.actionA g
    #fa0 = 0.5*f2 - (1.5*vol).float
    #t0 = 0.5*p2 - (16*vol).float
    fa0 = 0.5*f2
    t0 = 0.5*p2
    h0 = ga0 + fa0 + t0
  result = h0
  toc("end getH")
  if m.state == 0:
    inc m.state
    echo "Begin H: ",h0,"  Sg: ",ga0,"  Sf: ",fa0,"  T: ",t0
    tape.setTrack 1
    tape.run
    echo &"      H: {hv.obj}  Sg: {gav.obj}  T: {p2v.obj}"
    tape.setTrack 0
  else:
    echo "End H: ",h0,"  Sg: ",ga0,"  Sf: ",fa0,"  T: ",t0
    tape.setTrack 2
    tape.run
    echo &"      H: {hv.obj}  Sg: {gav.obj}  T: {p2v.obj}"
    tape.setTrack 0

proc getCost(m: Met): seq[float] =
  result.newSeq(0)
  var cost = vtau.obj * m.avgPAccept
  result.add cost
  for i in 0..<params.len:
    var costg = (if i==0: 1.0 else: 0.0)
    if m.hNew > m.hOld:
      costg = m.pAccept * (costg - vtau.obj * params[i].grad)
    result.add costg

proc checkGrad(m: Met) =
  #let tg = vtau.grad
  #let f0 = vtau.obj * exp(m.hOld-m.hNew)
  let eps = 1e-7
  var gs = newSeq[float](0)
  for i in 0..<params.len: gs.add params[i].grad
  for i in 0..<params.len:
    let t = params[i].obj
    params[i].obj += eps
    update()
    let h = m.getH
    echo "params[",i,"] grad: ", gs[i]
    echo "params[",i,"] diff: ", (h-m.hNew)/eps
    params[i].obj = t
  update()

# w = sum_i w_i
# sum = sum_i w_i x_i
# sum2 = sum_i w_i x_i^2
# X = sum/w
# var = sum_i w_i (x_i - X)^2 = sum_i w_i x_i^2 - w X^2
type DecayStat = object
  n: int
  fac: float
  w: float
  sum: float
  sum2: float
proc newDecayStat(f: float): DecayStat =
  result.fac = f
proc push(d: var DecayStat, x: float) =
  inc d.n
  let f = d.fac
  let o = 1.0 - f
  d.w = o*d.w + f
  d.sum = o*d.sum + f*x
  d.sum2 = o*d.sum2 + f*x*x
proc mean(d: DecayStat): float =
  result = d.sum / d.w
proc variance(d: DecayStat): float =
  let m = d.sum / d.w
  let m2 = d.sum2 / d.w
  result = sqrt(m2 - m*m)
#var decay = 0.1
var decay = 0.2
var cgstat = newSeq[RunningStat](0)
#var cgstat = newSeq[DecayStat](0)

proc getgrad(m: Met) =
  hv.grad = 1.0
  tape.setTrack 2
  tape.grad
  tape.setTrack 0
  #gc.gaugeDeriv2(gaugevs[^1].obj, gaugevs[^1].grad)
  #for mu in 0..<g.len:
  #  momvs[^1].grad[mu] := momvs[^1].obj[mu]
  tape.grad
  #let tg = vtau.grad
  let cst = getCost(m)
  echo "cost:  ", cst[0]
  if cgstat.len < params.len: cgstat.setLen(params.len)
  #if cgstat.len < params.len:
  #  for i in 0..<params.len: cgstat.add newDecayStat(decay)
  for i in 0..<cgstat.len:
    cgstat[i].push cst[i+1]
    let m = cgstat[i].mean
    let v = cgstat[i].variance
    #echo "costg: ", cgstat[i].mean, " ", cgstat[i].standardDeviationS
    #echo "costg: ", cgstat[i].mean, " ", cst[i+1]
    echo &"costg: {m/v:8.6f} {m:10.6f} {cst[i+1]}"

proc updateParams(rate: float) =
  let n = cgstat.len
  var s = 0.0
  var p = newSeq[float](n)
  var g = newSeq[float](n)
  for i in 0..<n:
    let m = cgstat[i].mean
    let v = cgstat[i].variance
    let t = m / v
    g[i] = t
    s += t*t
    p[i] = params[i].obj
  echo "Params: ", p
  s = rate*sqrt(n/s)
  for i in 0..<n:
    #let m = cgstat[i].mean
    let d = s*g[i]
    params[i].obj += d
    p[i] = params[i].obj
  echo "Params: ", p
  for i in 0..<n:
    clear cgstat[i]

proc finish*(m: var Met) =
  discard
  #for mu in 0..<g.len:
  #  g[mu] := g0[mu]
  #  p[mu] := p0[mu]
  #discard m.getH
  #mdt(0.5*tau)
  #echo "mdt: ", g.norm2, "  ", vg2.obj.norm2
  #mdv(tau)
  #echo "mdv: ", p.norm2, "  ", vp2.obj.norm2
  #mdt(0.5*tau)
  #echo "mdt: ", g.norm2, "  ", vg6.obj.norm2
  # pacc = max(1, exp(hOld-hNew))
  #if m.hNew <= m.hold:
  #gc.gaugeDeriv2(gaugevs[^1].obj, gaugevs[^1].grad)
  #for mu in 0..<g.len:
  #  momvs[^1].grad[mu] := momvs[^1].obj[mu]
  #tape.grad
  ##let tg = vtau.grad
  #let cst = getCost(m)
  #echo "cost:  ", cst[0]
  ##if cgstat.len < params.len: cgstat.setLen(params.len)
  #if cgstat.len < params.len:
  #  for i in 0..<params.len: cgstat.add newDecayStat(decay)
  #for i in 0..<cgstat.len:
  #  cgstat[i].push cst[i+1]
  #  #echo "costg: ", cgstat[i].mean, " ", cgstat[i].standardDeviationS
  #  echo "costg: ", cgstat[i].mean, " ", cst[i+1]
  #checkGrad(m)

proc disp(m: Met) =
  let p = g.plaq
  echo "plaq: ", 6.0*p
  echo "tplaq: ", p.sum
  echo "gauge force: ", sqrt(gf2s/ngf2)
  gf2s = 0.0
  ngf2 = 0
  echo "rmsDeltaH: ", sqrt(m.avgDeltaH2)
  echo "avgPAccept: ", m.avgPAccept

proc accept*(m: var Met) =
  disp(m)

proc reject*(m: var Met) =
  threads:
    for i in 0..<g.len:
      g[i] := g0[i]
  disp(m)

var alwaysAccept = false
proc globalRand*(m: Met): float =
  if not alwaysAccept:
    result = R.uniform

proc generate*(m: var Met) =
  update()

var m: Met
m.init

var src = lo.ColorVector
var prop = lo.ColorVector
let nt = lat[3]
let nt2 = (nt+2) div 2
#var picorr = newSeq[float](lat[3])
var picorr = newSeq[RunningStat](nt)
var picorrf = newSeq[RunningStat](nt2)
var sp = newSolverParams()
sp.r2req = 1e-20
block:
  var v: evalType(src{0})
  v := 0
  v[0] = 1
  threads:
    src.wallSource(0, v)
  #echo src.norm2slice(3)

proc measure =
  #threads:
  #  g.rephase
  #stag.solve(prop, src, mass, sp)
  #threads:
  #  g.rephase
  #picorr += prop.norm2slice(3)
  #echo picorr
  #let pic = prop.norm2slice(3)
  #for i in 0..<nt:
  #  picorr[i].push pic[i]
  #  let k = min(i, nt-i)
  #  picorrf[k].push pic[i]
  #for i in 0..<nt:
  #  echo i, " ", picorr[i].mean, " ", picorr[i].standardDeviationS
  #for i in 0..<nt2:
  #  let k = nt2 - 1 - i
  #  echo i, " ", picorr[k].mean, " ", picorr[k].standardDeviationS
  discard

#setupMD2()
#setupMD2p()
setupMD2g()
#setupMD3()
#setupMD3g()
#setupMD4()
#setupMD4g()
#setupMD6()

setupAction()

if nwarm > 0:
  echo "Starting warmups"
  #gsteps = 60
  #fsteps = 40
  #setupMDx()
  alwaysAccept = true
  for n in 1..nwarm:
    m.update
  m.clearStats

echo "Starting HMC"
#gsteps = gsteps0
#fsteps = fsteps0
#setupMD5()
alwaysAccept = false
gutime = 0.0
gftime = 0.0
#fftime = 0.0
block:
  tic()
  for n in 1..trajs:
    echo "Starting trajectory: ", n
    tic()
    m.update
    getGrad(m)
    if upit > 0:
      if n mod upit == 0:
        updateParams(0.001)
    let tup = getElapsedTime()
    measure()
    let ttot = getElapsedTime()
    echo "End trajectory update: ", tup, "  measure: ", ttot-tup, "  total: ", ttot
  let et = getElapsedTime()
  toc()
  echo "HMC time: ", et
  #let at = gutime + gftime + fftime
  #echo &"gu: {gutime}  gf: {gftime}  ff: {fftime}  ot: {et-at}  tt: {et}"

if outfn != "":
  echo "Saving gauge field to file: ", outfn
  let err = g.saveGauge outfn

#echoTimers()
qexfinalize()
