import base/globals
setDefaultNc(1)
setVLENmax(4)

import qex
import gauge, physics/qcdTypes
import algorithms/integrator, maths/special, utils/resample
import os, strutils, times
import numericalnim
from numericalnim/ode import IntegratorProc
type IntegratorProc*[T] = ode.IntegratorProc[T]

proc mean(xs: Ensemble[seq[float]]): float =
  var m = 0.0
  for i in 0..<xs.len:
    let x = xs[i]
    m += (x - m)/float(i+1)
  m

proc naiveIntAutocorr(xs: Ensemble[seq[float]], maxlen: int): float =
  ## Very rough integrated autocorrelation: sum of correlation out to some cut.
  let N = xs.len
  let m = mean(xs)
  # c(0)
  var c0 = 0.0
  for i in 0..<xs.len:
    let d = xs[i]-m
    c0 += d*d
  c0 /= float(N)
  if c0 <= 1e-14: return 1.0

  var sumRho = 1.0
  for lag in 1..(min(maxlen, N-1)):
    var cLag = 0.0
    for i in 0..<(N-lag):
      cLag += (xs[i] - m)*(xs[i+lag] - m)
    cLag /= float(N-lag)
    let rho = cLag/c0
    if rho < 0.0:
      break
    sumRho += 2.0*rho
  sumRho

proc tunnelingRate(xs: Ensemble[seq[float]]): float =
  var r = 0.0
  for i in 1..<xs.len:
    let d = xs[i]-xs[i-1]
    r += (d*d-r)/float(i)
  sqrt(r)

proc topo2DU1(g:auto):float =
  tic()
  const nc = g[0][0].nrows
  let
    lo = g[0].l
    nd = lo.nDim
    t = newTransporters(g, g[0], 1)
  var p = 0.0
  toc("topo2DU1 setup")
  threads:
    tic()
    var tp:type(atan2(g[0][0][0,0].im, g[0][0][0,0].re))
    for mu in 1..<nd:
      for nu in 0..<mu:
        let tpl = (t[mu]^*g[nu]) * (t[nu]^*g[mu]).adj
        for i in tpl:
          tp += atan2(tpl[i][0,0].im, tpl[i][0,0].re)
    var v = tp.simdSum
    v.threadRankSum
    threadSingle: p += v
    toc("topo2DU1 work")
  toc("topo2DU1 threads")
  p/TAU

proc infVolPlaq(beta: float): float = besselI1(beta)/besselI0(beta)
proc infVolChiQ(beta: float, tolerance = 1.0e-15): float =
  let i0beta = besselI0(beta)
  var
    term = PI * PI / 12.0
    sumVal = term
    n = 1
  while abs(term) > tolerance*abs(sumVal):
    let inbeta = besselIn(n, beta)
    # term = ((-1)^n / n^2) * (I_n(beta) / I_0(beta))
    term = inbeta / (i0beta * float(n*n))
    if n mod 2 == 0:
      sumVal += term
    else:
      sumVal -= term
    inc n
  sumVal / (PI * PI)

proc pnorm2(p:auto):float =
  tic()
  var p2 = 0.0
  threads:
    var p2t = 0.0
    for i in 0..<p.len:
      p2t += p[i].norm2
    threadMaster: p2 = p2t
  toc("pnorm2")
  p2

proc hamiltonian(gc:GaugeActionCoeffs, g:auto, p:auto):auto =
  tic()
  let
    vol = float(g[0].l.physVol)
    ga = gc.plaq*vol + g.gaugeAction2 gc
    p2 = p.pnorm2
    t = 0.5*p2 - vol
    h = ga + t
  toc("hamiltonian")
  (ga, t, h)

type MDAlgo = enum hamilton, nosehoover, ghmc, gv, test1
converter toMDAlgo(s:string):MDAlgo = parseEnum[MDAlgo](s)
let gAlgs = [ghmc, gv]

qexinit(verb=2)

tic()

let seed0 = defaultComm.broadcast(uint(1000*epochTime()))

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
      @[64,64]
  beta = 5.0
  bg = beta
  tau = 2.0
  ntraj = 128
  ntrajThermo = 64
  ntrajThermoAcc = ntrajThermo div 8
  jkBlockSize = max(1, ntraj div 64)
  seed = seed0
  ## hamilton | nosehoover [H + xi^2/(2 gamma); dU = dt P, dP = - dt S'(U) - xi P, dxi = gamma dt (p^2 - dof)]
  mdalgo:MDAlgo = "hamilton"
  ## extra params in other mdalgo
  gamma = 1.0
  ## 2MN,0.21 | 4MN3F1GP,0.27 | 4MN5F2GP
  gintalg:integrator.IntegratorProc = "2MN,0.21"
  omfLambda = 0.21
  gsteps = 20
  intsteps = 100
  relTol = 1e-12
  reduceT = true
  alwaysAccept:bool = 0
  revCheckFreq = ntraj
  verboseGCStats:bool = 0
  verboseTimer:bool = 0

installStandardParams()
echoParams()
#qexLog "rank ", myRank, "/", nRanks
#threads: qexLog "thread ", threadNum, "/", numThreads
processHelpParam()

VerboseGCStats = verboseGCStats
VerboseTimer = verboseTimer

let
  lo = lat.newLayout
  nd = lo.nDim
  gc = GaugeActionCoeffs(plaq:beta)
  vol = lo.physVol
  dof = float(2*vol)
  useG = mdAlgo in gAlgs

var
  g = lo.newGauge
  r = lo.newRNGField(RngMilc6, seed)
  R:RngMilc6  # global RNG
R.seed(seed, 987654321)

g.random r

qexLog "Initial plaq: ",g.plaq3

####  modify to add new HMC variants  ###
var
  p = lo.newgauge
  f = lo.newgauge
  gg = lo.newgauge  # FG backup gauge
  g0 = lo.newgauge
  p0 = lo.newgauge
  xi = 0.0
  xi1 = xi
  lnJ = 0.0

proc revBegin =
  case mdalgo
  of nosehoover:
    xi1 = xi
    xi *= -1
  else:
    discard

proc revEnd =
  case mdalgo
  of nosehoover:
    xi = xi1
  else:
    discard

proc extraInit =
  lnJ = 0.0
  case mdalgo
  of nosehoover:
    xi = R.gaussian * sqrt(gamma)
  else:
    discard

proc extraH: float =
  case mdalgo
  of nosehoover:
    result = xi * xi / (2.0 * gamma)
  else:
    discard

var
  gf = lo.newgauge
  gb: array[2,typeof(gf)]
  gbb = lo.newgauge
for i in 0..<nd: gb[i] = lo.newgauge()

proc initgV =
  let tf = newShifters(g[0], 1)
  let tb = newShifters(g[0], -1)
  threads:
    for i in 0..<nd:
      gf[i] := tf[1-i] ^* g[i]
      for j in 0..<nd:
        gb[i][j] := tb[i] ^* g[j]
        threadBarrier()
      gbb[i] := tb[i] ^* gf[i]

#[
proc gfunV(x,mu: auto): auto =
  let nu = 1 - mu
  let su = g[nu][x] * gf[mu][x] * gf[nu][x].adj
  let sd = gb[nu][nu][x].adj * gb[nu][mu][x] * gbb[nu][x]
  let s = su + sd
  let p = g[mu][x] * s.adj
  let a = beta * (2 - p.re)
  result = exp(a)
]#

proc gfunderivV(x,mu: auto): auto =
  let nu = 1 - mu
  let su = g[nu][x] * gf[mu][x] * gf[nu][x].adj
  let sd = gb[nu][nu][x].adj * gb[nu][mu][x] * gbb[nu][x]
  let s = su + sd
  let p = g[mu][x] * s.adj
  let pt = trace(p)
  let a = bg * (1 - pt.re)
  let f = exp(a)
  let d = (0.5*bg*f)*(p - p.adj)
  result = (f,d)

proc gpotV(x,mu: auto): auto =
  let nu = 1 - mu
  let su = g[nu][x] * gf[mu][x] * gf[nu][x].adj
  let sd = gb[nu][nu][x].adj * gb[nu][mu][x] * gbb[nu][x]
  let s = su + sd
  let p = g[mu][x] * s.adj
  let pt = trace(p)
  let a = exp(bg)
  let b = -bg * pt
  result = (a,b)

proc initg =
  case mdalgo
  of gv:
    initgV()
  else:
    discard

#[
proc gfun(x,mu: auto): auto =
  case mdalgo
  of ghmc:
    result = 1
  of gv:
    result = gfunV(x, mu)
]#

#[
proc gderiv(x,mu: auto): auto =
  case mdalgo
  of ghmc:
    result = 0
  of gv:
    result = gderivV(x, mu)
]#

proc gfunderiv(x,mu: auto): auto =
  var r: typeof(gfunderivV(x, mu))
  case mdalgo
  of ghmc:
    r[0] := 1
    r[1] := 0
  of gv:
    r = gfunderivV(x, mu)
  else:
    discard
  result = r

proc gpot(x,mu: auto): auto =
  var r: typeof(gpotV(x, mu))
  case mdalgo
  of ghmc:
    r[0] := 1
    r[1] := 0
  of gv:
    r = gpotV(x, mu)
  else:
    discard
  result = r

proc gupdate0(x,mu: int, t: float): auto =
  let dt = t / intsteps
  var lnJt: evalType(g[mu][x][0,0].re)
  for s in 0..<intsteps:
    let gd = gfunderiv(x,mu)
    let tg = dt*gd[0]
    let te = exp(tg*p[mu][x])
    let etpg = te*g[mu][x]
    let j = 1 + (dt*gd[1][0,0].im)*p[mu][x][0,0].im
    lnJt -= ln(j)
    g[mu][x] := etpg
  result = lnJt

proc gupdate1(x,mu: int, t: float): auto =
  var lnJt: evalType(g[mu][x][0,0].re)
  let dt = t / intsteps
  let gp = gpot(x,mu)
  let b = gp[0] * p[mu][x][0,0].im
  #var s: typeof(lnJt)
  var ae = gp[1]
  for i in 1..intsteps:
    let acc = ae.re
    let asc = ae.im
    let s1 = exp(acc)
    #let s0 = s1*dt
    let s2 = -0.5*b*asc*s1*s1
    let s0 = s1*dt + s2*dt*dt
    #let s3 = (1.0/6.0)*(2*asc*asc-acc)*b*b*s1*s1*s1
    #let s0 = s1*dt + s2*dt*dt + s3*dt*dt*dt
    #let s4 = (1.0/24.0)*(-6*asc*asc*asc+7*acc*asc+asc)*b*b*b*s1*s1*s1*s1
    #let s0 = s1*dt + s2*dt*dt + s3*dt*dt*dt + s4*dt*dt*dt*dt
    let f = exp(s0*gp[0]*p[mu][x])
    let u = f * g[mu][x]
    #let g0 = gp[0] * s1
    #let g1 = gp[0] * exp(re(ae*f[0,0]))
    #s += s0
    ae *= f[0,0]
    g[mu][x] := u
    #lnJt -= ln(g1/g0)
    lnJt -= re(ae) - acc
    #let j = 1 + re(ae)/acc
    #lnJt -= ln(j)
  result = lnJt

proc gupdate2(x,mu: int, t: float): auto =
  let gp = gpot(x,mu)
  let b = abs(gp[0] * p[mu][x][0,0].im)
  let s1 = 2*PI/b  # period in s
  let ds = s1 / intsteps
  let df = exp(ds*gp[0]*p[mu][x][0,0])
  var ae = gp[1]
  var dt = 0.5*exp(-ae.re)
  for i in 1..<intsteps:
    ae *= df
    dt += exp(-ae.re)
  ae *= df
  dt += 0.5*exp(-ae.re)
  dt /= intsteps  # period in t
  let nf = t / dt
  let r = dt*(nf - trunc(nf)) # the remainder of t modulo period dt
  ae = gp[1]
  let fac = 0.5/intsteps
  var s,r0,r1: typeof(b)
  var eae = fac * exp(-ae.re)
  for i in 1..intsteps:
    r0 = r1
    r1 += eae
    ae *= df
    eae = fac * exp(-ae.re)
    r1 += eae
    let dr = max(min(r,r1),r0)-r0
    s += dr*ds/(r1-r0)
  let f = exp(s*gp[0]*p[mu][x])
  let u = f * g[mu][x]
  g[mu][x] := u
  ae = gp[1]
  ae = gp[1]*f[0,0]
  result = gp[1].re - ae.re

const rkorder = 4
proc gupdaterk(x,mu: int, t: float): auto =
  let dt = t / intsteps
  let gp = gpot(x,mu)
  let lam = gp[0]*p[mu][x][0,0]
  var s: typeof(lam.re)
  for i in 1..intsteps:
    let ae = gp[1] * exp(s*lam)
    let k1 = dt * exp(ae.re)
    case rkorder
    of 1:
      let ds = k1  # first order
      s += ds
    of 2:
      let ae2 = ae * exp(k1*lam)
      let k2 = dt * exp(ae2.re)
      let ds = 0.5*(k1+k2)  # second order
      s += ds
    of 3:
      let ae2 = ae * exp(0.5*k1*lam)
      let k2 = dt * exp(ae2.re)
      let ae3 = ae * exp((2*k2-k1)*lam)
      let k3 = dt * exp(ae3.re)
      let ds = (1.0/6.0)*(k1+4*k2+k3)  # third order
      s += ds
    of 4:
      let ae2 = ae * exp(0.5*k1*lam)
      let k2 = dt * exp(ae2.re)
      let ae3 = ae * exp(0.5*k2*lam)
      let k3 = dt * exp(ae3.re)
      let ae4 = ae * exp(k3*lam)
      let k4 = dt * exp(ae4.re)
      let ds = (1.0/6.0)*(k1+2*k2+2*k3+k4)  # fourth order
      s += ds
    else: discard
  let f = exp(s*gp[0]*p[mu][x])
  let u = f * g[mu][x]
  g[mu][x] := u
  let ae = gp[1] * exp(s*lam)
  result = gp[1].re - ae.re

proc rkint(t: float, f: auto): float =
  let dt = t / intsteps
  var s = 0.0
  for i in 1..intsteps:
    let k1 = dt * f(s)
    case rkorder
    of 1:
      let ds = k1  # first order
      s += ds
    of 2:
      let k2 = dt * f(s+k1)
      let ds = 0.5*(k1+k2)  # second order
      s += ds
    of 3:
      let k2 = dt * f(s+0.5*k1)
      let k3 = dt * f(s+2*k2-k1)
      let ds = (1.0/6.0)*(k1+4*k2+k3)  # third order
      s += ds
    of 4:
      let k2 = dt * f(s+0.5*k1)
      let k3 = dt * f(s+0.5*k2)
      let k4 = dt * f(s+k3)
      let ds = (1.0/6.0)*(k1+2*k2+2*k3+k4)  # fourth order
      s += ds
    else: discard
  result = s

mixin ode.IntegratorProc
proc gupA(lam: auto, a: auto, t: float): float =
  # s ~ t / I0(|a|)
  let i0a = besselI0(sqrt(a.norm2))
  let smag = t / i0a
  var t = t
  if reduceT and lam.im != 0.0:
    let tp = i0a * 2 * PI / abs(lam.im)
    t -= tp * trunc(t/tp)
  #proc dsdt(y: float): float =
  proc dsdt(t: float, y: float, ctx: NumContext[float, float]): float =
    let ae = a * exp(y*lam)
    result = exp(ae.re)
  #result = rkint(t, dsdt)
  let s0 = 0.0
  let tspan = [t]
  #let odeOptions = newODEoptions(dtMin=0.0,absTol=1e-6)
  #let odeOptions = newODEoptions(dtMin=1e-14,dtMax=1e-3,absTol=1e-18,relTol=1e-20)
  #let odeOptions = newODEoptions(dtMin=1e-16,dtMax=1e-3,absTol=1e-12*smag,relTol=0.0)
  #let odeOptions = newODEoptions(dtMin=1e-16,dtMax=1e-3,absTol=0.0,relTol=1e-16/smag)
  let abst = 0.0
  let relt = relTol / smag
  let odeOptions = newODEoptions(dtMin=1e-16,dtMax=1e-3,absTol=abst,relTol=relt)
  #let intg = "rk21"
  let intg = "dopri54"
  #let intg = "tsit54"
  #let intg = "vern65"
  let (t1, y1) = solveOde(dsdt, s0, tspan, odeOptions, integrator=intg)
  result = y1[^1]

proc gupdateA(x,mu: int, t: float): auto =
  var s: evalType(g[mu][x][0,0].re)
  const vl = simdLength(s)
  let gp = gpot(x,mu)
  let lam = gp[0]*p[mu][x][0,0]
  for i in 0..<vl:
    let li = eval(lam[asSimd(i)])
    let ai = eval(gp[1][asSimd(i)])
    let si = gupA(li, ai, t)
    s[i] = si
  let f = exp(s*gp[0]*p[mu][x])
  let u = f * g[mu][x]
  g[mu][x] := u
  let ae = gp[1] * exp(s*lam)
  result = gp[1].re - ae.re

proc mdtfb(t:float, fb:int) =
  if useG:
    let eosub = [lo.getSubset("even"), lo.getSubset("odd")]
    let dt = t
    for mux in [0,1]:
      let mu = if fb==0: mux else: 1-mux
      for eox in [0,1]:
        let eo = if fb==0: eox else: 1-eox
        initg()
        threads:
          var lnJt: evalType(g[0][0][0,0].re)
          for x in eosub[eo]:
            lnJt += gupdateA(x,mu,dt)
          var lnJs = simdSum(lnJt)
          threadRankSum(lnJs)
          threadSingle: lnJ += lnJs
  else:
    threads:
      for i in 0..<g.len:
        for e in g[i]:
          let etpg = exp(t*p[i][e])*g[i][e]
          g[i][e] := etpg
    if mdalgo == nosehoover:
      xi += t * gamma * (p.pnorm2 - dof)

proc mdtf(t:float) =
  #let t = 0.25*t
  mdtfb(t, 0)
  #mdtfb(t, 1)
  #mdtfb(t, 0)
  #mdtfb(t, 1)

proc mdtb(t:float) =
  #let t = 0.25*t
  #mdtfb(t, 0)
  mdtfb(t, 1)
  #mdtfb(t, 0)
  #mdtfb(t, 1)

proc mdt(t:float) =
  if useG:
    let t2 = 0.5*t
    mdtf t2
    mdtb t2
    #mdtf t
  else:
    mdtf t

proc mdv(t:float) =
  if useG:
    gc.gaugeforce2(g, f)
    initg()
    #var f2 = 0.0
    threads:
      #var f2t: typeof(norm2(f[0][0]))
      for i in 0..<g.len:
        for e in g[i]:
          let gd = gfunderiv(e,i)
          let tg = t*gd[0]
          var ff = -tg * f[i][e]  # - to correct for sign of f
          ff += t * gd[1]
          p[i][e] += ff
          #f2t += ff.norm2
          #var f2s = simdSum(f2t)
          #threadRankSum(f2s)
          #threadSingle: f2 = f2s
    #echo "F2: ", f2 / g[0].l.physVol
  else:
    gc.gaugeforce2(g, f)
    let etxi = exp(-0.5*t*xi)
    if mdalgo == nosehoover:
      threads:
        for i in 0..<p.len:
          p[i] *= etxi
    threads:
      for i in 0..<f.len:
        for e in f[i]:
          let tf = (-t)*f[i][e]
          p[i][e] += tf
    if mdalgo == nosehoover:
      threads:
        for i in 0..<p.len:
          p[i] *= etxi
      lnJ += dof*t*xi

####  end of area to modify to add new HMC variants  ###

proc updatefga(ts:openarray[float]) =
  tic("updatefga")
  mdv ts[0]
  toc("mdv")

proc revCheck(evo:auto; h0,ga0,t0,eh0:float, g0,p0:auto) =
  tic("reversibility")
  var
    g1 = lo.newgauge
    p1 = lo.newgauge
    #xi1 = xi
    lnJ1 = lnJ
  threads:
    for i in 0..<g1.len:
      g1[i] := g[i]
      p1[i] := p[i]
      p[i] := -1*p[i]
  #xi = -xi
  revBegin()
  evo.evolve tau
  evo.finish
  let
    (ga1, t1, h1) = hamiltonian(gc,g,p)
    #eh1 = xi * xi / (2.0 * gamma)
    eh1 = extraH()
    dH = h1-h0+eh1-eh0+lnJ
    dSg = ga1-ga0
    dT = t1-t0
    deh = eh1-eh0
  var dg2 = 0.0
  var dp2 = 0.0
  threads:
    for i in 0..<g.len:
      g[i] -= g0[i]
      let tg2 = g[i].norm2
      threadSingle: dg2 += tg2
      p[i] += p0[i]
      let tp2 = p[i].norm2
      threadSingle: dp2 += tp2
      g[i] := g1[i]
      p[i] := p1[i]
  dg2 = sqrt(dg2/float(lo.nDim * lo.physVol))
  dp2 = sqrt(dp2/float(lo.nDim * lo.physVol))
  qexLog "Reversed H: ",h1,"  Sg: ",ga1,"  T: ",t1,"  extH: ",eh1,"  lnJ: ",lnJ
  #qexLog "Reversibility: dH: ",dH,"  dSg: ",dSg,"  dT: ",dT,"  dextH: ",deh
  qexLog "Reversibility: dg2: ",dg2,"  dp2: ",dp2,"  dSg: ",dSg,"  dT: ",dT
  #proc epf(d,x:float):float = abs(d/x)/dof
  #if epf(dH,h0+eh0)>1e-14 or epf(dSg,ga0)>1e-14 or epf(dT,t0)>1e-14 or abs(lnJ/dof)>1e-14:
  #  qexWarn "broken reversibility in error/volume: dH: ",epf(dH,h0+eh0),"  dSg: ",epf(dSg,ga0),"  dT: ",epf(dT,t0),"  lnJ: ",abs(lnJ/dof)
  #xi = xi1
  revEnd()
  lnJ = lnJ1
  toc("done")

proc obstat(Hvals, Jvals, Avals, Pvals, Qvals:seq[float]) =
  proc rootmean(xs: Ensemble[seq[float]]):float = sqrt(mean(xs))
  let dHrms = Hvals.jackknife(jkBlockSize, rootmean)
  let lnJ = Jvals.jackknife(jkBlockSize, mean)
  var J2vals = newseq[float](Jvals.len)
  for i in 0..<Jvals.len:
    let a = Jvals[i]
    J2vals[i] = a*a
  let lnJrms = J2vals.jackknife(jkBlockSize, rootmean)
  let expmdh = Avals.jackknife(jkBlockSize, mean)
  var APvals = newseq[float](Avals.len)
  for i in 0..<Avals.len:
    let a = Avals[i]
    if a>1.0:
      APvals[i] = 1.0
    else:
      APvals[i] = a
  let pacc = APvals.jackknife(jkBlockSize, mean)
  let Pmean = Pvals.jackknife(jkBlockSize, mean)
  let Qmean = Qvals.jackknife(jkBlockSize, mean)
  let Qac = Qvals.jackknife(jkBlockSize, naiveIntAutocorr, 200)
  let dQrms = Qvals.jackknife(jkBlockSize, tunnelingRate)

  var qmax = 0
  for qv in Qvals:
    let q = int(round(qv))
    if q>qmax:
      qmax = q
    elif -q>qmax:
      qmax = -q
  var qdist = newseq[JackknifeStat[float]]()
  proc count(xs:Ensemble[seq[float]], x:int):float =
    for i in 0..<xs.len:
      if x == int(abs(round(xs[i]))):
        result += 1.0
    if x!=0:
      result *= 0.5
  for i in 0..qmax:
    qdist.add Qvals.jackknife(jkBlockSize, count, i)

  var Q2vals = newseq[float](Qvals.len)
  for i in 0..<Qvals.len:
    let a = Qvals[i]
    Q2vals[i] = a*a
  let Q2 = Q2vals.jackknife(jkBlockSize, mean)
  let Q2ac = Q2vals.jackknife(jkBlockSize, naiveIntAutocorr, 200)

  echo "lnJ = ", lnJ.mean, " ± ", lnJ.stdev
  echo "lnJrms = ", lnJrms.mean, " ± ", lnJrms.stdev
  echo "dHrms = ", dHrms.mean, " ± ", dHrms.stdev
  echo "exp(-dH) = ", expmdh.mean, " ± ", expmdh.stdev
  echo "Pacc = ", pacc.mean, " ± ", pacc.stdev
  echo "Pmean = ", Pmean.mean, " ± ", Pmean.stdev, "  dP = ",Pmean.mean-infVolPlaq(beta)
  echo "Qmean = ", Qmean.mean, " ± ", Qmean.stdev
  echo "Tau_Q = ", Qac.mean, " ± ", Qac.stdev
  echo "dQrms = ", dQrms.mean, " ± ", dQrms.stdev
  echo "Q2/V = ", Q2.mean/float(vol), " ± ", Q2.stdev/float(vol), "  dQ2/V = ", Q2.mean/float(vol)-infVolChiQ(beta)
  echo "Tau_Q2 = ", Q2ac.mean, " ± ", Q2ac.stdev
  for i in 0..qmax:
    echo "P(Q=",i,") = ",qdist[i].mean/float(ntraj), " ± ", qdist[i].stdev/float(ntraj)

#let
#  (V,T) = newIntegratorPair(updatefga, mdt)
#  md = gintalg(T = T, V = V[0], steps = gsteps)

type Md = object
var md: Md
proc evolve0(md: Md, t: float) =
  let eps = t / gsteps
  let eps2 = 0.5 * eps
  mdv eps2
  mdt eps
  for i in 1..<gsteps:
    mdv eps
    mdt eps
  mdv eps2
proc evolve(md: Md, t: float) =
  let eps = t / gsteps
  let epsa1 = omfLambda * eps
  let epsa1x2 = 2*epsa1
  let epsa2 = eps - epsa1x2
  let epsb1 = 0.5 * eps
  mdv epsa1
  mdtf epsb1
  mdv epsa2
  mdtb epsb1
  for i in 1..<gsteps:
    mdv epsa1x2
    mdtf epsb1
    mdv epsa2
    mdtb epsb1
  mdv epsa1

proc finish(md: Md) = discard

proc mc =
  tic("mc")
  var
    Hvals = newSeq[float](ntraj)
    Jvals = newSeq[float](ntraj)
    Avals = newSeq[float](ntraj)
    Pvals = newSeq[float](ntraj)
    Qvals = newSeq[float](ntraj)

  for n in 1-ntrajThermo..ntraj:
    tic("traj")
    qexLog "Begin traj: ",n
    threads:
      p.randomTAH r
      for i in 0..<p.len:
        g0[i] := g[i]
        p0[i] := p[i]
    #if mdalgo == nosehoover:
    #  xi = R.gaussian * sqrt(gamma)
    #  lnJ = 0.0
    extraInit()
    let (ga0, t0, h0) = hamiltonian(gc,g,p)
    #let eh0 = xi * xi / (2.0 * gamma)
    let eh0 = extraH()
    qexLog "Begin H: ",h0,"  Sg: ",ga0,"  T: ",t0,"  extH: ",eh0
    md.evolve tau
    md.finish
    let (ga1, t1, h1) = hamiltonian(gc,g,p)
    #let eh1 = xi * xi / (2.0 * gamma)
    let eh1 = extraH()
    qexLog "End H: ",h1,"  Sg: ",ga1,"  T: ",t1,"  extH: ",eh1,"  lnJ: ",lnJ

    if revCheckFreq > 0 and n>=0 and n mod revCheckFreq == 0:
      md.revCheck(h0,ga0,t0,eh0,g0,p0)

    let
      dH = h1 - h0 + eh1 - eh0 + lnJ
      acc = exp(-dH)
      accr = R.uniform
    if n <= ntrajThermoAcc-ntrajThermo:
      qexLog "ACCEPT(FORCE):  dH: ",dH,"  exp(-dH): ",acc,"  r: ",accr
    elif accr <= acc:  # accept
      qexLog "ACCEPT:  dH: ",dH,"  exp(-dH): ",acc,"  r: ",accr
    else:  # reject
      qexLog "REJECT:  dH: ",dH,"  exp(-dH): ",acc,"  r: ",accr
      threads:
        for i in 0..<g.len:
          g[i] := g0[i]

    let pl = g.plaq3
    if n>0:
      Hvals[n-1] = dH*dH
      Avals[n-1] = acc
      Jvals[n-1] = lnJ
      Pvals[n-1] = pl.re
      Qvals[n-1] = g.topo2DU1
    qexLog "plaq: ",pl.re," ",pl.im," topo: ",Qvals[n-1]
    toc("done")

  obstat(Hvals, Jvals, Avals, Pvals, Qvals)
  toc("done")

toc("prep")
mc()
toc("hmc")

echoProf()
qexGC()
qexfinalize()
