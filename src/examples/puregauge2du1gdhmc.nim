import qex
import gauge, physics/qcdTypes
import algorithms/integrator, maths/special, utils/resample
import os, strutils

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

type MDAlgo = enum hamilton, nosehoover
converter toMDAlgo(s:string):MDAlgo = parseEnum[MDAlgo](s)

qexinit()

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
      @[64,64]
  beta = 5.0
  tau = 2.0
  ntraj = 128
  ntrajThermo = 64
  ntrajThermoAcc = ntrajThermo div 8
  jkBlockSize = max(1, ntraj div 64)
  seed:uint64 = 11u^11
  ## hamilton | nosehoover [H + xi^2/(2 gamma); dU = dt P, dP = - dt S'(U) - xi P, dxi = gamma dt (p^2 - dof)]
  mdalgo:MDAlgo = "hamilton"
  ## extra params in other mdalgo
  gamma = 1.0
  ## 2MN,0.21 | 4MN3F1GP,0.27 | 4MN5F2GP
  gintalg:IntegratorProc = "2MN,0.21"
  gsteps = 20
  alwaysAccept:bool = 0
  revCheckFreq = ntraj
  verboseGCStats:bool = 0
  verboseTimer:bool = 0

installStandardParams()
echoParams()
qexLog "rank ", myRank, "/", nRanks
threads: qexLog "thread ", threadNum, "/", numThreads
processHelpParam()

VerboseGCStats = verboseGCStats
VerboseTimer = verboseTimer

let
  lo = lat.newLayout
  gc = GaugeActionCoeffs(plaq:beta)
  vol = lo.physVol
  dof = float(2*vol)

var
  g = lo.newGauge
  r = lo.newRNGField(RngMilc6, seed)
  R:RngMilc6  # global RNG
R.seed(seed, 987654321)

g.random r

qexLog "Initial plaq: ",g.plaq3

var
  p = lo.newgauge
  f = lo.newgauge
  gg = lo.newgauge  # FG backup gauge
  g0 = lo.newgauge
  xi = 0.0
  lnJ = 0.0

proc mdt(t:float) =
  threads:
    for i in 0..<g.len:
      for e in g[i]:
        let etpg = exp(t*p[i][e])*g[i][e]
        g[i][e] := etpg
  if mdalgo == nosehoover:
    xi += t * gamma * (p.pnorm2 - dof)
proc mdv(t:float) =
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

# For force gradient update
proc fgv(t:float) =
  if mdalgo == nosehoover:
    qexError "force gradient update for nosehoover unimplemented."
  gc.gaugeforce2(g, f)
  threads:
    for i in 0..<g.len:
      for e in g[i]:
        let etfg = exp((-t)*f[i][e])*g[i][e]
        g[i][e] := etfg
proc fgsave =
  threads:
    for i in 0..<g.len:
      gg[i] := g[i]
proc fgload =
  threads:
    for i in 0..<g.len:
      g[i] := gg[i]
proc updatefga(ts,gs:openarray[float]) =
  tic("updatefga")
  let
    t = ts[0]
    g = gs[0]
  if g != 0:
    if t != 0:
      tic("FG")
      # Approximate the force gradient update with a Taylor expansion.
      let (tf,tg) = approximateFGcoeff(t,g)
      fgsave()
      toc("FG fgv")
      fgv tg
      toc("FG mdv")
      mdv tf
      toc("FG load")
      fgload()
      toc("done")
    else:
      quit("Force gradient without the force update.")
  elif t != 0:
    tic("MD mdv")
    mdv t
    toc("done")
  else:
    quit("No updates required.")
  toc("done")

proc revCheck(evo:auto; h0,ga0,t0,eh0:float) =
  tic("reversibility")
  var
    g1 = lo.newgauge
    p1 = lo.newgauge
    xi1 = xi
    lnJ1 = lnJ
  threads:
    for i in 0..<g1.len:
      g1[i] := g[i]
      p1[i] := p[i]
      p[i] := -1*p[i]
  xi = -xi
  evo.evolve tau
  evo.finish
  let
    (ga1, t1, h1) = hamiltonian(gc,g,p)
    eh1 = xi * xi / (2.0 * gamma)
    dH = h1-h0+eh1-eh0+lnJ
    dSg = ga1-ga0
    dT = t1-t0
    deh = eh1-eh0
  qexLog "Reversed H: ",h1,"  Sg: ",ga1,"  T: ",t1,"  extH: ",eh1,"  lnJ: ",lnJ
  qexLog "Reversibility: dH: ",dH,"  dSg: ",dSg,"  dT: ",dT,"  dextH: ",deh
  proc epf(d,x:float):float = abs(d/x)/dof
  if epf(dH,h0+eh0)>1e-14 or epf(dSg,ga0)>1e-14 or epf(dT,t0)>1e-14 or abs(lnJ/dof)>1e-14:
    qexWarn "broken reversibility in error/volume: dH: ",epf(dH,h0+eh0),"  dSg: ",epf(dSg,ga0),"  dT: ",epf(dT,t0),"  lnJ: ",abs(lnJ/dof)
  for i in 0..<g1.len:
    g[i] := g1[i]
    p[i] := p1[i]
  xi = xi1
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

let
  (V,T) = newIntegratorPair(updatefga, mdt)
  md = gintalg(T = T, V = V[0], steps = gsteps)

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
      for i in 0..<p.len:
        g0[i] := g[i]
      p.randomTAH r
    if mdalgo == nosehoover:
      xi = R.gaussian * sqrt(gamma)
      lnJ = 0.0
    let (ga0, t0, h0) = hamiltonian(gc,g,p)
    let eh0 = xi * xi / (2.0 * gamma)
    qexLog "Begin H: ",h0,"  Sg: ",ga0,"  T: ",t0,"  extH: ",eh0
    md.evolve tau
    md.finish
    let (ga1, t1, h1) = hamiltonian(gc,g,p)
    let eh1 = xi * xi / (2.0 * gamma)
    qexLog "End H: ",h1,"  Sg: ",ga1,"  T: ",t1,"  extH: ",eh1,"  lnJ: ",lnJ

    if revCheckFreq > 0 and n>=0 and n mod revCheckFreq == 0:
      md.revCheck(h0,ga0,t0,eh0)

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
