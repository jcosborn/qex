import qex, gauge, physics/qcdTypes, algorithms/integrator
import math, os, sequtils, strformat, strutils, times

type GaugeActType = enum ActWilson, ActAdjoint, ActRect, ActSymanzik, ActIwasaki, ActDBW2
converter toGaugeActType(s:string):GaugeActType = parseEnum[GaugeActType](s)

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

proc mplaq(g:auto) =
  tic()
  let
    pl = g.plaq
    nl = pl.len div 2
    ps = pl[0..<nl].sum * 2.0
    pt = pl[nl..^1].sum * 2.0
  echo "MEASplaq ss: ",ps,"  st: ",pt,"  tot: ",0.5*(ps+pt)
  toc("plaq")

proc ploop(g:auto) =
  tic()
  let pg = g[0].l.physGeom
  var pl = newseq[typeof(g.wline @[1])](pg.len)
  for i in 0..<pg.len:
    pl[i] = g.wline repeat(i+1, pg[i])
  let
    pls = pl[0..^2].sum / float(pl.len-1)
    plt = pl[^1]
  echo "MEASploop spatial: ",pls.re," ",pls.im," temporal: ",plt.re," ",plt.im
  toc("ploop")

qexinit()
tic()

letParam:
  gaugefile = ""
  savefile = "config"
  savefreq = 10
  lat =
    if fileExists(gaugefile):
      getFileLattice gaugefile
    else:
      if gaugefile.len > 0:
        qexWarn "Nonexistent gauge file: ", gaugefile
      @[8,8,8,8]
  gact:GaugeActType = "ActAdjoint"
  beta = 6.0
  adjFac = -0.25
  rectFac = -1.4088
  tau = 2.0
  inittraj = 0
  trajs = 10
  seed:uint64 = int(1000*epochTime())
  gintalg:IntegratorProc = "4MN5F2GP"
  gsteps = 4
  alwaysAccept:bool = 0
  revCheckFreq = savefreq
  useFG2:bool = 0
  showTimers:bool = 1
  timerWasteRatio = 0.05
  timerEchoDropped:bool = 0
  timerExpandRatio = 0.05
  verboseGCStats:bool = 0
  verboseTimer:bool = 0

echoParams()
echo "rank ", myRank, "/", nRanks
threads: echo "thread ", threadNum, "/", numThreads

DropWasteTimerRatio = timerWasteRatio
VerboseGCStats = verboseGCStats
VerboseTimer = verboseTimer

let
  gc = case gact
    of ActWilson: GaugeActionCoeffs(plaq: beta)
    of ActAdjoint: GaugeActionCoeffs(plaq: beta, adjplaq: beta*adjFac)
    of ActRect: gaugeActRect(beta, rectFac)
    of ActSymanzik: Symanzik(beta)
    of ActIwasaki: Iwasaki(beta)
    of ActDBW2: DBW2(beta)
  lo = lat.newLayout
  vol = lo.physVol

echo gc

var r = lo.newRNGField(MRG32k3a, seed)
var R:MRG32k3a  # global RNG
R.seed(seed, 987654321)

var
  g = lo.newgauge
  p = lo.newgauge
  f = lo.newgauge
  g0 = lo.newgauge
  gg = lo.newgauge  # FG backup gauge

proc fgsave =
  threads:
    for mu in 0..<g.len:
      gg[mu] := g[mu]
proc fgload =
  threads:
    for mu in 0..<g.len:
      g[mu] := gg[mu]

template pnorm2(p2:float) =
  threads:
    var p2t = 0.0
    for i in 0..<p.len:
      p2t += p[i].norm2
    threadMaster: p2 = p2t

proc gaction(g:auto, p2:float):auto =
  tic()
  let
    ga = if gact==ActAdjoint: gc.actionA g else: gc.gaugeAction1 g
    t = 0.5*p2 - float(16*vol)
    h = ga + t
  toc("gaction")
  (ga, t, h)

proc mdt(t: float) =
  tic()
  threads:
    for mu in 0..<g.len:
      for s in g[mu]:
        g[mu][s] := exp(t*p[mu][s])*g[mu][s]
  toc("mdt")
proc mdv(t: float) =
  tic()
  if gact==ActAdjoint:
    gc.forceA(g, f)
  else:
    gc.gaugeForce(g, f)
  qexGC "mdv forceA"
  threads:
    for mu in 0..<f.len:
      p[mu] -= t*f[mu]
  toc("mdv")
# FG update g from backup, gg
proc fgv(t: float) =
  tic()
  if gact==ActAdjoint:
    gc.forceA(gg, f)
  else:
    gc.gaugeForce(gg, f)
  qexGC "fgv forceA"
  threads:
    for mu in 0..<g.len:
      for s in g[mu]:
        g[mu][s] := exp((-t)*f[mu][s])*g[mu][s]
  toc("fgv")
# Combined update for sharing computations
proc mdvAllfga(ts,gs:openarray[float]) =
  tic("mdvAllfga")
  var
    updateG = false
    updateGG = false
  let approxOrder = if useFG2: 2 else: 1
  type GGstep = tuple[t,g:float]
  var ggs: array[2,GGstep]

  # gauge
  if gs[0] != 0:
    updateGG = true
    if ts[0] == 0:
      qexError "Force gradient without the force update."
    if useFG2:
      let (tf,tg) = approximateFGcoeff2(ts[0],gs[0])
      for o in 0..1:
        ggs[o] = (t:tf[o], g:tg[o])
    else:
      let (tf,tg) = approximateFGcoeff(ts[0],gs[0])
      ggs[0] = (t:tf, g:tg)
  elif ts[0] != 0:
    updateG = true

  if updateGG:
    tic()
    fgsave()
    toc("FG save done")

  # MD
  if updateG:
    tic("MD mdv")
    mdv ts[0]
    toc("done")

  # FG
  if updateGG:
    tic("updateFG")
    for o in 0..<approxOrder:
      tic()
      fgv ggs[o].g
      toc("fgv")
      mdv ggs[o].t
      toc("FG mdv")
      fgload()
      toc("load")
    toc("done")
  toc("done")

proc revCheck(evo:auto; h0,ga0,t0:float) =
  tic("reversibility")
  var
    g1 = lo.newgauge
    p1 = lo.newgauge
    p2 = 0.0
  threads:
    for i in 0..<g1.len:
      g1[i] := g[i]
      p1[i] := p[i]
      p[i] := -1*p[i]
  evo.evolve tau
  evo.finish
  p2.pnorm2
  toc("p norm2 2")
  let (ga1,t1,h1) = g.gaction p2
  qexLog "Reversed H: ",h1,"  Sg: ",ga1,"  T: ",t1
  echo "Reversibility: dH: ",h1-h0,"  dSg: ",ga1-ga0,"  dT: ",t1-t0
  for i in 0..<g1.len:
    g[i] := g1[i]
    p[i] := p1[i]
  qexGC "revCheck done"
  toc("done")

let
  (V,T) = newIntegratorPair(mdvAllfga, mdt)
  H = gintalg(T = T, V = V[0], steps = gsteps)

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

g.mplaq

echo H

toc("prep")

for n in inittraj+1..inittraj+trajs:
  tic("traj")
  var p2 = 0.0
  threads:
    p.randomTAH r
    for i in 0..<g.len:
      g0[i] := g[i]
  toc("p refresh, save g")
  p2.pnorm2
  toc("p norm2 1")
  let (ga0,t0,h0) = g0.gaction p2
  qexGC "init action"
  toc("init gauge action")
  qexLog "Begin H: ",h0,"  Sg: ",ga0,"  T: ",t0

  H.evolve tau
  H.finish
  toc("evolve")

  p2.pnorm2
  toc("p norm2 2")
  let (ga1,t1,h1) = g.gaction p2
  qexGC "end action"
  toc("final gauge action")
  qexLog "End H: ",h1,"  Sg: ",ga1,"  T: ",t1
  toc("end evolve")

  if revCheckFreq > 0 and n mod revCheckFreq == 0:
    H.revCheck(h0,ga0,t0)

  let
    dH = h1 - h0
    acc = exp(-dH)
    accr = R.uniform
  if accr <= acc or alwaysAccept:  # accept
    echo "ACCEPT:  dH: ",dH,"  exp(-dH): ",acc,"  r: ",accr,(if alwaysAccept:" (ignored)" else:"")
    g.reunit
  else:  # reject
    echo "REJECT:  dH: ",dH,"  exp(-dH): ",acc,"  r: ",accr
    threads:
      for i in 0..<g.len:
        g[i] := g0[i]
  qexGC "traj done"

  g.mplaq
  g.ploop
  qexGC "measure done"
  toc("measure")

  if savefreq > 0 and n mod savefreq == 0:
    tic("save")
    let fn = savefile & &".{n:05}.lime"
    if 0 != g.saveGauge(fn):
      qexError "Failed to save gauge to file: ",fn
    qexLog "saved gauge to file: ",fn," secs: ",getElapsedTime()
    toc("done")

  qexLog "traj ",n," secs: ",getElapsedTime()
  toc("traj end")

toc("hmc")

if showTimers: echoTimers(timerExpandRatio, timerEchoDropped)
qexfinalize()
