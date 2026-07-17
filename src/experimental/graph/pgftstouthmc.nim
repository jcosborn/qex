## 4D SU(3) stout FTHMC:
##   U = f(V),  S_eff(V) = S(U) - log det f'(V).
## The per-link log-Jacobian uses the SU(3) adjoint representation.
## Measure and save U; invert loaded U to V before sampling.

import qex
import sequtils
import utils/resample

import core
import scalar
import gauge
import hmcgauge/config
import hmcgauge/trajectory
import hmcgauge/params
import hmcgauge/rng
import hmcgauge/gauge_io
import hmcgauge/ftstout
import hmcgauge/stats

proc avgPlaq(g: auto): float =
  ## <Re tr U_plaq/N>.
  g.plaq.sum

proc ploops(g: auto): tuple[spatial, temporal: float] =
  ## Polyakov-loop magnitudes: spatial-average and temporal (real part of trace).
  let pg = g[0].l.physGeom
  var lines = newSeq[seq[int]](pg.len)
  for i in 0..<pg.len:
    lines[i] = repeat(i+1, pg[i])
  let pl = g.wilsonLines(lines)
  var pls = pl[0]
  for i in 1..<pl.len-1:
    pls += pl[i]
  (spatial: (pls / float(pl.len-1)).re, temporal: pl[^1].re)

proc obstat(Hvals, Avals, Pvals, Lvals, Jvals: seq[float]; mdvals: seq[MdForceStats]; jkBlockSize: int) =
  let
    dHrms = Hvals.jackknife(jkBlockSize, rms)
    expmdh = Avals.jackknife(jkBlockSize, mean)
    pmean = Pvals.jackknife(jkBlockSize, mean)
    lmean = Lvals.jackknife(jkBlockSize, mean)
    lnDet = Jvals.jackknife(jkBlockSize, mean)
    he = extrema(Hvals)
  proc err(x: JackknifeStat[float]): string =
    if x.hasStdev: $x.stdev else: "n/a"
  echo "ntraj: ", Hvals.len
  echo "dHrms: ", dHrms.mean, " ± ", err(dHrms)
  echo "dH min/max: ", he.lo, " / ", he.hi
  echo "exp(-dH): ", expmdh.mean, " ± ", err(expmdh)
  echo "Plaq: ", pmean.mean, " ± ", err(pmean)
  echo "Ploop(temporal): ", lmean.mean, " ± ", err(lmean)
  echo "lnDet: ", lnDet.mean, " ± ", err(lnDet)
  echoMdStats(mdvals, jkBlockSize)

qexInit()
echo "rank ", myRank, "/", nRanks
threads:
  echo "thread ", threadNum, "/", numThreads

let
  gp = readGaugeInputs(GaugeParams(
    lat: @[8, 8, 8, 16], beta: 5.4, seed: 1234567891, rng: rkPhilox4x64,
    savefile: "config", dt: 0.025, gsteps: 4, intalg: "2MN", trajs: 50))
  sp = readStoutParams(StoutParams(rho: 0.1, nsmear: 1))
  runConfig = gp.toRunConfig
runConfig.validateRunConfig
if gp.lat.len != 4:
  raiseValueError("pgftstouthmc requires a 4D lattice, got " & $gp.lat.len & " dimensions")
for n in gp.lat:
  if (n and 1) != 0:
    raiseValueError("pgftstouthmc requires even lattice extents")
if sp.nsmear < 0:
  raiseValueError("nsmear must be >= 0, got " & $sp.nsmear)
# One-link contraction: 8|epsilon| < 1, epsilon = rho/3.
if abs(sp.rho) >= 3.0/8.0:
  qexWarn "abs(rho) >= 3/8; the 4D SU(3) stout Jacobian is not guaranteed positive definite: rho ", sp.rho
installStandardParams()
echoParams()
processHelpParam()

proc runStoutHmc[T](R: typedesc[T]) =
  tic()
  let
    lo = gp.lat.newLayout
    grt = initGraphRuntime()
    gc = actWilson(scalar.toGvalue(grt, gp.beta))
    sa = stoutAction(gc, sp.rho, sp.nsmear)

  var
    randomField = lo.newRNGField(R, gp.seed)
    acceptRandom: T
    g = lo.newgauge
    p = lo.newgauge
    Uloaded: type(g)                 # physical config as loaded (empty on a hot start)
  acceptRandom.seed(gp.seed, 987654321)
  if runConfig.gaugefile.len == 0:
    g.random randomField          # hot start (fundamental field V)
  else:
    g.loadOrInitGauge runConfig.gaugefile   # the saved file holds the physical U = f(V)
    Uloaded = lo.newgauge
    threads:
      for mu in 0..<g.len: Uloaded[mu] := g[mu]
    let inv = invertStoutFlow(g, sp.rho, sp.nsmear)   # map back to the fundamental V
    echo "stout inverse U->V: iter ", inv.iter, "  rdf2 ", inv.rdf2

  let
    graph = buildTrajectoryGraph(grt, g, p, sa.action, runConfig, buildTraining = false)
    measure = sa.flow(graph.initialState.gauge)

  block:
    discard measure.smeared.eval
    let us = measure.smeared.gaugeSnapshot
    echo "Initial smeared plaq: ", us.avgPlaq
    if Uloaded.len > 0:   # f(f^-1(U)) must reproduce the loaded physical config
      echo "load round-trip |f(f^-1(U)) - U|_max^2: ", maxGaugeDiff2(us, Uloaded)

  var
    Hvals = newSeq[float](runConfig.trajs)   # dH
    Avals = newSeq[float](runConfig.trajs)   # exp(-dH)
    Pvals = newSeq[float](runConfig.trajs)   # average plaquette of U = f(V)
    Lvals = newSeq[float](runConfig.trajs)   # temporal Polyakov loop of U = f(V)
    Jvals = newSeq[float](runConfig.trajs)   # ln det f'(V)
    mdvals = newSeq[MdForceStats](runConfig.trajs)

  # Measure on the physical field U = f(V) at the committed configuration.
  proc measureTraj(traj: int; dH, acc: float; accepted: bool; forceStats: MdForceStats) =
    discard measure.smeared.eval
    let
      u = measure.smeared.gaugeSnapshot
      lndetCur = measure.lndet.eval.sval
      pl = u.avgPlaq
      lp = u.ploops
    echo "plaq: ", pl, "  ploop: ", lp.spatial, " ", lp.temporal, "  lnDet: ", lndetCur
    if traj > runConfig.trajsThermo:
      let i = traj - runConfig.trajsThermo - 1
      Hvals[i] = dH
      Avals[i] = acc
      Pvals[i] = pl
      Lvals[i] = lp.temporal
      Jvals[i] = lndetCur
      mdvals[i] = forceStats
      u.maybeSaveGauge(runConfig, traj)   # save the physical field U = f(V)

  runHmc(graph, runConfig, randomField, acceptRandom, measureTraj)

  if Hvals.len > 0:
    obstat(Hvals, Avals, Pvals, Lvals, Jvals, mdvals, gp.jkBlockSize)

  toc()
  if gp.showTimers: echoTimers()
  if gp.showRunStats: grt.echoRunStats

withRng(gp.rng, R):
  runStoutHmc(R)

processSaveParams()
writeParamFile()
qexFinalize()
