## 2D U(1) stout FTHMC:
##   U = f(V),  S_eff(V) = S(U) - log det f'(V).
## Measure and save U; invert loaded U to V before sampling.

import base/globals
setDefaultNc(1)
setVLENmax(4)

import qex

import core
import scalar
import gauge
import hmcgauge/config
import hmcgauge/trajectory
import hmcgauge/params
import hmcgauge/rng
import hmcgauge/gauge_io
import hmcgauge/ftstout
import hmcgauge/measure2du1

# The field transformation contributes log det f'(V); obstat reports its
# lnDet/lnDetrms when handed the per-trajectory Jvals.

qexInit()
echo "rank ", myRank, "/", nRanks
threads:
  echo "thread ", threadNum, "/", numThreads

let
  gp = readGaugeInputs(GaugeParams(
    lat: @[64, 64], beta: 6.5, seed: 1234567891, rng: rkPhilox4x64,
    savefile: "config", dt: 0.1, gsteps: 10, intalg: "2MN",
    trajsThermo: 64, trajs: 128, revCheckFreq: 0))
  sp = readStoutParams(StoutParams(rho: 0.1, nsmear: 1))
  runConfig = gp.toRunConfig
runConfig.validateRunConfig
if gp.lat.len != 2:
  raiseValueError("pg2du1ftstouthmc requires a 2D lattice, got " & $gp.lat.len & " dimensions")
for n in gp.lat:
  if (n and 1) != 0:
    raiseValueError("pg2du1ftstouthmc requires even lattice extents")
if sp.nsmear < 0:
  raiseValueError("nsmear must be >= 0, got " & $sp.nsmear)
# J = 1 - rho(cos p+ + cos p-) > 0 for every link.
if abs(sp.rho) >= 0.5:
  qexWarn "abs(rho) >= 1/2; the 2D U(1) stout Jacobian is not guaranteed positive: rho ", sp.rho
installStandardParams()
echoParams()
processHelpParam()

let
  lo = gp.lat.newLayout
  grt = initGraphRuntime()
  gc = actWilson(scalar.toGvalue(grt, gp.beta))
  sa = stoutAction(gc, sp.rho, sp.nsmear)

withRng(gp.rng, R):
  var
    random = lo.newRNGField(R, gp.seed)
    acceptRandom: R
  acceptRandom.seed(gp.seed, 987654321)

  var
    g = lo.newgauge
    p = lo.newgauge
    Uloaded: type(g)                 # physical config as loaded (empty on a hot start)
  if runConfig.gaugefile.len == 0:
    g.random random          # hot start (fundamental field V)
  else:
    g.loadOrInitGauge runConfig.gaugefile   # the saved file holds the physical U = f(V)
    Uloaded = lo.newgauge
    threads:
      for mu in 0..<g.len: Uloaded[mu] := g[mu]
    let inv = invertStoutFlow(g, sp.rho, sp.nsmear)   # map back to the fundamental V
    echo "stout inverse U->V: iter ", inv.iter, "  rdf2 ", inv.rdf2

  let
    graph = buildTrajectoryGraph(grt, g, p, sa.action, runConfig, buildTraining = false)
    # physical field U = f(V) and its log-Jacobian, for measurement on the V leaf
    measure = sa.flow(graph.initialState.gauge)
    measureLd = logDetJ(measure, graph.initialState.gauge)
    # proposed physical field U = f(V_prop) from the end-of-MD state, for the
    # proposed-configuration monitor (dQ, maxP)
    measureProp = sa.flow(graph.finalState.gauge)

  var
    prevQ = 0.0   # committed topological charge from the previous trajectory
    proposalQ = 0.0
  block:
    discard measure.eval
    let us = measure.gaugeSnapshot
    echo "Initial smeared plaq: ", us.plaq3
    prevQ = us.topo2DU1
    if Uloaded.len > 0:   # f(f^-1(U)) must reproduce the loaded physical config
      echo "load round-trip |f(f^-1(U)) - U|_max^2: ", maxGaugeDiff2(us, Uloaded)

  var
    Hvals = newSeq[float](runConfig.trajs)   # dH of each measured trajectory
    Jvals = newSeq[float](runConfig.trajs)   # ln det f'(V)
    Avals = newSeq[float](runConfig.trajs)   # exp(-dH)
    Pvals = newSeq[float](runConfig.trajs)   # average plaquette of U = f(V)
    Qvals = newSeq[float](runConfig.trajs)   # topological charge of U = f(V)
    dQchanged = newSeq[bool](runConfig.trajs)  # proposal changed Q? (per production traj)
    mdvals = newSeq[MdForceStats](runConfig.trajs)

  # Trajectories up to trajsThermo are thermalization (discarded); the first
  # trajsForceAcc are force-accepted to escape the hot start.
  # Measure on the physical field U = f(V) at the committed configuration.
  # Proposed-configuration monitor (pre-commit): dQ = Q(f(V_prop)) − Q(committed) and
  # max|plaquette angle| of the end-of-MD proposal, independent of accept/reject.
  proc proposalMon(traj: int; dH, acc: float) =
    discard measureProp.eval
    let tm = measureProp.gaugeSnapshot.topoMaxP2DU1
    let dq = int(round(tm.topo - prevQ))
    proposalQ = tm.topo
    echo "proposal: dQ ", dq, "  maxP ", tm.maxP
    if traj > runConfig.trajsThermo:
      dQchanged[traj - runConfig.trajsThermo - 1] = dq != 0

  proc measureTraj(traj: int; dH, acc: float; accepted: bool; forceStats: MdForceStats) =
    discard measure.eval
    let
      u = measure.gaugeSnapshot
      lndetCur = measureLd.eval.sval
      pl = u.plaq3
      q = if accepted: proposalQ else: prevQ
    echo "plaq: ", pl.re, "  topo: ", q, "  lnDet: ", lndetCur
    prevQ = q
    if traj > runConfig.trajsThermo:
      let i = traj - runConfig.trajsThermo - 1
      Hvals[i] = dH
      Avals[i] = acc
      Jvals[i] = lndetCur
      Pvals[i] = pl.re
      Qvals[i] = q
      mdvals[i] = forceStats
      u.maybeSaveGauge(runConfig, traj)   # save the physical field U = f(V)

  runHmc(graph, runConfig, random, acceptRandom, measureTraj, proposalMon)

  if Hvals.len > 0:
    obstat(Hvals, Avals, Pvals, Qvals, gp.beta, lo.physVol, runConfig.trajs, gp.jkBlockSize, Jvals, mdvals)
    statsByQ(Hvals, Avals, dQchanged, gp.jkBlockSize, Jvals, mdvals)

if gp.showTimers: echoTimers()
if gp.showRunStats: grt.echoRunStats
processSaveParams()
writeParamFile()
qexFinalize()
