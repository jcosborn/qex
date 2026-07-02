## 2D U(1) pure-gauge HMC using graph-based evolution.

import base/globals
setDefaultNc(1)
setVLENmax(4)

import qex

import core
import scalar
import gauge
import hmcgauge/config
import hmcgauge/trajectory
import hmcgauge/rng
import hmcgauge/gauge_io
import hmcgauge/params
import hmcgauge/measure2du1

qexInit()
echo "rank ", myRank, "/", nRanks
threads:
  echo "thread ", threadNum, "/", numThreads

let
  gp = readGaugeInputs(GaugeParams(
    lat: @[64, 64], beta: 6.5, seed: 1234567891, rng: rkPhilox4x64,
    savefile: "config", dt: 0.2, gsteps: 10, intalg: "2MN",
    trajsThermo: 64, trajs: 128, revCheckFreq: 0))
  runConfig = gp.toRunConfig
runConfig.validateRunConfig
installStandardParams()
echoParams()
processHelpParam()

let
  lo = gp.lat.newLayout
  grt = initGraphRuntime()
  gc = actWilson(scalar.toGvalue(grt, gp.beta))

withRng(gp.rng, R):
  var
    random = lo.newRNGField(R, gp.seed)
    acceptRandom: R
  acceptRandom.seed(gp.seed, 987654321)

  var
    g = lo.newgauge
    p = lo.newgauge
  if runConfig.gaugefile.len == 0:
    g.random random          # hot start
  else:
    g.loadOrInitGauge runConfig.gaugefile

  echo "Initial plaq: ", g.plaq3

  let
    action = proc(x: Ggauge): Gscalar = gaugeAction(gc, x)
    force = proc(x: Ggauge): Ggauge = gaugeForce(gc, x)
    graph = buildTrajectoryGraph(grt, g, p, action, runConfig, buildTraining = false, force = force)

  var
    Hvals = newSeq[float](runConfig.trajs)   # dH of each measured trajectory
    Avals = newSeq[float](runConfig.trajs)   # exp(-dH)
    Pvals = newSeq[float](runConfig.trajs)   # average plaquette (real part)
    Qvals = newSeq[float](runConfig.trajs)   # topological charge
    dQchanged = newSeq[bool](runConfig.trajs)  # proposal changed Q? (per production traj)
    mdvals = newSeq[MdForceStats](runConfig.trajs)

  # Track topology on the committed configuration.
  var
    prevQ = g.topo2DU1   # committed topological charge from the previous trajectory
    proposalQ = prevQ

  # dQ = Q(proposal) - Q(committed), before accept/reject.
  proc proposalMon(traj: int; dH, acc: float) =
    discard graph.finalState.gauge.eval
    let tm = graph.finalState.gauge.gaugeSnapshot.topoMaxP2DU1
    let dq = int(round(tm.topo - prevQ))
    proposalQ = tm.topo
    echo "proposal: dQ ", dq, "  maxP ", tm.maxP
    if traj > runConfig.trajsThermo:
      dQchanged[traj - runConfig.trajsThermo - 1] = dq != 0

  proc measureTraj(traj: int; dH, acc: float; accepted: bool; forceStats: MdForceStats) =
    let
      cur = graph.initialState.gauge.gaugeSnapshot
      pl = cur.plaq3
      q = if accepted: proposalQ else: prevQ
    echo "plaq: ", pl.re, "  topo: ", q
    prevQ = q
    if traj > runConfig.trajsThermo:
      let i = traj - runConfig.trajsThermo - 1
      Hvals[i] = dH
      Avals[i] = acc
      Pvals[i] = pl.re
      Qvals[i] = q
      mdvals[i] = forceStats
      cur.maybeSaveGauge(runConfig, traj)

  runHmc(graph, runConfig, random, acceptRandom, measureTraj, proposalMon)

  obstat(Hvals, Avals, Pvals, Qvals, gp.beta, lo.physVol, runConfig.trajs, gp.jkBlockSize, mdvals = mdvals)
  statsByQ(Hvals, Avals, dQchanged, gp.jkBlockSize, mdvals = mdvals)

if gp.showTimers: echoTimers()
if gp.showRunStats: grt.echoRunStats
processSaveParams()
writeParamFile()
qexFinalize()
