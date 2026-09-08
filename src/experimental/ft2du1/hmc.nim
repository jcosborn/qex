## Exact finite change-of-variables HMC for two-dimensional U(1).
## The sampled field is V with S_eff(V)=S(f(V))-ln det f'(V); measurements and
## saved configurations use the physical field U=f(V).

import base/globals
setDefaultNc(1)
setVLENmax(4)

import qex
import std/strformat

import scan

import ../graph/core
import ../graph/scalar
import ../graph/gauge
import ../graph/gauge/shared
import ../graph/gauge/action/ops
from ../graph/hmcgauge/config import validateRunConfig
import ../graph/hmcgauge/trajectory
import ../graph/hmcgauge/params
import ../graph/hmcgauge/rng
import ../graph/hmcgauge/gauge_io
import ../graph/hmcgauge/measure2du1

when isMainModule:
  qexInit()
  echo "rank ", myRank, "/", nRanks
  threads:
    echo "thread ", threadNum, "/", numThreads

  let
    gp = readGaugeInputs(GaugeParams(
      lat: @[16, 16], beta: 6.0, seed: 1234567891, rng: rkPhilox4x64,
      savefile: "config", dt: 0.1, gsteps: 10, intalg: "2MN",
      trajsThermo: 64, trajs: 128, revCheckFreq: 0))
    mp = readMapParams(MapParams(
      geometry: "plaq4", construction: "scalar", basis: "bspline",
      ctxBasis: "fourier",
      mapDepth: 1, flowDepth: 1, mapStrengths: @[0.7],
      mapFloor: 0.15, mapEpsilon: 0.05, mapRho: 0.1, mapOrder: 8,
      ctxOrder: 1, mapKnots: 8, ctxKnots: 64,
      mapCoeffs: @[], mapControls: @[],
      mapParamFile: "", fejerOrders: @[8], fejerCenters: @[],
      fejerWeights: @[], mapStageOrder: "0123",
      mapDirs: @[0], mapParities: @[0], mapOffsets: @[0, 0],
      mapStride: 0, mapInvTol: 2e-14, mapScanStep: 2e-5, mapInvIter: 80,
      mapScan: 96, mapDump: "", mapFunctionScan: 1024,
      mapFunctionDump: "", mapFunctionContexts: @[], monitorEvery: 1,
      checkMap: false, startWidth: -1.0))
    runConfig = gp.toRunConfig
  runConfig.validateRunConfig
  installStandardParams()
  echoParams()
  processHelpParam()
  if gp.lat.len != 2 or gp.lat[0] < 4 or gp.lat[1] < 4 or
      (gp.lat[0] and 1) != 0 or (gp.lat[1] and 1) != 0:
    raiseValueError("change-of-variables HMC needs an even two-dimensional lattice")

  let
    lo = gp.lat.newLayout
    grt = initGraphRuntime()
    gc = actWilson(scalar.toGvalue(grt, gp.beta))
    spec = buildMapSpec(mp, gp.beta)
  withRng(gp.rng, R):
    var
      random = lo.newRNGField(R, gp.seed)
      acceptRandom: R
      g = lo.newgauge
      p = lo.newgauge
      loaded: type(g)
    acceptRandom.seed(gp.seed, 987654321)
    let layout = buildMapLayout(g[0], mp, spec)

    if myRank == 0:
      let s = mapScan(spec, mp.mapScan, mp.mapDump)
      echoMapSummary(spec, s, mp.mapScan)
      if mp.mapFunctionDump.len > 0:
        dumpMapFunction(spec, mp.mapFunctionScan, mp.mapFunctionContexts,
          mp.mapFunctionDump)
        echo "map functions: grid=", mp.mapFunctionScan,
          " intervals path=", mp.mapFunctionDump
      echo "lattice map: circle masks=", layout.circleMasks.len,
        " pair layers=", layout.pairLayers.len, " flowDepth=", spec.flowDepth
    commsBarrier()

    if runConfig.gaugefile.len == 0:
      if mp.startWidth >= 0.0: g.warm(mp.startWidth, random)
      else: g.random random
    else:
      g.loadOrInitGauge runConfig.gaugefile
      loaded = lo.newgauge
      threads:
        for mu in 0..<g.len: loaded[mu] := g[mu]
      invertMapFlow(g, spec, layout)

    let
      action = mapAction(gc, spec, layout)
      graph = buildTrajectoryGraph(grt, g, p, action, runConfig,
        buildTraining = false)

    var
      prevQ = 0.0
      proposalQ = 0.0
      proposalQValid = false
    block:
      let u0 = mapHost(graph.initialState.gauge.gaugeSnapshot, spec, layout).u
      let m = u0.topoMaxP2DU1
      echo "Initial physical plaq: ", u0.plaq3, "  topo: ", m.topo
      prevQ = m.topo
      if loaded.len > 0:
        echo "load round-trip maxdiff2: ", maxGaugeDiff2(u0, loaded)

    var
      Hvals = newSeq[float](runConfig.trajs)
      Jvals = newSeq[float](runConfig.trajs)
      Avals = newSeq[float](runConfig.trajs)
      Pvals = newSeq[float](runConfig.trajs)
      Qvals = newSeq[float](runConfig.trajs)
      dQchanged = newSeq[bool](runConfig.trajs)
      mdvals = newSeq[MdForceStats](runConfig.trajs)

    proc proposalMon(traj: int; dH, acc: float) =
      proposalQValid = false
      let
        show = mp.monitorEvery > 0 and (traj-1) mod mp.monitorEvery == 0
        prod = traj > runConfig.trajsThermo
      if not (show or prod): return
      let
        u = mapHost(graph.finalState.gauge.gaugeSnapshot, spec, layout).u
        m = u.topoMaxP2DU1
        dq = int(round(m.topo-prevQ))
      proposalQ = m.topo
      proposalQValid = true
      if show: echo &"proposal: dQ {dq}  maxP {m.maxP:.4f}"
      if prod: dQchanged[traj-runConfig.trajsThermo-1] = dq != 0

    proc measureTraj(traj: int; dH, acc: float; accepted: bool;
                     forceStats: MdForceStats) =
      let
        sm = mapHost(graph.initialState.gauge.gaugeSnapshot, spec, layout)
        pl = sm.u.plaq3
        q =
          if not accepted: prevQ
          elif proposalQValid: proposalQ
          else: sm.u.topo2DU1
        prod = traj > runConfig.trajsThermo
      echo "plaq: ", pl.re, "  topo: ", q, "  lnDet: ", sm.lndet
      prevQ = q
      proposalQValid = false
      if prod:
        let i = traj-runConfig.trajsThermo-1
        Hvals[i] = dH
        Jvals[i] = sm.lndet
        Avals[i] = acc
        Pvals[i] = pl.re
        Qvals[i] = q
        mdvals[i] = forceStats
        sm.u.maybeSaveGauge(runConfig, traj)

    if mp.checkMap:
      let
        v0 = graph.initialState.gauge.gaugeSnapshot
        sm = mapHost(v0, spec, layout)
      var back = lo.newgauge
      threads:
        for mu in 0..<back.len: back[mu] := sm.u[mu]
      invertMapFlow(back, spec, layout)
      let again = mapHost(back, spec, layout).u
      echo "checkMap inverse maxdiff2: ", maxGaugeDiff2(back, v0)
      echo "checkMap forward maxdiff2: ", maxGaugeDiff2(again, sm.u)

    runHmc(graph, runConfig, random, acceptRandom, measureTraj, proposalMon)
    if Hvals.len > 0:
      obstat(Hvals, Avals, Pvals, Qvals, gp.beta, lo.physVol,
        runConfig.trajs, gp.jkBlockSize, Jvals, mdvals)
      statsByQ(Hvals, Avals, dQchanged, gp.jkBlockSize, Jvals, mdvals)

  if gp.showTimers: echoTimers()
  if gp.showRunStats: grt.echoRunStats
  processSaveParams()
  writeParamFile()
  qexFinalize()
