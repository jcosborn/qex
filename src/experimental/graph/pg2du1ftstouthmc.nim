## 2D U(1) subset-stout HMC, measuring and saving U = f(V).
import base/globals
setDefaultNc(1)
setVLENmax(4)

import qex
import core, scalar, gauge
import gauge/types as graphGauge
import hmcgauge/[config, trajectory, params, rng, gauge_io, ftstout, measure2du1]

qexInit()
echo "rank ", myRank, "/", nRanks
threads:
  echo "thread ", threadNum, "/", numThreads

let
  gp = readGaugeInputs(GaugeParams(
    lat: @[64, 64], beta: 6.5, seed: 1234567891'u, rng: rkPhilox4x64,
    savefile: "config", dt: 0.1, gsteps: 10, intalg: "2MN",
    trajsThermo: 64, trajs: 128, revCheckFreq: 0))
  sp = readStoutParams(StoutParams(rho: 0.1, nsmear: 1))
  runConfig = gp.toRunConfig
installStandardParams()
echoParams()
processHelpParam()
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
    g = lo.newGauge
    p = lo.newGauge
    loaded: graphGauge.Gauge
  if runConfig.gaugefile.len == 0:
    g.random random
  else:
    g.loadOrInitGauge runConfig.gaugefile
    loaded = lo.newGauge
    threads:
      for mu in 0..<g.len: loaded[mu] := g[mu]
    let inv = invertStoutFlow(g, sp.rho, sp.nsmear)
    echo "stout inverse U->V: iter ", inv.iter, "  rdf2 ", inv.rdf2
  let graph = buildTrajectoryGraph(grt, g, p, sa.action, runConfig, buildTraining = false)
  runFlowHmc(graph, sa.flow, runConfig, random, acceptRandom, gp.beta, gp.jkBlockSize, loaded = loaded)

if gp.showTimers: echoTimers()
if gp.showRunStats: grt.echoRunStats
processSaveParams()
writeParamFile()
qexFinalize()
