## Learned 2D U(1) HMC with cold or imported latent initialization.
import base/globals
setDefaultNc(1)
setVLENmax(4)

import qex
import ../graph/[core, scalar, gauge]
import ../graph/hmcgauge/[config, trajectory, params, rng, measure2du1]
import model, graph, io

qexInit()
echo "rank ", myRank, "/", nRanks
threads:
  echo "thread ", threadNum, "/", numThreads

let
  checkpoint = strParam("checkpoint", "")
  latent = strParam("latent", "")
  precision = strParam("precision", "single")
  gp = readGaugeInputs(GaugeParams(
    lat: @[8, 8], beta: 3.0, seed: 1029'u, rng: rkPhilox4x64,
    savefile: "config", dt: 0.35, gsteps: 10, intalg: "2MNp",
    trajsThermo: 0, trajs: 1, revCheckFreq: 0))
  runConfig = gp.toRunConfig
installStandardParams()
echoParams()
processHelpParam()
runConfig.validateRunConfig
if checkpoint.len == 0:
  raiseValueError("learned flow requires -checkpoint:<manifest directory>")
if precision notin ["single", "double"]:
  raiseValueError("precision must be single or double")
if gp.lat.len != 2:
  raiseValueError("learned HMC requires a 2D lattice")
for n in gp.lat:
  if n < 4 or n mod 2 != 0:
    raiseValueError("learned HMC requires even lattice extents of at least four")
if gp.gaugefile.len > 0:
  raiseValueError("learned HMC requires latent initialization; select -latent:<array name> in the checkpoint directory")

let
  lo = gp.lat.newLayout
  grt = initGraphRuntime()
  gc = actWilson(scalar.toGvalue(grt, gp.beta))
  sa =
    if precision == "double": learnedAction(gc, toNnftModel(grt, loadNnft(checkpoint, float64)))
    else: learnedAction(gc, toNnftModel(grt, loadNnft(checkpoint, float32)))

withRng(gp.rng, R):
  var
    random = lo.newRNGField(R, gp.seed)
    acceptRandom: R
  acceptRandom.seed(gp.seed, 987654321)
  let g = lo.newGauge
  let p = lo.newGauge
  threads:
    for mu in 0..<g.len: g[mu] := 1
  if latent.len > 0: loadLatent(g, checkpoint, latent)
  echo "Learned checkpoint: ", checkpoint, "  precision: ", precision,
    "  latent: ", (if latent.len > 0: latent else: "cold")
  let graph = buildTrajectoryGraph(grt, g, p, sa.action, runConfig, buildTraining = false)
  runFlowHmc(graph, sa.flow, runConfig, random, acceptRandom, gp.beta, gp.jkBlockSize)

if gp.showTimers: echoTimers()
if gp.showRunStats: grt.echoRunStats
processSaveParams()
writeParamFile()
qexFinalize()
