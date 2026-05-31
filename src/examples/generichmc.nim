import qex, macros, strutils, sequtils
#setGlobal("StaggeredSmearing","HISQ")
when getGlobal("DefaultLat","") != "":
  let defaultLat = split(getGlobal("DefaultLat",""),',').map(parseInt)

import physics/[stagD]
import physics/[stagSolve]

import hmc/[hmcAction]

qexInit(verb=2)

# specify parameters; overrwide with -<option>:<value> on command line
letParam:
  gaugeFilename = "checkpoint.lat"

  lattice = when declared(defaultLat): defaultLat else: @[8, 8, 8, 8]

  beta = 7.5
  mass = 0.005
  massH = 0.6
  massPV = 0.75

  numSt = 1 # number of staggered fermions
  numPV = 8 # number of staggered Pauli-Villars

  start = "cold"

  parallelSeed = 987654321
  serialSeed = 987654321

  outerIntegrator = "2MN"
  innerIntegrator = "2MN"
  outerSteps = 10
  innerSteps = 10

  startTraj = 0
  numTraj = 1
  trajectoryLength = 1.0

  actionMaxIter = 10000
  forceMaxIter = 10000
  actionTol = 1e-20
  forceTol = 1e-13

installStandardParams()
echoParams()
#echo "rank ", myRank, "/", nRanks
#threads: echo "thread ", threadNum, "/", numThreads
processHelpParam()

# set up lattice layout
let lo = lattice.newLayout()

# set up random number generator
type RngType = RngMilc6
var r = lo.newRNGField(RngType, parallelSeed.uint64)
var s: RngType
s.seed(serialSeed, 987654321)

# set up gauge field and gauge configuration
var
  g = lo.newGauge()
  uc = newGaugeConfiguration(g)

# initialize gauge configuration
if start == "cold": uc.cold()
elif start == "hot": uc.u.random(r)
elif start == "read": uc.read(gaugeFilename)
else: qexError "invalid start type"

#[ build staggered Dirac operator ]#

#when defined(HypSmearing):
when StaggeredSmearing == "HYP":
  let stag = newStag(uc.su)
#elif defined(StoutSmearing):
#  qexError "Stout smearing for HMC not yet implemented"
#elif defined(HisqSmearing):
elif StaggeredSmearing == "HISQ":
  let stag = newStag3(uc.su, uc.sul)
else:
  let stag = newStag(uc.u)

var spa = initSolverParams()
var spf = initSolverParams()

spa.r2req = actionTol
spf.r2req = forceTol

spa.maxits = actionMaxIter
spf.maxits = forceMaxIter

#spa.verbosity = 0
#spf.verbosity = 0

#[ build action in coordination with integrator levels ]#

var hmc* = uc.newHmcAction(s, r, trajectoryLength)

var gc = GaugeActionCoeffs(plaq: beta)
#var ga = gc.newGaugeAction(uc)
var ga = hmc.newGaugeAction(gc, uc)
var gaugeLevel = newActionLevel(multiplier = innerSteps, integrator = innerIntegrator)
gaugeLevel.add ga
hmc.add gaugeLevel   # inner level

var fermionLevel = newActionLevel(multiplier = outerSteps, integrator = outerIntegrator)

for i in 0..<numSt:
  fermionLevel.add newStaggeredFermionAction(stag, massH, spa, spf, r)
  fermionLevel.add newStaggeredRatioAction(stag, stag, mass, massH, spa, spf, r)

for i in 0..<numPV:
  fermionLevel.add newStaggeredPauliVillarsAction(stag, massPV, spa, spf, r)

hmc.add fermionLevel # outer level

#[ do HMC trajectory ]#

echo "== Action ===================="
echo hmc.description
echo "=============================="

hmc.init()
hmc.verbosity = 1
for traj in startTraj..<startTraj + numTraj: hmc.run()

#when declared(genericHmcTest):
#  genericHmcTest(hmc)

proc finalize* =
  processSaveParams()
  writeParamFile()
  qexFinalize()

when isMainModule:
  finalize()
