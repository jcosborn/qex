import qex

import physics/[stagD]
import physics/[stagSolve]

import hmc/[hmcAction]

qexInit()

# specify parameters; overrwide with -<option>:<value> on command line
letParam:
  gaugeFilename = "checkpoint.lat"
  
  lattice = @[8, 8, 8, 8]

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
echo "rank ", myRank, "/", nRanks
threads: echo "thread ", threadNum, "/", numThreads
processHelpParam()

# set up lattice layout
let lo = lattice.newLayout()

# set up random number generator
var r = lo.newRNGField(RngMilc6, parallelSeed.uint64)
var s: RngMilc6
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

when defined(HypSmearing):
  let stag = newStag(uc.su)
elif defined(StoutSmearing):
  qexError "Stout smearing for HMC not yet implemented"
elif defined(HisqSmearing): 
  let stag = newStag3(uc.su, uc.sul)
else: 
  let stag = newStag(uc.u)

var spa = initSolverParams()
var spf = initSolverParams()

spa.r2req = actionTol
spf.r2req = forceTol

spa.maxits = actionMaxIter
spf.maxits = forceMaxIter

spa.verbosity = 1
spf.verbosity = 1

#[ build action in coordination with integrator levels ]#

var gc = GaugeActionCoeffs(plaq: beta)
var ga = gc.newGaugeAction()

var fermionLevel = newActionLevel(multiplier = outerSteps, integrator = outerIntegrator)

for i in 0..<numSt:
  fermionLevel.add newStaggeredFermionAction(stag, massH, spa, spf)
  fermionLevel.add newStaggeredRatioAction(stag, stag, mass, massH, spa, spf)

for i in 0..<numPV:
  fermionLevel.add newStaggeredPauliVillarsAction(stag, massPV, spa, spf)

var gaugeLevel = newActionLevel(multiplier = innerSteps, integrator = innerIntegrator)
gaugeLevel.add ga

var hmc = uc.newHmcAction(s, r, trajectoryLength)
hmc.add gaugeLevel   # inner level
hmc.add fermionLevel # outer level

#[ do HMC trajectory ]#

hmc.init()
hmc.verbosity = 1
for traj in startTraj..<startTraj + numTraj: hmc.update()

qexFinalize()