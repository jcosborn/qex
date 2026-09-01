## Production script enabling HMC for eight flavor symmetric mass generation 
## project by the Lattice Strong Dynamics collaboration.
## 
## Designed to be backwards-compatible with stag_pv_hmc code under staghmc-devel
## branch of QEX. Hence, some of the code is reflective of conventions that I 
## would not normally use these days (such as the inputs not using camel case
## and some of the strange choices for the formatting of the input XML file). 
## 
## If you're using this as an example for your own project, please make better
## choices for your own naming conventions -- don't take them from here. 
## 
## Author: Curtis Taylor Peterson <curtistaylorpetersonwork@gmail.com>

import std/[os]
import std/[strutils]
import std/[sequtils]

import qex

setGlobal("StaggeredSmearing", "HYP") # set HYP smearing
setGlobal("FUELCompat", "1") # enable FUEL compatibility (RNG)

import physics/[stagD]
import physics/[stagSolve]
import hmc/[hmcAction]

qexInit(verb = 2)

echo "rank ", myRank, "/", nRanks
threads: echo "thread ", threadNum, "/", numThreads

#[ baseline setup ]#

# read command line inputs
letParam:
  start_config = 0 # starting configuration
  end_config   = 1 # ending configuration
  config_space = 1 # number of trajectories (units) between configurations
  save_freq    = 0 # frequency of saving configurations (in configuration units)
  
  path = "./"          # global path to IO directory
  filename = "ckpoint" # base filename for checkpointing 
  xml = currentSourcePath.parentDir() & "/eightFlavorSMG.xml" # path to input XML
  verbosity = 0 # verbosity of information printed by HMC engine

  rank_geom = newSeq[int](0) # <-+- if empty, guesses layout
  simd_geom = newSeq[int](0) # <-+

# read xml inputs
letXml xml:
  hmc:
    tau = 1.0    # trajectory length 
    f_steps = 10 # number of fermion (outer) fermion/PV integrator steps
    g_steps = 3  # number of gauge (inner) integrator steps *PER* outer gauge update
    ferm_int_alg =  "2MN" # integrator for fermion/PV level
    gauge_int_alg = "2MN" # integrator for gauge level
    no_metropolis_until = 0 # number of trajectories to run w/o Metropolis test

  config_opts:
    start = "cold" # starting configuration: cold, hot, or read

  rng:
    parallel_seed = 987654321 # seed for the parallel RNG
    serial_seed   = 987654321 # seed for the serial (Metropolis) RNG

  action:
    geom:
      Ns = 8  # spatial extent
      Nt = 16 # temporal extent
      bc = "aaaa" # boundary conditions (a = antiperiodic, p = periodic)

    gauge:
      beta = 6.0 # bare gauge coupling
    
    ferm:
      Nf = 2      # number of staggered species (e.g., Nf = 1 for 4 Dirac fermions)
      mass = 0.05 # bare staggered fermion mass
      Nh = 1       # number of Hasenbusch per staggered species
      mass_h = 0.3 # bare Hasenbusch mass
    
    pv:
      num_pv = 8     # overall number of staggered Pauli-Villars species
      mass_pv = 0.75 # bare staggered Pauli-Villars mass

  smearing:
    nhyp_smearing: # nHYP smearing parameters
      alpha_1 = 0.4
      alpha_2 = 0.5
      alpha_3 = 0.5

  solver:
    a_tol = 1e-24    # action solver tolerance
    a_maxits = 10000 # action solver max iterations
    f_tol = 1e-16    # force solver tolerance
    f_maxits = 10000 # force solver max iterations
  
  basic_meas:
    plaq:
      plaq_freq = 1 # frequency of plaquette measurement (trajectory units)
    
    ploop_freq = 1  # frequency of Polyakov loop measurement (trajectory units)
    s4_freq = 1     # frequency of "s4 order parameter" measurement (trajectory units)
    
    hmc_checks:
      rev_check_freq = 0 # frequency of reversibility checks (trajectory units)

installStandardParams()
processHelpParam()

echoParams()
echoXml()

if Nh != 0 and Nh != 1: qexError "eight flavor project only uses one Hasenbusch field"

# set lattice up
let lattice = @[Ns, Ns, Ns, Nt]
let lo = case rank_geom.len:
  of 0: 
    case simd_geom.len:
      of 0: lattice.newLayout()
      else: 
        let comm = getDefaultComm()
        comm.newLayoutX(lattice, VLEN, @[], simd_geom)
  else: 
    case simd_geom.len:
      of 0: lattice.newLayout(rank_geom)
      else: lattice.newLayout(VLEN, rank_geom, simd_geom)

# set RNG up
type RngType = RngMilc6
var r = lo.newRNGField(RngType, parallel_seed.uint64)
var s: RngType
s.seed(serial_seed, 987654321)

# set gauge field up
var
  g = lo.newGauge()
  uc = newGaugeConfiguration(g, alpha1 = alpha_1, alpha2 = alpha_2, alpha3 = alpha_3)

# set up staggered Dirac operator & solver
let stag = newStag(uc.su)
var spa = initSolverParams()
spa.r2req = a_tol
spa.maxits = a_maxits
var spf = initSolverParams()
spf.r2req = f_tol
spf.maxits = f_maxits

#[ build action/integrator from inputs ]#

var hmc = uc.newHmcAction(s, r, tau, revCheckFreq = rev_check_freq)

var fermionLevel = newActionLevel(multiplier = f_steps, integrator = ferm_int_alg)

for _ in 0..<Nf:
  if Nh == 0: fermionLevel.add newStaggeredFermionAction(stag, mass, spa, spf, r, bc = bc)
  else:
    fermionLevel.add newStaggeredFermionAction(stag, mass_h, spa, spf, r, bc = bc)
    fermionLevel.add newStaggeredRatioAction(stag, stag, mass, mass_h, spa, spf, r, bc = bc)

for _ in 0..<num_pv:
  fermionLevel.add newStaggeredPauliVillarsAction(stag, mass_pv, spa, spf, r, bc = bc)

hmc.add fermionLevel # outermost level

var gc = GaugeActionCoeffs(plaq: beta, adjplaq: -0.25*beta)
var ga = hmc.newGaugeAction(gc, uc)
var gaugeLevel = newActionLevel(multiplier = g_steps, integrator = gauge_int_alg)

gaugeLevel.add ga
hmc.add gaugeLevel # innermost level

#[ do HMC ]#

echo "== Action ===================="
echo hmc.description
echo "=============================="

# initialize gauge configuration
if start == "cold": hmc.cold()
elif start == "hot": hmc.hot()
elif start == "read": 
  hmc.read(
    readParallelRNG = true,
    readSerialRNG = true,
    gaugeFilename = path & "/" & filename & "_" & $(start_config) & ".lat",
    parallelRNGFilename = path & "/" & filename & "_" & $(start_config) & ".parallelRNG",
    serialRNGFilename = path & "/" & filename & "_" & $(start_config) & ".serialRNG"
  )
else: qexError "invalid start type"

# run HMC trajectories
hmc.init(verbosity = verbosity)
for config in start_config..<end_config:
  for traj in 0..<config_space:
    # run HMC
    hmc.run(forceAccept = config < no_metropolis_until)

    echo "== measurements =========="
    if (plaq_freq > 0) and ((traj + 1) mod plaq_freq == 0): hmc.measurePlaquette()
    if (s4_freq > 0) and ((traj + 1) mod s4_freq == 0): hmc.measurePlaquetteS4()
    if (ploop_freq > 0) and ((traj + 1) mod ploop_freq == 0): hmc.measurePolyakovLoop()
    echo "=========================="

    # because Anna will most certainly ask about it at some point, if you want to
    # make a measurement of some observable at this stage of the HMC, just copy-and-
    # paste the procedure for it in this code and call it at this location. You can
    # access the gauge field of the HMC object with hmc.getGauge(). That is, you 
    # don't need to add a procedure to hmcAction --- you can just measure it here.
    # So if Anna asks, "can it measure psi-bar-psi?", you can say "yes, I just need
    # to copy it into the HMC driver code" ;)
    
  # save configuration
  if (save_freq > 0) and ((config + 1) mod save_freq == 0): 
    hmc.write(
      writeParallelRNG = true,
      writeSerialRNG = true,
      gaugeFilename = path & "/" & filename & "_" & $(config + 1) & ".lat",
      parallelRNGFilename = path & "/" & filename & "_" & $(config + 1) & ".parallelRNG",
      serialRNGFilename = path & "/" & filename & "_" & $(config + 1) & ".serialRNG"
    )

qexFinalize()