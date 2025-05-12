import qex
import mcmc/mcmc
import ../actions/latticeActionUtils

qexInit()

#[ Simple reversibility check ]#

# Get information from JSON file
let 
  jsonInfo = %* {
    "hmc": {
        "lattice-geometry": [4,4,4,4],
        "mpi-geometry": [1,1,1,1],
        "monte-carlo-algorithm": "hmc",
        "trajectory-length": 1.0,
        "serial-random-number-seed": 987654321,
        "parallel-random-number-seed": 987654321,
        "serial-random-number-generator": "milc",
        "parallel-random-number-generator": "milc",
        "start": "cold"
    },
    "action": {
        "smearing": "nhyp",
        "smearing-coefficients": [0.4,0.5,0.5],
        "boundary-conditions": "aaaa"
    },
    "gauge": {
        "group": "su",
        "action": "Wilson",
        "beta": 9.0,
        "steps": 10,
        "integrator": "2MN"
    },
    "pv": {
        "mass": 0.75,
        "steps": 10,
        "integrator": "2MN",
        "number-staggered-pv-fields": 8
    }
  }
  tau = jsonInfo["hmc"]["trajectory-length"].getFloat()
  nPV = jsonInfo["pv"]["number-staggered-pv-fields"].getInt()

# Construct lattice field theory
var 
  hmc = newLatticeFieldTheory(jsonInfo["hmc"]):
    fieldTheory.addGaugeMatterAction(jsonInfo["action"]):
      action.addGaugeField(jsonInfo["gauge"])
      for _ in 0..<nPV: action.addStaggeredBoson(jsonInfo["pv"])

# Run HMC
hmc.runHamiltonianMonteCarlo(1,tau):
  echo hf-hi, " ", hi, " ", hf
  hi = hf
  runMolecularDynamics(-tau)
  for action in hmc.actions: action.smearGauge
  hf = hamiltonian()
  echo hf-hi, " ", hi, " ", hf

#[ Reversibility check on more complicated MD ]#

qexFinalize()