import qex
import mcmc

qexInit()

# HMC parameters
var
  nTrajectories = 10
  trajectoryLength = 1.0

# JSON data types for specifying simulation information.
# JSON data types in Nim are basically Nim's version
# of Python's "dictionary" data type
var 
  fieldTheoryInfo = %* { # Field theory info/parameters
    "lattice-geometry": @[4,4,4,4],
    "mpi-geometry": @[1,1,1,1],
    "monte-carlo-algorithm": "hmc",
    "serial-random-number-seed": 987654321,
    "parallel-random-number-seed": 987654321,
    "serial-random-number-generator": "milc",
    "parallel-random-number-generator": "milc",
    "start": "cold"
  }
  actionInfo = %* { # Action info/parameters
    "smearing": "nhyp",
    "smearing-coefficients": @[0.4,0.5,0.5],
    "boundary-conditions": "aaaa"
  }
  gaugeFieldInfo = %* { # Gauge field info/parameters
    "group": "su",
    "action": "Wilson",
    "beta": 9.0,
    "steps": 10,
    "integrator": "2MN"
  }
  fermionInfo = %* { # Fermion field info/parameters
    "mass": 0.0,
    "steps": 10,
    "integrator": "2MN",
  }
  rootedFermionInfo = %* { # Rooted fermion info/parameters
    "mass": 0.0,
    "steps": 10,
    "integrator": "2MN",
    "nf": 2 # flavor #; can be 1-4
  }
  pauliVillarsInfo = %* { # Pauli-Villars info/parameters
    "mass": 0.75,
    "steps": 10,
    "integrator": "2MN",
  }

# Specify number of unrooted staggered fermions & Pauli-Villars fields
var 
  nStag = 2 # Staggered fermions (not rooted)
  nPV = 12 # Pauli-Villars

# Construct full lattice field theory
var hmc = newLatticeFieldTheory(fieldTheoryInfo):
  fieldTheory.addGaugeMatterAction(actionInfo):
    action.addGaugeField(gaugeFieldInfo) # Add gauge field
    for _ in 0..<nStag: # Add "nStag=2" (not rooted) stag. species
        action.addStaggeredFermion(fermionInfo) 
    action.addRootedStaggeredFermion(rootedFermionInfo) # Add 2 more flavors
    for _ in 0..<nPV: # Add "nPV" Pauli-Villars
        action.addStaggeredBoson(pauliVillarsInfo) 

# Run HMC simulation!
hmc.runHamiltonianMonteCarlo(nTrajectories,trajectoryLength):
  #[ 
    This code block will be executed after each trajectory. 
    This allows you to do whatever you want w/ the gauge field
    after each trajectory & print out important information.
  ]#

  # For example, you can print out information about this trajectory
  echo sample, " ", accepted, " ", hf-hi, " ", hi, " ", hf

  # Or you could call a function that measures some observable, for example:
  # plaquette(u)

  # Or you can save the gauge field. It's up to you!
  hmc.write("<location>",onlyGauge=false)

qexFinalize()
