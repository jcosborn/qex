import qex
import mcmc

qexInit()

# HMC parameters
var
  nTrajectories = 10
  trajectoryLength = 1.0

# Specify field theory information
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
    "steps": 20,
    "integrator": "4MN5FP",
  }
  bosonInfo = %* { # Pauli-Villars info/parameters
    "mass": 0.75,
    "steps": 15,
    "integrator": "2MN",
  }

# Create a crazy action with many levels of nested 
# and parallel integration
var hmc = newLatticeFieldTheory(fieldTheoryInfo):
  fieldTheory.addGaugeAction(gaugeInfo)
  fieldTheory.addMatterAction(actionInfo):
    action.addStaggeredFermion(fermionInfo):
      subAction.addStaggeredFermion(fermionInfo)
      subAction.addStaggeredFermion(fermionInfo):
        subAction.addStaggeredBoson(bosonInfo):
          subAction.addStaggeredBoson(bosonInfo)
          subAction.addStaggeredBoson(bosonInfo)
          subAction.addStaggeredBoson(bosonInfo)
      subAction.addStaggeredBoson(bosonInfo):
        subAction.addStaggeredBoson(bosonInfo)
        subAction.addStaggeredBoson(bosonInfo)
        subAction.addStaggeredBoson(bosonInfo)
    action.addStaggeredFermion(fermionInfo)
    action.addStaggeredBoson(bosonInfo):
      subAction.addStaggeredBoson(bosonInfo)
      subAction.addStaggeredFermion(fermionInfo)
      subAction.addStaggeredFermion(fermionInfo)
      subAction.addStaggeredFermion(fermionInfo)

# Run HMC simulation
hmc.runHamiltonianMonteCarlo(nTrajectories,trajectoryLength):
  echo sample, " ", accepted, " ", hf-hi, " ", hi, " ", hf

qexFinalize()