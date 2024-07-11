import ../abstractFieldTheory

qexInit()

var
  fieldTheoryParams = %* {
    "lattice-geometry": @[4,4,4,4],
    "rank-geometry": @[1,1,1,1],
    "monte-carlo-algorithm": "hamiltonian-monte-carlo",
    "trajectory-length": 1.0,
    "serial-random-number-generator": "MILC",
    "parallel-random-number-generator": "MILC",
    "serial-random-number-seed": 987654321,
    "parallel-random-number-seed": 987654321,
    "start": "cold"
  }
  gaugeActionParams = %* {
    "action": "Wilson",
    "beta": 6.0,
    "steps": 10,
    "integrator": "MN2"
  }
  matterActionParams = %* {
    "smearing": "nHYP",
    "smearing-coefficients": @[0.4, 0.5, 0.5],
    "boundary-conditions": "aaaa"
  }
  fermionParams = %* {
    "mass": 0.0,
    "steps": 10,
    "integrator": "MN2",
  }
  bosonParams = %* {
    "mass": 0.75,
    "steps": 5,
    "integrator": "MN2"
  }
  ntrajs = 500
  ns = 2 # Number of staggered fermions (nf = 4*ns)
  nspv = 8 # Number of staggered Pauli-Villars bosons

# Create abstract field theory
var fieldTheory = newFieldTheory(fieldTheoryParams)

# Set gauge action
fieldTheory.setGaugeAction(gaugeActionParams)

# Add matter action w/ Nf = 8 (2 stag. fields) & 8 Pauli-Villars per fermion
fieldTheory.addMatterAction(matterActionParams, name = "PauliVillarsAction")
fieldTheory.addStaggeredFermion(fermionParams, ns): "PauliVillarsAction"
fieldTheory.addStaggeredBoson(bosonParams, nspv): "PauliVillarsAction"

# Compile action
fieldTheory.compile
  
# Run Hamiltonian Monte Carlo
for traj in 0..<ntrajs:
  fieldTheory.hamiltonian(heatbath = true)
  fieldTheory.evolve
  fieldTheory.hamiltonian
  let
    acc = fieldTheory.metropolis
    dH = fieldTheory.hf["hamiltonian"] - fieldTheory.hi["hamiltonian"]
  echo traj, " ", acc, " ", dH
  #if step == int(nsteps/2):
  #  fieldTheory.write("test")
  #  fieldTheory.read("test")

qexFinalize()