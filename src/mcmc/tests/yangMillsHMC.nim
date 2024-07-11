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
    "parallel-random-number-seed": 987654321
  }
  gaugeActionParams = %* {
    "action": "Wilson",
    "beta": 6.0,
    "steps": 10,
    "integrator": "MN2"
  }
  ntrajs = 500

var fieldTheory = newFieldTheory(fieldTheoryParams)
fieldTheory.setGaugeAction(gaugeActionParams)
fieldTheory.compile
  
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