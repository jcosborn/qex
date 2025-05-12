import qex
import mcmc/mcmc

qexInit()

let 
  # Get command line information
  cmdInfo = readCMD()
  ntraj = cmdInfo["ntraj"].getInt()
  saveInterval = cmdInfo["saveInterval"].getInt()
  ioPath = cmdInfo["ioPath"].getStr()
  jsonPath = cmdInfo["jsonPath"].getStr() & "gaugeBosonInput.json"

  # Get information from JSON file
  jsonInfo = readJSON(jsonPath)
  tau = jsonInfo["hmc"]["trajectory-length"].getFloat()
  nPV = jsonInfo["pv"]["number-staggered-pv-fields"].getInt()
var cfg = cmdInfo["cfg"].getInt()

# Construct lattice field theory
var 
  hmc = newLatticeFieldTheory(jsonInfo["hmc"]):
    fieldTheory.addGaugeMatterAction(jsonInfo["action"]):
      action.addGaugeField(jsonInfo["gauge"])
      for _ in 0..<nPV: action.addStaggeredBoson(jsonInfo["pv"])

# Read in gauge field
if cfg != 0: hmc.read(ioPath & "checkpoint_" & $(cfg), onlyGauge = false)

# Run HMC
hmc.runHamiltonianMonteCarlo(ntraj,tau):
  # Print out information about this trajectory
  echo trajectory, " ", accepted, " ", hf-hi, " ", hi, " ", hf

  # Simple measurements
  u.plaquette
  u.polyakov

  # Save gauge field
  if (trajectory+1) mod saveInterval == 0: 
    hmc.write(ioPath & "checkpoint_" & $(cfg+1), onlyGauge = false)
    inc cfg

qexFinalize()