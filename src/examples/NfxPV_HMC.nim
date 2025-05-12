import qex
import mcmc/mcmc

qexInit()

let 
  # Get command line information
  cmdInfo = readCMD()
  ntraj = cmdInfo["ntraj"].getInt()
  saveInterval = cmdInfo["saveInterval"].getInt()
  ioPath = cmdInfo["ioPath"].getStr()
  fn = cmdInfo["fn"].getStr()
  jsonPath = cmdInfo["jsonPath"].getStr()

  # Get information from JSON file
  jsonInfo = readJSON(jsonPath)
  tau = jsonInfo["hmc"]["trajectory-length"].getFloat()
  nPV = jsonInfo["pv"]["number-staggered-pv-fields"].getInt()
  nf = jsonInfo["Nf"]["number-staggered-Nf-fields"].getInt()
  nfr = jsonInfo["Nfr"]["number-staggered-Nfr-fields"].getInt()
var cfg = cmdInfo["cfg"].getInt()

# Construct lattice field theory
var 
  hmc = newLatticeFieldTheory(jsonInfo["hmc"]):
    fieldTheory.addGaugeMatterAction(jsonInfo["action"]):
      action.addGaugeField(jsonInfo["gauge"])
      for _ in 0..<nPV: action.addStaggeredBoson(jsonInfo["pv"])
      for _ in 0..<nf: action.addStaggeredFermion(jsonInfo["Nf"])
      for _ in 0..<nfr: action.addRootedStaggeredFermion(jsonInfo["Nfr"])

# Read in gauge field
if cfg != 0: hmc.read(ioPath & fn & "_" & $(cfg), onlyGauge = false)

# Run HMC
hmc.runHamiltonianMonteCarlo(ntraj,tau):
  # Print out information about this trajectory
  var acc=case accepted 
    of true: "ACC"
    of false: "REJ"     
  echo trajectory, " ", acc, " dH ", hf-hi, " ", hi, " ", hf

  # Simple measurements
  u.plaquette
  u.polyakov

  # Save gauge field
  if (trajectory+1) mod saveInterval == 0: 
    hmc.write(ioPath & fn & "_" & $(cfg+trajectory+1), onlyGauge = false)
    #inc cfg

qexFinalize()
