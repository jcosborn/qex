import qex
import mcmc/mcmc

import strutils

qexInit()

let 
  # Get command line information
  cmdInfo = readCMD()
  ntraj = cmdInfo["ntraj"].getInt()
  saveInterval = cmdInfo["saveInterval"].getInt()
  ioPath = cmdInfo["ioPath"].getStr()
  fn = cmdInfo["fn"].getStr()
  jsonPath = cmdInfo["jsonPath"].getStr()
  flowString = cmdInfo["flow"].getStr()
  flow = parseBool(flowString)
  flowInterval = case flow
    of true: cmdInfo["flowInterval"].getInt()
    of false: 0
  
  # Get information from JSON file
  jsonInfo = readJSON(jsonPath)
  tau = jsonInfo["hmc"]["trajectory-length"].getFloat()
  nspv = jsonInfo["staggered-pauli-villars"]["species"].getInt()
  nsf = jsonInfo["staggered-fermions"]["species"].getInt()
  nrsf = jsonInfo["rooted-staggered-fermions"]["species"].getInt()
var 
  cfg = cmdInfo["cfg"].getInt()
  flowInfo = case flow 
    of true: jsonInfo["flow"]
    of false: parseJson("{}")

# Construct lattice field theory
var 
  hmc = newLatticeFieldTheory(jsonInfo["hmc"]):
    fieldTheory.addGaugeMatterAction(jsonInfo["action"]):
      action.addGaugeField(jsonInfo["gauge"])
      for _ in 0..<nspv: 
        action.addStaggeredBoson(jsonInfo["staggered-pauli-villars"])
      for _ in 0..<nsf: 
        action.addStaggeredFermion(jsonInfo["staggered-fermions"])
      for _ in 0..<nrsf: 
        action.addRootedStaggeredFermion(jsonInfo["rooted-staggered-fermions"])

# Read in gauge field
if cfg != 0: hmc.read(ioPath & fn & "_" & $(cfg), onlyGauge = false)

# Run HMC
hmc.runHamiltonianMonteCarlo(ntraj,tau):
  # Print out information about this trajectory
  var acc = case accepted 
    of true: "ACC"
    of false: "REJ"     
  echo trajectory, " ", acc, " dH ", hf-hi, " ", hi, " ", hf

  # Simple measurements
  u.plaquette
  u.polyakov

  # Gradient flow
  if flow:
    if (trajectory+1) mod flowInterval == 0:
      for flowName in flowInfo.keys(): 
        let filename = flowName & "_" & $(trajectory+1) & ".log"
        flowInfo[flowName].setFilename(filename)
      u.gradientFlow(flowInfo):
        let output = measurements.formatMeasurements(style="KS_nHYP_FA")
        f.write(output & "\n")

  # Save gauge field
  if (trajectory+1) mod saveInterval == 0: 
    hmc.write(ioPath & fn & "_" & $(cfg+trajectory+1), onlyGauge = false)

qexFinalize()
