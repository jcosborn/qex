import qex
import examples/[hisqhmc_h]
import fermionflowutils
import json

import gauge/[gaugefix, hisqsmear]
import physics/[stagSolve]
import observables/[nlStagMeson]

import std/[os, times, sugar, tables] 
import std/[parseopt, strutils, strformat]
import std/[json, sequtils]
import std/[math]
import std/[re]

# import quda/quda

const
  DefaultCoulombGaugeFixEpsilon = 1e-8
  DefaultCoulombGaugeFixRelax = 1.75

const
  DefaultConjugateGradientTolerance = 1e-20
  DefaultConjugateGradientMaximumIterations = 10000
  DefaultConjugateGradientVerbosty = 1

const Corners = 8
const Masses = 2

# const
#   mass = [0.01, 0.02, 0.03]
proc extractMFloats(inputStr: string): seq[float] =
  ## Extracts all numbers following 'm' in a string and returns them as 0.xxxx floats.
  ## Example:
  ##   let floats = extractMFloats("l3248f211b580m002426m06730m8447a")
  ##   # Returns @[0.002426, 0.0673, 0.8447]
  result = @[]  # Initialize an empty sequence

  # Find all patterns where 'm' is followed by digits
  for m in inputStr.findAll(re"m(\d+)"):
    echo m[1 .. ^1]
    let decimalStr = "0." & m[1 .. ^1]  # Prepend "0." to the matched digits
    result.add(decimalStr.parseFloat)  # Convert to float and add to result

proc rephase[U](u: seq[U]) = 
  threads:
    u.setBC()
    threadBarrier()
    u.stagPhase()

const 
  logStyle = "KS_nHYP_FA"
  banner = """
|---------------------------------------------------------------|
 Quantum EXpressions (QEX)

 QEX authors: James Osborn & Xiao-Yong Jin
 QEX gradient flow authors: 
   - James Osborn (Argonne National Laboratory)
   - Curtis Taylor Peterson [C.T.P.] (Michigan State University)
 QEX fermion flow authors:
   - Mingwei Dai (University of Illinois Urbana-Champaign)
 QEX GitHub: https://github.com/jcosborn/qex
 Gauge flow GitHub: https://github.com/ctpeterson/qex
 C.T.P. email: curtistaylorpetersonwork@gmail.com
 cite: Proceedings of Science (PoS) LATTICE2016 (2017) 271
|---------------------------------------------------------------|
"""

qexInit()

var 
  #reads command line arguments into a json object
  cmd = readCMD() 
  #e.g. command: ./fermionflow --flow-json:flow-info.json --lattice-json:lattice-info.json --configuration:100 --base-filename:config_name_stem --output:fflow_corr
  flowInfo = case cmd.hasKey("flow-json")
    of true: readJSON(cmd["flow-json"].getStr())
    of false:
      qexError "json file for flow information not specified"
      parseJson("{}")
  latInfo = case cmd.hasKey("lattice-json")
    of true: readJSON(cmd["lattice-json"].getStr())
    of false:
      qexError "json file for lattice information not specified"
      parseJson("{}")
  #[
  #08/11/2025: contractInfo not used yet 
  contractInfo = case cmd.hasKey("contract-json") 
    of true: readJSON(cmd["contract-json"].getStr())
    of false:
      echo "json file for contract information not specified"
      parseJson("{}")
  ]#
  cfg = case cmd.hasKey("configuration")
    of true: $cmd["configuration"].getInt()
    of false:
      qexError "configuration number not specified"
      "0"
  # name of the configuration without the extension
  filename = case cmd.hasKey("base-filename")
    of true: cmd["base-filename"].getStr()
    of false: "checkpoint"
  latLayout = case latInfo.hasKey("lattice-geometry")
    of true: latInfo["lattice-geometry"].getIntSeq()
    of false: 
      qexError "must specify lattice-geometry in lattice input file"
      @[8,8,8,8]
  lo = case latInfo.hasKey("rank-geometry")
    of true:
      case latInfo.hasKey("simd-geometry"):
        of true:
          newLayout(
            latLayout,
            VLEN,
            latInfo["rank-geometry"].getIntSeq(),
            latInfo["simd-geometry"].getIntSeq()
          )
        of false: newLayout(latLayout,latInfo["rank-geometry"].getIntSeq())
    of false: newLayout(latLayout)
  u = lo.newGauge()
# Set mass from filename
let mass = extractMFloats(filename)[0..Masses]
#path to save the correlation functions
let corrFile = case cmd.hasKey("output")
  of true: cmd["output"].getStr()
  of false: "__NOOUTPUT__"
echo "corrFile: ", corrFile
#read configuration
echo "Start leading configuration. "
u.readGauge(filename & "." & cfg & ".ildg")
# u.random
echo "Loading complete. "
#initialize stagSolve
let
    # conjugate gradient parameters
  r2req = DefaultConjugateGradientTolerance
  maxits = DefaultConjugateGradientMaximumIterations
    # gauge fixing parameters
  fixGauge: bool = true
  gfstop = DefaultCoulombGaugeFixEpsilon
  gforf = DefaultCoulombGaugeFixRelax
    # nhyp smearing parameters
  smearGauge: bool = true
  # nhyp = HypCoefs(alpha1: 0.4, alpha2: 0.5, alpha3: 0.5)
  hisq = newHisq() 
    # source time
  srcT = 0
    # time extent
  nt = latInfo["lattice-geometry"].getIntSeq[^1]
threads:
  block:
    let d = u.checkSU
    echo "unitary deviation avg: ", d.avg, "  max: ", d.max
  threadBarrier()
  u.projectSU
  threadBarrier()
  block:
    let d = u.checkSU
    echo "new unitary deviation avg: ", d.avg, "  max: ", d.max
# u.echoPlaq
# Coulomb gauge fix
if fixGauge:
  tic "Coulomb gauge fixing"
  var tmat = lo.ColorMatrix()
  threads: tmat := 1
  getGaugeFixTransform(tmat, u, @[0,1,2], gfstop, gforf, verb=0)
  var fg = lo.newGauge
  fg.gaugeTransform(u, tmat)
  threads: u := fg
  qexLog "Coulomb gauge fixing done. Time: ", getElapsedTime()
# smear and rephase
u.rephase()
var su = lo.newGauge()
var sul = lo.newGauge()
if smearGauge:
  tic "HISQ gauge smearing"
  discard hisq.smearGetForce(u, su, sul)
  # threads: u := su ???????????
  qexLog "nHISQ gauge smearing done. Time: ", getElapsedTime()
# u.rephase()

let stag = newStag3(su,sul)
# instantiate solver parameters
var sp = initSolverParams()
sp.r2req = r2req
sp.maxits = maxits

#[ spectrum ]#
# build meson operators & select mesons as subset
let 
  mesonOps = buildMesonOps()
  mesons = [
    # mesonOps[0],  # γ_0 γ_5 ⊗ γ_0 γ_5 (P = +σ) I     π_05  <--+ 
      mesonOps[1],  # γ_5 ⊗ ξ_5         (P = -σ) II    π_5      |
      # mesonOps[6],  # γ_5 ⊗ ξ_μ ξ_5     (P = +σ) VII   π_i5     |
      mesonOps[7],  # γ_0 γ_5 ⊗ ξ_μ ξ_ν (P = -σ) VIII  π_ij     | scalar &
      # mesonOps[12], # γ_0 γ_5 ⊗ ξ_μ     (P = -σ) XIII  π_0i     | pseudoscalar
      # mesonOps[13], # γ_5 ⊗ ξ_μ         (P = -σ) XIV   π_i      |
      # mesonOps[16], # γ_0 γ_5 ⊗ I       (P = -σ) XVI   π_I      |
      # mesonOps[17], # γ_5 ⊗ ξ_0         (P = +σ) XVII  π_0   <--+
      # mesonOps[2],  # γ_0 γ_k ⊗ ξ_0 ξ_μ (P = +σ) XIII  ρ_0i  <--+
      # mesonOps[3],  # γ_μ ⊗ ξ_μ         (P = -σ) XIV   ρ_i      |
      # mesonOps[4],  # γ_μ ⊗ I           (P = +σ) XV    ρ_I      |
      # mesonOps[5],  # γ_0 γ_μ ⊗ ξ_0     (P = +σ) VI    ρ_0      |
      # mesonOps[8],  # γ_0 γ_μ ⊗ ξ_μ ξ_5 (P = +σ) XIV   ρ_i5     |
      # mesonOps[9],  # γ_μ ⊗ ξ_μ ξ_ν     (P = +σ) X     ρ_ij     | vector &
      # mesonOps[10], # γ_0 γ_μ ⊗ ξ_0 ξ_5 (P = -σ) XI    ρ_05     | pseudovector
      # mesonOps[11], # γ_0 γ_μ ⊗ ξ_5     (P = -σ) XII   ρ_5      |
      # mesonOps[14], # γ_μ ⊗ ξ_μ         (P = -σ) XIV   ρ_i*     |
      # mesonOps[15], # γ_0 γ_k ⊗ ξ_0 ξ_μ (P = -σ) XIII  ρ_0i*    |
      # mesonOps[18], # γ_0 γ_μ ⊗ ξ_μ ξ_5 (P = -σ) XVIII ρ_i5*    |
      # mesonOps[19], # γ_μ ⊗ ξ_μ ξ_ν     (P = +σ) XX    ρ_ij* <--+
    ]
  
  massname = [
    "ll", "ss", #"cc",# "ls", "lc", "sc", "sl"
  ]
  mesonname = [
    "g5g5", 
    "g0g5gigj",
  ]
# initialize propagators and correlators
const colors = u[0][0].nrows
var 
  props: array[colors, array[Corners, array[Masses, typeof(lo.ColorVector)]]]
  mprops = [[0,0]] #,[1,1],[2,2],[0,1],[0,2],[1,2],[1,0]
  corrs = newSeqWith(mprops.len, newSeqWith(mesons.len, newSeq[float]()))#array[mprops.len, newSeq[seq[float]](mesons.len)]
for color in 0..<colors:
  for massIdx in 0..<Masses:
    for corner in 0..<Corners: props[color][corner][massIdx] = lo.ColorVector()
for propindex in 0..<mprops.len:
  for pionIdx in 0..<mesons.len: corrs[propindex][pionIdx].newSeq(nt)
tic "HISQ Dirac inversion"
# get spectrum and do it at each color
for color in 0..<colors:
  for corner in 0..<Corners:
    # generate color cource
    var src = lo.ColorVector()
    threads:
      let cc = lo.coords
      src := 0
      threadBarrier()
      for i in src.sites:
        if cc[3][i] == srcT and
          (cc[0][i] and 1) == px(corner) and
          (cc[1][i] and 1) == py(corner) and
          (cc[2][i] and 1) == pz(corner):
          src{i}[color] := 1.0
    # invert color source to get propagator
    for massIdx in 0..<Masses:
      threads: props[color][corner][massIdx] := 0
      stag.solve(props[color][corner][massIdx], src, mass[massIdx], sp)
# echo "Saving the correlatos at tau=0.00"
qexLog "HISQ propagators done. Time: ", getElapsedTime()
var corrJSON = %* {}
if corrFile != "__NOOUTPUT__":
  let tauStr = formatFloat(0, ffDecimal, 2)
  corrJSON[tauStr] = %* {}
  qexLog "Measuring the propagators at tau=0.0"
  for color in 0..<colors:
    for pionIdx in 0..<mesons.len:
      for termIdx, term in mesons[pionIdx]:
        for forward in 0..<Corners:
          let 
            backward = shiftIdx(forward, term.delta)
            base = [px(forward), py(forward), pz(forward), srcT]
          var tmp = lo.ColorVector()
  
              # Build shifted + rephased copy of `props[backward]`.
          for propindex in 0..<mprops.len:
            shiftAndRephase(tmp, props[color][backward][mprops[propindex][0]], term.delta, 15 - term.phase, base)
      
                # Correlator for this pair, binned by timeslice
            let corr = dotByTimeslice(props[color][forward][mprops[propindex][1]], tmp)
            for t in 0..<nt: corrs[propindex][pionIdx][t] -= term.factor.float * corr[t]
  for propindex in 0..<mprops.len:
    for pionIdx in 0..<mesons.len:
      let combinedName = mesonname[pionIdx] & "_" & massname[propindex]
      corrJSON["0.00"][combinedName] = %* corrs[propindex][pionIdx]
  # let corrStr = pretty(corrJSON)
  # echo "Saving correlation functions at tau=0.0"
  # writeFile(corrFile & "_" & formatFloat(0.0, ffDecimal, 2), corrStr)
  
for flow in flowInfo.keys():
  flowInfo[flow]["filename"] = %* (flow & "_" & cfg & ".log")
echo "Starting fermion flow"
u.fermionFlow(props,flowInfo):
# save gauge field
  f.write(measurements.formatMeasurements(style = logStyle) & "\n")
  if abs(tau.round(1) - tau) < 1e-10: # if tau is 0.1, 0.2 ...
    if corrFile != "__NOOUTPUT__":
      echo "Measuring the propagators at tau= ", tau
      let tauStr = formatFloat(tau, ffDecimal, 2)
      corrJSON[tauStr] = %* {}
      for color in 0..<colors:
        for pionIdx in 0..<mesons.len:
          for termIdx, term in mesons[pionIdx]:
            for forward in 0..<Corners:
              let 
                backward = shiftIdx(forward, term.delta)
                base = [px(forward), py(forward), pz(forward), srcT]
              var tmp = lo.ColorVector()
      
              # Build shifted + rephased copy of `props[backward]`.
              for propindex in 0..<mprops.len:
                shiftAndRephase(tmp, props[color][backward][mprops[propindex][0]], term.delta, 15 - term.phase, base)
      
                # Correlator for this pair, binned by timeslice
                let corr = dotByTimeslice(props[color][forward][mprops[propindex][1]], tmp)
                for t in 0..<nt: corrs[propindex][pionIdx][t] -= term.factor.float * corr[t]
      
      # save output to convenient JSON format
      # echo "Saving correlation functions at tau= ", tau
      for propindex in 0..<mprops.len:
        for pionIdx in 0..<mesons.len:
          let combinedName = mesonname[pionIdx] & "_" & massname[propindex]
          corrJSON[tauStr][combinedName] = %* corrs[propindex][pionIdx]
echo "Saving correlation functions"
let corrStr = pretty(corrJSON)
writeFile(corrFile, corrStr)

qexFinalize()