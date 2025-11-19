## HISQ-smeared staggered pion spectrum
## Author: Xiao-Yong Jin
## Author: Curtis Taylor Peterson

import qex
import gauge/[gaugefix, hisqsmear]
import physics/[stagSolve]
import maths/[matproject]
import observables/[nlStagMeson]

import std/[os, times, sugar, tables] 
import std/[parseopt, strutils, strformat]
import std/[json, sequtils]
import std/[math]

const
  DefaultCoulombGaugeFixEpsilon = 1e-8
  DefaultCoulombGaugeFixRelax = 1.75

const
  DefaultConjugateGradientTolerance = 1e-20
  DefaultConjugateGradientMaximumIterations = 10000
  DefaultConjugateGradientVerbosty = 1

const Corners = 8

let
  defaultInputs = %* {
    "lattice-geometry": [8,8,8,16], # lattice geometry
    #"rank-geometry": [1,1,2,2], # MPI rank geometry -- if not specified, guessed
    #"simd-geometry": [2,2,1,2], # SIMD (local) geometry -- if not specified, guessed
    "solver": {
      "minimum-squared-residual": DefaultConjugateGradientTolerance,
      "maximum-iterations": DefaultConjugateGradientMaximumIterations,
      "verbosity": DefaultConjugateGradientVerbosty
    },
    "coulomb-gauge-fix": { # specify gauge fixing information
      "eps": DefaultCoulombGaugeFixEpsilon, # stopping criterion for gauge fixing
      "relax": DefaultCoulombGaugeFixRelax # gauge fixing relaxation factor
    },
    "spectrum": {
      "mass": 0.1,
      "source-time": 0
    }
  }

# rephases links
proc rephase[U](u: seq[U]) = 
  threads:
    u.setBC()
    threadBarrier()
    u.stagPhase()

# Reads command line information
proc readCMD: JsonNode = 
  var cmd = initOptParser()
  result = parseJson("{}")
  while true:
    cmd.next()
    case cmd.kind:
      of cmdShortOption,cmdLongOption,cmdArgument:
        try: result[cmd.key] = %* parseInt(cmd.val)
        except ValueError:
          try: result[cmd.key] = %* parseFloat(cmd.val)
          except ValueError: result[cmd.key] = %* cmd.val
      of cmdEnd: break

# gets float or integer sequence from JsonNode object
proc getSeq[T](input: JsonNode, t: typedesc[T]): seq[T] = 
  result = newSeq[T]()
  for elem in input.getElems(): 
    let telem = elem.getFloat().T
    result.add telem

# reads JSON file
proc readJSON(fn: string): JsonNode = fn.parseFile

# reads gauge configuration
proc readGauge[U](u: seq[U]; config: string) =
  if fileExists(config):
    if 0 != u.loadGauge(config): qexError "unable to read " & config
    else: discard
  else: qexError config & " does not exist"

# constructs lattice layout
proc newLayout(info: JsonNode): auto =
  let latLayout = info["lattice-geometry"].getSeq(int)
  result = case info.hasKey("rank-geometry")
    of true:
      case info.hasKey("simd-geometry"):
        of true:
          let rgeom = info["rank-geometry"].getSeq(int)
          let sgeom = info["simd-geometry"].getSeq(int)
          newLayout(latLayout, VLEN, rgeom, sgeom)
        of false: newLayout(latLayout, info["rank-geometry"].getSeq(int))
    of false: newLayout(latLayout)
  assert(latLayout.len == 4)

# Helper to map (operator n, local index i) → global meson index
proc mesonIndex[M](mesons: M; n, i: int):int =
    var off = 0
    for k in 0..<n: off += mesons[k].len
    result = off + i

when isMainModule:
  qexInit()

  # print timing information & rank/thread info
  echo "rank ", myRank, "/", nRanks
  threads: echo "thread ", threadNum, "/", numThreads

  #[ setup ]#

  # command line information
  let cmd = readCMD()

  # information from stored JSON file on disk
  let info = case cmd.hasKey("json")
    of true: readJSON(cmd["json"].getStr())
    of false: defaultInputs

  # location of saved JSON file
  let corrFile = case cmd.hasKey("output")
    of true: cmd["output"].getStr()
    of false: "__NOOUTPUT__"

  let
    # conjugate gradient parameters
    r2req = info["solver"]["minimum-squared-residual"].getFloat()
    maxits = info["solver"]["maximum-iterations"].getInt()

    # gauge fixing parameters
    fixGauge: bool = true
    gfstop = info["coulomb-gauge-fix"]["eps"].getFloat()
    gforf = info["coulomb-gauge-fix"]["relax"].getFloat()

    # nhyp smearing parameters
    smearGauge: bool = true
    hisq = newHISQ(0.0, 1.0)

    # lattice layout
    lo = info.newLayout()

    # source time
    srcT = info["spectrum"]["source-time"].getInt()

    # time extent
    nt = info["lattice-geometry"].getSeq(int)[^1]

  # read in gauge configuration or fill with random links
  var g = lo.newGauge()
  case cmd.hasKey("configuration"):
    of true: g.readGauge(cmd["configuration"].getStr())
    of false: g.random 

  # reunitarize
  threads:
    block:
      let d = g.checkSU
      echo "unitary deviation avg: ", d.avg, "  max: ", d.max
    threadBarrier()
    g.projectSU
    threadBarrier()
    block:
      let d = g.checkSU
      echo "new unitary deviation avg: ", d.avg, "  max: ", d.max
  g.echoPlaq

  # Coulomb gauge fix
  if fixGauge:
    tic "Coulomb gauge fixing"

    var tmat = lo.ColorMatrix()
    threads: tmat := 1
    getGaugeFixTransform(tmat, g, @[0,1,2], gfstop, gforf, verb=0)

    var fg = lo.newGauge
    fg.gaugeTransform(g, tmat)
    threads: g := fg

    qexLog "Coulomb gauge fixing done. Time: ", getElapsedTime()

  # rephase and smear
  var sg = lo.newGauge()
  var sgl = lo.newGauge()
  g.rephase()
  if smearGauge:
    tic "HISQ gauge smearing"
    discard hisq.smearGetForce(g, sg, sgl, displayPerformance = true)
    qexLog "HISQ gauge smearing done. Time: ", getElapsedTime()

  # instantiate staggered Dirac operator
  let stag = newStag3(sg, sgl)

  # instantiate solver parameters
  var sp = initSolverParams()
  sp.r2req = r2req
  sp.maxits = maxits

  #[ spectrum ]#

  # build meson operators & select mesons as subset
  let 
    mesonOps = buildMesonOps()
    mesons = [
      mesonOps[0],  # γ_0 γ_5 ⊗ γ_0 γ_5 (P = +σ) I     π_05  <--+ 
      mesonOps[1],  # γ_5 ⊗ ξ_5         (P = -σ) II    π_5      |
      mesonOps[6],  # γ_5 ⊗ ξ_μ ξ_5     (P = +σ) VII   π_i5     |
      mesonOps[7],  # γ_0 γ_5 ⊗ ξ_μ ξ_ν (P = -σ) VIII  π_ij     | scalar &
      mesonOps[12], # γ_0 γ_5 ⊗ ξ_μ     (P = -σ) XIII  π_0i     | pseudoscalar
      mesonOps[13], # γ_5 ⊗ ξ_μ         (P = -σ) XIV   π_i      |
      mesonOps[16], # γ_0 γ_5 ⊗ I       (P = -σ) XVI   π_I      |
      mesonOps[17], # γ_5 ⊗ ξ_0         (P = +σ) XVII  π_0   <--+
      mesonOps[2],  # γ_0 γ_k ⊗ ξ_0 ξ_μ (P = +σ) XIII  ρ_0i  <--+
      mesonOps[3],  # γ_μ ⊗ ξ_μ         (P = -σ) XIV   ρ_i      |
      mesonOps[4],  # γ_μ ⊗ I           (P = +σ) XV    ρ_I      |
      mesonOps[5],  # γ_0 γ_μ ⊗ ξ_0     (P = +σ) VI    ρ_0      |
      mesonOps[8],  # γ_0 γ_μ ⊗ ξ_μ ξ_5 (P = +σ) XIV   ρ_i5     |
      mesonOps[9],  # γ_μ ⊗ ξ_μ ξ_ν     (P = +σ) X     ρ_ij     | vector &
      mesonOps[10], # γ_0 γ_μ ⊗ ξ_0 ξ_5 (P = -σ) XI    ρ_05     | pseudovector
      mesonOps[11], # γ_0 γ_μ ⊗ ξ_5     (P = -σ) XII   ρ_5      |
      mesonOps[14], # γ_μ ⊗ ξ_μ         (P = -σ) XIV   ρ_i*     |
      mesonOps[15], # γ_0 γ_k ⊗ ξ_0 ξ_μ (P = -σ) XIII  ρ_0i*    |
      mesonOps[18], # γ_0 γ_μ ⊗ ξ_μ ξ_5 (P = -σ) XVIII ρ_i5*    |
      mesonOps[19], # γ_μ ⊗ ξ_μ ξ_ν     (P = +σ) XX    ρ_ij* <--+
    ]

  # initialize propagators and correlators
  var 
    props: array[Corners, typeof(lo.ColorVector)]
    corrs = newSeq[seq[float]](mesons.len)
  for corner in 0..<Corners: props[corner] = lo.ColorVector()
  for pionIdx in 0..<mesons.len: corrs[pionIdx].newSeq(nt)

  # get spectrum and do it at each color
  const colors = g[0][0].nrows
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
      threads: props[corner] := 0
      stag.solve(props[corner], src, info["spectrum"]["mass"].getFloat(), sp)

    # do contractions
    for pionIdx in 0..<mesons.len:
      for termIdx, term in mesons[pionIdx]:
        for forward in 0..<Corners:
          let 
            backward = shiftIdx(forward, term.delta)
            base = [px(forward), py(forward), pz(forward), srcT]
          var tmp = lo.ColorVector()

          # Build shifted + rephased copy of `props[backward]`.
          shiftAndRephase(tmp, props[backward], term.delta, 15 - term.phase, base)

          # Correlator for this pair, binned by timeslice
          let corr = dotByTimeslice(props[forward], tmp)
          for t in 0..<nt: corrs[pionIdx][t] -= term.factor.float * corr[t]

  # save output to convenient JSON format
  if corrFile != "__NOOUTPUT__":
    let 
      corrJSON = %* {
        "pi_05": corrs[0],
        "pi_5":  corrs[1],
        "pi_i5": corrs[2],
        "pi_ij": corrs[3],
        "pi_0i": corrs[4],
        "pi_i":  corrs[5],
        "pi_I":  corrs[6],
        "pi_0":  corrs[7],
        "rho_0i":  corrs[8],
        "rho_i":   corrs[9],
        "rho_I":   corrs[10],
        "rho_0":   corrs[11],
        "rho_i5":  corrs[12],
        "rho_ij":  corrs[13],
        "rho_05":  corrs[14],
        "rho_5":   corrs[15],
        "rho_i*":  corrs[16],
        "rho_0i*": corrs[17],
        "rho_i5*": corrs[18],
        "rho_ij*": corrs[19]
      }
      corrStr = pretty(corrJSON)
    writeFile(corrFile, corrStr)

  qexFinalize()