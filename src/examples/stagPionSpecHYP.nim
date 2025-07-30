## nHYP-smeared staggered pion spectrum
## Author: Xiao-Yong Jin
## Author: Curtis Taylor Peterson

import qex
import gauge/[gaugefix, hypsmear]
import physics/[stagSolve]
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
    nhyp = HypCoefs(alpha1: 0.4, alpha2: 0.5, alpha3: 0.5)

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

  # smear and rephase
  if smearGauge:
    tic "nHYP gauge smearing"
    var sg = lo.newGauge()
    nhyp.smear(g, sg)
    threads: g := sg
    qexLog "nHYP gauge smearing done. Time: ", getElapsedTime()
  g.rephase()

  # instantiate staggered Dirac operator
  let stag = g.newStag

  # instantiate solver parameters
  var sp = initSolverParams()
  sp.r2req = r2req
  sp.maxits = maxits

  #[ spectrum ]#

  # build meson operators & select pions as subset
  let 
    mesons = buildMesonOps()
    pions = [
      mesons[0],  # I ⊗ I         (P = +σ) I     π_I
      mesons[1],  # γ_5 ⊗ ξ_5     (P = -σ) II    π_5
      mesons[6],  # γ_5 ⊗ ξ_μ ξ_5 (P = +σ) VII   π_i5
      mesons[7],  # I ⊗ ξ_μ       (P = -σ) VIII  π_0i
      mesons[12], # I ⊗ ξ_μ ξ_ν   (P = -σ) XIII  π_ij
      mesons[13], # γ_5 ⊗ ξ_μ     (P = -σ) XIV   π_i
      mesons[16], # I ⊗ ξ_0 ξ_5   (P = -σ) XVII  π_05
      mesons[17], # γ_5 ⊗ ξ_0     (P = +σ) XVIII π_0
    ]

  # initialize propagators and correlators
  var 
    props: array[Corners, typeof(lo.ColorVector)]
    corrs = newSeq[seq[float]](pions.len)
  for corner in 0..<Corners: props[corner] = lo.ColorVector()
  for pionIdx in 0..<pions.len: corrs[pionIdx].newSeq(nt)

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
    for pionIdx in 0..<pions.len:
      for termIdx, term in pions[pionIdx]:
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
        "pi_I":  corrs[0],
        "pi_5":  corrs[1],
        "pi_i5": corrs[2],
        "pi_0i": corrs[3],
        "pi_ij": corrs[4],
        "pi_i":  corrs[5],
        "pi_05": corrs[6],
        "pi_0":  corrs[7]
      }
      corrStr = pretty(corrJSON)
    writeFile(corrFile, corrStr)

  qexFinalize()