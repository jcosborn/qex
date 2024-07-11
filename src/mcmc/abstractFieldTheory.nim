import qex
import layout

import json
import options
import os

import abstractActions
import monteCarlo
import typeUtilities
import randomNumberGeneration

export qex, json

type
  AbstractFieldTheory*[L:static[int],S,T,U,V,W] = object of MonteCarlo[S]
    l*: Layout[L]
    actions*: Table[string, AbstractAction[L,S,T,U,V,W]]
    nMatterActions*: int
    start*: string

proc u(self: AbstractFieldTheory): auto = 
  result = self.actions["gauge"].fields["gauge"].u

proc checkJSON(info: JsonNode): JsonNode =
  result = checkMonteCarloAlgorithm(info)
  if not result.hasKey("lattice-geometry"): qexError "Must specify lattice geometry"
  if not result.hasKey("rank-geometry"): qexError "Must specify rank geometry"
  case result["monte-carlo-algorithm"].getStr():
    of "HamiltonianMonteCarlo":
      if not result.hasKey("trajectory-length"): 
        qexError "Must specify trajectory length"
      if not result.hasKey("serial-random-number-generator"):
        result["serial-random-number-generator"] = %* "MILC"
      if not result.hasKey("serial-random-number-seed"):
        result["serial-random-number-generator"] = %* 987654321
      if not result.hasKey("parallel-random-number-generator"):
        result["parallel-random-number-generator"] = %* "MILC"
      if not result.hasKey("parallel-random-number-seed"):
        result["parallel-random-number-generator"] = %* 987654321
      if not result.hasKey("start"): result["start"] = %* "none"
    of "HeatbathOverrelax": discard
    else: discard

proc checkJSON(self: AbstractFieldTheory, info: JsonNode): JsonNode =
  result = parseJson("{}")
  for key in info.keys(): result[key] = info[key]
  case self.algorithm:
    of HamiltonianMonteCarlo: 
      result["monte-carlo-algorithm"] = %* "hamiltonian-monte-carlo"
    of HeatbathOverrelax: 
      result["monte-carlo-algorithm"] = %* "heatbath-overrelax"

proc newAbstractFieldTheory[L:static[int],S,T,U,V,W](
    l: Layout[L]; 
    info: JsonNode;
    s: typedesc[S];
    t: typedesc[T];
    u: typedesc[U];
    v: typedesc[V];
    w: typedesc[W]
  ): AbstractFieldTheory[L,S,T,U,V,W] =
  let algorithm = toMonteCarloType(info["monte-carlo-algorithm"].getStr())
  result = AbstractFieldTheory[L,S,T,U,V,W](algorithm: algorithm)
  result.l = l
  case result.algorithm:
    of HamiltonianMonteCarlo:
      result.bu = result.l.newGauge()
      result.p = result.l.newGauge()
      result.f = result.l.newGauge()
      result.serialRNG = newSerialRNG(
        info["serial-random-number-generator"].getStr(),
        (info["serial-random-number-seed"].getInt()).uint64
      )   
      result.parallelRNG = result.l.newParallelRNG(
        info["parallel-random-number-generator"].getStr(),
        (info["parallel-random-number-seed"].getInt()).uint64
      )
      result.stepsT = newSeq[float]()
      result.stepsV = initTable[string,Table[string,seq[seq[float]]]]()
      result.tau = info["trajectory-length"].getFloat()
      result.start = info["start"].getStr()
    of HeatbathOverrelax: discard
  result.actions = initTable[string, AbstractAction[L,S,T,U,V,W]]()
  result.nMatterActions = 0

proc newFieldTheory*(fieldTheoryInfo: JsonNode): auto = 
  var 
    info = checkJSON(fieldTheoryInfo)
    latticeGeometry = newSeq[int]()
    rankGeometry = newSeq[int]()

  for el in info["lattice-geometry"].getElems(): latticeGeometry.add el.getInt()
  for el in info["rank-geometry"].getElems(): rankGeometry.add el.getInt()

  let l = newLayout(latticeGeometry, rankGeometry)

  result = l.newAbstractFieldTheory( 
    info, l.typeS, l.typeT, l.typeU, l.typeV, l.typeW
  )

proc setGaugeAction*(self: var AbstractFieldTheory; gaugeInfo: JsonNode) =
  var info = self.checkJSON(gaugeInfo)
  self.actions["gauge"] = self.l.newGaugeAction(info)
  case self.start:
    of "cold": unit(self.u)
    of "hot": self.parallelRNG.warm(self.u)
    else: discard

proc addMatterAction*(
    self: var AbstractFieldTheory; 
    info: JsonNode;
    name: string = "matter"
  ) =
  let defaultName = "matter" & $(self.nMatterActions)
  case name:
    of "matter": self.actions[defaultName] = self.l.newMatterAction(info)
    else: self.actions[name] = self.l.newMatterAction(info)
  self.nMatterActions = self.nMatterActions + 1

proc addStaggeredFermion*(
    self: var AbstractFieldTheory; 
    fermionInfo: JsonNode;
    nf: int;
    action: string
  ) =
  var info = self.checkJSON(fermionInfo)
  for _ in 0..<nf: self.actions[action].addStaggeredFermion(info)

proc addStaggeredBoson*(
    self: var AbstractFieldTheory; 
    bosonInfo: JsonNode;
    nb: int;
    action: string
  ) =
  var info = self.checkJSON(bosonInfo)
  for _ in 0..<nb: self.actions[action].addStaggeredBoson(info)

proc addStaggeredHasenbuschFermion*(
    self: var AbstractFieldTheory; 
    hasenbuschInfo: seq[JsonNode];
    nh: int;
    action: string
  ) =
  var info = newSeq[JsonNode]()
  for index in 0..<hasenbuschInfo.len: 
    info.add self.checkJSON(hasenbuschInfo[index])
  for _ in 0..<nh: self.actions[action].addStaggeredHasenbuschFermion(info)

proc compile*(self: var AbstractFieldTheory) =
  case self.algorithm:
    of HamiltonianMonteCarlo:
      # Initialize Hamiltonian Monte Carlo steps
      for action, _ in self.actions:
        self.stepsV[action] = initTable[string,seq[seq[float]]]()
        for field, _ in self.actions[action].fields:
          self.actions[action].fields[field].initHamiltonianMonteCarlo(self.tau)
          self.stepsV[action][field] = newSeq[seq[float]]()
          for index in 0..<self.actions[action].fields[field].stepsV.len:
            self.stepsV[action][field].add newSeq[float]()

      # Fill in steps for all fields in each action
      while true:
        # Get smallest value of u update step
        var dtaus = newSeq[float]()
        for action, _ in self.actions:
          for field, _ in self.actions[action].fields:
            for index in 0..<self.actions[action].fields[field].stepsT.len:
              dtaus.add self.actions[action].fields[field].stepsT[index][0]
        let dtau = min(dtaus)
        self.stepsT.add dtau

        # Decide whether or not to add momentum update step or insert "do nothing"
        for action, _ in self.actions:
          for field, _ in self.actions[action].fields:
            for index in 0..<self.actions[action].fields[field].stepsV.len:
              let stepT = self.actions[action].fields[field].stepsT[index][0]
              if abs(dtau - stepT) < epsilon(float):
                let stepV = self.actions[action].fields[field].stepsV[index][0]
                self.stepsV[action][field][index].add stepV
                self.actions[action].fields[field].stepsT[index].delete(0)
                self.actions[action].fields[field].stepsV[index].delete(0)
              else:
                self.stepsV[action][field][index].add 0.0
                self.actions[action].fields[field].stepsT[index][0] -= dtau

        # Finish when there are no more updates left to check
        var totalSteps = 0
        for action, _ in self.actions:
          for field, _ in self.actions[action].fields:
            for index in 0..<self.actions[action].fields[field].stepsV.len:
              let fieldSteps = self.actions[action].fields[field].stepsV[index].len
              totalSteps = totalSteps + fieldSteps
        if totalSteps == 0: break
    of HeatbathOverrelax: discard

proc pdagp(p: auto): float = 
  var p2: float
  threads:
    var p2t = 0.0
    for mu in 0..<p.len: p2t += p[mu].norm2
    threadBarrier()
    threadMaster: p2 = p2t
  result = p2

proc hamiltonian*(
    self: var AbstractFieldTheory; 
    heatbath: bool = false
  ) =
  # Initialize table
  var result = initTable[string,float]()

  # Momentum heatbath
  if heatbath: self.parallelRNG.randomTAHGaussian(self.p)

  # Fermion heatbath + calculation of initial action
  result["kinetic"] = 0.5 * pdagp(self.p) - 16.0 * self.l.physVol
  result["hamiltonian"] = result["kinetic"]
  for action, _ in self.actions:
    let h = case self.actions[action].action
      of PureGauge: {"gauge": self.actions[action].action()}.toTable()
      of PureMatter:
        if heatbath: 
          self.actions[action].action(self.u, parallelRNG = some(self.parallelRNG))
        else: self.actions[action].action(self.u)
    for ak, av in h: 
      result[action & "_" & ak] = av
      result["hamiltonian"] = result["hamiltonian"] + av

  # Save result of calculation
  if heatbath: self.hi = result
  else: self.hf = result

proc setCoordinate(u: auto; v: auto) =
  threads:
    for mu in 0..<u.len:
      for s in u[mu]: u[mu][s] := v[mu][s]

proc updateCoordinate(u: auto; p: auto; dtau: float) =
  threads:
    for mu in 0..<u.len:
      for s in u[mu]: u[mu][s] := exp(dtau*p[mu][s])*u[mu][s]

proc evolve*(self: var AbstractFieldTheory) =
  case self.algorithm:
    of HamiltonianMonteCarlo:
      # Save gauge field before evolution
      self.bu.setCoordinate(self.u)

      # Do molecular dynamics update
      for index in 0..<self.stepsT.len:
        # Coordinate update
        let dtauT = self.stepsT[index]
        if dtauT != 0.0: updateCoordinate(self.u, self.p, dtauT)

        # Set momentum time steps
        for action, _ in self.actions:
          for field, _ in self.actions[action].fields:
            for idx in 0..<self.actions[action].fields[field].stepsV.len:
              let dtauV = self.stepsV[action][field][idx][index]
              self.actions[action].fields[field].dtau[idx] = dtauV

        # Momentum update
        for action, _ in self.actions:
          case self.actions[action].action:
            of PureGauge: self.actions[action].updateMomentum(self.f, self.p)
            of PureMatter: 
              self.actions[action].updateMomentum(self.u, self.f, self.p)
    of HeatbathOverrelax: discard

proc reunit*(u: auto) =
  threads: u.projectSU

proc metropolis*(self: var AbstractFieldTheory): bool =
  result = false
  case self.algorithm:
    of HamiltonianMonteCarlo:
      let dH = self.hf["hamiltonian"] - self.hi["hamiltonian"]
      if self.serialRNG.uniform <= exp(-dH): 
        reunit(self.u)
        result = true
      else: self.u.setCoordinate(self.bu)
    of HeatbathOverrelax: discard

proc read*(
    self: var AbstractFieldTheory; 
    fn: string;
    onlyGauge: bool = false
  ) =
  let
    gfn = fn & ".lat"
    pfn = fn & ".rng"
    sfn = fn & ".global_rng"
  if fileExists(gfn):
    if 0 != self.u.loadGauge(gfn): qexError "unable to read " & gfn
    reunit(self.u)
  case self.algorithm:
    of HamiltonianMonteCarlo:
      if not onlyGauge:
        if fileExists(pfn): self.parallelRNG.readRNG(pfn)
        else: qexError "unable to read " & pfn
        if fileExists(sfn): 
          self.serialRNG.readRNG(sfn)
          echo "read" & sfn
        else: qexError "unable to read " & sfn
    of HeatbathOverrelax: discard

proc write*(
    self: var AbstractFieldTheory;
    fn: string; 
    onlyGauge: bool = false
  ) =
  let
    gfn = fn & ".lat"
    pfn = fn & ".rng"
    sfn = fn & ".global_rng"
  if 0 != self.u.saveGauge(gfn): qexError "unable to write " & gfn
  case self.algorithm:
    of HamiltonianMonteCarlo:
      if not onlyGauge:
        self.parallelRNG.writeRNG(pfn)
        self.serialRNG.writeRNG(sfn)
    of HeatbathOverrelax: discard

if isMainModule:
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
    matterActionParams = %* {
      "smearing": "nHYP",
      "smearing-coefficients": @[0.4, 0.5, 0.5],
      "boundary-conditions": "aaaa"
    }
    fermionParams = %* {
      "mass": 0.1,
      "steps": 10,
      "integrator": "MN2",
    }
    bosonParams = %* {
      "mass": 0.75,
      "steps": 5,
      "integrator": "MN2"
    }
    hasenbuschParams = @[fermionParams, bosonParams, bosonParams]
      
  # Creat field theory
  var fieldTheory = newFieldTheory(fieldTheoryParams)

  # Set gauge action, then add matter action
  fieldTheory.setGaugeAction(gaugeActionParams)
  #fieldTheory.addMatterAction(matterActionParams, name = "fermion/boson")

  # Add fermions & bosons to matter action
  #fieldTheory.addStaggeredFermion(fermionParams, 1): "fermion/boson"
  #fieldTheory.addStaggeredBoson(bosonParams, 8): "fermion/boson"
  #fieldTheory.addStaggeredHasenbuschFermion(hasenbuschParams, 1): "fermion/boson"

  # Compile field theory for use in Monte Carlo
  fieldTheory.compile

  #echo fieldTheory.stepsT
  #echo fieldTheory.stepsV

  for step in 0..<100000:
    # Get initial hamiltonian
    fieldTheory.hamiltonian(heatbath = true)
    #echo fieldTheory.hi["hamiltonian"]

    # Evolve gauge field with molecular dynamics
    fieldTheory.evolve

    # Get final hamiltonian
    fieldTheory.hamiltonian
    #echo fieldTheory.hf["hamiltonian"]

    # Do metropolis accept/reject
    let 
      acc = fieldTheory.metropolis
      dH = fieldTheory.hf["hamiltonian"] - fieldTheory.hi["hamiltonian"]
    echo step, " ", acc, " ", dH

  qexFinalize()