import qex

import randomNumberGeneration

import strutils
import json

type
  IntegratorAtomType* = enum
    CustomIntegrator, LF, LF1, LF2, LF3, MN0, MN2, MN4FP4, MN4FP5, MN4FV5
  MonteCarloType* = enum HamiltonianMonteCarlo, HeatbathOverrelax

type
  MonteCarloAtom* {.inheritable.} = object
    case algorithm*: MonteCarloType
      of HamiltonianMonteCarlo:
        lmbda*, rho*, theta*, vartheta*: float
        integrator*: IntegratorAtomType
        steps*: int
        dtau*: seq[float]
        stepsT*: seq[seq[float]]
        stepsV*: seq[seq[float]]
      of HeatbathOverrelax: discard

  MonteCarlo*[S] {.inheritable.} = object
    case algorithm*: MonteCarloType
      of HamiltonianMonteCarlo:
        tau*: float
        hi*,hf*: Table[string,float]
        stepsT*: seq[float]
        stepsV*: Table[string,Table[string,seq[seq[float]]]]
        bu*,p*,f*: seq[S]
        serialRNG*: SerialRNG
        parallelRNG*: ParallelRNG
      of HeatbathOverrelax: discard

converter toIntegratorAtomType*(s: string):
  IntegratorAtomType = parseEnum[IntegratorAtomType](s)

converter toMonteCarloType*(s: string): 
  MonteCarloType = parseEnum[MonteCarloType](s)

proc checkMonteCarloAlgorithm*(info: JsonNode): JsonNode = 
  result = parseJson("{}")
  for key in info.keys(): result[key] = info[key]
  if not result.hasKey("monte-carlo-algorithm"):
    result["monte-carlo-algorithm"] = %* "HamiltonianMonteCarlo"
  else:
    case result["monte-carlo-algorithm"].getStr():
      of "hmc", "HMC", "hamiltonian-monte-carlo", "hybrid-monte-carlo":
        result["monte-carlo-algorithm"] = %* "HamiltonianMonteCarlo"
      of "heatbath", "heatbath-overrelax", "heatbath+overrelax":
        result["monte-carlo-algorithm"] = %* "HeatbathOverrelax"
      else:
        var algorithm = result["monte-carlo-algorithm"].getStr()
        qexError algorithm & " is not a valid algorithm"

proc initHamiltonianMonteCarlo*(self: var MonteCarloAtom; tau: float) =
  # Set steps for default integrators if not using "custom integrator"
  case self.integrator:
    of LF, LF1, MN0:
      self.stepsT = @[@[0.5, 0.5]]
      self.stepsV = @[@[1.0, 0.0]]
    of LF2:
      self.stepsT = @[@[0.25, 0.25, 0.25]]
      self.stepsV = @[@[0.5, 0.5, 0.0]]
    of LF3:
      self.stepsT = @[@[1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0]]
      self.stepsV = @[@[1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0]]
    of MN2:
      self.stepsT = @[@[self.lmbda, 1.0 - 2.0 * self.lmbda, self.lmbda]]
      self.stepsV = @[@[0.5, 0.5, 0.0]]
    of MN4FP4:
      self.stepsT = @[@[
        self.rho, self.theta,
        1.0 - 2.0 * self.theta - 2.0 * self.rho,
        self.theta, self.rho
      ]]
      self.stepsV = @[@[
        self.lmbda, 0.5 - self.lmbda,
        0.5 - self.lmbda,
        self.lmbda, 0.0
      ]]
    of MN4FP5:
      self.stepsT = @[@[
        self.rho, self.theta,
        0.5 - self.theta - self.rho,
        0.5 - self.theta - self.rho,
        self.theta, self.rho
      ]]
      self.stepsV = @[@[
        self.vartheta, self.lmbda,
        1.0 - 2.0 * self.lmbda - 2.0 * self.vartheta,
        self.lmbda,
        self.vartheta, 0.0
      ]]
    of MN4FV5:
      self.stepsT = @[@[
        0.0, self.rho, self.theta,
        1.0 - 2.0 * self.theta - 2.0 * self.rho, self.theta, self.rho
      ]]
      self.stepsV = @[@[
        self.vartheta, self.lmbda, 0.5 - self.lmbda - self.vartheta,
        0.5 - self.lmbda - self.vartheta, self.lmbda, self.vartheta
      ]]
    else: discard

  # Collect full set of steps
  let dtau = tau/float(self.steps)
  for index in 0..<self.stepsT.len:
    for indexT in 0..<self.stepsT[index].len:
      self.stepsT[index][indexT] = self.stepsT[index][indexT] * dtau
  for index in 0..<self.stepsV.len:
    for indexV in 0..<self.stepsV[index].len:
      self.stepsV[index][indexV] = self.stepsV[index][indexV] * dtau

  var
    stepsT: seq[seq[float]]
    stepsV: seq[seq[float]]
  for _ in self.stepsT: stepsT.add newSeq[float]()
  for _ in self.stepsV: stepsV.add newSeq[float]()

  for step in 0..<self.steps:
    for index in 0..<self.stepsT.len:
      for indexT in 0..<self.stepsT[index].len:  
        stepsT[index].add self.stepsT[index][indexT]
    for index in 0..<self.stepsV.len:
      for indexV in 0..<self.stepsV[index].len: 
        stepsV[index].add self.stepsV[index][indexV]
  
  # Save integrator across all steps while also cleaning up redundancies
  for index in 0..<stepsT.len:
    self.stepsT[index] = newSeq[float]()
    self.stepsV[index] = newSeq[float]()
    for indexT in 0..<stepsT[index].len-1:
      let 
        dT = abs(stepsT[index][indexT] - stepsT[index][indexT+1])
        dV = abs(stepsV[index][indexT] - stepsV[index][indexT+1])
      if (dT < epsilon(float)) and (stepsV[index][indexT] == 0.0):
        self.stepsT[index].add 2.0 * stepsT[index][indexT]
        stepsT[index][indexT+1] = 0.0
      else: self.stepsT[index].add stepsT[index][indexT]
      if (dV < epsilon(float)) and (stepsT[index][indexT+1] == 0.0):
        self.stepsV[index].add 2.0 * stepsV[index][indexT]
        stepsV[index][indexT+1] = 0.0
      else: self.stepsV[index].add stepsV[index][indexT]
    self.stepsT[index].add stepsT[index][^1]
    self.stepsV[index].add stepsV[index][^1]