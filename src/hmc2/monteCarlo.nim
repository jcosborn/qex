import qex

import strformat, strutils
import tables

import randomNumberGeneration

type
  IntegratorAtomType* = enum
    LF, LF1, LF2, LF3, MN0, MN2, MN4FP4, MN4FP5, MN4FV5
  AlgorithmType* = enum HamiltonianMonteCarlo

type
  IntegratorAtom* {.inheritable.} = object
    lambda*: float
    rho*, theta*: float
    case integrator*: IntegratorAtomType
      of LF, LF1, MN0, LF2, LF3: discard
      of MN2, MN4FP4: discard
      of MN4FP5, MN4FV5:
        vartheta*: float

    steps*: int
    b*: float

    stepsT*: seq[float]
    stepsV*: seq[float]

  MonteCarlo*[S] {.inheritable.} = object
    case algorithm*: AlgorithmType
      of HamiltonianMonteCarlo:
        hi*: float
        hf*: float

        c*: float
        tau*: float
        steps*: int

        serialRNG*: SerialRNG
        parallelRNG*: ParallelRNG

        p*: seq[S]
        f*: seq[S]
        bu*: seq[S]

        stepsT*: seq[float]
        stepsV*: Table[string, seq[float]]

converter toIntegratorAtomType*(s: string):
  IntegratorAtomType = parseEnum[IntegratorAtomType](s)

converter toAlgorithmType*(s: string): AlgorithmType = parseEnum[AlgorithmType](s)

proc constructIntegratorAtom*(self: var IntegratorAtom; steps: int) =
  self.steps = steps
  case self.integrator:
    of LF, LF1, MN0, LF2, LF3: discard
    of MN2: self.lambda = 0.1931833275037836
    of MN4FP4:
      self.rho = 0.1786178958448091
      self.theta = -0.06626458266981843
      self.lambda = 0.7123418310626056
    of MN4FP5:
      self.rho = 0.2750081212332419
      self.theta = -0.1347950099106792
      self.vartheta = -0.08442961950707149
      self.lambda = 0.3549000571574260
    of MN4FV5:
      self.rho = 0.2539785108410595
      self.theta = -0.03230286765269967
      self.vartheta = 0.08398315262876693
      self.lambda = 0.6822365335719091

proc setIntegratorAtom*(self: var IntegratorAtom; tau: float) =
  case self.integrator:
    of LF, LF1, MN0:
      self.stepsT = @[0.5, 0.5]
      self.stepsV = @[1.0, 0.0]
    of LF2:
      self.stepsT = @[0.25, 0.25, 0.25]
      self.stepsV = @[0.5, 0.5, 0.0]
    of LF3:
      self.stepsT = @[1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0]
      self.stepsV = @[1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0]
    of MN2:
      self.stepsT = @[self.lambda, 1.0 - 2.0 * self.lambda, self.lambda]
      self.stepsV = @[0.5, 0.5, 0.0]
    of MN4FP4:
      self.stepsT = @[
        self.rho, self.theta,
        1.0 - 2.0 * self.theta - 2.0 * self.rho,
        self.theta, self.rho
      ]
      self.stepsV = @[
        self.lambda, 0.5 - self.lambda,
        0.5 - self.lambda,
        self.lambda, 0.0
      ]
    of MN4FP5:
      self.stepsT = @[
        self.rho, self.theta,
        0.5 - self.theta - self.rho,
        0.5 - self.theta - self.rho,
        self.theta, self.rho
      ]
      self.stepsV = @[
        self.vartheta, self.lambda,
        1.0 - 2.0 * self.lambda - 2.0 * self.vartheta,
        self.lambda,
        self.vartheta, 0.0
      ]
    of MN4FV5:
      self.stepsT = @[
        0.0, self.rho, self.theta,
        1.0 - 2.0 * self.theta - 2.0 * self.rho, self.theta, self.rho
      ]
      self.stepsV = @[
        self.vartheta, self.lambda, 0.5 - self.lambda - self.vartheta,
        0.5 - self.lambda - self.vartheta, self.lambda, self.vartheta
      ]
    else: discard

  for index in 0..<self.stepsT.len:
    self.stepsT[index] = self.stepsT[index] * tau / float(self.steps)
    self.stepsV[index] = self.stepsV[index] * tau / float(self.steps)

  var
    stepsT = newSeq[float]()
    stepsV = newSeq[float]()

  for step in 0..<self.steps:
    for index in 0..<self.stepsT.len:
      stepsT.add self.stepsT[index]
      stepsV.add self.stepsV[index]
  
  self.stepsT = stepsT
  self.stepsV = stepsV