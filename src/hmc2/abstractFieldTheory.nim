import qex
import layout
import gauge/hypsmear
import gauge/stoutsmear

import tables
import json
import strutils
import math

import abstractFields
import randomNumberGeneration

export abstractFields

type
  FieldTheoryType* = enum PureGauge, StaggeredGaugeFermion, WilsonGaugeFermion
  SmearingType* = enum None, nHYP, Stout

type
  AbstractFieldTheory*[L:static[int],S,T,U,V,W] {.inheritable.} = object of MonteCarlo[S]
    l*: Layout[L]

    fields*: Table[string, AbstractField[S,T,U]]

    nFerm*, nBos*: int
    matterBoundaryConditions*: string
    case fieldTheory*: FieldTheoryType
      of PureGauge: discard
      of StaggeredGaugeFermion:
        sPsi*: T
        sD*: DiracOperator[S,T,V,W]
      of WilsonGaugeFermion:
        wPsi*: U
        wD*: DiracOperator[S,T,V,W]

    smearedForce*: proc(f, chain: seq[S])
    case smearing*: SmearingType
      of nHYP:
        uHYP*: seq[S]
        hyp*: HypCoefs
        info*: PerfInfo
      of Stout:
        uStout*: seq[S]
        stout*: StoutSmear[seq[S]]
      of None: discard

proc u*(self: AbstractFieldTheory): auto = self.fields["gauge"].u

proc su*(self: AbstractFieldTheory): auto =
  case self.smearing:
    of nHYP: result = self.uHYP
    else: result = self.uStout

proc newDiracOperator*(self: AbstractFieldTheory; fieldType: FieldType): auto =
  case self.smearing:
    of nHYP, Stout: result = newDiracOperator(self.su, fieldType)
    else: result = newDiracOperator(self.u, fieldType)

proc new(
    self: var AbstractFieldTheory;
    l: Layout;
    gauge: AbstractField;
    smearingCoeffs: seq[float];
    matterBoundaryConditions: string;
    sGen, pGen: string;
    sSeed, pSeed: uint64
  ) =
  self.fields = initTable[string, AbstractField]()
  self.fields["gauge"] = gauge

  self.l = l

  case self.smearing:
    of nHYP:
      self.uHYP = self.l.newGauge()
      self.hyp = HypCoefs(
        alpha1: smearingCoeffs[0],
        alpha2: smearingCoeffs[1],
        alpha3: smearingCoeffs[2]
      )
    of Stout:
      self.uStout = self.l.newGauge()
      self.stout = self.l.newStoutSmear(smearingCoeffs[0])
    else: discard

  self.nFerm = 0
  self.nBos = 0
  self.matterBoundaryConditions = matterBoundaryConditions
  case self.fieldTheory:
    of StaggeredGaugeFermion:
      self.sPsi = self.l.ColorVector()
      self.sD = self.newDiracOperator(StaggeredFermion)
      for mu in 0..<self.u.len: self.sD.stagShifter[mu] = newShifter(self.sPsi, mu, 1)
    of WilsonGaugeFermion:
      self.wPsi = self.l.DiracFermion()
      self.wD = self.newDiracOperator(WilsonFermion)
    of PureGauge: discard

  case self.algorithm:
    of HamiltonianMonteCarlo:
      self.p = self.l.newGauge()
      self.f = self.l.newGauge()
      self.bu = self.l.newGauge()
      self.serialRNG = newSerialRNG(sGen, sSeed)
      self.parallelRNG = self.l.newParallelRNG(pGen, pSeed)

  case self.fields["gauge"].start:
    of "cold": unit(self.u)
    of "hot": self.parallelRNG.warm(self.u)
    else: discard

proc newAbstractFieldTheory[L:static[int],S,T,U,V,W](
    l: Layout[L];
    gauge: AbstractField;
    fieldTheory: FieldTheoryType;
    smearing: SmearingType;
    smearingCoeffs: seq[float];
    matterBoundaryConditions: string;
    algorithm: AlgorithmType;
    sGen, pGen: string;
    sSeed, pSeed: uint64;
    s: typedesc[S];
    t: typedesc[T];
    u: typedesc[U];
    v: typedesc[V];
    w: typedesc[W]
  ): AbstractFieldTheory[L,S,T,U,V,W] =
  result = AbstractFieldTheory[L,S,T,U,V,W](
    fieldTheory: fieldTheory,
    smearing: smearing,
    algorithm: algorithm
  )
  result.new(
    l, gauge, smearingCoeffs,
    matterBoundaryConditions,
    sGen, pGen, sSeed, pSeed
  )

template newAbstractFieldTheory*(info: JsonNode; gaugeConstructor: untyped): auto =
  var
    gauge = gaugeConstructor
    fieldTheory: FieldTheoryType
    smearing: SmearingType
    monteCarloAlgorithm: AlgorithmType

  if not info.hasKey("theory"): info["theory"] = %* "gauge"
  case info["theory"].getStr():
    of "yang-mills", "Yang-Mills", "gauge", "Gauge", "pure":
      fieldTheory = PureGauge
    of "gauge-fermion", "Gauge-Fermion":
      if not info.hasKey("matter-discretization"):
        info["matter-discretization"] = %* "staggered"
      case info["matter-discretization"].getStr():
        of "staggered", "Staggered": fieldTheory = StaggeredGaugeFermion
        of "wilson", "Wilson": fieldTheory = WilsonGaugeFermion
        else: discard
    else: discard

  if not info.hasKey("smearing"): info["smearing"] = %* "none"
  case info["smearing"].getStr():
    of "nHYP", "nhyp": smearing = nHYP
    of "stout", "Stout": smearing = Stout
    else: smearing = None

  if not info.hasKey("smearing-coefficients"):
    case smearing:
      of nHYP: info["smearing-coefficients"] = %* @[0.4, 0.5, 0.5]
      of Stout: info["smearing-coefficients"] = %* @[0.1]
      else: info["smearing-coefficients"] = %* @[]
  var smearingCoefficients = newSeq[float]()
  for coeff in info["smearing-coefficients"].getElems():
    smearingCoefficients.add coeff.getFloat()

  if not info.hasKey("matter-boundary-conditions"):
    var bcs = ""
    for mu in 0..<gauge.u[0].l.nDim:
      if mu != gauge.u[0].l.nDim-1: bcs = bcs & "p"
      else: bcs = bcs & "a"
    info["matter-boundary-conditions"] = %* bcs
  let matterBoundaryConditions = info["matter-boundary-conditions"].getStr()

  if not info.hasKey("monte-carlo-algorithm"):
    monteCarloAlgorithm = HamiltonianMonteCarlo
  else:
    case info["monte-carlo-algorithm"].getStr():
      of "hamiltonian-monte-carlo", "Hamiltonian-Monte-Carlo", "hmc", "HMC":
        monteCarloAlgorithm = HamiltonianMonteCarlo
  if not info.hasKey("serial-random-number-generator"):
    info["serial-random-number-generator"] = %* "MILC"
  if not info.hasKey("parallel-random-number-generator"):
    info["parallel-random-number-generator"] = %* "MILC"
  if not info.hasKey("serial-random-number-seed"):
    info["serial-random-number-seed"] = %* 987654321
  if not info.hasKey("parallel-random-number-seed"):
    info["parallel-random-number-seed"] = %* 987654321

  template st: untyped =
    type(ColorMatrix(gauge.u[0].l))
  template tt: untyped =
    type(ColorVector(gauge.u[0].l))
  template ut: untyped =
    type(DiracFermion(gauge.u[0].l))
  template vt: untyped =
    type(ColorVector(gauge.u[0].l)[0])
  template wt: untyped =
    type(spproj1p(DiracFermion(gauge.u[0].l)[0]))

  let l = gauge.u[0].l
  l.newAbstractFieldTheory(
    gauge,
    fieldTheory,
    smearing,
    smearingCoefficients,
    matterBoundaryConditions,
    monteCarloAlgorithm,
    info["serial-random-number-generator"].getStr(),
    info["parallel-random-number-generator"].getStr(),
    uint64(info["serial-random-number-seed"].getInt()),
    uint64(info["parallel-random-number-seed"].getInt()),
    st, tt, ut, vt, wt
  )

proc addFermion*(
    self: var AbstractFieldTheory;
    fermionParams: JsonNode;
    tag: string = "DefaultTag"
  ) =
  case self.fieldTheory:
    of StaggeredGaugeFermion:
      var ferm = self.l.newMatterField(StaggeredFermion, fermionParams)
      if tag == "DefaultTag": self.fields["fermion" & $(self.nFerm)] = ferm
      else: self.fields[tag] = ferm
    of WilsonGaugeFermion:
      var ferm = self.l.newMatterField(WilsonFermion, fermionParams)
      if tag == "DefaultTag": self.fields["fermion" & $(self.nFerm)] = ferm
      else: self.fields[tag] = ferm
    else: qexError "cannot add fermions to this type of field theory"
  self.nFerm = self.nFerm + 1

proc addBoson*(
    self: var AbstractFieldTheory;
    bosonParams: JsonNode;
    tag: string = "DefaultTag"
  ) =
  case self.fieldTheory:
    of StaggeredGaugeFermion:
      var bos = self.l.newMatterField(StaggeredBoson, bosonParams)
      if tag == "DefaultTag": self.fields["boson" & $(self.nBos)] = bos
      else: self.fields[tag] = bos
    of WilsonGaugeFermion:
      var bos = self.l.newMatterField(WilsonBoson, bosonParams)
      if tag == "DefaultTag": self.fields["boson" & $(self.nBos)] = bos
      else: self.fields[tag] = bos
    else: qexError "cannot add bosons to this type of field theory"
  self.nBos = self.nBos + 1

#[ Hamiltonian Monte Carlo ]#

proc setMatterBoundaryConditions(self: AbstractFieldTheory; u: auto) =
  for mu in 0..<u.len:
    if $self.matterBoundaryConditions[mu] == $"a":
      tfor i, 0..<u[mu].l.nSites:
        if u[mu].l.coords[mu][i] == u[mu].l.physGeom[mu]-1:
          u[mu]{i} *= -1.0

proc rephase(self: AbstractFieldTheory; u: auto) =
  threads:
    self.setMatterBoundaryConditions(u)
    threadBarrier()
    u.stagPhase

proc staggeredSmearedOneLinkForce(self: AbstractFieldTheory; f: auto; u: auto) =
  self.rephase(f)
  threads:
    for mu in 0..<f.len:
      for s in f[mu].odd: f[mu][s] *= -1
  self.smearedForce(f, f)
  threads:
    for mu in 0..<f.len:
      for s in f[mu]:
        var temp {.noinit.}: typeof(f[0][0])
        temp := f[mu][s]*u[mu][s].adj
        projectTAH(f[mu][s], temp)

proc smearForce(self: AbstractFieldTheory; f: auto) =
  case self.fieldTheory:
    of Gauge: discard
    of StaggeredGaugeFermion: self.staggeredSmearedOneLinkForce(f)
    of WilsonGaugeFermion: discard

proc hamiltonian*(self: var AbstractFieldTheory; hif: string): seq[float] =
  result = newSeq[float]()

  # Calculate kinetic contribution
  if hif == "hi": self.parallelRNG.randomTAHGaussian(self.p)
  result.add 0.5*pdagp(self.p)

  # Gauge action
  result.add self.fields["gauge"].action()

  # Smearing
  case self.smearing:
    of nHYP:
      self.hyp.smear(self.u, self.su, self.info)
      self.rephase(self.su)
    of Stout: discard # not yet implemented
    else: self.rephase(self.u)

  # Generate pseudofermions, if necessary
  if hif == "hi": 
    for key, _ in self.fields:
      case self.fields[key].field:
        of StaggeredFermion, StaggeredBoson:
          self.parallelRNG.randomComplexGaussian(self.sPsi)
        else: discard

      case self.fields[key].field:
        of StaggeredFermion:
          Ddag(self.sD.stag, self.fields[key].sPhi, self.sPsi, self.fields[key].mass)
        of StaggeredBoson:
          solve(
            self.fields[key],
            self.sD,
            self.fields[key].sPhi,
            self.sPsi,
            self.fields[key].actionSolverParameters
          )
        else: discard

      case self.fields[key].field:
        of StaggeredFermion, StaggeredBoson: zeroOdd(self.fields[key].sPhi)
        else: discard

  # Calculate contributions to action
  for key, _ in self.fields:
    case self.fields[key].field:
      of Gauge: discard
      of StaggeredFermion, StaggeredBoson:
        result.add self.fields[key].action(self.sD, self.sPsi)
      else: discard # Not yet implemented

  # Restore phase on u if necessary
  if self.smearing == None: self.rephase(self.u)

  # Save hi, hf
  case hif:
    of "hi": self.hi = sum(result)
    of "hf": self.hf = sum(result)
    else: discard

proc constructIntegrator*(
    self: var AbstractFieldTheory;
    tau: float
  ) =
  var
    dtaus = newSeq[float]()
    dtau = 0.0
    totalSteps: int

  self.tau = tau

  self.stepsT = newSeq[float]()
  for key, _ in self.fields:
    self.fields[key].setIntegratorAtom(self.tau)
    self.stepsV[key] = newSeq[float]()

  while true:
    dtaus.setLen(0)
    for key, _ in self.fields: dtaus.add self.fields[key].stepsT[0]
    dtau = min(dtaus)

    self.stepsT.add dtau

    for key, _ in self.fields:
      if abs(dtau - self.fields[key].stepsT[0]) < epsilon(float):
        self.stepsV[key].add self.fields[key].stepsV[0]
        self.fields[key].stepsV.delete(0)
        self.fields[key].stepsT.delete(0)
      else:
        self.stepsV[key].add 0.0
        self.fields[key].stepsT[0] -= dtau

    totalSteps = 0
    for key, _ in self.fields: totalSteps += self.fields[key].stepsV.len
    if totalSteps == 0: break

proc runMolecularDynamics*(self: var AbstractFieldTheory) =
  for index in 0..<self.stepsT.len:
    # Flag for whether or not the fields were smeared
    var fieldsSmeared = false

    # Update gauge field
    if self.stepsT[index] != 0.0: updateU(self.u, self.p, self.stepsT[index])

    # Set momentum update step & check if smearing to be performed
    for key, _ in self.fields:
      self.fields[key].b = self.stepsV[key][index]
      if (key != "gauge") and (self.fields[key].b != 0.0):
        case self.smearing:
          of nHYP:
            self.smearedForce = self.hyp.smearGetForce(self.u, self.su, self.info)
            self.rephase(self.su)
            fieldsSmeared = true
          of Stout: discard # Not implemented
          else: self.rephase(self.u)

    # Calculate force from gauge field & update momentum
    self.fields["gauge"].force(self.f)
    if self.fields["gauge"].b != 0: subtractForce(self.p, self.f)
    zero(self.f)

    # Take care of fermions/bosons
    for key, _ in self.fields:
      # Calculate force from fermions/bosons & update momentum
      case self.fields[key].field:
        of Gauge: discard 
        of StaggeredFermion, StaggeredBoson:
          self.fields[key].force(self.sD, self.sPsi, self.u, self.f)
        of WilsonFermion, WilsonBoson: discard # Not implemented

    # Smear force, if necessary
    if fieldsSmeared:
      case self.smearing:
        of nHYP, Stout: self.staggeredSmearedOneLinkForce(self.f, self.u)
        else: self.rephase(self.u)

    # Update momentum
    subtractForce(self.p, self.f)

    # Zero out force
    zero(self.f)

proc backupU(self: var AbstractFieldTheory) = setU(self.bu, self.u)

proc metropolis(self: var AbstractFieldTheory): bool = 
  if self.serialRNG.uniform <= exp(-(self.hf - self.hi)):
    reunit(self.u)
    result = true
  else:
    setU(self.u, self.bu)
    result = false

if isMainModule:
  qexinit()

  #[ Set lattice & parameters up ]#

  let
    lat = intSeqParam("lat", @[4, 4, 4, 4])
    lo = lat.newLayout(@[1, 1, 1, 1])
    tau = 1.0
    nTestConfigs = 10

    fieldTheoryParams = %* {
      "theory": "gauge-fermion",
      "matter-discretization": "staggered",
      "smearing": "nhyp",
      "smearing-coefficients": @[0.4, 0.5, 0.5],
      "matter-boundary-conditions": "aaaa",
      "algorithm": "hamiltonian-monte-carlo",
      "serial-random-number-generator": "MILC",
      "parallel-random-number-generator": "MILC",
      "serial-random-number-seed": 987654321,
      "parallel-random-number-seed": 987654321
    }

    gaugeParams = %* {
      "action": "Wilson",
      "beta": 6.0,
      "steps": 10,
      "integrator": "MN2",
      "start": "cold"
    }

    ns = 2
    fermionParams = %* {
      "mass": 0.0,
      "steps": 5,
      "integrator": "MN2"
    }

    nPV = 8
    pauliVillarsParams = %* {
      "mass": 0.75,
      "steps": 5,
      "integrator": "MN2"
    }

  # Construct field theory object
  var fieldTheory = newAbstractFieldTheory(fieldTheoryParams): lo.newGaugeField(gaugeParams)

  # Construct integrator
  fieldTheory.constructIntegrator(tau)

  # Pure gauge HMC test
  for config in 0..<nTestConfigs:
    fieldTheory.backupU
    discard fieldTheory.hamiltonian("hi")
    fieldTheory.runMolecularDynamics
    discard fieldTheory.hamiltonian("hf")
    let
      accept = fieldTheory.metropolis
      dH = fieldTheory.hf - fieldTheory.hi
    echo config, ": ", fieldTheory.hi, " ", fieldTheory.hf, " ", dH, " ", accept

  # Add two massless staggered fields
  for _ in 0..<ns: fieldTheory.addFermion(fermionParams)

  # Add eight massive Pauli-Villars fields
  for _ in 0..<nPV: fieldTheory.addBoson(pauliVillarsParams)

  # Construct integrator with these new fields
  fieldTheory.constructIntegrator(tau)

  # Dynamical HMC w/ fermions & PV fields test
  for config in 0..<nTestConfigs:
    fieldTheory.backupU
    discard fieldTheory.hamiltonian("hi")
    fieldTheory.runMolecularDynamics
    discard fieldTheory.hamiltonian("hf")
    let
      accept = fieldTheory.metropolis
      dH = fieldTheory.hf - fieldTheory.hi
    echo config, ": ", fieldTheory.hi, " ", fieldTheory.hf, " ", dH, " ", accept

  # Add one more massive fermion for the heck of it
  fermionParams["mass"] = %* 0.5
  fieldTheory.addFermion(fermionParams)

  # Construct integrator with new (massive) fermion in it
  fieldTheory.constructIntegrator(tau)

  # Dynamical HMC w/ (massless/massive) fermions & PV fields test
  for config in 0..<nTestConfigs:
    fieldTheory.backupU
    discard fieldTheory.hamiltonian("hi")
    fieldTheory.runMolecularDynamics
    discard fieldTheory.hamiltonian("hf")
    let
      accept = fieldTheory.metropolis
      dH = fieldTheory.hf - fieldTheory.hi
    echo config, ": ", fieldTheory.hi, " ", fieldTheory.hf, " ", dH, " ", accept

  qexfinalize()