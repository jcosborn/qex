import qex
import layout
import gauge/[hypsmear, stoutsmear]

import json
import options

import abstractFields
import gaugeFields
import staggeredFields
import wilsonFields
import diracOperator
import typeUtilities
import randomNumberGeneration

export abstractFields

type
  SmearingType* = enum Hypercubic, Stout, NoSmearing
  ActionType* = enum PureGauge, PureMatter

type 
  Smearing[S] = object
    su*: seq[S]
    smearedForce*: proc(f, chain: seq[S])
    case smearing*: SmearingType
      of Hypercubic:
        nhyp*: HypCoefs
        nHYPInfo*: PerfInfo
      of Stout: 
        stout*: StoutSmear[seq[S]]
      of NoSmearing: discard

  AbstractAction*[L:static[int],S,T,U,V,W] = object
    fields*: Table[string, AbstractField[S,T,U]]
    case action*: ActionType
      of PureGauge: discard
      of PureMatter:
        l*: Layout[L]
        nf*, nb*: int
        boundaryConditions*: string
        D*: Table[FieldType, DiracOperator[S,T,U,V,W]]
        case smearing*: SmearingType
          of Hypercubic, Stout:
            smear*: Smearing[S]
          of NoSmearing: discard

proc stagPsi(self: AbstractAction): auto =
  result = self.D[StaggeredMatterField].stagPsi

proc wilsPsi(self: AbstractAction): auto =
  result = self.D[WilsonMatterField].wilsPsi

proc stagD(self: AbstractAction): auto = 
  result = self.D[StaggeredMatterField]

proc wilsD(self: AbstractAction): auto = 
  result = self.D[WilsonMatterField]

proc su(self: AbstractAction): auto = self.smear.su

# Procs for creating actions & adding fields to actions

proc newSmearing[S](
    l: Layout; 
    smearing: SmearingType; 
    coeffs: seq[float];
    s: typedesc[S]
  ): Smearing[S] =
  result = Smearing[S](smearing: smearing)
  case result.smearing:
    of Hypercubic:
      result.nhyp = HypCoefs(
        alpha1: coeffs[0],
        alpha2: coeffs[1],
        alpha3: coeffs[2]
      )
    of Stout: result.stout = l.newStoutSmear(coeffs[0])
    of NoSmearing: discard
  case result.smearing:
    of Hypercubic, Stout: result.su = l.newGauge()
    of NoSmearing: discard

proc newGaugeAction[L:static[int],S,T,U,V,W](
    l: Layout[L];
    info: JsonNode;
    s: typedesc[S];
    t: typedesc[T];
    u: typedesc[U];
    v: typedesc[V];
    w: typedesc[W]
  ): AbstractAction[L,S,T,U,V,W] =
  result = AbstractAction[L,S,T,U,V,W](action: PureGauge)
  result.fields["gauge"] = l.newGaugeField(info)

proc newGaugeAction*(l: Layout; info: JsonNode): auto = 
  result = l.newGaugeAction(info, l.typeS, l.typeT, l.typeU, l.typeV, l.typeW)

proc newMatterAction[L:static[int],S,T,U,V,W](
    l: Layout[L];
    info: JsonNode;
    s: typedesc[S];
    t: typedesc[T];
    u: typedesc[U];
    v: typedesc[V];
    w: typedesc[W]
  ): AbstractAction[L,S,T,U,V,W] =
  
  var smearing = NoSmearing
  case info["smearing"].getStr():
    of "nHYP", "nhyp", "hypercubic", "Hypercubic": smearing = Hypercubic
    of "stout", "Stout": smearing = Stout

  result = AbstractAction[L,S,T,U,V,W](action: PureMatter, smearing: smearing)

  if info.hasKey("boundary-conditions"): 
    result.boundaryConditions = info["boundary-conditions"].getStr()
  else: # Defaults to ppp...pa
    result.boundaryConditions = ""
    for _ in 0..<l.nDim-1: 
      result.boundaryConditions = result.boundaryConditions & "p"
    result.boundaryConditions = result.boundaryConditions & "a"

  var coeffs: seq[float]
  case result.smearing:
    of Hypercubic, Stout:
      if not info.hasKey("smearing-coefficients"): 
        qexError "Must specify smearing coefficient(s)"
      else: 
        coeffs = newSeq[float]()
        for coeff in info["smearing-coefficients"].getElems(): 
          coeffs.add coeff.getFloat()
      result.smear = l.newSmearing(smearing, coeffs, l.typeS)
    of NoSmearing: discard

  result.D = initTable[FieldType, DiracOperator[S,T,U,V,W]]()

  result.nf = 0
  result.nb = 0

  result.l = l

proc newMatterAction*(l: Layout; info: JsonNode): auto = 
  result = l.newMatterAction(info, l.typeS, l.typeT, l.typeU, l.typeV, l.typeW)

proc initUnsmearedDiracOperator*(
    self: var AbstractAction; 
    u: auto; 
    disc: FieldType
  ) =
  self.D[disc] = u.newDiracOperator(disc)

proc initSmearedDiracOperator*(self: var AbstractAction; disc: FieldType) = 
  self.D[disc] = newDiracOperator(self.su, disc)

proc addStaggeredFermion*(self: var AbstractAction; info: JsonNode) = 
  let disc = StaggeredMatterField
  self.fields["staggeredFermion" & $(self.nf)] = self.l.newStaggeredFermion(info)
  case self.smearing:
    of Hypercubic, Stout:
      if not self.D.hasKey(disc): 
        self.D[disc] = newDiracOperator(self.su, disc)
    of NoSmearing: discard
  self.nf = self.nf + 1

proc addStaggeredBoson*(self: var AbstractAction; info: JsonNode) = 
  let disc = StaggeredMatterField
  self.fields["staggeredBoson" & $(self.nb)] = self.l.newStaggeredBoson(info)
  case self.smearing:
    of Hypercubic, Stout:
      if not self.D.hasKey(disc): 
        self.D[disc] = newDiracOperator(self.su, disc)
    of NoSmearing: discard
  self.nb = self.nb + 1

proc addStaggeredHasenbuschFermion*(
    self: var AbstractAction; 
    infos: seq[JsonNode]
  ) = 
  let disc = StaggeredMatterField
  for index in 0..<infos.len-1:
    infos[index]["mass1"] = %* infos[index+1]["mass"].getFloat()
    infos[index]["mass2"] = %* infos[index]["mass"].getFloat()
    let tag = "staggeredFermion" & $(self.nf) & "_hasenbusch" & $(index)
    self.fields[tag] = self.l.newStaggeredHasenbuschFermion(infos[index])
  self.fields["staggeredFermion" & $(self.nf)] = self.l.newStaggeredFermion(infos[^1])
  case self.smearing:
    of Hypercubic, Stout:
      if not self.D.hasKey(disc): 
        self.D[disc] = newDiracOperator(self.su, disc)
    of NoSmearing: discard
  self.nf = self.nf + 1

# Helper procs for action/force calculation

proc rephase(u: auto) =
  threads: u.stagPhase

proc setboundaryConditions(self: AbstractAction; u: auto) =
  let bcs = self.boundaryConditions
  threads:
    for mu in 0..<u.len:
      if $bcs[mu] == $"a":
        tfor i, 0..<u[mu].l.nSites:
          if u[mu].l.coords[mu][i] == u[mu].l.physGeom[mu]-1:
            u[mu]{i} *= -1.0

# Procs for calculating action & forces

proc action*(self: var AbstractAction): float = 
  result = self.fields["gauge"].gaugeAction()

proc action*[S](
    self: var AbstractAction; 
    u: seq[S]; 
    parallelRNG: Option[ParallelRNG] = none(ParallelRNG)
  ): Table[string,float] =
  result = initTable[string,float]()
  case self.action:
    of PureGauge: result["gauge"] = self.fields["gauge"].gaugeAction()
    of PureMatter:
      # Switches for matter fields
      var rephased = false

      # Take care of smearing
      case self.smearing:
        of Hypercubic: 
          self.smear.nhyp.smear(u, self.su, self.smear.nHYPInfo)
          self.setboundaryConditions(self.su)
        of Stout: discard
        of NoSmearing: self.setboundaryConditions(u)

      # Calculate contribution from staggered matter field
      for key, _ in self.fields:
        if (self.fields[key].field == StaggeredMatterField):
          # Rephase
          if not rephased:
            case self.smearing:
              of Hypercubic, Stout: rephase(self.su)
              of NoSmearing: rephase(u)
            rephased = true

          # Optional fermion heatbath
          if not parallelRNG.isNone:
            randomComplexGaussian(parallelRNG.get, self.stagPsi)
            self.fields[key].getStaggeredField(self.stagD, self.stagPsi)

          # Calculate & append contribution for full action
          result[key] = self.fields[key].staggeredAction(self.stagD, self.stagPsi)
      if (self.smearing == NoSmearing) and (rephased): rephase(u)

      #[Wilson matter field - not implemented]#

      # Undo boundary condition rephase if no smearing
      case self.smearing:
        of Hypercubic, Stout: discard
        of NoSmearing: self.setboundaryConditions(u)

proc updateMomentum*[S](self: var AbstractAction; f,p: seq[S]) =
  # Switches for momentum update
  var updateMomentum = false

  # Momentum update from gauge field
  let dtau =  self.fields["gauge"].dtau[0]
  if dtau != 0.0: updateMomentum = true
  if updateMomentum: 
    self.fields["gauge"].gaugeForce(f)
    threads:
      for mu in 0..<f.len:
        for s in f[mu]: p[mu][s] -= dtau*f[mu][s]

proc updateMomentum*[S](self: var AbstractAction; u: seq[S]; f,p: seq[S]) =
  # Switches for momentum update
  var updateMomentum = false

  # Calculate force
  case self.action:
    of PureGauge:
      let dtau =  self.fields["gauge"].dtau[0]
      if dtau != 0.0: updateMomentum = true
      if updateMomentum: 
        self.fields["gauge"].gaugeForce(f)
        threads:
          for mu in 0..<f.len:
            for s in f[mu]: f[mu][s] := dtau*f[mu][s]
    of PureMatter:
      # Switches for matter fields
      var fieldsRephased = false

      # Check if smearing & momentum update need to be done
      for key, _ in self.fields:
        if self.fields[key].dtau[0] != 0.0: updateMomentum = true

      # Calculate force
      if updateMomentum:
        # Zero force
        threads:
          for mu in 0..<f.len: 
            for s in f[mu]: f[mu][s] := 0

        # Smear fields
        case self.smearing:
          of Hypercubic:
            self.smear.smearedForce = smearGetForce(
              self.smear.nhyp, u, self.su, self.smear.nHYPInfo
            )
            self.setboundaryConditions(self.su)
          of Stout: discard
          of NoSmearing: self.setboundaryConditions(u)

        # Take care of contribution to force from staggered fields
        for key, _ in self.fields:
          let dtau = self.fields[key].dtau[0]
          let field = self.fields[key].field
          if (dtau != 0.0) and (field == StaggeredMatterField):
            if not fieldsRephased: # Rephase gauge fields - done only once 
              case self.smearing:
                of Hypercubic: rephase(self.su)
                of Stout: discard
                of NoSmearing: rephase(u)
              fieldsRephased = true
            self.fields[key].staggeredPartialForce(
              self.D[StaggeredMatterField], self.stagPsi, f
            )
        threads: # Rephase force
          f.stagPhase
          threadBarrier()
          for mu in 0..<f.len:
            for s in f[mu].odd: f[mu][s] *= -1
        case self.smearing: # Undo rephase if no smearing
          of Hypercubic, Stout: discard
          of NoSmearing:
            if fieldsRephased: rephase(u)

        #[Update from Wilson fields - not implemeted, will need to undo su rephase]#

        # Set boundary conditions, smear, and project
        self.setboundaryConditions(f)
        case self.smearing: 
          of Hypercubic:
            self.smear.smearedForce(f,f)
            threads:
              for mu in 0..<f.len:
                for s in f[mu]:
                  var temp {.noinit.}: typeof(f[0][0])
                  temp := f[mu][s]*u[mu][s].adj
                  projectTAH(f[mu][s], temp)
          of Stout: discard
          of NoSmearing: self.setboundaryConditions(u) # Undoes setting matter bc

  #[Optional force-gradient - not implemented]#
        
  # Update momentum
  if updateMomentum:
    threads:
      for mu in 0..<p.len:
        for s in p[mu]: p[mu][s] -= f[mu][s]

if isMainModule:
  qexInit()

  var 
    lat = intSeqParam("lat", @[4, 4, 4, 4])
    lo = lat.newLayout(@[1, 1, 1, 1])

    gaugeParams = %* {
      "action": "Wilson",
      "beta": 6.0,
      "steps": 10,
      "integrator": "MN2",
      "monte-carlo-algorithm": "hamiltonian-monte-carlo"
    }
    matterActionParams = %* {
      "smearing": "nHYP",
      "smearing-coefficients": @[0.4, 0.5, 0.5],
      "boundary-conditions": "aaaa"
    }
    fermionParams = %* {
      "mass": 0.0,
      "steps": 5,
      "integrator": "MN2",
      "monte-carlo-algorithm": "hamiltonian-monte-carlo"
    }
    bosonParams = %* {
      "mass": 0.1,
      "steps": 5,
      "integrator": "MN2",
      "monte-carlo-algorithm": "hamiltonian-monte-carlo"
    }

    rng = lo.newParallelRNG("MILC", 987654321)
    u = lo.newGauge()
    f = lo.newGauge()
    p = lo.newGauge()

  var gauge = lo.newGaugeAction(gaugeParams)
  var matter = lo.newMatterAction(matterActionParams)

  matter.addStaggeredFermion(fermionParams)
  matter.addStaggeredBoson(bosonParams)
  matter.addStaggeredHasenbuschFermion(@[fermionParams, bosonParams, bosonParams])

  unit(u)
  unit(p)

  echo gauge.action(u)
  echo gauge.action()
  gauge.fields["gauge"].dtau = @[1.0]
  gauge.updateMomentum(u, f, p)

  echo matter.action(u, parallelRNG = some(rng))
  echo matter.action(u)
  for key, _ in matter.fields: matter.fields[key].dtau = @[1.0]
  matter.updateMomentum(u, f, p)

  qexFinalize()