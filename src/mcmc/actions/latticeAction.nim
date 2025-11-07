
import ../mcmcTypes
import ../fields/gaugeFields
import ../fields/staggeredFields

proc checkJSON(
    info: JsonNode;
    actionType: ActionType;
    algorithmType: MonteCarloAlgorithm  
  ): JsonNode =
  let action = case actionType:
    of PureGauge: "gauge"
    of PureMatter: "matter"
    of GaugeMatter: "gauge-matter"
  let algorithm = case algorithmType:
    of HamiltonianMonteCarlo: "hmc"
    of HeatbathOverrelax: "hb"
  result = parseJson("{}")
  for k,kv in info: result[k] = kv
  result["action-type"] = %* action
  result["monte-carlo-algorithm"] = %* algorithm

template diracOperator(self: LatticeAction): untyped = self.D

template diracOperator(self: LatticeSubAction): untyped = self.D[]

template su(self: LatticeAction): untyped = self.smear.su

template su(self: LatticeSubAction): untyped = self.smear[].su

proc initDiracOperator*(self: var LatticeAction; disc: FieldType) = 
  self.D[disc] = self.l[].newDiracOperator(self.smear.su, disc)

proc initDiracOperator*(self: var LatticeSubAction; disc: FieldType) = 
  self.D[][disc] = self.l[].newDiracOperator(self.smear[].su, disc)

proc newStaggeredField[A](
    self: A; 
    sat: StaggeredActionType; 
    solo: bool;
    matterInfo: JsonNode
  ): auto =
  result = self.newLatticeSubAction(matterInfo, solo)
  result.pField = case sat
    of StaggeredFermion: self.l[].newStaggeredFermion(matterInfo)
    of RootedStaggeredFermion: 
      self.l[].newRootedStaggeredFermion(matterInfo)
    of StaggeredBoson: self.l[].newStaggeredBoson(matterInfo)
    of StaggeredHasenbuschFermion: 
      self.l[].newStaggeredHasenbuschFermion(matterInfo)
  if not hasKey(result.diracOperator, StaggeredMatterField):
    result.initDiracOperator(StaggeredMatterField)

template newStaggeredField[A](
    self: var A; 
    sat: StaggeredActionType;
    matterInfo: JsonNode;
    construction: untyped
  ): auto = 
  block:
    var subAction {.inject.} = self.newStaggeredField(sat, false, matterInfo)
    construction
    subAction

proc addStaggeredFermion*[A](self: var A; info: JsonNode) =
  let sat = StaggeredFermion
  self.subActions.add self.newStaggeredField(sat,true,info)

proc addRootedStaggeredFermion*[A](self: var A; info: JsonNode) =
  let sat = RootedStaggeredFermion
  self.subActions.add self.newStaggeredField(sat,true,info)

proc addStaggeredBoson*[A](self: var A; info: JsonNode) =
  let sat = StaggeredBoson
  self.subActions.add self.newStaggeredField(sat,true,info)

proc addStaggeredHasenbuschFermion*[A](self: var A; info: JsonNode) =
  let sat = StaggeredHasenbuschFermion
  self.subActions.add self.newStaggeredField(sat,true,info)

template addStaggeredFermion*[A](
    self: var A;
    info: JsonNode;
    construction: untyped
  ) =
  let sat = StaggeredFermion
  self.subActions.add self.newStaggeredField(sat, info, construction)

template addRootedStaggeredFermion*[A](
    self: var A;
    info: JsonNode;
    construction: untyped
  ) =
  let sat = RootedStaggeredFermion
  self.subActions.add self.newStaggeredField(sat, info, construction)

template addStaggeredBoson*[A](
    self: var A;
    info: JsonNode;
    construction: untyped
  ) =
  let sat = StaggeredBoson
  self.subActions.add self.newStaggeredField(sat, info, construction)

template addStaggeredHasenbuschFermion*[A](
    self: var A;
    info: JsonNode;
    construction: untyped
  ) =
  let sat = StaggeredHasenbuschFermion
  self.subActions.add self.newStaggeredField(sat, info, construction)

proc addGaugeField*[A](self: var A; gaugeInfo: JsonNode) =
  var action = self.newLatticeSubAction(gaugeInfo, true)
  action.pField = self.l[].newGaugeField(gaugeInfo)
  self.subActions.add action

template addGaugeField*[A](
    self: var A; 
    gaugeInfo: JsonNode; 
    construction: untyped
  ) = 
  block:
    var subAction {.inject.} = self.newLatticeSubAction(gaugeInfo, false)
    subAction.pField = self.l[].newGaugeField(gaugeInfo)
    construction
    self.subActions.add subAction

template newLatticeAction(
    l: Layout;
    actionInfo: JsonNode;
    actionType: ActionType;
    algorithmType: MonteCarloAlgorithm;
    construction: untyped
  ): auto =
  block:
    var
      info = checkJSON(actionInfo, actionType, algorithmType)
      action {.inject.} = l.newLatticeAction(info)
    construction
    action

proc addGaugeAction*(
    self: var LatticeFieldTheory;
    gInfo: JsonNode
  ) = 
  let
    alg = self.algorithm
    act = PureGauge
    info = parseJson("{}")
  for k,kv in gInfo: info[k] = kv
  info["smearing"] = %* "none"
  let gaugeAction = self.l.newLatticeAction(info, act, alg): 
      action.addGaugeField(info)
  self.actions.add gaugeAction

template addMatterAction*(
    self: var LatticeFieldTheory;
    info: JsonNode;
    construction: untyped
  ) = 
  let
    alg = self.algorithm
    act = PureMatter
  self.actions.add self.l.newLatticeAction(info, act, alg, construction)

template addGaugeMatterAction*(
    self: var LatticeFieldTheory;
    info: JsonNode;
    construction: untyped
  ) = 
  let
    alg = self.algorithm
    act = GaugeMatter
  self.actions.add self.l.newLatticeAction(info, act, alg, construction)

if isMainModule:
  qexInit()

  let
    ftInfo = %* {
      "lattice-geometry": @[4,4,4,4],
      "mpi-geometry": @[1,1,1,1],
      "monte-carlo-algorithm": "hmc",
      "trajectory-length": 1.0,
      "serial-random-number-seed": 987654321,
      "parallel-random-number-seed": 987654321,
      "serial-random-number-generator": "milc",
      "parallel-random-number-generator": "milc"
    }
    aInfo = %* {
      "smearing": "nhyp",
      "smearing-coefficients": @[0.4,0.5,0.5],
      "boundary-conditions": "pppa"
    }
    gInfo = %* {
      "action": "Wilson",
      "beta": 6.0,
      "steps": 10,
      "integrator": "MN2"
    }
    fInfo = %* {
      "mass": 0.0,
      "integrator": "2MN",
      "steps": 10
    }
    bInfo = %* {
      "mass": 0.75,
      "integrator": "2MN",
      "steps": 10
    }
    hInfo = %* {
      "mass1": 0.2,
      "mass2": 0.3,
      "integrator": "2MN",
      "steps": 10
    }

  var fieldTheory = newLatticeFieldTheory(ftInfo):
    fieldTheory.addGaugeAction(gInfo)
    fieldTheory.addMatterAction(aInfo):
      action.addStaggeredFermion(fInfo)
      action.addStaggeredBoson(bInfo)
      action.addStaggeredHasenbuschFermion(hInfo)
      action.addGaugeField(gInfo)
    fieldTheory.addGaugeMatterAction(aInfo):
      action.addStaggeredBoson(bInfo):
        subAction.addStaggeredFermion(fInfo)
      action.addStaggeredHasenbuschFermion(hInfo):
        subAction.addStaggeredFermion(fInfo)
      action.addGaugeField(gInfo):
        subAction.addStaggeredFermion(fInfo)
      action.addStaggeredFermion(fInfo):
        subAction.addStaggeredFermion(fInfo):
          subAction.addStaggeredFermion(fInfo)
        subAction.addStaggeredBoson(bInfo):
          subAction.addStaggeredFermion(fInfo)
        subAction.addStaggeredHasenbuschFermion(hInfo):
          subAction.addStaggeredFermion(fInfo)
        subAction.addGaugeField(gInfo):
          subAction.addStaggeredFermion(fInfo)
  

  qexFinalize()