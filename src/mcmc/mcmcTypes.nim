import qex
import layout
import gauge/[hypsmear, stoutsmear]
import physics/[stagD, wilsonD]

import strutils
import json
import options

import utilities/stream

export qex
export json
export stream

const banner = """

|---------------------------------------------------------------|
 Quantum EXpressions (QEX) Markov chain Monte Carlo

 QEX authors: James Osborn & Xiao-Yong Jin
 Monte Carlo authors: 
   - Curtis Taylor Peterson [C.T.P.] (Michigan State University)
   - James Osborn (Argonne National Laboratory)
   - Xiao-Yong Jin (Argonne National Laboratory) 
   - Anna Hasenfratz (University of Colorado Boulder)
 QEX GitHub: https://github.com/jcosborn/qex
 QEX Monte Carlo GitHub: https://github.com/ctpeterson/qex
 C.T.P. email: curtistaylorpetersonwork@gmail.com
 cite: Proceedings of Science (PoS) LATTICE2016 (2017) 271
|---------------------------------------------------------------|
"""

let 
  Large64* = 1.0/epsilon(float)
  Small64* = epsilon(float)

type
  MonteCarloAlgorithm* = enum 
    HamiltonianMonteCarlo, 
    HeatbathOverrelax
  RandomNumberGeneratorType* = enum 
    MILC, 
    MRG
  MolecularDynamicsUpdateType* = enum
    UpdateV,
    UpdateT,
    UpdateVTV
  Integrator = enum
    LeapFrog,
    Om2MN,
    Om4MN4FP,
    Om4MN5FV,
    Om4MN5FP,
    CustomIntegrator,
    UpdaterIntegrator
  GaugeStart* = enum 
    UnitGauge, 
    RandomGauge,
    ReadGauge
  SmearingType* = enum 
    Hypercubic, 
    Stout, 
    NoSmearing
  ActionType* = enum 
    PureGauge, 
    PureMatter, 
    GaugeMatter
  FieldType* = enum
    GaugeField,
    StaggeredMatterField,
    WilsonMatterField,
    DummyField
  GaugeActionType* = enum 
    Wilson, 
    Adjoint, 
    Rectangle, 
    Symanzik, 
    Iwasaki, 
    DBW2
  StaggeredActionType* = enum
    StaggeredFermion,
    RootedStaggeredFermion,
    StaggeredHasenbuschFermion,
    StaggeredBoson
  WilsonActionType* = enum
    WilsonFermion,
    WilsonBoson

type
  MolecularDynamicsUpdate = object
    scale*,dtau*: float
    case update*: MolecularDynamicsUpdateType
      of UpdateV: discard
      of UpdateT:
        absdtau*: float
      of UpdateVTV: discard

  MonteCarloAtom[S] {.inheritable.} = object
    case algorithm*: MonteCarloAlgorithm
      of HamiltonianMonteCarlo:
          timeStep*,spaceStep*,iStep*,steps*: int
          running*: bool
          dtau*,vdtau*,absScale*,cumSum*,absCumSum*: float
          updates*: seq[MolecularDynamicsUpdate]
          case parallel: bool
            of true:
              tdtau*: float
              p*,f*,u*: ref seq[S]
              nested*,includeInStep*: bool
            of false: discard
          case integrator*: Integrator
            of LeapFrog: discard
            of Om2MN,Om4MN4FP,Om4MN5FV,Om4MN5FP:
              lmbda,theta,rho,vartheta,xi,mu: float
            of CustomIntegrator: discard
            of UpdaterIntegrator: discard
          pRNG*: ref ParallelRNG
          currentAction*: float
      of HeatbathOverrelax: discard

  MonteCarlo[S,V] {.inheritable.} = object
    start*: GaugeStart
    case algorithm*: MonteCarloAlgorithm
      of HamiltonianMonteCarlo:
        tau*,tdtau*: float
        bu*,p*,f*: seq[S]
        u*: ref seq[S]
        sRNG*: SerialRNG
        pRNG*: ParallelRNG
        nFields*: int
        nested*: bool
      of HeatbathOverrelax:
        sf: seq[Shifter[S,V]]

  RandomNumberGenerator {.inheritable.} = object
    generator*: RandomNumberGeneratorType
    seed: uint64

  ParallelRNG* = object of RandomNumberGenerator
    milc*: typeof(Field[1, RngMilc6])
    mrg*: typeof(Field[1, MRG32k3a])

  SerialRNG* = object of RandomNumberGenerator
    milc*: RngMilc6
    mrg*: MRG32k3a

  RemezCoefficients* = object
    nTerms*: int
    f0*,if0*: float
    alpha*,ialpha*: seq[float]
    beta*,ibeta*: seq[float]

  LatticeField*[S,T,U] = object
    id: string
    case field*: FieldType
      of GaugeField:
        gaugeAction*: GaugeActionType
        gaugeActionCoefficients*: GaugeActionCoeffs
        u*: seq[S]
      of StaggeredMatterField:
        staggeredFields*: seq[T]
        staggeredMasses*: seq[float]
        stagActionSolverParams*: SolverParams
        stagForceSolverParams*: SolverParams
        case staggeredAction*: StaggeredActionType
          of RootedStaggeredFermion:
            remez*: RemezCoefficients
            rPhi*: T
          else: discard
      of WilsonMatterField:
        wilsonFields*: seq[U] 
        wilsonMasses*: seq[float] 
        wilsonAction: WilsonActionType
        wilsubActionsolverParams: SolverParams
        wilsForceSolverParams: SolverParams
      of DummyField:
        du*: seq[S]
        ds*: T

  Smearing*[S] = ref object
    su*: seq[S]
    u*: ref seq[S]
    sf*: proc(f, chain: seq[S])
    rephased*: bool
    smeared*: bool
    bc*: string
    case smearing*: SmearingType
      of Hypercubic:
        nhyp*: HypCoefs
        nhypInfo*: PerfInfo
      of Stout: 
        stout*: StoutSmear[seq[S]]
      of NoSmearing: discard

  DiracOperator*[S,T,U,W,X] = object
    case discretization*: FieldType
      of StaggeredMatterField:
        stag*: Staggered[S,W]
        stagShifter*: seq[Shifter[T,W]]
        stagPsi*: T
      of WilsonMatterField:
        wils*: Wilson[S,W]
        wilsShifter*: seq[Shifter[U,X]]
        wilsPsi*: U
      else: discard

  LatticeSubAction*[L:static[int],S,T,U,W,X] = ref object of MonteCarloAtom[S]
    l*: ref Layout[L]
    pField*: LatticeField[S,T,U]
    case solo*: bool
      of false:
        subActions*: seq[LatticeSubAction[L,S,T,U,W,X]]
      of true: discard
    case action*: ActionType:
      of PureMatter,GaugeMatter:
        D*: ref OrderedTable[FieldType, DiracOperator[S,T,U,W,X]]
      of PureGauge: discard
    smearing*: SmearingType
    smear*: ref Smearing[S]

  LatticeAction*[L:static[int],S,T,U,W,X] = object
    l*: ref Layout[L]
    subActions*: seq[LatticeSubAction[L,S,T,U,W,X]]
    algorithm*: MonteCarloAlgorithm
    case action*: ActionType
      of PureGauge: discard
      of PureMatter,GaugeMatter:
        nf*, nb*: int
        bc*: string
        D*: OrderedTable[FieldType, DiracOperator[S,T,U,W,X]]
        smearing*: SmearingType
        smear*: Smearing[S]

  LatticeFieldTheory*[L:static[int],S,T,U,V,W,X] = object of MonteCarlo[S,V]
    l*: Layout[L]
    actions*: seq[LatticeAction[L,S,T,U,W,X]]

converter toRNGType*(s: string):
  RandomNumberGeneratorType = parseEnum[RandomNumberGeneratorType](s)

template SS*(l: Layout): untyped =
  type(l.ColorMatrix())

template TT*(l: Layout): untyped = 
  type(l.ColorVector())

template UU*(l: Layout): untyped = 
  type(l.DiracFermion())

template VV*(l: Layout): untyped =
  type(l.ColorMatrix()[0])

template WW*(l: Layout): untyped =
  type(l.ColorVector()[0])

template XX*(l: Layout): untyped =
  type(spproj1p(l.DiracFermion()[0]))

proc stagPsi*(self: LatticeAction): auto =
  result = self.D[StaggeredMatterField].stagPsi

proc wilsPsi*(self: LatticeAction): auto =
  result = self.D[WilsonMatterField].wilsPsi

proc stagD*(self: LatticeAction): auto = 
  result = self.D[StaggeredMatterField]

proc wilsD*(self: LatticeAction): auto = 
  result = self.D[WilsonMatterField]

proc new(self: var RandomNumberGenerator; generator: string; seed: uint64) =
  case generator:
    of "MILC", "RngMilc6": self.generator = MILC 
    of "MRG", "MRG32k3a": self.generator = MRG
    else: qexError generator, " not supported"
  self.seed = seed

proc newParallelRNG*(
    l: Layout;
    generator: string;
    seed: uint64;
  ): ParallelRNG =
  new(result, generator, seed)
  case result.generator:
    of MILC: result.milc = l.newRNGField(RngMilc6, result.seed)
    of MRG: result.mrg = l.newRNGField(MRG32k3a, result.seed)

proc seed*(self: var SerialRNG) =
  case self.generator:
    of MILC: self.milc.seed(self.seed, 987654321)
    of MRG: self.mrg.seed(self.seed, 987654321)

proc newSerialRNG*(generator: string; seed: uint64): SerialRng =
  new(result, generator, seed)
  result.seed()

proc mkLeapFrog(self: LatticeSubAction): auto = 
  result = @[
    (UpdateT,0.5),
    (UpdateV,1.0),
    (UpdateT,0.5),
    (UpdateV,0.0)
  ]

proc mkOm2MN(self: LatticeSubAction): auto =
  result = @[
    (UpdateT, self.lmbda),
    (UpdateV, 0.5),
    (UpdateT, 1.0 - 2.0*self.lmbda),
    (UpdateV, 0.5),
    (UpdateT, self.lmbda),
    (UpdateV, 0.0)
  ]

proc mkOm4MN4FP(self: LatticeSubAction): auto =
  result = @[
    (UpdateT, self.rho),
    (UpdateV, self.lmbda),
    (UpdateT, self.theta),
    (UpdateV, 0.5 - self.lmbda),
    (UpdateT, 1.0 - 2.0*self.theta - 2.0*self.rho),
    (UpdateV, 0.5 - self.lmbda),
    (UpdateT, self.theta),
    (UpdateV, self.lmbda),
    (UpdateT, self.rho),
    (UpdateV, 0.0)
  ]

proc mkOm4MN5FV(self: LatticeSubAction): auto =
  result = @[
    (UpdateT, 0.0),
    (UpdateV, self.vartheta),
    (UpdateT, self.rho),
    (UpdateV, self.lmbda),
    (UpdateT, self.theta),
    (UpdateV, 0.5 - self.lmbda - self.vartheta),
    (UpdateT, 1.0 - 2.0*self.theta - 2.0*self.rho),
    (UpdateV, 0.5 - self.lmbda - self.vartheta),
    (UpdateT, self.theta),
    (UpdateV, self.lmbda),
    (UpdateT, self.rho),
    (UpdateV, self.vartheta)
  ]

proc mkOm4MN5FP(self: LatticeSubAction): auto =
  result = @[
    (UpdateT, self.rho),
    (UpdateV, self.vartheta),
    (UpdateT, self.theta),
    (UpdateV, self.lmbda),
    (UpdateT, 0.5 - self.theta - self.rho),
    (UpdateV, 1.0 - 2.0*self.lmbda - 2.0*self.vartheta),
    (UpdateT, 0.5 - self.theta - self.rho),
    (UpdateV, self.lmbda),
    (UpdateT, self.theta),
    (UpdateV, self.vartheta),
    (UpdateT, self.rho),
    (UpdateV, 0.0)
  ]

proc newLatticeField[S,T,U](
    info: JsonNode;
    s: typedesc[S];
    t: typedesc[T];
    u: typedesc[U]
  ): LatticeField[S,T,U] =
  let
    fieldInfo = case info.hasKey("field-type")
      of true: info["field-type"].getStr()
      of false: ""
    id = case info.hasKey("id")
      of true: info["id"].getStr()
      of false: ""
  var field: FieldType

  case fieldInfo:
    of "gauge": field = GaugeField
    of "staggered": field = StaggeredMatterField
    of "wilson", "Wilson": field = WilsonMatterField
    of "dummy": field = DummyField
    of "": qexError "Must specify field type for " & id
    else: qexError fieldInfo & " is not a valid field type for " & id

  result = LatticeField[S,T,U](field: field)
  result.id = id

proc newLatticeField*(l: Layout; info: JsonNode): auto = 
  result = info.newLatticeField(l.SS, l.TT, l.UU)

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
  result.rephased = false
  result.smeared = false
  new(result.u)

proc newDiracOperator[S,T,U,W,X](
    l: Layout;
    g: auto;
    discretization: FieldType;
    s: typedesc[S];
    t: typedesc[T];
    u: typedesc[U];
    w: typedesc[W];
    x: typedesc[X]
  ): DiracOperator[S,T,U,W,X] =
  result = DiracOperator[S,T,U,W,X](discretization: discretization)
  case result.discretization:
    of StaggeredMatterField:
      result.stagPsi = l.ColorVector()
      result.stag = newStag(g)
      result.stagShifter = newSeq[Shifter[T,W]](g[0].l.nDim)
      for mu in 0..<g.len: 
        result.stagShifter[mu] = newShifter(result.stagPsi, mu, 1)
    of WilsonMatterField: 
      result.wilsPsi = l.DiracFermion()
      #[
      result.wils = newWilson(g)
      result.wilsShifter = newSeq[Shifter[U,W]](g[0].l.nDim)
      for mu in 0..<g.len: 
        result.wilsShifter[mu] = newShifter(result.wilsPsi, mu, 1)
      ]#
    else: discard

proc newDiracOperator*[S](l: Layout; u: seq[S]; disc: FieldType): auto = 
  result = l.newDiracOperator(u, disc, S, l.TT, l.UU, l.WW, l.XX)

proc checkJSON[A](
    self: A;
    info: JsonNode
  ): JsonNode =
  var 
    algorithm = case self.algorithm
      of HamiltonianMonteCarlo: "hmc"
      of HeatbathOverrelax: "hb"
    action = case self.action:
      of PureGauge: "gauge"
      of PureMatter: "matter"
      of GaugeMatter: "gauge-matter"
    smearing = case self.smearing:
      of Hypercubic: "nhyp"
      of Stout: "stout"
      of NoSmearing: "none"
  result = parseJson("{}")
  for k,kv in info: result[k] = kv
  result["monte-carlo-algorithm"] = %* algorithm
  result["action-type"] = %* action
  result["smearing"] = %* smearing

proc setSubActionDiracOperator[S,T,U,W,X](
    self: var LatticeSubAction;
    D: OrderedTable[FieldType, DiracOperator[S,T,U,W,X]]
  ) =
  self.D[] = D

proc setSubActionDiracOperator[S,T,U,W,X](
    self: var LatticeSubAction;
    D: ref OrderedTable[FieldType, DiracOperator[S,T,U,W,X]]
  ) =
  self.D[] = D[]

proc setSubActionSmearing[S](self: var LatticeSubAction; smear: Smearing[S]) =
  self.smear[] = smear

proc setSubActionSmearing[S](self: var LatticeSubAction; smear: ref Smearing[S]) =
  self.smear[] = smear[]

proc setSubActionDiracOperator[A](self: var LatticeSubAction; parent: A) = 
  case parent.action:
    of PureMatter,GaugeMatter: 
      self.setSubActionDiracOperator(parent.D)
    of PureGauge: discard  

proc setSubActionSmearing[A](self: var LatticeSubAction; parent: A) = 
  case parent.action:
    of PureMatter,GaugeMatter: 
      self.setSubActionSmearing(parent.smear)
    of PureGauge: discard

proc newLatticeSubAction[L:static[int],S,T,U,W,X](
    l: Layout[L];
    stream: var MCStream;
    info: JsonNode;
    solo: bool;
    s: typedesc[S];
    t: typedesc[T];
    u: typedesc[U];
    w: typedesc[W];
    x: typedesc[X]
  ): LatticeSubAction[L,S,T,U,W,X] = 
  let
    algInfo = case info.hasKey("monte-carlo-algorithm")
      of true: info["monte-carlo-algorithm"].getStr()
      of false: "hmc"
  var 
    alg: MonteCarloAlgorithm
    action: ActionType
    smearing: SmearingType

  case algInfo
    of "hmc", "HMC": alg = HamiltonianMonteCarlo
    of "hb": qexError "QEX MCMC does not yet support heatbath + overrelax"
    else: qexError "QEX MCMC does not support " & algInfo
    
  case info.hasKey("action-type"):
    of true:
      case info["action-type"].getStr():
        of "gauge": action = PureGauge
        of "matter": action = PureMatter
        of "gauge-matter": action = GaugeMatter
    else: qexError "Action type not specified for sub-action"

  case info.hasKey("smearing"):
    of true:
      case info["smearing"].getStr():
        of "nhyp": smearing = Hypercubic
        of "stout": smearing = Stout
        of "none": smearing = NoSmearing
    of false: qexError "Action smearing not specified for sub-action"

  result = LatticeSubAction[L,S,T,U,W,X](
    algorithm: alg, 
    solo: solo, 
    action: action, 
    smearing: smearing
  )
  new(result.l)
  result.l[] = l
  new(result.pRNG)
  new(result.smear)

  case result.action:
    of PureGauge: discard
    of PureMatter,GaugeMatter: new(result.D)

  if not solo:
    result.subActions = newSeq[LatticeSubAction[L,S,T,U,W,X]]()

  case alg:
    of HamiltonianMonteCarlo:
      var steps: int
      var updates: seq[tuple[first:MolecularDynamicsUpdateType,second:float]]
      result.parallel = not solo
      case info.hasKey("steps"):
        of true: result.steps = info["steps"].getInt()
        of false: qexError "Must specify number of integration steps"
      case info.hasKey("integrator"):
        of true: 
          stream.add "  integrator = " & info["integrator"].getStr()
          stream.add "  integrator steps = " & $(info["steps"].getInt())
          case info["integrator"].getStr()
            of "lf", "LF", "leapfrog", "leap-frog", "LeapFrog":
              result.integrator = LeapFrog
            of "sexton-weingarten", "SextonWeingarten":
              result.integrator = Om2MN
              result.lmbda = 1.0/6.0
            of "MinimumNorm2", "2mn", "2MN":
              result.integrator = Om2MN
              result.lmbda = 0.1931833275037836
            of "4mn4fp", "4MN4FP":
              result.integrator = Om4MN4FP
              result.rho = 0.1786178958448091
              result.theta = -0.06626458266981843
              result.lmbda = 0.7123418310626056
            of "4mn5fv", "4MN5FV":
              result.integrator = Om4MN5FV
              result.rho = 0.2539785108410595
              result.theta = -0.03230286765269967
              result.vartheta = 0.08398315262876693
              result.lmbda = 0.6822365335719091
            of "4mn5fp", "4MN5FP":
              result.integrator = Om4MN5FP
              result.rho = 0.2750081212332419
              result.theta = -0.1347950099106792
              result.vartheta = -0.08442961950707149
              result.lmbda = 0.3549000571574260
            else:
              var msg = info["integrator"].getStr() 
              msg = msg & " is not a valid integrator" 
              qexError msg
          for key in info.keys():
            case key:
              of "lambda": result.lmbda = info[key].getFloat()
              of "rho": result.rho = info[key].getFloat()
              of "theta": result.theta = info[key].getFloat()
              of "vartheta": result.vartheta = info[key].getFloat()
              of "xi": result.xi = info[key].getFloat()
              of "mu": result.mu = info[key].getFloat()
              else: discard
          case result.integrator:
            of LeapFrog: updates = result.mkLeapFrog
            of Om2MN: 
              stream.add "  lambda = " & $(result.lmbda)
              updates = result.mkOm2MN
            of Om4MN4FP: 
              stream.add "  lambda = " & $(result.lmbda)
              stream.add "  theta = " & $(result.theta)
              stream.add "  rho = " & $(result.rho)
              updates = result.mkOm4MN4FP
            of Om4MN5FV,Om4MN5FP:
              stream.add "  lambda = " & $(result.lmbda)
              stream.add "  vartheta = " & $(result.vartheta)
              stream.add "  theta = " & $(result.theta)
              stream.add "  rho = " & $(result.rho)
              case result.integrator:
                of Om4MN5FV: updates = result.mkOm4MN5FV
                of Om4MN5FP: updates = result.mkOm4MN5FP
                else: discard
            of CustomIntegrator: discard
            of UpdaterIntegrator: discard
          result.updates = newSeq[MolecularDynamicsUpdate]()
          for u in updates: 
            result.updates.add MolecularDynamicsUpdate(update:u[0],scale:u[1])
        of false: qexError "Must specify integrator"
    of HeatbathOverrelax: discard

proc newLatticeSubAction(l: Layout; info: JsonNode; solo: bool): auto = 
  var stream = newMCStream("new subaction")
  stream.add "  nested = " & $(solo)
  result = l.newLatticeSubAction(stream, info, solo, l.SS, l.TT, l.UU, l.WW, l.XX)
  stream.finishStream

proc newLatticeSubAction*[A](self: A; matterInfo: JsonNode; solo: bool): auto =
  let info = self.checkJSON(matterInfo)
  result = self.l[].newLatticeSubAction(info, solo)
  result.setSubActionDiracOperator(self)
  result.setSubActionSmearing(self)

proc newLatticeAction[L:static[int],S,T,U,W,X](
    l: Layout[L];
    info: JsonNode;
    s: typedesc[S];
    t: typedesc[T];
    u: typedesc[U];
    w: typedesc[W];
    x: typedesc[X]  
  ): LatticeAction[L,S,T,U,W,X] =
  var 
    action: ActionType
    stream = newMCStream("new action")
  case info["action-type"].getStr():
    of "gauge": action = PureGauge
    of "matter": action = PureMatter
    of "gauge-matter", "matter-gauge": action = GaugeMatter
    else: qexError "QEX MCMC does not support " & $(action)
  stream.add "  action type = " & info["action-type"].getStr()

  result = LatticeAction[L,S,T,U,W,X](action: action)
  result.subActions = newSeq[LatticeSubAction[L,S,T,U,W,X]]()
  
  new(result.l)
  result.l[] = l

  # Set up action for having matter fields
  case result.action:
    of GaugeMatter,PureMatter:
      var 
        bc = ""
        coeffs = newSeq[float]()
      
      # Matter boundary conditions
      if info["boundary-conditions"].getStr().len > l.nDim:
        qexError "Number of boundary conditions does not match dimension of lattice"
      else: bc = info["boundary-conditions"].getStr()
      stream.add "  boundary conditions = " & info["boundary-conditions"].getStr()

      # Set smearing
      case info["smearing"].getStr():
        of "hyp", "HYP", "nhyp", "nHYP", "hypercubic", "Hypercubic":
          for el in info["smearing-coefficients"].getElems():
            coeffs.add el.getFloat()
          result.smearing = Hypercubic
        of "stout", "Stout": 
          let rho = info["smearing-coefficients"].getElems()
          coeffs.add rho[0].getFloat()
          result.smearing = Stout
          qexError "QEX MCMC does not support stout smearing yet"
        of "none", "None": result.smearing = NoSmearing
        else: qexError "QEX MCMC does not support " & info["smearing"].getStr()
      case result.smearing:
        of Hypercubic, Stout: 
          result.smear = result.l[].newSmearing(
            result.smearing, coeffs, result.l[].SS
          )
          stream.add "  smearing = " & info["smearing"].getStr()
          for idx,el in info["smearing-coefficients"].getElems():
            stream.add "  smearing coefficient " & $(idx) & " = " & $(el.getFloat()) 
        of NoSmearing:
          result.smear = result.l[].newSmearing(result.smearing, @[], result.l[].SS)
      result.smear.bc = bc
    of PureGauge: discard

  stream.finishStream

proc newLatticeAction*(l: Layout; info: JsonNode): auto = 
  result = l.newLatticeAction(info, l.SS, l.TT, l.UU, l.WW, l.XX)

proc initGaugeField(self: var LatticeFieldTheory) =
  var gaugeFieldFound = false

  proc search(self: LatticeFieldTheory; subAction: LatticeSubAction) =
    case subAction.pField.field:
      of GaugeField:
        case gaugeFieldFound:
          of false:
            self.u[] = subAction.pField.u
            gaugeFieldFound = true
          of true: qexError "Found multiple gauge fields in action"
      else: discard
    case subAction.solo:
      of true: discard
      of false:
        for subSubAction in subAction.subActions:
          self.search(subSubAction)

  for action in self.actions:
    for sAction in action.subActions: self.search(sAction)

  if not gaugeFieldFound:
    qexError "No gauge field to do Monte Carlo with"

  for action in self.actions: 
    case action.action:
      of PureGauge: discard
      of GaugeMatter,PureMatter: action.smear.u[] = self.u[]

proc prepHMC(self: LatticeSubAction; pRNG: ParallelRNG) =
  self.running = true
  self.timeStep = 0
  self.spaceStep = 0
  self.iStep = 0
  self.dtau = 1.0 / float(self.steps)
  self.absScale = 0.0

  for step in 0..<self.steps:
    for mdIdx in 0..<self.updates.len:
      case self.updates[mdIdx].update:
        of UpdateT: 
          if step == 0: 
            self.updates[mdIdx].dtau = self.updates[mdIdx].scale*self.dtau
            if mdIdx == 0:
              self.cumSum = self.updates[mdIdx].dtau
              self.absCumSum = abs(self.cumSum)
          self.absScale += abs(self.updates[mdIdx].dtau) 
        of UpdateV: discard
        of UpdateVTV: discard
  for mdIdx in 0..<self.updates.len:
    case self.updates[mdIdx].update:
      of UpdateT: 
        self.updates[mdIdx].absdtau = abs(self.updates[mdIdx].dtau/self.absScale)
      of UpdateV: 
        self.updates[mdIdx].dtau = self.updates[mdIdx].scale*self.dtau
      of UpdateVTV: discard
  
  self.pRNG[] = pRNG

  if not self.solo:
    new(self.p)
    new(self.f)
    new(self.u)
    for saIdx in 0..<self.subActions.len: # Sub-sub-actions are nested
      self.subActions[saIdx].prepHMC(self.pRNG[])

proc setReferences(
    self: LatticeFieldTheory; 
    action: LatticeSubAction; 
    u: auto
  ) =
  var nested = false
  action.p[] = self.p
  action.f[] = self.f
  action.u[] = u
  for sAction in action.subActions:
    if not sAction.solo: 
      self.setReferences(sAction, self.u[])
      nested = true
  action.nested = nested

proc prepHMC(self: var LatticeFieldTheory) =
  var 
    nFields = 0
    nested = false
  for action in self.actions:
    for sAction in action.subActions: 
      sAction.prepHMC(self.pRNG)
      nFields += 1
      if not sAction.solo: 
        self.setReferences(sAction, self.u[])
        nested = true
  self.nested = true
  self.nFields = nFields

proc newLatticeFieldTheory[L:static[int],S,T,U,V,W,X](
    l: Layout[L];
    info: JsonNode;
    stream: var MCStream;
    s: typedesc[S];
    t: typedesc[T];
    u: typedesc[U];
    v: typedesc[V];
    w: typedesc[W];
    x: typedesc[X]
  ): LatticeFieldTheory[L,S,T,U,V,W,X] = 
  var algorithm: MonteCarloAlgorithm

  case info.hasKey("monte-carlo-algorithm"):
    of true: 
      case info["monte-carlo-algorithm"].getStr():
        of "hmc", "HMC", "hybrid-monte-carlo", "hamiltonian-monte-carlo":
          algorithm = HamiltonianMonteCarlo
        of "hb", "HB", "heatbath", "heatbath-overrelax", "heatbath+overrelax":
          algorithm = HeatbathOverrelax
    of false: qexError "Must specify monte carlo algorithm"
  stream.add "  Monte Carlo algorithm = " & info["monte-carlo-algorithm"].getStr()

  result = LatticeFieldTheory[L,S,T,U,V,W,X](algorithm: algorithm)
  result.l = l
  result.actions = newSeq[LatticeAction[L,S,T,U,W,X]]()

  case info.hasKey("start"):
    of true:
      case info["start"].getStr():
        of "cold","unit": result.start = UnitGauge
        of "hot","random": result.start = RandomGauge
        of "read","input","in": result.start = ReadGauge
    of false: result.start = UnitGauge

  case result.algorithm:
    of HamiltonianMonteCarlo:
      var
        sSeed,pSeed: uint64
        sGen,pGen: string

      result.bu = result.l.newGauge()
      result.p = result.l.newGauge()
      result.f = result.l.newGauge()
      new(result.u)

      sSeed = case info.hasKey("serial-random-number-seed")
        of true: (info["serial-random-number-seed"].getInt()).uint64
        of false: 987654321.uint64
      pSeed = case info.hasKey("parallel-random-number-seed")
        of true: (info["parallel-random-number-seed"].getInt()).uint64
        of false: 987654321.uint64

      case info.hasKey("serial-random-number-generator"):
        of true:
          let gen = info["serial-random-number-generator"].getStr()
          case gen:
            of "milc", "MILC", "RngMilc6": sGen = "RngMilc6"
            of "mrg", "MRG", "MRG32k3a": sGen = "MRG32k3a"
            else: qexError gen & " not supported for serial RNG"
        of false: sGen = "RngMilc6"

      case info.hasKey("parallel-random-number-generator"):
        of true:
          let gen = info["parallel-random-number-generator"].getStr()
          case gen:
            of "milc", "MILC", "RngMilc6": pGen = "RngMilc6"
            of "mrg", "MRG", "MRG32k3a": pGen = "MRG32k3a"
            else: qexError gen & " not supported for parallel RNG"
        of false: pGen = "RngMilc6"
      
      result.sRNG = newSerialRNG(sGen, sSeed)
      result.pRNG = result.l.newParallelRNG(pGen, pSeed)

      stream.add "  serial random number generator = " & sGen
      stream.add "  serial random number seed = " & $(sSeed)
      stream.add "  parallel random number generator = " & pGen
      stream.add "  parallel random number seed = " & $(pSeed)
    of HeatbathOverrelax:
      qexError "QEX MCMC does not yet support heatbath + overrelax"

proc newLatticeFieldTheory(info: JsonNode): auto =
  var 
    stream = newMCStream("new lattice field theory", start = true)
    latticeGeometry = newSeq[int]()
    mpiGeometry = newSeq[int]()
    simdGeometry = newSeq[int]()
    (mpiGeomSpecified,simdGeomSpecified) = (false,false)

  # "physical" box geometry
  if not info.hasKey("lattice-geometry"): qexError "Must specify lattice geometry"
  else:
    for idx,el in info["lattice-geometry"].getElems(): latticeGeometry.add el.getInt()

  # rank geometry
  if info.hasKey("mpi-geometry"):
    for idx,el in info["mpi-geometry"].getElems(): mpiGeometry.add el.getInt()
    mpiGeomSpecified = true
  
  #[
  # simd (vector) geometry
  if info.hasKey("inner-geometry"):
    for idx,el in info["inner-geometry"].getElems(): simdGeometry.add el.getInt()
    simdGeomSpecified = true
  ]#

  # set lattice layout & instantiate LatticeFieldTheory object
  let l = case mpiGeomSpecified
    of true: newLayout(latticeGeometry, mpiGeometry)
    of false: newLayout(latticeGeometry)
  result = l.newLatticeFieldTheory(info, stream, l.SS, l.TT, l.UU, l.VV, l.WW, l.XX)
  
  stream.finishStream

template newLatticeFieldTheory*(
    info: JsonNode;
    construction: untyped
  ): auto =
  echo banner
  echo "# ranks: ", nRanks
  threads: echo "# threads: ", numThreads
  block:
    var fieldTheory {.inject.} = newLatticeFieldTheory(info)
    construction
    fieldTheory.initGaugeField
    case fieldTheory.algorithm:
      of HamiltonianMonteCarlo: fieldTheory.prepHMC
      of HeatbathOverrelax: discard
    fieldTheory

if isMainModule:
  qexInit()

  var 
    lat = intSeqParam("lat", @[4, 4, 4, 4])
    lo = lat.newLayout(@[1, 1, 1, 1])
    integrators = @[
      "LF", "2MN", "sexton-weingarten", 
      "4MN4FP", "4MN5FV", "4MN5FP"
    ]
    fieldTypes = @["gauge", "staggered", "wilson"]
    fInfo = %* {"id": "test"}
    saInfo = %* {
      "steps": 10,
      "monte-carlo-algorithm": "hmc"
    }
    actionTypes = @["gauge", "matter", "gauge-matter"]
    aInfo = %* {
      "smearing": "nhyp",
      "smearing-coefficients": @[0.4,0.5,0.5],
      "boundary-conditions": "pppa"
    }
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

  for ft in fieldTypes:
    fInfo["field-type"] = %* ft
    discard lo.newLatticeField(fInfo)

  for integrator in integrators:
    saInfo["integrator"] = %* integrator
    for solo in [true,false]:
      discard lo.newLatticeSubAction(saInfo,solo)

  for at in actionTypes:
    aInfo["action-type"] = %* at
    discard lo.newLatticeAction(aInfo)

  discard newLatticeFieldTheory(ftInfo)

  qexFinalize()