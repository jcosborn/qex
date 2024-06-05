import qex
import base
import layout
import field
import gauge
import physics/[qcdTypes, stagD, stagSolve, wilsonD, wilsonSolve]

import tables
import math
import typetraits
import strformat, strutils
import json

import monteCarlo
import fieldUtilities

export monteCarlo
export fieldUtilities

#[ Definitions of abstract data types ]#

type
  GaugeActionType = enum Wilson, Adjoint, Rectangle, Symanzik, Iwasaki, DBW2
  FieldType* = enum Gauge, StaggeredFermion, WilsonFermion, StaggeredBoson, WilsonBoson

type
  DiracOperator*[S,U,V,W] = object
    case discretization: FieldType
      of StaggeredFermion:
        stag*: Staggered[S,V]
        stagShifter*: seq[Shifter[U,V]]
      of WilsonFermion:
        wils*: Wilson[S,W]
      else: discard

  AbstractField*[S,T,U] = object of IntegratorAtom
    
    actionSolverParameters*: SolverParams
    forceSolverParameters*: SolverParams

    case field*: FieldType
      of Gauge:
        u*: seq[S]
        action*: GaugeActionType
        actionCoefficients*: GaugeActionCoeffs
        start*: string
      of StaggeredFermion, StaggeredBoson:
        sPhi*: T # Boson field
        sMass*: float # Mass of boson field
      of WilsonFermion, WilsonBoson:
        wPhi*: U # Boson field
        wMass*: float # Mass of boson field

const
  # Default value for beta_A/beta_F; see gauge/gaugeAction.nim
  BetaAOverBetaF = -0.25

  # Default "c1" rectangle coefficients; see gauge/gaugeAction.nim
  C1Symanzik = -1.0/12.0
  C1Iwasaki = -0.331
  C1DBW2 = -1.4088

  # Default settings for conjugate gradient
  ActionCGTol = 1e-20
  ForceCGTol = 1e-12
  ActionMaxCGIter = 10000
  ForceMaxCGIter = 10000

converter toGaugeActionType(s: string):
  GaugeActionType = parseEnum[GaugeActionType](s)

proc mass*(self: AbstractField): float =
  result = case self.field:
    of StaggeredFermion, StaggeredBoson: self.sMass
    of WilsonFermion, WilsonBoson: self.wMass
    else: 0.0

#[ AbstractField constructors ]#

proc new(
    self: var AbstractField;
    l: Layout;
    action: string;
    bareGaugeCoupling, adjointRatio, rectangleCoefficient: float;
    start: string;
  ) =  
  self.u = l.newGauge()
  self.action = toGaugeActionType(action)
  self.actionCoefficients = case self.action
    of Wilson: GaugeActionCoeffs(plaq: bareGaugeCoupling)
    of Adjoint:
      GaugeActionCoeffs(
        plaq: bareGaugeCoupling,
        adjplaq: bareGaugeCoupling * adjointRatio
      )
    of Rectangle: gaugeActRect(bareGaugeCoupling, rectangleCoefficient)
    of Symanzik: gaugeActRect(bareGaugeCoupling, C1Symanzik)
    of Iwasaki: gaugeActRect(bareGaugeCoupling, C1Iwasaki)
    of DBW2: gaugeActRect(bareGaugeCoupling, C1DBW2)

  self.start = start

proc new(
    self: var AbstractField;
    l: Layout;
    mass: float;
    aTol, fTol: float;
    aMaxCG, fMaxCG: int;
  ) =
  case self.field:
    of StaggeredFermion, StaggeredBoson:
      self.sPhi = l.ColorVector()
      self.sMass = mass
    of WilsonFermion, WilsonBoson:
      self.wPhi = l.DiracFermion()
      self.wMass = mass
    else: discard

  self.actionSolverParameters = initSolverParams()
  self.actionSolverParameters.r2req = aTol
  self.actionSolverParameters.maxits = aMaxCG

  self.forceSolverParameters = initSolverParams()
  self.forceSolverParameters.r2req = fTol
  self.forceSolverParameters.maxits = fMaxCG

proc newDiracOperator[S,U,V,W](
    g: auto;
    discretization: FieldType;
    s: typedesc[S];
    u: typedesc[U];
    v: typedesc[V];
    w: typedesc[W]
  ): DiracOperator[S,U,V,W] =
  result = DiracOperator[S,U,V,W](discretization: discretization)
  case result.discretization:
    of Gauge, StaggeredBoson, WilsonBoson: discard
    of StaggeredFermion:
      result.stag = newStag(g)
      result.stagShifter = newSeq[Shifter[U,V]](g[0].l.nDim)
    of WilsonFermion: result.wils = newWilson(g)

proc newDiracOperator*(g: auto; discretization: FieldType): auto =
  let l = g[0].l

  template s: untyped =
    type(l.ColorMatrix())

  template u: untyped = 
    type(l.ColorVector())

  template v: untyped =
    type(l.ColorVector()[0])

  template w: untyped =
    type(spproj1p(l.DiracFermion()[0]))

  result = g.newDiracOperator(discretization, s, u, v, w)

proc newGaugeField[S,T,U](
    l: Layout;
    action: string;
    bareGaugeCoupling, adjointRatio, rectangleCoefficient: float;
    integrator: IntegratorAtomType;
    steps: int;
    start: string;
    s: typedesc[S];
    t: typedesc[T];
    u: typedesc[U]
  ): AbstractField[S,T,U] =
  result = AbstractField[S,T,U](field: Gauge, integrator: integrator)
  result.constructIntegratorAtom(steps)
  result.new(l, action, bareGaugeCoupling, adjointRatio, rectangleCoefficient, start)

proc newGaugeField*(
    l: Layout;
    action: string;
    bareGaugeCoupling: float;
    steps: int;
    start: string;
    adjointRatio: float = BetaAOverBetaF,
    rectangleCoefficient: float = C1Symanzik;
    integrator: string = "MN2"
  ): auto = l.newGaugeField(
    action, bareGaugeCoupling, adjointRatio, rectangleCoefficient,
    toIntegratorAtomType(integrator), steps, start,
    type(l.ColorMatrix()),
    type(l.ColorVector()),
    type(l.DiracFermion())
  )

proc newGaugeField*(l: Layout; info: JsonNode): auto =
  if not info.hasKey("adjointRatio"):
    info["adjointRatio"] = %* BetaAOverBetaF
  if not info.hasKey("rectangleCoefficient"):
    info["rectangleCoefficient"] = %* C1Symanzik
  if not info.hasKey("integrator"): info["integrator"] = %* "MN2"
  if not info.hasKey("start"): info["start"] = %* "cold"

  result = l.newGaugeField(
    info["action"].getStr(),
    info["beta"].getFloat(),
    info["adjointRatio"].getFloat(),
    info["rectangleCoefficient"].getFloat(),
    toIntegratorAtomType(info["integrator"].getStr()),
    info["steps"].getInt(),
    info["start"].getStr(),
    type(l.ColorMatrix()),
    type(l.ColorVector()),
    type(l.DiracFermion())
  )

proc newMatterField[S,T,U](
    l: Layout;
    field: FieldType;
    mass: float;
    integrator: IntegratorAtomType;
    steps: int;
    aTol, fTol: float;
    aMaxCG, fMaxCG: int;
    s: typedesc[S];
    t: typedesc[T];
    u: typedesc[U]
  ): AbstractField[S,T,U] =
  result = AbstractField[S,T,U](field: field, integrator: integrator)
  result.constructIntegratorAtom(steps)
  result.new(l, mass, aTol, fTol, aMaxCG, fMaxCG)

proc newMatterField*(
    l: Layout;
    field: FieldType;
    mass: float;
    steps: int;
    integrator: string = "MN2";
    actionTolCG: float = ActionCGTol;
    forceTolCG: float = ForceCGTol;
    actionMaxItnCG: int = ActionMaxCGIter;
    forceMaxItnCG: int = ForceMaxCGIter
  ): auto = l.newMatterField(
    field, mass,
    toIntegratorAtomType(integrator), steps,
    actionTolCG, forceTolCG,
    actionMaxItnCG, forceMaxItnCG,
    type(l.ColorMatrix()),
    type(l.ColorVector()),
    type(l.DiracFermion())
  )

proc newMatterField*(
    l: Layout;
    field: FieldType;
    info: JsonNode
  ): auto =

  if not info.hasKey("mass"): qexError "must specify fermion mass"
  if not info.hasKey("steps"): qexError "must specify number of fermion integration steps"
  
  if not info.hasKey("integrator"): info["integrator"] = %* "MN2"
  if not info.hasKey("cg-tolerance-action"): info["cg-tolerance-action"] = %* ActionCGTol
  if not info.hasKey("cg-maxits-action"): info["cg-maxits-action"] = %* ActionMaxCGIter
  if not info.hasKey("cg-tolerance-force"): info["cg-tolerance-force"] = %* ForceCGTol
  if not info.hasKey("cg-maxits-force"): info["cg-maxits-force"] = %* ForceMaxCGIter

  result = l.newMatterField(
    field, 
    info["mass"].getFloat(),
    toIntegratorAtomType(info["integrator"].getStr()),
    info["steps"].getInt(),
    info["cg-tolerance-action"].getFloat(),
    info["cg-tolerance-force"].getFloat(),
    info["cg-maxits-action"].getInt(),
    info["cg-maxits-force"].getInt(),
    type(l.ColorMatrix()),
    type(l.ColorVector()),
    type(l.DiracFermion())
  )

  for key in info.keys():
    case key:
      of "lambda": result.lambda = info[key].getFloat()
      of "rho": result.rho = info[key].getFloat()
      of "theta": result.theta = info[key].getFloat()
      of "vartheta": result.vartheta = info[key].getFloat()

#[ Gauge AbstractField methods ]#

proc force*(self: AbstractField; f: auto) =
  if self.b != 0:
    case self.action:
      of Adjoint: self.actionCoefficients.forceA(self.u, f)
      else: self.actionCoefficients.gaugeForce(self.u, f)
    threads:
      for mu in 0..<f.len:
        for s in f[mu]: f[mu][s] := self.b*f[mu][s]

proc action*(self: AbstractField): float =
  result = case self.action
    of Adjoint: self.actionCoefficients.actionA(self.u)
    else: self.actionCoefficients.gaugeAction1(self.u)

#[ Fermion AbstractField methods ]#

proc applyD(
    self: AbstractField;
    D: DiracOperator;
    psi: auto;
    phi: auto
  ) =
  let mass = self.mass
  case self.field:
    of StaggeredFermion, StaggeredBoson:
      threads: D(D.stag, psi, phi, mass)
    of WilsonFermion, WilsonBoson: discard
    else: discard

proc applyMasslesDDag1(
    self: AbstractField;
    D: DiracOperator;
    x: auto;
    b: auto
  ) =
  case self.field:
    of StaggeredFermion, StaggeredBoson:
      threads:
        stagD2(D.stag.so, x, D.stag.g, b, 0, 0)
        threadBarrier()
        x.odd := -0.5*x
        x.even := 0
    of WilsonFermion, WilsonBoson: discard
    else: discard

proc applyMasslesDDag2(
    self: AbstractField;
    D: DiracOperator;
    x: auto;
    b: auto
  ) =
  case self.field:
    of StaggeredFermion, StaggeredBoson:
      threads:
        stagD2(D.stag.so, x, D.stag.g, b, 0, 0)
        threadBarrier()
        x.even := b
    of WilsonFermion, WilsonBoson: discard
    else: discard

proc solve*(
    self: AbstractField;
    D: DiracOperator;
    psi: auto;
    phi: auto;
    sp0: var SolverParams;
    mFac: float = 1;
  ) = 
  threads: psi := 0
  case self.field:
    of Gauge: discard
    of StaggeredFermion, StaggeredBoson:
      if self.mass != 0: solve(D.stag, psi, phi, mFac*self.mass, sp0)
      else:
        var sp = sp0
        sp.resetStats()
        sp.verbosity = sp0.verbosity
        sp.usePrevSoln = false

        threads: psi := 0
        solveEE(D.stag, psi, phi, 0, sp)
        threads: psi.even := 4*psi

        sp.calls = 1
        sp0.addStats(sp)
    of WilsonFermion, WilsonBoson: discard # Not yet implemented

proc norm(psi: auto): auto =
  var term = 0.0
  threads:
    let psi2 = psi.norm2()
    threadBarrier()
    threadMaster: term = psi2
  result = term

proc action*(
  self: var AbstractField;
  D: DiracOperator;
  psi: auto
  ): float =
  var action = newSeq[float]()
  case self.field:
    of StaggeredFermion:
      threads: psi := 0
      self.solve(D, psi, self.sPhi, self.actionSolverParameters, mFac = -1)
      if self.mass == 0: self.applyMasslesDDag1(D, psi, psi)
      action.add 0.5*psi.norm
    of StaggeredBoson:
      threads: psi := 0
      self.applyD(D, psi, self.sPhi)
      action.add 0.5*psi.norm
    of WilsonFermion: discard # Not implemented
    of WilsonBoson: discard # Not implemented
    else: discard

  # Result action as sum over contributions
  result = sum(action)

proc stagForce(
    self: AbstractField;
    psi: auto;
    shifter: auto;
    f: auto
  ) =
  var sc = -0.5 * self.b
  let n = psi[0].len
  
  case self.field:
    of StaggeredFermion:
      if self.mass != 0: sc = sc / self.mass
      else: sc = -0.5 * sc
    of StaggeredBoson: sc = 0.5 * sc
    else: discard

  threads:
    for mu in 0..<f.len:
      for s in f[mu]:
        forO a, 0, n-1:
          forO b, 0, n-1:
            f[mu][s][a,b] += sc*psi[s][a]*shifter[mu].field[s][b].adj

proc force*(
    self: var AbstractField;
    D: var DiracOperator;
    psi: auto;
    u: auto;
    f: auto
  ) =
  if self.b != 0:
    case self.field:
      of StaggeredFermion, StaggeredBoson:
        threads: psi := 0
      of WilsonFermion, WilsonBoson: discard # Not implemented
      else: discard

    case self.field:
      of StaggeredFermion:
        self.solve(D, psi, self.sPhi, self.forceSolverParameters)
        if self.mass == 0: self.applyMasslesDDag2(D, psi, psi)
      of StaggeredBoson: self.applyMasslesDDag2(D, psi, self.sPhi)
      of WilsonFermion: discard # Not implemented
      of WilsonBoson: discard # Not implemented
      else: discard

    case self.field:
      of StaggeredFermion, StaggeredBoson:
        for mu in 0..<f.len: discard D.stagShifter[mu] ^* psi
        self.stagForce(psi, D.stagShifter, f)
      of WilsonFermion, WilsonBoson: discard # Not implemented
      else: discard

if isMainModule:
  qexinit()

  let
    lat = intSeqParam("lat", @[4, 4, 4, 4])
    lo = lat.newLayout(@[1, 1, 1, 1])
    u = lo.newGauge()
    f = lo.newGauge()
    sPsi = lo.ColorVector()
    hPsi = lo.DiracFermion()

    beta = 6.0

    massF1 = 0.0
    massF2 = 0.5

    integrator = "MN2"
    steps = 10

    numberTestPV = 4
    numberTestH = 4

    massPV = 0.75
    massH = 0.5

    gaugeFieldParams = %* {
      "action": "Wilson",
      "beta": beta,
      "steps": steps,
      "integrator": integrator,
      "start": "cold"
    }

  # Test gauge
  var gauge = lo.newGaugeField("Wilson", beta, steps, "cold", integrator = integrator)
  gauge.force(f)
  discard gauge.action()

  # Equivalent way to create gauge field
  gauge = lo.newGaugeField(gaugeFieldParams)

  # Test creation of Staggered fermion fields
  var
    staggered1 = lo.newMatterField(
      StaggeredFermion, massF1, steps, integrator = integrator
    )
    staggered2 = lo.newMatterField(
      StaggeredFermion, massF2, steps,
      integrator = integrator
    )

  # Test adding staggered Pauli-Villars & Hasenbusch
  var sPauliVillars = lo.newMatterField(StaggeredBoson, massPV, steps)

  # Test creation of Wilson fermions
  var
    wilson1 = lo.newMatterField(WilsonFermion, massF1, steps)
    wilson2 = lo.newMatterField(
      WilsonFermion, massF2, steps,
      integrator = integrator
    )

  # Test adding Wilson Pauli-Villars/Hasenbusch
  var wPauliVillars = lo.newMatterField(WilsonBoson, massPV, steps)

  # Test out Staggered Dirac Operator
  var sD = u.newDiracOperator(StaggeredFermion)
  echo "staggered action fermiom: ", staggered1.action(sD, sPsi)
  echo "staggered Pauli-Villars boson: ", sPauliVillars.action(sD, sPsi)

  qexfinalize()