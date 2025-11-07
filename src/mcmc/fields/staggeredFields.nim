import ../mcmcTypes
import ../utilities/rational

import layout
import physics/[qcdTypes, stagSolve]

import strutils
import sequtils
import json

const
  ActionCGTol = 1e-20
  ForceCGTol = 1e-12
  ActionMaxCGIter = 10000
  ForceMaxCGIter = 10000

# Procs for creation of different kinds of staggered matter fields

proc checkJSON(info: JsonNode): JsonNode = 
  result = parseJson("{}")
  for key, keyVal in info: result[key] = keyVal
  if not result.hasKey("cg-tolerance-action"): result["cg-tolerance-action"] = %* ActionCGTol
  if not result.hasKey("cg-maxits-action"): result["cg-maxits-action"] = %* ActionMaxCGIter
  if not result.hasKey("cg-tolerance-force"): result["cg-tolerance-force"] = %* ForceCGTol
  if not result.hasKey("cg-maxits-force"): result["cg-maxits-force"] = %* ForceMaxCGIter
  result["field-type"] = %* "staggered"

proc newStaggeredField(self: var LatticeField; info: JsonNode) = 
  self.stagActionSolverParams = initSolverParams()
  self.stagForceSolverParams = initSolverParams()
  self.stagActionSolverParams.r2req = info["cg-tolerance-action"].getFloat()
  self.stagActionSolverParams.maxits = info["cg-maxits-action"].getInt()
  self.stagForceSolverParams.r2req = info["cg-tolerance-force"].getFloat()
  self.stagForceSolverParams.maxits = info["cg-maxits-force"].getInt()
  self.staggeredMasses = newSeq[float]()

proc newStaggeredFermion*(l: Layout; staggeredInformation: JsonNode): auto = 
  var 
    info = checkJSON(staggeredInformation)
    stream = newMCStream("new staggered fermion")
  if not info.hasKey("mass"): qexError "Must specify fermion mass"

  result = l.newLatticeField(info)
  result.newStaggeredField(info)
  result.staggeredFields = newSeq[l.TT]()

  result.staggeredAction = StaggeredFermion
  result.staggeredFields.add l.ColorVector()
  result.staggeredMasses.add info["mass"].getFloat()

  stream.add "  mass = " & $(info["mass"].getFloat())
  stream.finishStream

proc readJSONArray(input: JsonNode): seq[float] =
  result = newSeq[float]()
  for el in input.getElems(): result.add el.getFloat()

proc newRootedStaggeredFermion*(lo: Layout; staggeredInformation: JsonNode): auto =
  var 
    remezOrder,ratio: string
    info = checkJSON(staggeredInformation)
  let nf = info["nf"].getInt()
  
  if not info.hasKey("mass"): qexError "Must specify fermion mass"
  if not info.hasKey("nf"): qexError "Must specify nf for rooted staggered fermion"
  case $nf:
    of "1","2","3": discard
    else: qexError "Only nf=1,2,3 supported for rooted staggered fermions"
  remezOrder = case info.hasKey("remez-order")
    of true: $info["remez-order"].getInt()
    of false: "15"

  result = lo.newLatticeField(info)
  result.newStaggeredField(info)
  result.staggeredFields = newSeq[lo.TT]()

  result.staggeredAction = RootedStaggeredFermion
  for _ in 0..<parseInt(remezOrder)+1: result.staggeredFields.add lo.ColorVector()
  result.rPhi = lo.ColorVector()
  result.staggeredMasses.add info["mass"].getFloat()

  ratio = case $nf
    of "1": "(1,4)"
    of "2": "(1,2)"
    of "3": "(3,4)"
    else: "(1,4)"
  result.remez = RemezCoefficients(
    nTerms: parseInt(remezOrder),
    f0: rationalCoefficients[ratio][remezOrder]["f0"].getFloat(),
    if0: rationalCoefficients[ratio][remezOrder]["if0"].getFloat(),
    alpha: rationalCoefficients[ratio][remezOrder]["alpha"].readJSONArray(),
    ialpha: rationalCoefficients[ratio][remezOrder]["ialpha"].readJSONArray(),
    beta: rationalCoefficients[ratio][remezOrder]["beta"].readJSONArray(),
    ibeta: rationalCoefficients[ratio][remezOrder]["ibeta"].readJSONArray()
  )

proc newStaggeredBoson*(l: Layout; staggeredInformation: JsonNode): auto = 
  var 
    info = checkJSON(staggeredInformation)
    stream = newMCStream("new staggered boson")
  if not info.hasKey("mass"): qexError "Must specify boson mass"

  result = l.newLatticeField(info)
  result.newStaggeredField(info)
  result.staggeredFields = newSeq[l.TT]()

  result.staggeredAction = StaggeredBoson
  result.staggeredFields.add l.ColorVector()
  result.staggeredMasses.add info["mass"].getFloat()

  stream.add "  mass = " & $(info["mass"].getFloat())
  stream.finishStream

proc newStaggeredHasenbuschFermion*(
    l: Layout; 
    staggeredInformation: JsonNode
  ): auto = 
  var 
    info = checkJSON(staggeredInformation)
    stream = newMCStream("new staggered Hasenbusch fermion")
  if not info.hasKey("mass1"): qexError "Must specify first Hasenbusch mass"
  if not info.hasKey("mass2"): qexError "Must specify second Hasenbusch mass"
  let masses = @[info["mass1"].getFloat(), info["mass2"].getFloat()]

  result = l.newLatticeField(info)
  result.newStaggeredField(info)
  result.staggeredFields = newSeq[l.TT]()

  result.staggeredAction = StaggeredHasenbuschFermion
  for mass in masses:
    result.staggeredFields.add l.ColorVector()
    result.staggeredMasses.add mass

  stream.add "  mass1 = " & $(info["mass1"].getFloat())
  stream.add "  mass2 = " & $(info["mass2"].getFloat())
  stream.finishStream

# Generic procs for operations w/ staggered matter fields

proc phi(self: LatticeField): auto = 
  result = case self.staggeredAction:
    of RootedStaggeredFermion: self.rPhi
    else: self.staggeredFields[0]

proc phi1(self: LatticeField): auto = self.staggeredFields[0]

proc phi2(self: LatticeField): auto = self.staggeredFields[1]

proc mass(self: LatticeField): float = self.staggeredMasses[0]

proc mass1(self: LatticeField): float = self.staggeredMasses[0]

proc mass2(self: LatticeField): float = self.staggeredMasses[1]

proc zero(phi: auto) =
  threads: phi := 0

proc zero[T](phis: seq[T]) =
  threads: 
    for phi in phis: phi := 0

proc zeroOdd(phi: auto) =
  threads: phi.odd := 0

proc sq(x: float): float = x*x

proc normSquared(psi: auto): float =
  var nrm2 = 0.0
  threads:
    let psi2 = psi.norm2
    threadBarrier()
    threadMaster: nrm2 = psi2
  result = nrm2

proc applyD(
    D: DiracOperator;
    psi: auto;
    phi: auto;
    mass: float;
  ) = 
  threads: D(D.stag, psi, phi, mass)

proc applyDdag(
    D: DiracOperator;
    psi: auto;
    phi: auto;
    mass: float;
  ) =
  threads: D(D.stag, psi, phi, -mass)

proc rationalApprox[T](
    phi: auto; 
    psi: T; 
    phis: seq[T]; 
    alpha0: float; 
    alpha: seq[float]
  ) = 
  threads:
    phi := alpha0*psi
    threadBarrier()
    for idx in 0..<alpha.len: 
      phi += alpha[idx]*phis[idx+1]

proc applyD[T](
    D: DiracOperator;
    phi: T;
    psi: T;
    phis: seq[T];
    mass,f: float;
    alpha,beta: seq[float];
    sp0: var SolverParams
  ) =
  var shifts = newSeq[float](beta.len+1)
  for idx in 0..<shifts.len: 
    shifts[idx] = case idx == 0 
      of true: mass
      of false: mass + beta[idx-1]
  solve(D.stag, phis, psi, shifts, sp0)
  phi.rationalApprox(psi,phis,f,alpha)

proc applyDdag[T](
    D: DiracOperator;
    phi: T;
    psi: T;
    phis: seq[T];
    mass,f: float;
    alpha,beta: seq[float],
    sp0: var SolverParams;
    rescale: bool = false
  ) = D.applyD(phi,psi,phis,-mass,f,alpha,beta,sp0)

proc applyNegDdagOdd(
    D: DiracOperator;
    psi: auto;
    phi: auto;
  ) =
  threads:
    stagD2(D.stag.so, psi, D.stag.g, phi, 0, 0)
    threadBarrier()
    psi.odd := -0.5*psi
    psi.even := 0

proc applyDdag2OddAndReplaceEven(
    D: DiracOperator;
    psi: auto;
    phi: auto;
  ) =
  threads:
    stagD2(D.stag.so, psi, D.stag.g, phi, 0, 0)
    threadBarrier()
    psi.even := phi

proc solveD*(
    D: DiracOperator;
    psi: auto;
    phi: auto;
    mass: float;
    sp0: var SolverParams
  ) =
  psi.zero
  if mass > Small64: solve(D.stag, psi, phi, mass, sp0)
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

proc solveDdag*(
    D: DiracOperator;
    psi: auto;
    phi: auto;
    mass: float;
    sp0: var SolverParams
  ) =
  psi.zero
  if mass > epsilon(float): solve(D.stag, psi, phi, -mass, sp0)
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

proc solveD[T](
    D: DiracOperator;
    phis: seq[T];
    psi: T;
    mass: float;
    beta: seq[float];
    sp0: var SolverParams
  ) =
  var shifts = newSeq[float](beta.len+1)
  for idx in 0..<shifts.len: 
    shifts[idx] = case idx == 0 
      of true: mass
      of false: mass + beta[idx-1]
  solve(D.stag, phis, psi, shifts, sp0)

proc solveD[T](
    D: DiracOperator;
    phis: seq[T];
    psi: T;
    masses: seq[float];
    sp0: var SolverParams
  ) = solve(D.stag, phis, psi, masses, sp0)

proc solveDdag[T](
    D: DiracOperator;
    phis: seq[T];
    psi: T;
    mass: float;
    beta: seq[float];
    sp0: var SolverParams
  ) = D.solveD(phis,psi,-mass,beta,sp0)

proc solveDdag[T](
    D: DiracOperator;
    phis: seq[T];
    psi: T;
    masses: seq[float];
    sp0: var SolverParams
  ) = 
  var nmasses = newSeq[float](masses.len)
  for idx in 0..<nmasses.len: nmasses[idx] = -masses[idx]
  solve(D.stag, phis, psi, nmasses, sp0)

proc outer(f: auto; psi: auto; shifter: auto; dtau: float) =
    let n = psi[0].len
    threads:
      for mu in 0..<f.len:
        for s in f[mu]:
          forO a, 0, n-1:
            forO b, 0, n-1:
              f[mu][s][a,b]+=dtau*psi[s][a]*shifter[mu].field[s][b].adj

#[ 
Methods for 
1.) getting phi, 
2.) calculating action, 
3.) calculating (partial) force
]#

proc getStaggeredField*(self: var LatticeField; D: DiracOperator) =
  case self.staggeredAction:
    of StaggeredFermion,StaggeredBoson: zero(self.phi)
    of RootedStaggeredFermion:
      zero(self.phi)
      for f in 0..<self.remez.nTerms: zero(self.staggeredFields[f])
    of StaggeredHasenbuschFermion:
      zero(self.phi1)
      zero(self.phi2)
  case self.staggeredAction:
    of StaggeredFermion: D.applyDdag(self.phi, D.stagPsi, self.mass)
    of RootedStaggeredFermion:
      D.applyDdag(
        self.phi, 
        D.stagPsi, 
        self.staggeredFields,
        self.mass,
        self.remez.f0,
        self.remez.alpha,
        self.remez.beta,
        self.stagActionSolverParams
      ) # Applies D^{+}; still a solve when rooted
    of StaggeredHasenbuschFermion: 
      D.applyDdag(self.phi2, D.stagPsi, self.mass1)
      D.solveDdag(self.phi1, self.phi2, self.mass2, self.stagActionSolverParams)
    of StaggeredBoson: 
      D.solveDdag(self.phi, D.stagPsi, self.mass, self.stagActionSolverParams)
  zeroOdd(self.phi)

proc staggeredAction*(self: var LatticeField; D: var DiracOperator): float =
  case self.staggeredAction:
    of StaggeredFermion,StaggeredHasenbuschFermion,StaggeredBoson: zero(D.stagPsi)
    of RootedStaggeredFermion: zero(self.staggeredFields)
  case self.staggeredAction:
    of StaggeredFermion: 
      D.solveDdag(D.stagPsi, self.phi, self.mass, self.stagActionSolverParams)
      if self.mass <= Small64: D.applyNegDdagOdd(D.stagPsi, D.stagPsi)
    of RootedStaggeredFermion:
      D.applyDdag( 
        D.stagPsi,
        self.phi, 
        self.staggeredFields,
        self.mass,
        self.remez.if0,
        self.remez.ialpha,
        self.remez.ibeta,
        self.stagActionSolverParams
      ) # Applies (rational) inverse of D^{+}
    of StaggeredHasenbuschFermion:
      zero(D.stagPsi)
      D.applyD(self.phi2, self.phi1, self.mass2)
      D.solveD(D.stagPsi, self.phi2, self.mass1, self.stagActionSolverParams)
    of StaggeredBoson: D.applyDdag(D.stagPsi, self.phi, self.mass)
  case self.staggeredAction:
    of StaggeredHasenbuschFermion,StaggeredBoson: 
      result = 0.5*D.stagPsi.normSquared
    of StaggeredFermion, RootedStaggeredFermion:
      result = 0.5*D.stagPsi.normSquared

proc stagForce*(
    self: var LatticeField;
    D: var DiracOperator;
    fdtau: float;
    f: auto
  ) = 
  var dtau = case self.staggeredAction
    of StaggeredFermion,StaggeredHasenbuschFermion: -0.5*fdtau
    of StaggeredBoson: -0.25*fdtau
    of RootedStaggeredFermion: 0.0
  case self.staggeredAction:
    of StaggeredFermion,StaggeredHasenbuschFermion: D.stagPsi.zero
    of StaggeredBoson: D.stagPsi.zero
    of RootedStaggeredFermion: discard
  case self.staggeredAction:
    of StaggeredFermion:
      D.solveD(D.stagPsi, self.phi, self.mass, self.stagForceSolverParams)
      if self.mass < Small64: D.applyDdag2OddAndReplaceEven(D.stagPsi, D.stagPsi)
    of RootedStaggeredFermion:
      var effMasses = newSeq[float](self.staggeredFields.len)
      for idx in 0..<self.staggeredFields.len:
        effMasses[idx] = case idx == 0
          of true: self.mass
          of false: sqrt(self.mass*self.mass+self.remez.ibeta[idx-1])
      D.solveD(
        self.staggeredFields,
        self.phi,
        effMasses,
        self.stagForceSolverParams
      ) # Multimass solve
    of StaggeredHasenbuschFermion: 
      D.solveD(D.stagPsi, self.phi1, self.mass1, self.stagForceSolverParams)
    of StaggeredBoson: D.applyDdag2OddAndReplaceEven(D.stagPsi, self.phi)
  case self.staggeredAction:
    of StaggeredFermion, StaggeredHasenbuschFermion, StaggeredBoson:
      for mu in 0..<f.len: discard D.stagShifter[mu] ^* D.stagPsi
    of RootedStaggeredFermion: discard
  case self.staggeredAction:
    of StaggeredFermion:
      if self.mass > Small64: dtau = dtau/self.mass
      else: dtau = -0.5*dtau
    of StaggeredHasenbuschFermion: 
      dtau = dtau*(self.mass2.sq-self.mass1.sq)/self.mass1
    of StaggeredBoson,RootedStaggeredFermion: discard
  case self.staggeredAction:
    of StaggeredFermion,StaggeredBoson: f.outer(D.stagPsi, D.stagShifter, dtau)
    of StaggeredHasenbuschFermion: f.outer(D.stagPsi, D.stagShifter, dtau)
    of RootedStaggeredFermion:
      for idx in 0..<self.remez.ialpha.len:
        dtau = -0.5*fdtau*self.remez.ialpha[idx]
        dtau = dtau/sqrt(self.mass*self.mass+self.remez.ibeta[idx])
        for mu in 0..<f.len: discard D.stagShifter[mu] ^* self.staggeredFields[idx+1]
        f.outer(self.staggeredFields[idx+1], D.stagShifter, dtau)