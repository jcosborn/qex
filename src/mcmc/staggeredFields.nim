import qex
import layout
import physics/[qcdTypes, stagSolve]

import strutils
import json

import abstractFields
import diracOperator
import typeUtilities

const
  ActionCGTol = 1e-20
  ForceCGTol = 1e-12
  ActionMaxCGIter = 10000
  ForceMaxCGIter = 10000

# Procs for creation of different kinds of staggered matter fields

proc checkJSON(info: JsonNode): JsonNode = 
  result = checkMonteCarloAlgorithm(info)
  if not result.hasKey("cg-tolerance-action"): 
    result["cg-tolerance-action"] = %* ActionCGTol
  if not result.hasKey("cg-maxits-action"): 
    result["cg-maxits-action"] = %* ActionMaxCGIter
  if not result.hasKey("cg-tolerance-force"): 
    result["cg-tolerance-force"] = %* ForceCGTol
  if not result.hasKey("cg-maxits-force"): 
    result["cg-maxits-force"] = %* ForceMaxCGIter
  if not result.hasKey("steps"): qexError "Missing # integrator steps for fermion"

proc newStaggeredField*[S,T,U](
    staggeredInformation: JsonNode;
    s: typedesc[S];
    t: typedesc[T];
    u: typedesc[U];
  ): AbstractField[S,T,U] = 
  
  # Check standard JSON info
  var info = checkJSON(staggeredInformation)
  let algorithm = toMonteCarloType(info["monte-carlo-algorithm"].getStr())

  # Standard stuff that is necessary for all abstract fields
  result = AbstractField[S,T,U](field: StaggeredMatterField, algorithm: algorithm)
  result.newAbstractField(info)

  # Initialize staggered field parameters
  result.staggeredActionSolverParameters = initSolverParams()
  result.staggeredForceSolverParameters = initSolverParams()
  result.staggeredActionSolverParameters.r2req = info["cg-tolerance-action"].getFloat()
  result.staggeredActionSolverParameters.maxits = info["cg-maxits-action"].getInt()
  result.staggeredForceSolverParameters.r2req = info["cg-tolerance-force"].getFloat()
  result.staggeredForceSolverParameters.maxits = info["cg-maxits-force"].getInt()

  # Initialize array of fields & masses
  result.staggeredFields = newSeq[t]()
  result.staggeredMasses = newSeq[float]()

proc newStaggeredFermion*(l: Layout; info: JsonNode): auto = 

  # Check Json information
  if not info.hasKey("mass"): qexError "Must specify fermion mass"

  # Create abstract field object
  result = newStaggeredField(info, l.typeS, l.typeT, l.typeU)

  # Create data for staggered fermion
  result.staggeredAction = StaggeredFermion
  result.staggeredFields.add l.ColorVector()
  result.staggeredMasses.add info["mass"].getFloat()

proc newStaggeredBoson*(l: Layout; info: JsonNode): auto = 

  # Check Json information
  if not info.hasKey("mass"): qexError "Must specify boson mass"

  # Create abstract field object
  result = newStaggeredField(info, l.typeS, l.typeT, l.typeU)

  # Create data for staggered boson
  result.staggeredAction = StaggeredBoson
  result.staggeredFields.add l.ColorVector()
  result.staggeredMasses.add info["mass"].getFloat()

proc newStaggeredHasenbuschFermion*(l: Layout; info: JsonNode): auto = 

  # Check Json information
  if not info.hasKey("mass1"): qexError "Must specify first Hasenbusch mass"
  if not info.hasKey("mass2"): qexError "Must specify second Hasenbusch mass"
  let masses = @[info["mass1"].getFloat(), info["mass2"].getFloat()]

  # Create abstract field object
  result = newStaggeredField(info, l.typeS, l.typeT, l.typeU)

  # Create data for staggered Hasenbusch fermion
  result.staggeredAction = StaggeredFermion
  for mass in masses:
    result.staggeredFields.add l.ColorVector()
    result.staggeredMasses.add mass

# Generic procs for operations w/ staggered matter fields

proc phi(self: AbstractField): auto = self.staggeredFields[0]

proc phi1(self: AbstractField): auto = self.staggeredFields[0]

proc phi2(self: AbstractField): auto = self.staggeredFields[1]

proc mass(self: AbstractField): float = self.staggeredMasses[0]

proc mass1(self: AbstractField): float = self.staggeredMasses[0]

proc mass2(self: AbstractField): float = self.staggeredMasses[1]

proc zero(phi: auto) =
  threads: phi := 0

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
  threads: Ddag(D.stag, psi, phi, mass)

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
  if mass != 0: solve(D.stag, psi, phi, mass, sp0)
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
  if mass != 0: solve(D.stag, psi, phi, -mass, sp0)
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

proc outer(f: auto; psi: auto; shifter: auto; dtau: float) =
    let n = psi[0].len
    threads:
      for mu in 0..<f.len:
        for s in f[mu]:
          forO a, 0, n-1:
            forO b, 0, n-1:
              f[mu][s][a,b] += dtau * psi[s][a] * shifter[mu].field[s][b].adj

#[ 
Methods for 
1.) getting phi, 
2.) calculating action, 
3.) calculating (partial) force
]#

proc getStaggeredField*(
    self: var AbstractField;
    D: DiracOperator;
    psi: auto
  ) =
  case self.staggeredAction:
    of StaggeredFermion: D.applyDdag(self.phi, psi, self.mass)
    of StaggeredHasenbuschFermion: 
      D.solveDdag(self.phi2, psi, self.mass2, self.staggeredActionSolverParameters)
      D.applyDdag(self.phi1, self.phi2, self.mass1)
    of StaggeredBoson: 
      D.solveDdag(self.phi, psi, self.mass, self.staggeredActionSolverParameters)
  zeroOdd(self.phi)

proc staggeredAction*(
    self: var AbstractField;
    D: DiracOperator;
    psi: auto
  ): float =
  psi.zero
  case self.staggeredAction:
    of StaggeredFermion: 
      D.solveDdag(psi, self.phi, self.mass, self.staggeredActionSolverParameters)
      if self.mass == 0.0: D.applyNegDdagOdd(psi, psi)
    of StaggeredHasenbuschFermion:
      zero(self.phi2)
      D.applyDdag(self.phi2, self.phi1, self.mass2)
      D.solveDdag(psi, self.phi2, self.mass1, self.staggeredActionSolverParameters)
    of StaggeredBoson: D.applyDdag(psi, self.phi, self.mass)
  result = 0.5 * psi.normSquared

proc staggeredPartialForce*(
    self: var AbstractField;
    D: var DiracOperator;
    psi: auto;
    f: auto
  ) = 
  var dtau = -0.5 * self.dtau[0]
  psi.zero
  case self.staggeredAction:
    of StaggeredFermion:
      D.solveD(psi, self.phi, self.mass, self.staggeredForceSolverParameters)
      if self.mass == 0.0: D.applyDdag2OddAndReplaceEven(psi, psi)
    of StaggeredHasenbuschFermion: 
      D.solveD(psi, self.phi1, self.mass1, self.staggeredForceSolverParameters)
    of StaggeredBoson: D.applyDdag2OddAndReplaceEven(psi, self.phi)
  case self.staggeredAction:
    of StaggeredFermion, StaggeredHasenbuschFermion, StaggeredBoson:
      for mu in 0..<f.len: discard D.stagShifter[mu] ^* psi
  case self.staggeredAction:
    of StaggeredFermion:
      if self.mass != 0.0: dtau = dtau / self.mass
      else: dtau = -0.5 * dtau
    of StaggeredHasenbuschFermion: 
      dtau = dtau * (self.mass2.sq - self.mass1.sq) / self.mass1
    of StaggeredBoson: dtau = 0.5 * dtau
  case self.staggeredAction:
    of StaggeredFermion, StaggeredHasenbuschFermion, StaggeredBoson:
      f.outer(psi, D.stagShifter, dtau)