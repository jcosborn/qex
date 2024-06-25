import qex
import base
import layout
import field
import gauge
import physics/[qcdTypes, stagSolve, wilsonSolve]

import strutils
import json

import monteCarlo
export monteCarlo

type
  FieldType* = enum
    GaugeField,
    StaggeredMatterField,
    WilsonMatterField
  GaugeActionType* = enum 
    Wilson, 
    Adjoint, 
    Rectangle, 
    Symanzik, 
    Iwasaki, 
    DBW2
  StaggeredActionType* = enum
    StaggeredFermion,
    StaggeredHasenbuschFermion,
    StaggeredBoson
  WilsonActionType* = enum
    WilsonFermion,
    WilsonBoson

type
  AbstractField*[S,T,U] = object of MonteCarloAtom    
    case field*: FieldType
      of GaugeField:
        gaugeAction*: GaugeActionType # Gauge action
        gaugeActionCoefficients*: GaugeActionCoeffs # Action coefficients
        u*: seq[S]
      of StaggeredMatterField:
        staggeredFields*: seq[T] # Pseudofermion/boson field
        staggeredMasses*: seq[float] # Masses
        staggeredAction*: StaggeredActionType # Staggered action
        staggeredActionSolverParameters*: SolverParams # Solver parameters
        staggeredForceSolverParameters*: SolverParams # Solver parameters
      of WilsonMatterField:
        wilsonFields*: seq[U] # Pseuodofermion/boson field
        wilsonMasses*: seq[float] # Masses
        wilsonAction*: WilsonActionType # Wilson action
        wilsonActionSolverParameters*: SolverParams # Solver parameters
        wilsonForceSolverParameters*: SolverParams # Solver parameters

proc newAbstractField*(self: var AbstractField; info: JsonNode) = 
  case self.algorithm:
    of HamiltonianMonteCarlo:
      if not info.hasKey("integrator"): self.integrator = MN2
      else: self.integrator = toIntegratorAtomType(info["integrator"].getStr())
      if not info.hasKey("steps"): qexError "Must specify number of integ. steps"
      else: self.steps = info["steps"].getInt()
      
      case self.integrator:
        of LF, LF1, MN0, LF2, LF3: discard
        of MN2: self.lmbda = 0.1931833275037836
        of MN4FP4:
          self.rho = 0.1786178958448091
          self.theta = -0.06626458266981843
          self.lmbda = 0.7123418310626056
        of MN4FP5:
          self.rho = 0.2750081212332419
          self.theta = -0.1347950099106792
          self.vartheta = -0.08442961950707149
          self.lmbda = 0.3549000571574260
        of MN4FV5:
          self.rho = 0.2539785108410595
          self.theta = -0.03230286765269967
          self.vartheta = 0.08398315262876693
          self.lmbda = 0.6822365335719091
        of CustomIntegrator: discard # Need to add in option

      case self.integrator:
        of CustomIntegrator: discard
        else:
          for key in info.keys():
            case key:
              of "lambda": self.lmbda = info[key].getFloat()
              of "rho": self.rho = info[key].getFloat()
              of "theta": self.theta = info[key].getFloat()
              of "vartheta": self.vartheta = info[key].getFloat()

      self.dtau = newSeq[float](self.stepsT.len)  
    of HeatbathOverrelax:
      case self.field:
        of GaugeField: discard
        else: qexError "Heatbath only supported/possible for gauge fields"
  
