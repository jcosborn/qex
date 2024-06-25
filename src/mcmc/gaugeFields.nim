import qex
import gauge

import strutils
import json

import abstractFields
import typeUtilities

const
  BetaAOverBetaF = -1.0/4.0
  C1Symanzik = -1.0/12.0
  C1Iwasaki = -0.331
  C1DBW2 = -1.4088

converter toGaugeActionType(s: string):
  GaugeActionType = parseEnum[GaugeActionType](s)

proc checkJSON(info: JsonNode): JsonNode = 
  result = checkMonteCarloAlgorithm(info)
  if not result.hasKey("adjoint-ratio"):
    result["adjoint-ratio"] = %* BetaAOverBetaF
  if not result.hasKey("rectangle-coefficient"):
    result["rectangle-coefficient"] = %* C1Symanzik
  if not result.hasKey("action"): result["action"] = %* "Wilson"
  if not result.hasKey("beta"): qexError "beta not specified for gauge field"

proc newGaugeField[S,T,U](
    l: Layout;
    action: GaugeActionType;
    beta, adjRat, rectCoeff: float;
    algorithm: MonteCarloType;
    s: typedesc[S];
    t: typedesc[T];
    u: typedesc[U];
  ): AbstractField[S,T,U] =

  # Standard stuff that is necessary for all abstract fields
  result = AbstractField[S,T,U](field: GaugeField, algorithm: algorithm)
  
  # Initialization specific to gauge fields
  result.gaugeAction = action
  result.gaugeActionCoefficients = case result.gaugeAction
    of Wilson: GaugeActionCoeffs(plaq: beta)
    of Adjoint: GaugeActionCoeffs(plaq: beta, adjplaq: beta*adjRat)
    of Rectangle: gaugeActRect(beta, rectCoeff)
    else: # Annoying Nim 2.0 compiler workaround
      if result.gaugeAction == Symanzik: gaugeActRect(beta, C1Symanzik)
      elif result.gaugeAction == Iwasaki: gaugeActRect(beta, C1Iwasaki)
      else: gaugeActRect(beta, C1DBW2)

  result.u = l.newGauge()

proc newGaugeField*(l: Layout; gaugeInformation: JsonNode): auto =
  
  # Check JSON keys
  let info = checkJSON(gaugeInformation)

  # Create new gauge field
  result = l.newGaugeField(
    toGaugeActionType(info["action"].getStr()),
    info["beta"].getFloat(),
    info["adjoint-ratio"].getFloat(),
    info["rectangle-coefficient"].getFloat(),
    toMonteCarloType(info["monte-carlo-algorithm"].getStr()),
    l.typeS, l.typeT, l.typeU
  )
  
  # Set abstract field information
  result.newAbstractField(info)

proc gaugeAction*(self: AbstractField): float =
  result = case self.gaugeAction
    of Adjoint: self.gaugeActionCoefficients.actionA(self.u)
    else: self.gaugeActionCoefficients.gaugeAction1(self.u)

proc gaugeForce*[S](self: AbstractField; f: seq[S]) =
  case self.gaugeAction:
    of Adjoint: self.gaugeActionCoefficients.forceA(self.u, f)
    else: self.gaugeActionCoefficients.gaugeForce(self.u, f)

if isMainModule:
  qexInit()

  var 
    lat = intSeqParam("lat", @[4, 4, 4, 4])
    lo = lat.newLayout(@[1, 1, 1, 1])
    f = lo.newGauge()
    u = lo.newGauge()

  for action in ["Wilson", "Adjoint", "Rectangle", "Symanzik", "Iwasaki", "DBW2"]:
    var params = %* {
      "action": action,
      "beta": 6.0,
      "steps": 10,
      "integrator": "MN2",
      "monte-carlo-algorithm": "hamiltonian-monte-carlo"
    }

    case action:
      of "Adjoint": params["adjoint-ratio"] = %* BetaAOverBetaF
      of "Rectangle": params["rectangle-coefficient"] = %* C1Symanzik
      else: discard

    var gauge = lo.newGaugeField(params)

    unit(u)
    gauge.dtau = @[1.0]
    gauge.gaugeForce(f)
    var f2: float
    threads:
      var f2t = 0.0
      for mu in 0..<f.len: f2t += f[mu].norm2
      threadBarrier()
      threadMaster: f2 = f2t
    echo action, ": ", gauge.gaugeAction(), " ", f2

  qexFinalize()