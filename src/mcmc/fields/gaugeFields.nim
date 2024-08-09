import ../mcmcTypes

import gauge

import strutils
import json

const
  BetaAOverBetaF = -1.0/4.0
  C1Symanzik = -1.0/12.0
  C1Iwasaki = -0.331
  C1DBW2 = -1.4088

converter toGaugeActionType(s: string):
  GaugeActionType = parseEnum[GaugeActionType](s)

proc checkJSON(info: JsonNode): JsonNode = 
  result = parseJson("{}")
  for key, keyVal in info: result[key] = keyVal
  if not result.hasKey("adjoint-ratio"):
    result["adjoint-ratio"] = %* BetaAOverBetaF
  if not result.hasKey("rectangle-coefficient"):
    result["rectangle-coefficient"] = %* C1Symanzik
  if not result.hasKey("action"): result["action"] = %* "Wilson"
  if not result.hasKey("beta"): qexError "beta not specified for gauge field"
  result["field-type"] = %* "gauge"

proc newGaugeField(
    self: var LatticeField;
    action: GaugeActionType;
    beta, adjRat, rectCoeff: float
  ) =
  # Create new stream
  var stream = newMCStream("new gauge field")
  
  # Initialization specific to gauge fields
  self.gaugeAction = action
  self.gaugeActionCoefficients = case self.gaugeAction
    of Wilson: GaugeActionCoeffs(plaq: beta)
    of Adjoint: GaugeActionCoeffs(plaq: beta, adjplaq: beta*adjRat)
    of Rectangle: gaugeActRect(beta, rectCoeff)
    else: # Annoying Nim 2.0 compiler workaround
      if self.gaugeAction == Symanzik: gaugeActRect(beta, C1Symanzik)
      elif self.gaugeAction == Iwasaki: gaugeActRect(beta, C1Iwasaki)
      else: gaugeActRect(beta, C1DBW2)

  stream.add "  action = " & $(action)
  stream.add "  beta = " & $(beta)
  case self.gaugeAction:
    of Adjoint: stream.add "  beta_a/beta_f = " & $(adjRat)
    of Rectangle: stream.add "  c1 = " & $(rectCoeff)
    else: discard
  stream.finishStream

proc newGaugeField*(l: Layout; gaugeInformation: JsonNode): auto =
  let info = checkJSON(gaugeInformation)
  result = l.newLatticeField(info)
  result.u = l.newGauge()
  result.newGaugeField(
    toGaugeActionType(info["action"].getStr()),
    info["beta"].getFloat(),
    info["adjoint-ratio"].getFloat(),
    info["rectangle-coefficient"].getFloat(),
  )

proc gaugeAction*(self: LatticeField; u: auto): float =
  result = case self.gaugeAction
    of Adjoint: self.gaugeActionCoefficients.actionA(u)
    else: self.gaugeActionCoefficients.gaugeAction1(u)

proc gaugeForce*[S](self: LatticeField; u,f: seq[S]) =
  case self.gaugeAction:
    of Adjoint: self.gaugeActionCoefficients.forceA(u,f)
    else: self.gaugeActionCoefficients.gaugeForce(u,f)

proc gaugeAction*(self: LatticeField): float =
  result = case self.gaugeAction
    of Adjoint: self.gaugeActionCoefficients.actionA(self.u)
    else: self.gaugeActionCoefficients.gaugeAction1(self.u)

proc gaugeForce*[S](self: LatticeField; f: seq[S]) =
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
      "monte-carlo-algorithm": "hmc"
    }

    case action:
      of "Adjoint": params["adjoint-ratio"] = %* BetaAOverBetaF
      of "Rectangle": params["rectangle-coefficient"] = %* C1Symanzik
      else: discard

    var gauge = lo.newGaugeField(params)

    unit(u)
    gauge.gaugeForce(f)
    var f2: float
    threads:
      var f2t = 0.0
      for mu in 0..<f.len: f2t += f[mu].norm2
      threadBarrier()
      threadMaster: f2 = f2t
    echo action, ": ", gauge.gaugeAction(), " ", f2

  qexFinalize()