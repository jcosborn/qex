from ../../core/base import raiseValueError
from ../../support/op import raiseUnsupportedPath
import ../../../../layout, ../../../../gauge
import ../../../../physics/qcdTypes
import ../shared

type
  GaugeActionFamily = enum
    gafGaugeAction1, gafActionA

const
  C1Symanzik* = -1.0/12.0
  C1Iwasaki* = -0.331
  C1DBW2* = -1.4088

proc raiseUnsupportedGaugeCoeff(gc: GaugeActionCoeffs) {.noreturn.} =
  raiseValueError("Gauge coefficient unsupported: " & $gc)

proc gaugeActionFamily(gc: GaugeActionCoeffs): GaugeActionFamily =
  if gc.adjplaq == 0:
    return gafGaugeAction1
  if gc.rect == 0 and gc.pgm == 0:
    return gafActionA
  raiseUnsupportedGaugeCoeff(gc)

proc negatedGaugeCoeffs(gc: GaugeActionCoeffs): GaugeActionCoeffs =
  result = gc
  for f in result.fields:
    f = -f

proc evalGaugeActionValue*(gc: GaugeActionCoeffs, g: shared.Gauge): float =
  case gc.gaugeActionFamily
  of gafGaugeAction1:
    gc.gaugeAction1 g
  of gafActionA:
    gc.actionA g

proc evalGaugeForceValue*(gc: GaugeActionCoeffs,
                          g: shared.Gauge,
                          outg: shared.Gauge) =
  let coeffs = gc.negatedGaugeCoeffs
  case coeffs.gaugeActionFamily
  of gafGaugeAction1:
    coeffs.gaugeActionDeriv(g, outg)
  of gafActionA:
    coeffs.gaugeADeriv(g, outg)

proc evalGaugeForceJacobian*(b: shared.Gauge,
                             gc: GaugeActionCoeffs,
                             g: shared.Gauge,
                             outg: shared.Gauge) =
  case gc.gaugeActionFamily
  of gafGaugeAction1:
    threads:
      for mu in 0..<outg.len:
        outg[mu] := 0.0
    gc.gaugeDerivDeriv2(g, b, outg)
  of gafActionA:
    raiseUnsupportedPath("evalGaugeForceJacobian", "ActionA-family second derivatives")
