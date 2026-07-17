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

proc evalProjectedGaugeForceValue*(gc: GaugeActionCoeffs, g, outg: shared.Gauge) =
  case gc.gaugeActionFamily
  of gafGaugeAction1:
    gc.gaugeForce(g, outg)
  of gafActionA:
    gc.forceA(g, outg)

proc evalGaugeForceSubset*(gc: GaugeActionCoeffs, g, outg: shared.Gauge, sd, sf, sb: auto, parity, dir: int) =
  ## Subset Wilson derivative; reject non-plaquette coefficients.
  if gc.rect != 0 or gc.pgm != 0 or gc.adjplaq != 0:
    raiseUnsupportedGaugeCoeff(gc)
  gc.gaugeDeriv2SubsetWork(g, outg, sd, sf, sb, parity, dir, clear=false)

proc evalGaugeForceJacobian*(b: shared.Gauge,
                             gc: GaugeActionCoeffs,
                             g: shared.Gauge,
                             outg: shared.Gauge) =
  case gc.gaugeActionFamily
  of gafGaugeAction1:
    outg.zeroGaugeStorage
    gc.gaugeDerivDeriv2(g, b, outg)
  of gafActionA:
    raiseUnsupportedPath("evalGaugeForceJacobian", "ActionA-family second derivatives")

proc evalGaugeForceJacobianSubset*(b: shared.Gauge, gc: GaugeActionCoeffs, g, outg: shared.Gauge, parity, dir: int) =
  ## outg = H_g(b[dir]|parity); writes all neighbouring links. Plaquette-only.
  if gc.rect != 0 or gc.pgm != 0 or gc.adjplaq != 0:
    raiseUnsupportedGaugeCoeff(gc)
  gc.gaugeDerivDeriv2Subset(g, b, outg, parity, dir)

proc evalGaugeForceJacobianSubsetSum*(b: seq[DLatticeColorMatrixV], w: DLatticeColorMatrixV, gc: GaugeActionCoeffs, g, outg: shared.Gauge, parity, dir: int) =
  if gc.rect != 0 or gc.pgm != 0 or gc.adjplaq != 0:
    raiseUnsupportedGaugeCoeff(gc)
  gc.gaugeDerivDeriv2SubsetSum(g, b, w, outg, parity, dir)
