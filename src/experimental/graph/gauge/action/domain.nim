from ../../core/base import raiseValueError
from ../../support/op import raiseUnsupportedPath
import ../../../../layout, ../../../../gauge
import ../../../../physics/qcdTypes
import ../types

type
  GaugeActionFamily = enum
    gafGaugeAction1, gafActionA

const
  C1Symanzik* = -1.0/12.0
  C1Iwasaki* = -0.331
  C1DBW2* = -1.4088

proc raiseUnsupportedGaugeCoeff(gc: GaugeActionCoeffs) {.noreturn.} =
  raiseValueError("Gauge coefficient unsupported: " & $gc)

proc isPlaqOnly*(gc: GaugeActionCoeffs): bool =
  gc.rect == 0 and gc.pgm == 0 and gc.adjplaq == 0

proc isPlaqRect*(gc: GaugeActionCoeffs): bool =
  gc.pgm == 0 and gc.adjplaq == 0

proc isAdjPlaq*(gc: GaugeActionCoeffs): bool =
  ## The kernel family with adjoint plaquettes (actionA), selected by adjplaq != 0.
  gc.rect == 0 and gc.pgm == 0 and gc.adjplaq != 0

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

proc evalGaugeActionValue*(gc: GaugeActionCoeffs, g: types.Gauge): float =
  case gc.gaugeActionFamily
  of gafGaugeAction1:
    # gaugeAction1 and gaugeActionDeriv carry no parallelogram terms.
    if gc.pgm != 0:
      raiseUnsupportedPath("gaugeAction", "parallelogram coefficients")
    gc.gaugeAction1 g
  of gafActionA:
    gc.actionA g

proc evalGaugeForceValue*(gc: GaugeActionCoeffs,
                          g: types.Gauge,
                          outg: types.Gauge) =
  let coeffs = gc.negatedGaugeCoeffs
  case coeffs.gaugeActionFamily
  of gafGaugeAction1:
    if gc.pgm != 0:
      raiseUnsupportedPath("gaugeActionDeriv", "parallelogram coefficients")
    coeffs.gaugeActionDeriv(g, outg)
  of gafActionA:
    coeffs.gaugeADeriv(g, outg)

proc evalProjectedGaugeForceValue*(gc: GaugeActionCoeffs, g, outg: types.Gauge) =
  case gc.gaugeActionFamily
  of gafGaugeAction1:
    gc.gaugeForce(g, outg)
  of gafActionA:
    gc.forceA(g, outg)

proc evalGaugeForceSubset*(gc: GaugeActionCoeffs, g, outg: types.Gauge, sd, sf, sb: auto, parity, dir: int) =
  ## Subset Wilson derivative; reject non-plaquette coefficients.
  if not gc.isPlaqOnly:
    raiseUnsupportedGaugeCoeff(gc)
  gc.gaugeDeriv2SubsetWork(g, outg, sd, sf, sb, parity, dir, clear=false)

proc evalGaugeForceJacobian*(b: types.Gauge,
                             gc: GaugeActionCoeffs,
                             g: types.Gauge,
                             outg: types.Gauge) =
  case gc.gaugeActionFamily
  of gafGaugeAction1:
    # gaugeDerivDeriv2 has plaquette terms only; gaugeActionGraph differentiates
    # the rectangle family through basic ops instead.
    if not gc.isPlaqOnly:
      raiseUnsupportedPath("evalGaugeForceJacobian", "rectangle and parallelogram Hessians")
    outg.zeroGaugeStorage
    gc.gaugeDerivDeriv2(g, b, outg)
  of gafActionA:
    raiseUnsupportedPath("evalGaugeForceJacobian", "ActionA-family second derivatives")

proc evalGaugeForceJacobianSubset*(b: types.Gauge, gc: GaugeActionCoeffs, g, outg: types.Gauge, parity, dir: int) =
  ## outg = H_g(b[dir]|parity); writes all neighbouring links. Plaquette-only.
  if not gc.isPlaqOnly:
    raiseUnsupportedGaugeCoeff(gc)
  gc.gaugeDerivDeriv2Subset(g, b, outg, parity, dir)

proc evalGaugeForceJacobianSubsetSum*(b: seq[DLatticeColorMatrixV], w: DLatticeColorMatrixV, gc: GaugeActionCoeffs, g, outg: types.Gauge, parity, dir: int) =
  if not gc.isPlaqOnly:
    raiseUnsupportedGaugeCoeff(gc)
  gc.gaugeDerivDeriv2SubsetSum(g, b, w, outg, parity, dir)
