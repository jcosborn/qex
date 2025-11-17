import qex
import grid/[Grid]
import gauge/[gaugeAction]

{.pragma: gauge, header: "<Grid/qcd/action/gauge/OneLoopGaugeAction.h>".}

const ONELOOPACTION = "Grid::OneLoopGaugeAction<Grid::PeriodicGimplR>"

type OneLoopGaugeAction* {.importcpp: ONELOOPACTION, gauge.} = object

proc newOneLoopGaugeAction(
  grid: ptr GridCartesian, 
  beta, cp, cr, cpr: cdouble
): OneLoopGaugeAction 
  {.importcpp: ONELOOPACTION & "(#, #, #, #, #)", constructor, gauge.}

proc S(action: OneLoopGaugeAction; u: GridLatticeGaugeField): cdouble 
  {.importcpp: "#.S(#)", gauge.}

proc gaugeActionOneLoopHISQ*[U](gc: GaugeActionCoeffs; u: openArray[U]): cdouble =
  # prepare grid
  let lo = u[0].l
  let
    lat = lo.physGeom
    latSize = newCoordinate(lat)
  let
    simdLayout = GridDefaultSimd(len(lat), Nsimd(GridVComplex))
    mpiLayout = newCoordinate(lo.rankGeom)
  let grid = latSize.newGridCartesian(simdLayout, mpiLayout)

  # get grid link and action object
  var g = grid.gauge()
  let action = gc.newOneLoopGaugeAction(addr grid, gc.plaq, gc.plaq, gc.rect, gc.pgm)

  # do action calculation and return
  g := u
  return action.S(g)
  
