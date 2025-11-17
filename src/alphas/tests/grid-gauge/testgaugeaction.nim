import qex
import grid/[Grid]
import gauge/[gaugeAction]

{.pragma: grid, header: "<Grid/Grid.h>".}
{.pragma: gauge, header: "<Grid/qcd/action/gauge/OneLoopGaugeAction.h>".}
{.pragma: altgauge, header: "<Grid/qcd/action/gauge/PlaqPlusRectangleAction.h>".}

const GRIDONELOOPACTION = "Grid::OneLoopGaugeAction<Grid::PeriodicGimplR>"
const ALTGAUGEACTION = "Grid::PlaqPlusRectangleAction<Grid::PeriodicGimplR>"

type
  GridOneLoopGaugeAction* {.importcpp: GRIDONELOOPACTION, gauge.} = object

type
  GridPlaqPlusRectangleAction* {.importcpp: ALTGAUGEACTION, altgauge.} = object

proc newGridOneLoopGaugeAction*(
  grid: ptr GridCartesian, 
  beta, cp, cr, cpr: cdouble
): GridOneLoopGaugeAction
  {.importcpp: GRIDONELOOPACTION & "(#, #, #, #, #)", constructor, gauge.}

proc newGridPlaqPlusRectangleAction*(cp, cr: cdouble): GridPlaqPlusRectangleAction
  {.importcpp: ALTGAUGEACTION & "(#, #)", constructor, altgauge.}

proc gaugeAction2*(action: GridOneLoopGaugeAction; u: GridLatticeGaugeField): cdouble 
  {.importcpp: "#.S(#)", gauge.}

proc gaugeAction2*(action: GridPlaqPlusRectangleAction; u: GridLatticeGaugeField): cdouble 
  {.importcpp: "#.S(#)", altgauge.}

when isMainModule:
  qexInit()

  defaultSetup()

  let (beta, cp, cr, cpr) = (1.0, 1.0, 1.0, 1.0)

  proc testGrid() =
    let latSize = newCoordinate(lat)
    let
      simdLayout = GridDefaultSimd(len(lat), Nsimd(GridVComplex))
      mpiLayout = newCoordinate(lo.rankGeom)
    let grid = latSize.newGridCartesian(simdLayout, mpiLayout)

    var qexGauge = lo.newGauge()
    var gridGauge = grid.gauge()
    qexGauge.random()
    gridGauge := qexGauge

    var unit = lo.newGauge()
    var gridUnit = grid.gauge()
    for mu in 0..<unit.len: unit[mu] := 1
    gridUnit := unit

    let gridAction = newGridOneLoopGaugeAction(addr grid, beta, cp, cr, cpr)
    let gridResult = gridAction.gaugeAction2(gridGauge)

    discard gridAction.gaugeAction2(gridUnit)

    let altGridAction = newGridPlaqPlusRectangleAction(beta*cp, beta*cr)
    let altGridResult = altGridAction.gaugeAction2(gridGauge)

    let qexAction = GaugeActionCoeffs(plaq: beta*cp, rect: beta*cr, pgm: beta*cpr)
    let qexResult = qexAction.gaugeAction2(qexGauge) - qexAction.gaugeAction2(unit)

    echo "Grid: ", gridResult
    echo "Grid (alt, no parallelogram): ", altGridResult
    echo "QEX:  ", qexResult
    echo "ratio: ", gridResult / qexResult
    echo "diff: ", gridResult - qexResult

  testGrid()

  qexFinalize()