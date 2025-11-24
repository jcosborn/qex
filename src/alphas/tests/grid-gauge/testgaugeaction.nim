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

proc gaugeAction2*(
  action: GridPlaqPlusRectangleAction; 
  u: GridLatticeGaugeField
): cdouble {.importcpp: "#.S(#)", altgauge.}

proc gaugeForce2*(
  action: GridOneLoopGaugeAction; 
  u: GridLatticeGaugeField; 
  force: GridLatticeGaugeField
) {.importcpp: "#.deriv(#, #)", gauge.}

proc gaugeForce2*(
  action: GridPlaqPlusRectangleAction; 
  u: GridLatticeGaugeField; 
  force: GridLatticeGaugeField
) {.importcpp: "#.deriv(#, #)", gauge.}

proc checkerboard(x: ptr GridLatticeGaugeField): int 
  {.importcpp: "#->Checkerboard()", grid.}

proc `:=`(r0: seq[Field], x0: GridLatticeGaugeField) =
  type GridScalarObject = GridLatticeGaugeField.scalarObj
  let 
    r = addr r0
    x = addr x0
  let
    lo = r0[0].l
    nd = lo.nDim
    nc = r0[0][0].ncols
    nSites = lo.nSites
  var subset = lo.getSubset("all")
  let glSites = x0.Grid.lSites
  var c0 = newSeq[cint](nd)

  lo.coord(c0, lo.myrank, 0)
  if glSites != nSites:
    subset = case x.checkerboard == 0
      of true: lo.getSubset("even")
      of false: lo.getSubset("odd")

  threads:
    {.emit: "using namespace Grid;".}
    {.emit: "Coordinate c(`nd`);".}
    {.emit: ["autoView(dst, ", x[], ", CpuRead);"].}
    
    var t: GridScalarObject
    
    for s in subset.singleSites:
      # set coordinate
      for mu in 0..<nd: {.emit: ["c[", mu, "] = ", lo.coords[mu][s].cint - c0[mu], ";"].}
      {.emit: "peekLocalSite(`t`, dst, c);".}

      # set link values
      for mu in 0..<nd:
        for a in 0..<nc:
          for b in 0..<nc:
            var tr, ti {.noinit.}: float
            {.emit: "`tr` = `t`._internal[`mu`]._internal._internal[`a`][`b`].real();".}
            {.emit: "`ti` = `t`._internal[`mu`]._internal._internal[`a`][`b`].imag();".}
            r[][mu]{s}[a,b] := newComplex(tr, ti)

when isMainModule:
  qexInit()

  defaultSetup()

  let (beta, cp, cr, cpr) = (1.0, 1.0, 1.0, 0.0)

  type TestType = enum ActionTest, ForceTest

  proc testGrid(test: TestType) =
    let latSize = newCoordinate(lat)
    let
      simdLayout = GridDefaultSimd(len(lat), Nsimd(GridVComplex))
      mpiLayout = newCoordinate(lo.rankGeom)
    let grid = latSize.newGridCartesian(simdLayout, mpiLayout)

    var qexGauge = lo.newGauge()
    var gridGauge = grid.gauge()

    let gridAction = newGridOneLoopGaugeAction(addr grid, beta, cp, cr, cpr)
    let altGridAction = newGridPlaqPlusRectangleAction(beta*cp, beta*cr)
    let qexAction = GaugeActionCoeffs(plaq: beta*cp, rect: beta*cr, pgm: beta*cpr)

    case test:
    of ActionTest:
      qexGauge.random()
      gridGauge := qexGauge

      var unit = lo.newGauge()
      var gridUnit = grid.gauge()
      for mu in 0..<unit.len: unit[mu] := 1
      gridUnit := unit

      let gridResult = gridAction.gaugeAction2(gridGauge)
      let altGridResult = altGridAction.gaugeAction2(gridGauge)
      let qexResult = qexAction.gaugeAction2(qexGauge) - qexAction.gaugeAction2(unit)

      echo "Grid: ", gridResult
      echo "Grid (alt, no parallelogram): ", altGridResult
      echo "QEX:  ", qexResult
      echo "ratio: ", gridResult / qexResult
      echo "diff: ", gridResult - qexResult
    of ForceTest:
      qexGauge.random()
      gridGauge := qexGauge

      var gridForce = grid.gauge()
      var gridForce2 = grid.gauge()
      var gridForceQEX = lo.newGauge()
      var gridForceQEX2 = lo.newGauge()
      var qexForce = lo.newGauge()
      gridAction.gaugeForce2(gridGauge, gridForce)
      altGridAction.gaugeForce2(gridGauge, gridForce2)
      qexAction.gaugeForce2(qexGauge, qexForce)

      gridForceQEX := gridForce
      gridForceQEX2 := gridForce2

      let gridResult = qexAction.gaugeAction2(gridForceQEX)
      let altGridResult = qexAction.gaugeAction2(gridForceQEX2)
      let qexResult = qexAction.gaugeAction2(qexForce)

      echo "Grid: ", gridResult
      #echo "Grid (alt, no parallelogram): ", altGridResult
      echo "QEX:  ", qexResult
      echo "ratio: ", gridResult / qexResult
      echo "diff: ", gridResult - qexResult

  testGrid(ForceTest)

  qexFinalize()