## Brief: HISQ using Grid
## 
## Author: Curtis Taylor Peterson <curtistaylorpetersonwork@gmail.com>
##
## Details:
## This module wraps Curtis Taylor Peterson's implementation of HISQ in Grid

import qex
import grid/[Grid]
import alphaslinks
import alphasproject

import gauge/[gaugeAction] # for debugging

export alphaslinks
export alphasproject

{.pragma: grid, header: "<Grid/Grid.h>".}
{.pragma: hisq, header: "<Grid/qcd/utils/HighlyImprovedStaggeredFermionImpl.h>".}

const HISQIMPL = "Grid::HighlyImprovedStaggeredFermionImpl"
const GIMPL = "Grid::PeriodicGimplR"

# some Grid/C++ types

type vector*[T] {.importcpp: "std::vector<'*0>", header:"<vector>".} = object

type GridHISQ {.importcpp: HISQIMPL & "<" & GIMPL & ">", hisq.} = object

proc newVector*[T](x,y: ptr T): vector[T] 
  {.importcpp: "std::vector<'*0>(#, #)", constructor.}
proc toVector*[T](x: openArray[T]): vector[T] =
  let x0 = unsafeAddr x[0]
  return newVector(x0, x0 + x.len)

proc newGridHISQ(grid: ptr GridCartesian, calcStagPhases: bool = true): GridHISQ 
  {.importcpp: HISQIMPL & "<" & GIMPL & ">(#, #)", constructor, hisq.}

# Helper Grid procedures

proc checkerboard(x: ptr GridLatticeGaugeField): int 
  {.importcpp: "#->Checkerboard()", grid.}

#{.emit:["SchurRedBlackStaggeredSolve<FermionField> SchurSolver(CG);"].}
#{.emit:["SchurSolver(Ds,",gsrc,",",gsoln,");"].}

proc toQEX(r0: seq[Field], x0: GridLatticeGaugeField) =
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

proc toGrid(r: var GridLatticeGaugeField, x0: openArray[Field]) = r := x0

# Grid HISQ method wrappers

proc rephase(self: GridHISQ; X: GridLatticeGaugeField; W: GridLatticeGaugeField) 
  {.importcpp: "#.rephase(#, #)", hisq.}

proc smear(
  self: GridHISQ;
  X: GridLatticeGaugeField; 
  WWW: GridLatticeGaugeField;
  W: GridLatticeGaugeField
) {.importcpp: "#.smear(#, #, #)", hisq.}

proc smear(
  self: GridHISQ;
  X: GridLatticeGaugeField; 
  W: GridLatticeGaugeField
) {.importcpp: "#.smear(#, #)", hisq.}

proc project(
  self: GridHISQ;
  V: GridLatticeGaugeField; 
  U: GridLatticeGaugeField;
  forDerivative: bool
) {.importcpp: "#.project(#, #, #)", hisq.}

proc smearDerivative(
  self: GridHISQ;
  dXdU: GridLatticeGaugeField; 
  dXdW: GridLatticeGaugeField; 
  W: GridLatticeGaugeField
) {.importcpp: "#.smearDerivative(#, #, #)", hisq.} 

proc smearDerivative(
  self: GridHISQ;
  dXdU: GridLatticeGaugeField; 
  dXdW: GridLatticeGaugeField;
  dXdWWW: GridLatticeGaugeField;
  W: GridLatticeGaugeField
) {.importcpp: "#.smearDerivative(#, #, #, #)", hisq.} 

proc projectionDerivative(
  self: GridHISQ;
  dVdU: GridLatticeGaugeField; 
  dZdV: GridLatticeGaugeField; 
  V: GridLatticeGaugeField;
  U: GridLatticeGaugeField
) {.importcpp: "#.projectionDerivative(#, #, #, #)", hisq.}

# primary smearing and force derivative procedure

proc smearGetForce*[T](
  self: alphaslinks.HisqCoefs; 
  u: T; 
  su, sul: T;
  displayPerformance: bool = false,
  regulate: bool = false
): proc(dsdu: var T; dsdsu, dsdsul: T) =
  var lo = u[0].l
  var 
    g = lo.newGauge()
    v = lo.newGauge()
    w = lo.newGauge()

  block:
    # grid prep
    let
      lat = lo.physGeom
      latSize = newCoordinate(lat)
    let
      simdLayout = GridDefaultSimd(len(lat), Nsimd(GridVComplex))
      mpiLayout = newCoordinate(lo.rankGeom)
    let grid = latSize.newGridCartesian(simdLayout, mpiLayout)
    var hisq = newGridHISQ(addr grid)

    # force
    var
      gu = grid.gauge()
      gw = grid.gauge()
      gv = grid.gauge()

    # output
    var gsu = grid.gauge()
    var gsul = grid.gauge()

    # load inputs & rephase
    gu.toGrid(u)
    hisq.rephase(gu, gu)

    # smear
    hisq.smear(gv, gu)
    hisq.project(gw, gv, regulate)
    hisq.smear(gsu, gsul, gw)

    # save for force
    g.toQEX(gu)
    v.toQEX(gv)
    w.toQEX(gw)

    # save output
    su.toQEX(gsu)
    sul.toQEX(gsul)

    #[]
    for mu in 0..<lo.nDim:
      {.emit: "using namespace Grid;".}
      {.emit: "std::cout << \"RESULT OF SMEAR: \" << `mu` << \" \" << sum(trace(PeekIndex<LorentzIndex>(`gu`, `mu`))) << std::endl;".}
      echo "RESULT OF SMEAR: ", mu, " ", simdSum(trace(g[mu]))
      echo ""
      {.emit: "std::cout << \"RESULT OF SMEAR: \" << `mu` << \" \" << sum(trace(PeekIndex<LorentzIndex>(`gv`, `mu`))) << std::endl;".}
      echo "RESULT OF SMEAR: ", mu, " ", simdSum(trace(v[mu]))
      echo ""
      {.emit: "std::cout << \"RESULT OF SMEAR: \" << `mu` << \" \" << sum(trace(PeekIndex<LorentzIndex>(`gw`, `mu`))) << std::endl;".}
      echo "RESULT OF SMEAR: ", mu, " ", simdSum(trace(w[mu]))
      echo ""
      {.emit: "std::cout << \"RESULT OF SMEAR: \" << `mu` << \" \" << sum(trace(PeekIndex<LorentzIndex>(`gsu`, `mu`))) << std::endl;".}
      echo "RESULT OF SMEAR: ", mu, " ", simdSum(trace(su[mu]))
      echo ""
      {.emit: "std::cout << \"RESULT OF SMEAR: \" << `mu` << \" \" << sum(trace(PeekIndex<LorentzIndex>(`gsul`, `mu`))) << std::endl;".}
      echo "RESULT OF SMEAR: ", mu, " ", simdSum(trace(sul[mu]))
      echo ""
    ]#

    #[
    var gc = GaugeActionCoeffs(plaq: 1.0, rect: 1.0)
    echo "action from smeared link: ", gc.gaugeAction2(su)
    for mu in 0..<lo.nDim:
      {.emit: "using namespace Grid;".}
      {.emit: "std::cout << \"trace sum from smeared link (Grid): \" << `mu` << \" \" << sum(trace(PeekIndex<LorentzIndex>(`gsu`, `mu`))) << std::endl;".}
      echo "trace sum from smeared link (QEX): ", simdSum(trace(su[mu]))
    echo "action from long link: ", gc.gaugeAction2(sul)
    for mu in 0..<lo.nDim:
      {.emit: "using namespace Grid;".}
      {.emit: "std::cout << \"trace sum from long link (Grid): \" << `mu` << \" \" << sum(trace(PeekIndex<LorentzIndex>(`gsul`, `mu`))) << std::endl;".}
      echo "trace sum from long link (QEX): ", simdSum(trace(sul[mu]))
  ]#

  return proc(dsdu: var T; dsdsu, dsdsul: T) =
    # grid prep
    let
      lat = lo.physGeom
      latSize = newCoordinate(lat)
    let
      simdLayout = GridDefaultSimd(len(lat), Nsimd(GridVComplex))
      mpiLayout = newCoordinate(lo.rankGeom)
    let grid = latSize.newGridCartesian(simdLayout, mpiLayout)
    var hisq = newGridHISQ(addr grid)

    # various smearing levels
    var
      gu = grid.gauge()
      gv = grid.gauge()
      gw = grid.gauge()

    # chain rule
    var 
      gdsdu = grid.gauge()
      gdsdsu = grid.gauge()
      gdsdsul = grid.gauge()
    
     # temporary
    var gt = grid.gauge()

    # convert to Grid layout
    gu.toGrid(g)
    gv.toGrid(v)
    gw.toGrid(w)
    gdsdsu.toGrid(dsdsu)
    gdsdsul.toGrid(dsdsul)

    # derivative
    hisq.smearDerivative(gt, gdsdsu, gdsdsul, gw)
    hisq.projectionDerivative(gt, gt, gw, gv)
    hisq.smearDerivative(gdsdu, gt, gu)

    # save result
    dsdu.toQEX(gdsdu)

when isMainModule:
  qexInit()

  defaultSetup()

  let hisq = newHISQ()

  g.random()
  var 
    sg = g.newOneOf
    sgl = g.newOneOf
  var
    c1 = g.newOneOf
    c2 = g.newOneOf
    f = g.newOneOf
  c1.random()
  c2.random()

  var force = hisq.smearGetForce(g, sg, sgl)
  force(f, c1, c2)

  qexFinalize()