## Brief: HISQ using Grid
## 
## Author: Curtis Taylor Peterson <curtistaylorpetersonwork@gmail.com>
##
## Details:
## This module wraps Curtis Taylor Peterson's implementation of HISQ in Grid

import qex
import grid/[Grid]
import alphaslinks

export alphaslinks

{.pragma: grid, header: "<Grid/Grid.h>".}
{.pragma: hisq, header: "HighlyImprovedStaggeredFermionImpl.h".}

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
    var t: GridScalarObject
    {.emit: "Coordinate c(`nd`);".}
    {.emit: ["autoView(dst, ", x[], ", CpuRead);"].}
    for s in subset.singleSites:
      for mu in 0..<nd: 
        let smu = lo.coords[mu][s].cint - c0[mu]
        {.emit: "c[`mu`] = `smu`;".}
      {.emit: "peekLocalSite(`t`, dst, c);".}
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
  W: GridLatticeGaugeField; 
  coeffs: vector[float]; 
  lepage, naik: float
) {.importcpp: "#.smear(#, #, #, #, #, #)", hisq.}

proc smear(
  self: GridHISQ;
  X: GridLatticeGaugeField; 
  W: GridLatticeGaugeField; 
  coeffs: vector[float]; 
  lepage: float
) {.importcpp: "#.smear(#, #, #, #)", hisq.}

proc project(
  self: GridHISQ;
  V: GridLatticeGaugeField; 
  U: GridLatticeGaugeField;
  derivativeCutoff: float;
  backupSVD: bool;
  svdTolerance: float
) {.importcpp: "#.project(#, #, #, #, #)", hisq.}

proc smearDerivative(
  self: GridHISQ;
  dXdU: GridLatticeGaugeField; 
  dXdW: GridLatticeGaugeField; 
  W: GridLatticeGaugeField; 
  coeffs: vector[float]; 
  lepage: float
) {.importcpp: "#.smearDerivative(#, #, #, #, #)", hisq.} 

proc smearDerivative(
  self: GridHISQ;
  dXdU: GridLatticeGaugeField; 
  dXdW: GridLatticeGaugeField;
  dXdWWW: GridLatticeGaugeField;
  W: GridLatticeGaugeField; 
  coeffs: vector[float];
  lepage, naik: float;
) {.importcpp: "#.smearDerivative(#, #, #, #, #, #, #)", hisq.} 

proc projectionDerivative(
  self: GridHISQ;
  dVdU: GridLatticeGaugeField; 
  dZdV: GridLatticeGaugeField; 
  V: GridLatticeGaugeField;
  U: GridLatticeGaugeField
) {.importcpp: "#.projectionDerivative(#, #, #, #)", hisq.}

proc smearGetForce*[T](
  self: alphaslinks.HisqCoefs; 
  u: T; 
  su, sul: T;
  displayPerformance: bool = false,
  regulate: bool = false
): proc(dsdu: var T; dsdsu, dsdsul: T) =
  let
    f7l1 = @[
      self.fat7first.oneLink,
      self.fat7first.threeStaple,
      self.fat7first.fiveStaple,
      self.fat7first.sevenStaple
    ].toVector()
    f7l2 = @[
      self.fat7second.oneLink,
      self.fat7second.threeStaple,
      self.fat7second.fiveStaple,
      self.fat7second.sevenStaple
    ].toVector()
    lpl1 = self.fat7first.lepage
    lpl2 = self.fat7second.lepage
    naik = self.naik
  var lo = u[0].l
  var 
    g = lo.newGauge()
    v = lo.newGauge()
    w = lo.newGauge()

  block:
    let
      lat = lo.physGeom
      latSize = newCoordinate(lat)
      simdLayout = GridDefaultSimd(len(lat), Nsimd(GridVComplex))
      mpiLayout = newCoordinate(lo.rankGeom)
    let grid = latSize.newGridCartesian(simdLayout, mpiLayout)
    var hisq = newGridHISQ(addr grid)
    var
      gu = grid.gauge()
      gw = grid.gauge()
      gv = grid.gauge()
    var
      gsu = grid.gauge()
      gsul = grid.gauge()

    gu.toGrid(u)
    hisq.rephase(gu, gu)

    # smear
    hisq.smear(gv, gu, f7l1, lpl1)
    case regulate: # ensure that regulation only done for force action
      of true: hisq.project(gw, gv, 5e-5, true, 1e-8)
      of false: hisq.project(gw, gv, 1e-20, true, 1e-8)
    hisq.smear(gsu, gsul, gw, f7l2, lpl2, naik)

    # save results
    g.toQEX(gu)
    v.toQEX(gv)
    w.toQEX(gw)
    su.toQEX(gsu)
    sul.toQEX(gsul)

  proc smearedForce(dsdu: var T; dsdsu, dsdsul: T) =
    let
      lat = lo.physGeom
      latSize = newCoordinate(lat)
      simdLayout = GridDefaultSimd(len(lat), Nsimd(GridVComplex))
      mpiLayout = newCoordinate(lo.rankGeom)
    let grid = latSize.newGridCartesian(simdLayout, mpiLayout)
    var hisq = newGridHISQ(addr grid)
    var
      gu = grid.gauge()
      gv = grid.gauge()
      gw = grid.gauge()
    var 
      gdsdu = grid.gauge()
      gdsdsu = grid.gauge()
      gdsdsul = grid.gauge()
    var gt = grid.gauge()

    gu.toGrid(g)
    gv.toGrid(v)
    gw.toGrid(w)
    gdsdsu.toGrid(dsdsu)
    gdsdsul.toGrid(dsdsul)

    hisq.smearDerivative(gt, gdsdsu, gdsdsul, gw, f7l2, lpl2, naik)
    hisq.projectionDerivative(gt, gt, gw, gv)
    hisq.smearDerivative(gdsdu, gt, gu, f7l1, lpl1)

    dsdu.toQEX(gdsdu)

    qexError "intentional exit"

  return smearedForce

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