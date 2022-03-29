import qex
import physics/stagSolve
import grid/gridImpl

template getGrid(g: Field): untyped =
  let lo = g.l
  let latt_size = newCoordinate(lo.physGeom)
  #let simd_layout = newCoordinate(lo.innerGeom)
  #let simd_layout = GridDefaultSimd(lo.nDim, GridVComplex.Nsimd);
  let simd_layout = GridDefaultSimd(lo.nDim, Nsimd(GridVComplex));
  let mpi_layout = newCoordinate(lo.rankGeom)
  let grid = newGridCartesian(latt_size,simd_layout,mpi_layout)

proc gridSolveEE*(s:Staggered; r,t:Field; m:SomeNumber; sp: var SolverParams) =
  let lo = r.l
  let latt_size = newCoordinate(lo.physGeom)
  #let simd_layout = newCoordinate(lo.innerGeom)
  #let simd_layout = GridDefaultSimd(lo.nDim, GridVComplex.Nsimd);
  let simd_layout = GridDefaultSimd(lo.nDim, Nsimd(GridVComplex));
  let mpi_layout = newCoordinate(lo.rankGeom)
  let grid = newGridCartesian(latt_size,simd_layout,mpi_layout)
  #r.getGrid
  #let grid = r.getGrid()
  #let rbgrid = newGridRedBlackCartesian(grid)
  var gfl = grid.gauge()
  var gll = grid.gauge()
  if s.g.len == 4: # plain staggered
    gfl := s.g[0..3]
    {.emit:"gll = Grid::Zero();".}
  elif s.g.len == 8: # Naik staggered
    gfl := @[s.g[0],s.g[2],s.g[4],s.g[6]]
    gll := @[s.g[1],s.g[3],s.g[5],s.g[7]]
  else:
    qexError "unknown s.g.len: ", s.g.len

  type ferm = GridImprovedStaggeredFermionR
  let rbgrid = newGridRedBlackCartesian(grid)
  var gsrc = rbgrid.fermion(ferm)
  var gsoln = rbgrid.fermion(ferm)
  gsrc.even
  gsoln.even
  gsrc := t
  var mass = m
  var res = sqrt sp.r2req
  var maxit = sp.maxits
  {.emit:"using namespace Grid;".}
  {.emit:"gsoln = Zero();".}
  {.emit:"using ImpStag = ImprovedStaggeredFermionR;".}
  {.emit:"using FermionField = ImpStag::FermionField;".}
  {.emit:"ImpStag Ds(grid,rbgrid,2.*mass,2.,2.,1.);".}
  {.emit:"Ds.ImportGaugeSimple(gll,gfl);".}
  {.emit:"SchurStaggeredOperator<ImpStag,FermionField> HermOp(Ds);".}
  {.emit:"ConjugateGradient<FermionField> CG(res, maxit, false);".}
  {.emit:"CG(HermOp, gsrc, gsoln);".}
  var rr = r
  rr := gsoln
  #sp.iterations = iters.int
  {.emit:"sp.iterations = CG.IterationsToComplete;".}
  #[
    let t0 = getTics()
    let t1 = getTics()
    echo "Grid time: ", (t1-t0).seconds
    #soln2 := gsrc
    soln2 := gsoln
  ]#
