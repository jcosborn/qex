import qex
import physics/stagSolve
import grid/gridImpl

#[
template getGrid(g: Field): untyped =
  let lo = g.l
  let latt_size = newCoordinate(lo.physGeom)
  #let simd_layout = newCoordinate(lo.innerGeom)
  #let simd_layout = GridDefaultSimd(lo.nDim, GridVComplex.Nsimd);
  let simd_layout = GridDefaultSimd(lo.nDim, Nsimd(GridVComplex));
  let mpi_layout = newCoordinate(lo.rankGeom)
  let grid = newGridCartesian(latt_size,simd_layout,mpi_layout)
]#

proc gridSolveXX*(s:Staggered; r,t:Field; m:SomeNumber; sp: var SolverParams;
                  parEven = true) =
  let lo = r.l
  let latt_size = newCoordinate(lo.physGeom)
  #let simd_layout = newCoordinate(lo.innerGeom)
  #let simd_layout = GridDefaultSimd(lo.nDim, GridVComplex.Nsimd);
  let simd_layout = GridDefaultSimd(lo.nDim, Nsimd(GridVComplex));
  let mpi_layout = newCoordinate(lo.rankGeom)
  let grid = newGridCartesian(latt_size,simd_layout,mpi_layout)
  let rbgrid = newGridRedBlackCartesian(grid)
  var mass = m
  var res = sqrt sp.r2req
  var maxit = sp.maxits
  var gfl = grid.gauge()

  if s.g.len == 4: # plain staggered
    type ferm = GridNaiveStaggeredFermionR
    var gsrc = rbgrid.fermion(ferm)
    var gsoln = rbgrid.fermion(ferm)
    if parEven:
      gsrc.even
      gsoln.even
    else:
      gsrc.odd
      gsoln.odd
    gsrc := t
    {.emit:"using namespace Grid;".}
    {.emit:[gsoln," = Zero();"].}
    s.g.stagPhase([0,1,3,7])
    gfl := s.g[0..3]
    s.g.stagPhase([0,1,3,7])
    {.emit:["using Stag = ",GridNaiveStaggeredFermionR,";"].}
    {.emit:"using FermionField = Stag::FermionField;".}
    {.emit:["Stag Ds(grid,rbgrid,2.*",mass,",2.,1.);"].}
    {.emit:"Ds.ImportGauge(gfl);".}
    {.emit:"SchurStaggeredOperator<Stag,FermionField> HermOp(Ds);".}
    {.emit:["ConjugateGradient<FermionField> CG(",res,", ",maxit,", false);"].}
    {.emit:["CG(HermOp, ",gsrc,", ",gsoln,");"].}
    {.emit:[sp,".iterations = CG.IterationsToComplete;"].}
    var rr = r
    rr := gsoln
  elif s.g.len == 8: # Naik staggered
    type ferm = GridImprovedStaggeredFermionR
    var gsrc = rbgrid.fermion(ferm)
    var gsoln = rbgrid.fermion(ferm)
    if parEven:
      gsrc.even
      gsoln.even
    else:
      gsrc.odd
      gsoln.odd
    gsrc := t
    {.emit:"using namespace Grid;".}
    {.emit:[gsoln," = Zero();"].}
    var gll = grid.gauge()
    gfl := @[s.g[0],s.g[2],s.g[4],s.g[6]]
    gll := @[s.g[1],s.g[3],s.g[5],s.g[7]]
    {.emit:["using ImpStag = ",GridImprovedStaggeredFermionR,";"].}
    {.emit:"using FermionField = ImpStag::FermionField;".}
    {.emit:"ImpStag Ds(grid,rbgrid,2.*mass,2.,2.,1.);".}
    {.emit:"Ds.ImportGaugeSimple(gll,gfl);".}
    {.emit:"SchurStaggeredOperator<ImpStag,FermionField> HermOp(Ds);".}
    {.emit:"ConjugateGradient<FermionField> CG(res, maxit, false);".}
    {.emit:["CG(HermOp, ",gsrc,", ",gsoln,");"].}
    {.emit:[sp,".iterations = CG.IterationsToComplete;"].}
    var rr = r
    rr := gsoln
  else:
    qexError "unknown s.g.len: ", s.g.len

  #sp.iterations = iters.int
  #[
    let t0 = getTics()
    let t1 = getTics()
    echo "Grid time: ", (t1-t0).seconds
    #soln2 := gsrc
    soln2 := gsoln
  ]#

proc gridSolveEE*(s:Staggered; r,t:Field; m:SomeNumber; sp: var SolverParams) =
  gridSolveXX(s, r, t, m, sp, parEven=true)

proc gridSolveOO*(s:Staggered; r,t:Field; m:SomeNumber; sp: var SolverParams) =
  gridSolveXX(s, r, t, m, sp, parEven=false)
