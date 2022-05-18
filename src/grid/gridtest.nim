import grid/gridImpl

when isMainModule:
  import qex
  import physics/stagSolve
  import grid/Grid

  proc pointSource(r: Field; c: openArray[int]; ic: int) =
    let (ptRank,ptIndex) = r.l.rankIndex(c)
    if myRank==ptRank:
      r{ptIndex}[ic] := 1

  template noSimd[V,T](x: typedesc[Field[V,T]]): untyped =
    mixin noSimd
    Field[V,noSimd(type T)]
  template noSimd[T](x: typedesc[Color[T]]): untyped =
    mixin noSimd
    Color[noSimd(type T)]
  template noSimd[T](x: typedesc[AsVector[T]]): untyped =
    mixin noSimd
    AsVector[noSimd(type T)]
  template noSimd[N,T](x: typedesc[VectorArrayObj[N,T]]): untyped =
    mixin noSimd
    VectorArrayObj[N,noSimd(type T)]
  template noSimd[T](x: typedesc[ComplexType[T]]): untyped =
    mixin noSimd
    ComplexType[noSimd(type T)]
  template noSimd[T](x: typedesc[Simd[T]]): untyped =
    mixin noSimd
    numberType(type T)
  proc comp(a: Field, b: Field) =
    var ax,bx: a[0].type.noSimd
    var crd: array[4,cint]
    var ne = 0
    for i in a.l.singleSites:
      ax := a{i}
      bx := b{i}
      let d = norm2(ax-bx)
      if d>1e-10 and ne<10:
        inc ne
        #echo ne
        a.l.coord(crd,i)
        echo i, " ", crd, " ", d
        echo "  ", ax
        echo "  ", bx

  proc checkPlaq(g: seq[Field], gg: GridLatticeGaugeField) =
    let qp = g.plaq.sum
    let gp = GridWilsonLoops[GridPeriodicGimplR].avgPlaquette(gg)
    echo "QEX plaq:  ", qp
    echo "Grid plaq: ", gp
    echo " rel diff: ", abs(qp-gp)/qp

  proc checkNaive(g: seq[Field], gg: GridLatticeGaugeField) =
    let lo = g[0].l
    g.setBC
    g.stagPhase
    var s = newStag(g)
    var sp = initSolverParams()
    sp.backend = sbQex
    sp.r2req = 1e-12
    sp.verbosity = 1
    var src = lo.ColorVector()
    var src2 = lo.ColorVector()
    var soln = lo.ColorVector()
    var soln2 = lo.ColorVector()
    threads:
      src := 0
      soln := 0
      soln2 := 0
    pointSource(src, [0,0,0,0], 0)
    #pointSource(src, [0,0,0,2], 0)
    #pointSource(src, [0,0,2,0], 0)
    #pointSource(src, [0,0,0,1], 0)
    #pointSource(src, [0,0,1,0], 0)
    #pointSource(src, [0,1,0,0], 0)
    #pointSource(src, [1,0,0,0], 0)
    var mass = floatParam("mass", 1.0)
    echo "mass: ", mass
    #s.solve(soln, src, mass, sp)
    s.solveEE(soln, src, mass, sp)
    threads:
      soln := 0
    s.solveEE(soln, src, mass, sp)
    #src2 := 0
    #s.D(src2,src,mass)
    #s.D(soln,src,mass)

    type ferm = GridNaiveStaggeredFermionR
    let grid = cast[ptr GridCartesian](gg.Grid())
    var gfl = grid[].gauge()
    g.stagPhase([0,1,3,7])
    gfl := g
    g.stagPhase([0,1,3,7])
    let rbgrid = newGridRedBlackCartesian(grid)
    var gsrc = rbgrid.fermion(ferm)
    var gsoln = rbgrid.fermion(ferm)
    gsrc.even
    gsoln.even
    #src2 := 0
    #s.D(src2,src,mass)
    #s.Ddag(soln,src2,mass)
    gsrc := src
    #gsoln := soln
    {.emit:"using namespace Grid;".}
    {.emit:"gsoln = Zero();".}
    {.emit:"using ImpStag = NaiveStaggeredFermionR;".}
    {.emit:"using FermionField = ImpStag::FermionField;".}
    #{.emit:"ImpStag Ds(gll,gfl,*grid,rbgrid,mass,1.,1.,1.);".}
    {.emit:"ImpStag Ds(*grid,rbgrid,2.*mass,2.,1.);".}
    {.emit:"Ds.ImportGauge(gfl);".}
    #{.emit:"MdagMLinearOperator<ImprovedStaggeredFermionR,FermionField> HermOp(Ds);".}
    {.emit:"SchurStaggeredOperator<ImpStag,FermionField> HermOp(Ds);".}
    #{.emit:"HermOp.Op(gsrc,gsoln);".}
    #{.emit:"Ds.M(gsrc,gsoln);".}
    #{.emit:"gsrc = Zero();".}
    #{.emit:"Ds.M(gsoln,gsrc);".}
    {.emit:"ConjugateGradient<FermionField> CG(1e-6, 400, false);".}
    let t0 = getTics()
    {.emit:"CG(HermOp, gsrc, gsoln);".}
    let t1 = getTics()
    echo "Grid time: ", (t1-t0).seconds
    #soln2 := gsrc
    soln2 := gsoln
    #soln2 *= 0.25
    #soln.odd := 0
    echo norm2(soln-soln2)
    echo norm2(soln)
    echo norm2(soln2)
    comp(soln, soln2)
    #echo soln[0]
    #echo soln2[0]
    #echo soln[[0,0,1,0]]
    #echo soln2[lo.nEvenOuter]

  proc checkHisq(g: seq[Field], gg: GridLatticeGaugeField) =
    let lo = g[0].l
    var coef: HisqCoefs
    g.setBC
    g.stagPhase
    coef.init()
    var fl = lo.newGauge()
    var ll = lo.newGauge()
    coef.smear(g, fl, ll)
    #for i in 0..3: ll[i] := 0
    #var s = newStag3(g,g)
    var s = newStag3(fl, ll)
    var sp = initSolverParams()
    sp.backend = sbQex
    sp.r2req = 1e-12
    sp.verbosity = 1
    var src = lo.ColorVector()
    var src2 = lo.ColorVector()
    var soln = lo.ColorVector()
    var soln2 = lo.ColorVector()
    threads:
      src := 0
      soln := 0
      soln2 := 0
    pointSource(src, [0,0,0,0], 0)
    #pointSource(src, [0,0,0,2], 0)
    #pointSource(src, [0,0,2,0], 0)
    #pointSource(src, [0,0,0,1], 0)
    #pointSource(src, [0,0,1,0], 0)
    #pointSource(src, [0,1,0,0], 0)
    #pointSource(src, [1,0,0,0], 0)
    var mass = floatParam("mass", 1.0)
    echo "mass: ", mass
    #s.solve(soln, src, mass, sp)
    s.solveEE(soln, src, mass, sp)
    threads:
      soln := 0
    s.solveEE(soln, src, mass, sp)
    #src2 := 0
    #s.D(src2,src,mass)
    #s.D(soln,src,mass)

    type ferm = GridImprovedStaggeredFermionR
    let grid = cast[ptr GridCartesian](gg.Grid())
    var gfl = grid[].gauge()
    var gll = grid[].gauge()
    #fl.stagPhase([0,1,3,7])
    gfl := fl
    gll := ll
    #fl.stagPhase([0,1,3,7])
    let rbgrid = newGridRedBlackCartesian(grid)
    var gsrc = rbgrid.fermion(ferm)
    var gsoln = rbgrid.fermion(ferm)
    gsrc.even
    gsoln.even
    #src2 := 0
    #s.D(src2,src,mass)
    #s.Ddag(soln,src2,mass)
    gsrc := src
    #gsoln := soln
    {.emit:"using namespace Grid;".}
    {.emit:"gsoln = Zero();".}
    {.emit:"using ImpStag = ImprovedStaggeredFermionR;".}
    {.emit:"using FermionField = ImpStag::FermionField;".}
    #{.emit:"ImpStag Ds(gll,gfl,*grid,rbgrid,mass,1.,1.,1.);".}
    {.emit:"ImpStag Ds(*grid,rbgrid,2.*mass,2.,2.,1.);".}
    {.emit:"Ds.ImportGaugeSimple(gll,gfl);".}
    #{.emit:"MdagMLinearOperator<ImprovedStaggeredFermionR,FermionField> HermOp(Ds);".}
    {.emit:"SchurStaggeredOperator<ImpStag,FermionField> HermOp(Ds);".}
    #{.emit:"HermOp.Op(gsrc,gsoln);".}
    #{.emit:"Ds.M(gsrc,gsoln);".}
    #{.emit:"gsrc = Zero();".}
    #{.emit:"Ds.M(gsoln,gsrc);".}
    {.emit:"ConjugateGradient<FermionField> CG(1e-6, 400, false);".}
    let t0 = getTics()
    {.emit:"CG(HermOp, gsrc, gsoln);".}
    let t1 = getTics()
    echo "Grid time: ", (t1-t0).seconds
    #soln2 := gsrc
    soln2 := gsoln
    #soln2 *= 0.25
    #soln.odd := 0
    echo norm2(soln-soln2)
    echo norm2(soln)
    echo norm2(soln2)
    comp(soln, soln2)
    #echo soln[0]
    #echo soln2[0]
    #echo soln[[0,0,1,0]]
    #echo soln2[lo.nEvenOuter]

  proc testHisq(g: seq[Field]) =
    let lo = g[0].l
    var coef: HisqCoefs
    g.setBC
    g.stagPhase
    coef.init()
    var fl = lo.newGauge()
    var ll = lo.newGauge()
    coef.smear(g, fl, ll)
    for i in 0..3: ll[i] := 0
    #var s = newStag3(g,g)
    var s = newStag3(fl, ll)
    var sp = initSolverParams()
    sp.backend = sbQex
    sp.r2req = 1e-12
    sp.verbosity = 1
    var src = lo.ColorVector()
    var src2 = lo.ColorVector()
    var soln = lo.ColorVector()
    var soln2 = lo.ColorVector()
    threads:
      src := 0
      soln := 0
      soln2 := 0
    pointSource(src, [0,0,0,0], 0)
    #pointSource(src, [0,0,0,2], 0)
    #pointSource(src, [0,0,2,0], 0)
    #pointSource(src, [0,0,0,1], 0)
    #pointSource(src, [0,0,1,0], 0)
    #pointSource(src, [0,1,0,0], 0)
    #pointSource(src, [1,0,0,0], 0)
    var mass = floatParam("mass", 1.0)
    echo "mass: ", mass
    #s.solve(soln, src, mass, sp)
    s.solveEE(soln, src, mass, sp)
    threads:
      soln := 0
    s.solveEE(soln, src, mass, sp)
    #src2 := 0
    #s.D(src2,src,mass)
    #s.D(soln,src,mass)

  proc test() =
    defaultSetup()
    g.random
    #testHisq(g)

    let latt_size = newCoordinate(lat)
    #let simd_layout = newCoordinate(lo.innerGeom)
    let simd_layout = GridDefaultSimd(lat.len, GridVComplex.Nsimd);
    let mpi_layout = newCoordinate(lo.rankGeom)
    let grid = newGridCartesian(latt_size,simd_layout,mpi_layout)
    #let rbgrid = newGridRedBlackCartesian(grid)

    var gg = grid.gauge()
    #GridSU[3].ColdConfiguration(gg)
    gg := g

    checkPlaq(g,gg)
    checkNaive(g,gg)
    #checkHisq(g,gg)

  #[
  {.emit:"/*INCLUDESECTION*/\n#include <Grid/Grid.h>".}
  proc testgrid =
    #let latt_size = newCoordinate(lat)
    {.emit:"""
    using namespace Grid;
    typedef typename ImprovedStaggeredFermionR::FermionField FermionField;
    typename ImprovedStaggeredFermionR::ImplParams params;

    Coordinate latt_size   = GridDefaultLatt();
    Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();
    GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
    GridRedBlackCartesian     RBGrid(&Grid);

    //std::vector<int> seeds({1,2,3,4});
    //GridParallelRNG          pRNG(&Grid);  pRNG.SeedFixedIntegers(seeds);

    LatticeGaugeField Umu(&Grid);
    //SU<Nc>::HotConfiguration(pRNG,Umu);
    SU<Nc>::ColdConfiguration(Umu);
    auto gp = WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu);
    std::cout<<gp<<std::endl;

    /*
    FermionField    src(&Grid); random(pRNG,src);
    FermionField result(&Grid); result=Zero();
    FermionField  resid(&Grid);

    RealD mass=0.1;
    RealD c1=9.0/8.0;
    RealD c2=-1.0/24.0;
    RealD u0=1.0;
    ImprovedStaggeredFermionR Ds(Umu,Umu,Grid,RBGrid,mass,c1,c2,u0);

    ConjugateGradient<FermionField> CG(1.0e-8,10000);
    SchurRedBlackStaggeredSolve<FermionField> SchurSolver(CG);

    double volume=1.0;
    for(int mu=0;mu<Nd;mu++){
      volume=volume*latt_size[mu];
    }
    double t1=usecond();
    SchurSolver(Ds,src,result);
    double t2=usecond();

    // Schur solver: uses DeoDoe => volume * 1146
    double ncall=CG.IterationsToComplete;
    double flops=(16*(3*(6+8+8)) + 15*3*2)*volume*ncall; // == 66*16 +  == 1146

    std::cout<<GridLogMessage << "usec    =   "<< (t2-t1)<<std::endl;
    std::cout<<GridLogMessage << "flop/s  =   "<< flops<<std::endl;
    std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t2-t1)<<std::endl;
    */
    """.}
  ]#

  #Grid_init()
  qexInit()
  #qexSetFinalizeComms(false)
  test()
  #testgrid()
  #echoTimers()
  qexFinalize()
  #Grid_finalize()
