import grid/gridImpl

when isMainModule:
  import qex
  import grid/Grid
  import os
  {.emit:"/*INCLUDESECTION*/\n#include <Grid/Grid.h>".}
  {.emit:"/*INCLUDESECTION*/\nusing namespace Grid;".}
  discard paramCount()
  var argc {.importc: "cmdCount", global.}: cint
  var argv {.importc: "cmdLine", global.}: cstringArray

  proc unitPlaq =
    {.emit:"Grid_init(&cmdCount,&cmdLine);".}
    {.emit:"""
  Coordinate latt_size   = GridDefaultLatt();
  //Coordinate latt_size(16,16,16,16);
  Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();
  GridCartesian Grid(latt_size,simd_layout,mpi_layout);
  LatticeGaugeField Umu(&Grid);
  SU<Nc>::ColdConfiguration(Umu);
  auto gp = WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu);
  std::cout<<gp<<std::endl;
""".}
    {.emit:"Grid_finalize();".}

  #qexInit()
  #var argc {.importc: "cmdCount", global.}: cint
  #var argv {.importc: "cmdLine", global.}: cstringArray
  #{.emit:"Grid_init(&argc,&argv);
  unitPlaq()
  #qexFinalize()
