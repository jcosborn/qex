import ospaths
import io/qio

const gridDir {.strdefine.} = getHomeDir() & "/lqcd/install/grid"
const gridPassC = "-I" & gridDir / "include"
const gridPassL = "-L" & gridDir & "/lib -lGrid -lz"
{.passC: gridPassC.}
{.passL: gridPassL.}
#{.passC: "-diag-disable=469".}
static:
  echo "Using Grid: ", gridDir
  echo "Grid compile flags: ", gridPassC
  echo "Grid link flags: ", gridPassL

{.pragma: gh, header:"Grid/Grid.h".}

type
  siginfo_t* {.importcpp,header:"signal.h".} = object
  stdvector*[T] {.importcpp:"std::vector",header:"<vector>".} = object
  GridVComplex* {.importcpp:"Grid::vComplex",gh.} = object
  Coordinate* {.importcpp:"Grid::Coordinate",gh.} = object
  GridBase* {.importcpp:"Grid::GridBase",gh,inheritable.} = object
  GridCartesian* {.importcpp:"Grid::GridCartesian",gh.} = object of GridBase
  GridRedBlackCartesian* {.importcpp:"Grid::GridRedBlackCartesian",gh.} = object of GridBase
  GridLatticeGaugeField* {.importcpp:"Grid::LatticeGaugeField",gh,byref.} = object
  GridVectorObject*[T] {.importcpp:"'0::vector_object",gh.} = object
  GridScalarObject*[T] {.importcpp:"'0::scalar_object",gh.} = object
  GridSU*[T:static int] {.importcpp:"Grid::SU",gh.} = object
  GridPeriodicGimplR* {.importcpp:"Grid::PeriodicGimplR",gh.} = object
  GridWilsonLoops*[T] {.importcpp:"Grid::WilsonLoops",gh.} = object
  GridFermion*[T] {.importcpp:"'0::FermionField",gh,byref.} = object
  GridNaiveStaggeredFermionR* {.
    importcpp:"Grid::NaiveStaggeredFermionR",gh.} = object
  GridImprovedStaggeredFermionR* {.
    importcpp:"Grid::ImprovedStaggeredFermionR",gh.} = object

type
  GridParity* = enum
    gpEven, gpOdd
#var
#  GridEven {.importc:"Grid::Even",gh.}
#  GridOdd {.importc:"Grid::Odd",gh.}

proc newStdVector*[T](n: int): stdvector[T] {.
  constructor,importcpp:"'0(#)".}

proc newStdVector*[T](a,b: ptr T): stdvector[T] {.
  constructor,importcpp:"'0(#,#)".}

proc newStdVector*[R,T](a,b: ptr T): stdvector[R] {.
  constructor,importcpp:"'0(#,#)".}

template `+`*[T](a: ptr T, n: int): untyped =
  cast[ptr T](cast[ByteAddress](a) + n*sizeof(T))

proc newStdVector*[T](a: openArray[T]): stdvector[T] =
  let a0 = unsafeaddr a[0]
  let a1 = a0 + a.len
  newStdVector(a0, a1)

proc newStdVector*[R,T](t: typedesc[R]; a: openArray[T]): stdvector[R] =
  let a0 = unsafeaddr a[0]
  let a1 = a0 + a.len
  newStdVector[R,T](a0, a1)

proc Nsimd*(t: typedesc[GridVComplex]): int {.importcpp:"Grid::vComplex::Nsimd()".}

proc newCoordinate*(v: stdvector[cint]): Coordinate {.
  constructor,importcpp:"Grid::Coordinate(#)".}

proc newCoordinate*(v: seq[cint]): Coordinate {.
  constructor,importcpp:"Grid::Coordinate(#)".}

proc newCoordinate*(s: seq[int]): Coordinate =
  let v = newStdVector(cint,s)
  newCoordinate(v)

proc GridDefaultSimd*(ndim,nsimd: int): Coordinate {.
  importc:"Grid::GridDefaultSimd",gh.}

proc newGridCartesian*(dimensions,simd_layout,processor_grid:Coordinate):
  GridCartesian {.importcpp,constructor.}

proc newGridRedBlackCartesian*(x: ptr GridBase):
  GridRedBlackCartesian {.importcpp,constructor.}
proc newGridRedBlackCartesian*(x: ptr GridCartesian):
  GridRedBlackCartesian {.importcpp,constructor.}
template newGridRedBlackCartesian*(x: GridCartesian): untyped =
  newGridRedBlackCartesian(unsafeaddr x)

proc newGridLatticeGaugeField*(x: ptr GridCartesian):
  GridLatticeGaugeField {.importcpp,constructor.}
template gauge*(x: GridCartesian): untyped =
  newGridLatticeGaugeField(unsafeaddr x)

proc lSites*(x: ptr GridBase): cint {.importcpp,gh.}

proc Grid*(x: GridLatticeGaugeField): ptr GridBase {.importcpp.}
proc Grid*(x: GridFermion): ptr GridBase {.importcpp.}

template vector_obj*(x: typedesc[GridLatticeGaugeField]): untyped =
  GridVectorObject[GridLatticeGaugeField]
template scalar_obj*(x: typedesc[GridLatticeGaugeField]): untyped =
  GridScalarObject[GridVectorObject[GridLatticeGaugeField]]

template vector_obj*[T:GridFermion](x: T): untyped =
  GridVectorObject[T]
template scalar_obj*[T:GridFermion](x: T): untyped =
  GridScalarObject[T]

proc vectorizeFromLexOrdArray*[T](x: stdVector[T], r: any) {.
  importc,gh.}
proc vectorizeFromLexOrdArray*[T](x: stdVector[T], r: GridLatticeGaugeField) {.
  importcpp:"vectorizeFromLexOrdArray(#,#)",gh.}

proc ColdConfiguration*[T](t: typedesc[T], g: GridLatticeGaugeField) {.
  importcpp:"'1::ColdConfiguration(# #)".}

proc avgPlaquette*[T](t: typedesc[T], g: GridLatticeGaugeField):
  float {.importcpp:"'1::avgPlaquette(# #)".}

proc newGridFermion*[T](x: ptr GridBase):
  GridFermion[T] {.importcpp,constructor.}
template fermion*[T](x: GridBase, y: typedesc[T]): untyped =
  newGridFermion[T](unsafeaddr x)
template fermion*[T](x: ptr GridBase, y: typedesc[T]): untyped =
  newGridFermion[T](x)

proc checkerboard*(x: ptr GridFermion): int =
  var r = 0
  #{.emit:"auto t = x->Checkerboard();".}
  #{.emit:"if(t==Grid::Odd){r=1;}".}
  {.emit:"r = x->Checkerboard();".}
  r

proc checkerboard*(x: ptr GridFermion, y: GridParity) =
  {.emit:"x->Checkerboard() = y;".}
template checkerboard*(x: var GridFermion, y: GridParity) =
  checkerboard(addr x, y)
template even*(x: var GridFermion) = x.checkerboard(gpEven)
template odd*(x: var GridFermion) = x.checkerboard(gpOdd)
