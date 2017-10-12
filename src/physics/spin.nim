import base/wrapperTypes
import maths/types
import maths

type
  Spin*[T] = object
    v*: T
  Spin2*[T] = Spin[T]
  Spin3*[T] = Spin[T]
  SpinMatrix[I,J:static[int],T] = Spin[MatrixArray[I,J,T]]

template asSpin*(xx: untyped): untyped =
  staticTraceBegin: asSpin
  let x_asSpin = xx
  staticTraceEnd: asSpin
  Spin[type(x_asSpin)](v: x_asSpin)

template isWrapper*(x: Spin): untyped = true
template asWrapper*(x: Spin, y: typed): untyped =
  asSpin(y)
template asVarWrapper*(x: Spin, y: typed): untyped =
  asVar(asSpin(y))
makeDeref(Spin, x.T)
template `[]`*(x: Spin, i: any): untyped = x[][i]
template `[]`*(x: Spin, i: any, j: any): untyped = x[][i,j]
template `[]`*(x: Spin, i: any, j: any, y: any): untyped =
  x[][i,j] = y
forwardFunc(Spin, len)
forwardFunc(Spin, nrows)
forwardFunc(Spin, ncols)
forwardFunc(Spin, numberType)
forwardFunc(Spin, nVectors)
forwardFunc(Spin, simdType)
forwardFunc(Spin, simdLength)
template row*(x: Spin, i: any): untyped =
  mixin row
  asSpin(row(x[],i))
template setRow*(r: Spin; x: Spin2; i: int): untyped =
  setRow(r[], x[], i)

template binDDRet(fn,wr,T1,T2) =
  template fn*(x: T1, y: T2): untyped =
    #echoUntyped: fn
    wr(fn(x[], y[]))

binDDRet(`+`, asSpin, Spin, Spin2)
binDDRet(`-`, asSpin, Spin, Spin2)
binDDRet(`*`, asSpin, Spin, Spin2)
binDDRet(`/`, asSpin, Spin, Spin2)

template numberType*[T](x: typedesc[Spin[T]]): untyped = numberType(type(T))
#template numNuumbers*[T](x: typedesc[Spin[T]]): untyped = numberType(T)
template numNumbers*(x: Spin): untyped = numNumbers(x[])
template load1*(x: Spin): untyped = asSpin(load1(x[]))
template assign*(r: var Spin, x: SomeNumber) =
  assign(r[], x)
template assign*(r: var Spin, x: AsComplex) =
  assign(r[], x)
template assign*(r: var Spin, x: Spin2) =
  assign(r[], x[])
template `:=`*(r: var Spin, x: SomeNumber) =
  `:=`(r[], x)
template `:=`*(r: var Spin, x: Spin2) =
  r[] := x[]
template `+=`*(r: var Spin, x: Spin2) =
  staticTraceBegin: peqSpinSpin
  r[] += x[]
  staticTraceEnd: peqSpinSpin
template `-=`*(r: var Spin, x: Spin2) =
  r[] -= x[]
template `*=`*(r: var Spin, x: SomeNumber) =
  `*=`(r[], x)
template iadd*(r: var Spin, x: AsComplex) =
  iadd(r[], x)
template iadd*(r: var Spin, x: Spin2) =
  iadd(r[], x[])
template isub*(r: var Spin, x: Spin2) =
  isub(r[], x[])
template imul*(r: var Spin, x: SomeNumber) =
  imul(r[], x)
template imadd*(r: var Spin, x: Spin2, y: Spin3) =
  imadd(r[], x[], y[])
template imsub*(r: var Spin, x: Spin2, y: Spin3) =
  imsub(r[], x[], y[])
template `*`*(x: Spin, y: SomeNumber): untyped =
  asSpin(x[] * y)
template `*`*(x: SomeNumber, y: Spin2): untyped =
  asSpin(x * y[])
template `*`*(x: AsComplex, y: Spin2): untyped =
  asSpin(x * y[])
template mul*(r: var Spin, x: Spin2, y: Spin3) =
  mul(r[], x[], y[])
template random*(x: var Spin) =
  gaussian(x[], r)
template gaussian*(x: var Spin, r: var untyped) =
  gaussian(x[], r)
template projectU*(r: var Spin, x: Spin2) =
  projectU(r[], x[])
template norm2*(x: Spin): untyped = norm2(x[])
template inorm2*(r: var any, x: Spin2) = inorm2(r, x[])
template dot*(x: Spin, y: Spin2): untyped =
  dot(x[], y[])
template idot*(r: var any, x: Spin2, y: Spin3) = idot(r, x[], y[])
template redot*(x: Spin, y: Spin2): untyped =
  redot(x[], y[])
template trace*(x: Spin): untyped = trace(x[])

#template spinMatrix*[I,J:static[int],T](a: untyped): untyped =
#  Spin[MatrixArray[I,J,T]](v: MatrixArray[I,J,T](v: MatrixArrayObj[I,J,T](mat: a)))
template spinMatrix*(I,J,T,a: untyped): untyped =
  Spin[MatrixArray[I,J,T]](v: MatrixArray[I,J,T](v: MatrixArrayObj[I,J,T](mat: a)))

const z0 = ComplexType[float](v: ComplexObj[float](re: 0.0, im: 0.0))
const z1 = ComplexType[float](v: ComplexObj[float](re: 1.0, im: 0.0))
const zi = ComplexType[float](v: ComplexObj[float](re: 0.0, im: 1.0))
var gamma0* = spinMatrix(4,4,ComplexType[float],[[z1,z0,z0,z0],[z0,z1,z0,z0],[z0,z0,z1,z0],[z0,z0,z0,z1]])
var gamma4* = spinMatrix(4,4,ComplexType[float],[[z0,z0,z1,z0],[z0,z0,z0,z1],[z1,z0,z0,z0],[z0,z1,z0,z0]])

const spprojmat1p* = spinMatrix(2,4,ComplexType[float],
                                [[z1,z0,z0,z0],[z0,z1,z0,z0]])

const spreconmat1p* = spinMatrix(4,2,ComplexType[float],
                                 [[z1,z0],[z0,z0],[z0,z1],[z0,z0]])


proc spproj1p*(r: var any, x: any) =
  ## r: HalfFermion
  ## x: DiracFermion
  #let nc = x[0].len
  #for i in 0..<nc:
  #  r[0][i] = x[0][i] + x[2][i]
  #  r[1][i] = x[1][i] + x[3][i]
  r := spprojmat1p * x

template spproj1p*(x: any): untyped = spprojmat1p * x
template spproj2p*(x: any): untyped = spprojmat1p * x
template spproj3p*(x: any): untyped = spprojmat1p * x
template spproj4p*(x: any): untyped = spprojmat1p * x
template spproj1m*(x: any): untyped = spprojmat1p * x
template spproj2m*(x: any): untyped = spprojmat1p * x
template spproj3m*(x: any): untyped = spprojmat1p * x
template spproj4m*(x: any): untyped = spprojmat1p * x

proc sprecon1p*(r: var any, x: any) =
  ## r: DiracFermion
  ## x: HalfFermion
  let nc = x[0].len
  for i in 0..<nc:
    r[0][i] = x[0][i]
    r[1][i] = x[1][i]
    r[2][i] = x[0][i]
    r[3][i] = x[1][i]

template sprecon1p*(x: any): untyped = spreconmat1p * x
template sprecon2p*(x: any): untyped = spreconmat1p * x
template sprecon3p*(x: any): untyped = spreconmat1p * x
template sprecon4p*(x: any): untyped = spreconmat1p * x
template sprecon1m*(x: any): untyped = spreconmat1p * x
template sprecon2m*(x: any): untyped = spreconmat1p * x
template sprecon3m*(x: any): untyped = spreconmat1p * x
template sprecon4m*(x: any): untyped = spreconmat1p * x

when isMainModule:
  echo gamma0[0,0]
  let g2 = gamma0 + gamma4
  #echo g2[0,0]
