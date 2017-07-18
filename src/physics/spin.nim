import base/wrapperTypes
import maths/types
import maths

type
  Spin*[T] = object
    v*: T
  Spin2*[T] = Spin[T]
  Spin3*[T] = Spin[T]

template asSpin*(xx: typed): untyped =
  #static: echo "asSpin typed"
  #dumpTree: xx
  let x_asSpin = xx
  Spin[type(x_asSpin)](v: x_asSpin)

template isWrapper*(x: Spin): untyped = true
template asWrapper*(x: Spin, y: typed): untyped =
  #static: echo "asWrapper Spin"
  #dumpTree: y
  asSpin(y)
template asVarWrapper*(x: Spin, y: typed): untyped =
  #static: echo "asVarWrapper Spin"
  #var cy = asSpin(y)
  #cy
  asVar(asSpin(y))
#template `[]`*(x: Spin): untyped = x.v
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
    wr(fn(x[], y[]))

binDDRet(`+`, asSpin, Spin, Spin2)
binDDRet(`-`, asSpin, Spin, Spin2)
binDDRet(`*`, asSpin, Spin, Spin2)
binDDRet(`/`, asSpin, Spin, Spin2)

template numberType*[T](x: typedesc[Spin[T]]): untyped = numberType(T)
template numNuumbers*[T](x: typedesc[Spin[T]]): untyped = numberType(T)
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
  r[] += x[]
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
