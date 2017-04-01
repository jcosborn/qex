import macros
import complexProxy
export complexProxy

type
  Imag*[T] = ImagProxy[T]
  ComplexObj*[TR,TI] = object
    reX*: TR
    imX*: TI
  Complex*[TR,TI] = ComplexProxy[ComplexObj[TR,TI]]
  Complex2*[TR,TI] = Complex[TR,TI]
  Complex3*[TR,TI] = Complex[TR,TI]
  AsComplex*[T] = ComplexProxy[T]

template newRealImpl*(x: typed): untyped = x
template newImagImpl*(x: typed): untyped = newImagProxy(x)
template newComplexImpl*(x,y: typed): untyped =
  newComplexProxy(ComplexObj[type(x),type(y)](reX: x, imX: y))
template newReal*(x: typed): untyped = newRealImpl(x)
template newImag*(x: typed): untyped = newImagImpl(x)
template newComplex*(x,y: typed): untyped = newComplexImpl(x,y)

template re*(x: ComplexObj): untyped = x.reX
macro re*(x: ComplexObj{nkObjConstr}): auto =
  #echo x.treerepr
  result = x[1][1]
  #echo result.treerepr
template im*(x: ComplexObj): untyped = x.imX
macro im*(x: ComplexObj{nkObjConstr}): auto =
  #echo x.treerepr
  result = x[2][1]
  #echo result.treerepr
template `re=`*(x: ComplexObj, y: untyped): untyped =
  x.reX = y
template `im=`*(x: ComplexObj, y: untyped): untyped =
  x.imX = y

overloadAsReal(SomeNumber)
template I*(x: SomeNumber): untyped = newImag(x)

template numberType*[T](x: ComplexProxy[T]): untyped = numberType(T)
template numberType*[T](x: typedesc[ComplexProxy[T]]): untyped =
  mixin numberType
  numberType(T)
template numberType*[T](x: ComplexObj[T,T]): untyped = numberType(T)
#template nVectors*[T](x: Complex[T,T]): untyped =
#  mixin nVectors
#  nVectors(T)
template numNumbers*(x: ComplexProxy): untyped =
  mixin numNumbers
  2*numNumbers(x.re)
template simdSum*(x: ComplexObj): untyped =
  newComplex(simdSum(x.re),simdSum(x.im))

template isComplex*(x: ComplexProxy): untyped = true
template asComplex*(x: untyped): untyped = newComplexProxy(x)
template asVarComplex*(x: untyped): untyped = newComplexProxy(x)
template imaddCRC*(r: untyped, x: untyped, y: untyped) =
  r.re += x * y.re
  r.im += x * y.im
template imaddCIC*(r: untyped, x: untyped, y: untyped) =
  r.re -= x * y.im
  r.im += x * y.re
template imaddCCR*(r: untyped, x: untyped, y: untyped) =
  r.re += x.re * y
  r.im += x.im * y
template imaddCCI*(r: untyped, x: untyped, y: untyped) =
  r.re -= x.im * y
  r.im += x.re * y
template imadd*(r: ComplexProxy, x: ComplexProxy2, y: ComplexProxy3) =
  r += x*y

when isMainModule:
  template pos(x: SomeNumber): untyped = x
  template neg(x: SomeNumber): untyped = -x
  template conj(x: SomeNumber): untyped = x
  template adj(x: SomeNumber): untyped = x
  template transpose(x: SomeNumber): untyped = x
  template trace(x: SomeNumber): untyped = x
  template norm2(x: SomeNumber): untyped = x*x
  template inv(x: SomeNumber): untyped = ((type(x))(1))/x

  proc testadd(a,b: float) =
    var z0 = newComplex(a,b)
    var z1 = newComplex(a,b)
    var z2 = newComplex(a,b)
    var z3 = newComplex(a,b)
    var z4 = newComplex(a,b)
    var z6 = z0+z1+z2+z3+z4
    z6 := z1+z2+z3
    z6 += z2+z3+z4
    echo z6
  testadd(1.0, 2.0)

  proc testmul(a,b: float) =
    var z0 = newComplex(a,b)
    var z1 = newComplex(a+1,b+1)
    var z2 = newComplex(a+2,b+2)
    var z3 = newComplex(a,b)
    var z4 = newComplex(a,b)
    var z5 = newComplex(a,b)
    z5 += z0*z1*z2*z3*z4
    echo z5
    #var z6 = z0*z1*z2*z3*z4
    #echo z6
  testmul(1.0, 2.0)

  proc testdiv(a,b: float) =
    var z0 = newComplex(a,b)
    var z1 = newComplex(a+1,b+1)
    var z2 = newComplex(a+2,b+2)
    var z3 = newComplex(a,b)
    var z4 = newComplex(a,b)
    var z5 = newComplex(a,b)
    z5 += z0*z1*z2*z3*z4/1.0
    z5 += z0*z1*z2*z3*z4/2.0.I
    z5 += z0*z1*z2*z3*z4/z4
    z5 *= 2.0.I + 3.0
    echo z5
    #var z6 = z0*z1*z2*z3*z4
    #echo z6
  testdiv(1.0, 2.0)

  var a = 1.0 + 2.0.I
  echo a

  echo a.pos
  echo +a

  echo a.neg
  echo -a

  echo a.conj
  echo ~a

  echo a.adj
  echo @a

  echo a.transpose
  echo <>a

  echo a.trace
  echo %a

  echo a.norm2
  echo |a

  echo a.inv
  echo /a

  echo divd(1.0,a)
  echo 1.0/a
