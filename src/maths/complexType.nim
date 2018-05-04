import macros
import base
import complexProxy
export complexProxy

type
  AsReal*[T] = RealProxy[T]
  AsImag*[T] = ImagProxy[T]
  ComplexObj*[TR,TI] = object
    reX*: TR
    imX*: TI
  ComplexObj2*[TR,TI] = ComplexObj[TR,TI]
  Complex*[TR,TI] = ComplexProxy[ComplexObj[TR,TI]]
  Complex2*[TR,TI] = Complex[TR,TI]
  Complex3*[TR,TI] = Complex[TR,TI]
  AsComplex*[T] = ComplexProxy[T]
  ComplexType*[T] = Complex[T,T]

template newRealImpl*(x: typed): untyped = x
template newImagImpl*(x: typed): untyped = newImagProxy(x)
template newComplexImplU*(x,y: typed): untyped =
  newComplexProxy(ComplexObj[type(x),type(y)](reX: x, imX: y))
template newComplexImpl*(x,y: typed): untyped =
  flattenCallArgs(newComplexImplU, x, y)
template newReal*(x: typed): untyped = newRealImpl(x)
template newImag*(x: typed): untyped = newImagImpl(x)
template newComplex*(x,y: typed): untyped = newComplexImpl(x,y)
template asReal*(x: untyped): untyped = newRealProxy(x)

template isWrapper*(x: ComplexObj): untyped = false
template isWrapper*(x: ComplexProxy): untyped = true
template asWrapper*(x: ComplexProxy, y: typed): untyped =
  newComplexProxy(y)
template asVarWrapper*(x: ComplexProxy, y: typed): untyped =
  asVar(newComplexProxy(y))

template re*(x: ComplexObj): untyped = x.reX
macro re*(x: ComplexObj{nkObjConstr}): untyped =
  #echo x.treerepr
  result = x[1][1]
  #echo result.treerepr
template im*(x: ComplexObj): untyped = x.imX
macro im*(x: ComplexObj{nkObjConstr}): untyped =
  #echo x.treerepr
  result = x[2][1]
  #echo result.treerepr

template `re=`*(x: ComplexObj, y: untyped): untyped =
  #dumpTree: x
  mixin `:=`
  x.reX := y
  #assign(x.reX, y)
#template `re=`*[TR,TI](x: ComplexObj[TR,TI], y: TR): untyped =
#  #x.reX = y
#  assign(x.reX, y)
template `im=`*(x: ComplexObj, y: untyped): untyped =
  x.imX := y
  #assign(x.imX, y)
#template `im=`*[TR,TI](x: ComplexObj[TR,TI], y: TI): untyped =
#  #x.imX = y
#  assign(x.imX, y)

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
template simdType*[T](x: ComplexProxy[T]): untyped = simdType(T)
template simdType*[TR,TI](x: ComplexObj[TR,TI]): untyped =
  mixin simdType
  simdType(TR)
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

template load1*(x: ComplexProxy): untyped = x
template load1*(x: RealProxy): untyped = x
template load1*(x: ImagProxy): untyped = x
template eval*(x: ComplexProxy): untyped = newComplexProxy(eval(x[]))
template eval*(x: ComplexObj): untyped =
  mixin eval
  let er = eval(x.re)
  let ei = eval(x.im)
  ComplexObj2[type(er),type(ei)](reX: er, imX: ei)
#template eval*(x: Complex): untyped =
#  newComplex(eval(x.re),eval(x.im))

template map*(x: ComplexObj; f: untyped): untyped =
  let fr = f(x.re)
  let fi = f(x.im)
  ComplexObj2[type(fr),type(fi)](reX: fr, imX: fi)

template toDoubleImpl*(xx: ComplexObj): untyped =
  let x = xx
  let tdiR = toDouble(x.re)
  let tdiI = toDouble(x.im)
  ComplexObj2[type(tdiR),type(tdiI)](reX: tdiR, imX: tdiI)


template add*(r: var ComplexProxy, x: ComplexProxy2, y: ComplexProxy3):
         untyped =  assign(r,x+y)
template add*(r: ComplexProxy, x: SomeNumber, y: ComplexProxy3): untyped =
  r := x + y

template sub*(r: ComplexProxy, x: SomeNumber, y: ComplexProxy3): untyped =
  r := x - y
template sub*(r: ComplexProxy, x: ComplexProxy2, y: SomeNumber): untyped =
  r := x - y
template sub*(r: ComplexProxy, x: ComplexProxy2, y: ComplexProxy3): untyped =
  r := x - y

template neg*(r: ComplexProxy, x: ComplexProxy2): untyped =
  r := neg(x)
template dot*(x: ComplexProxy, y: ComplexProxy2): untyped =
  trace( x.adj * y )
template mulCCR*(r: var ComplexProxy, y: ComplexProxy2, x: untyped):
         untyped =  assign(r,x*y)
template mul*(r: ComplexProxy, x: ComplexProxy2, y: SomeNumber): untyped =
  r := x * y
template mul*(r: ComplexProxy, x: SomeNumber, y: ComplexProxy3): untyped =
  r := x * y
template mul*(r: ComplexProxy, x: ImagProxy, y: ComplexProxy3): untyped =
  r := x * y
template mul*(r: var ComplexProxy, x: ComplexProxy2, y: ComplexProxy3):
         untyped =  assign(r,x*y)
template imsub*(r: var ComplexProxy, x: ComplexProxy2, y: ComplexProxy3):
         untyped =  r -= x*y

template norm2*(r: any, x: ComplexProxy2): untyped =
  r = x.norm2
template inorm2*(r: any, x: ComplexProxy2): untyped =
  r += x.norm2

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
