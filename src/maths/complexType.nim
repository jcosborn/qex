import macros
import base/basicOps
import complexProxy
export complexProxy

type
  AsReal*[T] = RealProxy[T]
  AsImag*[T] = ImagProxy[T]
  ComplexObj*[TR,TI] = object
    reX*: TR
    imX*: TI
  ComplexObj2*[TR,TI] = ComplexObj[TR,TI]
  ComplexTT*[TR,TI] = ComplexProxy[ComplexObj[TR,TI]]
  ComplexTT2*[TR,TI] = ComplexTT[TR,TI]
  ComplexTT3*[TR,TI] = ComplexTT[TR,TI]
  AsComplex*[T] = ComplexProxy[T]
  AsComplex2*[T] = ComplexProxy[T]
  ComplexT*[T] = ComplexTT[T,T]
  ComplexType*[T] = ComplexT[T]

template complexObj*[TR,TI](x: TR, y: TI): untyped =
  ComplexObj[typeof(TR),typeof(TI)](reX: x, imX: y)
template complexObj*[TR,TI](x: typedesc[TR], y: typedesc[TI]): typedesc =
  ComplexObj[typeof(TR),typeof(TI)]
template newComplexObj*[TR,TI](x: TR, y: TI): untyped = complexObj(x, y)

template newRealImpl*(x: typed): untyped = x
template newImagImpl*(x: typed): untyped = newImagProxy(x)
#template newComplexImplU*(x,y: typed): untyped =
#  newComplexProxy(ComplexObj[type(x),type(y)](reX: x, imX: y))
#template newComplexImplU*[T,U](x: T, y: U): untyped =
#  newComplexProxy(ComplexObj[type(T),type(U)](reX: x, imX: y))
template newComplexImpl*[TR,TI](x: TR, y: TI): untyped =
  #flattenCallArgs(newComplexImplU, x, y)
  newComplexProxy(ComplexObj[type(TR),type(TI)](reX: x, imX: y))
template newReal*(x: typed): untyped = newRealImpl(x)
template newImag*(x: typed): untyped = newImagImpl(x)
template newComplex*(x,y: typed): untyped = newComplexImpl(x,y)
template asReal*(x: typed): untyped = newRealProxy(x)
template asImag*(x: typed): untyped = newImagProxy(x)
template asComplex*(x: typed): auto = newComplexProxy(x)
template asComplex*(x: typedesc): typedesc = newComplexProxy(x)

template isWrapper*(x: ComplexObj): untyped = false

template `[]`*[T](x: AsComplex, i: T): untyped =
  when T is AsComplex:
    x[][i[]]
  elif T.isWrapper:
    indexed(x, i)
    #asVar(asComplex(x[][i]))
  else:
    x[][i]

template index*[TR,TI,I](x: typedesc[ComplexObj[TR,TI]], i: typedesc[I]): typedesc =
  when I.isWrapper:
    ComplexObj[index(TR.type,I.type),index(TI.type,I.type)]
  else:
    {.error.} #index(X[], I)

template index*[X:AsComplex,I](x: typedesc[X], i: typedesc[I]): typedesc =
  when I is AsComplex:
    index(X[], I[])
  elif I.isWrapper:
    asComplex(index(X.type[], I.type))
  else:
    index(X[], I)

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

template `re=`*(x: var ComplexObj, y: typed): untyped =
  #static: echo "co re="
  #debugType: x
  #debugType: y
  mixin assign
  #debugCall:
  assign(x.reX, y)
template `im=`*(x: ComplexObj, y: typed): untyped =
  #x.imX := y
  assign(x.imX, y)

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
template simdType*[T](x: type ComplexProxy[T]): untyped = simdType(T)
template simdType*[TR,TI](x: ComplexObj[TR,TI]): untyped =
  mixin simdType
  simdType(TR)
template simdLength*[TR,TI](x: ComplexObj[TR,TI]): untyped =
  mixin simdLength
  simdLength(TR)
template simdLength*[T](x: ComplexProxy[T]): untyped = simdLength(T)
template simdLength*[T](x: type ComplexProxy[T]): untyped = simdLength(T)
template simdSum*(x: ComplexObj): untyped =
  newComplexObj(simdSum(x.re),simdSum(x.im))
template getNc*(x: ComplexProxy): untyped = 1
template getNs*(x: ComplexProxy): untyped = 1

template toSingle*[TR,TI](x: typedesc[ComplexObj[TR,TI]]): untyped =
  ComplexObj[toSingle(type(TR)),toSingle(type(TI))]
template toSingle*[T](x: typedesc[ComplexProxy[T]]): untyped =
  ComplexProxy[toSingle(type(T))]

template load1*(x: ComplexProxy): untyped = x
template load1*(x: RealProxy): untyped = x
template load1*(x: ImagProxy): untyped = x
#template eval*(x: ComplexProxy): untyped = newComplexProxy(eval(x[]))
#template eval*(x: ComplexObj): untyped =
#  mixin eval
#  let er = eval(x.re)
#  let ei = eval(x.im)
#  ComplexObj2[type(er),type(ei)](reX: er, imX: ei)
template eval*(xx: ComplexProxy): untyped =
  let x = xx[]
  newComplex(eval(x.re),eval(x.im))
template eval*[TR,TI](x: typedesc[ComplexObj[TR,TI]]): typedesc =
  mixin eval
  complexObj(eval(typeof TR), eval(typeof TI))

template map*(xx: ComplexProxy; f: untyped): untyped =
  #let fr = f(x.re)
  #let fi = f(x.im)
  #ComplexObj2[type(fr),type(fi)](reX: fr, imX: fi)
  let x = xx[]
  newComplex(f(x.re),f(x.im))

#template toDoubleImpl*(xx: ComplexProxy): untyped =
  #let x = xx
  #let tdiR = toDouble(x.re)
  #let tdiI = toDouble(x.im)
  #ComplexObj2[type(tdiR),type(tdiI)](reX: tdiR, imX: tdiI)
  #let x = xx[]
  #newComplex(toDouble(x.re),toDouble(x.im))
  #mixin toDoubleX
  #toDoubleX(toDerefPtr xx)

template toSingleImpl*(xx: ComplexObj): untyped =
  let x = xx
  newComplexObj(toSingle(x.re),toSingle(x.im))
  #mixin toSingleX
  #toSingleX(toDerefPtr xx)

template toDoubleImpl*(xx: ComplexObj): untyped =
  let x = xx
  newComplexObj(toDouble(x.re),toDouble(x.im))
  #mixin toDoubleX
  #toDoubleX(toDerefPtr xx)

template add*(r: ComplexProxy, x: ComplexProxy2, y: ComplexProxy3):
         untyped =  assign(r,x+y)
template add*(r: ComplexProxy, x: SomeNumber, y: ComplexProxy3): untyped =
  r := x + y
template add*(r: ComplexProxy, x: ComplexProxy2, y: SomeNumber): untyped =
  r := x + y
template add*(r: ComplexProxy, x: RealProxy, y: ComplexProxy3): untyped =
  r := x + y
template add*(r: ComplexProxy, x: ComplexProxy2, y: RealProxy): untyped =
  r := x + y

template sub*(r: ComplexProxy, x: SomeNumber, y: ComplexProxy3): untyped =
  r := x - y
template sub*(r: ComplexProxy, x: RealProxy, y: ComplexProxy3): untyped =
  r := x - y
template sub*(r: ComplexProxy, x: ComplexProxy2, y: SomeNumber): untyped =
  r := x - y
template sub*(r: ComplexProxy, x: ComplexProxy2, y: ComplexProxy3): untyped =
  r := x - y

template neg*(r: ComplexProxy, x: ComplexProxy2): untyped =
  r := neg(x)
#template dot*(x: ComplexProxy, y: ComplexProxy2): untyped =
#  trace( x.adj * y )
template mulCCR*(r: ComplexProxy, y: ComplexProxy2, x: untyped):
         untyped =  assign(r,x*y)
template mul*(r: ComplexProxy, x: ComplexProxy2, y: SomeNumber): untyped =
  r := x * y
template mul*(r: ComplexProxy, x: SomeNumber, y: ComplexProxy3): untyped =
  r := x * y
template mul*(r: ComplexProxy, x: RealProxy, y: ComplexProxy3): untyped =
  r := x * y
template mul*(r: ComplexProxy, x: ImagProxy, y: ComplexProxy3): untyped =
  r := x * y
template mul*(r: ComplexProxy, x: ComplexProxy2, y: ComplexProxy3):
         untyped =  assign(r,x*y)
template mul*(r: SomeNumber, x: ImagProxy2, y: ImagProxy3): untyped =
  r := x * y
template mul*(r: ImagProxy, x: ImagProxy2, y: SomeNumber): untyped =
  r := x * y

template imadd*(r: SomeNumber, x: ImagProxy2, y: ImagProxy3):
         untyped =  r -= x*y
template imadd*(r: ImagProxy, x: ImagProxy2, y: SomeNumber):
         untyped =  r -= x*y
template imadd*(r: ComplexProxy, x: ComplexProxy2, y: ComplexProxy3) =
  r += x*y
template imaddCRC*(r: typed, x: typed, y: typed) =
  r.re += x * y.re
  r.im += x * y.im
template imaddCIC*(r: typed, x: typed, y: typed) =
  r.re -= x * y.im
  r.im += x * y.re
template imaddCCR*(r: typed, x: typed, y: typed) =
  r.re += x.re * y
  r.im += x.im * y
template imaddCCI*(r: typed, x: typed, y: typed) =
  r.re -= x.im * y
  r.im += x.re * y

template imsub*(r: ComplexProxy, x: ComplexProxy2, y: ComplexProxy3):
         untyped =  r -= x*y

template norm2*(r: auto, x: ComplexProxy2): untyped =
  r = x.norm2
template inorm2*(r: auto, x: ComplexProxy2): untyped =
  r += x.norm2


import simd

overloadAsReal(Simd)
template add*(r: AsComplex, x: Simd, y: AsComplex2) =
  r := add(x,y)
template sub*(r: AsComplex, x: Simd, y: AsComplex2) =
  r := sub(x,y)
template mul*(r: AsComplex, x: Simd, y: AsComplex2) =
  r := mul(x,y)


when isMainModule:
  template pos(x: SomeNumber): untyped = x
  template neg(x: SomeNumber): untyped = -x
  template conj(x: SomeNumber): untyped = x
  template adj(x: SomeNumber): untyped = x
  template transpose(x: SomeNumber): untyped = x
  template trace(x: SomeNumber): untyped = x
  template norm2(x: SomeNumber): untyped = x*x
  template inv[T: SomeNumber](x: T): untyped = ((T)1)/x

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
