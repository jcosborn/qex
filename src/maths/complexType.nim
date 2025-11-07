import macros
import base/[basicOps, wrapperTypes]
import complexProxy
export complexProxy

getOptimPragmas()

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

template complexObj*[TR,TI](x: typedesc[TR], y: typedesc[TI]): typedesc =
  ComplexObj[typeof(TR),typeof(TI)]
template newComplexObj*[TR,TI](x: TR, y: TI): untyped =
  ComplexObj[typeof(TR),typeof(TI)](reX: x, imX: y)

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
template asVarComplex*(x: typed): auto = newComplexProxy(x)

template isComplex*(x: ComplexProxy): auto = true
template isWrapper*(x: ComplexObj): bool = false
template isWrapper*(x: typedesc[ComplexObj]): bool = false

template has*[R,I](x: typedesc[ComplexObj[R,I]], y: typedesc): bool =
  mixin has
  has(R.type, y) or has(I.type, y)

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

template re*(x: ComplexObj): auto = x.reX
macro re*(x: ComplexObj{nkObjConstr}): auto = x[1][1]

template im*(x: ComplexObj): auto = x.imX
macro im*(x: ComplexObj{nkObjConstr}): auto = x[2][1]

template `re=`*(x: ComplexObj, y: typed) =
  mixin assign
  assign(x.reX, y)
template `im=`*(x: ComplexObj, y: typed) =
  mixin assign
  assign(x.imX, y)

overloadAsReal(SomeNumber)
template I*(x: SomeNumber): auto = newImag(x)

template numberType*[T](x: ComplexProxy[T]): typedesc = numberType(T)
template numberType*[T](x: typedesc[ComplexProxy[T]]): typedesc =
  mixin numberType
  numberType(T)
template numberType*[T](x: ComplexObj[T,T]): typedesc = numberType(T)
template numberType*[T](x: type ComplexObj[T,T]): typedesc = numberType(T)

template numNumbers*[TR,TI](x: ComplexObj[TR,TI]): auto =
  mixin numNumbers
  numNumbers(TR)+numNumbers(TI)
template numNumbers*[TR,TI](x: type ComplexObj[TR,TI]): auto =
  mixin numNumbers
  numNumbers(TR)+numNumbers(TI)
template numNumbers*[T](x: ComplexProxy[T]): auto =
  mixin numNumbers
  numNumbers(T)
template numNumbers*[T](x: type ComplexProxy[T]): auto =
  mixin numNumbers
  numNumbers(T)

template simdType*[T](x: ComplexObj[T,T]): typedesc =
  mixin simdType
  simdType(T)
template simdType*[T](x: type ComplexObj[T,T]): typedesc =
  mixin simdType
  simdType(T)
template simdType*[T](x: ComplexProxy[T]): typedesc = simdType(T)
template simdType*[T](x: type ComplexProxy[T]): typedesc = simdType(T)

template simdLength*[T](x: ComplexObj[T,T]): auto =
  mixin simdLength
  simdLength(T)
template simdLength*[T](x: type ComplexObj[T,T]): auto =
  mixin simdLength
  simdLength(T)
template simdLength*[T](x: ComplexProxy[T]): auto = simdLength(T)
template simdLength*[T](x: type ComplexProxy[T]): auto = simdLength(T)

#template simdSum*(x: ComplexObj): auto =
#  newComplexObj(simdSum(x.re),simdSum(x.im))
proc simdSum*(x: ComplexObj): auto {.alwaysInline.} =
  newComplexObj(simdSum(x.re),simdSum(x.im))
template simdSum*(x: ComplexProxy): auto = asComplex(simdSum(x[]))

template getNc*(x: ComplexProxy): auto = 1
template getNs*(x: ComplexProxy): auto = 1

template toSingle*[TR,TI](x: typedesc[ComplexObj[TR,TI]]): typedesc =
  mixin toSingle
  ComplexObj[toSingle(type(TR)),toSingle(type(TI))]
template toSingle*[T](x: typedesc[ComplexProxy[T]]): typedesc =
  mixin toSingle
  ComplexProxy[toSingle(type(T))]
template toSingleImpl*(xx: ComplexObj): auto =
  let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
  newComplexObj(toSingle(x.re),toSingle(x.im))
  #mixin toSingleX
  #toSingleX(toDerefPtr xx)

template toDouble*[TR,TI](x: typedesc[ComplexObj[TR,TI]]): typedesc =
  mixin toDouble
  ComplexObj[toDouble(type(TR)),toDouble(type(TI))]
template toDouble*[T](x: typedesc[ComplexProxy[T]]): typedesc =
  mixin toDouble
  ComplexProxy[toDouble(type(T))]
#template toDoubleImpl*(xx: ComplexObj): auto =
#  let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
#  newComplexObj(toDouble(x.re),toDouble(x.im))
#  #mixin toDoubleX
#  #toDoubleX(toDerefPtr xx)
proc toDoubleImpl*(x: ComplexObj): auto {.alwaysInline,noInit.} =
  newComplexObj(toDouble(x.re),toDouble(x.im))

template load1*(x: ComplexProxy): auto = x
template load1*(x: RealProxy): auto = x
template load1*(x: ImagProxy): auto = x

template eval*[TR,TI](x: typedesc[ComplexObj[TR,TI]]): typedesc =
  mixin eval
  complexObj(eval(typeof TR), eval(typeof TI))

template map*(xx: ComplexProxy; f: untyped): auto =
  let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
  newComplex(f(x.re),f(x.im))

template add*(r: ComplexProxy, x: SomeNumber, y: ComplexProxy3) =
  r := x + y
template add*(r: ComplexProxy, x: ComplexProxy2, y: SomeNumber) =
  r := x + y
template add*(r: ComplexProxy, x: RealProxy, y: ComplexProxy3) =
  r := x + y
template add*(r: ComplexProxy, x: ComplexProxy2, y: RealProxy) =
  r := x + y
proc add*[R,X,Y:ComplexProxy](r: var R, x: X, y: Y) {.alwaysInline.} =
#template add*[R,X,Y:ComplexProxy](rr: R, xx: X, yy: Y) =
#  let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
#  let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
#  let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
  mixin add
  add(r.re, x.re, y.re)
  add(r.im, x.im, y.im)

template sub*(r: ComplexProxy, x: SomeNumber, y: ComplexProxy3) =
  r := x - y
template sub*(r: ComplexProxy, x: RealProxy, y: ComplexProxy3) =
  r := x - y
template sub*(r: ComplexProxy, x: ComplexProxy2, y: SomeNumber) =
  r := x - y
template sub*(r: ComplexProxy, x: ComplexProxy2, y: ComplexProxy3) =
  r := x - y

template neg*(r: ComplexProxy, x: ComplexProxy2) =
  r := neg(x)


template mul*(r: SomeNumber, x: ImagProxy2, y: ImagProxy3) =
  r := x * y
template mul*(r: ImagProxy, x: ImagProxy2, y: SomeNumber) =
  r := x * y

template mul*[R,X:ComplexProxy;Y:SomeNumber](rr: R, xx: X, yy: Y) =
  mixin mul
  let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
  let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
  let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
  mul(r.re, x.re, y)
  mul(r.im, x.im, y)

proc mul*[R,Y:ComplexProxy;X:SomeNumber](r: var R, x: X, y: Y) {.alwaysInline.} =
#template mul*[R,Y:ComplexProxy;X:SomeNumber](rr: R, xx: X, yy: Y) =
#  let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
#  let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
#  let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
  mixin mul
  mul(r.re, x, y.re)
  mul(r.im, x, y.im)

template mul*[R,X:ComplexProxy;Y:RealProxy](rr: R, xx: X, yy: Y) =
  mixin mul
  let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
  let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
  let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
  mul(r.re, x.re, y[])
  mul(r.im, x.im, y[])

template mul*[R,Y:ComplexProxy;X:RealProxy](rr: R, xx: X, yy: Y) =
  mixin mul
  let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
  let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
  let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
  mul(r.re, x[], y.re)
  mul(r.im, x[], y.im)

template mul*[R,X:ComplexProxy;Y:ImagProxy](rr: R, xx: X, yy: Y) =
  mixin mul
  let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
  let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
  let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
  mul(r.re, x.im, -y[])
  mul(r.im, x.re,  y[])

template mul*[R,Y:ComplexProxy;X:ImagProxy](rr: R, xx: X, yy: Y) =
  mixin mul
  let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
  let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
  let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
  mul(r.re, -x[], y.im)
  mul(r.im,  x[], y.re)

proc mul*[R,X,Y:ComplexProxy](r: var R, x: X, y: Y) {.alwaysInline.} =
#template mul*[R,X,Y:ComplexProxy](r: R, xx: X, yy: Y) =
#  #let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
#  let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
#  let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
  #mixin mul, imadd
  #mul(  r, asReal(unsafeaddr x.re), y)
  #imadd(r, asImag(unsafeaddr x.im), y)
  #r.re = x.re*y.re - x.im*y.im
  #r.im = x.re*y.im + x.im*y.re
  mixin mul, imadd, imsub
  mul(r.re, x.re, y.re)
  mul(r.im, x.re, y.im)
  imsub(r.re, x.im, y.im)
  imadd(r.im, x.im, y.re)

template imadd*(r: SomeNumber, x: ImagProxy2, y: ImagProxy3) =  r -= x*y
template imadd*(r: ImagProxy, x: ImagProxy2, y: SomeNumber) =  r += x*y

proc imadd*[R,Y:ComplexProxy;X:RealProxy](r: var R, x: X, y: Y) {.alwaysInline.} =
#template imadd*[R,Y:ComplexProxy;X:RealProxy](r: R, xx: X, yy: Y) =
#  let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
#  let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
  mixin imsub, imadd
  imadd(r.re, x[], y.re)
  imadd(r.im, x[], y.im)

proc imadd*[R,Y:ComplexProxy;X:ImagProxy](r: var R, x: X, y: Y) {.alwaysInline.} =
#template imadd*[R,Y:ComplexProxy;X:ImagProxy](r: R, xx: X, yy: Y) =
#  let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
#  let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
  mixin imsub, imadd
  imsub(r.re, x[], y.im)
  imadd(r.im, x[], y.re)

proc imadd*[R,X,Y:ComplexProxy](r: var R, x: X, y: Y) {.alwaysInline.} =
#template imadd*[R,X,Y:ComplexProxy](r: R, xx: X, yy: Y) =
#  let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
#  let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
  #mixin imadd
  #imadd(r, asReal(unsafeaddr x.re), y)
  #imadd(r, asImag(unsafeaddr x.im), y)
  #r.re += x.re*y.re - x.im*y.im
  #r.im += x.re*y.im + x.im*y.re
  mixin imadd, imsub
  imadd(r.re, x.re, y.re)
  imadd(r.im, x.re, y.im)
  imsub(r.re, x.im, y.im)
  imadd(r.im, x.im, y.re)

#template imsub*(r: ComplexProxy, x: ComplexProxy2, y: ComplexProxy3) =  r -= x*y
proc imsub*[R,X,Y:ComplexProxy](r: var R, x: X, y: Y) {.alwaysInline.} =
  mixin imadd, imsub
  imsub(r.re, x.re, y.re)
  imsub(r.im, x.re, y.im)
  imadd(r.re, x.im, y.im)
  imsub(r.im, x.im, y.re)


template norm2*(r: auto, x: ComplexProxy2) =
  r = x.norm2
template inorm2*(r: auto, x: ComplexProxy2) =
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
  template pos(x: SomeNumber): auto = x
  template neg(x: SomeNumber): auto = -x
  template conj(x: SomeNumber): auto = x
  template adj(x: SomeNumber): auto = x
  template transpose(x: SomeNumber): auto = x
  template trace(x: SomeNumber): auto = x
  template norm2(x: SomeNumber): auto = x*x
  template inv[T: SomeNumber](x: T): auto = ((T)1)/x

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
