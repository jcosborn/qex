import base/basicOps
import ../base/wrapperTypes
import ../maths/types
#import ../maths/complexNumbers
#export complexNumbers

#makeWrapperType(Simd)
makeWrapperF({wfDeref}, Simd)
type
  Simd2*[T] = Simd[T]
  Simd3*[T] = Simd[T]

#template asSimd*(x: Indexed): untyped = x[]
template asSimd*(x: Indexed): untyped =
  #asReal(x)
  asScalar(x)

template `[]`*[T](x: Simd[T]): untyped =
  when T is (RefWrap|ptr):
    derefXX(x)[]
  else:
    derefXX(x)
template `[]`*[T](x: typedesc[Simd[T]]): typedesc =
  #static: echo "x: ", x, " []: ", T.type
  when T is (RefWrap|ptr):
    T.type[]
  else:
    T.type
#[
template `[]`*[T](x: typedesc[Simd[T]]): untyped =
    T.type
template `[]`*[T](x: typedesc[Simd[RefWrap[T]]]): untyped =
    T.type[]
]#

# assume these don't have any lazy evaluations   # FIXME could be toDouble, etc.
template eval*[T](x: typedesc[Simd[T]]): typedesc =
  asSimd(typeof(T))

template `[]`*(x: Simd, i: typed): untyped =
  when i is Simd:
    when i[] is AsView:
      #indexed(x, i[][])
      var t = indexed(x, i[][])
      t
    else:
      #x[][i[]]
      var t = indexed(x, i[])
      t
  else:
    when i is AsView:
      indexed(x, i[])
    else:
      #x[][i]
      indexed(x, i)

template index*[T,I](x: typedesc[Simd[T]], i: typedesc[I]): typedesc =
  numberType(T)

#template doIndexed[T](x: T): untyped =
#  when T is Indexed:
#    x[]
#  else:
#    x
#template doIndexed[T](x: T): untyped = x
#template doIndexed(x: Indexed): untyped = doIndexed x[]
#template doIndexed(x: Scalar): untyped = doIndexed x[]
template doIndexed[T](x: T): untyped = eval(x)

template toPrec*(x: Simd, y: typedesc[float32]): untyped = toSingle(x)
template toPrec*(x: Simd, y: typedesc[float64]): untyped = toDouble(x)

template stripSimdAsView*[T](x: T): untyped = x
template stripSimdAsView*(x: AsView): untyped = stripSimdAsView x[]
template stripSimdAsView*(x: Simd): untyped = stripSimdAsView x[]
template `[]=`*(x: Simd, i: typed, y: typed): untyped =
  when y is Simd:
    x[][stripSimdAsView i] = doIndexed(y[])
  else:
    x[][stripSimdAsView i] = doIndexed(y)

template attrib(att: untyped): untyped {.dirty.} =
  template att*[T](x: typedesc[Simd[T]]): untyped =
    mixin att
    att(T)
  template att*[T](x: Simd[T]): untyped =
    mixin att
    att(T)

attrib(numberType)
attrib(numNumbers)
attrib(simdType)
attrib(simdLength)

template noSimd*[T](x: typedesc[Simd[T]]): untyped =
  numberType(type T)

# no return value
template p2(f: untyped) {.dirty.} =
  #template f*(x: var Simd, y: Simd2) =
  template f*[T1,T2](x: var Simd[T1], y: Simd[T2]) =
    #static: echo "f Simd Simd"
    #static: echo "  ", x.type
    #static: echo "  ", y.type
    mixin f
    #f(x[], y[])
    when numberType(T1) is numberType(T2):
      f(x[], doIndexed(y[]))
    else:
      f(x[], y.toPrec(numberType(T1))[])

template p2s(f: untyped) {.dirty.} =
  template f*(x: var Simd, y: SomeNumber) =
    #static: echo "f Simd Somenumber"
    #static: echo "  ", x.type
    #static: echo "  ", y.type
    mixin f
    f(x[], y)
  p2(f)

template p3(f: untyped) {.dirty.} =
  template f*[T1,T2,T3](x: var Simd[T1], y: Simd[T2], z: Simd[T3]) =
    mixin f
    when numberType(T1) is numberType(T2):
      when numberType(T1) is numberType(T3):
        f(x[], doIndexed(y[]), doIndexed(z[]))
      else:
        f(x[], doIndexed(y[]), z.toPrec(numberType(T1))[])
    else:
      when numberType(T1) is numberType(T3):
        f(x[], y.toPrec(numberType(T1))[], doIndexed(z[]))
      else:
        f(x[], y.toPrec(numberType(T1))[], z.toPrec(numberType(T1))[])

template p3s(f: untyped) {.dirty.} =
  template f*(x: var Simd, y: SomeNumber, z: Simd3) =
    mixin f
    f(x[], y, z[])
  template f*(x: var Simd, y: Simd2, z: SomeNumber) =
    mixin f
    f(x[], y[], z)
  p3(f)

p2(rsqrt)
p2(norm2)
p2s(assign)
p2s(`:=`)
p2s(`+=`)
p2s(`-=`)
p2s(`*=`)
p2s(`/=`)
p2s(iadd)
p2s(isub)
p2s(inorm2)
p3(imadd)
p3(imsub)
p3s(add)
p3s(sub)
p3s(mul)
p3s(divd)


# with return value
template f1(f: untyped): untyped {.dirty.} =
  template f*(x: Simd): untyped =
    mixin f
    asSimd(f(x[]))

template f2(f: untyped): untyped {.dirty.} =
  template f*[T1,T2](x: Simd[T1], y: Simd[T2]): untyped =
    mixin f
    #static: echo numberType(T1), " ", numberType(T2)
    when numberType(T1) is numberType(T2):
      asSimd(f(x[], y[]))
    elif numberType(T1) is float64:
      asSimd(f(x[], y[].toDoubleImpl))
    else:
      asSimd(f(x[].toDoubleImpl, y[]))

template f2s(f: untyped): untyped {.dirty.} =
  template f*(x: Simd, y: SomeNumber): untyped =
    mixin f
    asSimd(f(x[], y))
  template f*(x: SomeNumber, y: Simd): untyped =
    mixin f
    asSimd(f(x, y[]))
  f2(f)

f1(`-`)
f1(abs)
f1(inv)
f1(sqrt)
f1(rsqrt)
f1(sin)
f1(cos)
f1(acos)
f1(load1)
f2(atan2)
f2s(`+`)
f2s(`-`)
#f2s(`*`)
f2s(`/`)
f2s(add)
f2s(sub)
f2s(mul)
f2s(min)
f2s(max)

setBinop(`*`,mul,Simd,Simd2,asSimd(type(x[]*y[])))
setBinop(`*`,mul,Simd,SomeNumber,asSimd(type(x[]*y)))
setBinop(`*`,mul,SomeNumber,Simd2,asSimd(type(x*y[])))

# special cases

template getNc*(x: Simd): untyped = 1
template getNs*(x: Simd): untyped = 1

template to*[T](x: SomeNumber, y: typedesc[Simd[T]]): untyped =
  asSimd(x.to(type(T)))

template `+`*(x: Simd): untyped = x
template adj*(x: Simd): untyped = x

template norm2*[T](x: Simd[T]): untyped =
  mixin norm2
  when T is Masked:
    norm2(x[])
  elif T is Indexed:
    norm2(x[][])
  else:
    asSimd(norm2(x[]))

template trace*(x: Simd): untyped = x

template simdReduce*(x: Simd): untyped =
  mixin simdReduce
  simdReduce(x[])
template simdMaxReduce*(x: Simd): untyped =
  mixin simdMaxReduce
  simdMaxReduce(x[])
template simdMinReduce*(x: Simd): untyped =
  mixin simdMinReduce
  simdMinReduce(x[])
template simdSum*(x: Simd): untyped = simdReduce(x)
template simdMax*(x: Simd): untyped = simdMaxReduce(x)
template simdMin*(x: Simd): untyped = simdMinReduce(x)

template re*(x: Simd): untyped = x
template im*[T:Simd](x: T): untyped = 0.to(type(T))

template assign*(x: Simd, y: array) = assign(x[], y)
template `:=`*(x: Simd, y: array) = assign(x[], y)

template assign*(x: SomeNumber, y: Simd) = assign(x, y[])  # Masked
template `:=`*(x: SomeNumber, y: Simd) = assign(x, y[])  # Masked
template `+=`*(x: SomeNumber, y: Simd) =
  when y.simdLength == 1:
    x += y[0]
  else:
    x += y[]  # Masked

template exp*(xx: Simd[Indexed]): untyped =
  let x = xx
  exp(x[][x.indexedIdx])
