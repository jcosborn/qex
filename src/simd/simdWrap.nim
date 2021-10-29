import ../base/wrapperTypes
import ../maths/types
#import ../maths/complexNumbers
#export complexNumbers

makeWrapperType(Simd)
type
  Simd2*[T] = Simd[T]
  Simd3*[T] = Simd[T]

template `[]`*(x: Simd, i: typed): untyped =
  when i is Simd:
    #x[][i[]]
    indexed(x, i[])
  else:
    x[][i]
    #indexed(x, i)

template doIndexed[T](x: T): untyped =
  when T is Indexed:
    x[]
  else:
    x

template `[]=`*(x: Simd, i: typed, y: typed): untyped =
  when i is Simd:
    when y is Simd:
      x[][i[]] = doIndexed(y[])
    else:
      x[][i[]] = y
  else:
    when y is Simd:
      x[][i] = doIndexed(y[])
    else:
      x[][i] = y

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
  template f*[T](x: var Simd[T], y: Simd2) =
    #static: echo "f Simd Simd"
    #static: echo "  ", x.type
    #static: echo "  ", y.type
    mixin f
    #f(x[], y[])
    when numberType(T) is float64:
      f(x[], y[].toDoubleImpl)
    elif numberType(T) is float32:
      f(x[], y[].toSingleImpl)
    else:
      f(x[], y[])

template p2s(f: untyped) {.dirty.} =
  template f*(x: var Simd, y: SomeNumber) =
    #static: echo "f Simd Somenumber"
    #static: echo "  ", x.type
    #static: echo "  ", y.type
    mixin f
    f(x[], y)
  p2(f)

template p3(f: untyped) {.dirty.} =
  template f*(x: var Simd, y: Simd2, z: Simd3) =
    mixin f
    f(x[], y[], z[])

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
  template f*(x: Simd, y: Simd2): untyped =
    mixin f
    asSimd(f(x[], y[]))

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
f2s(`*`)
f2s(`/`)
f2s(add)
f2s(sub)
f2s(mul)
f2s(min)
f2s(max)


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

template simdSum*(x: Simd): untyped = simdSum(x[])
template simdMax*(x: Simd): untyped = simdMax(x[])
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
