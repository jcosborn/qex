import ../base/wrapperTypes

makeWrapperType(Simd)
type
  Simd2*[T] = Simd[T]
  Simd3*[T] = Simd[T]

template `[]`*(x: Simd, i: typed): untyped = x[][i]
template `[]=`*(x: Simd, i: typed, y: typed): untyped = x[][i] = y

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

# no return value
template p2(f: untyped) {.dirty.} =
  template f*(x: var Simd, y: Simd2) =
    #static: echo "f Simd Simd"
    #static: echo "  ", x.type
    #static: echo "  ", y.type
    mixin f
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

p2(rsqrt)
p2s(assign)
p2s(`:=`)
p2s(`+=`)
p2s(`-=`)
p2s(iadd)
p2s(inorm2)
p3(imadd)
p3(imsub)


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

template `+`*(x: Simd): untyped = x

template norm2*[T](x: Simd[T]): untyped =
  when T is Masked:
    norm2(x[])
  else:
    asSimd(norm2(x[]))

template trace*(x: Simd): untyped = x

template simdSum*(x: Simd): untyped = simdSum(x[])
template simdMax*(x: Simd): untyped = simdMax(x[])

template assign*(x: Simd, y: array) = assign(x[], y)
template `:=`*(x: Simd, y: array) = assign(x[], y)

template assign*(x: SomeNumber, y: Simd) = assign(x, y[])  # Masked
template `:=`*(x: SomeNumber, y: Simd) = assign(x, y[])  # Masked
template `+=`*(x: SomeNumber, y: Simd) = x += y[]  # Masked

