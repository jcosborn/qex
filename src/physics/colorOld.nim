import base/basicOps
import base/wrapperTypes
export wrapperTypes
import maths/types
import maths
import simd/simdWrap

makeWrapperType(Color):
  ## wrapper type for colored objects

type
  Color2*[T] = Color[T]
  Color3*[T] = Color[T]
  Color4*[T] = Color[T]

template asVarWrapper*(x: Color, y: typed): untyped =
  #static: echo "asVarWrapper Color"
  #var cy = asColor(y)
  #cy
  asVar(asColor(y))

template index*[T,I](x: typedesc[Color[T]], i: typedesc[I]): typedesc =
  when I is Color:
    index(T.type, I.type[])
  elif I.isWrapper:
    Color[index(T.type, I.type)]
  else:
    index(T.type, I.type)

template `[]`*[T](x: Color, i: T): untyped =
  when T is Color:
    x[][i[]]
  elif T.isWrapper:
    #indexed(x, i)
    var tColorBracket = asColor(x[][i])
    tColorBracket
  else:
    x[][i]
template `[]`*(x: Color, i,j: SomeInteger): auto = x[][i,j]

#[
template `[]=`*[T](x: Color, i:T, y: typed): untyped =
  when T.isWrapper and T isnot Color:
    x[][asScalar(i)] = y
  else:
    x[][i] = y
]#
template `[]=`*[T](x: Color, i: T; y: auto) =
  when T is Color2:
    x[][i[]] = y
  elif T.isWrapper:
    #indexed(x, i)
    var tColorBracket = asColor(x[][i])
    tColorBracket := y
  else:
    x[][i] = y
template `[]=`*(x: Color, i,j: SomeInteger, y: auto) =
  x[][i,j] = y

# forward from value to value
template forwardVV(f: untyped) {.dirty.} =
  template f*(x: Color): auto =
    mixin f
    f(x[])
# forward from value to type
#template forwardVT(f: untyped) {.dirty.} =
#  template f*[T](x: Color[T]): untyped =
#    mixin f
#    f(type T)
# forward from type to type
template forwardTT(f: untyped) {.dirty.} =
  template f*[T](x: typedesc[Color[T]]): auto =
    mixin f
    f(type T)
# forward from type to type and wrap
template forwardTTW(f: untyped) {.dirty.} =
  template f*[T](x: typedesc[Color[T]]): auto =
    mixin f
    Color[f(type T)]

forwardVV(len)
forwardVV(nrows)
forwardVV(ncols)
forwardVV(nVectors)
forwardVV(simdType)
#forwardVV(simdLength)
forwardVV(getNs)
forwardVV(numNumbers)

forwardTT(len)
forwardTT(nrows)
forwardTT(ncols)
forwardTT(nVectors)
forwardTT(simdType)
forwardTT(simdLength)
forwardTT(getNs)
forwardTT(numberType)

forwardTTW(toSingle)
forwardTTW(toDouble)

#template eval*[T](x: typedesc[Color[T]]): typedesc = asColor(eval(type T))
template eval*[T:Color](x: typedesc[T]): typedesc = asColor(eval((type T)[]))

template has*[T:Color](x: typedesc[T], y: typedesc): bool =
  mixin has
  when y is Color: true
  else: has(T.type[], y)

template row*(x: Color, i: untyped): untyped =
  mixin row
  asColor(row(x[],i))
template setRow*(r: Color; x: Color2; i: untyped): untyped =
  setRow(r[], x[], i)

template getNc*[T](x: Color[T]): untyped =
  when T is Mat1:
    x[].nrows
  elif T is Vec1:
    x[].len
  else:
    static:
      echo "error: unknown Nc"
      echo x.repr
      echo type(x).name
      qexExit 1
    0

template binDDRet(fn,wr,T1,T2) =
  template fn*(x: T1, y: T2): untyped =
    wr(fn(x[], y[]))
  #template fn*(x: T1, y: T2): auto =
  #  var tmp {.noInit.}: wr(evalType(fn(x[],y[]))
  #  wr(fn(x[], y[]))

#binDDRet(`+`, asColor, Color, Color2)
binDDRet(`-`, asColor, Color, Color2)
#binDDRet(`*`, asColor, Color, Color2)
binDDRet(`/`, asColor, Color, Color2)

setBinop(`+`, add, Color, Color2, asColor(evalType(x[]+y[])))
setBinop(`*`, mul, Color, Color2, asColor(evalType(x[]*y[])))

template load1*(x: Color): untyped = asColor(load1(x[]))
template `-`*(x: Color): untyped = asColor(-(x[]))
template assign*(r: var Color, x: SomeNumber) =
  assign(r[], x)
template assign*(r: var Color, x: AsComplex) =
  assign(r[], x)
template assign*(r: var Color, x: Color2) =
  assign(r[], x[])
template `:=`*(r: var Color, x: SomeNumber) =
  r[] := x
template `:=`*(r: var Color, x: AsComplex) =
  r[] := x
template `:=`*(r: var Color, x: Color2) =
  r[] := x[]
template `+=`*(r: var Color, x: Color2) =
  r[] += x[]
template `-=`*(r: var Color, x: Color2) =
  r[] -= x[]
template `*=`*(r: var Color, x: SomeNumber) =
  r[] *= x
template iadd*(r: var Color, x: SomeNumber) =
  iadd(r[], x)
template iadd*(r: var Color, x: AsComplex) =
  iadd(r[], x)
template iadd*(r: var Color, x: Color2) =
  iadd(r[], x[])
template isub*(r: var Color, x: SomeNumber) =
  isub(r[], x)
template isub*(r: var Color, x: Color2) =
  isub(r[], x[])
template imul*(r: var Color, x: SomeNumber) =
  imul(r[], x)
template imadd*(r: var Color, x: Color2, y: Color3) =
  imadd(r[], x[], y[])
template imadd*(r: var Color, x: AsComplex, y: Color3) =
  imadd(r[], x, y[])
template imsub*(r: var Color, x: Color2, y: Color3) =
  imsub(r[], x[], y[])
template peqOuter*(r: var Color, x: Color2, y: Color3) =
  peqOuter(r[], x[], y[])
template meqOuter*(r: var Color, x: Color2, y: Color3) =
  meqOuter(r[], x[], y[])
template `+`*(r: Color, x: SomeNumber): untyped =
  asColor(r[] + x)
template `+`*(r: SomeNumber, x: Color): untyped =
  asColor(r + x[])
template `-`*(r: Color, x: SomeNumber): untyped =
  asColor(r[] - x)
template `+`*(r: Color, x: AsComplex): untyped =
  asColor(r[] + x)
template `-`*(r: Color, x: AsComplex): untyped =
  asColor(r[] - x)
template add*(r: var Color, x: Color2, y: Color3) =
  add(r[], x[], y[])
template sub*(r: var Color, x: Color2, y: Color3) =
  sub(r[], x[], y[])
template `*`*(x: Color, y: SomeNumber): untyped =
  asColor(x[] * y)
#template `*`*(x: SomeNumber, y: Color): auto =
  #asColor(x * y[])
template `*`*[X:SomeNumber,Y:Color](x: X, y: Y): auto =
  #static: echo "SomeNumber * Color"
  #var tmp {.noInit.}: asColor(type(X)*type(Y)[])
  var tmp {.noInit.}: asColor(evalType(x*y[]))
  #static: echo $type(tmp)
  mul(tmp[], x, y[])
  tmp
template `*`*(x: Simd, y: Color2): untyped =
  asColor(x * y[])
#template `*`*(x: AsReal, y: Color2): untyped =
#  asColor(x * y[])
template `*`*[X:AsReal,Y:Color](x: X, y: Y): auto =
  #var tmp {.noInit.}: asColor(type(X)*type(Y)[])
  var tmp {.noInit.}: asColor(evalType(x*y[]))
  mul(tmp[], x, y[])
  tmp
template `*`*(x: AsImag, y: Color2): untyped =
  asColor(x * y[])
template `*`*(x: AsComplex, y: Color2): untyped =
  asColor(x * y[])
template `*`*(x: Color, y: AsComplex): untyped =
  asColor(x[] * y)
#template mul*(x: SomeNumber, y: Color2): untyped =
#  asColor(`*`(x, y[]))
template mul*[X:SomeNumber,Y:Color](x: X, y: Y): auto =
  #var tmp {.noInit.}: asColor(type(X)*type(Y)[])
  var tmp {.noInit.}: asColor(evalType(x*y[]))
  mul(tmp[], x, y[])
  tmp
template mul*(x: Color, y: Color2): untyped =
  asColor(`*`(x[], y[]))
template mul*(r: var Color, x: Color2, y: Color3) =
  mul(r[], x[], y[])
template mul*(r: var Color, x: SomeNumber, y: Color3) =
  mul(r[], x, y[])
template mul*(r: var Color, x: Color2, y: SomeNumber) =
  mul(r[], x[], y)
template mul*(r: var Color, x: AsComplex, y: Color3) =
  mul(r[], x, y[])
template mul*(x: AsComplex, y: Color2): untyped =
  asColor(mul(x, y[]))
template random*(x: var Color) =
  gaussian(x[], r)
template gaussian*(x: Color, r: untyped) =
  gaussian(x[], r)
template uniform*(x: var Color, r: var untyped) =
  uniform(x[], r)
template z4*(x: var Color, r: var untyped) =
  z4(x[], r)
template z2*(x: var Color, r: var untyped) =
  z2(x[], r)
template u1*(x: var Color, r: var untyped) =
  u1(x[], r)
template projectU*(r: var Color) =
  projectU(r[])
template projectU*(r: var Color, x: Color2) =
  projectU(r[], x[])
template projectUderiv*(r: var Color, u: Color2, x: Color3, chain: Color4) =
  projectUderiv(r[], u[], x[], chain[])
template projectUderiv*(r: var Color, x: Color3, chain: Color4) =
  projectUderiv(r[], x[], chain[])
template projectSU*(r: var Color) =
  projectSU(r[])
template projectSU*(r: var Color, x: Color2) =
  projectSU(r[], x[])
template projectTAH*(r: var Color) =
  projectTAH(r[])
template projectTAH*(r: var Color, x: Color2) =
  projectTAH(r[], x[])
template checkU*(x: Color):untyped = checkU(x[])
template checkSU*(x: Color):untyped = checkSU(x[])
template norm2*(x: Color): untyped = norm2(x[])
template norm2*(r: var auto, x: Color): untyped = norm2(r, x[])
template inorm2*(r: var auto, x: Color2) = inorm2(r, x[])
template dot*(x: Color, y: Color2): untyped =
  dot(x[], y[])
template idot*(r: var auto, x: Color2, y: Color3) = idot(r, x[], y[])
template redot*(x: Color, y: Color2): untyped =
  redot(x[], y[])
template trace*(x: Color): untyped = trace(x[])
template simdSum*(x: Color): untyped = asColor(simdSum(x[]))
template re*(x: Color): untyped = asColor(re(x[]))
template im*(x: Color): untyped = asColor(im(x[]))
template exp*(x: Color): untyped = asColor(exp(x[]))
template expDeriv*(x: Color, c: Color2): untyped = asColor(expDeriv(x[], c[]))
template ln*(x: Color): untyped = asColor(ln(x[]))
