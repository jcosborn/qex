import base/wrapperTypes
import maths/types
import maths

type
  Color*[T] = object
    v*: T
  Color2*[T] = Color[T]
  Color3*[T] = Color[T]

#template asColorX*(xx: typed): untyped =
#  #lets(x,xx):
#  #static: echo "asColor typed"
#  #dumpTree: xx
#  let x_asColor = xx
#  Color[type(x_asColor)](v: x_asColor)
#template asColorX*(xx: typed{nkObjConstr}): untyped =
#  static: echo "asColor typed{nkObjConstr}"
#  dumpTree: xx
#  lets(x,xx):
#    Color[type(x)](v: x)
#template asColor*(xx: typed): untyped = asColorX(normalizeAst(xx))
template asColor*(xx: typed): untyped =
  #lets(x,xx):
  #static: echo "asColor typed"
  #dumpTree: xx
  let x_asColor = xx
  Color[type(x_asColor)](v: x_asColor)

template isWrapper*(x: Color): untyped = true
template asWrapper*(x: Color, y: typed): untyped =
  #static: echo "asWrapper Color"
  #dumpTree: y
  asColor(y)
template asVarWrapper*(x: Color, y: typed): untyped =
  #static: echo "asVarWrapper Color"
  #var cy = asColor(y)
  #cy
  asVar(asColor(y))
#template `[]`*(x: Color): untyped = x.v
makeDeref(Color, x.T)
template `[]`*(x: Color, i: any): untyped = x[][i]
template `[]`*(x: Color, i: any, j: any): untyped = x[][i,j]
template `[]`*(x: Color, i: any, j: any, y: any): untyped =
  x[][i,j] = y
forwardFunc(Color, len)
forwardFunc(Color, nrows)
forwardFunc(Color, ncols)
forwardFunc(Color, numberType)
forwardFunc(Color, nVectors)
forwardFunc(Color, simdType)
forwardFunc(Color, simdLength)
template row*(x: Color, i: any): untyped =
  mixin row
  asColor(row(x[],i))
template setRow*(r: Color; x: Color2; i: int): untyped =
  setRow(r[], x[], i)

template binDDRet(fn,wr,T1,T2) =
  template fn*(x: T1, y: T2): untyped =
    wr(fn(x[], y[]))

binDDRet(`+`, asColor, Color, Color2)
binDDRet(`-`, asColor, Color, Color2)
binDDRet(`*`, asColor, Color, Color2)
binDDRet(`/`, asColor, Color, Color2)

template numberType*[T](x: typedesc[Color[T]]): untyped = numberType(T)
#template numNumbers*[T](x: typedesc[Color[T]]): untyped = numberType(T)
template numNumbers*(x: Color): untyped = numNumbers(x[])
template load1*(x: Color): untyped = asColor(load1(x[]))
template assign*(r: var Color, x: SomeNumber) =
  assign(r[], x)
template assign*(r: var Color, x: AsComplex) =
  assign(r[], x)
template assign*(r: var Color, x: Color2) =
  assign(r[], x[])
template `:=`*(r: var Color, x: SomeNumber) =
  `:=`(r[], x)
template `:=`*(r: var Color, x: Color2) =
  r[] := x[]
template `+=`*(r: var Color, x: Color2) =
  r[] += x[]
template `*=`*(r: var Color, x: SomeNumber) =
  `*=`(r[], x)
template iadd*(r: var Color, x: AsComplex) =
  iadd(r[], x)
template iadd*(r: var Color, x: Color2) =
  iadd(r[], x[])
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
template `*`*(x: Color, y: SomeNumber): untyped =
  asColor(x[] * y)
template `*`*(x: SomeNumber, y: Color2): untyped =
  asColor(x * y[])
template `*`*(x: AsComplex, y: Color2): untyped =
  asColor(x * y[])
template mul*(r: var Color, x: Color2, y: Color3) =
  mul(r[], x[], y[])
template mul*(r: var Color, x: SomeNumber, y: Color3) =
  mul(r[], x, y[])
template mul*(r: var Color, x: AsComplex, y: Color3) =
  mul(r[], x, y[])
template random*(x: var Color) =
  gaussian(x[], r)
template gaussian*(x: var Color, r: var untyped) =
  gaussian(x[], r)
template uniform*(x: var Color, r: var untyped) =
  uniform(x[], r)
template z4*(x: var Color, r: var untyped) =
  z4(x[], r)
template z2*(x: var Color, r: var untyped) =
  z2(x[], r)
template u1*(x: var Color, r: var untyped) =
  u1(x[], r)
template projectU*(r: var Color, x: Color2) =
  projectU(r[], x[])
template projectSU*(r: var Color, x: Color2) =
  projectSU(r[], x[])
template checkSU*(x: Color):untyped = checkSU(x[])
template norm2*(x: Color): untyped = norm2(x[])
template norm2*(r: var any, x: Color): untyped = norm2(r, x[])
template inorm2*(r: var any, x: Color2) = inorm2(r, x[])
template dot*(x: Color, y: Color2): untyped =
  dot(x[], y[])
template idot*(r: var any, x: Color2, y: Color3) = idot(r, x[], y[])
template redot*(x: Color, y: Color2): untyped =
  redot(x[], y[])
template trace*(x: Color): untyped = trace(x[])
