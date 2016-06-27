# unary ops: assign(=),neg(-),iadd(+=),isub(-=),imul(*=),idiv(/=)
# binary ops: add(+),sub(-),mul(*),divd(/),imadd(+=*),imsub(-=*)
# ternary ops: madd(*+),msub(*-),nmadd(-*+),nmsub(-*-)
# sqrt, rsqrt
# trace, norm2

# imadd(r=x*y+r),imsub(r=x*y-r),inmadd(r=-x*y+r),imsub(r=-x*y-r)

import math
export math

{.passL:"-lm".}

type
  SomeInteger2* = int|int8|int16|int32|int64
  SomeInteger3* = int|int8|int16|int32|int64
  SomeInteger4* = int|int8|int16|int32|int64
  SomeReal2* = float|float32|float64
  SomeReal3* = float|float32|float64
  SomeReal4* = float|float32|float64
  SomeNumber2* = SomeInteger2 | SomeReal2
  SomeNumber3* = SomeInteger3 | SomeReal3
  SomeNumber4* = SomeInteger4 | SomeReal4

var FLT_EPSILON*{.importC,header:"<float.h>".}:float32
var DBL_EPSILON*{.importC,header:"<float.h>".}:float64
template epsilon*(x:float32):untyped = FLT_EPSILON
template epsilon*(x:float64):untyped = DBL_EPSILON
template basicNumberDefines(T,N,F) {.dirty.} =
  template numberType*(x:T):expr = F
  template numberType*(x:typedesc[T]):expr = F
  template numNumbers*(x:T):expr = N
  template numNumbers*(x:typedesc[T]):expr = N
basicNumberDefines(float32, 1, float32)
basicNumberDefines(float64, 1, float64)
template `[]`*(x:SomeNumber; i:SomeInteger):expr = x

template cnvrt(r,x):expr = ((type(r))(x))
template to*(x:any; t:typedesc[SomeNumber]):expr =
  when x.type is t:
    x
  else:
    #var r{.noInit.}:t
    #assign(r, x)
    #r
    t(x)
template to*(t:typedesc[SomeNumber]; x:any):expr =
  when x.type is t:
    x
  else:
    t(x)

template assign*(x:var SomeNumber; y:ptr SomeNumber2):untyped =
  x = cnvrt(x,y[])

#template assign*(r:var SomeNumber, x:SomeNumber2):untyped =
proc assign*(r:var SomeNumber, x:SomeNumber2) {.inline.} =
  r = cnvrt(r,x)
template neg*(r:var SomeNumber, x:SomeNumber2):untyped =
  r = cnvrt(r,-x)
template iadd*(r:var SomeNumber, x:SomeNumber2):untyped =
  r += cnvrt(r,x)
template isub*(r:var SomeNumber, x:SomeNumber2):untyped =
  r -= cnvrt(r,x)
template imul*(r:var SomeNumber, x:SomeNumber2):untyped =
  r *= cnvrt(r,x)
template idiv*(r:var SomeNumber, x:SomeNumber2):untyped =
  r /= cnvrt(r,x)
template add*(r:var SomeNumber, x:SomeNumber2, y:SomeNumber3):untyped =
  r = cnvrt(r,x) + cnvrt(r,y)
template sub*(r:var SomeNumber, x:SomeNumber2, y:SomeNumber3):untyped =
  r = cnvrt(r,x) - cnvrt(r,y)
template mul*(r:var SomeNumber, x:SomeNumber2, y:SomeNumber3):untyped =
  r = cnvrt(r,x) * cnvrt(r,y)
template divd*(r:var SomeNumber, x:SomeNumber2, y:SomeNumber3):untyped =
  r = cnvrt(r,x) / cnvrt(r,y)
template imadd*(r:var SomeNumber, x:SomeNumber2, y:SomeNumber3):untyped =
  r += cnvrt(r,x) * cnvrt(r,y)
template imsub*(r:var SomeNumber, x:SomeNumber2, y:SomeNumber3):untyped =
  r -= cnvrt(r,x) * cnvrt(r,y)
template madd*(r:var SomeNumber, x:SomeNumber2,
               y:SomeNumber3, z:SomeNumber4):untyped =
  r = (cnvrt(r,x) * cnvrt(r,y)) + cnvrt(r,z)
template msub*(r:var SomeNumber, x:SomeNumber2,
               y:SomeNumber3, z:SomeNumber4):untyped =
  r = (cnvrt(r,x) * cnvrt(r,y)) - cnvrt(r,z)
template nmadd*(r:var SomeNumber, x:SomeNumber2,
                y:SomeNumber3, z:SomeNumber4):untyped =
  r = cnvrt(r,z) - (cnvrt(r,x) * cnvrt(r,y))
template nmsub*(r:var SomeNumber, x:SomeNumber2,
                y:SomeNumber3, z:SomeNumber4):untyped =
  r = cnvrt(r,-z) - (cnvrt(r,x) * cnvrt(r,y))
template conj*(r:var SomeNumber, x:SomeNumber2):untyped = assign(r, x)
template adj*(r:var SomeNumber, x:SomeNumber2):untyped = assign(r, x)
template trace*(x:SomeNumber):untyped = x
template norm2*(r:var SomeNumber, x:SomeNumber2):untyped = mul(r, x, x)
template norm2*(x:SomeNumber):expr = x*x
template inorm2*(r:var SomeNumber; x:SomeNumber2):untyped = imadd(r, x, x)
template dot*(x:SomeNumber; y:SomeNumber2):expr = x*y
template idot*(r:var SomeNumber; x:SomeNumber2;y:SomeNumber3):expr =
  imadd(r,x,y)
template simdSum*(x:SomeNumber):expr = x
template simdSum*(r:var SomeNumber; x:SomeNumber2):untyped =
 r = (type(r))(x)
template simdReduce*(x:SomeNumber):expr = x
template perm1*(r:var SomeNumber; x:SomeNumber2):untyped =
 r = (type(r))(x)
template perm2*(r:var SomeNumber; x:SomeNumber2):untyped =
 r = (type(r))(x)
template perm4*(r:var SomeNumber; x:SomeNumber2):untyped =
 r = (type(r))(x)
template perm8*(r:var SomeNumber; x:SomeNumber2):untyped =
 r = (type(r))(x)
#proc sqrt*(x:float32):float32 {.importC:"sqrtf",header:"math.h".}
#proc sqrt*(x:float64):float64 {.importC:"sqrt",header:"math.h".}
proc acos*(x:float64):float64 {.importC:"acos",header:"math.h".}
template rsqrt*(r:var SomeNumber; x:SomeNumber) =
  r = cnvrt(r,1)/sqrt(cnvrt(r,x))

template load1*(x:SomeNumber):expr = x

template tmpvar*(r:untyped; x:untyped):untyped =
  mixin load1
  var r{.noInit.}:type(load1(x))
template load2*(r:untyped, x:untyped):untyped =
  mixin load1,assign
  #tmpvar(r, x)
  var r{.noInit.}:type(load1(x))
  assign(r, x)
template store*(r:var untyped, x:untyped):untyped =
  mixin assign
  assign(r, x)

template load*(x:untyped):expr =
  mixin load1
  load1(x)
template load*(r:untyped, x:untyped):untyped =
  mixin load2
  load2(r, x)

template `:=`*(x:var SomeNumber; y:SomeNumber2):untyped = assign(x,y)
template `+`*(x:SomeReal; y:SomeInteger):auto = x + cnvrt(x,y)
template `+`*(x:SomeInteger; y:SomeReal):auto = cnvrt(y,x) + y
template `-`*(x:SomeReal; y:SomeInteger):auto = x - cnvrt(x,y)
template `-`*(x:SomeInteger; y:SomeReal):auto = cnvrt(y,x) - y
template `*`*(x:SomeInteger; y:SomeReal):auto = cnvrt(y,x) * y

when isMainModule:
  var d1,d2:float
  var s1,s2:float32
  var i1,i2:int
  assign(d1,s1)
  assign(d1,i1)
  imadd(d1, s1, i1)
  load(t, d1)
  madd(t, s2, i1, i2)
  store(d1, t)
  echo d1
