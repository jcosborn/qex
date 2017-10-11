# unary ops: assign(=),neg(-),iadd(+=),isub(-=),imul(*=),idiv(/=)
# binary ops: add(+),sub(-),mul(*),divd(/),imadd(+=*),imsub(-=*)
# ternary ops: madd(*+),msub(*-),nmadd(-*+),nmsub(-*-)
# sqrt, rsqrt
# trace, norm2

# imadd(r=x*y+r),imsub(r=x*y-r),inmadd(r=-x*y+r),imsub(r=-x*y-r)

import globals
import math
export math
import macros

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
  template numberType*(x:T):untyped = F
  template numberType*(x:typedesc[T]):untyped = F
  template numNumbers*(x:T):untyped = N
  template numNumbers*(x:typedesc[T]):untyped = N
basicNumberDefines(float32, 1, float32)
basicNumberDefines(float64, 1, float64)

template numberType*[T](x:tuple[re,im:T]):untyped = numberType(T)
template numberType*[T](x:typedesc[tuple[re,im:T]]):untyped = numberType(T)
template numberType*[I,T](x:array[I,T]):untyped = numberType(type(T))
template numberType*[I,T](x:typedesc[array[I,T]]):untyped = numberType(T)
#template numberType*(x:not typedesc):untyped = numberType(type(x))
template `[]`*(x:SomeNumber; i:SomeInteger):untyped = x

template cnvrt(r,x):untyped = ((type(r))(x))
template to*(x:any; t:typedesc[SomeNumber]):untyped =
  when x.type is t:
    x
  else:
    #var r{.noInit.}:t
    #assign(r, x)
    #r
    t(x)
template to*(t:typedesc[SomeNumber]; x:any):untyped =
  when x.type is t:
    x
  else:
    t(x)
template toDoubleImpl*(x:SomeNumber):untyped =
  when type(x) is float64:
    x
  else:
    float64(x)

template assign*(x:var SomeNumber; y:ptr SomeNumber2):untyped =
  x = cnvrt(x,y[])
template assign*(r:var SomeNumber, x:SomeNumber2):untyped =
#proc assign*(r:var SomeNumber, x:SomeNumber2) {.inline.} =
  r = cnvrt(r,x)

#template adj*(x: SomeNumber): untyped = x
template inv*(x: SomeNumber): untyped = (type(x)(1))/x

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
template norm2*(x:SomeNumber):untyped = x*x
template inorm2*(r:var SomeNumber; x:SomeNumber2):untyped = imadd(r, x, x)
template dot*(x:SomeNumber; y:SomeNumber2):untyped = x*y
template idot*(r:var SomeNumber; x:SomeNumber2;y:SomeNumber3):untyped =
  imadd(r,x,y)
template redot*(x:SomeNumber; y:SomeNumber2):untyped = x*y
template redotinc*(r:var SomeNumber; x:SomeNumber2; y:SomeNumber3):untyped =
  r += x*y
template simdSum*(x:SomeNumber):untyped = x
template simdSum*(r:var SomeNumber; x:SomeNumber2):untyped =
 r = (type(r))(x)
template simdReduce*(x:SomeNumber):untyped = x
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
proc atan2*(x,y:float64):float64 {.importC:"atan2",header:"math.h".}
proc atan2*(x,y:float32):float32 {.importC:"atan2f",header:"math.h".}
template rsqrt*(r:var SomeNumber; x:SomeNumber) =
  r = cnvrt(r,1)/sqrt(cnvrt(r,x))

template load1*(x:SomeNumber):untyped = x

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

template load*(x:untyped):untyped =
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

template setUnopP*(op,fun,t1,t2: untyped): untyped {.dirty.} =
  proc op*(x: t1): auto {.inline,noInit.} =
    var r{.noInit.}: t2
    fun(r, x)
    r
template setUnopT*(op,fun,t1,t2: untyped): untyped {.dirty.} =
  template op*(xx: t1): untyped =
    #subst(xt,xx,r,_):
    #  lets(x,xt):
    subst(r,_):
      lets(x,xx):
        var r{.noInit.}: t2
        fun(r, x)
        r

template setBinopP*(op,fun,t1,t2,t3: untyped): untyped {.dirty.} =
  proc op*(x: t1; y: t2): auto {.inline,noInit.} =
    var r{.noInit.}: t3
    fun(r, x, y)
    r
template setBinopT*(op,fun,t1,t2,t3: untyped): untyped {.dirty.} =
  subst(r_setBinopT,_):
    template op*(xx: t1; yy: t2): untyped =
      #dumpTree: setBinopT op
      #echoType: xx
      #echoType: yy
      lets(x,xx,y,yy):
        var r_setBinopT{.noInit.}: t3
        fun(r_setBinopT, x, y)
        r_setBinopT

when forceInline:
  template setUnop*(op,fun,t1,t2: untyped): untyped {.dirty.} =
    setUnopT(op, fun, t1, t2)
  template setBinop*(op,fun,t1,t2,t3: untyped): untyped {.dirty.} =
    setBinopT(op, fun, t1, t2, t3)
else:
  template setUnop*(op,fun,t1,t2: untyped): untyped {.dirty.} =
    setUnopP(op, fun, t1, t2)
  template setBinop*(op,fun,t1,t2,t3: untyped): untyped {.dirty.} =
    setBinopP(op, fun, t1, t2, t3)

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
