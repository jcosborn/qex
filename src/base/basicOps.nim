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

template getOptimPragmas* =
  {.pragma: alwaysInline, inline,
    codegenDecl: "inline __attribute__((always_inline)) $# $#$#".}
getOptimPragmas()

type
  SomeInteger2* = int|int8|int16|int32|int64
  SomeInteger3* = int|int8|int16|int32|int64
  SomeInteger4* = int|int8|int16|int32|int64
  SomeFloat2* = float|float32|float64
  SomeFloat3* = float|float32|float64
  SomeFloat4* = float|float32|float64
  SomeNumber2* = SomeInteger2 | SomeFloat2
  SomeNumber3* = SomeInteger3 | SomeFloat3
  SomeNumber4* = SomeInteger4 | SomeFloat4

var FLT_EPSILON*{.importC,header:"<float.h>".}:float32
var DBL_EPSILON*{.importC,header:"<float.h>".}:float64
template epsilon*(x:float32):untyped = FLT_EPSILON
template epsilon*(x:typedesc[float32]):untyped = FLT_EPSILON
template epsilon*(x:float64):untyped = DBL_EPSILON
template epsilon*(x:typedesc[float64]):untyped = DBL_EPSILON
template basicNumberDefines(T,N,F) {.dirty.} =
  template numberType*(x:T):untyped = F
  template numberType*(x:typedesc[T]):untyped = F
  template numNumbers*(x:T):untyped = N
  template numNumbers*(x:typedesc[T]):untyped = N
basicNumberDefines(float32, 1, float32)
basicNumberDefines(float64, 1, float64)

template numberType*[T](x:ptr T):untyped = numberType(T)
template numberType*[T](x:tuple[re,im:T]):untyped = numberType(T)
template numberType*[T](x:typedesc[tuple[re,im:T]]):untyped = numberType(T)
template numberType*[I,T](x:array[I,T]):untyped = numberType(type(T))
template numberType*[I,T](x:typedesc[array[I,T]]):untyped = numberType(type(T))
#template numberType*(x:not typedesc):untyped = numberType(type(x))
template `[]`*[T](x:typedesc[ptr T]):untyped = T
template `[]`*(x:SomeNumber; i:SomeInteger):untyped = x
template isWrapper*(x: SomeNumber): untyped = false
template isWrapper*(x: typedesc[SomeNumber]): untyped = false
template eval*[T:SomeNumber](x: typedesc[T]): untyped = typeof(T)

template cnvrt(r,x):untyped = ((type(r))(x))
template to*(x:auto; t:typedesc[SomeNumber]):untyped =
  when x.type is t:
    x
  else:
    #var r{.noInit.}:t
    #assign(r, x)
    #r
    t(x)
template to*(t:typedesc[SomeNumber]; x:auto):untyped =
  when x.type is t:
    x
  else:
    t(x)
template toDoubleImpl*(x:SomeNumber):untyped =
  when type(x) is float64:
    x
  else:
    float64(x)

template assign*[R,X:SomeNumber](r: R; x: ptr X) =
  r = R(x[])
template assign*[R,X:SomeNumber](r: R, x: X) =
  r = R(x)
template `:=`*[R,X:SomeNumber](r: R; x: X) =
  r = R(x)
template `:=`*[R,X:SomeNumber](r: R; x: ptr X) =
  r = R(x[])
#proc `+=`*[R,X:SomeNumber](r: var R; x: X) {.alwaysInline.} =
#  r = r + R(x)
#proc `-=`*[R,X:SomeNumber](r: var R; x: X) {.alwaysInline.} =
#  r = r - R(x)

template adj*(x: SomeNumber): untyped = x
template inv*[T:SomeNumber](x: T): untyped = (T(1))/x

template neg*[R,X:SomeNumber](r: R, x: X) =
  r := -x
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
template re*(x:SomeNumber): untyped = x
template im*[T:SomeNumber](x:T): untyped = T(0)
template conj*(r:var SomeNumber, x:SomeNumber2) = assign(r, x)
template adj*(r:var SomeNumber, x:SomeNumber2) = assign(r, x)
template trace*(x:SomeNumber):untyped = x
template norm2*(r:var SomeNumber, x:SomeNumber2):untyped =
  let tNorm2VNum = x
  mul(r, tNorm2VNum, tNorm2VNum)
template norm2*(x:SomeNumber):untyped =
  let tNorm2Num = x
  tNorm2Num * tNorm2Num
template inorm2*(r:var SomeNumber; x:SomeNumber2):untyped =
  let tInorm2Num = x
  imadd(r, tInorm2Num, tInorm2Num)
template dot*(x:SomeNumber; y:SomeNumber2):untyped = x*y
template idot*(r:var SomeNumber; x:SomeNumber2;y:SomeNumber3):untyped =
  imadd(r,x,y)
template redot*(x:SomeNumber; y:SomeNumber2):untyped = x*y
template redotinc*(r:var SomeNumber; x:SomeNumber2; y:SomeNumber3):untyped =
  r += x*y
template simdLength*(x:SomeNumber):untyped = 1
template simdLength*(x:typedesc[SomeNumber]):untyped = 1
template simdReduce*(x:SomeNumber):untyped = x
template simdReduce*[T,X:SomeNumber](r:var T; x:X) =
  r = T(x)
template simdMaxReduce*(x:SomeNumber):untyped = x
template simdMinReduce*(x:SomeNumber):untyped = x
template simdSum*(x:SomeNumber):untyped = simdReduce(x)
template simdSum*[T,X:SomeNumber](r:var T; x:X) = simdReduce(r,x)
template simdMax*(x:SomeNumber):untyped = simdMaxReduce(x)
template simdMin*(x:SomeNumber):untyped = simdMinReduce(x)
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
template rsqrt*[R,X:SomeNumber](r:var R; x:X) =
  r = R(1)/sqrt(R(x))
template rsqrt*(x: SomeNumber): untyped = 1/sqrt(x)
template select*(c: bool, a,b: typed): untyped =
  if c: a else: b

template load1*(x:SomeNumber):untyped = x

template tmpvar*(r:untyped; x:untyped):untyped =
  mixin load1
  var r{.noInit.}:evalType(load1(x))
template load2*(r:untyped, x:untyped):untyped =
  mixin load1,assign
  #tmpvar(r, x)
  var r{.noInit.}:evalType(load1(x))
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

template `+`*(x:SomeFloat; y:SomeInteger):auto = x + cnvrt(x,y)
template `+`*(x:SomeInteger; y:SomeFloat):auto = cnvrt(y,x) + y
template `-`*(x:SomeFloat; y:SomeInteger):auto = x - cnvrt(x,y)
template `-`*(x:SomeInteger; y:SomeFloat):auto = cnvrt(y,x) - y
template `*`*[T:SomeFloat](x:SomeInteger; y:T):auto = (T(x)) * y
template `*`*[T:SomeFloat](x:T; y:SomeInteger):auto = x * (T(y))
template `/`*[T:SomeFloat](x:SomeInteger,y:T):auto = (T(x)) / y
template `/`*[T:SomeFloat](x:T,y:SomeInteger):auto = x / (T(y))

template `:=`*[T](x: SomeNumber; y: array[1,T]) = assign(x,y[0])

template setUnopP*(op,fun,t1,t2: untyped): untyped {.dirty.} =
  proc op*(x: t1): auto {.alwaysInline,noInit.} =
    var r{.noInit.}: t2
    fun(r, x)
    r
template setUnopT*(op,fun,t1,t2: untyped): untyped {.dirty.} =
  template op*(x: t1): untyped =
    var rSetUnopT{.noInit.}: t2
    fun(rSetUnopT, x)
    rSetUnopT

template setBinopP*(op,fun,t1,t2,t3: untyped): untyped {.dirty.} =
  #template op*(x: typedesc[t1]; y: typedesc[t2]): typedesc = t3
  proc op*(x: t1; y: t2): auto {.alwaysInline,noInit.} =
    var r{.noInit.}: t3
    fun(r, x, y)
    r
template setBinopT*(op,fun,t1,t2,t3: untyped) {.dirty.} =
  #template op*(x: typedesc[t1]; y: typedesc[t2]): typedesc = t3
  template op*(x: t1; y: t2): untyped =
    var rSetBinopT{.noInit.}: t3
    fun(rSetBinopT, x, y)
    rSetBinopT

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

import numberWrap
export numberWrap

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
