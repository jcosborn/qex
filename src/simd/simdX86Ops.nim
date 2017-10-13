{. deadCodeElim: on .}

import simdX86Types
import simdSse
import simdAvx
import simdAvx512
#import ../basicOps
import base
import math
import macros

template binaryMixed(T,op1,op2:untyped):untyped =
  template op1*(x:SomeNumber; y:T):T = op2(x.to(T),y)
  template op1*(x:T; y:SomeNumber):T = op2(x,y.to(T))
template unaryMixedVar(T,op1,op2:untyped):untyped =
  template op1*(r: T; x:SomeNumber) = op2(r,x.to(T))
template binaryMixedVar(T,op1,op2:untyped):untyped =
  template op1*(r: T; x:SomeNumber; y:T) = op2(r,x.to(T),y)
  template op1*(r: T; x:T; y:SomeNumber) = op2(r,x,y.to(T))
template trinaryMixedVar(T,op1,op2:untyped):untyped =
  template op1*(r: T; x:SomeNumber; y:T; z:T) = op2(r,x.to(T),y,z)
  template op1*(r: T; x:T; y:SomeNumber; z:T) = op2(r,x,y.to(T),z)
  template op1*(r: T; x:T; y:T; z:SomeNumber) = op2(r,x,y,z.to(T))
template map1(T,N,op:untyped):untyped {.dirty.} =
  proc op*(x:T):T {.inline,noInit.} =
    let t = x.toArray
    var r{.noInit.}:type(t)
    for i in 0..<N:
      r[i] = op(t[i])
    assign(result, r)
template map2(T,N,op:untyped):untyped {.dirty.} =
  proc op*(x:T, y:T):T {.inline,noInit.} =
    let
      t = x.toArray
      u = y.toArray
    var r{.noInit.}:type(t)
    for i in 0..<N:
      r[i] = op(t[i], u[i])
    assign(result, r)

template basicDefs(T,F,N,P,S:untyped):untyped {.dirty.} =
  template numberType*(x:typedesc[T]):typedesc = F
  template numberType*(x:T):typedesc = F
  template numNumbers*(x:typedesc[T]):untyped = N
  template numNumbers*(x:T):untyped = N
  template simdType*(x:typedesc[T]):typedesc = T
  template simdType*(x:T):typedesc = T
  template simdLength*(x:T):untyped = N
  template simdLength*(x:typedesc[T]):untyped = N
  template load1*(x:T):untyped = x
  proc assign*(r:ptr F; x:T) {.inline.} =
    `P "_storeu_" S`(r, x)
  proc assign*(r:var T; x:ptr SomeNumber) {.inline.} =
    when x[] is F:
      r = `P "_loadu_" S`(x)
    else:
      let y = cast[ptr array[N,type(x[])]](x)
      var t{.noInit.}:array[N,F]
      for i in 0..<N: t[i] = F(y[][i])
      assign(r, t)
  template toSimd*(x:array[N,F]):untyped =
    `P "_loadu_" S`(cast[ptr F](unsafeAddr(x)))
  proc toArray*(x:T):array[N,F] {.inline,noInit.} =
    `P "_storeu_" S`(cast[ptr F](result.addr), x)
  template to*(x:SomeNumber; t:typedesc[T]):untyped =
    `P "_set1_" S`(F(x))
  template to*(x:array[N,F]; t:typedesc[T]):untyped =
    toSimd(x)
  proc to*(x:T; t:typedesc[array[N,F]]):array[N,F] {.inline,noInit.} =
    `P "_storeu_" S`(cast[ptr F](result.addr), x)
  proc assign1*(r:var T; x:SomeNumber) {.inline.} =
    r = `P "_set1_" S`(F(x))
  template setX:untyped = `P "_setr_" S`()
  template setF(x):untyped = F(x)
  macro assign*(r:var T; x:varargs[SomeNumber]):auto =
    if x.len==1:
      result = newCall(!"assign1", r, x[0])
    else:
      result = newStmtList()
      var call = getAst(setX())[0]
      for i in 0..<N:
        call.add getAst(setF(x[i mod x.len]))
      template asgn(r,c:untyped):untyped = r = c
      result = getAst(asgn(r, call))
    #echo result.treerepr
  proc assign*(r:var T; x:array[N,SomeNumber]) {.inline.} =
    when x[0] is F:
      r = `P "_loadu_" S`(cast[ptr F](unsafeAddr(x)))
    else:
      var t{.noInit.}:array[N,F]
      for i in 0..<N: t[i] = F(x[i])
      assign(r, t)
  #proc assign*(r:var T; x:T) {.inline.} =
  #  r = x
  #template `=`*(r: var T; x: T) = {.emit: [r, " = ", x].}
  template assign*(r: T; x:T) =
    r = x
  #proc assign*(r:var array[N,F]; x:T) {.inline.} =
  #  assign(cast[ptr F](r.addr), x)
  proc assign*(r:var array[N,F]; x:T) {.inline.} =
    assign(r[0].addr, x)
  proc `[]`*(x:T; i:SomeInteger):F {.inline,noInit.} =
    toArray(x)[i]
  proc `[]=`*(r:var T; i:SomeInteger; x:SomeNumber) {.inline,noInit.} =
    var a = toArray(r)
    a[i] = F(x)
    assign(r, a)
  proc `$`*(x:T):string =
    result = "[" & $x[0]
    for i in 1..<N:
      result &= "," & $x[i]
    result &= "]"
  proc prefetch*(x:ptr T) {.inline.} =
    mm_prefetch(cast[cstring](x), 3)

  template add*(x,y:T):T = `P "_add_" S`(x,y)
  template sub*(x,y:T):T = `P "_sub_" S`(x,y)
  template mul*(x,y:T):T = `P "_mul_" S`(x,y)
  template divd*(x,y:T):T = `P "_div_" S`(x,y)
  template neg*(x:T):T = sub(`P "_setzero_" S`(), x)

  binaryMixed(T, add, add)
  binaryMixed(T, sub, sub)
  binaryMixed(T, mul, mul)
  binaryMixed(T, divd, divd)

  template neg*(r:var T; x:T) = r = neg(x)
  template add*(r: T; x,y:T) = r = add(x,y)
  template sub*(r: T; x,y:T) = r = sub(x,y)
  template mul*(r: T; x,y:T) = r = mul(x,y)
  template divd*(r: T; x,y:T) = r = divd(x,y)

  binaryMixedVar(T, add, add)
  binaryMixedVar(T, sub, sub)
  binaryMixedVar(T, mul, mul)
  binaryMixedVar(T, divd, divd)

  template iadd*(r: T; x:T) = add(r,r,x)
  template isub*(r: T; x:T) = sub(r,r,x)
  template imul*(r: T; x:T) = mul(r,r,x)
  template idiv*(r: T; x:T) = divd(r,r,x)
  template imadd*(r: T; x,y: T) = iadd(r,mul(x,y))
  template imsub*(r: T; x,y:T) = isub(r,mul(x,y))
  template madd*(r: T; x,y,z:T) = add(r,mul(x,y),z)
  template msub*(r: T; x,y,z:T) = sub(r,mul(x,y),z)

  unaryMixedVar(T, iadd, iadd)
  unaryMixedVar(T, isub, isub)
  unaryMixedVar(T, imul, imul)
  unaryMixedVar(T, idiv, idiv)
  binaryMixedVar(T, imadd, imadd)
  binaryMixedVar(T, imsub, imsub)
  trinaryMixedVar(T, madd, madd)
  trinaryMixedVar(T, msub, msub)

  template `-`*(x:T):T = neg(x)
  template `+`*(x,y:T):T = add(x,y)
  template `-`*(x,y:T):T = sub(x,y)
  template `*`*(x,y:T):T = mul(x,y)
  template `/`*(x,y:T):T = divd(x,y)

  binaryMixed(T, `+`, add)
  binaryMixed(T, `-`, sub)
  binaryMixed(T, `*`, mul)
  binaryMixed(T, `/`, divd)

  template `:=`*(r: T; x:T) = assign(r,x)
  template `+=`*(r: T, x:T) = iadd(r,x)
  template `-=`*(r: T, x:T) = isub(r,x)
  template `*=`*(r: T, x:T) = imul(r,x)
  template `/=`*(r: T, x:T) = idiv(r,x)

  unaryMixedVar(T, `:=`, assign)
  template `:=`*(r: T; x:openArray[SomeNumber]) = assign(r,x)
  unaryMixedVar(T, `+=`, iadd)
  unaryMixedVar(T, `-=`, isub)
  unaryMixedVar(T, `*=`, imul)
  #unaryMixedVar(T, `/=`, idiv)
  template `/=`*(r: T; x:SomeNumber) = idiv(r,x.to(T))

  proc trace*(x:T):T {.inline,noInit.}= x
  proc norm2*(x:T):T {.inline,noInit.} = mul(x,x)
  proc norm2*(r:var T; x:T) {.inline.} = mul(r,x,x)
  proc inorm2*(r:var T; x:T) {.inline.} = imadd(r,x,x)
  proc max*(x,y:T):T {.inline,noInit.} = `P "_max_" S`(x,y)
  proc min*(x,y:T):T {.inline,noInit.} = `P "_min_" S`(x,y)
  proc abs*(x:T):T {.inline,noInit.} = max(x,neg(x))
  proc sqrt*(x:T):T {.inline,noInit.} = `P "_sqrt_" S`(x)
  proc rsqrt*(x:T):T {.inline,noInit.} = divd(sqrt(x),x)
  proc rsqrt*(r:var T; x:T) {.inline.} = r = rsqrt(x)
  map1(T,N, sin)
  map1(T,N, cos)
  map1(T,N, acos)
  map2(T,N, atan2)

basicDefs(m128,  float32,  4, mm, ps)
basicDefs(m128d, float64,  2, mm, pd)
basicDefs(m256,  float32,  8, mm256, ps)
basicDefs(m256d, float64,  4, mm256, pd)
basicDefs(m512,  float32, 16, mm512, ps)
basicDefs(m512d, float64,  8, mm512, pd)


proc simdReduce*(r:var SomeNumber; x:m128) {.inline.} =
  let y = mm_hadd_ps(x, x)
  let z = mm_hadd_ps(y, y)
  var t{.noInit.}:float32
  mm_store_ss(t.addr, z)
  r = (type(r))(t)
proc simdReduce*(r:var SomeNumber; x:m256) {.inline.} =
  let y = mm256_hadd_ps(x, mm256_permute2f128_ps(x, x, 1))
  let z = mm256_hadd_ps(y, y)
  let w = mm256_hadd_ps(z, z)
  r = (type(r))(w[0])
proc simdReduce*(r:var SomeNumber; x:m256d) {.inline.} =
  let y = mm256_hadd_pd(x, mm256_permute2f128_pd(x, x, 1))
  let z = mm256_hadd_pd(y, y)
  r = (type(r))(z[0])
proc simdReduce*(r:var SomeNumber; x:m512) {.inline.} =
  #r = (type(r))(mm512_reduce_add_ps(x))
  #let t = mm512_shuffle_f32x4(x, x, BASE4(1,0,3,2))
  #let t2 = add(x, t)
  r = x[0]
  for i in 1..<16:
    r += x[i]
proc simdReduce*(r:var SomeNumber; x:m512d) {.inline.} =
  #r = (type(r))(mm512_reduce_add_pd(x))
  r = x[0]
  for i in 1..<8:
    r += x[i]
proc simdReduce*(x:m128):float32 {.inline,noInit.} = simdReduce(result, x)
proc simdReduce*(x:m256):float32 {.inline,noInit.} = simdReduce(result, x)
proc simdReduce*(x:m256d):float64 {.inline,noInit.} = simdReduce(result, x)
proc simdReduce*(x:m512):float32 {.inline,noInit.} = simdReduce(result, x)
proc simdReduce*(x:m512d):float64 {.inline,noInit.} = simdReduce(result, x)

template simdSum*(x:m128):untyped = simdReduce(x)
template simdSum*(x:m256):untyped = simdReduce(x)
template simdSum*(x:m256d):untyped = simdReduce(x)
template simdSum*(x:m512):untyped = simdReduce(x)
template simdSum*(x:m512d):untyped = simdReduce(x)
template simdSum*(r:var SomeNumber; x:m128) = simdReduce(r, x)
template simdSum*(r:var SomeNumber; x:m256) = simdReduce(r, x)
template simdSum*(r:var SomeNumber; x:m256d) = simdReduce(r, x)
template simdSum*(r:var SomeNumber; x:m512) = simdReduce(r, x)
template simdSum*(r:var SomeNumber; x:m512d) = simdReduce(r, x)

proc simdMaxReduce*(r:var SomeNumber; x:m128) {.inline.} =
  r = x[0]
  for i in 1..<4:
    if r < x[i]: r = (type(r))(x[i])
proc simdMaxReduce*(r:var SomeNumber; x:m256) {.inline.} =
  r = x[0]
  for i in 1..<8:
    if r < x[i]: r = (type(r))(x[i])
proc simdMaxReduce*(r:var SomeNumber; x:m256d) {.inline.} =
  r = x[0]
  for i in 1..<4:
    if r < x[i]: r = (type(r))(x[i])
proc simdMaxReduce*(r:var SomeNumber; x:m512) {.inline.} =
  r = x[0]
  for i in 1..<16:
    if r < x[i]: r = (type(r))(x[i])
proc simdMaxReduce*(r:var SomeNumber; x:m512d) {.inline.} =
  r = x[0]
  for i in 1..<8:
    if r < x[i]: r = (type(r))(x[i])
proc simdMaxReduce*(x:m128):float32 {.inline,noInit.} = simdMaxReduce(result, x)
proc simdMaxReduce*(x:m256):float32 {.inline,noInit.} = simdMaxReduce(result, x)
proc simdMaxReduce*(x:m256d):float64 {.inline,noInit.} = simdMaxReduce(result, x)
proc simdMaxReduce*(x:m512):float32 {.inline,noInit.} = simdMaxReduce(result, x)
proc simdMaxReduce*(x:m512d):float64 {.inline,noInit.} = simdMaxReduce(result, x)

template simdMax*(x:m128):untyped = simdMaxReduce(x)
template simdMax*(x:m256):untyped = simdMaxReduce(x)
template simdMax*(x:m256d):untyped = simdMaxReduce(x)
template simdMax*(x:m512):untyped = simdMaxReduce(x)
template simdMax*(x:m512d):untyped = simdMaxReduce(x)
template simdMax*(r:var SomeNumber; x:m128) = simdMaxReduce(r, x)
template simdMax*(r:var SomeNumber; x:m256) = simdMaxReduce(r, x)
template simdMax*(r:var SomeNumber; x:m256d) = simdMaxReduce(r, x)
template simdMax*(r:var SomeNumber; x:m512) = simdMaxReduce(r, x)
template simdMax*(r:var SomeNumber; x:m512d) = simdMaxReduce(r, x)

# include perm, pack and blend
include simdX86Ops1

when defined(AVX):
  when defined(AVX512):
    proc toDouble*(x:SimdS8):SimdD8 {.inline,noInit.} =
      result = mm512_cvtps_pd(x)
  else:
    proc toDoubleA*(x:SimdS8):array[2,SimdD4] {.inline,noInit.} =
      result[0] = mm256_cvtps_pd(mm256_extractf128_ps(x,0))
      result[1] = mm256_cvtps_pd(mm256_extractf128_ps(x,1))
      #for i in 0..3: result[0][i] = x[i]
      #for i in 0..3: result[1][i] = x[4+i]

when defined(AVX512):
  proc toDoubleA*(x:SimdS16):array[2,SimdD8] {.inline,noInit.} =
    result[0] = mm512_cvtps_pd(mm512_castps512_ps256(x))
    var y{.noInit.}:SimdS16
    perm8(y, x)
    result[1] = mm512_cvtps_pd(mm512_castps512_ps256(y))

when defined(SimdS4):
  proc mm_cvtph_ps(x:m128i):m128
    {.importC:"_mm_cvtph_ps",header:"f16cintrin.h".}
  proc mm_cvtps_ph(x:m128,y:cint):m128i
    {.importC:"_mm_cvtps_ph",header:"f16cintrin.h".}
  template toHalf(x:SimdS4):SimdH4 = SimdH4(mm_cvtps_ph(x))
  template toSingle(x:SimdH4):SimdS4 = mm_cvtph_ps(x)
when defined(SimdS8):
  proc mm256_cvtph_ps(x:m128i):m256
    {.importC:"_mm256_cvtph_ps",header:"f16cintrin.h".}
  proc mm256_cvtps_ph(x:m256,y:cint):m128i
    {.importC:"_mm256_cvtps_ph",header:"f16cintrin.h".}
  template toHalf(x:SimdS8):SimdH8 = SimdH8(mm256_cvtps_ph(x,0))
  template toSingle(x:SimdH8):SimdS8 = mm256_cvtph_ps(m128i(x))
when defined(SimdS16):
  template toHalf(x:SimdS16):SimdH16 = SimdH16(mm512_cvtps_ph(x,0))
  template toSingle(x:SimdH16):SimdS16 = mm512_cvtph_ps(m256i(x))

# toSingle, toDouble, to(x,float32), to(x,float64)
discard """
template `lid`(x:untyped):untyped = to(x,`id`)
"""

when isMainModule:
  var x,y,z:m256d
  var d:float64
  var a = [1.0,2.0,3.0,4.0]
  assign(x, 1)
  assign(y, 2)
  assign(z, 0)

  echo z
  z = x+y
  echo z
  simdReduce(d, z)
  echo d
  assign(z, a)
  echo z
  simdReduce(d, z)
  echo d
  perm1(y, z)
  echo y
  echo z
  perm2(y, z)
  echo y

  assign(x, 1)
  echo x
  assign(x, 1, 2)
  echo x
  assign(x, 1, 2, 3)
  echo x
  assign(x, 1, 2, 3, 4)
  echo x
  assign(x, 1, 2, 3, 4, 5)
  echo x

  assign(x, a[0], a[1], a[2], a[3])
  echo x

  var s:m256
  assign(s, a[0], a[1], a[2], a[3])
  echo s

  var s8:SimdS8
  assign(s8, [0,1,2,3,4,5,6,7])
  var d8 = toDoubleA(s8)
  echo d8[0]
  echo d8[1]

  var h:SimdH8
  s = toSingle(h)
  h = toHalf(s)
  assign(s,[1,2,3,4,5,6,7,8])
  h = toHalf(s)
  s8 = toSingle(h)
  echo s8

  when declared(SimdS16):
    var s16:SimdS16
    assign(s16, [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])
    var h16 = toHalf(s16)
    var t16 = toSingle(h16)
    echo t16
