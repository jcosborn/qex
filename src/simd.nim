#import simdGcc
#export simdGcc
import base/metaUtils
import math

when defined(SSE) or defined(AVX) or defined(AVX512):
  import simd/simdX86
  export simdX86
elif defined(QPX):
  import simd/simdQpx
  export simdQpx
else:
  import base/globals
  import base/basicOps
  import simd/simdArray
  export simdArray
  when VLEN==1:
    makeSimdArray(SimdS1, 1, float32)
    makeSimdArray(SimdD1, 1, float64)
  when VLEN==2:
    makeSimdArray(SimdS2, 2, float32)
    makeSimdArray(SimdD2, 2, float64)
  when VLEN==4:
    makeSimdArray(SimdS4, 4, float32)
    makeSimdArray(SimdD4, 4, float64)
  when VLEN==8:
    makeSimdArray(SimdS8, 8, float32)
    makeSimdArray(SimdD8, 8, float64)
    proc toDouble*(x: SimdS8): SimdD8 {.inline,noInit.} =
      for i in 0..<8:
        result[i] = x[i]

#import simd/simdGeneric
#export simdGeneric

template mapSimd*(t,f: untyped) {.dirty.} =
  proc f*(x: t): t {.inline,noInit.} =
    forStatic i, 0, x.numNumbers-1:
      result[i] = f(x[i])

when declared(SimdD1):
  template eval*(x: SimdD1): untyped = x
  template toSingleImpl*(x: SimdD1): untyped = toSingle(x)
  template toDoubleImpl*(x: SimdD1): untyped = x
  mapSimd(SimdD1, exp)

when declared(SimdD2):
  template eval*(x: SimdD2): untyped = x
  template toSingleImpl*(x: SimdD2): untyped = toSingle(x)
  template toDoubleImpl*(x: SimdD2): untyped = x

when declared(SimdS4):
  template eval*(x: SimdS4): untyped = x
  template toSingleImpl*(x: SimdS4): untyped = x
  template toDoubleImpl*(x: SimdS4): untyped = toDouble(x)
  template inv*(x: SimdS4): untyped = 1.0/x

when declared(SimdD4):
  template assign*(r: array[4,float32], x: SimdD4): untyped =
    assign(r, toSingle(x))
  template eval*(x: SimdD4): untyped = x
  template toSingleImpl*(x: SimdD4): untyped = toSingle(x)
  template toDoubleImpl*(x: SimdD4): untyped = x
  mapSimd(SimdD4, exp)
  template inv*(x: SimdD4): untyped = 1.0/x

when declared(SimdS4) and declared(SimdD4):
  proc toSingle*(x: SimdD4): SimdS4 {.inline,noInit.} =
    for i in 0..<4:
      result[i] = x[i]
  proc toDouble*(x: SimdS4): SimdD4 {.inline,noInit.} =
    for i in 0..<4:
      result[i] = x[i]
  template assign*(r: SimdS4, x: SimdD4): untyped =
    r = toSingle(x)
  template assign*(r: SimdD4, x: SimdS4): untyped =
    r = toDouble(x)
  converter promote*(x: SimdS4): SimdD4 {.inline,noInit.} =
    assign(result, x)
  proc inorm2*(r:var SimdD4; x:SimdS4) {.inline.} =
    let y = toDouble(x)
    inorm2(r, y)

when declared(SimdS8):
  template toSingleImpl*(x: SimdS8): untyped = x
  template toSingleImpl*(x: SimdD8): untyped = toSingle(x)
  template toDoubleImpl*(x: SimdS8): untyped = toDouble(x)
  template toDoubleImpl*(x: SimdD8): untyped = x
  mapSimd(SimdS8, exp)
  template inv*(x: SimdS8): untyped = 1.0/x

when declared(SimdD8):
  template eval*(x: SimdD8): untyped = x
  mapSimd(SimdD8, exp)
  template inv*(x: SimdD8): untyped = 1.0/x

when declared(SimdD8) and declared(SimdS8):
  #template toDouble*(x: SimdD8): untyped = x
  proc toSingle*(x: SimdD8): SimdS8 {.inline,noInit.} =
    for i in 0..<8:
      result[i] = x[i]
  #proc toDouble*(x: SimdS8): SimdD8 {.inline,noInit.} =
  #  for i in 0..<8:
  #    result[i] = x[i]
  template assign*(r: SimdS8, x: SimdD8): untyped =
    r := toSingle(x)
  template `:=`*(r: SimdS8, x: SimdD8): untyped =
    assign(r, x)
  template assign*(r: SimdD8, x: SimdS8): untyped =
    r = toDouble(x)
  converter promote*(x: SimdS8): SimdD8 =
    assign(result, x)
  template isub*(r: SimdD8, x: SimdS8): untyped =
    isub(r, toDouble(x))
  template imadd*(r: SimdD8, x: SimdD8, y: SimdS8): untyped =
    imadd(r, x, toDouble(y))
  template imsub*(r: SimdD8, x: SimdD8, y: SimdS8): untyped =
    imsub(r, x, toDouble(y))
  template `-=`*(r: SimdD8, x: SimdD8): untyped = isub(r, x)


when declared(SimdS16):
  proc toSingle*(x: SimdD16): SimdS16 {.inline,noInit.} =
    for i in 0..<16:
      result[i] = x[i]
  template assign*(r: SimdS16, x: SimdD16): untyped =
    r = toSingle(x)
  template assign*(r: SimdD16, x: SimdS16): untyped =
    r = toDouble(x)
  converter promote*(x: SimdS16): SimdD16 =
    assign(result, x)
  template toSingleImpl*(x: SimdS16): untyped = x
  template toSingleImpl*(x: SimdD16): untyped = toSingle(x)
  template toDoubleImpl*(x: SimdS16): untyped = toDouble(x)
  template toDoubleImpl*(x: SimdD16): untyped = x

when declared(SimdD16):
  template eval*(x: SimdD16): auto = x

