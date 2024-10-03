#import base/globals
#import simdGcc
#export simdGcc
import base/metaUtils
import math

import simd/simdWrap
export simdWrap

import base/basicOps
import base/stdUtils
import simd/simdArray
export simdArray

template msa(T,N,F: untyped) {.dirty,used.} =
  makeSimdArray(`T Obj`, N, F)
  type T* = Simd[`T Obj`]
  #template `T Array` = discard
  type `T Array`* = `T Obj`
  #static: echo "made type", $T

when defined(SSE) or defined(AVX) or defined(AVX512):
  import simd/simdX86
  export simdX86

#when VLEN>=16:
when true:
  when not declared(SimdS16):
    when declared(SimdS8):
      msa(SimdS16, 2, SimdS8[])
    elif declared(SimdS4):
      msa(SimdS16, 4, SimdS4[])
    elif declared(SimdS2):
      msa(SimdS16, 8, SimdS2[])
    else:
      msa(SimdS16, 16, float32)
  when not declared(SimdD16):
    when declared(SimdD8):
      msa(SimdD16, 2, `[]`SimdD8)
    elif declared(SimdD4):
      msa(SimdD16, 4, SimdD4[])
    elif declared(SimdD2):
      msa(SimdD16, 8, SimdD2[])
    else:
      msa(SimdD16, 16, float64)
  when not declared(SimdS16Obj):
    type SimdS16Obj* = `[]`(SimdS16)
  when not declared(SimdD16Obj):
    type SimdD16Obj* = `[]`(SimdD16)

#when VLEN>=8:
when true:
  when not declared(SimdS8):
    when declared(SimdS4):
      msa(SimdS8, 2, SimdS4[])
    elif declared(SimdS2):
      msa(SimdS8, 4, SimdS2[])
    else:
      msa(SimdS8, 8, float32)
  when not declared(SimdD8):
    when declared(SimdD4):
      msa(SimdD8, 2, SimdD4[])
    elif declared(SimdD2):
      msa(SimdD8, 4, SimdD2[])
    else:
      msa(SimdD8, 8, float64)
  when not declared(SimdS8Obj):
    type SimdS8Obj* = `[]`(SimdS8)
  when not declared(SimdD8Obj):
    type SimdD8Obj* = `[]`(SimdD8)

#when VLEN>=4:
when true:
  when not declared(SimdS4):
    when declared(SimdS2):
      msa(SimdS4, 2, SimdS2[])
    else:
      msa(SimdS4, 4, float32)
  when not declared(SimdD4):
    when declared(SimdD2):
      msa(SimdD4, 2, SimdD2[])
    else:
      msa(SimdD4, 4, float64)
  when not declared(SimdS4Obj):
    type SimdS4Obj* = `[]`(SimdS4)
  when not declared(SimdD4Obj):
    type SimdD4Obj* = `[]`(SimdD4)

#when VLEN>=2:
when true:
  when not declared(SimdS2):
    msa(SimdS2, 2, float32)
  when not declared(SimdD2):
    msa(SimdD2, 2, float64)
  when not declared(SimdD2Obj):
    type SimdD2Obj* = `[]`(SimdD2)

#when VLEN>=1:
when true:
  msa(SimdS1, 1, float32)
  msa(SimdD1, 1, float64)


## mixed precision assignment

when declared(SimdS1Array):
  when declared(SimdD1Array):
    template assign*(r: var SimdS1Obj, x: SimdD1Obj) = assign(r[], x[])
    template assign*(r: var SimdD1Obj, x: SimdS1Obj) = assign(r[], x[])
  else:
    template assign*(r: var SimdS1Obj, x: SimdD1Obj) = assign(r[], x)
    template assign*(r: var SimdD1Obj, x: SimdS1Obj) = assign(r, x[])
else:
  when declared(SimdD1Array):
    template assign*(r: var SimdS1Obj, x: SimdD1Obj) = assign(r, x[])
    template assign*(r: var SimdD1Obj, x: SimdS1Obj) = assign(r[], x)

when declared(SimdS2Array):
  when declared(SimdD2Array):
    template assign*(r: var SimdS2Obj, x: SimdD2Obj) = assign(r[], x[])
    template assign*(r: var SimdD2Obj, x: SimdS2Obj) = assign(r[], x[])
  else:
    template assign*(r: var SimdS2Obj, x: SimdD2Obj) = assign(r[], x)
    template assign*(r: var SimdD2Obj, x: SimdS2Obj) = assign(r, x[])
else:
  when declared(SimdD2Array):
    template assign*(r: var SimdS2Obj, x: SimdD2Obj) = assign(r, x[])
    template assign*(r: var SimdD2Obj, x: SimdS2Obj) = assign(r[], x)

when declared(SimdS4Array):
  when declared(SimdD4Array):
    template assign*(r: var SimdS4Obj, x: SimdD4Obj) = assign(r[], x[])
    template assign*(r: var SimdD4Obj, x: SimdS4Obj) = assign(r[], x[])
  else:
    template assign*(r: var SimdS4Obj, x: SimdD4Obj) = assign(r[], x)
    template assign*(r: var SimdD4Obj, x: SimdS4Obj) = assign(r, x[])
else:
  when declared(SimdD4Array):
    template assign*(r: var SimdS4Obj, x: SimdD4Obj) = assign(r, x[])
    template assign*(r: var SimdD4Obj, x: SimdS4Obj) = assign(r[], x)

when declared(SimdS8Array):
  when declared(SimdD8Array):
    template assign*(r: var SimdS8Obj, x: SimdD8Obj) = assign(r[], x[])
    template assign*(r: var SimdD8Obj, x: SimdS8Obj) = assign(r[], x[])
  else:
    template assign*(r: var SimdS8Obj, x: SimdD8Obj) = assign(r[], x)
    template assign*(r: var SimdD8Obj, x: SimdS8Obj) = assign(r, x[])
else:
  when declared(SimdD8Array):
    template assign*(r: var SimdS8Obj, x: SimdD8Obj) = assign(r, x[])
    template assign*(r: var SimdD8Obj, x: SimdS8Obj) = assign(r[], x)

when declared(SimdS16Array):
  when declared(SimdD16Array):
    template assign*(r: var SimdS16Obj, x: SimdD16Obj) = assign(r[], x[])
    template assign*(r: var SimdD16Obj, x: SimdS16Obj) = assign(r[], x[])
  else:
    template assign*(r: var SimdS16Obj, x: SimdD16Obj) = assign(r[], x)
    template assign*(r: var SimdD16Obj, x: SimdS16Obj) = assign(r, x[])
else:
  when declared(SimdD16Array):
    template assign*(r: var SimdS16Obj, x: SimdD16Obj) = assign(r, x[])
    template assign*(r: var SimdD16Obj, x: SimdS16Obj) = assign(r[], x)

## other generic helpers

template convert(x: typed, T: typedesc): untyped =
  var r {.noInit.}: T
  assign(r, x)
  r

template mapSimd*(t,f: untyped) {.dirty.} =
  proc f*(x: t): t {.inline,noInit.} =
    forStatic i, 0, x.numNumbers-1:
      result[i] = f(x[i])

template makeBinaryMixed(S,D,op) =
  template op*(x: S, y: D): untyped =
    op(toDouble(x),y)
  template op*(x: D, y: S): untyped =
    op(x,toDouble(y))

when declared(SimdS1):
  template eval*(x: SimdS1): untyped = x
  template toSingle*(x: typedesc[SimdS1Obj]): typedesc = SimdS1Obj
  template toDouble*(x: typedesc[SimdS1Obj]): typedesc = SimdD1Obj
  template toSingleImpl*(x: SimdS1Obj): untyped = x
  template toDoubleImpl*(x: SimdS1Obj): untyped = convert(x, SimdD1Obj)
  mapSimd(SimdS1, exp)
  mapSimd(SimdS1, ln)

when declared(SimdD1):
  template eval*(x: SimdD1): untyped = x
  template toSingle*(x: typedesc[SimdD1Obj]): typedesc = SimdS1Obj
  template toDouble*(x: typedesc[SimdD1Obj]): typedesc = SimdD1Obj
  template toSingleImpl*(x: SimdD1Obj): untyped = convert(x, SimdS1Obj)
  template toDoubleImpl*(x: SimdD1Obj): untyped = x
  mapSimd(SimdD1, exp)
  mapSimd(SimdD1, ln)

when declared(SimdS2):
  template eval*(x: SimdS2): untyped = x
  template toSingle*(x: typedesc[SimdS2Obj]): typedesc = SimdS2Obj
  template toDouble*(x: typedesc[SimdS2Obj]): typedesc = SimdD2Obj
  template toSingleImpl*(x: SimdS2Obj): untyped = x
  template toDoubleImpl*(x: SimdS2Obj): untyped = convert(x, SimdD2Obj)
  mapSimd(SimdS2, exp)
  mapSimd(SimdS2, ln)

when declared(SimdD2):
  template eval*(x: SimdD2): untyped = x
  template toSingle*(x: typedesc[SimdD2Obj]): typedesc = SimdS2Obj
  template toDouble*(x: typedesc[SimdD2Obj]): typedesc = SimdD2Obj
  template toSingleImpl*(x: SimdD2Obj): untyped = convert(x, SimdS2Obj)
  template toDoubleImpl*(x: SimdD2Obj): untyped = x
  mapSimd(SimdD2, exp)
  mapSimd(SimdD2, ln)

when declared(SimdS4):
  template eval*(x: SimdS4): untyped = x
  template toSingleImpl*(x: SimdS4Obj): untyped = x
  template toDoubleImpl*(x: SimdS4Obj): untyped = convert(x, SimdD4Obj)
  mapSimd(SimdS4, exp)
  mapSimd(SimdS4, ln)

when declared(SimdD4):
  #template assign*(r: array[4,float32], x: SimdD4): untyped =
  #  assign(r, toSingle(x))
  template eval*(x: SimdD4): untyped = x
  template toSingleImpl*(x: SimdD4Obj): untyped = convert(x, SimdS4Obj)
  template toDoubleImpl*(x: SimdD4Obj): untyped = x
  #template min*(x: SomeNumber, y: SimdD4): untyped = min(x.to(SimdD4), y)
  #template max*(x: SomeNumber, y: SimdD4): untyped = max(x.to(SimdD4), y)
  mapSimd(SimdD4, exp)
  mapSimd(SimdD4, ln)

when declared(SimdS4) and declared(SimdD4):
  #proc toSingle*(x: SimdD4): SimdS4 {.inline,noInit.} =
  #  for i in 0..<4:
  #    result[i] = x[i]
  #proc toDouble*(x: SimdS4): SimdD4 {.inline,noInit.} =
  #  for i in 0..<4:
  #    result[i] = x[i]
  #template assign*(r: SimdS4, x: SimdD4): untyped =
  #  r = toSingle(x)
  #template assign*(r: SimdD4, x: SimdS4): untyped =
  #  r = toDouble(x)
  #converter promote*(x: SimdS4): SimdD4 {.inline,noInit.} =
  #  assign(result, x)
  #template assign*(x: SimdS4; y: SimdD4): untyped =
  #  assign(x, toSingle(y))
  makeBinaryMixed(SimdS4, SimdD4, `+`)
  makeBinaryMixed(SimdS4, SimdD4, `-`)
  makeBinaryMixed(SimdS4, SimdD4, `*`)
  proc inorm2*(r:var SimdD4; x:SimdS4) {.inline.} =
    let y = toDouble(x)
    inorm2(r, y)

when declared(SimdS8):
  template toSingleImpl*(x: SimdS8Obj): untyped = x
  template toDoubleImpl*(x: SimdS8Obj): untyped = convert(x, SimdD8Obj)
  mapSimd(SimdS8, exp)
  mapSimd(SimdS8, ln)

when declared(SimdD8):
  template eval*(x: SimdD8): untyped = x
  template toSingleImpl*(x: SimdD8Obj): untyped = convert(x, SimdS8Obj)
  template toDoubleImpl*(x: SimdD8Obj): untyped = x
  mapSimd(SimdD8, exp)
  mapSimd(SimdD8, ln)

when declared(SimdD8) and declared(SimdS8):
  #template toDouble*(x: SimdD8): untyped = x
  #proc toSingle*(x: SimdD8): SimdS8 {.inline,noInit.} =
  #  for i in 0..<8:
  #    result[i] = x[i]
  #template toSingle*(x: typedesc[SimdD8]): untyped = SimdS8
  #proc toDouble*(x: SimdS8): SimdD8 {.inline,noInit.} =
  #  for i in 0..<8:
  #    result[i] = x[i]
  #template assign*(r: SimdS8, x: SimdD8): untyped =
  #  r := toSingle(x)
  template `:=`*(r: SimdS8, x: SimdD8): untyped =
    assign(r, x)
  #template assign*(r: SimdD8, x: SimdS8): untyped =
  #  r = toDouble(x)
  #converter promote*(x: SimdS8): SimdD8 {.inline,noInit.} =
  #  assign(result, x)
  template isub*(r: SimdD8, x: SimdS8): untyped =
    isub(r, toDouble(x))
  template imadd*(r: SimdD8, x: SimdD8, y: SimdS8): untyped =
    imadd(r, x, toDouble(y))
  template imsub*(r: SimdD8, x: SimdD8, y: SimdS8): untyped =
    imsub(r, x, toDouble(y))
  #template `-=`*(r: SimdD8, x: SimdD8): untyped = isub(r, x)


when declared(SimdS16):
  template toSingleImpl*(x: SimdS16Obj): untyped = x
  template toDoubleImpl*(x: SimdS16Obj): untyped = convert(x, SimdD16Obj)
  mapSimd(SimdS16, exp)
  mapSimd(SimdS16, ln)

when declared(SimdD16):
  template eval*(x: SimdD16): auto = x
  template toSingleImpl*(x: SimdD16Obj): untyped = convert(x, SimdS16Obj)
  template toDoubleImpl*(x: SimdD16Obj): untyped = x
  mapSimd(SimdD16, exp)
  mapSimd(SimdD16, ln)

template assignX*(x: var Simd, y: SomeNumber) =
  static: echo "assignX Simd SomeNumber"
  debugType: x
  debugType: y
  assign(x[], y)

template assignX*(x: var Simd, y: Simd2) =
  static: echo "assignX Simd Simd"
  debugType: x
  debugType: y
  assign(x[], y[])


