#import simdGcc
#export simdGcc

when defined(QPX):
  import simd/simdQpx
  export simdQpx
else:
  import simd/simdX86
  export simdX86


proc toSingle*(x: SimdD4): SimdS4 {.inline,noInit.} =
  for i in 0..<4:
    result[i] = x[i]
proc toSingle*(x: SimdD8): SimdS8 {.inline,noInit.} =
  for i in 0..<8:
    result[i] = x[i]

when declared(SimdS4):
  proc toDouble*(x: SimdS4): SimdD4 {.inline,noInit.} =
    for i in 0..<4:
      result[i] = x[i]
  proc inorm2*(r:var SimdD4; x:SimdS4) {.inline.} =
    let y = toDouble(x)
    inorm2(r, y)

template assign*(r: array[4,float32], x: SimdD4): untyped =
  assign(r, toSingle(x))
template assign*(r: SimdS4, x: SimdD4): untyped =
  r = toSingle(x)
template assign*(r: SimdS8, x: SimdD8): untyped =
  r = toSingle(x)

template isub*(r: SimdD8, x: SimdS8): untyped =
  isub(r, toDouble(x))
template imadd*(r: SimdD8, x: SimdD8, y: SimdS8): untyped =
  imadd(r, x, toDouble(y))
template imsub*(r: SimdD8, x: SimdD8, y: SimdS8): untyped =
  imsub(r, x, toDouble(y))
