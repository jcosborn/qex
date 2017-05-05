import macros
#import ../metaUtils
#import ../basicOps
import base
import simdX86Types
export simdX86Types

import simdX86Ops
export simdX86Ops

import simdArray

template tryArray(T,L,B:untyped):untyped =
  when (not declared(T)) and declared(B):
    makeSimdArray(T, L, B)
macro makeArray(P,N:untyped):auto =
  let n = N.intVal
  let t = ident("Simd" & $P & $n)
  var m = n div 2
  result = newStmtList()
  while m>0:
    let b = ident("Simd" & $P & $m)
    let l = n div m
    result.add getAst(tryArray(t,newLit(l),b))
    m = m div 2
  #echo result.repr

makeArray(D, 16)
makeArray(D,  8)
makeArray(D,  4)

#when defined(SSE):
#proc toDoubleA*(x:SimdS4):array[2,SimdD2] {.inline,noInit.} =
#  result[0] = mm_cvtps_pd(x)
#  var y{.noInit.}:SimdS4
#  perm2(y, x)
#  result[1] = mm_cvtps_pd(y)

when defined(AVX):
  when not defined(AVX512):
    proc toDouble*(x:SimdS8):SimdD8 {.inline,noInit.} =
      #result = SimdD8(toDoubleA(x))
      result := toDoubleA(x)

#when declared(SimdS4):
#  proc toDouble*(x:SimdS4):SimdD4 {.inline,noInit.} =
#    result = SimdD4(toDoubleA(x))
#  proc inorm2*(r:var SimdD4; x:SimdS4) {.inline.} =
#    let y = toDouble(x)
#    inorm2(r, y)

when declared(SimdS8):
  #proc toDouble*(x:SimdS8):SimdD8 {.inline,noInit.} =
  #  result = SimdD8(toDoubleA(x))
  proc inorm2*(r:var SimdD8; x:SimdS8) {.inline.} =
    var xx{.noInit.} = toDouble(x)
    inorm2(r, xx)
  proc imadd*(r:var SimdD8; x,y:SimdS8) {.inline.} =
    var xx{.noInit.} = toDouble(x)
    var yy{.noInit.} = toDouble(y)
    imadd(r, xx, yy)
  proc imsub*(r:var SimdD8; x,y:SimdS8) {.inline.} =
    let xd = toDouble(x)
    let yd = toDouble(y)
    imsub(r, xd, yd)

when declared(SimdS16):
  proc toDouble*(x:SimdS16):SimdD16 {.inline,noInit.} =
    #for i in 0..15: result[i] = float64(x[i])
    result = SimdD16(v: toDoubleA(x))
  proc inorm2*(r:var SimdD16; x:SimdS16) {.inline.} = inorm2(r, toDouble(x))
  proc imadd*(r:var SimdD16; x,y:SimdS16) {.inline.} =
    var xx{.noInit.} = toDouble(x)
    var yy{.noInit.} = toDouble(y)
    imadd(r, xx, yy)
  proc imsub*(r:var SimdD16; x,y:SimdS16) {.inline.} =
    let xd = toDouble(x)
    let yd = toDouble(y)
    imsub(r, xd, yd)

when declared(SimdD4):
  template toDouble*(x:SimdD4):untyped = x
when declared(SimdD8):
  template toDouble*(x:SimdD8):untyped = x
when declared(SimdD16):
  template toDouble*(x:SimdD16):untyped = x

when isMainModule:
  var s8:SimdS8
  assign(s8, [0,1,2,3,4,5,6,7])
  var d8 = toDouble(s8)
  echo d8
  inorm2(d8, s8)
  echo d8
