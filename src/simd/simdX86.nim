import macros
import ../metaUtils
import ../basicOps
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
