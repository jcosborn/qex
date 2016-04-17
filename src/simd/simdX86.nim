import ../metaUtils
import ../basicOps
import simdX86Types
export simdX86Types

import simdX86Ops
export simdX86Ops

# move to simdArray
type SimdD8* = distinct array[2,SimdD4]
template `[]`*(x:SimdD8):expr = (array[2,SimdD4])(x)
template `[]`*(x:SimdD8; i:SomeInteger):expr = x[][i div 4][i mod 4]
template `[]=`*(x:SimdD8; i:SomeInteger; y:any) = x[][i div 4][i mod 4] = y
template numberType*(x:typedesc[SimdD8]):typedesc = float64
template numberType*(x:SimdD8):typedesc = float64
template numNumbers*(x:typedesc[SimdD8]):untyped = 8
template numNumbers*(x:SimdD8):untyped = 8
template simdType*(x:typedesc[SimdD8]):typedesc = SimdD8
template simdType*(x:SimdD8):typedesc = SimdD8
template simdLength*(x:SimdD8):untyped = 8
template simdLength*(x:typedesc[SimdD8]):untyped = 8

proc to*(x:SomeNumber; y:typedesc[SimdD8]):SimdD8 {.inline,noInit.} =
  #{.emit:"#define memset(a,b,c)".}
  assign(result[][0], x)
  assign(result[][1], x)

proc simdReduce*(x:SimdD8):float64 {.inline,noInit.} =
  result = (type(result))(simdReduce(x[][0])+simdReduce(x[][1]))
proc simdReduce*(r:var SomeNumber; x:SimdD8) {.inline.} =
  r = (type(r))(simdReduce(x[][0])+simdReduce(x[][1]))
template simdSum*(x:SimdD8):expr = simdReduce(x)
template simdSum*(r:var SomeNumber; x:SimdD8) = simdReduce(r, x)
proc `-`*(x:SimdD8):SimdD8 {.inline,noInit.} =
  neg(result[][0], x[][0])
  neg(result[][1], x[][1])

template map011(T,n,op1,op2:untyped):untyped =
  proc op1*(x:T):T {.inline,noInit.} =
    result[][0] = op2(x[][0])
    result[][1] = op2(x[][1])
template map021(T,n,op1,op2:untyped):untyped =
  proc op1*(x,y:T):T {.inline,noInit.} =
    result[][0] = op2(x[][0], y[][0])
    result[][1] = op2(x[][1], y[][1])
template map110(T,n,op1,op2:untyped):untyped =
  proc op1*(r:var T; x:T) {.inline.} =
    op2(r[][0], x[][0])
    op2(r[][1], x[][1])
template map120(T,n,op1,op2:untyped):untyped =
  proc op1*(r:var T; x,y:T) {.inline.} =
    op2(r[][0], x[][0], y[][0])
    op2(r[][1], x[][1], y[][1])
template map130(T,n,op1,op2:untyped):untyped =
  proc op1*(r:var T; x,y,z:T) {.inline.} =
    op2(r[][0], x[][0], y[][0], z[][0])
    op2(r[][1], x[][1], y[][1], z[][1])

map011(SimdD8, 2, abs, abs)
map011(SimdD8, 2, trace, trace)
map011(SimdD8, 2, norm2, norm2)
map011(SimdD8, 2, sqrt, sqrt)
map011(SimdD8, 2, sin, sin)
map011(SimdD8, 2, cos, cos)
map011(SimdD8, 2, acos, acos)

map021(SimdD8, 2, add, add)
map021(SimdD8, 2, sub, sub)
map021(SimdD8, 2, mul, mul)
map021(SimdD8, 2, divd, divd)
map021(SimdD8, 2, `+`, add)
map021(SimdD8, 2, `-`, sub)
map021(SimdD8, 2, `*`, mul)
map021(SimdD8, 2, `/`, divd)

map110(SimdD8, 2, assign, assign)
map110(SimdD8, 2, neg, neg)
map110(SimdD8, 2, iadd, iadd)
map110(SimdD8, 2, isub, isub)
map110(SimdD8, 2, imul, imul)
map110(SimdD8, 2, norm2, norm2)
map110(SimdD8, 2, inorm2, inorm2)
map110(SimdD8, 2, rsqrt, rsqrt)

map120(SimdD8, 2, add, add)
map120(SimdD8, 2, sub, sub)
map120(SimdD8, 2, mul, mul)
map120(SimdD8, 2, divd, divd)
map120(SimdD8, 2, imadd, imadd)
map120(SimdD8, 2, imsub, imsub)

map130(SimdD8, 2, msub, msub)

#proc assign*(r:var SimdD8; x:SomeNumber) {.inline,neverInit.} =
proc assign*(r:var SimdD8; x:SomeNumber) {.inline.} =
  #{.emit:"#define memset(a,b,c)".}
  assign(r, x.to(SimdD8))
  #assign(r[][0], x)
  #assign(r[][1], r[][0])
proc assign*(r:var SimdD8; x:array[8,SomeNumber]) {.inline.} =
  assign(r[][0], cast[ptr array[4,type(x[0])]](unsafeAddr(x[0]))[])
  assign(r[][1], cast[ptr array[4,type(x[4])]](unsafeAddr(x[4]))[])
template add*(r:var SimdD8; x:SomeNumber; y:SimdD8) = add(r, x.to(SimdD8), y)
template sub*(r:var SimdD8; x:SomeNumber; y:SimdD8) = sub(r, x.to(SimdD8), y)
template mul*(r:var SimdD8; x:SomeNumber; y:SimdD8) = mul(r, x.to(SimdD8), y)
template mul*(r:var SimdD8; x:SimdD8; y:SomeNumber) = mul(r, x, y.to(SimdD8))
template imsub*(r:var SimdD8; x:SomeNumber; y:SimdD8) = imsub(r, x.to(SimdD8), y)
template divd*(r:var SimdD8; x:SomeNumber; y:SimdD8) = divd(r, x.to(SimdD8), y)
template imul*(r:var SimdD8; x:SomeNumber) = imul(r, x.to(SimdD8))
template msub*(r:var SimdD8; x:SomeNumber; y,z:SimdD8) = msub(r, x.to(SimdD8), y, z)
template `:=`*(r:var SimdD8; x:array[8,SomeNumber]) = assign(r, x)
template `+`*(x:SomeNumber; y:SimdD8):SimdD8 = add(x.to(SimdD8), y)
template `-`*(x:SomeNumber; y:SimdD8):SimdD8 = sub(x.to(SimdD8), y)
template `*`*(x:SomeNumber; y:SimdD8):SimdD8 = mul(x.to(SimdD8), y)
template `*`*(x:SimdD8; y:SomeNumber):SimdD8 = mul(x, y.to(SimdD8))
template `/`*(x:SomeNumber; y:SimdD8):SimdD8 = divd(x.to(SimdD8), y)
template `/`*(x:SimdD8; y:SomeNumber):SimdD8 = divd(x, y.to(SimdD8))

proc `$`*(x:SimdD8):string =
  result = "[" & $x[0]
  for i in 1..<8:
    result &= "," & $x[i]
  result &= "]"

map110(SimdD8, 2, perm1, perm1)
map110(SimdD8, 2, perm2, perm2)
proc perm4*(r:var SimdD8; x:SimdD8) {.inline.} =
  assign(r[][0], x[][1])
  assign(r[][1], x[][0])

template pck(op:untyped) =
  proc op*(r:var openArray[SomeNumber], x:SimdD8,
           l:var openarray[SomeNumber]) {.inline.} =
    let ra = cast[ptr array[2,type(r[0])]](r[2].addr)
    let la = cast[ptr array[2,type(l[0])]](l[2].addr)
    op(r, x[][0], l)
    op(ra[], x[][1], la[])
pck(packp1)
pck(packm1)
pck(packp2)
pck(packm2)
proc packp4*(r:var openArray[SomeNumber], x:SimdD8,
             l:var openarray[SomeNumber]) {.inline.} =
  assign(cast[ptr float64](l.addr), x[][0])
  assign(cast[ptr float64](r.addr), x[][1])
proc packm4*(r:var openArray[SomeNumber], x:SimdD8,
             l:var openarray[SomeNumber]) {.inline.} =
  assign(cast[ptr float64](r.addr), x[][0])
  assign(cast[ptr float64](l.addr), x[][1])

template blnd(op:untyped) =
  proc op*(x:var SimdD8; r:openArray[SomeNumber];
           l:openArray[SomeNumber]) {.inline.} =
    let ra = cast[ptr array[2,type(r[0])]](unsafeAddr(r[2]))
    let la = cast[ptr array[2,type(l[0])]](unsafeAddr(l[2]))
    op(x[][0], r, l)
    op(x[][1], ra[], la[])
blnd(blendp1)
blnd(blendm1)
blnd(blendp2)
blnd(blendm2)
proc blendp4*(x:var SimdD8, r:openArray[SomeNumber],
             l:openarray[SomeNumber]) {.inline.} =
  assign(x[][0], cast[ptr float64](unsafeAddr(l)))
  assign(x[][1], cast[ptr float64](unsafeAddr(r)))
proc blendm4*(x:var SimdD8, r:openArray[SomeNumber],
             l:openarray[SomeNumber]) {.inline.} =
  assign(x[][0], cast[ptr float64](unsafeAddr(r)))
  assign(x[][1], cast[ptr float64](unsafeAddr(l)))
