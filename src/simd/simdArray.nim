import macros
import base
#import ../metaUtils
#import ../basicOps

template map011(T,L,op1,op2:untyped):untyped {.dirty.} =
  proc op1*(x:T):T {.inline,noInit.} =
    forStatic i, 0, L-1:
      result[][i] = op2(x[][i])
template map021(T,L,op1,op2:untyped):untyped {.dirty.} =
  proc op1*(x,y:T):T {.inline,noInit.} =
    forStatic i, 0, L-1:
      result[][i] = op2(x[][i], y[][i])
template map110(T,L,op1,op2:untyped):untyped {.dirty.} =
  proc op1*(r:var T; x:T) {.inline.} =
    forStatic i, 0, L-1:
      op2(r[][i], x[][i])
template map120(T,L,op1,op2:untyped):untyped {.dirty.} =
  proc op1*(r:var T; x,y:T) {.inline.} =
    forStatic i, 0, L-1:
      op2(r[][i], x[][i], y[][i])
template map130(T,L,op1,op2:untyped):untyped {.dirty.} =
  proc op1*(r:var T; x,y,z:T) {.inline.} =
    forStatic i, 0, L-1:
      op2(r[][i], x[][i], y[][i], z[][i])
template makePermX(F,P,T,L,N0) {.dirty.} =
  when N0>P:
    bind map110
    map110(T, L, F, F)
  else:
    proc F*(r:var T; x:T) {.inline.} =
      const b = (P div N0) and (L-1)
      forStatic i, 0, L-1:
        assign(r[][i], x[][i xor b])
template makePerm(P,T,L,N0) {.dirty.} =
  bind makePermX
  evalBacktic:
    makePermX(`"perm" P`,P,T,L,N0)
template makePackPX(F,P,T,L,N0) {.dirty.} =
  when N0>P:
    proc F*(r:var openArray[SomeNumber], x:T,
            l:var openarray[SomeNumber]) {.inline.} =
      const N02 = N0 div 2
      let ra = cast[ptr array[L,array[N02,type(r[0])]]](r[0].addr)
      let la = cast[ptr array[L,array[N02,type(l[0])]]](l[0].addr)
      forStatic i, 0, L-1:
        F(ra[][i], x[][i], la[][i])
  else:
    proc F*(r:var openArray[SomeNumber], x:T,
            l:var openarray[SomeNumber]) {.inline.} =
      const L2 = L div 2
      let ra = cast[ptr array[L2,array[N0,type(r[0])]]](r[0].addr)
      let la = cast[ptr array[L2,array[N0,type(l[0])]]](l[0].addr)
      const b = (P div N0) and (L-1)
      var ir,il = 0
      forStatic i, 0, L-1:
        if (i and b) == 0:
          assign(la[][il], x[][i])
          inc il
        else:
          assign(ra[][ir], x[][i])
          inc ir
template makePackMX(F,P,T,L,N0) {.dirty.} =
  # P: perm, T: vector type, L: outer vec len, N0: inner vec len
  when N0>P:
    proc F*(r:var openArray[SomeNumber], x:T,
            l:var openarray[SomeNumber]) {.inline.} =
      const N02 = N0 div 2
      let ra = cast[ptr array[L,array[N02,type(r[0])]]](r[0].addr)
      let la = cast[ptr array[L,array[N02,type(l[0])]]](l[0].addr)
      forStatic i, 0, L-1:
        F(ra[][i], x[][i], la[][i])
  else:
    proc F*(r:var openArray[SomeNumber], x:T,
            l:var openarray[SomeNumber]) {.inline.} =
      const L2 = L div 2
      let ra = cast[ptr array[L2,array[N0,type(r[0])]]](r[0].addr)
      let la = cast[ptr array[L2,array[N0,type(l[0])]]](l[0].addr)
      const b = (P div N0) and (L-1)
      var ir,il = 0
      forStatic i, 0, L-1:
        if (i and b) == 0:
          assign(ra[][ir], x[][i])
          inc ir
        else:
          assign(la[][il], x[][i])
          inc il
template makePackP(P,T,L,N0) {.dirty.} =
  bind makePackPX
  evalBacktic:
    makePackPX(`"packp" P`,P,T,L,N0)
template makePackM(P,T,L,N0) {.dirty.} =
  bind makePackMX
  evalBacktic:
    makePackMX(`"packm" P`,P,T,L,N0)
template makeBlendPX(F,P,T,L,N0) {.dirty.} =
  when N0>P:
    proc F*(x:var T; r:openArray[SomeNumber];
            l:openArray[SomeNumber]) {.inline.} =
      const N02 = N0 div 2
      let ra = cast[ptr array[L,array[N02,type(r[0])]]](unsafeAddr(r[0]))
      let la = cast[ptr array[L,array[N02,type(l[0])]]](unsafeAddr(l[0]))
      forStatic i, 0, L-1:
        F(x[][i], ra[][i], la[][i])
  else:
    proc F*(x:var T; r:openArray[SomeNumber];
            l:openArray[SomeNumber]) {.inline.} =
      const L2 = L div 2
      let ra = cast[ptr array[L2,array[N0,type(r[0])]]](unsafeAddr(r[0]))
      let la = cast[ptr array[L2,array[N0,type(l[0])]]](unsafeAddr(l[0]))
      const b = (P div N0) and (L-1)
      var ir,il = 0
      forStatic i, 0, L-1:
        if (i and b) == 0:
          assign(x[][i], la[][il])
          inc il
        else:
          assign(x[][i], ra[][ir])
          inc ir
template makeBlendMX(F,P,T,L,N0) {.dirty.} =
  when N0>P:
    proc F*(x:var T; r:openArray[SomeNumber];
            l:openArray[SomeNumber]) {.inline.} =
      const N02 = N0 div 2
      let ra = cast[ptr array[L,array[N02,type(r[0])]]](unsafeAddr(r[0]))
      let la = cast[ptr array[L,array[N02,type(l[0])]]](unsafeAddr(l[0]))
      forStatic i, 0, L-1:
        F(x[][i], ra[][i], la[][i])
  else:
    proc F*(x:var T; r:openArray[SomeNumber];
            l:openArray[SomeNumber]) {.inline.} =
      const L2 = L div 2
      let ra = cast[ptr array[L2,array[N0,type(r[0])]]](unsafeAddr(r[0]))
      let la = cast[ptr array[L2,array[N0,type(l[0])]]](unsafeAddr(l[0]))
      const b = (P div N0) and (L-1)
      var ir,il = 0
      forStatic i, 0, L-1:
        if (i and b) == 0:
          assign(x[][i], ra[][ir])
          inc ir
        else:
          assign(x[][i], la[][il])
          inc il
template makeBlendP(P,T,L,N0) {.dirty.} =
  bind makeBlendPX
  evalBacktic:
    makeBlendPX(`"blendp" P`,P,T,L,N0)
template makeBlendM(P,T,L,N0) {.dirty.} =
  bind makeBlendMX
  evalBacktic:
    makeBlendMX(`"blendm" P`,P,T,L,N0)

discard """
  proc op*(x:var T; r:openArray[SomeNumber];
           l:openArray[SomeNumber]) {.inline.} =
    let ra = cast[ptr array[2,type(r[0])]](unsafeAddr(r[2]))
    let la = cast[ptr array[2,type(l[0])]](unsafeAddr(l[2]))
    op(x[][0], r, l)
    op(x[][1], ra[], la[])
"""


# T = Simd{S,D}{L} = array[L,B]
# B ~ array[N0,F]
template makeSimdArray*(T:untyped;L,B:typed):untyped {.dirty.} =
  makeSimdArray2(T, L, B, numberType(B), numNumbers(B), L*numNumbers(B))
template makeSimdArray2*(T:untyped;L,B,F,N0,N:typed):untyped {.dirty.} =
  bind map011, map021, map110, map120, map130
  bind makePerm, makePackP, makePackM, makeBlendP, makeBlendM
  #type T* = distinct array[L,B]
  type T* = object
    v*: array[L,B]
  template numberType*(x:typedesc[T]):typedesc = F
  template numberType*(x:T):typedesc = F
  template numNumbers*(x:typedesc[T]):untyped = N
  template numNumbers*(x:T):untyped = N
  template simdType*(x:typedesc[T]):typedesc = T
  template simdType*(x:T):typedesc = T
  template simdLength*(x:T):untyped = N
  template simdLength*(x:typedesc[T]):untyped = N
  #template `[]`*(x:T):untyped = (array[L,B])(x)
  template `[]`*(x:T):untyped = x.v
  template `[]`*(x:T; i:SomeInteger):untyped = x[][i div N0][i mod N0]
  template `[]=`*(x:T; i:SomeInteger; y:any) = x[][i div N0][i mod N0] = y
  template load1*(x:T):untyped = x
  proc to*(x:SomeNumber; y:typedesc[T]):T {.inline,noInit.} =
    subst(i,_):
      forStatic i, 0, L-1:
        assign(result[][i], x)
  proc simdReduce*(r: var SomeNumber; x: T) {.inline.} =
    var y = add(x[][0], x[][1])
    forStatic i, 2, L-1:
      iadd(y, x[][i])
    r = (type(r))(simdReduce(y))
  proc simdReduce*(x:T):F {.noInit,inline.} = simdReduce(result, x)
  template simdSum*(r:var SomeNumber; x:T) = simdReduce(r, x)
  template simdSum*(x:T):untyped = simdReduce(x)
  proc simdMaxReduce*(r:var SomeNumber; x:T) {.inline.} =
    mixin simdMaxReduce
    var y = x[][0]
    subst(i,_):
      forStatic i, 1, L-1:
        let c = x[][i]
        y = max(y, c)
    r = (type(r))(simdMaxReduce(y))
  proc simdMaxReduce*(x:T):F {.noinit,inline.} = simdMaxReduce(result, x)
  template simdMax*(r:var SomeNumber; x:T) = simdMaxReduce(r, x)
  template simdMax*(x:T):untyped = simdMaxReduce(x)
  proc `-`*(x:T):T {.inline,noInit.} =
    forStatic i, 0, L-1:
      neg(result[][i], x[][i])

  map011(T, L, abs, abs)
  map011(T, L, trace, trace)
  map011(T, L, norm2, norm2)
  map011(T, L, sqrt, sqrt)
  map011(T, L, sin, sin)
  map011(T, L, cos, cos)
  map011(T, L, acos, acos)

  map021(T, L, atan2, atan2)
  map021(T, L, min, min)
  map021(T, L, max, max)
  map021(T, L, add, add)
  map021(T, L, sub, sub)
  map021(T, L, mul, mul)
  map021(T, L, divd, divd)
  map021(T, L, `+`, add)
  map021(T, L, `-`, sub)
  map021(T, L, `*`, mul)
  map021(T, L, `/`, divd)

  map110(T, L, assign, assign)
  map110(T, L, neg, neg)
  map110(T, L, iadd, iadd)
  map110(T, L, isub, isub)
  map110(T, L, imul, imul)
  map110(T, L, norm2, norm2)
  map110(T, L, inorm2, inorm2)
  map110(T, L, rsqrt, rsqrt)
  map110(T, L, `+=`, iadd)

  map120(T, L, add, add)
  map120(T, L, sub, sub)
  map120(T, L, mul, mul)
  map120(T, L, divd, divd)
  map120(T, L, imadd, imadd)
  map120(T, L, imsub, imsub)

  map130(T, L, msub, msub)

  #template `assign`*(rr: T; x: T): untyped =
  #  #echotype: rr
  #  #echotype: x
  #  #echo rr
  #  #echo rr[][0]
  #  let xx = x
  #  subst(r,rr):
  #    forStatic i, 0, L-1:
  #      assign(r[][i], xx[][i])
  #proc `:=`*(r: var T; x: T) {.inline.} =
  #  forStatic i, 0, L-1:
  #    r[][i] = x[][i]
  template `:=`*(r: var T; x: T): untyped = assign(r, x)
  #template `:=`*(r: T; x: T) =
  #  let xx = x
  #  forStatic i, 0, L-1:
  #    r[][i] = xx[][i]
  template `:=`*(r: T; x: array[L,B]) =
    let xx = x
    forStatic i, 0, L-1:
      r[][i] = xx[i]
  #proc assign*(r:var T; x:SomeNumber) {.inline,neverInit.} =
  proc assign*(r:var T; x:SomeNumber) {.inline.} =
    #{.emit:"#define memset(a,b,c)".}
    assign(r, x.to(T))
    #assign(r[][0], x)
    #assign(r[][1], r[][0])
  proc assign*(r:var T; x:array[N,SomeNumber]) {.inline.} =
    when compiles(assign(r[][0], unsafeAddr(x[0]))):
      forStatic i, 0, L-1:
        assign(r[][i], unsafeAddr(x[i*N0]))
    else:
      var y = x
      forStatic i, 0, L-1:
        assign(r[][i], unsafeAddr(y[i*N0]))
  proc assign*(r:var array[N,SomeNumber], x:T) {.inline.} =
    subst(i,_):
      forStatic i, 0, L-1:
        assign(addr(r[i*N0]), x[][i])
  template add*(r:var T; x:SomeNumber; y:T) = add(r, x.to(T), y)
  template sub*(r:var T; x:SomeNumber; y:T) = sub(r, x.to(T), y)
  template sub*(r:var T; x:T; y:SomeNumber) = sub(r, x, y.to(T))
  template mul*(r:var T; x:SomeNumber; y:T) = mul(r, x.to(T), y)
  template mul*(r:var T; x:T; y:SomeNumber) = mul(r, x, y.to(T))
  template iadd*(r:var T; x:SomeNumber) = iadd(r, x.to(T))
  template imadd*(r:var T; x:SomeNumber; y:T) = imadd(r, x.to(T), y)
  template imsub*(r:var T; x:SomeNumber; y:T) = imsub(r, x.to(T), y)
  template divd*(r:var T; x:SomeNumber; y:T) = divd(r, x.to(T), y)
  template imul*(r:var T; x:SomeNumber) = imul(r, x.to(T))
  template msub*(r:var T; x:SomeNumber; y,z:T) = msub(r, x.to(T), y, z)
  template `:=`*(r:var T; x:array[N,SomeNumber]) = assign(r, x)
  template `+`*(x:SomeNumber; y:T):T = add(x.to(T), y)
  template `-`*(x:SomeNumber; y:T):T = sub(x.to(T), y)
  template `*`*(x:SomeNumber; y:T):T = mul(x.to(T), y)
  template `*`*(x:T; y:SomeNumber):T = mul(x, y.to(T))
  template `/`*(x:SomeNumber; y:T):T = divd(x.to(T), y)
  template `/`*(x:T; y:SomeNumber):T = divd(x, y.to(T))

  proc `$`*(x:T):string =
    result = "[" & $x[0]
    for i in 1..<N:
      result &= "," & $x[i]
    result &= "]"

  makePerm(1,T,L,N0)
  makePerm(2,T,L,N0)
  makePerm(4,T,L,N0)
  makePerm(8,T,L,N0)
  makePerm(16,T,L,N0)

  makePackP(1,T,L,N0)
  makePackP(2,T,L,N0)
  makePackP(4,T,L,N0)
  makePackP(8,T,L,N0)
  makePackP(16,T,L,N0)

  makePackM(1,T,L,N0)
  makePackM(2,T,L,N0)
  makePackM(4,T,L,N0)
  makePackM(8,T,L,N0)
  makePackM(16,T,L,N0)

  makeBlendP(1,T,L,N0)
  makeBlendP(2,T,L,N0)
  makeBlendP(4,T,L,N0)
  makeBlendP(8,T,L,N0)
  makeBlendP(16,T,L,N0)

  makeBlendM(1,T,L,N0)
  makeBlendM(2,T,L,N0)
  makeBlendM(4,T,L,N0)
  makeBlendM(8,T,L,N0)
  makeBlendM(16,T,L,N0)

  discard """
template pck(op:untyped) =
  proc op*(r:var openArray[SomeNumber], x:T,
           l:var openarray[SomeNumber]) {.inline.} =
    let ra = cast[ptr array[2,type(r[0])]](r[2].addr)
    let la = cast[ptr array[2,type(l[0])]](l[2].addr)
    op(r, x[][0], l)
    op(ra[], x[][1], la[])
pck(packp1)
pck(packm1)
pck(packp2)
pck(packm2)
proc packp4*(r:var openArray[SomeNumber], x:T,
             l:var openarray[SomeNumber]) {.inline.} =
  assign(cast[ptr float64](l.addr), x[][0])
  assign(cast[ptr float64](r.addr), x[][1])
proc packm4*(r:var openArray[SomeNumber], x:T,
             l:var openarray[SomeNumber]) {.inline.} =
  assign(cast[ptr float64](r.addr), x[][0])
  assign(cast[ptr float64](l.addr), x[][1])

template blnd(op:untyped) =
  proc op*(x:var T; r:openArray[SomeNumber];
           l:openArray[SomeNumber]) {.inline.} =
    let ra = cast[ptr array[2,type(r[0])]](unsafeAddr(r[2]))
    let la = cast[ptr array[2,type(l[0])]](unsafeAddr(l[2]))
    op(x[][0], r, l)
    op(x[][1], ra[], la[])
blnd(blendp1)
blnd(blendm1)
blnd(blendp2)
blnd(blendm2)
proc blendp4*(x:var T, r:openArray[SomeNumber],
             l:openarray[SomeNumber]) {.inline.} =
  assign(x[][0], cast[ptr float64](unsafeAddr(l)))
  assign(x[][1], cast[ptr float64](unsafeAddr(r)))
proc blendm4*(x:var T, r:openArray[SomeNumber],
             l:openarray[SomeNumber]) {.inline.} =
  assign(x[][0], cast[ptr float64](unsafeAddr(r)))
  assign(x[][1], cast[ptr float64](unsafeAddr(l)))
"""

when isMainModule:
  makeSimdArray(S2, 2, float32)

  var x,y,z:S2
  assign(x, [1'f32, 2])
  echo x
  y += x
  echo y
