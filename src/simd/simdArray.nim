import macros
import base
import base/metaUtils
import maths/types
export types
#import base/basicOps
getOptimPragmas()

template `[]`*(x: SomeNumber): untyped = x
template `[]`*(x: SomeNumber, i: typed): untyped = x

type
  SimdArrayObj*[N:static int,T] = object
    v*: array[N,T]
    #v* {.align(L*sizeof(B)).}: arrayType(L,B)
  #SimdArray*[N:static int,T] = Simd[SimdArrayObj[N,T]]
type NumberType[T] = numberType(T)

template liftX(f,R,v: untyped) {.dirty.} =
  template f*[N:static int,T](x: SimdArrayObj[N,T]): R = v
  template f*[N:static int,T](x: typedesc[SimdArrayObj[N,T]]): R = v
template liftT(f: untyped) {.dirty.} =
  template f*[N:static int,T](x: SimdArrayObj[N,T]): typedesc =
    mixin f
    f(type T)
  template f*[N:static int,T](x: typedesc[SimdArrayObj[N,T]]): typedesc =
    mixin f
    f(type T)
template liftV(f,v: untyped) {.dirty.} =
  template f*[N:static int,T](x: SimdArrayObj[N,T]): auto = v
  template f*[N:static int,T](x: typedesc[SimdArrayObj[N,T]]): auto = v

liftT(numberType)
liftV(isWrapper, false)
liftV(simdLength, N*simdLength(T))
liftV(numNumbers, N*numNumbers(T))
template simdType*[T:SimdArrayObj](x:T):typedesc = T
template simdType*[T:SimdArrayObj](x:typedesc[T]):typedesc = T

proc to*[T:SimdArrayObj](x: SomeNumber; y: typedesc[T]): T {.alwaysInline,noInit.} =
  forStatic i, 0, y.N-1:
    assign(result[][i], x)

template toSingleImpl*[N:static int,T](x: SimdArrayObj[N,T]): auto =
  when numberType(T) is float32: x
  else: {.error.}
template toDoubleImpl*[N:static int,T](x: SimdArrayObj[N,T]): auto =
  when numberType(T) is float64: x
  else:
    type D = simdObjType(N, toDouble(T))
    var r {.noInit.}: D
    assign(r, x)

template `[]`*(x: SimdArrayObj): auto = x.v
proc `[]`*[N:static int,T](x: SimdArrayObj[N,T], i: SomeInteger):
        NumberType[T] {.alwaysInline.} =
  when T is SomeNumber:
    result = x[][i]
  else:
    let N0 = numNumbers(T)
    let j = i div N0
    let k = i mod N0
    result = x[][j][k]
proc `[]=`*[N:static int,T](r: var SimdArrayObj[N,T], i: SomeInteger, x: auto)
 {.alwaysInline.} =
  when T is SomeNumber:
    r[][i] = T(x)
  else:
    let N0 = numNumbers(T)
    let j = i div N0
    let k = i mod N0
    r[][j][k] = x

proc assign*[N:static int,T:SomeNumber](r: var SimdArrayObj; x: array[N,T]) {.alwaysInline.} =
  when N != r.numNumbers: {.error.}
  const N0 = numNumbers(r.T)
  for i in 0 ..< r.N:
    assign(r[][i], unsafeAddr x[i*N0])

proc assign*[N:static int,T:SomeNumber](r: var array[N,T], x: SimdArrayObj) {.alwaysInline.} =
  when N != x.numNumbers: {.error.}
  const N0 = numNumbers(x.T)
  for i in 0 .. x.N-1:
    assign(addr r[i*N0], x[][i])

# no return value
template p2(f: untyped) {.dirty.} =
  proc f*[S:SimdArrayObj](r: var S, x: S) {.alwaysInline.} =
    mixin f
    for i in 0..<S.N:
      f(r[][i], x[][i])
template p2s(f: untyped) {.dirty.} =
  proc f*[S:SimdArrayObj](r: var S, x: Number) {.alwaysInline.} =
    mixin f
    let t = eval(x)
    for i in 0..<S.N:
      f(r[][i], t)
  p2(f)
template p3(f: untyped) {.dirty.} =
  proc f*[S:SimdArrayObj](r: var S, x,y: S) {.alwaysInline.} =
    mixin f
    for i in 0..<S.N:
      f(r[][i], x[][i], y[][i])
template p3s(f: untyped) {.dirty.} =
  proc f*[S:SimdArrayObj](r: var S, x: S, y: Number) {.alwaysInline.} =
    mixin f
    let t = eval(y)
    for i in 0..<S.N:
      f(r[][i], x[][i], t)
  proc f*[S:SimdArrayObj](r: var S, x: Number, y: S) {.alwaysInline.} =
    mixin f
    let t = eval(x)
    for i in 0..<S.N:
      f(r[][i], t, y[][i])
  p3(f)

p2(neg)
p2(rsqrt)
p2(norm2)
p2s(assign)
p2s(`:=`)
p2s(`+=`)
p2s(`-=`)
p2s(`*=`)
p2s(`/=`)
p2s(iadd)
p2s(isub)
p2s(imul)
p2s(idiv)
p2s(inorm2)
p3s(imadd)
p3s(imsub)
p3s(add)
p3s(sub)
p3s(mul)
p3s(divd)

# with return value
template f1(f: untyped) {.dirty.} =
  proc f*[S:SimdArrayObj](x: S): S {.noInit,alwaysInline.} =
    mixin f
    for i in 0..<S.N:
      result[][i] = f(x[][i])
template f2(f: untyped) {.dirty.} =
  proc f*[N:static int,T](x,y: SimdArrayObj[N,T]):
        SimdArrayObj[N,T] {.noInit,alwaysInline.} =
    mixin f
    for i in 0..<N:
      result[][i] = f(x[][i], y[][i])
template f2s(f: untyped) {.dirty.} =
  proc f*[S:SimdArrayObj](x: S, y: Number): S {.noInit,alwaysInline.} =
    mixin f
    let t = eval(y)
    for i in 0..<S.N:
      result[][i] = f(x[][i], t)
  proc f*[S:SimdArrayObj](x: Number, y: S): S {.noInit,alwaysInline.} =
    mixin f
    let t = eval(x)
    for i in 0..<S.N:
      result[][i] = f(t, y[][i])
  f2(f)

f1(`-`)
f1(abs)
f1(norm2)
f1(inv)
f1(sqrt)
f1(rsqrt)
f1(sin)
f1(cos)
f1(acos)
f1(load1)
f2(atan2)
f2s(`+`)
f2s(`-`)
f2s(`*`)
f2s(`/`)
f2s(add)
f2s(sub)
f2s(mul)
f2s(min)
f2s(max)


proc simdReduce*(r: var SomeNumber; x: SimdArrayObj) {.alwaysInline.} =
  mixin simdReduce
  var y = x[][0]
  forStatic i, 1, x.N-1:
    iadd(y, x[][i])
  r = (type(r))(simdReduce(y))
proc simdReduce*[N:static int,T](x: SimdArrayObj[N,T]):
               NumberType[T] {.noInit,alwaysInline.} =
  simdReduce(result, x)
template simdSum*(r: var SomeNumber; x: SimdArrayObj) = simdReduce(r, x)
template simdSum*(x: SimdArrayObj): auto = simdReduce(x)
proc simdMaxReduce*(r: var SomeNumber; x: SimdArrayObj) {.alwaysInline.} =
  mixin simdMaxReduce
  var y = x[][0]
  forStatic i, 1, x.N-1:
    y = max(y, x[][i])
  r = (type(r))(simdMaxReduce(y))
proc simdMaxReduce*[N:static int,T](x: SimdArrayObj[N,T]):
               NumberType[T] {.noInit,alwaysInline.} =
  simdMaxReduce(result, x)
template simdMax*(r: var SomeNumber; x: SimdArrayObj) = simdMaxReduce(r, x)
template simdMax*(x: SimdArrayObj): auto = simdMaxReduce(x)
proc simdMinReduce*(r: var SomeNumber; x: SimdArrayObj) {.alwaysInline.} =
  mixin simdMinReduce
  var y = x[][0]
  forStatic i, 1, x.N-1:
    y = max(y, x[][i])
  r = (type(r))(simdMinReduce(y))
proc simdMinReduce*[N:static int,T](x: SimdArrayObj[N,T]):
               NumberType[T] {.noInit,alwaysInline.} =
  simdMinReduce(result, x)
template simdMin*(r: var SomeNumber; x: SimdArrayObj) = simdMinReduce(r, x)
template simdMin*(x: SimdArrayObj): auto = simdMinReduce(x)

proc perm*[T:SimdArrayObj](r: var T; x: T, p: int) {.alwaysInline.} =
  mixin perm, assign
  const L = x.N
  const N0 = numNumbers(x.T)
  if N0 > p:
    when N0 > 1:
      forStatic i, 0, L-1:
        perm(r[][i], x[][i], p)
  else:
    let b = (p div N0) and (L-1)
    forStatic i, 0, L-1:
      assign(r[][i], x[][i xor b])
template perm1*[T:SimdArrayObj](r: var T; x: T) = perm(r, x, 1)
template perm2*[T:SimdArrayObj](r: var T; x: T) = perm(r, x, 2)
template perm4*[T:SimdArrayObj](r: var T; x: T) = perm(r, x, 4)
template perm8*[T:SimdArrayObj](r: var T; x: T) = perm(r, x, 8)

template assign[R,X:SomeNumber](r: array[1,R], x: X) = assign(r[0], x)
template assign[R,X:SomeNumber](r: R, x: array[1,X]) = assign(r, x[0])
proc packp*(r: var openArray[SomeNumber], x: SimdArrayObj,
            l: var openarray[SomeNumber], p: int) {.inline.} =
  mixin packp, assign
  const L = x.N
  const N0 = numNumbers(x.T)
  #static: echo N0
  if N0 > p:
    when N0 > 1:
      const N02 = N0 div 2
      let ra = cast[ptr array[L,array[N02,type(r[0])]]](r[0].addr)
      let la = cast[ptr array[L,array[N02,type(l[0])]]](l[0].addr)
      forStatic i, 0, L-1:
        packp(ra[][i], x[][i], la[][i], p)
  else:
    const L2 = L div 2
    let ra = cast[ptr array[L2,array[N0,type(r[0])]]](r[0].addr)
    let la = cast[ptr array[L2,array[N0,type(l[0])]]](l[0].addr)
    let b = (p div N0) and (L-1)
    var ir,il = 0
    forStatic i, 0, L-1:
      if (i and b) == 0:
        assign(la[][il], x[][i])
        inc il
      else:
        assign(ra[][ir], x[][i])
        inc ir
  #echo r, l
template packp1*[R,L:openArray[SomeNumber]](r: var R, x: SimdArrayObj, l: var L) =
  when numNumbers(x) > 1:
    packp(r, x, l, 1)
template packp2*[R,L:openArray[SomeNumber]](r: var R, x: SimdArrayObj, l: var L) =
  when numNumbers(x) > 2:
    packp(r, x, l, 2)
template packp4*[R,L:openArray[SomeNumber]](r: var R, x: SimdArrayObj, l: var L) =
  when numNumbers(x) > 4:
    packp(r, x, l, 4)
template packp8*[R,L:openArray[SomeNumber]](r: var R, x: SimdArrayObj, l: var L) =
  when numNumbers(x) > 8:
    packp(r, x, l, 8)

proc packm*(r: var openArray[SomeNumber], x: SimdArrayObj,
            l: var openarray[SomeNumber], p: int) {.inline.} =
  mixin packm
  const L = x.N
  const N0 = numNumbers(x.T)
  if N0 > p:
    when N0 > 1:
      const N02 = N0 div 2
      let ra = cast[ptr array[L,array[N02,type(r[0])]]](r[0].addr)
      let la = cast[ptr array[L,array[N02,type(l[0])]]](l[0].addr)
      forStatic i, 0, L-1:
        packm(ra[][i], x[][i], la[][i], p)
  else:
    const L2 = L div 2
    let ra = cast[ptr array[L2,array[N0,type(r[0])]]](r[0].addr)
    let la = cast[ptr array[L2,array[N0,type(l[0])]]](l[0].addr)
    let b = (p div N0) and (L-1)
    var ir,il = 0
    forStatic i, 0, L-1:
      if (i and b) == 0:
        assign(ra[][ir], x[][i])
        inc ir
      else:
        assign(la[][il], x[][i])
        inc il
template packm1*[R,L:openArray[SomeNumber]](r: var R, x: SimdArrayObj, l: var L) =
  when numNumbers(x) > 1:
    packm(r, x, l, 1)
template packm2*[R,L:openArray[SomeNumber]](r: var R, x: SimdArrayObj, l: var L) =
  when numNumbers(x) > 2:
    packm(r, x, l, 2)
template packm4*[R,L:openArray[SomeNumber]](r: var R, x: SimdArrayObj, l: var L) =
  when numNumbers(x) > 4:
    packm(r, x, l, 4)
template packm8*[R,L:openArray[SomeNumber]](r: var R, x: SimdArrayObj, l: var L) =
  when numNumbers(x) > 8:
    packm(r, x, l, 8)

proc blendp*(x: var SimdArrayObj; r: openArray[SomeNumber];
             l: openArray[SomeNumber], p: int) {.inline.} =
  mixin blendp
  const L = x.N
  const N0 = numNumbers(x.T)
  if N0 > p:
    when N0 > 1:
      const N02 = N0 div 2
      let ra = cast[ptr array[L,array[N02,type(r[0])]]](unsafeAddr(r[0]))
      let la = cast[ptr array[L,array[N02,type(l[0])]]](unsafeAddr(l[0]))
      forStatic i, 0, L-1:
        blendp(x[][i], ra[][i], la[][i], p)
  else:
    const L2 = L div 2
    let ra = cast[ptr array[L2,array[N0,type(r[0])]]](unsafeAddr(r[0]))
    let la = cast[ptr array[L2,array[N0,type(l[0])]]](unsafeAddr(l[0]))
    let b = (p div N0) and (L-1)
    var ir,il = 0
    forStatic i, 0, L-1:
      if (i and b) == 0:
        assign(x[][i], la[][il])
        inc il
      else:
        assign(x[][i], ra[][ir])
        inc ir
template blendp1*(x: var SimdArrayObj; r: openArray[SomeNumber];
                  l: openArray[SomeNumber]) =
  when numNumbers(x) > 1:
    blendp(x, r, l, 1)
template blendp2*(x: var SimdArrayObj; r: openArray[SomeNumber];
                  l: openArray[SomeNumber]) =
  when numNumbers(x) > 2:
    blendp(x, r, l, 2)
template blendp4*(x: var SimdArrayObj; r: openArray[SomeNumber];
                  l: openArray[SomeNumber]) =
  when numNumbers(x) > 4:
    blendp(x, r, l, 4)
template blendp8*(x: var SimdArrayObj; r: openArray[SomeNumber];
                  l: openArray[SomeNumber]) =
  when numNumbers(x) > 8:
    blendp(x, r, l, 8)

template blendm1*(x: var SimdArrayObj; r: openArray[SomeNumber];
                  l: openArray[SomeNumber]) = blendp1(x, l, r)
template blendm2*(x: var SimdArrayObj; r: openArray[SomeNumber];
                  l: openArray[SomeNumber]) = blendp2(x, l, r)
template blendm4*(x: var SimdArrayObj; r: openArray[SomeNumber];
                  l: openArray[SomeNumber]) = blendp4(x, l, r)
template blendm8*(x: var SimdArrayObj; r: openArray[SomeNumber];
                  l: openArray[SomeNumber]) = blendp8(x, l, r)


# map (var param) (param) (return)

template map011(T,L,op1,op2:untyped):untyped {.dirty.} =
  getOptimPragmas()
  proc op1*(x:T):T {.alwaysInline,noInit.} =
    bind forStatic
    #forStatic i, 0, L-1:
    for i in 0..<L:
      result[][i] = op2(x[][i])
template map021(T,L,op1,op2:untyped):untyped {.dirty.} =
  proc op1*(x,y:T):T {.alwaysInline,noInit.} =
    bind forStatic
    #forStatic i, 0, L-1:
    for i in 0..<L:
      result[][i] = op2(x[][i], y[][i])
template map021x(T1,T2,TR,L,op1,op2:untyped):untyped {.dirty.} =
  proc op1*(x:T1,y:T2):TR {.alwaysInline,noInit.} =
    mixin `[]`
    #bind forStatic
    #forStatic i, 0, L-1:
    for i in 0..<L:
      result[][i] = op2(x[][i], y[][i])
#template map110(T,L,op1,op2:untyped):untyped {.dirty.} =
#  proc op1*(r:var T; x:T) {.inline.} =
#  #template op1*(r: T; xx: T) =
#    #let x = xx
#    forStatic i, 0, L-1:
#      op2(r[][i], x[][i])
macro map110(T,L,op1,op2: untyped): untyped =
  #template tmpldef(f,r,xx,T,body: untyped) =
  #  template f*(r: T; xx: T) = body
  template tmpldef(f,r,x,T,body: untyped) =
    getOptimPragmas()
    proc f*(r: var T; x: T) {.alwaysInline.} = body
  let r = gensym(nskParam,"r")
  let x = gensym(nskParam,"x")
  var body = newStmtList()
  let n = L.intval.int
  for i in 0..<n:
    template bb(t: untyped): untyped =
      newCall(ident"[]",newCall(ident"[]",t),newLit(i))
    body.add newCall(op2, bb(r), bb(x))
  result = getAst(tmpldef(op1,r,x,T,body))
  #echo result.repr

template map120(T,L,op1,op2:untyped):untyped {.dirty.} =
  getOptimPragmas()
  proc op1*(r:var T; x,y:T) {.alwaysInline,noInit.} =
    #bind forStatic
    #forStatic i, 0, L-1:
    for i in 0..<L:
      op2(r[][i], x[][i], y[][i])
template map120x(T1,T2,TR,L,op1,op2:untyped):untyped {.dirty.} =
  proc op1*(r:var TR; x: T1, y: T2) {.inline,noInit.} =
    #bind forStatic
    #forStatic i, 0, L-1:
    for i in 0..<L:
      op2(r[][i], x[][i], y[][i])
template map130(T,L,op1,op2:untyped):untyped {.dirty.} =
  proc op1*(r:var T; x,y,z:T) {.inline,noInit.} =
    bind forStatic
    #forStatic i, 0, L-1:
    for i in 0..<L:
      op2(r[][i], x[][i], y[][i], z[][i])
template makePermX(F,P,T,L,N0) {.dirty.} =
  when N0>P:
    bind map110
    map110(T, L, F, F)
  else:
    proc F*(r:var T; x:T) {.inline.} =
      const b = (P div N0) and (L-1)
      bind forStatic
      forStatic i, 0, L-1:
        assign(r[][i], x[][i xor b])
template makePerm(P,T,L,N0) {.dirty.} =
  bind evalBacktic, makePermX
  evalBacktic:
    makePermX(`"perm" P`,P,T,L,N0)
template makePackPX(F,P,T,L,N0) {.dirty.} =
  # P: perm, T: vector type, L: outer vec len, N0: inner vec len
  when N0>P:
    proc F*(r:var openArray[SomeNumber], x:T,
            l:var openarray[SomeNumber]) {.inline.} =
      const N02 = N0 div 2
      let ra = cast[ptr array[L,array[N02,type(r[0])]]](r[0].addr)
      let la = cast[ptr array[L,array[N02,type(l[0])]]](l[0].addr)
      bind forStatic
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
      bind forStatic
      forStatic i, 0, L-1:
        if (i and b) == 0:
          when N0==1:
            la[][il][0] = x[][i]
          else:
            assign(la[][il], x[][i])
          inc il
        else:
          when N0==1:
            ra[][ir][0] = x[][i]
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
      bind forStatic
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
      bind forStatic
      forStatic i, 0, L-1:
        if (i and b) == 0:
          when N0==1:
            ra[][ir][0] = x[][i]
          else:
            assign(ra[][ir], x[][i])
          inc ir
        else:
          when N0==1:
            la[][il][0] = x[][i]
          else:
            assign(la[][il], x[][i])
          inc il
template makePackP(P,T,L,N0) {.dirty.} =
  bind evalBacktic, makePackPX
  evalBacktic:
    makePackPX(`"packp" P`,P,T,L,N0)
template makePackM(P,T,L,N0) {.dirty.} =
  bind evalBacktic, makePackMX
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
      bind forStatic
      forStatic i, 0, L-1:
        if (i and b) == 0:
          when N0==1:
            x[][i] = la[][il][0]
          else:
            assign(x[][i], la[][il])
          inc il
        else:
          when N0==1:
            x[][i] = ra[][ir][0]
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
      bind forStatic
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
      bind forStatic
      forStatic i, 0, L-1:
        if (i and b) == 0:
          when N0==1:
            x[][i] = ra[][ir][0]
          else:
            assign(x[][i], ra[][ir])
          inc ir
        else:
          when N0==1:
            x[][i] = la[][il][0]
          else:
            assign(x[][i], la[][il])
          inc il
template makeBlendP(P,T,L,N0) {.dirty.} =
  bind evalBacktic, makeBlendPX
  evalBacktic:
    makeBlendPX(`"blendp" P`,P,T,L,N0)
template makeBlendM(P,T,L,N0) {.dirty.} =
  bind evalBacktic, makeBlendMX
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

template arrayType[T](N: typed, x: typedesc[T]): untyped =
  type B {.gensym.} = T.type
  #static: echo "arrayType: ", N, " ", B
  array[N,B]
# T = Simd{S,D}{L} = array[L,B]
# B ~ array[N0,F]
template makeSimdArray*(T:untyped;L:typed;B:typedesc):untyped {.dirty.} =
#template makeSimdArray*(L:typed;B:typedesc;T:untyped) {.dirty.} =
  #static: echo "makeSimdArray: ", L, " ", $B.type
  #makeSimdArray2(T, L, type B, numberType(B), numNumbers(B), L*numNumbers(B))
  makeSimdArray2(L, B, numberType(B), numNumbers(B), L*numNumbers(B), T)
#template makeSimdArray2*(T:untyped;L,B,F,N0,N:typed):untyped {.dirty.} =
#template makeSimdArray2*(T:untyped;L:typed;BB,F:typedesc;N0,N:typed):untyped {.dirty.} =
template makeSimdArray2*(L:typed;B,F:typedesc;N0,N:typed,T:untyped) {.dirty.} =
  getOptimPragmas()
  #static: echo "makeSimdArray2: ", L, " ", $B
  bind map011, map021, map021x, map110, map120, map120x, map130, arrayType
  bind makePerm, makePackP, makePackM, makeBlendP, makeBlendM
  #type B {.gensym.} = typeof(BB)
  #type T* = distinct array[L,B]
  type T* = object
    #v*: array[L,B.type]
    v*: arrayType(L,B)
    #v* {.align(L*sizeof(B)).}: arrayType(L,B)
  #type T* = SimdArrayObj[L,B]
  template isWrapper*(x:typedesc[T]): bool = false
  template isWrapper*(x:T): bool = false
  template numberType*(x:typedesc[T]):typedesc = F
  template numberType*(x:T):typedesc = F
  template numNumbers*(x:typedesc[T]):untyped = N
  template numNumbers*(x:T):untyped = N
  template simdType*(x:typedesc[T]):typedesc = T
  template simdType*(x:T):typedesc = T
  template simdLength*(x:T):untyped = N
  template simdLength*(x:typedesc[T]):untyped = N
  template eval*(x:T): untyped = x
  template eval*(x:typedesc[T]): typedesc = T
  #template `[]`*(x:T):untyped = (array[L,B])(x)
  template `[]`*(x:T):untyped = x.v
  template `[]=`*(x: T; y: typed) = x.v = y
  when B is SomeNumber:
    template `[]`*(x:T; i:SomeInteger):untyped = x[][i div N0]
    template `[]=`*(x:T; i:SomeInteger; y:auto) = x[][i div N0] := y
  else:
    template `[]`*(x:T; i:SomeInteger):untyped = x[][i div N0][i mod N0]
    template `[]=`*(x:T; i:SomeInteger; y:auto) = x[][i div N0][i mod N0] = y
  template load1*(x:T):untyped = x
  proc to*(x:SomeNumber; y:typedesc[T]):T {.alwaysInline,noInit.} =
    bind forStatic
    forStatic i, 0, L-1:
      assign(result[][i], x)
  #echoAst: F
  #static: echo F.treerepr
  #when not(F is float32):
  # should be handled in simd.nim now
  #when F is float64:
  #  template toDoubleImpl*(x: T): untyped = x
  #else:
  #  template toSingleImpl*(x: T): untyped = x
  #when F is float32:
  #  template toSingleImpl*(x: T): untyped = x
  #template toDoubleImpl*(x: T): auto =
  #  mixin simdObjType, assign
  #  when F is float64: x
  #  else:
  #    type D = simdObjType(N, float64)
  #    var r {.noInit.}: D
  #    assign(r, x)
  proc simdReduce*(r: var SomeNumber; x: T) {.inline.} =
    #mixin add
    var y = x[][0]
    forStatic i, 1, L-1:
      iadd(y, x[][i])
    r = (type(r))(simdReduce(y))
  proc simdReduce*(x:T):F {.noInit,inline.} = simdReduce(result, x)
  template simdSum*(r:var SomeNumber; x:T) = simdReduce(r, x)
  template simdSum*(x:T):untyped = simdReduce(x)
  proc simdMaxReduce*(r:var SomeNumber; x:T) {.inline.} =
    mixin simdMaxReduce
    var y = x[][0]
    forStatic i, 1, L-1:
      let c = x[][i]
      y = max(y, c)
    r = (type(r))(simdMaxReduce(y))
  proc simdMaxReduce*(x:T):F {.noinit,inline.} = simdMaxReduce(result, x)
  template simdMax*(r:var SomeNumber; x:T) = simdMaxReduce(r, x)
  template simdMax*(x:T):untyped = simdMaxReduce(x)
  proc simdMinReduce*(r:var SomeNumber; x:T) {.inline.} =
    mixin simdMinReduce
    var y = x[][0]
    forStatic i, 1, L-1:
      let c = x[][i]
      y = min(y, c)
    r = (type(r))(simdMinReduce(y))
  proc simdMinReduce*(x:T):F {.noinit,inline.} = simdMinReduce(result, x)
  template simdMin*(r:var SomeNumber; x:T) = simdMinReduce(r, x)
  template simdMin*(x:T):untyped = simdMinReduce(x)
  proc `-`*(x:T):T {.inline,noInit.} =
    forStatic i, 0, L-1:
      neg(result[][i], x[][i])

  map011(T, L, abs, abs)
  map011(T, L, trace, trace)
  map011(T, L, norm2, norm2)
  map011(T, L, sqrt, sqrt)
  map011(T, L, rsqrt, rsqrt)
  map011(T, L, inv, inv)
  map011(T, L, sin, sin)
  map011(T, L, cos, cos)
  map011(T, L, acos, acos)

  map021(T, L, atan2, atan2)
  map021(T, L, min, min)
  map021(T, L, max, max)
  map021(T, L, add, `+`)
  map021(T, L, sub, `-`)
  map021(T, L, mul, `*`)
  map021(T, L, divd, `/`)
  map021(T, L, `+`, `+`)
  map021(T, L, `-`, `-`)
  map021(T, L, `*`, `*`)
  map021(T, L, `/`, `/`)

  map110(T, L, assign, assign)
  map110(T, L, neg, neg)
  map110(T, L, iadd, iadd)
  map110(T, L, isub, isub)
  map110(T, L, imul, imul)
  map110(T, L, idiv, idiv)
  map110(T, L, norm2, norm2)
  map110(T, L, inorm2, inorm2)
  map110(T, L, rsqrt, rsqrt)
  map110(T, L, `+=`, iadd)
  map110(T, L, `-=`, isub)
  map110(T, L, `*=`, imul)
  map110(T, L, `/=`, idiv)

  map120(T, L, add, add)
  map120(T, L, sub, sub)
  map120(T, L, mul, mul)
  map120(T, L, divd, divd)
  map120(T, L, imadd, imadd)
  map120(T, L, imsub, imsub)

  map130(T, L, msub, msub)

  template `:=`*(r: T; x: T): untyped = assign(r, x)
  template `:=`*(r: T; x: array[L,B]) =
    let xx = x
    forStatic i, 0, L-1:
      r[][i] = xx[i]
  when N==1:
    template assign*(r: SomeNumber, x: T): untyped = r = x[0]
  proc assign*(r: var T; x: SomeNumber) {.alwaysInline.} =
    #assign(r, x.to(T))
    forStatic i, 0, L-1:
      assign(r[][i], x)
  template `:=`*(r: T, x: SomeNumber) = assign(r, x)
  proc assign*(r:var T; x:array[N,SomeNumber]) {.alwaysInline.} =
    when compiles(assign(r[][0], unsafeAddr(x[0]))):
      forStatic i, 0, L-1:
        assign(r[][i], unsafeAddr(x[i*N0]))
    else:
      var y = x
      forStatic i, 0, L-1:
        assign(r[][i], unsafeAddr(y[i*N0]))
  proc assign*(r:var array[N,SomeNumber], x:T) {.alwaysInline.} =
    bind forStatic
    when B is SomeNumber:
      forStatic i, 0, L-1:
        r[i] = x[][i]
    else:
      forStatic i, 0, L-1:
        assign(addr(r[i*N0]), x[][i])
  proc assign*(r: var SimdArrayObj, x: T) {.alwaysInline.} =
    for i in 0..<N:
      r[i] = x[i]
  proc assign*(m: Masked[T], x: SomeNumber) {.inline.} =
    #static: echo "a mask"
    var i = 0
    var b = m.mask
    while b != 0:
      if (b and 1) != 0:
        m.pobj[][i] = x
      b = b shr 1
      i.inc
    #static: echo "end a mask"
  template add*(r:var T; x:SomeNumber; y:T) = add(r, x.to(type(T)), y)
  template add*(r:var T; x:T; y:SomeNumber) = add(r, x, y.to(type(T)))
  template sub*(r:var T; x:SomeNumber; y:T) = sub(r, x.to(type(T)), y)
  template sub*(r:var T; x:T; y:SomeNumber) = sub(r, x, y.to(type(T)))
  template mul*(r:var T; x:SomeNumber; y:T) = mul(r, x.to(type(T)), y)
  #map120x(SomeNumber,T,T,L,mul,mul)
  template mul*(r:var T; x:T; y:SomeNumber) = mul(r, x, y.to(type(T)))
  template iadd*(r:var T; x:SomeNumber) = iadd(r, x.to(type(T)))
  template isub*(r:var T; x:SomeNumber) = isub(r, x.to(type(T)))
  template imadd*(r:var T; x:SomeNumber; y:T) = imadd(r, x.to(type(T)), y)
  template imsub*(r:var T; x:SomeNumber; y:T) = imsub(r, x.to(type(T)), y)
  template divd*(r:var T; x:SomeNumber; y:T) = divd(r, x.to(type(T)), y)
  template imul*(r:var T; x:SomeNumber) = imul(r, x.to(type(T)))
  template idiv*(r:var T; x:SomeNumber) = idiv(r, x.to(type(T)))
  template msub*(r:var T; x:SomeNumber; y,z:T) = msub(r, x.to(type(T)), y, z)
  template `:=`*(r:var T; x:array[N,SomeNumber]) = assign(r, x)
  template `+`*(x:SomeNumber; y:T):T = add(x.to(type(T)), y)
  template `+`*(x:T; y:SomeNumber):T = add(x, y.to(type(T)))
  template `-`*(x:SomeNumber; y:T):T = sub(x.to(type(T)), y)
  template `-`*(x:T; y:SomeNumber):T = sub(x, y.to(type(T)))
  template `*`*(x:SomeNumber; y:T):T = mul(x.to(type(T)), y)
  #map021x(SomeNumber,T,T,L,`*`,`*`)
  template `*`*(x:T; y:SomeNumber):T = mul(x, y.to(type(T)))
  template `/`*(x:SomeNumber; y:T):T = divd(x.to(type(T)), y)
  template `/`*(x:T; y:SomeNumber):T = divd(x, y.to(type(T)))
  template `+=`*(r:var T; x:SomeNumber) = iadd(r, x.to(type(T)))
  template `-=`*(r:var T; x:SomeNumber) = isub(r, x.to(type(T)))
  template `*=`*(r:var T; x:SomeNumber) = imul(r, x.to(type(T)))
  template `/=`*(r:var T; x:SomeNumber) = idiv(r, x.to(type(T)))

  proc `$`*(x:T):string =
    result = "[" & $x[0]
    for i in 1..<N.int:
      result &= "," & $x[i]
    result &= "]"

  makePerm(1,T,L,N0)
  makePerm(2,T,L,N0)
  makePerm(4,T,L,N0)
  makePerm(8,T,L,N0)
  makePerm(16,T,L,N0)
  proc perm*(x: T, p: SomeNumber): T {.inline,noInit.} =
    let b = (p div N0) and (L-1)
    forStatic i, 0, L-1:
      assign(result[][i], perm(x[][i xor b],p))

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
