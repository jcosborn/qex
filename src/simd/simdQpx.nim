import base
import base/wrapperTypes

type
  #SimdS4* = distinct array[4,float32]
  SimdS4* = object
    v: array[4,float32]
  #SimdS4* = object
  #  v: ptr array[4,float32]
  SimdD4* {.importc:"__vector4double".} = object
  SimdSD4 * = SimdS4 | SimdD4
  ToSingle = object
  ToDouble = object
makeDeref(SimdS4, array[4,float32])
#makeWrapper(ToSingle, toSingleImpl)
#makeWrapper(ToDouble, toDoubleImpl)
#template toSingle*(x:SimdS4|ToSingle):untyped = x
#template toSingle*(x:SimdD4):untyped = toSingleImpl(x)
#template toSingle*(x:ToDouble):untyped = x[]
#template toDouble*(x:SimdD4|ToDouble):untyped = x
#template toDouble*(x:SimdS4):untyped = toDoubleImpl(x)
#template toDouble*(x:ToSingle):untyped = x[]
type SimdSAny* = SimdS4 | ToSingle
type SimdSAny2* = SimdS4 | ToSingle
type SimdSAny3* = SimdS4 | ToSingle
type SimdDAny* = SimdD4 | ToDouble
type SimdDAny2* = SimdD4 | ToDouble
type SimdDAny3* = SimdD4 | ToDouble
type SimdAny* = SimdSAny | SimdDAny
type SimdAny2* = SimdSAny2 | SimdDAny2
type SimdAny3* = SimdSAny3 | SimdDAny3

template numberType*(x:SimdS4):untyped = float32
template numberType*(x:SimdD4):untyped = float64
template numberType*(x:typedesc[SimdS4]):untyped = float32
template numberType*(x:typedesc[SimdD4]):untyped = float64
template simdType*(x:SimdS4):untyped = SimdS4
template simdType*(x:SimdD4):untyped = SimdD4
template simdLength*(x:SimdS4):untyped = 4
template simdLength*(x:SimdD4):untyped = 4
template simdLength*(x:typedesc[SimdS4]):untyped = 4
template simdLength*(x:typedesc[SimdD4]):untyped = 4
template numNumbers*(x:SimdS4|SimdD4):untyped = simdLength(x)
template numNumbers*(x:typedesc[SimdS4]|typedesc[SimdD4]):untyped = simdLength(x)

proc vecLd*(x:cint; y:ptr float32):SimdD4 {.importC:"vec_ld",noDecl.}
proc vecLd*(x:cint; y:ptr float64):SimdD4 {.importC:"vec_ld",noDecl.}
proc vecLd*(x:cint; y:array[4,float32]):SimdD4 {.importC:"vec_ld",noDecl.}
proc vecLd*(x:cint; y:array[4,float64]):SimdD4 {.importC:"vec_ld",noDecl.}
#proc vecLd*(x:cint; y:SimdD4):SimdD4 {.importC:"vec_ld",noDecl.}
proc vecLd2*(x:cint; y:ptr float32):SimdD4 {.importC:"vec_ld2",noDecl.}
proc vecLd2*(x:cint; y:ptr float64):SimdD4 {.importC:"vec_ld2",noDecl.}
proc vecSt*(x:SimdD4; y:cint; r:ptr float32) {.importC:"vec_st",noDecl.}
proc vecSt*(x:SimdD4; y:cint; r:ptr float64) {.importC:"vec_st",noDecl.}
proc vecSplats*(x:float64):SimdD4 {.importC:"vec_splats",noDecl.}
proc vecExtract*(x:SimdD4; y:cint):float64 {.importC:"vec_extract",noDecl.}
proc vecInsert*(x:float64; y:SimdD4; z:cint):SimdD4 {.importC:"vec_insert",noDecl.}
proc vecGpci*(x:cint):SimdD4 {.importC:"vec_gpci",noDecl.}
proc vecPerm*(x,y,z:SimdD4):SimdD4 {.importC:"vec_perm",noDecl.}
proc vecSt2*(x:SimdD4; y:cint; z:ptr cdouble) {.importC:"vec_st2",noDecl.}
proc vecRsp*(x:SimdD4):SimdD4 {.importC:"vec_rsp",noDecl.}
proc vecNeg*(x:SimdD4):SimdD4 {.importC:"vec_neg",noDecl.}
proc vecAdd*(x,y:SimdD4):SimdD4 {.importC:"vec_add",noDecl.}
proc vecSub*(x,y:SimdD4):SimdD4 {.importC:"vec_sub",noDecl.}
proc vecMul*(x,y:SimdD4):SimdD4 {.importC:"vec_mul",noDecl.}
proc vecDiv*(x,y:SimdD4):SimdD4 {.importC:"vec_swdiv",noDecl.}
proc vecMadd*(x,y,z:SimdD4):SimdD4 {.importC:"vec_madd",noDecl.} # x*y+z
proc vecMsub*(x,y,z:SimdD4):SimdD4 {.importC:"vec_msub",noDecl.} # x*y-z
proc vecNmadd*(x,y,z:SimdD4):SimdD4 {.importC:"vec_nmadd",noDecl.} # -x*y-z
proc vecNmsub*(x,y,z:SimdD4):SimdD4 {.importC:"vec_nmsub",noDecl.} # -x*y+z
proc vecAbs*(x:SimdD4):SimdD4 {.importC:"vec_abs",noDecl.}
proc vecSqrt*(x:SimdD4):SimdD4 {.importC:"vec_swsqrt",noDecl.}
template ld(x:SimdS4):untyped = vecLd(0i32, x[])
template ld(x:SimdD4):untyped = x #vecLd(0i32, x)
template ld(x:ToSingle):untyped = ld(x[])
template ld(x:ToDouble):untyped = ld(x[])
template ld(x:SomeNumber):untyped = vecSplats(x.float64)

template load1*(x:SimdS4):untyped =
  bind ld
  ld(x)
template load1*(x:SimdD4):untyped =
  bind ld
  ld(x)

template `[]`*(x:SimdS4; y:SomeInteger):float32 = x[][y]
template `[]`*(x:SimdD4; y:SomeInteger):float64 =
  vecExtract(x, y.cint)
#proc `[]`*(x:SimdD4; y:SomeInteger):float64 {.inline,noInit.} =
#  #var t{.noInit.}:array[4,float64]
#  #assign(t, x)
#  #vecSt(x, 0, t[0].addr)
#  let t = cast[ptr array[4,float64]](unsafeAddr(x))
#  t[y]
template `[]=`*(x:SimdS4; y:SomeInteger; z:any) =
  x[][y] = z.float32
template `[]=`*(x:SimdD4; y:SomeInteger; z:any) =
  x = vecInsert(z.float64, x, y.cint)
#proc `[]=`*(x:var SimdD4; y:SomeInteger; z:any) {.inline.} =
#  #var t{.noInit.}:array[4,float64]
#  #vecSt(x, 0, t[0].addr)
#  #t[y] = float64(z)
#  #x = vecLd(0, t[0].addr)
#  let t = cast[ptr array[4,float64]](unsafeAddr(x))
#  t[y] = float64(z)

proc `$`*(x:SimdS4):string =
  result = "SimdS4[" & $x[][0]
  for i in 1..3:
    result &= "," & $x[][i]
  result &= "]"
proc `$`*(x:SimdD4):string =
  result = "SimdD4[" & $x[0]
  for i in 1..3:
    result &= "," & $x[i]
  result &= "]"

#proc assign*(r:var SimdS4; x:SimdS4) {.inline.} =
template assign*(r:var SimdS4; x:SimdS4) =
  r[] = x[]
#proc assign*(r:var SimdS4; x:SimdD4) {.inline.} =
template assign*(r:var SimdS4; x:SimdD4) =
  vecSt(vecRsp(x), 0i32, r[][0].addr)
template assign*(r:ptr float32; x:SimdD4) =
  vecSt(vecRsp(x), 0i32, r)
proc assign*(r:var SimdS4; x:SomeNumber) {.inline.} =
  let t = x.float32
  r[] = [t,t,t,t]
#proc assign*(r:var SimdD4; x:SimdS4) {.inline.} =
template assign*(r:var SimdD4; x:SimdS4) =
  r = ld(x)
#proc assign*(r:var SimdD4; x:SimdD4) {.inline.} =
template assign*(r:var SimdD4; x:SimdD4) =
  r = x
#proc assign*(r:var SimdD4; x:SomeNumber) {.inline.} =
template assign*(r:var SimdD4; x:SomeNumber) =
  r = vecSplats(x.float64)
proc assign*(r:var array[4,float64]; x:SimdD4) {.inline.} =
  vecSt(x, 0, r[0].addr)
template assign*(r:var array[4,float32]; x:SimdD4):untyped = assign(r[0].addr, x)
template assign*(r:var array[4,float32]; x:SimdS4):untyped =
  r = x[]
proc assign*(r:var SimdS4; x:array[4,SomeNumber]) {.inline.} =
  for i in 0..3: r[i] = (float32)(x[i])
proc assign*(r:var SimdD4; x:array[4,SomeNumber]) {.inline.} =
  for i in 0..3: r[i] = (float64)(x[i])
template assign*(r:var SimdSD4, x:ToSingle):untyped = assign(r, x[])
template assign*(r:var SimdSD4, x:ToDouble):untyped = assign(r, x[])
template `:=`*(r:var SimdSD4; x:SimdAny|SomeNumber):untyped = assign(r, x)
template `:=`*(r:var SimdSD4; x:array[4,SomeNumber]):untyped = assign(r, x)
template toArray*(x:SimdS4):untyped = x[]
proc toArray*(x:SimdD4):array[4,float64] {.inline,noInit.} =
  assign(result, x)

#proc negImpl*(r:var SimdD4; x:SimdD4) {.inline.} =
template negImpl*(r:var SimdD4; x:SimdD4) =
  r = vecNeg(x)
#proc iaddImpl*(r:var SimdD4; x:SimdD4) {.inline.} =
template iaddImpl*(r:var SimdD4; x:SimdD4) =
  r = vecAdd(r, x)
#proc isubImpl*(r:var SimdD4; x:SimdD4) {.inline.} =
template isubImpl*(r:var SimdD4; x:SimdD4) =
  r = vecSub(r, x)
#proc imulImpl*(r:var SimdD4; x:SimdD4) {.inline.} =
template imulImpl*(r:var SimdD4; x:SimdD4) =
  r = vecMul(r, x)
#proc idivdImpl*(r:var SimdD4; x:SimdD4) {.inline.} =
template idivdImpl*(r:var SimdD4; x:SimdD4) =
  r = vecDiv(r, x)

#proc addImpl*(r:var SimdD4; x:SimdD4; y:SimdD4) {.inline.} =
template addImpl*(r:var SimdD4; x:SimdD4; y:SimdD4) =
  r = vecAdd(x, y)
#proc subImpl*(r:var SimdD4; x:SimdD4; y:SimdD4) {.inline.} =
template subImpl*(r:var SimdD4; x:SimdD4; y:SimdD4) =
  r = vecSub(x, y)
#proc mulImpl*(r:var SimdD4; x:SimdD4; y:SimdD4) {.inline.} =
template mulImpl*(r:var SimdD4; x:SimdD4; y:SimdD4) =
  r = vecMul(x, y)
#proc divdImpl*(r:var SimdD4; x:SimdD4; y:SimdD4) {.inline.} =
template divdImpl*(r:var SimdD4; x:SimdD4; y:SimdD4) =
  r = vecDiv(x, y)

#proc imaddImpl*(r:var SimdD4; x:SimdD4; y:SimdD4) {.inline.} =
template imaddImpl*(r:var SimdD4; x:SimdD4; y:SimdD4) =
  r = vecMadd(x, y, r)
#proc imsubImpl*(r:var SimdD4; x:SimdD4; y:SimdD4) {.inline.} =
template imsubImpl*(r:var SimdD4; x:SimdD4; y:SimdD4) =
  r = vecNMsub(x, y, r)
#proc inmaddImpl*(r:var SimdD4; x:SimdD4; y:SimdD4) {.inline.} =
#template inmaddImpl*(r:var SimdD4; x:SimdD4; y:SimdD4) =
#  r = vecNmsub(x, y, r)
#proc inmsubImpl*(r:var SimdD4; x:SimdD4; y:SimdD4) {.inline.} =
#  r = vecNmadd(x, y, r)

#proc maddImpl*(r:var SimdD4; x:SimdD4; y:SimdD4; z:SimdD4) {.inline.} =
template maddImpl*(r:var SimdD4; x:SimdD4; y:SimdD4; z:SimdD4) =
  r = vecMadd(x, y, z)
#proc msubImpl*(r:var SimdD4; x:SimdD4; y:SimdD4; z:SimdD4) {.inline.} =
template msubImpl*(r:var SimdD4; x:SimdD4; y:SimdD4; z:SimdD4) =
  r = vecMsub(x, y, z)
#proc nmaddImpl*(r:var SimdD4; x:SimdD4; y:SimdD4; z:SimdD4) {.inline.} =
#  r = vecNmsub(x, y, z)
#proc nmsubImpl*(r:var SimdD4; x:SimdD4; y:SimdD4; z:SimdD4) {.inline.} =
#  r = vecNmadd(x, y, z)

template makeUnary(name,op:untyped):untyped =
  proc name*(r:var SimdD4; x:any) {.inline.} =
  #template name*(r:var SimdD4; x:any):untyped =
    `name Impl`(r, ld(x))
  proc name*(r:var SimdS4; x:any) {.inline.} =
  #template name*(r:var SimdS4; x:any):untyped =
    tmpvar(rr, r)
    `name Impl`(rr, ld(x))
    store(r, rr)
  proc name*(x:SimdAny):SimdD4 {.inline,noInit.} =
    `name Impl`(result, ld(x))
  template op*(x:SimdAny):untyped = name(x)
makeUnary(neg, `-`)

template makeIUnary(name,op:untyped):untyped =
  proc name*(r:var SimdD4; x:any) {.inline.} =
  #template name*(r:var SimdD4; x:any):untyped =
    `name Impl`(r, ld(x))
  proc name*(r:var SimdS4; x:any) {.inline.} =
  #template name*(r:var SimdS4; x:any):untyped =
    load2(rr, r)
    `name Impl`(rr, ld(x))
    store(r, rr)
  template op*(r:var SimdSD4; x:any):untyped = name(r, x)
makeIUnary(iadd, `+=`)
makeIUnary(isub, `-=`)
makeIUnary(imul, `*=`)
makeIUnary(idivd, `/=`)

template makeBinary(name,op:untyped):untyped =
  proc name*(r:var SimdD4; x:SimdAny; y:SimdAny2) {.inline.} =
  #template name*(r:var SimdD4; x:SimdAny; y:SimdAny2):untyped =
    `name Impl`(r, ld(x), ld(y))
  proc name*(r:var SimdD4; x:SimdAny; y:SomeNumber) {.inline.} =
  #template name*(r:var SimdD4; x:SimdAny; y:SomeNumber):untyped =
    `name Impl`(r, ld(x), ld(y))
  proc name*(r:var SimdD4; x:SomeNumber; y:SimdAny) {.inline.} =
  #template name*(r:var SimdD4; x:SomeNumber; y:SimdAny):untyped =
    `name Impl`(r, ld(x), ld(y))
  proc name*(r:var SimdS4; x:any; y:any) {.inline.} =
  #template name*(r0:var SimdS4; x0:any; y0:any):untyped =
    #subst(r,r0,x,x0,y,y0,rr,_):
      tmpvar(rr, r)
      `name Impl`(rr, ld(x), ld(y))
      store(r, rr)
  proc name*(x:SimdAny; y:SimdAny2):SimdD4 {.inline,noInit.} =
    `name Impl`(result, ld(x), ld(y))
  proc name*(x:SimdAny; y:SomeNumber):SimdD4 {.inline,noInit.} =
    `name Impl`(result, ld(x), ld(y))
  proc name*(x:SomeNumber; y:SimdAny2):SimdD4 {.inline,noInit.} =
    `name Impl`(result, ld(x), ld(y))
  template op*(x:SimdAny; y:SimdAny2):untyped = name(x,y)
  template op*(x:SimdAny; y:SomeNumber):untyped = name(x,y)
  template op*(x:SomeNumber; y:SimdAny2):untyped = name(x,y)
makeBinary(add, `+`)
makeBinary(sub, `-`)
makeBinary(mul, `*`)
makeBinary(divd, `/`)

template makeIBinary(name):untyped =
  proc name*(r:var SimdD4; x:SimdAny; y:SimdAny2) {.inline.} =
  #template name*(r:var SimdD4; x:SimdAny; y:SimdAny2):untyped =
    `name Impl`(r, ld(x), ld(y))
  proc name*(r:var SimdD4; x:SimdAny; y:SomeNumber) {.inline.} =
  #template name*(r:var SimdD4; x:SimdAny; y:SomeNumber):untyped =
    `name Impl`(r, ld(x), ld(y))
  proc name*(r:var SimdD4; x:SomeNumber; y:SimdAny) {.inline.} =
  #template name*(r:var SimdD4; x:SomeNumber; y:SimdAny):untyped =
    `name Impl`(r, ld(x), ld(y))
  proc name*(r:var SimdS4; x:SimdAny; y:SimdAny2) {.inline.} =
  #template name*(r0:var SimdS4; x0:SimdAny; y0:SimdAny2):untyped =
    #subst(r,r0,x,x0,y,y0,rr,_):
      load2(rr, r)
      `name Impl`(rr, ld(x), ld(y))
      store(r, rr)
  proc name*(r:var SimdS4; x:SimdAny; y:SomeNumber) {.inline.} =
  #template name*(r0:var SimdS4; x0:SimdAny; y0:SomeNumber):untyped =
    #subst(r,r0,x,x0,y,y0,rr,_):
      load2(rr, r)
      `name Impl`(rr, ld(x), ld(y))
      store(r, rr)
  proc name*(r:var SimdS4; x:SomeNumber; y:SimdAny) {.inline.} =
  #template name*(r0:var SimdS4; x0:SomeNumber; y0:SimdAny):untyped =
    #subst(r,r0,x,x0,y,y0,rr,_):
      load2(rr, r)
      `name Impl`(rr, ld(x), ld(y))
      store(r, rr)
makeIBinary(imadd)
makeIBinary(imsub)

template makeTrinary(name):untyped =
  template name*(r:var SimdD4; x:SimdAny; y:SimdAny2; z:SimdAny3):untyped =
    `name Impl`(r, ld(x), ld(y), ld(z))
  template name*(r:var SimdD4; x:SomeNumber; y:SimdAny2; z:SimdAny3):untyped =
    `name Impl`(r, ld(x), ld(y), ld(z))
  template name*(r:var SimdD4; y:SimdAny; x:SomeNumber; z:SimdAny3):untyped =
    `name Impl`(r, ld(x), ld(y), ld(z))
  template name*(r:var SimdD4; x:SimdAny; y:SimdAny2; z:SomeNumber):untyped =
    `name Impl`(r, ld(x), ld(y), ld(z))
  template name*(rr:var SimdS4; xx:any; yy:any; zz:any):untyped =
    subst(r,rr,rt,_,x,xx,y,yy,z,zz):
      tmpvar(rt, r)
      `name Impl`(rt, ld(xx), ld(yy), ld(zz))
      store(r, rt)
makeTrinary(madd)
makeTrinary(msub)

template abs*(x: SimdSD4): untyped = vecAbs(ld(x))
template sqrt*(x: SimdSD4): untyped = vecSqrt(ld(x))

proc rsqrt*(x: SimdD4): SimdD4 {.inline,noInit.} =
  result = 1.0/sqrt(x)
template rsqrt*(x: SimdS4): untyped = rsqrt(ld(x))

template rsqrt*(r: var SimdSD4, x: SimdSD4) = r = rsqrt(x)

template map1(T,N,op: untyped): untyped {.dirty.} =
  proc op*(x: T): T {.inline,noInit.} =
    let t = x.toArray
    var r{.noInit.}: type(t)
    for i in 0..<N:
      r[i] = op(t[i])
    assign(result, r)

map1(SimdSD4, 4, sin)
map1(SimdSD4, 4, cos)
map1(SimdSD4, 4, acos)


proc trace*(x:SimdSD4):SimdSD4 {.inline,noInit.}= x
proc norm2*(x:SimdSD4):SimdD4 {.inline,noInit.} = mul(x,x)
proc norm2*(r:var SimdD4; x:SimdSD4) {.inline.} = mul(r,x,x)
proc inorm2*(r:var SimdD4; x:SimdSD4) {.inline.} = imadd(r,x,x)

proc simdSum*(r:var SomeNumber; x:SimdAny) {.inline.} =
  r = x[0]
  forStatic i, 1, 3:
    r += x[i]
#template simdSum*(x:SimdAny):untyped =
  #var r{.noInit.}:numberType(x)
  #simdSum(r, x)
  #r
proc simdSum*(x:SimdAny):auto {.inline,noInit.} =
  var r:numberType(x)
  simdSum(r, x)
  r

proc perm1*(r:var SimdD4; x:SimdD4) {.inline.} =
  r = vecPerm(x,x,vecGpci(0o1032.cint))
proc perm2*(r:var SimdD4; x:SimdD4) {.inline.} =
  r = vecPerm(x,x,vecGpci(0o2301.cint))
  #r = vecPerm(x,x,vecGpci(0o0123.cint))
proc perm4*(r:var SimdD4; x:SimdD4) {.inline.} =
  assert(false, "perm4 not valid for SimdD4")

template perm1*(r:var SimdD4; x:SimdS4) = perm1(r, ld(x))
template perm2*(r:var SimdD4; x:SimdS4) = perm2(r, ld(x))
template perm4*(r:var SimdD4; x:SimdS4) = perm4(r, ld(x))

proc perm1*(r:var SimdS4; x:SimdS4) {.inline.} =
  let t = x.toArray
  r[0] = t[1]
  r[1] = t[0]
  r[2] = t[3]
  r[3] = t[2]
proc perm2*(r:var SimdS4; x:SimdS4) {.inline.} =
  let t = x.toArray
  r[0] = t[2]
  r[1] = t[3]
  r[2] = t[0]
  r[3] = t[1]
proc perm4*(r:var SimdS4; x:SimdS4) {.inline.} =
  assert(false, "perm4 not valid for SimdS4")

proc perm*(x:SimdD4; p:int):SimdD4 {.noInit.} =
  case p
  of 0: assign(result, x)
  of 1: perm1(result, ld(x))
  of 2: perm2(result, ld(x))
  of 4: perm4(result, ld(x))
  else: discard
template perm*(x:SimdS4; p:int):SimdD4 = perm(ld(x), p)

proc packp1*(r:var openArray[SomeNumber]; x:SimdD4;
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  l[0] = t[0]
  r[0] = t[1]
  l[1] = t[2]
  r[1] = t[3]
  #vecSt2(vecPerm(x,x,vecGpci(0o1302.cint)),0.cint,r[0].addr)
  #vecSt2(vecPerm(x,x,vecGpci(0o0213.cint)),0.cint,l[0].addr)
proc packm1*(r:var openArray[SomeNumber]; x:SimdD4;
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  r[0] = t[0]
  l[0] = t[1]
  r[1] = t[2]
  l[1] = t[3]
  #vecSt2(vecPerm(x,x,vecGpci(0o0213.cint)),0.cint,r[0].addr)
  #vecSt2(vecPerm(x,x,vecGpci(0o1302.cint)),0.cint,l[0].addr)
proc packp2*(r:var openArray[SomeNumber]; x:SimdD4;
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  l[0] = t[0]
  l[1] = t[1]
  r[0] = t[2]
  r[1] = t[3]
  #vecSt2(vecPerm(x,x,vecGpci(0o2301.cint)),0.cint,r[0].addr)
  #vecSt2(x,0.cint,l[0].addr)
proc packm2*(r:var openArray[SomeNumber]; x:SimdD4;
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  r[0] = t[0]
  r[1] = t[1]
  l[0] = t[2]
  l[1] = t[3]
  #vecSt2(x,0.cint,r[0].addr)
  #vecSt2(vecPerm(x,x,vecGpci(0o2301.cint)),0.cint,l[0].addr)
proc packp4*(r:var openArray[SomeNumber]; x:SimdD4;
             l:var openArray[SomeNumber]) {.inline.} =
  assert(false, "packp4 not valid for SimdD4")
proc packm4*(r:var openArray[SomeNumber]; x:SimdD4;
             l:var openArray[SomeNumber]) {.inline.} =
  assert(false, "packm4 not valid for SimdD4")
proc packp8*(r:var openArray[SomeNumber]; x:SimdD4;
             l:var openArray[SomeNumber]) {.inline.} =
  assert(false, "packp8 not valid for SimdD4")
proc packm8*(r:var openArray[SomeNumber]; x:SimdD4;
             l:var openArray[SomeNumber]) {.inline.} =
  assert(false, "packm8 not valid for SimdD4")

proc packp1*(r:var openArray[SomeNumber]; x:SimdS4;
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  l[0] = t[0]
  r[0] = t[1]
  l[1] = t[2]
  r[1] = t[3]
proc packm1*(r:var openArray[SomeNumber]; x:SimdS4;
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  r[0] = t[0]
  l[0] = t[1]
  r[1] = t[2]
  l[1] = t[3]
proc packp2*(r:var openArray[SomeNumber]; x:SimdS4;
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  l[0] = t[0]
  l[1] = t[1]
  r[0] = t[2]
  r[1] = t[3]
proc packm2*(r:var openArray[SomeNumber]; x:SimdS4;
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  r[0] = t[0]
  r[1] = t[1]
  l[0] = t[2]
  l[1] = t[3]
proc packp4*(r:var openArray[SomeNumber]; x:SimdS4;
             l:var openArray[SomeNumber]) {.inline.} =
  assert(false, "packp4 not valid for SimdS4")
proc packm4*(r:var openArray[SomeNumber]; x:SimdS4;
             l:var openArray[SomeNumber]) {.inline.} =
  assert(false, "packm4 not valid for SimdS4")

proc blendp1*(x:var SimdSD4; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = l[0]
  t[1] = r[0]
  t[2] = l[1]
  t[3] = r[1]
  assign(x, t)
  #let tr = vecLd2(0.cint,unsafeAddr(r[0]))
  #let tl = vecLd2(0.cint,unsafeAddr(l[0]))
  #let y = vec_perm(tr,tl,vec_gpci(0o4051))
  #assign(x, y)
proc blendm1*(x:var SimdSD4; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = r[0]
  t[1] = l[0]
  t[2] = r[1]
  t[3] = l[1]
  assign(x, t)
proc blendp2*(x:var SimdSD4; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = l[0]
  t[1] = l[1]
  t[2] = r[0]
  t[3] = r[1]
  assign(x, t)
  #let tr = vecLd2(0.cint,unsafeAddr(r[0]))
  #let tl = vecLd2(0.cint,unsafeAddr(l[0]))
  #let y = vec_perm(tr,tl,vec_gpci(0o4501))
  #assign(x, y)
proc blendm2*(x:var SimdSD4; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = r[0]
  t[1] = r[1]
  t[2] = l[0]
  t[3] = l[1]
  assign(x, t)
  #let tr = vecLd2(0.cint,unsafeAddr(r[0]))
  #let tl = vecLd2(0.cint,unsafeAddr(l[0]))
  #let y = vec_perm(tr,tl,vec_gpci(0o0145))
  #assign(x, y)
proc blendp4*(x:var SimdSD4; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
  assert(false, "blendp4 not valid for SimdD4")
proc blendm4*(x:var SimdSD4; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
  assert(false, "blendm4 not valid for SimdD4")
proc blendp8*(x:var SimdSD4; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
  assert(false, "blendp8 not valid for SimdD4")
proc blendm8*(x:var SimdSD4; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
  assert(false, "blendm8 not valid for SimdD4")

when isMainModule:
  var s1,s2,s3:SimdS4
  var d1,d2,d3:SimdD4
  s1 := 1.0
  s2 := 2.0

  s3 := s1 + s2
  echo s3
  s3 := s1 - s2
  echo s3
  s3 := s1 + s2 * s1 - s2
  echo s3

  macro dup(x:untyped):auto =
    echo x.repr
    x
  proc foo1(r:var any; x,y:any) {.inline.} =
    dup:
      tmpvar(rr, r)
      load2(xx, x)
      load2(yy, y)
      rr = xx + yy
      store(r, rr)

  foo1(s3, s1, s2)
  echo s3
  foo1(d3, d1, d2)
  s3 := d3
  echo s3

  s3 := s1*d1 - s2*d2
  d3 := s1*d2 + s2*d1
  echo s3
  s3 := d3
  echo s3
