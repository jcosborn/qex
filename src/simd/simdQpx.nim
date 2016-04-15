import wrapperTypes

type
  SimdS4* = distinct array[4,float32]
  SimdD4* {.importc:"__vector4double".} = object
  SimdSD4 * = SimdS4 | SimdD4
makeDeref(SimdS4, array[4,float32])
makeWrapper(ToSingle, toSingleImpl)
makeWrapper(ToDouble, toDoubleImpl)
template toSingle*(x:SimdS4|ToSingle):expr = x
template toSingle*(x:SimdD4):expr = toSingleImpl(x)
template toSingle*(x:ToDouble):expr = x[]
template toDouble*(x:SimdD4|ToDouble):expr = x
template toDouble*(x:SimdS4):expr = toDoubleImpl(x)
template toDouble*(x:ToSingle):expr = x[]
type SimdSAny* = SimdS4 | ToSingle
type SimdSAny2* = SimdS4 | ToSingle
type SimdDAny* = SimdD4 | ToDouble
type SimdDAny2* = SimdD4 | ToDouble
type SimdAny* = SimdSAny | SimdDAny
type SimdAny2* = SimdSAny2 | SimdDAny2

proc vecLd*(x:SomeInteger; y:SimdS4):SimdD4 {.importC:"vec_ld",noDecl.}
proc vecSt*(x:SimdD4; y:SomeInteger; r:var SimdS4) {.importC:"vec_st",noDecl.}
proc vecSplats*(x:SomeNumber):SimdD4 {.importC:"vec_splats",noDecl.}
proc vecRsp*(x:SimdD4):SimdD4 {.importC:"vec_rsp",noDecl.}
proc vecNeg*(x:SimdD4):SimdD4 {.importC:"vec_neg",noDecl.}
proc vecAdd*(x,y:SimdD4):SimdD4 {.importC:"vec_add",noDecl.}
proc vecSub*(x,y:SimdD4):SimdD4 {.importC:"vec_sub",noDecl.}
proc vecMul*(x,y:SimdD4):SimdD4 {.importC:"vec_mul",noDecl.}
proc vecDiv*(x,y:SimdD4):SimdD4 {.importC:"vec_div",noDecl.}
proc vecMadd*(x,y,z:SimdD4):SimdD4 {.importC:"vec_madd".}
proc vecMsub*(x,y,z:SimdD4):SimdD4 {.importC:"vec_msub".}
proc vecNmadd*(x,y,z:SimdD4):SimdD4 {.importC:"vec_nmadd".}
proc vecNmsub*(x,y,z:SimdD4):SimdD4 {.importC:"vec_nmsub".}
template ld(x:SimdS4):expr = vecLd(0i32, x)
template ld(x:SimdD4):expr = x
template ld(x:ToSingle):expr = ld(x[])
template ld(x:ToDouble):expr = ld(x[])
template ld(x:SomeNumber):expr =
  var r{.noInit.}:SimdD4
  assign(r, x)
  r

template tmpvar*(r:untyped; x:SimdS4):untyped =
  var r{.noInit.}:SimdD4
template tmpvar*(r:untyped; x:SimdD4):untyped =
  var r{.noInit.}:SimdD4
template load*(r:untyped; x:SimdS4):untyped =
  var r{.noInit.}:SimdD4
  assign(r, x)
template load*(r:untyped; x:SimdD4):untyped =
  var r{.noInit.}:SimdD4
  assign(r, x)
template store*(r:var SimdS4; x:untyped):untyped =
  assign(r, x)
template store*(r:var SimdD4; x:untyped):untyped =
  assign(r, x)

proc `$`*(x:SimdS4):string =
  result = "SimdS4[" & $x[][0]
  for i in 1..3:
    result &= "," & $x[][i]
  result &= "]"
#proc `$`*(x:SimdD4):string =
#  result = $(x[])

proc assign*(r:var SimdS4; x:SimdS4) {.inline.} =
  r[] = x[]
proc assign*(r:var SimdS4; x:SimdD4) {.inline.} =
  vecSt(vecRsp(x), 0i32, r)
proc assign*(r:var SimdS4; x:SomeNumber) {.inline.} =
  let t = x.float32
  r[] = [t,t,t,t]
proc assign*(r:var SimdD4; x:SimdS4) {.inline.} =
  r = ld(x)
proc assign*(r:var SimdD4; x:SimdD4) {.inline.} =
  r = x
proc assign*(r:var SimdD4; x:SomeNumber) {.inline.} =
  r = vecSplats(x.double)
template assign*(r:var SimdSD4, x:ToSingle):untyped = assign(r, x[])
template assign*(r:var SimdSD4, x:ToDouble):untyped = assign(r, x[])
template `:=`*(r:var SimdSD4; x:SimdAny|SomeNumber):untyped = assign(r, x)

proc negImpl*(r:var SimdD4; x:SimdD4) {.inline.} =
  r = vecNeg(x)
proc iaddImpl*(r:var SimdD4; x:SimdD4) {.inline.} =
  r = vecAdd(r, x)
proc isubImpl*(r:var SimdD4; x:SimdD4) {.inline.} =
  r = vecSub(r, x)
proc imulImpl*(r:var SimdD4; x:SimdD4) {.inline.} =
  r = vecMul(r, x)
proc idivdImpl*(r:var SimdD4; x:SimdD4) {.inline.} =
  r = vecDiv(r, x)

proc addImpl*(r:var SimdD4; x:SimdD4; y:SimdD4) {.inline.} =
  r = vecAdd(x, y)
proc subImpl*(r:var SimdD4; x:SimdD4; y:SimdD4) {.inline.} =
  r = vecSub(x, y)
proc mulImpl*(r:var SimdD4; x:SimdD4; y:SimdD4) {.inline.} =
  r = vecMul(x, y)
proc divdImpl*(r:var SimdD4; x:SimdD4; y:SimdD4) {.inline.} =
  r = vecDiv(x, y)

proc imaddImpl*(r:var SimdD4; x:SimdD4; y:SimdD4) {.inline.} =
  r = vecMadd(x, y, r)
proc imsubImpl*(r:var SimdD4; x:SimdD4; y:SimdD4) {.inline.} =
  r = vecMsub(x, y, r)
proc inmaddImpl*(r:var SimdD4; x:SimdD4; y:SimdD4) {.inline.} =
  r = vecNmsub(x, y, r)
proc inmsubImpl*(r:var SimdD4; x:SimdD4; y:SimdD4) {.inline.} =
  r = vecNmadd(x, y, r)

proc maddImpl*(r:var SimdD4; x:SimdD4; y:SimdD4; z:SimdD4) {.inline.} =
  r = vecMadd(x, y, z)
proc msubImpl*(r:var SimdD4; x:SimdD4; y:SimdD4; z:SimdD4) {.inline.} =
  r = vecMsub(x, y, z)
proc nmaddImpl*(r:var SimdD4; x:SimdD4; y:SimdD4; z:SimdD4) {.inline.} =
  r = vecNmsub(x, y, z)
proc nmsubImpl*(r:var SimdD4; x:SimdD4; y:SimdD4; z:SimdD4) {.inline.} =
  r = vecNmadd(x, y, z)

template makeUnary(name,op:untyped):untyped =
  template name*(r:var SimdD4; x:any):untyped =
    `name Impl`(r, ld(x))
  template name*(r:var SimdS4; x:any):untyped =
    tmpvar(rr, r)
    `name Impl`(rr, ld(x))
    store(r, rr)
  proc name*(x:SimdAny):SimdD4 {.inline,noInit.} =
    `name Impl`(result, ld(x))
  template op*(x:SimdAny):untyped = name(x)
makeUnary(neg, `-`)

template makeIUnary(name,op:untyped):untyped =
  template name*(r:var SimdD4; x:any):untyped =
    `name Impl`(r, ld(x))
  template name*(r:var SimdS4; x:any):untyped =
    tmpvar(rr, r)
    `name Impl`(rr, ld(x))
    store(r, rr)
  template op*(r:var SimdSD4; x:any):untyped = name(r, x)
makeIUnary(iadd, `+=`)
makeIUnary(isub, `-=`)
makeIUnary(imul, `*=`)
makeIUnary(idivd, `/=`)

template makeBinary(name,op:untyped):untyped =
  template name*(r:var SimdD4; x:any; y:any):untyped =
    `name Impl`(r, ld(x), ld(y))
  template name*(r:var SimdS4; x:any; y:any):untyped =
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
      load(xx, x)
      load(yy, y)
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
