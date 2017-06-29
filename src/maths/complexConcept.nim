import macros
#import metaUtils
#import basicOps
import base
import base/wrapperTypes
import maths/types

# unary ops: unm(-),inc(+=),dec(-=),muleq(*=),diveq(/=)
# binary ops: add(+),sub(-),mul(*),divd(/),minc(+=*),mdec(-=*)
# ternary ops: madd(*+),msub(*-),nmadd(-*+),nmsub(-*-)
# conj,adj,trace,det,transpose,dot,outer,lu,qr,svd,eig,norm2,inv
# sqrt,exp,rsqrt,log,groupProject,groupCheck
# mapX(f,x,r),mapXY(f,x,y,r)

# funcImpl
# renameFunctions(prefix, suffix)
# RealType, ImagType, ComplexType

# types:
#  RealObjConcept, ImagObjConcept, ComplexObjConcept (.re,.im)
#  RealConcept, ImagConcept
#  AsRealObj, AsImagObj, AsComplexObj
#  AsReal, AsImag

proc isKnown*(x:any):bool =
  mixin isReal, isImag, isComplex
  isReal(x) or isImag(x) or isComplex(x)
proc isUnknown*(x:any):bool =
  mixin isKnown
  not isKnown(x)

makeWrapper(AsReal, asReal)
makeWrapper(AsVarReal, asVarReal)
makeWrapper(AsImag, asImag)
makeWrapper(AsVarImag, asVarImag)
makeWrapper(AsComplex, asComplex)
makeWrapper(AsVarComplex, asVarComplex)
type
  AsRI* = AsReal | AsImag
#template `[]`*(x:AsRI; i:SomeInteger):untyped = x[][i]
#template `[]`*(x:AsRI; i,j:SomeInteger):untyped = x[][j,i]
template re*(x:AsReal):untyped = x[]
template im*(x:AsReal):untyped = 0
template re*(x:AsVarReal):untyped = x[]
template im*(x:AsVarReal):untyped = 0
template re*(x:AsImag):untyped = 0
template im*(x:AsImag):untyped = x[]
template re*(x:AsVarImag):untyped = 0
template im*(x:AsVarImag):untyped = x[]
template re*(x:AsComplex):untyped = x[].re
template im*(x:AsComplex):untyped = x[].im
template `re=`*(x:AsComplex; y:any):untyped = x[].re = y
template `im=`*(x:AsComplex; y:any):untyped = x[].im = y
template re*(x:AsVarComplex):untyped = x[].re
template im*(x:AsVarComplex):untyped = x[].im
template `re=`*(x:AsVarComplex; y:any):untyped = x[].re = y
template `im=`*(x:AsVarComplex; y:any):untyped = x[].im = y
template len*(x:AsRI):untyped = x[].len
template nrows*(x:AsRI):untyped = x[].ncols
template ncols*(x:AsRI):untyped = x[].nrows
template numNumbers*(x:AsRI):untyped =
  mixin numNumbers
  numNumbers(x[])
template numNumbers*(x:AsComplex):untyped =
  mixin numNumbers
  2*numNumbers(x.re)
declareReal(AsReal)
declareImag(AsImag)
declareComplex(AsComplex)
declareComplex(AsVarComplex)
template mvLevel*(x:AsRI):untyped =
  mixin mvLevel
  mvLevel(x[])
template deref(x:typed):untyped =
  when isReal(x):
    x.re
  else:
    x

type
  #R1* = concept x
  #  isReal(x)
  #R2* = concept x
  #  isReal(x)
  #R3* = concept x
  #  isReal(x)
  R1*[T] = AsReal[T]
  R2*[T] = AsReal[T]
  R3*[T] = AsReal[T]
  #I1* = concept x
  #  isImag(x)
  #I2* = concept x
  #  isImag(x)
  #I3* = concept x
  #  isImag(x)
  I1*[T] = AsImag[T]
  I2*[T] = AsImag[T]
  I3*[T] = AsImag[T]
  #C1* = concept x
  #  isComplex(x)
  #C2* = concept x
  #  isComplex(x)
  #C3* = concept x
  #  isComplex(x)
  #C4* = concept x
  #  isComplex(x)
  C1*[T] = AsComplex[T]
  C2*[T] = AsComplex[T]
  C3*[T] = AsComplex[T]
  C4*[T] = AsComplex[T]
  RIC1* = R1 | I1 | C1
  RIC2* = R2 | I2 | C2
  RIC3* = R3 | I3 | C3
  VC1* = var C1 | AsVarComplex
  U1* = concept x
    x isnot RIC1
  U2* = concept x
    x isnot RIC2
  U3* = concept x
    x isnot RIC3
  #RIC1 = concept x
  #  isRIC(x)
  #RIC2* = concept x
  #  isRIC(x)
  #RIC3 = concept x
  #  isRIC(x)
  #Iconcept1* = I1
  #Cconcept1* = C1
  #IC1* = Iconcept1 | Cconcept1
  ComplexObj*[T] = object
    re*,im*: T
  ComplexObj2*[T] = ComplexObj[T]
  ComplexType*[T] = AsComplex[ComplexObj[T]]

template complexObj*(x,y: typed): untyped =
  ComplexObj[type(x)](re: x, im: y)
template toDoubleImpl*(x: ComplexObj): untyped = toDoubleX(x)
#proc `$`*(x:C1):string =
#  result = "(" & $x.re & "," & $x.im & ")"

template numberType*[T](x:ComplexObj[T]):untyped = numberType(type(T))
template numberType*[T](x:typedesc[ComplexObj[T]]):untyped = numberType(T)
template numberType*[T](x:AsComplex[T]):untyped = numberType(type(T))
template numberType*[T](x:typedesc[AsComplex[T]]):untyped = numberType(T)

template haveR(x:typed, body:untyped):untyped =
  when not x.isImag:
    block:
      body
template haveI(x:untyped, body:untyped):untyped =
  when not x.isReal:
    block:
      body
template haveR(body:untyped):untyped = haveR(result, body)
template haveI(body:untyped):untyped = haveI(result, body)

template load1*(x:R1):untyped =
  mixin load1
  var r{.noInit.}:AsReal[type(load1(x.re))]
  assign(r, x)
  r
template load1*(x:I1):untyped =
  var r{.noInit.}:AsImag[type(load1(x.im))]
  assign(r, x)
  r
subst(r,_):
  template load1*(x:C1):untyped =
  #subst(x,xx,r,_):
    var r{.noInit.}:ComplexType[type(load1(x.re))]
    assign(r, x)
    r
subst(r,_):
  template load1*(x: AsVarComplex):untyped =
  #subst(x,xx,r,_):
    var r{.noInit.}:ComplexType[type(load1(x.re))]
    assign(r, x)
    r

template map*(x: ComplexObj; f: untyped): untyped =
  let tx = x
  complexObj(f(tx.re), f(tx.im))
template map*(x: AsComplex; f: untyped): untyped = asComplex(f(x[]))
#template map*(x: ; f: untyped): untyped =
#  ComplexType(re:f(x.re),im:f(x.im))

template apply*(result:RIC1; f:untyped) =
  haveR: f(result.re)
  haveI: f(result.im)

#template map*(result:RIC1; f:untyped; x:RIC2):untyped =
#  haveR: `f`(result.re, x.re)
#  haveI: `f`(result.im, x.im)
#template map*(result:var RIC1; f:untyped; x:RIC2):untyped =
#template map*[T:RIC1](result:var T; f:untyped; x:RIC2):untyped =
#template map*(result:untyped; f:untyped; x:untyped):untyped =
#  #bind haveR, haveI
#  mixin re, im
#  haveR: f(result.re, x.re)
#  haveI: f(result.im, x.im)

#proc toDouble*(x:C1):auto {.inline,noInit.} =
#  var r{.noInit.}:ComplexType[type(toDouble(x.re))]
#  r.re = toDouble(x.re)
#  r.im = toDouble(x.im)
#  r

template isWrapper*(x: ComplexObj): untyped = false

template isWrapper*(x: AsComplex): untyped = true
template copyWrapper*(x: AsComplex, y: typed): untyped =
  #static: echo "copyWrapper AsComplex"
  asComplex(y)
template asWrapper*(x: AsComplex, y: typed): untyped =
  #static: echo "asWrapper AsComplex"
  asComplex(y)
template asVarWrapper*(x: AsComplex, y: typed): untyped =
  #static: echo "asVarWrapper AsComplex"
  #var t = asComplex(y)
  #t
  asVar(asComplex(y))

#template masked*(x: AsComplex; msk: int): untyped =
#  asVarComplex(masked(x[],msk))

#template eval*(x: ComplexObj): untyped =
#  let tx = x
#  #ComplexObj(re: eval(tx.re), im: eval(tx.im))
#  ComplexObj(re: tx.re, im: tx.im)
proc eval*(x: ComplexObj): auto =
  let tx = x
  ComplexObj2(re: eval(tx.re), im: eval(tx.im))
  #ComplexObj(re: tx.re, im: tx.im)
template eval*(x: AsComplex): untyped =
  #echoType: x
  #echoType: x[]
  asComplex(eval(x[]))

template trace*(r:var RIC1; x:RIC2):untyped = map(r, trace, x)
proc trace*(x: C1): auto {.inline.} =
  var r{.noInit.}: ComplexType[type(trace(x.re)+trace(x.im))]
  r.re = trace(x.re)
  r.im = trace(x.im)
  r

template norm2*(r:var any; x:I1) =
  norm2(r, x.im)
template norm2*(r:var any; x:C1) =
  mixin norm2, inorm2
  norm2(r, x.re)
  inorm2(r, x.im)
template norm2*(r:var any; x: AsVarComplex) =
  mixin norm2, inorm2
  norm2(r, x.re)
  inorm2(r, x.im)
template norm2*(x:C1):untyped =
  mixin norm2, inorm2
  var r:type(norm2(x.re))
  norm2(r, x.re)
  inorm2(r, x.im)
  r
template norm2*(x: AsVarComplex): untyped =
  mixin norm2, inorm2
  var r:type(norm2(x.re))
  norm2(r, x.re)
  inorm2(r, x.im)
  r
template inorm2*(r:var any; x:C1) =
  mixin inorm2
  inorm2(r, x.re)
  inorm2(r, x.im)

template dR(x:untyped):untyped =
  mixin declaredReal
  when compiles(declaredReal(x)):
    when declaredReal(x):
      x.re
    else:
      x
  else:
    x

template makeUnary(op:untyped):untyped {.dirty.} =
  template op*(r:var R1; x:R2) = op(r.re, x.re)
  template op*(r:var R1; x:C2) = op(r.re, x.re)
  template op*(r:var I1; x:I2) = op(r.im, x.im)
  template op*(r:var I1; x:C2) = op(r.im, x.im)
  template `op CU`*(r:untyped; x:untyped):untyped =
    op(r.re, x)
    op(r.im, 0)
  template op*(r:var C1; x:U2) = `op CU`(r, x)
  template op*(r:AsVarComplex; x:U2) = `op CU`(r, x)
  template op*(r:var C1; x:R2) =
    op(r.re, x.re)
    op(r.im, 0)
  template op*(r:var C1; x:I2) =
    op(r.re, 0)
    op(r.im, x.im)
  template op*(r:var C1; xx:C2) =
    lets(x,xx):
      op(r.re, x.re)
      op(r.im, x.im)
  template op*(r:var C1; xx:AsVarComplex) =
    lets(x,xx):
      op(r.re, x.re)
      op(r.im, x.im)
makeUnary(assign)
makeUnary(neg)
makeUnary(iadd)
makeUnary(isub)

template `:=`*(x:VC1; y:SomeNumber) = assign(x, y)
template `:=`*(x:VC1; y:C2) = assign(x, y)
template `+=`*(r:var RIC1; x:typed):untyped = iadd(r, x)
template `-=`*(r:var RIC1; x:typed):untyped = isub(r, x)


template makeConj(op:untyped):untyped {.dirty.} =
  template op*(r:var R1; x:R2) = op(r.re, x.re)
  template op*(r:var R1; x:C2) = op(r.re, x.re)
  template op*(r:var I1; x:I2) = op(r.im, -x.im)
  template op*(r:var I1; x:C2) = op(r.im, -x.im)
  template op*(r:var C1; x:R2) = op(r.re, x.re)
  template op*(r:var C1; x:I2) = op(r.im, -x.im)
  template op*(r:var C1; x:C2) =
    op(r.re, x.re)
    op(r.im, -x.im)
  template op*(x:RIC1):untyped =
    mixin op
    var r{.noInit.}:type(x)
    op(r, x)
    r
makeConj(conj)
#makeConj(adj)

template imulCU*(r:typed; x:typed):untyped =
  mixin imul
  imul(r.re, x)
  imul(r.im, x)
#proc imul*(r:VC1; x:U2) = imulCU(r, x)
proc imul*(r: var C1; x:U2) = imulCU(r, x)
proc imul*(r: AsVarComplex; x:U2) = imulCU(r, x)

template makeBinary(op:untyped):untyped =
  template op*(r:var R1; x:R2; y:R3) = op(r.re, x.re, y.re)
  template op*(r:var I1; x:I2; y:I3) = op(r.im, x.im, y.im)
  template op*(r:var C1; x:R2; y:I3) =
    op(r.re, x.re, 0)
    op(r.im, 0, y.im)
  template op*(r:var C1; x:I2; y:R3) =
    op(r.re, 0, y.re)
    op(r.im, x.im, 0)
  template op*(r:var C1; x:R2|U2; y:C3) =
    op(r.re, deref(x), y.re)
    op(r.im, 0, y.im)
  template op*(r:var C1; x:C2; y:SomeNumber) =
    op(r.re, x.re, y)
    op(r.im, x.im, 0)
  template op*(r:var C1; x:C2; y:R3) =
    op(r.re, x.re, y.re)
    op(r.im, x.im, 0)
  template op*(r:var C1; x:I2; y:C3) =
    op(r.re, 0, y.re)
    op(r.im, x.im, y.im)
  template op*(r:var C1; x:C2; y:I3) =
    op(r.re, x.re, 0)
    op(r.im, x.im, y.im)
  template op*(r:var C1; x:C2; y:C3) =
    op(r.re, x.re, y.re)
    op(r.im, x.im, y.im)
makeBinary(add)
makeBinary(sub)

setBinop(`+`,add,C1,U2,ComplexType[type(x.re+y)])

proc `+`*(x:U1; y:C2):auto {.inline.} =
  var r{.noInit.}:ComplexType[type(x+y.re)]
  add(r, x, y)
  r
proc `+`*(x:C1; y:C2):auto {.inline.} =
  var r{.noInit.}:ComplexType[type(x.re+y.re)]
  add(r, x, y)
  r
#proc `-`*(x:U1; y:C2):auto {.inline.} =
proc `-`*(x:SomeNumber; y:C2):auto {.inline.} =
  var r{.noInit.}:ComplexType[type(x-y.re)]
  sub(r, x, y)
  r
proc `-`*(x:C1; y:SomeNumber):auto {.inline.} =
  var r{.noInit.}:ComplexType[type(x.re-y)]
  sub(r, x, y)
  r
proc `-`*(x:C1; y:C2):auto {.inline.} =
  var r{.noInit.}:ComplexType[type(x.re-y.re)]
  sub(r, x, y)
  r

template mulCCR*(rr:typed; xx,yy:typed):untyped =
  # r.re = x.re * y
  # r.im = x.im * y
  mixin mul
  subst(r,rr,x,xx,y,yy):
    mul(r.re, x.re, y)
    mul(r.im, x.im, y)

template mulCRC*(rr:typed; xx,yy:typed):untyped =
  # r.re = x * y.re
  # r.im = x * y.im
  mixin mul
  subst(r,rr,x,xx,y,yy):
    mul(r.re, x, y.re)
    mul(r.im, x, y.im)
#template mulCIC*(rr:typed; xx,yy:typed):untyped =
#  # r.re = - x * y.im
#  # r.im =   x * y.re
#  mixin mul, nmul
#  subst(r,rr,x,xx,y,yy):
#    nmul(r.re, x, y.im)
#    mul(r.im, x, y.re)

#proc mul*(r:var RIC1; x:U2; y:RIC3) {.inline.} =
template mul*(r:var RIC1; x:U2; y:RIC3) =
  mixin mul
  mul(r.re, x, y.re)
  mul(r.im, x, y.im)

proc mul*(r:var U1; x:C2; y:C3) {.inline.} =
  # r = x.re*y.re - x.im*y.im
  mixin mul, imsub
  mul(r, x.re, y.re)
  imsub(r, x.im, y.im)
proc mul*(r:var C1; x:C2; y:U3) {.inline.} =
  # r.re = x.re * y.re
  # r.im = x.im * y.re
  mul(r.re, x.re, y)
  mul(r.im, x.im, y)
#proc mul*(r:var C1; x:C2; y:R3) {.inline.} =
template mul*(r:var C1; x:C2; y:R3) =
  # r.re = x.re * y.re
  # r.im = x.im * y.re
  mul(r.re, x.re, y.re)
  mul(r.im, x.im, y.re)
#proc mul*(r:var C1; x:C2; y:C3) {.inline.} =
template mul*(r:var C1; x:C2; y:C3) =
  # r.re = x.re*y.re - x.im*y.im
  # r.im = x.im*y.re + x.re*y.im
  mixin mul, imadd
  mul(r, x, y.re.asReal)
  imadd(r, x, y.im.asImag)

proc redot*(x:C2; y:C3):auto {.inline,noInit.} =
  # x.re*y.re + x.im*y.im
  mixin mul, imadd
  var r{.noInit.}:type(x.re*y.re)
  mul(r, x.re, y.re)
  imadd(r, x.im, y.im)
  r
proc redotinc*(r:var any; x:C2; y:C3) {.inline.} =
  # r += x.re*y.re + x.im*y.im
  mixin imadd
  imadd(r, x.re, y.re)
  imadd(r, x.im, y.im)

proc `*`*(x:SomeNumber; y:C2):auto {.inline.} =
  var r{.noInit.}:ComplexType[type(x*y.re)]
  mul(r, x, y)
  r
proc `*`*(x:C1; y:SomeNumber):auto {.inline.} =
  var r{.noInit.}:ComplexType[type(x.re*y)]
  mul(r, x, y)
  r
proc `*`*(x:C1; y:C2):auto {.inline.} =
  var r{.noInit.}:ComplexType[type(x.re*y.re)]
  mul(r, x, y)
  r

proc divd*(r:var C1; x:C2; y:R3) {.inline.} =
  divd(r.re, x.re, y.re)
  divd(r.im, x.im, y.re)
#proc divd*(r:var C1; x:C2; y:int) {.inline.} =
#  #mixin divd
#  divd(r.re, x.re, y)
#  divd(r.im, x.im, y)
proc divd*(r:var C1; x:C2; y:U3) {.inline.} =
  #mixin divd
  divd(r.re, x.re, y)
  divd(r.im, x.im, y)

proc divd*(r:var C1; x:C2; y:C3) {.inline.} =
  var d = 1/(y.re*y.re+y.im*y.im)
  r.re = d * (x.re*y.re + x.im*y.im)
  r.im = d * (x.im*y.re - x.re*y.im)

proc `/`*(x:C1; y:C2):auto {.inline.} =
  var r{.noInit.}:ComplexType[type(x.re*y.re)]
  divd(r, x, y)
  r

proc `/`*(x:C1; y:SomeNumber):auto =
  mixin divd
  var r{.noInit.}:type(x)
  #echoType: x
  divd(r, x, y)
  r

template imaddCCR*(rr:typed; xx,yy:typed):untyped =
  # r.re += x.re * y
  # r.im += x.im * y
  mixin imadd
  subst(r,rr,x,xx,y,yy):
    imadd(r.re, x.re, y)
    imadd(r.im, x.im, y)
template imaddCCI*(rr:typed; xx,yy:typed):untyped =
  # r.re -= x.im * y
  # r.im += x.re * y
  mixin imadd, imsub
  subst(r,rr,x,xx,y,yy):
    imsub(r.re, x.im, y)
    imadd(r.im, x.re, y)
template imaddCRC*(rr:typed; xx,yy:typed):untyped =
  # r.re += x * y.re
  # r.im += x * y.im
  mixin imadd
  subst(r,rr,x,xx,y,yy):
    imadd(r.re, x, y.re)
    imadd(r.im, x, y.im)
template imaddCIC*(rr:typed; xx,yy:typed):untyped =
  # r.re -= x * y.im
  # r.im += x * y.re
  mixin imadd, imsub
  subst(r,rr,x,xx,y,yy):
    imsub(r.re, x, y.im)
    imadd(r.im, x, y.re)

proc imadd*(r:var U1; x:C2; y:C3) {.inline.} =
  # r += x.re*y.re - x.im*y.im
  mixin imadd
  imadd(r, x.re, y.re)
  imadd(r, -x.im, y.im)
#proc imadd*(r:var C1; x:C2; y:R3) {.inline.} =
template imadd*(r:var C1; x:C2; y:R3) =
  mixin imadd
  imadd(r.re, x.re, y.re)
  imadd(r.im, x.im, y.re)
#proc imadd*(r:var C1; x:C2; y:I3) {.inline.} =
template imadd*(r:var C1; x:C2; y:I3) =
  # r.re -= x.im * y.im
  # r.im += x.re * y.im
  mixin imadd, imsub
  imsub(r.re, x.im, y.im)
  imadd(r.im, x.re, y.im)
#proc imadd*(r:var C1; x:C2; y:C3) {.inline.} =
template imadd*(r:var C1; x:C2; y:C3) =
  # r.re += x.re*y.re - x.im*y.im
  # r.im += x.re*y.im + x.im*y.re
  mixin imadd
  imadd(r, x, y.re.asReal)
  imadd(r, x, y.im.asImag)

template imsubCRC*(rr:typed; xx,yy:typed):untyped =
  # r.re -= x * y.re
  # r.im -= x * y.im
  mixin imsub
  subst(r,rr,x,xx,y,yy):
    imsub(r.re, x, y.re)
    imsub(r.im, x, y.im)
template imsubCIC*(rr:typed; xx,yy:typed):untyped =
  # r.re += x * y.im
  # r.im -= x * y.re
  mixin imadd, imsub
  subst(r,rr,x,xx,y,yy):
    imadd(r.re, x, y.im)
    imsub(r.im, x, y.re)

proc imsub*(r:var C1; x:C2; y:R3) {.inline.} =
  mixin imsub
  imsub(r.re, x.re, y.re)
  imsub(r.im, x.im, y.re)
proc imsub*(r:var C1; x:C2; y:I3) {.inline.} =
  # r.re += x.im * y.im
  # r.im -= x.re * y.im
  mixin imsub, imsub
  imadd(r.re, x.im, y.im)
  imsub(r.im, x.re, y.im)
proc imsub*(r:var C1; x:U2; y:C3) {.inline.} =
  imsub(r.re, x, y.re)
  imsub(r.im, x, y.im)
proc imsub*(r:var C1; x:C2; y:C3) {.inline.} =
  # r.re -= x.re*y.re - x.im*y.im
  # r.im -= x.re*y.im + x.im*y.re
  mixin imsub
  imsub(r, x, y.re.asReal)
  imsub(r, x, y.im.asImag)

proc msub*(r:var C1; x:U2; y:C3; z:C4) {.inline.} =
  mixin msub
  msub(r.re, x, y.re, z.re)
  msub(r.im, x, y.im, z.im)

discard """
proc mdec*(result:var RIC1; x:RIC2; y:RIC3) {.inline.} =
  # r.re -= x.re*y.re - x.im*y.im
  # r.im -= x.re*y.im + x.im*y.re
  haveR:
    mdec(result.re, x.re, y.re)
    minc(result.re, x.im, y.im)
  haveI:
    mdec(result.im, x.im, y.re)
    mdec(result.im, x.re, y.im)
#madd
#msub
#nmadd
#nmsub

"""

when isMainModule:
  import typetraits
  #type RealObj*[T] = object
  #  re: T
  type RealObj*[T] = T
  #type ImagObj*[T] = object
  #  im: T
  type ImagObj*[T] = T
  type ComplexObj*[T] = object
    re: T
    im: T
  type ComplexObj2*[TR,TI] = object
    re: TR
    im: TI
  type Real*[T] = AsReal[RealObj[T]]
  type Imag*[T] = AsImag[ImagObj[T]]
  type Complex*[T] = AsComplex[ComplexObj[T]]
  type Complex2*[TR,TI] = AsComplex[ComplexObj2[TR,TI]]
  template `[]`(x: Complex): untyped = x

  #template newReal*(x): untyped =
  #  Real[type(x)](v:RealObj[type(x)](re:x))
  template newReal*(x): untyped =
    Real[type(x)](v:RealObj[type(x)](x))
  template newImag*(x): untyped =
    Imag[type(x)](v:ImagObj[type(x)](x))
  template newComplex*(x,y): untyped =
    Complex[type(x)](v:ComplexObj[type(x)](re:x,im:y))

  #declareReal(Real)
  #declareImag(Imag)
  #declareComplex(Complex)

  #template add(r: var SomeNumber, x: SomeNumber2, y: SomeNumber3) =
  #  r = (type(r))(x) + (type(r))(y)

  proc `==`*(x:Complex,y:Complex):bool =
    (x.re == y.re) and (x.im == y.im)

  #template assign(r: var SomeNumber, x: RealObj) =
  #  r = x
  #template assign(r: var RealObj, x: RealObj) =
  #  r = x
  #template assign(r: var SomeNumber, x: ImagObj) =
  #  r = x.im
  #template assign(r: var ImagObj, x: ImagObj) =
  #  r.im = x.im
  template `:=`[T1,T2](r: var Real[T1], x: Real[T2]) = assign(r, x)
  template `:=`[T1,T2](r: var Imag[T1], x: Imag[T2]) = assign(r, x)
  template `:=`[T1,T2](r: var Complex[T1], x: Real[T2]) = assign(r, x)
  template `:=`[T1,T2](r: var Complex[T1], x: Imag[T2]) = assign(r, x)
  template `:=`[T1,T2](r: var Complex[T1], x: Complex[T2]) = assign(r, x)

  template set2(op,f,t1,r:untyped):untyped {.dirty.} =
    proc op*[T](x:t1[T]):r[T] {.noInit,inline.} = f(result,x)
  template set22(op,f,t1,r:untyped):untyped {.dirty.} =
    proc op*[T1,T2](x:t1[T1,T2]):r[T1,T2] {.noInit,inline.} = f(result,x)

  #template neg*(r: RealObj, x: RealObj) = neg(r, x)
  set2(`-`, neg, AsReal, AsReal)
  set2(`-`, neg, AsImag, AsImag)
  set2(`-`, neg, AsComplex, AsComplex)

  template set3(op,f,t1,t2,r:untyped):untyped =
    proc op*[T](x:t1[T],y:t2[T]):r[T] {.noInit,inline.} = f(result,x,y)

  template set32(op,f,T1,T2,t1,t2,r:untyped):untyped =
    proc op*[T1,T2](x:t1[T1],y:t2[T2]):r {.noInit,inline.} = f(result,x,y)

  set3(`+`, add, Real, Real, Real)
  set3(`+`, add, Real, Imag, Complex)
  set3(`+`, add, Imag, Real, Complex)
  set3(`+`, add, Imag, Imag, Imag)
  set32(`+`, add, T1,T2, AsReal, AsComplex, AsComplex[T2])
  set32(`+`, add, T1,T2, AsComplex, AsReal, AsComplex[T1])
  set32(`+`, add, T1,T2, AsImag, AsComplex, AsComplex[T2])
  set32(`+`, add, T1,T2, AsComplex, AsImag, AsComplex[T1])
  set3(`+`, add, AsComplex, AsComplex, AsComplex)

  set3(`*`, mul, Real, Real, Real)
  set3(`*`, mul, Real, Imag, Imag)
  set3(`*`, mul, Imag, Real, Imag)
  set3(`*`, mul, Imag, Imag, Real)
  set3(`*`, mul, Real, Complex, Complex)
  set3(`*`, mul, Complex, Real, Complex)
  set3(`*`, mul, Imag, Complex, Complex)
  set3(`*`, mul, Complex, Imag, Complex)
  set3(`*`, mul, AsComplex, AsComplex, AsComplex)

  var f0 = 0.0
  var r0 = newReal(0.0)
  var r1 = newReal(1.0)
  var r2 = newReal(2.0)
  var i0 = newImag(0.0)
  var i1 = newImag(1.0)
  var i2 = newImag(2.0)
  var c0 = newComplex(0.0, 0.0)
  var c1 = newComplex(1.0, 1.0)
  var c2 = newComplex(2.0, 2.0)
  var z1 = c0
  var z2 = c0

  echo(r0 is R1)
  #echo(r0 is RIC1)
  echo(declaredReal(r0))
  #template `:=`*(r:var typed; x:typed):untyped =
  #  echoTyped(x)
  #  mixin assign
  #  assign(r, x)
  static:
    echo r1.type.name
    echo r1.re.type.name
    echo((-r1).type.name)

  r0 := r1
  i0 := i1
  c0 := c1
  c0 := r1
  c0 := i1

  r0 := -r1
  i0 := -i1
  c0 := -c1
  c0 := -r1
  c0 := -i1

  echo isKnown(r0)
  echo isUnknown(r0)
  echo(r0 is RIC1)
  echo(r0 is U1)
  echo(f0 is RIC1)
  echo(f0 is U1)
  r0 += r1

  r0 := r1 + r2
  c0 := r1 + i2
  c0 := i1 + r2
  i0 := i1 + i2
  c0 := r1 + c2
  c0 := c1 + r2
  c0 := i1 + c2
  c0 := c1 + i2
  c0 := c1 + c2

  z1 := c1/c2
  z2 := c2*z1
  doAssert( z2 == c1 )

  z1 := toDouble(c1)

  discard """
  r0 := r1 * r2
  i0 := r1 * i2
  i0 := i1 * r2
  r0 := i1 * i2
  c0 := r1 * c2
  c0 := c1 * r2
  c0 := i1 * c2
  c0 := c1 * i2
  c0 := c1 * c2

  echo($r0.re)
  echo($r0.im)
  echo($i0.re)
  echo($i0.im)
  echo($c0.re)
  echo($c0.im)
"""
  discard """
  #template myBorrow1(op,st,rt:untyped) =
  #  proc op*[T](x:st):rt =
  #    result = cast[rt](op(cast[T](x)))
  template myBorrow2(op,st,rt:untyped) =
  proc op*[T](x,y:st):rt =
    result = cast[rt](op(cast[T](x),cast[T](y)))
#template imBorrow1(op:untyped) =
#  myBorrow1(op, Imag[T], Imag[T])


## unary operators

proc `$`*[T](x:Imag[T]):string =
  result = $(cast[T](x)) & "i"
proc `$`*[T](x:Complex[T]):string =
  result = "(" & $x.re & "," & $x.im & ")"

proc `-`*[T](x:Imag[T]):Imag[T] {.noInit.} =
  result = cast[type(result)](-cast[T](x))
proc `-`*[T](x:Complex[T]):Complex[T] {.noInit.} =
  result.re = -x.re
  result.im = -x.im

proc conj*[T](x:Real[T]):Real[T] {.noInit.} =
  result = cast[type(result)](cast[T](x))
proc conj*[T](x:Imag[T]):Imag[T] {.noInit.} =
  result = cast[type(result)](-cast[T](x))
proc conj*[T](x:Complex[T]):Complex[T] {.noInit.} =
  result.re = x.re
  result.im = -x.im

template `~`*[T](x:Real[T]):untyped = conj(x)
template `~`*[T](x:Imag[T]):untyped = conj(x)
template `~`*[T](x:Complex[T]):untyped = conj(x)


## binary comparisons

template imBorrowBool(op:untyped) =
  myBorrow2(op, Imag[T], bool)

imBorrowBool(`==`)
imBorrowBool(`<`)
imBorrowBool(`<=`)

proc `==`*[T](x,y:Complex[T]):bool =
  result = (x.re==y.re) and (x.im==y.im)


## binary operators
## ReIm, ImRe, ReCo, CoRe, ImIm, ImCo, CoIm, CoCo

proc `+`*[T](x:Real[T], y:Imag[T]):Complex[T] {.noInit,inline.} =
  result.re = x
  result.im = cast[T](y)
proc `+`*[T](x:Imag[T], y:Real[T]):Complex[T] {.noInit,inline.} =
  result.re = y
  result.im = cast[T](x)
proc `+`*[T](x:Real[T], y:Complex[T]):Complex[T] {.noInit,inline.} =
  result.re = x + y.re
  result.im = y.im
proc `+`*[T](x:Complex[T], y:Real[T]):Complex[T] {.noInit,inline.} =
  result.re = x.re + y
  result.im = x.im
proc `+`*[T](x:Imag[T], y:Imag[T]):Imag[T] {.noInit,inline.} =
  result = cast[type(result)](cast[T](x)+cast[T](y))
proc `+`*[T](x:Imag[T], y:Complex[T]):Complex[T] {.noInit,inline.} =
  result.re = y.re
  result.im = cast[T](x) + y.im
proc `+`*[T](x:Complex[T], y:Imag[T]):Complex[T] {.noInit,inline.} =
  result.re = x.re
  result.im = x.im + cast[T](y)
proc `+`*[T](x:Complex[T], y:Complex[T]):Complex[T] {.noInit,inline.} =
  result.re = x.re + y.re
  result.im = x.im + y.im

proc `-`*[T](x:Real[T], y:Imag[T]):Complex[T] {.noInit,inline.} =
  result.re = x
  result.im = -cast[T](y)
proc `-`*[T](x:Imag[T], y:Real[T]):Complex[T] {.noInit,inline.} =
  result.re = -y
  result.im = cast[T](x)
proc `-`*[T](x:Real[T], y:Complex[T]):Complex[T] {.noInit,inline.} =
  result.re = x - y.re
  result.im = -y.im
proc `-`*[T](x:Complex[T], y:Real[T]):Complex[T] {.noInit,inline.} =
  result.re = x.re - y
  result.im = x.im
proc `-`*[T](x:Imag[T], y:Imag[T]):Imag[T] {.noInit,inline.} =
  result = cast[type(result)](cast[T](x)-cast[T](y))
proc `-`*[T](x:Imag[T], y:Complex[T]):Complex[T] {.noInit,inline.} =
  result.re = -y.re
  result.im = cast[T](x) - y.im
proc `-`*[T](x:Complex[T], y:Imag[T]):Complex[T] {.noInit,inline.} =
  result.re = x.re
  result.im = x.im - cast[T](y)
proc `-`*[T](x:Complex[T], y:Complex[T]):Complex[T] {.noInit,inline.} =
  result.re = x.re - y.re
  result.im = x.im - y.im

proc `*`*[T](x:Real[T], y:Imag[T]):Imag[T] {.noInit,inline.} =
  result = cast[type(result)](cast[T](x)*cast[T](y))
proc `*`*[T](x:Imag[T], y:Real[T]):Imag[T] {.noInit,inline.} =
  result = cast[type(result)](cast[T](x)*cast[T](y))
proc `*`*[T](x:Real[T], y:Complex[T]):Complex[T] {.noInit,inline.} =
  result.re = x * y.re
  result.im = x * y.im
proc `*`*[T](x:Complex[T], y:Real[T]):Complex[T] {.noInit,inline.} =
  result.re = x.re * y
  result.im = x.im * y
proc `*`*[T](x:Imag[T], y:Imag[T]):Real[T] {.noInit,inline.} =
  result = cast[type(result)](-cast[T](x)*cast[T](y))
proc `*`*[T](x:Imag[T], y:Complex[T]):Complex[T] {.noInit,inline.} =
  result.re = - cast[T](x) * y.im
  result.im =   cast[T](x) * y.re
proc `*`*[T](x:Complex[T], y:Imag[T]):Complex[T] {.noInit,inline.} =
  result.re = - x.im * cast[T](y)
  result.im =   x.re * cast[T](y)
proc `*`*[T](x:Complex[T], y:Complex[T]):Complex[T] {.noInit,inline.} =
  result.re = x.re*y.re - x.im*y.im
  result.im = x.im*y.re + x.re*y.im

proc `/`*[T](x:Real[T], y:Imag[T]):auto {.noInit,inline.} =
  var r = - cast[T](x) / cast[T](y)
  result = cast[Imag[type(r)]](r)
proc `/`*[T](x:Imag[T], y:Real[T]):auto {.noInit,inline.} =
  var r = cast[T](x) / cast[T](y)
  result = cast[Imag[type(r)]](r)
proc `/`*[T](x:Real[T], y:Complex[T]):auto {.noInit,inline.} =
  var r = cast[T](x)/(y.re*y.re+y.im*y.im)
  var rr = r * (type(r))(y.re)
  var ri = r * (type(r))(-y.im)
  result = cast[Complex[type(r)]]((rr,ri))
proc `/`*[T](x:Complex[T], y:Real[T]):auto {.noInit,inline.} =
  var r = T(1) / cast[T](y)
  var rr = r * (type(r))(x.re)
  var ri = r * (type(r))(x.im)
  result = cast[Complex[type(r)]]((rr,ri))
proc `/`*[T](x:Imag[T], y:Imag[T]):auto {.noInit,inline.} =
  var r = cast[T](x) / cast[T](y)
  result = cast[Real[type(r)]](r)
proc `/`*[T](x:Imag[T], y:Complex[T]):auto {.noInit,inline.} =
  var r = cast[T](x)/(y.re*y.re+y.im*y.im)
  var rr = r * (type(r))(y.im)
  var ri = r * (type(r))(y.re)
  result = cast[Complex[type(r)]]((rr,ri))
proc `/`*[T](x:Complex[T], y:Imag[T]):auto {.noInit,inline.} =
  var r = T(1) / cast[T](y)
  var rr = r * (type(r))(x.im)
  var ri = r * (type(r))(-x.re)
  result = cast[Complex[type(r)]]((rr,ri))
proc `/`*[T](x,y:Complex[T]):auto {.noInit,inline.} =
  var r = T(1)/(y.re*y.re+y.im*y.im)
  var rr = r * (type(r))(x.re*y.re + x.im*y.im)
  var ri = r * (type(r))(x.im*y.re - x.re*y.im)
  result = cast[Complex[type(r)]]((rr,ri))

proc `+=`*[T](x:var Complex[T], y:Real[T]) {.inline.} =
  x.re += y
proc `+=`*[T](x:var Imag[T], y:Imag[T]) {.inline.} =
  x = cast[Imag[T]](cast[T](x)+cast[T](y))
proc `+=`*[T](x:var Complex[T], y:Imag[T]) {.inline.} =
  x.im += cast[T](y)
proc `+=`*[T](x:var Complex[T], y:Complex[T]) {.inline.} =
  x.re += y.re
  x.im += y.im

proc `-=`*[T](x:var Complex[T], y:Real[T]) {.inline.} =
  x.re -= y
proc `-=`*[T](x:var Imag[T], y:Imag[T]) {.inline.} =
  x = cast[Imag[T]](cast[T](x)-cast[T](y))
proc `-=`*[T](x:var Complex[T], y:Imag[T]) {.inline.} =
  x.im -= cast[T](y)
proc `-=`*[T](x:var Complex[T], y:Complex[T]) {.inline.} =
  x.re -= y.re
  x.im -= y.im

proc `*=`*[T](x:var Imag[T], y:Real[T]) {.inline.} =
  x = cast[Imag[T]](cast[T](x)*cast[T](y))
proc `*=`*[T](x:var Complex[T], y:Real[T]) {.inline.} =
  x.re *= y
  x.im *= y
proc `*=`*[T](x:var Complex[T], y:Imag[T]) {.inline.} =
  var r = x.re * cast[T](y)
  x.re = - x.im * cast[T](y)
  x.im = r
proc `*=`*[T](x:var Complex[T], y:Complex[T]) {.inline.} =
  var r = x.im*y.re + x.re*y.im
  x.re = x.re*y.re - x.im*y.im
  x.im = r

proc `/=`*[T](x:var Imag[T], y:Real[T]):auto {.inline.} =
  var r = cast[T](x) / cast[T](y)
  x = cast[Imag[type(r)]](r)
proc `/=`*[T](x:var Complex[T], y:Real[T]):auto {.inline.} =
  var r = T(1) / cast[T](y)
  var rr = r * (type(r))(x.re)
  var ri = r * (type(r))(x.im)
  x = cast[Complex[type(r)]]((rr,ri))
proc `/=`*[T](x:var Complex[T], y:Imag[T]):auto {.inline.} =
  var r = T(1) / cast[T](y)
  var rr = r * (type(r))(x.im)
  var ri = r * (type(r))(-x.re)
  x = cast[Complex[type(r)]]((rr,ri))
proc `/=`*[T](x:var Complex[T]; y:Complex[T]):auto {.inline.} =
  var r = T(1)/(y.re*y.re+y.im*y.im)
  var rr = r * (type(r))(x.re*y.re + x.im*y.im)
  var ri = r * (type(r))(x.im*y.re - x.re*y.im)
  x = cast[Complex[type(r)]]((rr,ri))

when isMainModule:
  template test(T:typedesc) =
    block:
      var zRe = Real(T(0))
      var oRe = Real(T(1))
      var tRe = Real(T(2))
      var xRe = Real(T(1))
      var zIm = Imag(T(0))
      var oIm = Imag(T(1))
      var tIm = Imag(T(2))
      var xIm = Imag(T(1))
      var z = Cmplx(T(0),T(0))
      var o = Cmplx(T(1),T(1))
      var t = Cmplx(T(2),T(2))
      var a = Cmplx(T(1),T(2))
      var b = Cmplx(T,-1,-2)

      assert( $oRe == $T(1) )
      assert( -oRe == -1 )
      assert( ~oRe == 1 )
      assert( oRe.conj == 1 )
      assert( oRe == 1 )
      assert( oRe == oRe )
      assert( oRe != tRe )
      assert( oRe < tRe )
      assert( oRe <= tRe )
      assert( oRe > zRe )
      assert( oRe >= zRe )
      assert( oRe+oRe == tRe )
      assert( oRe-tRe == -oRe )
      assert( oRe*tRe == tRe )
      assert( tRe/oRe == tRe/1 )
      xRe += oRe
      assert( xRe == tRe )
      xRe -= oRe
      assert( xRe == oRe )
      xRe *= oRe
      assert( xRe == oRe )
      when type(xRe) is type(xRe/oRe):
        #echo "testing /[", xRe.type.name, "]"
        xRe /= oRe
        assert( xRe == oRe )

      assert( $oIm == $T(1) & "i" )
      assert( -oIm == Imag[T](-1) )
      assert( ~oIm == Imag[T](-1) )
      assert( oIm.conj == Imag[T](-1) )
      assert( oIm == Imag[T](1) )
      assert( oIm == oIm )
      assert( oIm != tIm )
      assert( oIm < tIm )
      assert( oIm <= tIm )
      assert( oIm > zIm )
      assert( oIm >= zIm )
      assert( oIm+oIm == tIm )
      assert( oIm-tIm == -oIm )
      assert( oIm*tIm == -tRe )
      assert( tIm/oIm == tRe/1 )
      xIm += oIm
      assert( xIm == tIm )
      xIm -= oIm
      assert( xIm == oIm )

      assert( $a == "(" & $T(1) & "," & $T(2) & ")" )
      assert( -a == b )
      assert( ~a == Cmplx(T,1,-2) )
      assert( a.conj == Cmplx(T,1,-2) )
      assert( a == a )
      assert( b* ~a == Cmplx(T,-5,0) )
      assert( a != b )
      assert( a+b == z )
      assert( a-a == z )
      var r1 = a*b
      assert( r1 == Cmplx(T,3,-4) )
      var r2 = a*b*t
      assert( r2 == Cmplx(T,14,-2) )
      var r3 = a/b
      echo r3
      assert( r3 == Cmplx(T(-1)/T(1),T(0)/T(1)) )
      a += b
      assert( a == z )
      a -= b
      assert( a == -b )
      a *= b
      assert( a == -b*b )
      when type(xRe) is type(xRe/oRe):
        a /= b
        assert( a == -b )
      a = -b

      assert( oRe+tIm == a )
      assert( tIm+oRe == a )
      assert( oRe+a == t )
      assert( a+oRe == t )
      assert( oIm+a == Cmplx(T,1,3) )
      assert( a+oIm == Cmplx(T,1,3) )
      assert( b+oRe+tIm == z )

      assert( -oRe-tIm == b )
      assert( -tIm-oRe == b )
      assert( oRe-b == t )
      assert( b-oRe == -t )
      assert( oIm-a == -o )
      assert( a-oIm == o )
      assert( a-oRe-tIm == z )

      assert( oRe*tIm == tIm )
      assert( tIm*oRe == tIm )
      assert( oRe*a == a )
      assert( b*oRe == b )
      assert( oIm*a == Cmplx(T,-2,1) )
      assert( a*oIm == Cmplx(T,-2,1) )
      assert( a*oRe*tIm == Cmplx(T,-4,2) )

      assert( tRe/oIm == Imag(T(-2)/T(1)) )
      assert( tIm/oRe == Imag(T(2)/T(1)) )
      assert( oRe/o == Cmplx(T(1)/T(2),T(-1)/T(2)) )
      assert( a/tRe == Cmplx(T(1)/T(2),T(2)/T(2)) )
      assert( oIm/o == Cmplx(T(1)/T(2),T(1)/T(2)) )
      assert( a/tIm == Cmplx(T(2)/T(2),T(-1)/T(2)) )
      assert( a*oRe*tIm == Cmplx(T,-4,2) )

      a += oRe
      assert( a == t )
      a += oIm
      assert( a == o-b )

      a -= oIm
      assert( a == t )
      a -= oRe
      assert( a == -b )

      oIm *= tRe
      assert( oIm == tIm )
      oIm = Imag(T(1))
      a *= oRe
      assert( a == -b )
      a *= oIm
      assert( a == Cmplx(T,-2,1) )
      a = -b

      a = o/z
      b = 2*a*z
      assert( b == o )

      when type(xRe) is type(xRe/oRe):
        tIm /= oRe
        assert( tIm == oIm+oIm )
        tIm = Imag(T(2))
        a /= oRe
        assert( a == -b )
        a /= oIm
        assert( a == Cmplx(T,2,-1) )
      a = -b

      echo oRe
      echo oIm
      echo a

  test(int)
  test(float)
  test(float32)
  """
