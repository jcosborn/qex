import macros
#import safecall

type
  RealProxy*[T] = object
    v*: T
  RealProxy2*[T] = RealProxy[T]
  RealProxy3*[T] = RealProxy[T]
  RealProxy4*[T] = RealProxy[T]

  ImagProxy*[T] = object
    v*: T
  ImagProxy2*[T] = ImagProxy[T]
  ImagProxy3*[T] = ImagProxy[T]
  ImagProxy4*[T] = ImagProxy[T]

  ComplexProxy*[T] = object
    v*: T
  ComplexProxy2*[T] = ComplexProxy[T]
  ComplexProxy3*[T] = ComplexProxy[T]
  ComplexProxy4*[T] = ComplexProxy[T]

template `[]`*(x: RealProxy): untyped = x.v
macro `[]`*(x: RealProxy{nkObjConstr}): auto =
  #echo x.treerepr
  result = x[1][1]
  #echo result.treerepr
template `[]`*(x: ImagProxy): untyped = x.v
macro `[]`*(x: ImagProxy{nkObjConstr}): auto =
  #echo x.treerepr
  result = x[1][1]
  #echo result.treerepr
template `[]`*(x: ComplexProxy): untyped = x.v
macro `[]`*(x: ComplexProxy{nkObjConstr}): auto =
  #echo x.treerepr
  result = x[1][1]
  #echo result.treerepr

template `[]=`*(x: RealProxy, y: typed): untyped =
  #when x.T is type(y):
  #  x.v = y
  #else:
  mixin `:=`
  x.v := y
template `[]=`*(x: ImagProxy, y: typed): untyped =
  #when x.T is type(y):
  #  x.v = y
  #else:
    x.v := y
template `[]=`*(x: ComplexProxy, y: typed): untyped =
  #when x.T is type(y):
  #  x.v = y
  #else:
    x.v := y

proc `$`*(x: RealProxy): string =
  result = $x[]
proc `$`*(x: ImagProxy): string =
  result = $x[] & "I"
proc `$`*(x: ComplexProxy): string =
  #result = "(" & $x.reX & "," & $x.imX & ")"
  result = $x[].re
  var t = $x[].im & "I"
  if t[0]!='-': result.add "+"
  result.add t

template newRealProxy*(x: typed): untyped =
  RealProxy[type(x)](v: x)
template newRealProxy*(x: typed{call}): untyped =
  let t = x
  RealProxy[type(t)](v: t)

template newImagProxy*(x: typed): untyped =
  ImagProxy[type(x)](v: x)
template newImagProxy*(x: typed{call}): untyped =
  let t = x
  ImagProxy[type(t)](v: t)

template newComplexProxy*(x: typed): untyped =
  ComplexProxy[type(x)](v: x)
template newComplexProxy*(x: typed{call}): untyped =
  let t = x
  ComplexProxy[type(t)](v: t)

template newRealP*(x: typed): untyped =
  newRealImpl(x)
template newImagP*(x: typed): untyped =
  newImagImpl(x)
template newComplexP*(x,y: typed): untyped =
  newComplexImpl(x, y)

template re*(x: RealProxy): untyped = x[]
template im*(x: RealProxy): untyped = 0
template re*(x: ImagProxy): untyped = 0
template im*(x: ImagProxy): untyped = x[]
template re*(x: ComplexProxy): untyped = x[].re
template im*(x: ComplexProxy): untyped = x[].im

template `re=`*(x: RealProxy, y: untyped): untyped = x[] = y
template `im=`*(x: RealProxy, y: untyped): untyped = discard
template `re=`*(x: ImagProxy, y: untyped): untyped = discard
template `im=`*(x: ImagProxy, y: untyped): untyped = x[] = y
template `re=`*(x: ComplexProxy, y: untyped): untyped = x[].re = y
template `im=`*(x: ComplexProxy, y: untyped): untyped = x[].im = y

template assign*(x: RealProxy, y: RealProxy2): untyped =
  x[] = y[]
template assign*(x: ImagProxy, y: ImagProxy2): untyped =
  x[] = y[]
template assign*(x: ComplexProxy, y: RealProxy2): untyped =
  x[].re = y[]
  x[].im = 0
template assign*(x: ComplexProxy, y: ImagProxy2): untyped =
  x[].re = 0
  x[].im = y[]
template assign*(x: ComplexProxy, y: ComplexProxy2): untyped =
  let t = y
  x[].re = t[].re
  x[].im = t[].im
template assign*(x: ComplexProxy, y: ComplexProxy2{atom}): untyped =
  x[].re = y[].re
  x[].im = y[].im

template `:=`*(x: RealProxy, y: RealProxy2): untyped = assign(x,y)
template `:=`*(x: ImagProxy, y: ImagProxy2): untyped = assign(x,y)
template `:=`*(x: ComplexProxy, y: RealProxy2): untyped = assign(x,y)
template `:=`*(x: ComplexProxy, y: ImagProxy2): untyped = assign(x,y)
template `:=`*(x: var ComplexProxy, y: ComplexProxy2): untyped = assign(x,y)

# pos, neg, conj, adj, transpose, trace, norm2, inv

template unaryOverloads(op,fn,implR,implI: untyped) {.dirty.} =
  #template fn*(x: RealProxy): untyped = newRealP(implR(x,0))
  proc fn*(x: RealProxy): auto {.inline,noInit.} = newRealP(implR(x[],0))
  template op*(x: RealProxy): untyped = fn(x)

  #template fn*(x: ImagProxy): untyped = newImagP(implI(0,x))
  proc fn*(x: ImagProxy): auto {.inline,noInit.} = newImagP(implI(0,x[]))
  template op*(x: ImagProxy): untyped = fn(x)

  #template `fn U`*(x: ComplexProxy): untyped =
  #  newComplexP(implR(x.re,x.im), implI(x.re,x.im))
  #template fn*(x: ComplexProxy): untyped =
  #  safecall(`fn U`, x)
  #proc fn*(x: ComplexProxy): auto {.inline,noInit.} =
  #  newComplexP(implR(x.re,x.im), implI(x.re,x.im))
  template fn*(xx: ComplexProxy): untyped =
    let x = xx
    newComplexP(implR(x.re,x.im), implI(x.re,x.im))
  template op*(x: ComplexProxy): untyped = fn(x)

template unaryOverloadsR(op,fn,implR,implI: untyped) {.dirty.} =
  #template fn*(x: RealProxy): untyped = newRealP(implR(x,0))
  proc fn*(x: RealProxy): auto {.inline,noInit.} = newRealP(implR(x[],0))
  template op*(x: RealProxy): untyped = fn(x)

  #template fn*(x: ImagProxy): untyped = newRealP(implI(0,x))
  proc fn*(x: ImagProxy): auto {.inline,noInit.} = newRealP(implI(0,x[]))
  template op*(x: ImagProxy): untyped = fn(x)

  #template `fn U`*(x: ComplexProxy): untyped =
  #  newRealP(implR(x.re,x.im)+implI(x.re,x.im))
  #template fn*(x: ComplexProxy): untyped =
  #  safecall(`fn U`, x)
  proc fn*(x: ComplexProxy): auto {.inline,noInit.} =
    newRealP(implR(x.re,x.im)+implI(x.re,x.im))
  template op*(x: ComplexProxy): untyped = fn(x)

template posRealU(xr,xi: untyped): untyped = xr
template posImagU(xr,xi: untyped): untyped = xi
unaryOverloads(`+`, pos, posRealU, posImagU)

template negRealU(xr,xi: untyped): untyped = -xr
template negImagU(xr,xi: untyped): untyped = -xi
unaryOverloads(`-`, neg, negRealU, negImagU)

template conjRealU(xr,xi: untyped): untyped =
  mixin conj
  conj(xr)
template conjImagU(xr,xi: untyped): untyped =
  mixin conj
  -conj(xi)
unaryOverloads(`~`, conj, conjRealU, conjImagU)

template adjRealU(xr,xi: untyped): untyped =
  mixin adj
  adj(xr)
template adjImagU(xr,xi: untyped): untyped =
  mixin adj
  #-adj(xi)
  -adj(xi)
unaryOverloads(`@`, adj, adjRealU, adjImagU)

template transposeRealU(xr,xi: untyped): untyped =
  mixin transpose
  transpose(xr)
template transposeImagU(xr,xi: untyped): untyped =
  mixin transpose
  transpose(xi)
unaryOverloads(`<>`, transpose, transposeRealU, transposeImagU)

template traceRealU(xr,xi: untyped): untyped =
  mixin trace
  trace(xr)
template traceImagU(xr,xi: untyped): untyped =
  mixin trace
  trace(xi)
unaryOverloads(`%`, trace, traceRealU, traceImagU)

template norm2RealU(xr,xi: untyped): untyped =
  mixin norm2
  norm2(xr)
template norm2ImagU(xr,xi: untyped): untyped =
  mixin norm2
  norm2(xi)
unaryOverloadsR(`|`, norm2, norm2RealU, norm2ImagU)

#template inv*(x: RealProxy): untyped = newRealP(x[].inv)
proc inv*(x: RealProxy): auto {.inline,noInit.} =
  mixin inv
  newRealP(x[].inv)
template `/`*(x: RealProxy): untyped = inv(x)
#template inv*(x: ImagProxy): untyped = newImagP(-x[].inv)
proc inv*(x: ImagProxy): auto {.inline,noInit.} =
  mixin inv
  newImagP(-x[].inv)
template `/`*(x: ImagProxy): untyped = inv(x)
#template invComplexU(x: untyped): untyped = x.adj * x.norm2.inv
#template inv*(x: ComplexProxy): untyped = safecall(invComplexU, x)
proc inv*(x: ComplexProxy): auto {.inline,noInit.} =
  mixin inv
  x.adj * x.norm2.inv
template `/`*(x: ComplexProxy): untyped = inv(x)

# add, sub, mul, divd

template binaryOverloadsAddSub(op,fn: untyped) {.dirty.} =
  template fn*(x: RealProxy, y: RealProxy2): untyped = newRealP(op(x[],y[]))
  template op*(x: RealProxy, y: RealProxy2): untyped = fn(x,y)

  template fn*(x: ImagProxy, y: ImagProxy2): untyped = newImagP(op(x[],y[]))
  template op*(x: ImagProxy, y: ImagProxy2): untyped = fn(x,y)

  template fn*(x: RealProxy, y: ImagProxy2): untyped = newComplexP(x[],op(y[]))
  template op*(x: RealProxy, y: ImagProxy2): untyped = fn(x,y)
  template fn*(x: ImagProxy, y: RealProxy2): untyped = newComplexP(op(y[]),x[])
  template op*(x: ImagProxy, y: RealProxy2): untyped = fn(x,y)

  #template `fn RCU`*(x: RealProxy, y: ComplexProxy2): untyped =
  #  newComplexP(op(x[],y.re), op(y.im))
  #template fn*(x: RealProxy, y: ComplexProxy2): untyped =
  #  safecall(`fn RCU`, x, y)
  #proc fn*(x: RealProxy, y: ComplexProxy2): auto {.inline,noInit.} =
  #  newComplexP(op(x[],y.re), op(y.im))
  template fn*(x: RealProxy, yy: ComplexProxy2): untyped =
    let y = yy
    newComplexP(op(x[],y.re), op(y.im))
  template op*(x: RealProxy, y: ComplexProxy2): untyped = fn(x,y)
  #template `fn CRU`*(x: ComplexProxy, y: RealProxy2): untyped =
  #  newComplexP(op(x.re,y[]), x.im)
  #template fn*(x: ComplexProxy, y: RealProxy2): untyped =
  #  safecall(`fn CRU`, x, y)
  proc fn*(x: ComplexProxy, y: RealProxy2): auto {.inline,noInit.} =
    newComplexP(op(x.re,y[]), x.im)
  template op*(x: ComplexProxy, y: RealProxy2): untyped = fn(x,y)

  #template `fn ICU`*(x: ImagProxy, y: ComplexProxy2): untyped =
  #  newComplexP(op(y.re), op(x.im,y.im))
  #template fn*(x: ImagProxy, y: ComplexProxy2): untyped =
  #  safecall(`fn ICU`, x, y)
  proc fn*(x: ImagProxy, y: ComplexProxy2): auto {.inline,noInit.} =
    newComplexP(op(y.re), op(x.im,y.im))
  template op*(x: ImagProxy, y: ComplexProxy2): untyped = fn(x,y)
  #template `fn CIU`*(x: ComplexProxy, y: ImagProxy2): untyped =
  #  newComplexP(x.re, op(x.im,y[]))
  #template fn*(x: ComplexProxy, y: ImagProxy2): untyped =
  #  safecall(`fn CIU`, x, y)
  proc fn*(x: ComplexProxy, y: ImagProxy2): auto {.inline,noInit.} =
    newComplexP(x.re, op(x.im,y[]))
  template op*(x: ComplexProxy, y: ImagProxy2): untyped = fn(x,y)

  #template `fn CCU`*(x: ComplexProxy, y: ComplexProxy2): untyped =
  #  newComplexP(op(x.re,y.re), op(x.im,y.im))
  #template fn*(x: ComplexProxy, y: ComplexProxy2): untyped =
  #  safecall(`fn CCU`, x, y)
  #proc fn*(x: ComplexProxy, y: ComplexProxy2): auto {.inline,noInit.} =
  #  newComplexP(op(x.re,y.re), op(x.im,y.im))
  template fn*(xx: ComplexProxy, yy: ComplexProxy2): untyped =
    let x = xx
    let y = yy
    newComplexP(op(x.re,y.re), op(x.im,y.im))
  template op*(x: ComplexProxy, y: ComplexProxy2): untyped = fn(x,y)

binaryOverloadsAddSub(`+`, add)
binaryOverloadsAddSub(`-`, sub)

template binaryOverloadsMul(op,fn: untyped) {.dirty.} =
  template fn*(x: RealProxy, y: RealProxy2): untyped = newRealP(op(x[],y[]))
  template op*(x: RealProxy, y: RealProxy2): untyped = fn(x,y)

  template fn*(x: ImagProxy, y: ImagProxy2): untyped = newRealP(-op(x[],y[]))
  template op*(x: ImagProxy, y: ImagProxy2): untyped = fn(x,y)

  template fn*(x: RealProxy, y: ImagProxy2): untyped = newImagP(op(x[],y[]))
  template op*(x: RealProxy, y: ImagProxy2): untyped = fn(x,y)
  template fn*(x: ImagProxy, y: RealProxy2): untyped = newImagP(op(x[],y[]))
  template op*(x: ImagProxy, y: RealProxy2): untyped = fn(x,y)

  #template `fn RCU`*(x: RealProxy, y: ComplexProxy2): untyped =
  #  newComplexP(op(x[],y.re), op(x[],y.im))
  #template fn*(x: RealProxy, y: ComplexProxy2): untyped =
  #  safecall(`fn RCU`, x, y)
  #proc fn*(x: RealProxy, y: ComplexProxy2): auto {.inline,noInit.} =
  #  newComplexP(op(x[],y.re), op(x[],y.im))
  template fn*(xx: RealProxy, yy: ComplexProxy2): untyped =
    let x = xx[]
    let y = yy
    newComplexP(op(x,y.re), op(x,y.im))
  template op*(x: RealProxy, y: ComplexProxy2): untyped = fn(x,y)
  #template `fn CRU`*(x: ComplexProxy, y: RealProxy2): untyped =
  #  newComplexP(op(x.re,y[]), op(x.im,y[]))
  #template fn*(x: ComplexProxy, y: RealProxy2): untyped =
  #  safecall(`fn CRU`, x, y)
  #proc fn*(x: ComplexProxy, y: RealProxy2): auto {.inline,noInit.} =
  #  newComplexP(op(x.re,y[]), op(x.im,y[]))
  template fn*(xx: ComplexProxy, yy: RealProxy2): untyped =
    mixin op
    let x = xx
    let y = yy[]
    newComplexP(op(x.re,y), op(x.im,y))
  template op*(x: ComplexProxy, y: RealProxy2): untyped = fn(x,y)

  #template `fn ICU`*(x: ImagProxy, y: ComplexProxy2): untyped =
  #  newComplexP(-op(x[],y.im), op(x[],y.re))
  #template fn*(x: ImagProxy, y: ComplexProxy2): untyped =
  #  safecall(`fn ICU`, x, y)
  #proc fn*(x: ImagProxy, y: ComplexProxy2): auto {.inline,noInit.} =
  #  newComplexP(-op(x[],y.im), op(x[],y.re))
  template fn*(xx: ImagProxy, yy: ComplexProxy2): untyped =
    let x = xx[]
    let y = yy
    newComplexP(-op(x,y.im), op(x,y.re))
  template op*(x: ImagProxy, y: ComplexProxy2): untyped = fn(x,y)
  #template `fn CIU`*(x: ComplexProxy, y: ImagProxy2): untyped =
  #  newComplexP(-op(x.im,y[]), op(x.re,y[]))
  #template fn*(x: ComplexProxy, y: ImagProxy2): untyped =
  #  safecall(`fn CIU`, x, y)
  proc fn*(x: ComplexProxy, y: ImagProxy2): auto {.inline,noInit.} =
    newComplexP(-op(x.im,y[]), op(x.re,y[]))
  template op*(x: ComplexProxy, y: ImagProxy2): untyped = fn(x,y)

  #template `fn CCU`*(x: ComplexProxy, y: ComplexProxy2): untyped =
  #  newComplexP(op(x.re,y.re)-op(x.im,y.im),op(x.re,y.im)+op(x.im,y.re))
  #template fn*(x: ComplexProxy, y: ComplexProxy2): untyped =
  #  safecall(`fn CCU`, x, y)
  #proc fn*(x: ComplexProxy, y: ComplexProxy2): auto {.inline,noInit.} =
  #  newComplexP(op(x.re,y.re)-op(x.im,y.im),op(x.re,y.im)+op(x.im,y.re))
  template fn*(xx: ComplexProxy, yy: ComplexProxy2): untyped =
    let x = xx
    let y = yy
    newComplexP(op(x.re,y.re)-op(x.im,y.im),op(x.re,y.im)+op(x.im,y.re))
  template op*(x: ComplexProxy, y: ComplexProxy2): untyped = fn(x,y)

binaryOverloadsMul(`*`, mul)

template binaryOverloadsF(op,fn,impl: untyped) {.dirty.} =
  template fn*(x: RealProxy, y: RealProxy2): untyped = impl(x,y)
  template op*(x: RealProxy, y: RealProxy2): untyped = fn(x,y)

  template fn*(x: ImagProxy, y: ImagProxy2): untyped = impl(x,y)
  template op*(x: ImagProxy, y: ImagProxy2): untyped = fn(x,y)

  template fn*(x: RealProxy, y: ImagProxy2): untyped = impl(x,y)
  template op*(x: RealProxy, y: ImagProxy2): untyped = fn(x,y)
  template fn*(x: ImagProxy, y: RealProxy2): untyped = impl(x,y)
  template op*(x: ImagProxy, y: RealProxy2): untyped = fn(x,y)

  template fn*(x: RealProxy, y: ComplexProxy2): untyped = impl(x,y)
  template op*(x: RealProxy, y: ComplexProxy2): untyped = fn(x,y)
  template fn*(x: ComplexProxy, y: RealProxy2): untyped = impl(x,y)
  template op*(x: ComplexProxy, y: RealProxy2): untyped = fn(x,y)

  template fn*(x: ImagProxy, y: ComplexProxy2): untyped = impl(x,y)
  template op*(x: ImagProxy, y: ComplexProxy2): untyped = fn(x,y)
  template fn*(x: ComplexProxy, y: ImagProxy2): untyped = impl(x,y)
  template op*(x: ComplexProxy, y: ImagProxy2): untyped = fn(x,y)

  template fn*(x: ComplexProxy, y: ComplexProxy2): untyped = impl(x,y)
  template op*(x: ComplexProxy, y: ComplexProxy2): untyped = fn(x,y)

template divdComplexU*(x,y: untyped): untyped = x*inv(y)
binaryOverloadsF(`/`, divd, divdComplexU)

# iadd, isub, imul, idivd

template iBinaryOverloads(op,fn,impl: untyped) {.dirty.} =
  template fn*(x: ComplexProxy, y: RealProxy2) =    assign(x, impl(x,y))
  template fn*(x: ComplexProxy, y: ImagProxy2) =    assign(x, impl(x,y))
  template fn*(x: ComplexProxy, y: ComplexProxy2) = assign(x, impl(x,y))
  template op*(x: ComplexProxy, y: RealProxy2) =    fn(x,y)
  template op*(x: ComplexProxy, y: ImagProxy2) =    fn(x,y)
  template op*(x: ComplexProxy, y: ComplexProxy2) = fn(x,y)

iBinaryOverloads(`+=`, iadd, add)
iBinaryOverloads(`-=`, isub, sub)
iBinaryOverloads(`*=`, imul, mul)
iBinaryOverloads(`/=`, idivd, divd)

# inorm2, redot, iredot, dot, idot
# sqrt, rsqrt, exp, ...
# import complexFuncs
template redot*(x: ComplexProxy, y: ComplexProxy2): untyped =
  (x.adj * y).re
template redotinc*(r: RealProxy, x: ComplexProxy2, y: ComplexProxy3):
  untyped =  r += redot(x,y)


template overloadAsReal2(T: typedesc, op,fn: untyped) {.dirty.} =
  template fn*(x: T, y: RealProxy): untyped = fn(newRealProxy(x),y)
  template op*(x: T, y: RealProxy): untyped = fn(newRealProxy(x),y)
  template fn*(x: RealProxy, y: T): untyped = fn(x,newRealProxy(y))
  template op*(x: RealProxy, y: T): untyped = fn(x,newRealProxy(y))
  template fn*(x: T, y: ImagProxy): untyped = fn(newRealProxy(x),y)
  template op*(x: T, y: ImagProxy): untyped = fn(newRealProxy(x),y)
  template fn*(x: ImagProxy, y: T): untyped = fn(x,newRealProxy(y))
  template op*(x: ImagProxy, y: T): untyped = fn(x,newRealProxy(y))
  template fn*(x: T, y: ComplexProxy): untyped = fn(newRealProxy(x),y)
  template op*(x: T, y: ComplexProxy): untyped = fn(newRealProxy(x),y)
  template fn*(x: ComplexProxy, y: T): untyped = fn(x,newRealProxy(y))
  template op*(x: ComplexProxy, y: T): untyped = fn(x,newRealProxy(y))
template overloadAsReal*(T: typedesc) {.dirty.} =
  bind overloadAsReal2
  overloadAsReal2(T, `+`, add)
  overloadAsReal2(T, `-`, sub)
  overloadAsReal2(T, `*`, mul)
  overloadAsReal2(T, `/`, divd)
  overloadAsReal2(T, `:=`, assign)
  overloadAsReal2(T, `+=`, iadd)
  overloadAsReal2(T, `-=`, isub)
  overloadAsReal2(T, `*=`, imul)
  overloadAsReal2(T, `/=`, idivd)

when isMainModule:
  # complexType, complexProxy, complexFuncs
  #template neg(x: SomeNumber): untyped = -x
  template conj(x: SomeNumber): untyped = x
  template adj(x: SomeNumber): untyped = x
  template transpose(x: SomeNumber): untyped = x
  template trace(x: SomeNumber): untyped = x
  template norm2(x: SomeNumber): untyped = x*x
  template inv(x: SomeNumber): untyped = (type(x))(1)/x
  template `:=`(r: SomeNumber, x: int): untyped = r = (type(r))(x)
  template `:=`(r: SomeNumber, x: float): untyped = r = (type(r))(x)

  type
    Imag*[T] = ImagProxy[T]
    ComplexObj*[TR,TI] = object
      reX*: TR
      imX*: TI
    Complex*[TR,TI] = ComplexProxy[ComplexObj[TR,TI]]

  template newRealImpl(x: typed): untyped = x
  template newImagImpl(x: typed): untyped = newImagProxy(x)
  template newComplexImpl(x,y: typed): untyped =
    newComplexProxy(ComplexObj[type(x),type(y)](reX: x, imX: y))
  template newReal(x: typed): untyped = newRealImpl(x)
  template newImag(x: typed): untyped = newImagImpl(x)
  template newComplex(x,y: typed): untyped = newComplexImpl(x,y)

  template re*(x: ComplexObj): untyped = x.reX
  macro re*(x: ComplexObj{nkObjConstr}): auto =
    #echo x.treerepr
    result = x[1][1]
    #echo result.treerepr
  template im*(x: ComplexObj): untyped = x.imX
  macro im*(x: ComplexObj{nkObjConstr}): auto =
    #echo x.treerepr
    result = x[2][1]
    #echo result.treerepr
  template `re=`*(x: ComplexObj, y: untyped): untyped =
    x.reX := y
  template `im=`*(x: ComplexObj, y: untyped): untyped =
    x.imX := y

  overloadAsReal(SomeNumber)
  template I(x: SomeNumber): untyped = newImag(x)

  proc testadd(a,b: float) =
    var z0 = newComplex(a,b)
    var z1 = newComplex(a,b)
    var z2 = newComplex(a,b)
    var z3 = newComplex(a,b)
    var z4 = newComplex(a,b)
    var z6 = z0+z1+z2+z3+z4
    z6 := z1+z2+z3
    z6 += z2+z3+z4
    echo z6
  testadd(1.0, 2.0)

  proc testmul(a,b: float) =
    var z0 = newComplex(a,b)
    var z1 = newComplex(a+1,b+1)
    var z2 = newComplex(a+2,b+2)
    var z3 = newComplex(a,b)
    var z4 = newComplex(a,b)
    var z5 = newComplex(a,b)
    z5 += z0*z1*z2*z3*z4
    echo z5
    #var z6 = z0*z1*z2*z3*z4
    #echo z6
  testmul(1.0, 2.0)

  proc testdiv(a,b: float) =
    var z0 = newComplex(a,b)
    var z1 = newComplex(a+1,b+1)
    var z2 = newComplex(a+2,b+2)
    var z3 = newComplex(a,b)
    var z4 = newComplex(a,b)
    var z5 = newComplex(a,b)
    z5 += z0*z1*z2*z3*z4/1.0
    z5 += z0*z1*z2*z3*z4/2.0.I
    z5 += z0*z1*z2*z3*z4/z4
    z5 *= 2.0.I + 3.0
    echo z5
    #var z6 = z0*z1*z2*z3*z4
    #echo z6
  testdiv(1.0, 2.0)

  var a = 1.0 + 2.0.I
  echo a

  echo a.pos
  echo +a

  echo a.neg
  echo -a

  echo a.conj
  echo ~a

  echo a.adj
  echo @a

  echo a.transpose
  echo <>a

  echo a.trace
  echo %a

  echo a.norm2
  echo |a

  echo a.inv
  echo /a

  echo divd(1.0,a)
  echo 1.0/a

  a := 1
