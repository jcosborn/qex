import macros, math
import base/[basicOps,metaUtils]
getOptimPragmas()

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

# optimization convention:
#   flatten args if accessed more than once
#   calculate results in temporaries then return object with temps
#     newComplexP(x.re+y.re, x.im+y.im)
#   if accessing any element more than once, load it to a temporary
#     let xr = x.re
#     let xi = x.im
#     let yr = y.re
#     let yi = y.im
#     newComplexP(xr*yr-xi*yi,xr*yi+xi*yr)

template asRealProxy*[T](x: T): auto = RealProxy[typeof T](v: x)
template asRealProxy*[T](x: typedesc[T]): typedesc = RealProxy[typeof T]
template asImagProxy*[T](x: T): auto = ImagProxy[typeof T](v: x)
template asImagProxy*[T](x: typedesc[T]): typedesc = ImagProxy[typeof T]
template asComplexProxy*[T](x: T): auto = ComplexProxy[typeof T](v: x)
template asComplexProxy*[T](x: typedesc[T]): typedesc = ComplexProxy[typeof T]

template `[]`*[T](x: RealProxy[T]): auto =
  when T is ptr:
    x.v[]
  else:
    x.v

macro `[]`*[T](x: RealProxy[T]{nkObjConstr}): auto =
  when T is ptr:
    result = newBracketExpr(x[1][1])
  else:
    result = x[1][1]

template `[]`*[T](x: ImagProxy[T]): auto =
  when T is ptr:
    x.v[]
  else:
    x.v
macro `[]`*[T](x: ImagProxy[T]{nkObjConstr}): auto =
  when T is ptr:
    result = newBracketExpr(x[1][1])
  else:
    result = x[1][1]

template isWrapper*(x: RealProxy): auto = true
template asWrapper*(x: RealProxy, y: typed): auto =
  asRealProxy(y)
template asVarWrapper*(x: RealProxy, y: typed): auto =
  asVar(asRealProxy(y))
template isWrapper*(x: ImagProxy): auto = true
template asWrapper*(x: ImagProxy, y: typed): auto =
  asImagProxy(y)
template asVarWrapper*(x: ImagProxy, y: typed): auto =
  asVar(asImagProxy(y))
template isWrapper*(x: ComplexProxy): bool = true
template isWrapper*(x: typedesc[ComplexProxy]): bool = true
template asWrapper*(x: ComplexProxy, y: typed): auto =
  asComplexProxy(y)
template asWrapper*(x: typedesc[ComplexProxy], y: typed): auto =
  asComplexProxy(y)
template asVarWrapper*(x: ComplexProxy, y: typed): auto =
  asVar(asComplexProxy(y))

template has*[T:ComplexProxy](x: typedesc[T], y: typedesc): bool =
  mixin has
  when y is ComplexProxy: true
  else: has(T.type[], y)

template `[]`*[T](x: typedesc[ComplexProxy[T]]): typedesc =
  when T is ptr:
    T.type[]
  else:
    T.type

template `[]`*[T](x: ComplexProxy[T]): auto =
  when T is ptr:
    x.v[]
  else:
    x.v
macro `[]`*[T](x: ComplexProxy[T]{nkObjConstr}): auto =
  #echo x.treerepr
  when T is ptr:
    result = newBracketExpr(x[1][1])
  else:
    result = x[1][1]
  #echo result.treerepr
macro `[]`*[T](x: ComplexProxy[T]{nkStmtListExpr}): auto =
  #echo x.treerepr
  result = newNimNode(nnkStmtListExpr)
  for i in 0..(x.len-2):
    result.add x[i]
  result.add newCall(ident"[]", x[^1])

template `[]=`*(x: RealProxy, y: typed) =
  #when x.T is type(y):
  #  x.v = y
  #else:
  mixin `:=`
  x.v := y
template `[]=`*(x: ImagProxy, y: typed) =
  #when x.T is type(y):
  #  x.v = y
  #else:
    x.v := y
template `[]=`*(x: ComplexProxy, y: typed) =
  #when x.T is type(y):
  #  x.v = y
  #else:
    x.v := y
proc `[]=`*(x: var ComplexProxy, i: auto, y: ComplexProxy2) {.alwaysInline.} =
  mixin re, im, `[]=`
  when isWrapper i:
    x.re[i] = y.re
    x.im[i] = y.im
  else:
    {.error.}

template eval*[T](x: typedesc[ComplexProxy[T]]): typedesc =
  mixin eval
  asComplexProxy(eval typeof T)

proc `$`*(x: RealProxy): string =
  result = $x[]
proc `$`*(x: ImagProxy): string =
  result = $x[] & "I"
proc `$`*(x: ComplexProxy): string =
  mixin re, im
  #result = "(" & $x.reX & "," & $x.imX & ")"
  result = $x[].re
  var t = $x[].im & "I"
  if t[0]!='-': result.add "+"
  result.add t

proc `|`*(x: ComplexProxy, y: tuple): string =
  result = x[].re | y
  var t = (x[].im | y) & "I"
  if t[0]!='-': result.add "+"
  result.add t


template newRealProxy*[T](x: T): untyped =
  RealProxy[type(T)](v: x)
#template newRealProxy*(x: typed{call}): untyped =
#  let t = x
#  RealProxy[type(t)](v: t)

template newImagProxy*[T](x: T): untyped =
  ImagProxy[type(T)](v: x)
#template newImagProxy*(x: typed{call}): untyped =
#  let t = x
#  ImagProxy[type(t)](v: t)

#template newComplexProxyU*[T](x: T): untyped =
#  ComplexProxy[type(T)](v: x)
#template newComplexProxy*(x: typed): untyped =
#  flattenCallArgs(newComplexProxyU, x)
template newComplexProxy*[T](x: typedesc[T]): typedesc =
  ComplexProxy[type(T)]
template newComplexProxy*[T](x: T): untyped =
  ComplexProxy[type(T)](v: x)
#template newComplexProxy*(x: typed{call}): untyped =
#  let t = x
#  ComplexProxy[type(t)](v: t)

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

template `re=`*(x: RealProxy, y: typed) = x[] = y
template `im=`*(x: RealProxy, y: typed) = discard
template `re=`*(x: ImagProxy, y: typed) = discard
template `im=`*(x: ImagProxy, y: typed) = x[] = y
template `re=`*(x: ComplexProxy, y: typed) = x[].re = y
template `im=`*(x: ComplexProxy, y: typed) = x[].im = y

#template setU*(r: ComplexProxy, x: typed, y: typed) =
#  r[].re = x
#  r[].im = y
#template set*(r: ComplexProxy, x: typed, y: typed) =
#  flattenCallArgs(setU, r, x, y)
proc set*(r: var ComplexProxy, x: auto, y: auto) {.alwaysInline.} =
  r.re = x
  r.im = y

template assign*(x: RealProxy, y: RealProxy2): untyped =
  x[] = y[]
template assign*(x: ImagProxy, y: ImagProxy2): untyped =
  x[] = y[]
template assignU*(x: ComplexProxy, y: RealProxy2): untyped =
  x[].re = y[]
  x[].im = 0
#template assign*(x: ComplexProxy, y: RealProxy2): untyped =
#  #echoRepr: x
#  let assignCR = y
#  #assignU(x, yy)
#  flattenCallArgs(assignU, x, assignCR)
proc assign*(x: var ComplexProxy, y: RealProxy2) {.alwaysInline.} =
  x.re = y[]
  x.im = 0
template assign*(x: ComplexProxy, y: ImagProxy2): untyped =
  x[].re = 0
  x[].im = y[]
template assignU*(x: ComplexProxy, y: ComplexProxy2): untyped =
  #echoRepr: x
  #static: echo "cp assignU"
  x.re = y.re
  x.im = y.im
#template assignU*(x: ComplexProxy, yy: ComplexProxy2{call}): untyped =
#  let ya {.noInit.} = yy
#  #var ya {.noInit.} = yy
#  #var ya {.noInit.}: T2
#  #ya = yy
#  assignU(x, ya)
proc assign*[R,X:ComplexProxy](r: var R, x: X) {.alwaysInline.} =
#template assign*[R,X:ComplexProxy](rr: R, xx: X) =
  #echoRepr: x
  #let assignCC = xx
  #assignU(x, yy)
  #flattenCallArgs(assignU, r, assignCC)
  #let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
  #let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
  r.re = x.re
  r.im = x.im

template `:=`*(x: RealProxy, y: RealProxy2): untyped = assign(x,y)
template `:=`*(x: ImagProxy, y: ImagProxy2): untyped = assign(x,y)
template `:=`*(x: ComplexProxy, y: RealProxy2): untyped = assign(x,y)
template `:=`*(x: ComplexProxy, y: ImagProxy2): untyped = assign(x,y)
template `:=`*(x: ComplexProxy, y: ComplexProxy2): untyped = assign(x,y)

# pos, neg, conj, adj, transpose, trace
template unaryOverloads(op,fn,implR,implI: untyped) {.dirty.} =
  template fn*(x: RealProxy): untyped = newRealP(implR(x[],0))
  #proc fn*(x: RealProxy): auto {.inline,noInit.} = newRealP(implR(x[],0))
  template op*(x: RealProxy): untyped = fn(x)

  template fn*(x: ImagProxy): untyped = newImagP(implI(0,x[]))
  #proc fn*(x: ImagProxy): auto {.inline,noInit.} = newImagP(implI(0,x[]))
  template op*(x: ImagProxy): untyped = fn(x)

  template `fn U`*(x: ComplexProxy): untyped =
    newComplexP(implR(x.re,x.im), implI(x.re,x.im))
  #template fn*(x: ComplexProxy): untyped =
  #  flattenCallArgs(`fn U`, x)
  proc fn*(x: ComplexProxy): auto {.alwaysInline,noInit.} =
    newComplexP(implR(x.re,x.im), implI(x.re,x.im))
  #template fn*(xx: ComplexProxy): untyped =
  #  let x = xx
  #  newComplexP(implR(x.re,x.im), implI(x.re,x.im))
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

#template adjRealU(xr,xi: untyped): untyped =
#  mixin adj
#  adj(xr)
#template adjImagU(xr,xi: untyped): untyped =
#  mixin adj
#  -adj(xi)
#unaryOverloads(`@`, adj, adjRealU, adjImagU)

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

# norm2, abs
template unaryOverloadsR(op,fn,implR,implI: untyped) {.dirty.} =
  template fn*(x: RealProxy): untyped = newRealP(implR(x[],0))
  #proc fn*(x: RealProxy): auto {.inline,noInit.} = newRealP(implR(x[],0))
  template op*(x: RealProxy): untyped = fn(x)

  template fn*(x: ImagProxy): untyped = newRealP(implI(0,x[]))
  #proc fn*(x: ImagProxy): auto {.inline,noInit.} = newRealP(implI(0,x[]))
  template op*(x: ImagProxy): untyped = fn(x)

  template `fn U`*(x: ComplexProxy): untyped =
    newRealP(implR(x.re,x.im)+implI(x.re,x.im))
  #template fn*(x: ComplexProxy): untyped =
  #  flattenCallArgs(`fn U`, x)
  proc fn*(x: ComplexProxy): auto {.alwaysinline,noInit.} =
    newRealP(implR(x.re,x.im)+implI(x.re,x.im))
  template op*(x: ComplexProxy): untyped = fn(x)

template norm2RealU(xr,xi: untyped): untyped =
  mixin norm2
  norm2(xr)
template norm2ImagU(xr,xi: untyped): untyped =
  mixin norm2
  norm2(xi)
unaryOverloadsR(`|`, norm2, norm2RealU, norm2ImagU)

template abs*(x: ComplexProxy): untyped =
  sqrt(norm2(x))


#template inv*(x: RealProxy): untyped = newRealP(x[].inv)
proc inv*(x: RealProxy): auto {.alwaysInline,noInit.} =
  mixin inv
  newRealP(x[].inv)
template `/`*(x: RealProxy): untyped = inv(x)
#template inv*(x: ImagProxy): untyped = newImagP(-x[].inv)
proc inv*(x: ImagProxy): auto {.alwaysInline,noInit.} =
  mixin inv
  newImagP(-x[].inv)
template `/`*(x: ImagProxy): untyped = inv(x)
#template invComplexU(x: untyped): untyped = x.adj * x.norm2.inv
#template inv*(x: ComplexProxy): untyped = safecall(invComplexU, x)
proc inv*(x: ComplexProxy): auto {.alwaysInline,noInit.} =
  mixin inv
  x.adj * x.norm2.inv
template `/`*(x: ComplexProxy): untyped = inv(x)

proc exp*(x: ComplexProxy): auto {.inline,noInit.} =
  mixin exp, cos, sin
  let er = exp(x.re)
  let xi = x.im
  let ci = cos(xi)
  let si = sin(xi)
  newComplexP(er*ci, er*si)

proc ln*(x: ComplexProxy): auto {.inline,noInit.} =
  mixin ln, atan2
  let n = 0.5 * x.norm2.ln
  let xr = x.re
  let xi = x.im
  newComplexP(n, atan2(xi,xr))

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

  template `fn RCU`*(x: RealProxy, y: ComplexProxy2): untyped =
    newComplexP(op(x[],y.re), op(y.im))
  template fn*(x: RealProxy, y: ComplexProxy2): untyped =
    flattenCallArgs(`fn RCU`, x, y)
  #proc fn*(x: RealProxy, y: ComplexProxy2): auto {.inline,noInit.} =
  #  newComplexP(op(x[],y.re), op(y.im))
  #template fn*(x: RealProxy, yy: ComplexProxy2): untyped =
  #  let y = yy
  #  newComplexP(op(x[],y.re), op(y.im))
  template op*(x: RealProxy, y: ComplexProxy2): untyped = fn(x,y)
  template `fn CRU`*(x: ComplexProxy, y: RealProxy2): untyped =
    newComplexP(op(x.re,y[]), x.im)
  template fn*(x: ComplexProxy, y: RealProxy2): untyped =
    flattenCallArgs(`fn CRU`, x, y)
  #proc fn*(x: ComplexProxy, y: RealProxy2): auto {.inline,noInit.} =
  #  newComplexP(op(x.re,y[]), x.im)
  template op*(x: ComplexProxy, y: RealProxy2): untyped = fn(x,y)

  template `fn ICU`*(x: ImagProxy, y: ComplexProxy2): untyped =
    newComplexP(op(y.re), op(x.im,y.im))
  template fn*(x: ImagProxy, y: ComplexProxy2): untyped =
    flattenCallArgs(`fn ICU`, x, y)
  #proc fn*(x: ImagProxy, y: ComplexProxy2): auto {.inline,noInit.} =
  #  newComplexP(op(y.re), op(x.im,y.im))
  template op*(x: ImagProxy, y: ComplexProxy2): untyped = fn(x,y)
  template `fn CIU`*(x: ComplexProxy, y: ImagProxy2): untyped =
    newComplexP(x.re, op(x.im,y.im))
  template fn*(x: ComplexProxy, y: ImagProxy2): untyped =
    flattenCallArgs(`fn CIU`, x, y)
  #proc fn*(x: ComplexProxy, y: ImagProxy2): auto {.inline,noInit.} =
  #  newComplexP(x.re, op(x.im,y[]))
  template op*(x: ComplexProxy, y: ImagProxy2): untyped = fn(x,y)

  template `fn U`*(x: ComplexProxy, y: ComplexProxy2): untyped =
    newComplexP(op(x.re,y.re), op(x.im,y.im))
  template fn*(x: ComplexProxy, y: ComplexProxy2): untyped =
    flattenCallArgs(`fn U`, x, y)
  #proc fn*(x: ComplexProxy, y: ComplexProxy2): auto {.inline,noInit.} =
  #  newComplexP(op(x.re,y.re), op(x.im,y.im))
  #template fn*(xx: ComplexProxy, yy: ComplexProxy2): untyped =
  #  let x = xx
  #  let y = yy
  #  newComplexP(op(x.re,y.re), op(x.im,y.im))
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

  #template `fn U`*(x: RealProxy, y: ComplexProxy2): untyped =
  #  let `t fn R` = op(x[], y.re)
  #  let `t fn I` = op(x[], y.im)
  #  newComplexP(`t fn R`, `t fn I`)
  #template fn*(x: RealProxy, y: ComplexProxy2): untyped =
  #  flattenCallArgs(`fn U`, x, y)
  #proc fn*(x: RealProxy, y: ComplexProxy2): auto {.inline,noInit.} =
  #  newComplexP(op(x[],y.re), op(x[],y.im))
  template `fn RCU`*(xx: RealProxy, y: ComplexProxy2): untyped =
    let x = xx[]
    newComplexP(op(x,y.re), op(x,y.im))
  template fn*(x: RealProxy, y: ComplexProxy2): untyped =
    flattenCallArgs(`fn RCU`, x, y)
  template op*(x: RealProxy, y: ComplexProxy2): untyped = fn(x,y)

  template `fn CRU`*(x: ComplexProxy, yy: RealProxy2): untyped =
    mixin op
    let y = yy[]
    newComplexP(op(x.re,y), op(x.im,y))
  template fn*(x: ComplexProxy, y: RealProxy2): untyped =
    flattenCallArgs(`fn CRU`, x, y)
  #proc fn*(x: ComplexProxy, y: RealProxy2): auto {.inline,noInit.} =
  #  newComplexP(op(x.re,y[]), op(x.im,y[]))
  #template fn*(xx: ComplexProxy, yy: RealProxy2): untyped =
  #  mixin op
  #  let x = xx
  #  let y = yy[]
  #  newComplexP(op(x.re,y), op(x.im,y))
  template op*(x: ComplexProxy, y: RealProxy2): untyped = fn(x,y)

  template `fn ICU`*(xx: ImagProxy, y: ComplexProxy2): untyped =
    let x = xx[]
    newComplexP(-op(x,y.im), op(x,y.re))
  template fn*(x: ImagProxy, y: ComplexProxy2): untyped =
    flattenCallArgs(`fn ICU`, x, y)
  #proc fn*(x: ImagProxy, y: ComplexProxy2): auto {.inline,noInit.} =
  #  newComplexP(-op(x[],y.im), op(x[],y.re))
  #template fn*(xx: ImagProxy, yy: ComplexProxy2): untyped =
  #  let x = xx[]
  #  let y = yy
  #  newComplexP(-op(x,y.im), op(x,y.re))
  template op*(x: ImagProxy, y: ComplexProxy2): untyped = fn(x,y)

  template `fn CIU`*(x: ComplexProxy, yy: ImagProxy2): untyped =
    let y = yy[]
    newComplexP(-op(x.im,y), op(x.re,y))
  template fn*(x: ComplexProxy, y: ImagProxy2): untyped =
    flattenCallArgs(`fn CIU`, x, y)
  #proc fn*(x: ComplexProxy, y: ImagProxy2): auto {.inline,noInit.} =
  #  newComplexP(-op(x.im,y[]), op(x.re,y[]))
  template op*(x: ComplexProxy, y: ImagProxy2): untyped = fn(x,y)

  template `fn CCU`*(x: ComplexProxy, y: ComplexProxy2): untyped =
    #let `t fn R` = op(x.re,y.re)-op(x.im,y.im)
    #let `t fn I` = op(x.re,y.im)+op(x.im,y.re)
    #newComplexP(`t fn R`, `t fn I`)
    let xr = x.re
    let xi = x.im
    let yr = y.re
    let yi = y.im
    newComplexP(op(xr,yr)-op(xi,yi),op(xr,yi)+op(xi,yr))
  template fn*(x: ComplexProxy, yy: ComplexProxy2): untyped =
    let y = yy
    flattenCallArgs(`fn CCU`, x, y)
  #proc fn*(x: ComplexProxy, y: ComplexProxy2): auto {.inline,noInit.} =
  #  newComplexP(op(x.re,y.re)-op(x.im,y.im),op(x.re,y.im)+op(x.im,y.re))
  #template fn*(xx: ComplexProxy, yy: ComplexProxy2): untyped =
  #  let x = xx
  #  let y = yy
  #  newComplexP(op(x.re,y.re)-op(x.im,y.im),op(x.re,y.im)+op(x.im,y.re))
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

template divdComplexU*(x,y: untyped): untyped =
  let iy = inv(y)
  x*iy
binaryOverloadsF(`/`, divd, divdComplexU)

# iadd, isub, imul, idivd

template iBinaryOverloads(op,fn,impl: untyped) {.dirty.} =
  #template `fn U`*(x: ComplexProxy, y: RealProxy2) = assign(x, impl(x,y))
  #template fn*(x: ComplexProxy, y: RealProxy2) =
  #  flattenCallArgs(`fn U`, x, y)
  proc fn*(x: var ComplexProxy, y: RealProxy2) {.alwaysInline.} =
    assign(x, impl(x,y))
  template fn*(x: ComplexProxy, y: ImagProxy2) =    assign(x, impl(x,y))
  #template `fn U`*(x: ComplexProxy, y: ComplexProxy2) = assign(x, impl(x,y))
  #template fn*(x: ComplexProxy, y: ComplexProxy2) =
  #  flattenCallArgs(`fn U`, x, y)
  proc fn*(x: var ComplexProxy, y: ComplexProxy2) {.alwaysInline.} =
    assign(x, impl(x,y))
  template op*(x: ComplexProxy, y: RealProxy2) =    fn(x,y)
  template op*(x: ComplexProxy, y: ImagProxy2) =    fn(x,y)
  template op*(x: ComplexProxy, y: ComplexProxy2) = fn(x,y)
  template `fn U`*(x: ImagProxy, y: ImagProxy2) = assign(x, impl(x,y))
  template fn*(x: ImagProxy, y: ImagProxy2) =
    flattenCallArgs(`fn U`, x, y)
  template op*(x: ImagProxy, y: ImagProxy2) =    fn(x,y)

iBinaryOverloads(`+=`, iadd, add)
iBinaryOverloads(`-=`, isub, sub)
iBinaryOverloads(`*=`, imul, mul)
iBinaryOverloads(`/=`, idivd, divd)

#[
template iadd*[R,X:ComplexProxy](r: R, xx: X) =
  let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
  iadd(r.re, x.re)
  iadd(r.im, x.im)
template `+=`*[R,X:ComplexProxy](r: R, x: X) = iadd(r,x)
]#

proc sqrt*(x: ComplexProxy): auto =
  mixin copySign
  let n = sqrt(x.norm2)
  let r = sqrt(0.5*(n + x.re))
  #let i = select(x.im<0, -1, 1)*sqrt(0.5*(n - x.re))
  let i = copySign(sqrt(0.5*(n - x.re)), x.im)
  newComplexP(r, i)

# inorm2, redot, iredot, dot, idot
# sqrt, rsqrt, exp, ...
# import complexFuncs
template redot*(x: ComplexProxy, y: ComplexProxy2): untyped =
  #(x.adj * y).re
  # below needed to workaround C++ backend issue (duplicate variable name)
  let xx = x
  let yy = y
  xx.re*yy.re + xx.im*yy.im

template redotinc*(r: RealProxy, x: ComplexProxy2,
                   y: ComplexProxy3): untyped =
  r += redot(x,y)

template dot*(x: ComplexProxy, y: ComplexProxy2): untyped =
  x.adj * y


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
