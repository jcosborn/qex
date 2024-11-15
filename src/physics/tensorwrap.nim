import base/basicOps
import base/wrapperTypes  # for getPtr
export wrapperTypes
import maths
import simd/simdWrap
getOptimPragmas()

type
  TensorObj*[K,T] = object
    obj*: T
  SomeTensor*[K,T] = TensorObj[K,T]
  SomeTensor2*[K,T] = TensorObj[K,T]
  SomeTensor3*[K,T] = TensorObj[K,T]
  SomeTensor4*[K,T] = TensorObj[K,T]

template tensorObj*[K,T](k: typedesc[K], t: typedesc[T]): typedesc =
  TensorObj[K,T]
template tensorObj*[K,T](k: typedesc[K], x: T): auto =
  TensorObj[K,T](obj: x)

template kind*[K,T](x: typedesc[SomeTensor[K,T]]): typedesc = K
template kind*[K,T](x: SomeTensor[K,T]): typedesc = K

template `[]`*[K,T](x: typedesc[SomeTensor[K,T]]): typedesc = T
template `[]`*[K,T](x: TensorObj[K,T]): auto = x.obj
#proc `[]`*[K,T](x: var TensorObj[K,T]): var T {.alwaysInline.} = x.obj

template isWrapper*[K,T](x: typedesc[SomeTensor[K,T]]): bool = true
template isWrapper*[K,T](x: SomeTensor[K,T]): bool = true
template asWrapper*[K,T,Y](x: typedesc[SomeTensor[K,T]], y: typedesc[Y]): typedesc =
  tensorObj(K, Y)
template asWrapper*[K,T](x: typedesc[SomeTensor[K,T]], y: auto): auto =
  tensorObj(K, y)
template asWrapper*[K,T](x: SomeTensor[K,T], y: auto): auto =
  tensorObj(K, y)

template eval*[K,T](x: typedesc[TensorObj[K,T]]): typedesc =
  tensorObj(K,eval(T))

template has*[K,T](x: typedesc[SomeTensor[K,T]], y: typedesc): bool =
#template has*[K,T;Y:typedesc](x: typedesc[SomeTensor[K,T]], y: Y): bool =
  mixin has
  when y is SomeTensor[K,auto]: true
  else: has(T[], y)

template index*[K,T,I](x: typedesc[TensorObj[K,T]], i: typedesc[I]): typedesc =
  when I is SomeTensor[K,auto]:
    index(K, I[])
  elif I.isWrapper:
    tensorObj(K, index(T, I))
  else:
    index(T.type, I.type)

template index*[K,T](x: typedesc[TensorObj[K,T]], i,j: typedesc[int]): typedesc =
  index(T.type, int, int)

template `[]`*[K,T,I](x: TensorObj[K,T], i: I): auto =
  when I is SomeTensor[K,auto]:
    x[][i[]]
  elif I.isWrapper:
    #indexed(x, i)
    var tTensorObjBracket = tensorObj(K, x[][i])
    tTensorObjBracket  # need to return var type
  else:
    x[][i]

template `[]=`*[K,T,I](x: SomeTensor[K,T], i: I, y: auto) =
  mixin `:=`
  when I is SomeTensor[K,auto]:
    x[][i[]] = y
  elif y is SomeTensor[K,auto]:
    x[][i] = y[]
  else:
    x[][i] = y

template `[]`*[K,T,I](x: TensorObj[K,T], i,j: I): auto =
  x[][i,j]
template `[]=`*[K,T,I](x: TensorObj[K,T], i,j: I, y: auto) =
  x[][i,j] = y

# forward from value to value
template forwardVV(f: untyped) {.dirty.} =
  template f*(x: SomeTensor): auto =
    mixin f
    f(x[])
# forward from value to type
#template forwardVT(f: untyped) {.dirty.} =
#  template f*[T](x: SomeTensor[T]): auto =
#    mixin f
#    f(type T)
# forward from type to type
template forwardTT(f: untyped) {.dirty.} =
  template f*[K,T](x: typedesc[SomeTensor[K,T]]): auto =
    mixin f
    f(type T)
# forward from type to type and wrap
template forwardTTW(f: untyped) {.dirty.} =
  template f*[K,T](x: typedesc[TensorObj[K,T]]): auto =
    mixin f
    tensorObj(K, f(type T))

forwardVV(len)
forwardVV(nrows)
forwardVV(ncols)
forwardVV(nVectors)
forwardVV(simdType)
#forwardVV(simdLength)
forwardVV(getNs)
forwardVV(numberType)
forwardVV(numNumbers)

forwardTT(len)
forwardTT(nrows)
forwardTT(ncols)
forwardTT(nVectors)
forwardTT(simdType)
forwardTT(simdLength)
forwardTT(getNs)
forwardTT(numberType)

forwardTTW(toSingle)
forwardTTW(toDouble)

#** Unary assignment functions
template setUnaryAssignT(f) {.dirty.} =
  template f*[R:SomeTensor,X:SomeTensor2](r: var R, x: X) =
    when R.kind is X.kind:
      f(r[], x[])
    else:  # assume R has X
      f(r[], x)
template setUnaryAssignX(f,X) {.dirty.} =
  template f*[R:SomeTensor](r: var R, x: X) =
    f(r[], x)

setUnaryAssignT(assign)
setUnaryAssignX(assign, SomeNumber)
setUnaryAssignX(assign, AsComplex)

setUnaryAssignT(`:=`)
setUnaryAssignX(`:=`, SomeNumber)
setUnaryAssignX(`:=`, AsComplex)

setUnaryAssignT(`+=`)
setUnaryAssignX(`+=`, SomeNumber)
setUnaryAssignX(`+=`, AsComplex)

setUnaryAssignT(`-=`)
setUnaryAssignX(`-=`, SomeNumber)
setUnaryAssignX(`-=`, AsComplex)

setUnaryAssignT(`*=`)
setUnaryAssignX(`*=`, SomeNumber)
setUnaryAssignX(`*=`, AsComplex)

setUnaryAssignT(iadd)
setUnaryAssignX(iadd, SomeNumber)
setUnaryAssignX(iadd, AsComplex)

setUnaryAssignT(isub)
setUnaryAssignX(isub, SomeNumber)
setUnaryAssignX(isub, AsComplex)

setUnaryAssignT(imul)
setUnaryAssignX(imul, SomeNumber)
setUnaryAssignX(imul, AsComplex)

#** Unary functions
template setUnopT(f, g) {.dirty.} =
  template f*[X:SomeTensor](x: typedesc[X]): typedesc =
    tensorObj(X.kind, f(X[]))
  template f*[X:SomeTensor](x: X): auto =
    tensorObj(X.kind, f(x[]))

setUnopT(load1, load1)
setUnopT(`-`, neg)
setUnopT(re, re)
setUnopT(im, im)
setUnopT(simdSum, simdSum)

#** Binary assignment functions

template setBinopAssignTT(f) {.dirty.} =
  template f*[R:SomeTensor,X:SomeTensor2,Y:SomeTensor3](r: var R, x: X, y: Y) =
    when R.kind is X.kind:
      when R.kind is Y.kind:
        f(r[], x[], y[])
      else:
        f(r[], x[], y)
    else:
      when R.kind is Y.kind:
        f(r[], x, y[])
      else:
        f(r[], x, y)
template setBinopAssignTX(f,Y) {.dirty.} =
  template f*[R:SomeTensor,X:SomeTensor2](r: var R, x: X, y: Y) =
    when R.kind is X.kind:
      f(r[], x[], y)
    else:
      f(r[], x, y)
template setBinopAssignXT(f,X) {.dirty.} =
  template f*[R:SomeTensor,Y:SomeTensor3](r: var R, x: X, y: Y) =
    when R.kind is Y.kind:
      f(r[], x, y[])
    else:
      f(r[], x, y)

setBinopAssignTT(add)
setBinopAssignTT(sub)

setBinopAssignTT(mul)
setBinopAssignTX(mul, SomeNumber)
setBinopAssignXT(mul, SomeNumber)
setBinopAssignTX(mul, AsComplex)
setBinopAssignXT(mul, AsComplex)

setBinopAssignTT(imadd)
setBinopAssignTX(imadd, AsComplex)
setBinopAssignXT(imadd, AsComplex)

setBinopAssignTT(imsub)

setBinopAssignTT(peqOuter)
setBinopAssignTT(meqOuter)

#** Binary functions

template setBinopTT(f, g) {.dirty.} =
  template f*[X:SomeTensor,Y:SomeTensor2](x: typedesc[X], y: typedesc[Y]): typedesc =
    when X.kind is Y.kind:
      tensorObj(X.kind, f(X[],Y[]))
    elif X.kind.has Y.kind:
      f(X[], Y)
    else:  # assume Y.kind has X.kind
      f(X, Y[])
  template f*[X:SomeTensor,Y:SomeTensor2](x: X, y: Y): auto =
    when X.kind is Y.kind:
      tensorObj(X.kind, f(x[],y[]))
    elif X.kind.has Y.kind:
      f(x[], y)
    else:  # assume Y.kind has X.kind
      f(x, y[])
  #template f*(x: SomeTensor, y: SomeTensor2): auto =
  #  var r {.noInit.}: f(x.type, y.type)
  #  g(r, x, y)
  #  r
template setBinopTX(f,g,Y) {.dirty.} =
  template f*[X:SomeTensor](x: typedesc[X], y: typedesc[Y]): typedesc =
    tensorObj(X.kind, f(X[], Y))
  template f*[X:SomeTensor](x: X, y: Y): auto =
    tensorObj(X.kind, f(x[], y))
  #template f*(x: SomeTensor, y: X): auto =
  #  var r {.noInit.}: f(x.type, X)
  #  g(r, x, y)
  #  r
template setBinopXT(f,g,X) {.dirty.} =
  template f*[Y:SomeTensor](x: typedesc[X], y: typedesc[Y]): typedesc =
    tensorObj(Y.kind, f(X, Y[]))
  template f*[Y:SomeTensor](x: X, y: Y): auto =
    tensorObj(Y.kind, f(x, y[]))
  #template f*[Y:SomeTensor](x: X, y: Y): auto =
  #  var r {.noInit.}: f(X, Y)
  #  g(r, x, y)
  #  r

setBinopTT(`+`, add)
setBinopTX(`+`, add, SomeNumber)
setBinopXT(`+`, add, SomeNumber)

setBinopTT(`-`, sub)
setBinopTX(`-`, sub, SomeNumber)
setBinopXT(`-`, sub, SomeNumber)

setBinopTT(`*`, mul)
setBinopTX(`*`, mul, SomeNumber)
setBinopXT(`*`, mul, SomeNumber)
setBinopTX(`*`, mul, Simd)
setBinopXT(`*`, mul, Simd)
setBinopTX(`*`, mul, AsReal)
setBinopXT(`*`, mul, AsReal)
setBinopTX(`*`, mul, AsImag)
setBinopXT(`*`, mul, AsImag)
setBinopTX(`*`, mul, AsComplex)
setBinopXT(`*`, mul, AsComplex)

#[
#template mul*(x: SomeNumber, y: SomeTensor2): auto =
#  asSomeTensor(`*`(x, y[]))
template mul*[X:SomeNumber,Y:SomeTensor](x: X, y: Y): auto =
  #var tmp {.noInit.}: asSomeTensor(type(X)*type(Y)[])
  var tmp {.noInit.}: asSomeTensor(evalType(x*y[]))
  mul(tmp[], x, y[])
  tmp
template mul*(x: SomeTensor, y: SomeTensor2): auto =
  asSomeTensor(`*`(x[], y[]))
template mul*(x: AsComplex, y: SomeTensor2): auto =
  asSomeTensor(mul(x, y[]))
]#

#template row*[K,T](x: SomeTensor[K,T], i: auto): auto =
#  mixin row
#  tensorObj(x.kind, row(x[],i))
#template setRow*[K,T1,T2](r: SomeTensor[K,T1]; x: SomeTensor2[K,T2]; i: auto) =
#  setRow(r[], x[], i)
setBinopAssignTX(setRow, auto)
setBinopTX(row, setRow, auto)

#template random*(x: var SomeTensor) = gaussian(x[], r)
setUnaryAssignX(gaussian, auto)
setUnaryAssignX(uniform, auto)
setUnaryAssignX(z2, auto)
setUnaryAssignX(z4, auto)
setUnaryAssignX(u1, auto)

template projectU*(r: var SomeTensor) =  # in place
  projectU(r[])
#template projectU*(r: var SomeTensor, x: SomeTensor2) =
#  projectU(r[], x[])
setUnaryAssignT(projectU)

template projectUderiv*(r: var SomeTensor, u: SomeTensor2,
                        x: SomeTensor3, chain: SomeTensor4) =
  projectUderiv(r[], u[], x[], chain[])
template projectUderiv*(r: var SomeTensor, x: SomeTensor3, chain: SomeTensor4) =
  projectUderiv(r[], x[], chain[])
template projectSU*(r: var SomeTensor) =
  projectSU(r[])
template projectSU*(r: var SomeTensor, x: SomeTensor2) =
  projectSU(r[], x[])
template projectTAH*(r: var SomeTensor) =
  projectTAH(r[])
template projectTAH*(r: var SomeTensor, x: SomeTensor2) =
  projectTAH(r[], x[])
template checkU*(x: SomeTensor): auto = checkU(x[])
template checkSU*(x: SomeTensor): auto = checkSU(x[])

template norm2*(x: SomeTensor): auto = norm2(x[])
template norm2*(r: var auto, x: SomeTensor) = norm2(r, x[])
template inorm2*(r: var auto, x: SomeTensor) = inorm2(r, x[])
template dot*(x: SomeTensor, y: SomeTensor2): auto = dot(x[], y[])
template idot*(r: var auto, x: SomeTensor2, y: SomeTensor3) = idot(r, x[], y[])
template redot*(x: SomeTensor, y: SomeTensor2): auto = redot(x[], y[])
template trace*(x: SomeTensor): auto = trace(x[])

template exp*[X:SomeTensor](x: X): auto =
  tensorObj(X.kind, exp(x[]))
template expDeriv*[X:SomeTensor,C:SomeTensor2](x: X, c: C): auto =
  tensorObj(X.kind, expDeriv(x[], c[]))
template ln*[X:SomeTensor](x: X): auto =
  tensorObj(X.kind, ln(x[]))
