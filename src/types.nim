import stdUtils
import comms
import complexConcept
import matrixConcept
import metaUtils
import wrapperTypes
import macros

# converter adjMat*[T](x:Adjointed[AsMatrix[T]]):AsMatrix[Adjointed[T]]

type
  #Adjointed*{.borrow: `.`.}[T] = distinct T
  #Adjointed*[T] = distinct T
  Adjointed*[T] = object
    v*:T
  #Adjointed*[T] = object
  #  v*:ptr T
template adj*(xx:typed):expr =
  subst(x,xx):
    when isComplex(x):
      asComplex(adj(x[]))
    elif isVector(x):
      asVector(adj(x[]))
    elif isMatrix(x):
      asMatrix(adj(x[]))
    else:
      #when compiles(addr(x)):
      when compiles(unsafeAddr(x)):
        cast[ptr Adjointed[type(x)]](unsafeAddr(x))[]
        #cast[Adjointed[type(x)]](unsafeAddr(x))
      else:
        #(Adjointed[type(x)])(x)
        cast[Adjointed[type(x)]](x)
        #cast[Adjointed[type(x)]]((var t=x; t.addr))
#template `[]`*[T](x:Adjointed[T]):expr = cast[T](x)
makeDeref(Adjointed, x.T)
template `[]`*(x:Adjointed; i:SomeInteger):expr = x[][i].adj
template `[]`*(x:Adjointed; i,j:SomeInteger):expr = x[][j,i].adj
template len*(x:Adjointed):expr = x[].len
template nrows*(x:Adjointed):expr = x[].ncols
template ncols*(x:Adjointed):expr = x[].nrows
template isVector*(x:Adjointed):expr = isVector(x[])
template isMatrix*(x:Adjointed):expr = isMatrix(x[])
template re*(x:Adjointed):expr = x[].re
template im*(x:Adjointed):expr = -(x[].im)
#template mvLevel*(x:Adjointed):expr =
#  mixin mvLevel
#  mvLevel(x[])


type
  MaskedObj*[T] = object
    pobj*:ptr T
    mask*:int
  Masked*[T] = distinct MaskedObj[T]
template pobj*(x:Masked):expr = ((MaskedObj[x.T])(x)).pobj
template mask*(x:Masked):expr = ((MaskedObj[x.T])(x)).mask
template `pobj=`*(x:Masked;y:untyped):untyped = ((MaskedObj[x.T])(x)).pobj = y
template `mask=`*(x:Masked;y:untyped):untyped = ((MaskedObj[x.T])(x)).mask = y
template maskedObj*(x:typed; msk:int):expr =
  (MaskedObj[type(x)])(pobj:x.addr,mask:msk)
template masked*(x:typed; msk:int):expr =
  mixin isVector,isMatrix
  when isComplex(x):
    #ctrace()
    #asComplex((Masked[type(x)])(maskedObj(x,msk)))
    asVarComplex(masked(x[],msk))
  elif isVector(x):
    #ctrace()
    #asVarMatrix((Masked[type(x)])(maskedObj(x,msk)))
    asVarVector(masked(x[],msk))
  elif isMatrix(x):
    #ctrace()
    #asVarMatrix((Masked[type(x)])(maskedObj(x,msk)))
    asVarMatrix(masked(x[],msk))
  elif x is SomeNumber:
    #ctrace()
    #asVarMatrix((Masked[type(x)])(maskedObj(x,msk)))
    x
  else:
    #ctrace()
    (Masked[type(x)])(maskedObj(x,msk))
#template isRIC*(x:int):expr = true
#template isRIC*(m:Masked):expr = isRIC(m.pobj[])
#template isComplex*(m:Masked):expr = isComplex(m.pobj[])
template declaredComplex*(m:Masked):expr =
  mixin declaredComplex
  declaredComplex(m.pobj[])
template isVector*(m:Masked):expr =
  mixin isVector
  isVector(m.pobj[])
#template isMatrix*(m:Masked):expr =
#  mixin isMatrix
  #echo "isMatrix"
  #echo isMatrix(m.pobj[])
#  isMatrix(m.pobj[])
template mvLevel*(m:Masked):expr =
  mixin mvLevel
  mvLevel(m.pobj[])
template numNumbers*(m:Masked):expr = numNumbers(m[])
template len*(m:Masked):expr =
  mixin len
  len(m.pobj[])
template nrows*(m:Masked):expr =
  mixin nrows
  nrows(m.pobj[])
template ncols*(m:Masked):expr =
  mixin ncols
  ncols(m.pobj[])
template re*(m:Masked):expr =
  mixin re
  masked(m.pobj[].re, m.mask)
template im*(m:Masked):expr =
  mixin im
  masked(m.pobj[].im, m.mask)
#template `[]`*(x:Masked; i:int):expr = Masked(x:x.pobj[i],mask:x.mask)
#template `[]`*(m:Masked; i,j:int):untyped =
#  Masked(x:unsafeAddr(m.pobj[][i,j]), mask:m.mask)
proc `[]`*[T](m:Masked[T]):var T {.inline.} = m.pobj[]
template `[]`*(m:Masked; i:int):expr = masked(m[][i],m.mask)
#proc `[]`*[T](m:Masked[T]; i:int):auto =
#  var r:Masked[type(m.pobj[][i])]
#  #Masked(pobj:unsafeAddr(m.pobj[][i,j]), mask:m.mask)
#  r.pobj = addr(m.pobj[][i])
#  r.mask = m.mask
#  r
#proc `[]`*(m:Masked; i,j:int):var Masked[type(m.pobj[][i,j])] =
template `[]`*(m:Masked; i,j:int):expr = masked(m[][i,j],m.mask)
  #var t = m[].addr
  #var tij = t[][i,j].addr
  #let tm = masked(tij[],m.mask)
  #ctrace()
  #tm
  #var r:Masked[type(tij)]
  #var r:Masked[type(m[][i,j])]
#  #Masked(pobj:unsafeAddr(m.pobj[][i,j]), mask:m.mask)
  #r.pobj = m.pobj[][i,j].addr
  #r.mask = m.mask
  #0
#template `[]`*(m:Masked; i,j:int):untyped =
#  Masked[type((pobj:m[][i,j].addr,mask:m.mask)
#template `[]=`*(m:Masked; i,j:int; y:untyped):untyped =
#  set(Masked(pobj:m.pobj[i,j].addr,mask:x.mask), y)
#proc `:=`*(x:Masked; y:int) =
#  mixin assign
#  var t = x
#  assign(t, y)
#proc `*=`*(x:Masked; y:int) =
  #echo "*="
#  mixin mul
  #echoAll isMatrix(t)
  #echoAll isMatrix(x.pobj[])
  #echoAll isScalar(y)
#  mul(x, x[], y)
  #echo "*="
#proc `$`*(x:Masked):string =
#  result = $(x[])
