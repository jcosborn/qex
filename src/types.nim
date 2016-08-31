import stdUtils
import comms
#import complexConcept
#import matrixConcept
import metaUtils
import wrapperTypes
import macros

# converter adjMat*[T](x:Adjointed[AsMatrix[T]]):AsMatrix[Adjointed[T]]

type
  #Adjointed*{.borrow: `.`.}[T] = distinct T
  #Adjointed*[T] = distinct T
  Adjointed*[T] = object
    v*:T
template adj*(xx:typed):untyped =
  lets(x,xx):
    when isComplex(x):
      asComplex(adj(x[]))
    elif isVector(x):
      asVector(adj(x[]))
    elif isMatrix(x):
      asMatrix(adj(x[]))
    elif x is SomeNumber:
      x
    else:
      when compiles(addr(x)):
      #when compiles(unsafeAddr(x)):
        cast[ptr Adjointed[type(x)]](addr(x))[]
        #cast[Adjointed[type(x)]](x)
      else:
        #(Adjointed[type(x)])(x)
        cast[Adjointed[type(x)]](x)
        #cast[Adjointed[type(x)]]((var t=x; t.addr))
#template `[]`*[T](x:Adjointed[T]):untyped = cast[T](x)
makeDeref(Adjointed, x.T)
template `[]`*(x:Adjointed; i:SomeInteger):untyped = x[][i].adj
template `[]`*(x:Adjointed; i,j:SomeInteger):untyped = x[][j,i].adj
template len*(x:Adjointed):untyped = x[].len
template nrows*(x:Adjointed):untyped = x[].ncols
template ncols*(x:Adjointed):untyped = x[].nrows
template isVector*(x:Adjointed):untyped = isVector(x[])
template isMatrix*(x:Adjointed):untyped = isMatrix(x[])
template re*(x:Adjointed):untyped = x[].re
template im*(x:Adjointed):untyped = -(x[].im)
#template mvLevel*(x:Adjointed):untyped =
#  mixin mvLevel
#  mvLevel(x[])


type
  MaskedObj*[T] = object
    pobj*:ptr T
    mask*:int
  Masked*[T] = distinct MaskedObj[T]
template pobj*(x:Masked):untyped = ((MaskedObj[x.T])(x)).pobj
template mask*(x:Masked):untyped = ((MaskedObj[x.T])(x)).mask
template `pobj=`*(x:Masked;y:untyped):untyped = ((MaskedObj[x.T])(x)).pobj = y
template `mask=`*(x:Masked;y:untyped):untyped = ((MaskedObj[x.T])(x)).mask = y
template maskedObj*(x:typed; msk:int):untyped =
  (MaskedObj[type(x)])(pobj:x.addr,mask:msk)
template masked*(x:typed; msk:int):untyped =
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
#template isRIC*(x:int):untyped = true
#template isRIC*(m:Masked):untyped = isRIC(m.pobj[])
#template isComplex*(m:Masked):untyped = isComplex(m.pobj[])
template declaredComplex*(m:Masked):untyped =
  mixin declaredComplex
  declaredComplex(m.pobj[])
template isVector*(m:Masked):untyped =
  mixin isVector
  isVector(m.pobj[])
#template isMatrix*(m:Masked):untyped =
#  mixin isMatrix
  #echo "isMatrix"
  #echo isMatrix(m.pobj[])
#  isMatrix(m.pobj[])
template mvLevel*(m:Masked):untyped =
  mixin mvLevel
  mvLevel(m.pobj[])
template numNumbers*(m:Masked):untyped = numNumbers(m[])
template len*(m:Masked):untyped =
  mixin len
  len(m.pobj[])
template nrows*(m:Masked):untyped =
  mixin nrows
  nrows(m.pobj[])
template ncols*(m:Masked):untyped =
  mixin ncols
  ncols(m.pobj[])
template re*(m:Masked):untyped =
  mixin re
  masked(m.pobj[].re, m.mask)
template im*(m:Masked):untyped =
  mixin im
  masked(m.pobj[].im, m.mask)
template `re=`*(m:Masked; x:any):untyped =
  mixin re
  assign(masked(m.pobj[].re, m.mask), x)
template `im=`*(m:Masked; x:any):untyped =
  mixin im
  assign(masked(m.pobj[].im, m.mask), x)
#template `[]`*(x:Masked; i:int):untyped = Masked(x:x.pobj[i],mask:x.mask)
#template `[]`*(m:Masked; i,j:int):untyped =
#  Masked(x:unsafeAddr(m.pobj[][i,j]), mask:m.mask)
proc `[]`*[T](m:Masked[T]):var T {.inline.} = m.pobj[]
template `[]`*(m:Masked; i:int):untyped = masked(m[][i],m.mask)
#proc `[]`*[T](m:Masked[T]; i:int):auto =
#  var r:Masked[type(m.pobj[][i])]
#  #Masked(pobj:unsafeAddr(m.pobj[][i,j]), mask:m.mask)
#  r.pobj = addr(m.pobj[][i])
#  r.mask = m.mask
#  r
#proc `[]`*(m:Masked; i,j:int):var Masked[type(m.pobj[][i,j])] =
template `[]`*(m:Masked; i,j:int):untyped = masked(m[][i,j],m.mask)
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
