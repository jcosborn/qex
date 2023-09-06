import macros
import base/basicOps
import base/metaUtils
#import stdUtils
#import metaUtils
#import basicOps
import maths/types
import base/wrapperTypes
import complexNumbers
import simd/simdWrap
#import complexType
# opinc, opdec
# minc, mdec, redotinc, redotdec
# unary ops: assign(=),neg(-),iadd(+=),isub(-=),imul(*=),idiv(/=)
# binary ops: add(+),sub(-),mul(*),divd(/),imadd(+=*),imsub(-=*)
# iaddmul, isubmul
# ternary ops: madd(*+),msub(*-),nmadd(-*+),nmsub(-*-)
# assign: trace,dot,outer
# wrap: conj,adj,transpose
# norm2,det,lu,qr,svd,eig
# sqrt,rsqrt,exp,log,groupProject,groupCheck
# mapX(f,x,r),mapXY(f,x,y,r)
getOptimPragmas()

template createAsType2(t,c:untyped) =
  mixin `[]`, index
  #makeWrapper(t, c)
  makeWrapperType(t)
  template `[]`*[T](x: t; i: T): untyped =
    when T is t:
      x[][i[]]
    elif T.isWrapper:
      indexed(x, i)
    else:
      x[][i]
  template index*[X:t,T](x: typedesc[X]; i: typedesc[T]): typedesc =
    mixin index
    when T is t:
      index(X[], T[])
    elif T.isWrapper:
      c(index(X.type[], type T))
    else:
      index(X[], T)
  template `[]`*(x: t; i: Scalar): untyped = c(x[][i])
  template `[]`*(x:t; i,j:SomeInteger):untyped =
    #echoType: x
    #ctrace()
    x[][i,j]
    #let t1 = x[]
    #ctrace()
    #echoType: x[]
    #mixin `[]`
    #let t2 = t1[i,j]
    #let t2 = t1[i][j]
    #ctrace()
    #t2
  #template `[]`*(x: t; i,j: SomeInteger): untyped =
  #  x[][i,j]
  template `[]=`*[T](x: t; i: T; y: typed) =
    when T is t:
      x[][i[]] = y
    elif T.isWrapper:
      indexed(x, i) := y
    else:
      x[][i] = y
  template `[]=`*(x: t; i,j: SomeInteger; y: typed) =
    x[][i,j] = y
  template len*(x:t):untyped = getConst(x[].len)
  template nrows*(x:t):untyped = getConst(x[].nrows)
  template ncols*(x:t):untyped = getConst(x[].ncols)
  template len*(x:typedesc[t]):auto = getConst(x[].len)
  #template mvLevel*(x:t):untyped =
  #  mixin mvLevel
  #  mvLevel(x[])
#template createAsType(t:untyped):untyped = createAsType2(`As t`, `as t`)
macro createAsType(t: untyped): untyped =
  newCall(bindsym("createAsType2"), ident("As" & t.repr), ident("as" & t.repr))

#createAsType(Scalar)
#createAsType(VarScalar)
createAsType(Vector)
#createAsType(VarVector)
createAsType(Matrix)
#createAsType(VarMatrix)

template eval*[T](x: typedesc[AsVector[T]]): typedesc =
  mixin eval
  asVector(eval(typeof T))
template eval*[T](x: typedesc[AsMatrix[T]]): typedesc =
  mixin eval
  asMatrix(eval(typeof T))

template has*[T:AsVector](x: typedesc[T], y: typedesc): bool =
  mixin has
  when y is AsVector: true
  else: has(T.type[], y)
template has*[T:AsMatrix](x: typedesc[T], y: typedesc): bool =
  mixin has
  when y is AsMatrix: true
  else: has(T.type[], y)

#declareScalar(AsScalar)
#declareScalar(AsVarScalar)
#declareVector(AsVector)
#declareVector(AsVarVector)
#declareMatrix(AsMatrix)
#declareMatrix(AsVarMatrix)
#template deref(x:typed):untyped =
#  #when type(x) is AsScalar|AsVarScalar:
#  #when type(x) is AsScalar:
#  #  x[]
#  #else:
#    x

type
  #VectorArrayObj*[I:static[int],T] = array[I,T]
  VectorArrayObj*[I:static[int],T] = object
    vec*: array[I,T]
  VectorArray*[I:static[int],T] = AsVector[VectorArrayObj[I,T]]
  #VectorArray*[I:static[int],T] = AsVector[array[I,T]]
  MatrixArrayObj*[I,J:static[int],T] = object
    mat*: array[I,array[J,T]]
  MatrixArray*[I,J:static[int],T] = AsMatrix[MatrixArrayObj[I,J,T]]
  MatrixRowObj*[T] = object
    rw*: int
    mat*: ptr T
  MatrixRow*[T] = AsVector[MatrixRowObj[T]]
  #MatrixCol*[T] = tuple[col:int,mat:ptr T]
  #MatrixDiag*[T] = tuple[diag:int,mat:ptr T]

  #Vec1* = concept x
  #  #mixin isVector
  #  x.isVector
  #Vec2* = concept x
  #  mixin isVector
  #  x.isVector
  #Vec3* = concept x
  #  #mixin isVector
  #  x.isVector
  Vec1*[T] = AsVector[T]
  Vec2*[T] = AsVector[T]
  Vec3*[T] = AsVector[T]
  #Mat1* = concept x
  #  #mixin isMatrix
  #  x.isMatrix
  #Mat2* = concept x
  #  #mixin isMatrix
  #  x.isMatrix
  #Mat3* = concept x
  #  #mixin isMatrix
  #  x.isMatrix
  Mat1*[T] = AsMatrix[T]
  Mat2*[T] = AsMatrix[T]
  Mat3*[T] = AsMatrix[T]
  Mat4*[T] = AsMatrix[T]
  MV1* = Mat1 | Vec1
  MV2* = Mat2 | Vec2
  MV3* = Mat3 | Vec3
  MVconcept1* = Vec1 | Mat1
  MVconcept2* = Vec2 | Mat2
  MVconcept3* = Vec3 | Mat3
  Sca1* = concept x
    x isnot MVconcept1
  Sca2* = concept x
    x isnot MVconcept2
  Sca3* = concept x
    x isnot MVconcept3
  #VarSca1* = var Sca1
  #VarVec1* = var Vec1
  #VarMat1* = var Mat1
  #VarMV1* = AsVarMatrix | AsVarVector
  #VarAuto* = var auto #| AsVarMatrix
  #AsVarVector*[T] = AsVar[AsVector[T]]
  #AsVarMatrix*[T] = AsVar[AsMatrix[T]]

template isWrapper*(x: array): untyped = false

#template isWrapper*(x: AsVector): untyped = true
#template asWrapper*(x: AsVector, y: typed): untyped =
#  #static: echo "asWrapper AsVector"
#  asVector(y)
#template asVarWrapper*(x: AsVector, y: typed): untyped =
#  #static: echo "asVarWrapper AsVector"
#  asVar(asVector(y))

template `len`*(x:VectorArrayObj):untyped = x.I
template `len`*[I:static[int],T](x:typedesc[VectorArrayObj[I,T]]):auto = I
template `[]`*(x:VectorArrayObj):untyped = x.vec
template `[]`*(x:VectorArrayObj; i:int):untyped =
  #var vecIdx = toRef x.vec[i]
  #vecIdx
  x.vec[i]
#template `[]`*(x:var VectorArrayObj; i:int):untyped =
#  toRef x.vec[i]
template `[]=`*(x:VectorArrayObj; i:int, y:typed):untyped =
  x.vec[i] := y
template asVectorArray*[N:static[int],T](x: array[N,T]): untyped =
  #static: echo "asVectorArray"
  #let x_asVectorArray = xx
  #const n_asVectorArray = x_asVectorArray.len
  #static: echo "asVectorArray2"
  #asVector( VectorArrayObj[n_asVectorArray,
  #                         type(x_asVectorArray[0])](x_asVectorArray) )
  #asVector( VectorArrayObj[N,type(T)](x) )
  asVector( VectorArrayObj[N,type(T)](vec: x) )
  #let t1 = VectorArrayObj[n_asVectorArray,
  #                        type(x_asVectorArray[0])](vec: x_asVectorArray)
  #static: echo "asVectorArray1"
  #let t = asVector(t1)
  #static: echo "asVectorArray2"
  #t
template asVectorArray*[T](N:static[int], x: typedesc[T]): untyped =
  asVector( VectorArrayObj[N,type(T)] )

template index*[I:static[int],T,K](x: typedesc[VectorArrayObj[I,T]];
                                   k: typedesc[K]): typedesc =
  when K.isWrapper:
    VectorArrayObj[I,index(type T, type K)]
  else:
    false # error

template index*[I,J:static[int],T,K](x: typedesc[MatrixArrayObj[I,J,T]];
                                     k: typedesc[K]): typedesc =
  when K.isWrapper:
    MatrixArrayObj[I,J,index(type T, type K)]
  else:
    false # error

template `len`*(x:MatrixArrayObj):untyped = x.I
template nrows*(x:MatrixArrayObj):untyped = x.I
template ncols*(x:MatrixArrayObj):untyped = x.J
template `[]`*(x:MatrixArrayObj):untyped = x.mat
template `[]`*(x:MatrixArrayObj; i:Scalar):untyped = indexed(x, i[])
template `[]`*(x:MatrixArrayObj; i,j:int):untyped = x.mat[i][j]
template `[]`*(x:var MatrixArrayObj; i,j:int):untyped = x.mat[i][j]
template `[]=`*(x:MatrixArrayObj; i,j:int, y:typed):untyped = x.mat[i][j] = y
template numberType*[T](x:AsVector[T]):typedesc = numberType(type(T))
template numberType*[T](x:AsMatrix[T]):typedesc = numberType(type(T))
template numberType*[T](x:typedesc[AsVector[T]]):typedesc = numberType(type(T))
template numberType*[T](x:typedesc[AsMatrix[T]]):typedesc = numberType(type(T))
template numberType*[I,T](x:VectorArrayObj[I,T]):typedesc =
  numberType(type(T))
template numberType*[I,T](x:typedesc[VectorArrayObj[I,T]]):typedesc =
  numberType(type(T))
template numberType*[I,J,T](x:MatrixArrayObj[I,J,T]):typedesc =
  numberType(type(T))
template numberType*[I,J,T](x:typedesc[MatrixArrayObj[I,J,T]]):typedesc =
  numberType(type(T))
template numNumbers*(x:AsVector):untyped =
  mixin numNumbers
  x.len*numNumbers(x[0])
template numNumbers*(x:AsMatrix):untyped =
  x.nrows*x.ncols*numNumbers(x[0,0])
#template `[]`*(x:array; i,j:int):untyped = x[i][j]
#template `[]=`*(x:array; i,j:int, y:untyped):untyped = x[i][j] = y
template getNc*(x: AsVector): untyped = getNc(x[0])
template getNc*(x: AsMatrix): untyped = getNc(x[0,0])
template getNs*(x: AsVector): untyped = getNs(x[0])
template getNs*(x: AsMatrix): untyped = getNs(x[0,0])

template simdLength*[T](x:typedesc[AsVector[T]]):untyped = simdLength(type(T))
template simdLength*[I,T](x:typedesc[VectorArrayObj[I,T]]):untyped =
  simdLength(type(T))
template simdLength*[T](x:typedesc[AsMatrix[T]]):untyped = simdLength(type(T))
template simdLength*[I,J,T](x:typedesc[MatrixArrayObj[I,J,T]]):untyped =
  simdLength(type(T))

template len*(x:MatrixRowObj):untyped = x.mat[].ncols
template `[]`*(x:MatrixRowObj; i:int):untyped = x.mat[][x.rw,i]
template `[]=`*(x:MatrixRowObj; i:int; y:typed):untyped = x.mat[][x.rw,i] = y

template isWrapper*(x: VectorArrayObj): bool = false
template isWrapper*(x: type VectorArrayObj): bool = false
template isWrapper*(x: MatrixArrayObj): bool = false
template isWrapper*(x: type MatrixArrayObj): bool = false

template eval*[I,T](x: typedesc[VectorArrayObj[I,T]]): typedesc =
  mixin eval
  VectorArrayObj[I,eval(typeof T)]
template eval*[I,J,T](x: typedesc[MatrixArrayObj[I,J,T]]): typedesc =
  mixin eval
  MatrixArrayObj[I,J,eval(typeof T)]

template has*[I,T](x: typedesc[VectorArrayObj[I,T]], y: typedesc): bool =
  mixin has
  has(T.type, y)
template has*[I,J,T](x: typedesc[MatrixArrayObj[I,J,T]], y: typedesc): bool =
  mixin has
  has(T.type, y)

template isWrapper*(x: AsVar[AsMatrix]): untyped = true
template asWrapper*(x: AsVar[AsMatrix], y: typed): untyped =
  #static: echo "asWrapper AsVarMatrix"
  asVar(asMatrix(v: y))
template asVarWrapper*(x: AsVar[AsMatrix], y: typed): untyped =
  #static: echo "asVarWrapper AsVarMatrix"
  asVar(asMatrix(v: y))

template toSingle*[I,T](x: typedesc[VectorArrayObj[I,T]]): typedesc =
  mixin toSingle
  VectorArrayObj[I,toSingle(type(T))]
template toSingle*[T](x: typedesc[AsVector[T]]): typedesc =
  AsVector[toSingle(type(T))]
template toSingleImpl*(x: VectorArrayObj): untyped =
  mixin toSingleX
  toSingleX(toRef x)
template toSingle*[I,J,T](x: typedesc[MatrixArrayObj[I,J,T]]): typedesc =
  mixin toSingle
  MatrixArrayObj[I,J,toSingle(type(T))]
template toSingle*[T](x: typedesc[AsMatrix[T]]): typedesc =
  AsMatrix[toSingle(type(T))]
template toSingleImpl*(x: MatrixArrayObj): untyped =
  mixin toSingleX
  toSingleX(toRef x)

template toDouble*[I,T](x: typedesc[VectorArrayObj[I,T]]): untyped =
  mixin toDouble
  VectorArrayObj[I,toDouble(type(T))]
template toDouble*[T](x: typedesc[AsVector[T]]): untyped =
  AsVector[toDouble(type(T))]
template toDoubleImpl*(x: VectorArrayObj): untyped =
  mixin toDoubleX
  toDoubleX(toRef x)
template toDouble*[I,J,T](x: typedesc[MatrixArrayObj[I,J,T]]): typedesc =
  mixin toDouble
  MatrixArrayObj[I,J,toDouble(type(T))]
template toDouble*[T](x: typedesc[AsMatrix[T]]): typedesc =
  AsMatrix[toDouble(type(T))]
template toDoubleImpl*(x: MatrixArrayObj): untyped =
  mixin toDoubleX
  toDoubleX(toRef x)

#template masked*(x: AsMatrix, msk: typed): untyped =
#  static: echo "masked AsMatrix"
#  asVarMatrix(masked(x[],msk))
#template masked*(x: AsVarMatrix, msk: typed): untyped =
#  asVarMatrix(masked(x[],msk))

#template isVector(x:Row):untyped = true
#template isVector(x:Col):untyped = true
#template mvLevel(x:Sca1):untyped = -1

#template simpleAssign(
#template tmpluntyped*(x:typed):untyped = x
#template tmptype*(x:Vec1):untyped = VectorArray[x.len,type(load1(x[0]))]

template load1*(x: VectorArrayObj): untyped = x
  #mixin load1
  #let x = xx
  #var r_load1V{.noInit.}: VectorArray[x.len,type(load1(x[0]))]
  #assign(r_load1V, x)
  #r_load1V
  #asVector(load1(x[]))
template load1*(xx: Vec1): untyped =
  mixin load1
  let x = xx
  var r_load1V{.noInit.}: VectorArray[x.len,evalType(load1(x[0]))]
  assign(r_load1V, x)
  r_load1V
  #asVector(load1(x[]))
template load1*(xx: AsVar[AsVector]):untyped =
  mixin load1
  lets(x,xx):
    var r{.noInit.}:VectorArray[x.len,evalType(load1(x[0]))]
    assign(r, x)
    r
template load1*(xx: Mat1): untyped =
  lets(x,xx):
    var r_load1M{.noInit.}: MatrixArray[x.nrows,x.ncols,evalType(load1(x[0,0]))]
    assign(r_load1M, x)
    r_load1M
#template tmpvar1*(x:Vec1):untyped =
#  lets(x,xx):
#    var r{.noInit.}:VectorArray[x.len,type(load1(x[0]))]
#    #assign(r, x)
#    r

template row*(x:AsVector; i:int):untyped = x
#proc row*(x:AsMatrix; i:int):auto {.inline,noInit.} =
template row*(x:AsMatrix; i:int):untyped =
  #const nc = x.ncols
  #var r{.noInit.}:VectorArray[nc,type(load1(x[0,0]))]
  #for j in 0..<nc:
  #  assign(r[j], x[i,j])
  #r
  asVector(MatrixRowObj[type(x)](rw:i,mat:unsafeAddr(x)))
template setRow*(r:AsVector; x:AsVector; i:int):untyped =
  assign(r, x)
#proc setRow*(r:var AsMatrix; x:AsVector; i:int) {.inline.} =
template setRow*(r:AsMatrix; xx:AsVector; ii:int) =
  block:
    let xp = getPtr xx; template x:untyped = xp[]
    let ip = getPtr ii; template i:untyped = ip[]
    const nc = getConst(r.ncols)
    for j in 0..<nc:
      assign(r[i,j], x[j])
template column*(x:AsVector; i:int):untyped = x
proc column*(x:AsMatrix; i:int):auto {.inline,noInit.} =
  const nr = x.nrows
  var r{.noInit.}:VectorArray[nr,evalType(x[0,0])]
  for j in 0..<nr:
    assign(r[j], x[j,i])
  r
template setColumn*(r:AsVector; x:AsVector; i:int):untyped =
  assign(r, x)
proc setColumn*(r:var AsMatrix; x:AsVector; i:int) {.inline.} =
  const nr = r.nrows
  for j in 0..<nr:
    assign(r[j,i], x[j])

import matrixOps
export matrixOps

proc toString*(x:Vec1):string =
  mixin `$`
  result = $(x[0])
  forO i, 1, x.len-1:
    result.add "," & $(x[i])
proc `toString`*(x:Mat1):string =
  mixin `$`
  result = ""
  forO i, 0, x.nrows-1:
    result.add $(x[i,0])
    forO j, 1, x.ncols-1:
      result.add "," & $(x[i,j])
    if i<x.nrows-1: result.add "\n"

#template makeLevel1P(f,s1,t1,s2,t2:untyped):untyped {.dirty.} =
#  proc f*(r:t1, x:t2) {.inline.} =
#    `f s1 s2`(r, deref(x))
template makeLevel1T(f,s1,t1,s2,t2:untyped):untyped {.dirty.} =
  #template `f U`*(r: t1, x: t2): untyped =
  #  `f s1 s2`(r, x)
  template f*(r: t1, x: t2): untyped =
    #flattenCallArgs(`f U`, r, x)
    `f s1 s2`(r, x)
template makeLevel1(f,s1,t1,s2,t2:untyped):untyped =
  makeLevel1T(f,s1,t1,s2,t2)

#macro makeLevel2(f,s1,t1,s2,t2,s3,t3:untyped):auto =
#  let f123 = ident($f & $s1 & $s2 & $s3)
#  result = quote do:
#    proc `f`*(r:`t1`; x:`t2`; y:`t3`) {.inline.} =
#      `f123`(r, x, y)
#template func3(f,a,b,c: untyped): untyped = f(a,b,c)
#template makeLevel2P(f,s1,t1,s2,t2,s3,t3:untyped):untyped {.dirty.} =
#  proc f*(r:t1, x:t2, y:t3) {.inline.} =
#    #`f s1 s2 s3`(r, x, y)
#    func3(`f s1 s2 s3`, r, x, y)
template makeLevel2T(f,s1,t1,s2,t2,s3,t3: untyped): untyped {.dirty.} =
  #template `f U`*(r: t1, x: t2, y: t3): untyped =
  #  `f s1 s2 s3`(r, x, y)
  template f*(r: t1, x: t2, y: t3): untyped =
    #echoType: r
    #flattenCallArgs(`f U`, r, x, y)
    #`f U`(r, x, y)
    `f s1 s2 s3`(r, x, y)
template makeLevel2(f,s1,t1,s2,t2,s3,t3:untyped):untyped {.dirty.} =
  makeLevel2T(f,s1,t1,s2,t2,s3,t3)

template makeMap1(op:untyped):untyped =
  makeLevel1(op, S, var Sca1, V, Vec2)
  makeLevel1(op, S, var Sca1, M, Mat2)
  makeLevel1(op, V, var Vec1, S, Sca2)
  makeLevel1(op, V, var Vec1, V, Vec2)
  #makeLevel1(op, V, var Vec1, V, AsVarVector)
  #makeLevel1(op, V, AsVarVector, S, Sca2)
  #makeLevel1(op, V, AsVarVector, V, Vec2)
  makeLevel1(op, M, var Mat1, S, Sca2)
  makeLevel1(op, M, var Mat1, V, Vec2)
  makeLevel1(op, M, var Mat1, M, Mat2)
  #makeLevel1(op, M, AsVarMatrix, S, Sca2)
  #makeLevel1(op, M, AsVarMatrix, V, Vec2)
  #makeLevel1(op, M, AsVarMatrix, M, Mat2)

makeMap1(assign)
makeMap1(neg)
makeMap1(iadd)
makeMap1(isub)

setUnop(`-`,neg,Vec1,VectorArray[x.len,evalType(x[0])])
setUnop(`-`,neg,Mat1,MatrixArray[x.nrows,x.ncols,evalType(x[0,0])])

#template assign*(x:Mat1; y:SomeNumber) =
#  echo "test"

template `:=`*(x:var Vec1; y:SomeNumber) = assign(x, y)
#template `:=`*(x:AsVarVector; y:SomeNumber) = assign(x, y)
template `:=`*(x:var Vec1; y:Vec2): untyped = assign(x, y)
#template `:=`*(x:AsVarVector; y:Vec2): untyped = assign(x, y)
template `:=`*(x:var Mat1; y:SomeNumber) = assign(x, y)
template `:=`*(x:var Mat1; y:Simd) = assign(x, y)
template `:=`*(x:var Mat1; y:AsComplex) = assign(x, y)
#template `:=`*(x:AsVarMatrix; y:SomeNumber) = assign(x, y)
template `:=`*(x:var Mat1; y:Vec2) = assign(x, y)
#template `:=`*(x:AsVarMatrix; y:Vec2) = assign(x, y)
template `:=`*(x:var Mat1; y:Mat2) = assign(x, y)
#template `:=`*(x:AsVarMatrix; y:Mat2) = assign(x, y)

template `+=`*(x: var Vec1; y: Vec2) =
  staticTraceBegin: peqVV
  iadd(x, y)
  staticTraceEnd: peqVV
template `+=`*(x:var Mat1; y:Mat2) = iadd(x, y)
template `+=`*(x:var Mat1; y:Sca2) = iadd(x, y)

template `-=`*(x:var Vec1; y:Vec2) = isub(x, y)
template `-=`*(x:var Mat1; y:Sca2) = isub(x, y)
template `-=`*(x:var Mat1; y:Mat2) = isub(x, y)

makeLevel1(imul, V, var Vec1, S, Sca2)
makeLevel1(imul, M, var Mat1, S, Sca2)
#makeLevel1(imul, M, AsVarMatrix, S, Sca2)

template `*=`*(r:var MV1; x:Sca2) = imul(r, x)
template `/=`*(r:var MV1; x:Sca2) = imul(r, x)
#template `*=`*(r:AsVarMatrix; x:Sca2) = imul(r, x)

template makeMap2(op:untyped):untyped =
  makeLevel2(op, V, var Vec1, V, Vec2, S, Sca3)
  makeLevel2(op, V, var Vec1, S, Sca2, V, Vec3)
  makeLevel2(op, V, var Vec1, V, Vec2, V, Vec3)
  makeLevel2(op, M, var Mat1, S, Sca2, S, Sca3)
  makeLevel2(op, M, var Mat1, V, Vec2, S, Sca3)
  makeLevel2(op, M, var Mat1, S, Sca2, V, Vec3)
  makeLevel2(op, M, var Mat1, V, Vec2, V, Vec3)
  makeLevel2(op, M, var Mat1, M, Mat2, S, Sca3)
  makeLevel2(op, M, var Mat1, S, Sca2, M, Mat3)
  makeLevel2(op, M, var Mat1, M, Mat2, V, Vec3)
  makeLevel2(op, M, var Mat1, V, Vec2, M, Mat3)
  makeLevel2(op, M, var Mat1, M, Mat2, M, Mat3)

makeMap2(add)
makeMap2(sub)

setBinop(`+`,add,Vec1,Sca2,VectorArray[x.len,evalType(x[0]+y)])
setBinop(`-`,sub,Vec1,Sca2,VectorArray[x.len,evalType(x[0]-y)])

setBinop(`+`,add,Vec1,Vec2,VectorArray[x.len,evalType(x[0]+y[0])])
setBinop(`-`,sub,Vec1,Vec2,VectorArray[x.len,evalType(x[0]-y[0])])

setBinop(`+`,add,Sca1,Mat2,MatrixArray[y.nrows,y.ncols,evalType(x+y[0,0])])
setBinop(`-`,sub,Sca1,Mat2,MatrixArray[y.nrows,y.ncols,evalType(x-y[0,0])])

setBinop(`+`,add,Mat1,Sca2,MatrixArray[x.nrows,x.ncols,evalType(x[0,0]+y)])
setBinop(`-`,sub,Mat1,Sca2,MatrixArray[x.nrows,x.ncols,evalType(x[0,0]-y)])

setBinop(`+`,add,Mat1,Mat2,MatrixArray[x.nrows,x.ncols,evalType(x[0,0]+y[0,0])])
setBinop(`-`,sub,Mat1,Mat2,MatrixArray[x.nrows,x.ncols,evalType(x[0,0]-y[0,0])])

makeLevel2(mul, V, var Vec1, V, Vec2, S, Sca3)

makeLevel2(mul, V, var Vec1, S, Sca2, V, Vec3)
#makeLevel2(mul, V, var AsVector, S, Sca2, V, Vec3)

#makeLevel2(op, S, Sca1, V, Vec2, V, Vec3)
makeLevel2(mul, M, var Mat1, M, Mat2, S, Sca3)
makeLevel2(mul, M, var Mat1, S, Sca2, M, Mat3)
makeLevel2(mul, V, var Vec1, M, Mat2, V, Vec3)
#makeLevel2(op, V, Vec1, V, Vec2, M, Mat3)
makeLevel2(mul, M, var Mat1, M, Mat2, M, Mat3)
#makeLevel2(op, M, Mat1, S, Sca2, S, Sca3)
#makeLevel2(op, M, Mat1, V, Vec2, S, Sca3)
#makeLevel2(op, M, Mat1, S, Sca2, V, Vec3)
#makeLevel2(op, M, Mat1, V, Vec2, V, Vec3)
#makeLevel2(op, M, Mat1, V, Vec2, M, Mat3)
setBinop(mul,mul, Sca1,Vec2,VectorArray[getConst(y.len),evalType(x*y[0])])

#setBinop(`*`,mul, Sca1,AsVector,VectorArray[y.len,evalType(x*y[0])])
#setBinop(`*`,mul, float,Vec2,VectorArray[y.len,evalType(x*y[0])])
#setBinop(`*`,mul, AsScalar,Vec2,VectorArray[y.len,evalType(x*y[0])])
#template `*`*[X:Sca1,Y:Vec](x: typedesc[X], y: typedesc[Y]): typedesc =
#  VectorArray[len Y, X*indexed(Y,int)]
setBinop(`*`,mul, Sca1,Vec2,VectorArray[getConst(y.len),evalType(x*y[0])])
setBinop(`*`,mul, Vec1,Sca2,VectorArray[x.len,evalType(x[0]*y)])
setBinop(`*`,mul, Mat1,Vec2,VectorArray[x.nrows,evalType(x[0,0]*y[0])])

setBinop(`*`,mul, Sca1,Mat2,MatrixArray[y.nrows,y.ncols,evalType(x*y[0,0])])
setBinop(`*`,mul, Mat1,Sca2,MatrixArray[x.nrows,x.ncols,evalType(x[0,0]*y)])
setBinop(`*`,mul, Mat1,Mat2,MatrixArray[x.nrows,y.ncols,evalType(x[0,0]*y[0,0])])

makeLevel2(imadd, V, var Vec1, S, Sca2, V, Vec3)
#makeLevel2(imadd, V, AsVarVector, S, Sca2, V, Vec3)
makeLevel2(imadd, V, var Vec1, M, Mat2, V, Vec3)
#makeLevel2(imadd, V, AsVarVector, M, Mat2, V, Vec3)
makeLevel2(imadd, M, var Mat1, M, Mat2, M, Mat3)

makeLevel2(imsub, V, var Vec1, S, Sca2, V, Vec3)
makeLevel2(imsub, V, var Vec1, M, Mat2, V, Vec3)
#makeLevel2(imsub, V, AsVarVector, M, Mat2, V, Vec3)
makeLevel2(imsub, M, var Mat1, M, Mat2, M, Mat3)

#proc imadd*(r:VarVec1; x:Mat2; y:Vec3) {.inline.} = imaddVMV(r, x, y)
#proc imadd*(r:AsVarVector; x:Mat2; y:Vec3) {.inline.} = imaddVMV(r, x, y)
#proc imadd*(r:VarMat1; x:Mat2; y:Mat3) {.inline.} = imaddMMM(r, x, y)

#proc imsub*(r:var Vec1; x:Sca2; y:Vec3) {.inline.} = imsubVSV(r, x, y)
#proc imsub*(r:VarVec1; x:Mat2; y:Vec3) {.inline.} = imsubVMV(r, x, y)
#proc imsub*(r:AsVarVector; x:Mat2; y:Vec3) {.inline.} = imsubVMV(r, x, y)
#proc imsub*(r:VarMat1; x:Mat2; y:Mat3) {.inline.} = imsubMMM(r, x, y)

proc msub*(r:var Vec1; x:auto; y:Vec2; z:Vec3) {.inline.} = msubVSVV(r,x,y,z)

#proc trace*(r:var Sca1; x:Mat2) {.inline.} =
proc trace*(r: var auto; x: Mat2) {.alwaysInline.} =
  mixin nrows, ncols, trace, iadd
  let n = min(x.nrows, x.ncols)
  assign(r, 0)
  for i in 0..<n:
    let t = trace(x[i,i])
    iadd(r, t)
proc trace*(x: Mat1): auto {.alwaysInline,noInit.} =
  var t{.noInit.}: evalType(trace(x[0,0]))
  #static: echo "trace"
  trace(t, x)
  t

# FIXME: make generic reduce

#proc inorm2*(r:var auto; x:Vec2) {.inline.} =
#  mixin inorm2
#  for i in 0..<x.len:
#    #echo r
#    inorm2(r, x[i])

template inorm2*(r: var auto; x: Vec2) =
  mixin inorm2
  let xx = x
  for i in 0..<xx.len:
    #echo r
    inorm2(r, xx[i])


proc inorm2*(r:var auto; x:Mat2) {.inline.} =
  mixin nrows, ncols, inorm2
  for i in 0..<x.nrows:
    for j in 0..<x.ncols:
      inorm2(r, x[i,j])
proc norm2*(r:var auto; x:Vec2) {.inline.} =
  mixin norm2, iadd
  assign(r, 0)
  for i in 0..<x.len:
    var t{.noInit.}:evalType(r)
    norm2(t, x[i])
    iadd(r, t)
#proc norm2*(r:var auto; x: AsVarVector) {.inline.} =
#  mixin norm2, iadd
#  assign(r, 0)
#  for i in 0..<x.len:
#    var t{.noInit.}:evalType(r)
#    norm2(t, x[i])
#    iadd(r, t)
proc norm2*(r:var auto; x:Mat2) {.inline.} =
  mixin nrows, ncols, norm2, iadd
  assign(r, 0)
  for i in 0..<x.nrows:
    for j in 0..<x.ncols:
      var t{.noInit.}:evalType(r)
      norm2(t, x[i,j])
      iadd(r, t)
proc norm2*(x:Vec1):auto {.inline,noInit.} =
  var t{.noInit.}:evalType(norm2(x[0]))
  norm2(t, x)
  t
#proc norm2*(x: AsVarVector): auto {.inline,noInit.} =
#  var t{.noInit.}:evalType(norm2(x[0]))
#  norm2(t, x)
#  t
proc norm2*(x:Mat1):auto {.inline,noInit.} =
  var t{.noInit.}:evalType(norm2(x[0,0]))
  norm2(t, x)
  t
proc norm2X*(x:Vec1):auto {.inline,noInit.} =
  var t{.noInit.}:evalType(norm2X(x[0]))
  norm2(t, x)
  t
proc norm2X*(x:Mat1):auto {.inline,noInit.} =
  var t{.noInit.}:evalType(norm2X(x[0,0]))
  norm2(t, x)
  t

#proc idot*(r:var Sca1; x:Vec2; y:Vec3) {.inline.} =
template idot*(r: var Sca1; xx: Vec2; yy: Vec3) =
  let x = xx
  let y = yy
  for i in 0..<len(x):
    r += dot(x[i],y[i])
    #idot(r,x[i],y[i])
#proc dot*(r:var Sca1; x:Vec2; y:Vec3) {.inline.} =
template dot*(r: var Sca1; x: Vec2; y: Vec3) =
  r := 0
  idot(r, x, y)
setBinop(dot, dot, Vec1, Vec2, evalType(dot(x[0],y[0])))

proc dot*(x: Mat2; y: Mat3): auto {.inline,noInit.} =
  mixin dot
  result = dot(x[0,0],y[0,0])
  forO j, 1, x.len.pred:
    result += dot(x[0,j],y[0,j])
  forO i, 1, x.len.pred:
    forO j, 0, x.len.pred:
      result += dot(x[i,j],y[i,j])

proc redot*(x:Vec2; y:Vec3): auto {.inline,noInit.} =
  result = redot(x[0],y[0])
  forO i, 1, x.len.pred:
    result += redot(x[i],y[i])

proc redot*(x: Mat2; y: Mat3): auto {.inline,noInit.} =
  result = redot(x[0,0],y[0,0])
  forO j, 1, x.len.pred:
    result += redot(x[0,j],y[0,j])
  forO i, 1, x.len.pred:
    forO j, 0, x.len.pred:
      result += redot(x[i,j],y[i,j])

proc peqOuter*(r: var Mat1; x: Vec2; y: Vec3) {.alwaysInline.} =
  forO i, 0, r.nrows-1:
    forO j, 0, r.ncols-1:
      r[i,j] += x[i] * y[j].adj
proc meqOuter*(r: var Mat1; x: Vec2; y: Vec3) {.alwaysInline.} =
  forO i, 0, r.nrows-1:
    forO j, 0, r.ncols-1:
      r[i,j] -= x[i] * y[j].adj

proc simdSum*(x: Vec1): auto {.noInit.} =
  var r{.noInit.}: VectorArray[x.len,evalType(simdSum(x[0]))]
  forO i, 0, x.len.pred:
    r[i] := simdSum(x[i])
  r
proc simdSum*(x: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[x.ncols,x.nrows,evalType(simdSum(x[0,0]))]
  forO i, 0, x.ncols.pred:
    forO j, 0, x.nrows.pred:
      r[i,j] := simdSum(x[i,j])
  r

#[
proc simdSum*(r: var auto; x: Mat2) {.inline.} =
  mixin nrows, ncols, trace, iadd
  assign(r, 0)
  for i in 0..<r.nrows:
    for j in 0..<r.ncols:
      r[i,j] := simdSum(x[i,j])
proc simdSum*(x: Mat1): auto {.inline,noInit.} =
  var t{.noInit.}: MatrixArray[x.nrows,x.ncols,evalType(simdSum(x[0,0]))]
  #static: echo "trace"
  simdSum(t, x)
  t
]#

when isMainModule:
  const nc = 3
  const ns = 4
  type
    Color[T] = object
      v: T
    Color2[T] = Color[T]
    Spin[T] = object
      v: T
    Spin2[T] = Spin[T]
    CVec = Color[VectorArray[nc,float]]
    CMat = Color[MatrixArray[nc,nc,float]]
    SCVec = Spin[VectorArray[ns,CVec]]
    SCMat = Spin[MatrixArray[ns,ns,CMat]]
  template `[]`(x: Color): untyped = x.v
  template `[]`(x: Color, i: int): untyped = x.v[i]
  template `[]`(x: Color, i,j: int): untyped = x.v[i,j]
  template load1(x: Color): untyped = x
  template assign(x: Color, y: Color2): untyped = assign(x[], y[])
  template assign(x: Color, y: SomeNumber): untyped = assign(x[], y)
  template assign(x: SomeNumber, y: Color2): untyped = assign(x, y[])
  template redot(x: Color, y: Color2): untyped = redot(x[], y[])
  template `*`(x: Color, y: Color2): untyped = `*`(x[], y[])
  template `[]`(x: Spin): untyped = x.v
  template `[]`(x: Spin, i: int): untyped = x.v[i]
  template `[]`(x: Spin, i,j: int): untyped = x.v[i,j]
  template assign(x: Spin, y: Spin2): untyped = assign(x[], y[])
  template assign(x: Spin, y: Color2): untyped = assign(x[], y)
  template assign(x: Spin, y: SomeNumber): untyped = assign(x[], y)
  template assign(x: SomeNumber, y: Spin2): untyped = assign(x, y[])
  template adj(x: Color): untyped =
    Color2[type(adj(x[]))](v: adj(x[]))

  var cv1,cv2:CVec
  var cm1,cm2:CMat
  var scv1,scv2:SCVec
  var scm1,scm2:SCMat

  proc test =
    var s:type(cv1[0])
    assign(cv1, 1)
    assign(s, cv1)
    echo s
    assign(cv2, cv1)
    assign(cm1, 1)
    assign(s, cm1)
    echo s
    assign(cm1, cv1)
    assign(cm2, cm1)
    echo cv2[0]
    echo cm2[0,0]

    assign(scv1, 2)
    echo "scv1: ", scv1[3][0]
    assign(scv1, cv1)
    echo "scv1: ", scv1[3][0]
    assign(scv1, 2)
    echo "scv1: ", scv1[3][0]
    assign(scv1, cv1)
    echo "scv1: ", scv1[3][0]

    var rd = redot(cm1,cm2)
    var rd2 = cm1*cm2
    var rd3 = cm1.adj*cm2
    var rd4 = trace(cm1.adj*cm2)
    echo rd
    echo rd2

    #neg(scm1, cv1)
    #add(scm1, scv1, cm1)
    #add(scm1, scm1, cm1)
    #echo scm1[0,0][0,0]

    #var vv:array[nc,float]
    #echo vv[0]
    #assign(asVector(vv), cv1)
    #echo vv[0]

  test()

  discard """
  template isScalar(x:MyScalar):untyped = true
  template isScalar(x:MyScalar2):untyped = true
  template isScalar(x:SomeNumber):untyped = true
  template needScalarOps(x:float):untyped = true
  #template `[]`(x:MyScalar, i:SomeInteger):untyped = float(x)
  #template `[]=`(x:MyScalar, i:SomeInteger, y:untyped):untyped =
  #  x = MyScalar(y)
  #template neg(x:MyScalar, r:var MyScalar) = r = -x
  #template `-`*(x:Myscalar):MyScalar2 = MyScalar2(-float(x))

  template isVector(x:MyVector):untyped = true
  template len(x:MyVector):untyped = n
  #template `[]`(x:MyVector, i:SomeInteger):untyped = (array[n,float](x))[i]
  #template `[]=`(x:MyVector, i:SomeInteger, y:untyped):untyped =
  #  (array[n,float](x))[i] = float(y)
  template isVector(x:MyVector2):untyped = true
  template len(x:MyVector2):untyped = n
  template `[]`(x:MyVector2, i:SomeInteger):untyped = (array[n,float](x))[i]
  template `[]=`(x:MyVector2, i:SomeInteger, y:untyped):untyped =
    (array[n,float](x))[i] = float(y)

  template isMatrix(x:MyMatrix):untyped = true
  template nrows(x:MyMatrix):untyped = n
  template ncols(x:MyMatrix):untyped = n
  template `[]`(x:MyMatrix, i,j:SomeInteger):untyped = x[i][j]
  template `[]=`(x:MyMatrix, i,j:SomeInteger, y:untyped):untyped = x[i][j] = y
  template isMatrix(x:MyMatrix2):untyped = true
  template nrows(x:MyMatrix2):untyped = n
  template ncols(x:MyMatrix2):untyped = n
  template `[]`(x:MyMatrix2, i,j:SomeInteger):untyped =
    (array[n,MyVector](x))[i][j]
  template `[]=`(x:MyMatrix2, i,j:SomeInteger, y:untyped):untyped =
    (array[n,MyVector](x))[i][j] = y

  template negTypes(X,R:typedesc):untyped =
    proc `-`*(x:X):R {.noInit,inline.} = neg(result,x)
  #negTypes(MyScalar, MyScalar2)
  negTypes(MyVector, MyVector2)
  negTypes(MyMatrix, MyMatrix2)
  template addTypes(X,Y,R:typedesc):untyped =
    proc `+`*(x:X,y:Y):R {.noInit,inline.} = add(result,x,y)
    proc `-`*(x:X,y:Y):R {.noInit,inline.} = sub(result,x,y)
  addTypes(MyVector, MyScalar, MyVector2)
  addTypes(MyScalar, MyVector, MyVector2)
  addTypes(MyVector, MyVector, MyVector2)
  addTypes(MyMatrix, MyScalar, MyMatrix2)
  addTypes(MyScalar, MyMatrix, MyMatrix2)
  addTypes(MyMatrix, MyVector, MyMatrix2)
  addTypes(MyVector, MyMatrix, MyMatrix2)
  addTypes(MyMatrix, MyMatrix, MyMatrix2)
  template mulTypes(X,Y,R:typedesc):untyped =
    proc `*`*(x:X,y:Y):R {.noInit,inline.} = mul(result,x,y)
  mulTypes(MyVector, MyScalar, MyVector2)
  mulTypes(MyScalar, MyVector, MyVector2)
  mulTypes(MyMatrix, MyScalar, MyMatrix2)
  mulTypes(MyScalar, MyMatrix, MyMatrix2)
  mulTypes(MyMatrix, MyVector, MyVector2)
  mulTypes(MyMatrix, MyMatrix, MyMatrix2)

  var
    s0 = MyScalar2(0.0)
    s1 = MyScalar(1.0)
    s2 = MyScalar(2.0)
    v0 = MyVector2([s2,s1,s1])
    v1 = MyVector([s1,s2,s1])
    v2 = MyVector([s1,s1,s2])
    m0 = MyMatrix2([v2,v1,v1])
    m1 = MyMatrix([v1,v2,v1])
    m2 = MyMatrix([v1,v1,v2])

  s0 = -s1
  v0 = -v1
  m0 = -m1

  s0 = s1 + s2
  v0 = v1 + v2
  m0 = m1 + m2
  v0 = v1 + s2
  v0 = s1 + v2
  m0 = m1 + s2
  m0 = s1 + m2
  m0 = m1 + v2
  m0 = v1 + m2

  s0 = s1 - s2
  v0 = v1 - v2
  m0 = m1 - m2
  v0 = v1 - s2
  v0 = s1 - v2
  m0 = m1 - s2
  m0 = s1 - m2
  m0 = m1 - v2
  m0 = v1 - m2

  s0 = s1 * s2
  echo($s0)
  m0 = m1 * m2
  echo($m0[0,0])
  v0 = v1 * s2
  echo($v0[0])
  v0 = s1 * v2
  echo($v0[0])
  m0 = m1 * s2
  echo($m0[0,0])
  m0 = s1 * m2
  echo($m0[0,0])
  v0 = m1 * v2
  echo($v0[0])
"""
