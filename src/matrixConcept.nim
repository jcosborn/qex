import stdUtils
import metaUtils
import macros
import basicOps
import types
import wrapperTypes
import complexConcept
# unary ops: assign(=),neg(-),iadd(+=),isub(-=),imul(*=),idiv(/=)
# binary ops: add(+),sub(-),mul(*),divd(/),imadd(+=*),imsub(-=*)
# ternary ops: madd(*+),msub(*-),nmadd(-*+),nmsub(-*-)
# assign: trace,dot,outer
# wrap: conj,adj,transpose
# norm2,det,lu,qr,svd,eig
# sqrt,rsqrt,exp,log,groupProject,groupCheck
# mapX(f,x,r),mapXY(f,x,y,r)

template createAsType2(t,c:untyped):untyped =
  mixin `[]`
  makeWrapper(t, c)
  template `[]`*(x:t; i:SomeInteger):expr = x[][i]
  template `[]`*(x:t; i,j:SomeInteger):expr = 
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
  template `[]=`*(x:t; i,j:SomeInteger; y:untyped):expr =
    x[][i,j] = y
  template len*(x:t):expr = x[].len
  template nrows*(x:t):expr = x[].nrows
  template ncols*(x:t):expr = x[].ncols
  #template mvLevel*(x:t):expr =
  #  mixin mvLevel
  #  mvLevel(x[])
template createAsType(t:untyped):untyped = createAsType2(`As t`, `as t`)

createAsType(Scalar)
createAsType(VarScalar)
createAsType(Vector)
createAsType(VarVector)
createAsType(Matrix)
createAsType(VarMatrix)

template makeDeclare(s:untyped):untyped {.dirty.} =
  template `declare s`*(t:typedesc):untyped {.dirty.} =
    template `declared s`*(y:t):expr {.dirty.} = true
  template `is s`*(x:typed):expr {.dirty.} =
    when compiles(`declared s`(x)):
      `declared s`(x)
    else:
      false
makeDeclare(Scalar)
makeDeclare(Matrix)
makeDeclare(Vector)
declareScalar(AsScalar)
declareScalar(AsVarScalar)
declareVector(AsVector)
declareVector(AsVarVector)
declareMatrix(AsMatrix)
declareMatrix(AsVarMatrix)
template deref(x:typed):expr =
  when type(x) is AsScalar|AsVarScalar:
    x[]
  else:
    x

type
  Vec1* = concept x
    #mixin isVector
    x.isVector
  Vec2* = concept x
    #mixin isVector
    x.isVector
  Vec3* = concept x
    #mixin isVector
    x.isVector
  Mat1* = concept x
    #mixin isMatrix
    x.isMatrix
  Mat2* = concept x
    #mixin isMatrix
    x.isMatrix
  Mat3* = concept x
    #mixin isMatrix
    x.isMatrix
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
  VarSca1* = var Sca1
  VarVec1* = var Vec1
  VarMat1* = var Mat1
  VarMV1* = AsVarMatrix | AsVarVector
  VarAny* = var any #| AsVarMatrix
  VectorArray*[I:static[int],T] = AsVector[array[I,T]]
  MatrixArrayObj*[I,J:static[int],T] = array[I,array[J,T]]
  MatrixArray*[I,J:static[int],T] = AsMatrix[MatrixArrayObj[I,J,T]]
  MatrixRowObj*[T] = object
    row:int
    mat:ptr T
  MatrixRow*[T] = AsVector[MatrixRowObj[T]]
  #MatrixCol*[T] = tuple[col:int,mat:ptr T]
  #MatrixDiag*[T] = tuple[diag:int,mat:ptr T]

template nrows*(x:MatrixArrayObj):expr = x.I
template ncols*(x:MatrixArrayObj):expr = x.J
template `[]`*(x:MatrixArrayObj; i,j:int):expr = x[i][j]
template `[]=`*(x:MatrixArrayObj; i,j:int, y:untyped):untyped = x[i][j] = y
template numNumbers*(x:AsVector):expr = x.len*numNumbers(x[0])
template numNumbers*(x:AsMatrix):expr = x.nrows*x.ncols*numNumbers(x[0])
#template `[]`*(x:array; i,j:int):expr = x[i][j]
#template `[]=`*(x:array; i,j:int, y:untyped):untyped = x[i][j] = y

template len*(x:MatrixRowObj):expr = x.mat[].ncols
template `[]`*(x:MatrixRowObj; i:int):expr = x.mat[][x.row,i]
template `[]=`*(x:MatrixRowObj; i:int; y:untyped):expr = x.mat[][x.row,i] = y

#template isVector(x:Row):expr = true
#template isVector(x:Col):expr = true
#template mvLevel(x:Sca1):expr = -1

macro getConst(x:typed):auto =
  #echo x.treerepr
  #result = newLit(3)
  result = newLit(x.intVal)

template load1*(x:Vec1):expr =
  var r{.noInit.}:VectorArray[x.len,type(load1(x[0]))]
  assign(r, x)
  r
template load1*(x:Mat1):expr =
  var r{.noInit.}:MatrixArray[getConst(x.nrows),getConst(x.ncols),type(load1(x[0,0]))]
  assign(r, x)
  r

template cfor*(i,r0,r1,b:untyped):untyped =
  block:
    var i = r0
    while i <= r1:
      b
      inc(i)
#template forO*(i,r0,r1,b:untyped):untyped = #cfor(i,r0,r1,b)
#  var i:int
#  for ii{.gensym.} in r0..r1:
#    i = ii
#    b
macro forO*(i,r0,r1,b:untyped):auto =
  #echo b.repr
  result = quote do:
    for `i` in `r0`..`r1`:
      `b`
#template forO*(i,r0,r1,b:untyped):untyped = forStatic(i,r0,r1,b)

template row*(x:AsVector; i:int):expr = x
proc row*(x:AsMatrix; i:int):auto {.inline,noInit.} =
  const nc = x.ncols
  var r{.noInit.}:VectorArray[nc,type(load1(x[0,0]))]
  for j in 0..<nc:
    assign(r[j], x[i,j])
  r
  #return asVector(MatrixRowObj[type(x)](row:i,mat:unsafeAddr(x)))
template setRow*(r:AsVector; x:AsVector; i:int):untyped =
  assign(r, x)
proc setRow*(r:var AsMatrix; x:AsVector; i:int) {.inline.} =
  const nc = r.ncols
  for j in 0..<nc:
    assign(r[i,j], x[j])
template column*(x:AsVector; i:int):expr = x
proc column*(x:AsMatrix; i:int):auto {.inline,noInit.} =
  const nr = x.nrows
  var r{.noInit.}:VectorArray[nr,type(x[0,0])]
  for j in 0..<nr:
    assign(r[j], x[j,i])
  r
template setColumn*(r:AsVector; x:AsVector; i:int):untyped =
  assign(r, x)
proc setColumn*(r:var AsMatrix; x:AsVector; i:int) {.inline.} =
  const nr = r.nrows
  for j in 0..<nr:
    assign(r[j,i], x[j])

proc `$`*(x:Vec1):string =
  mixin `$`
  result = $(x[0])
  forO i, 1, x.len-1:
    result.add "," & $(x[i])
proc `$`*(x:Mat1):string =
  mixin `$`
  result = ""
  forO i, 0, x.nrows-1:
    result.add $(x[i,0])
    forO j, 1, x.ncols-1:
      result.add "," & $(x[i,j])
    if i<x.nrows-1: result.add "\n"

template makeLevel1(f,s1,t1,s2,t2:untyped):untyped =
  proc f*(r:t1, x:t2) {.inline.} =
    `f s1 s2`(r, x)
#macro makeLevel1(f,s1,t1,s2,t2:untyped):auto =
  #mixin mvLevel
#  let f12 = ident($f & $s1 & $s2)
#  result = quote do:
    #template `f`*(r:`t1`, x:`t2`):untyped =
#    proc `f`*(r:`t1`, x:`t2`) {.inline.} =
    #when mvLevel(r) > mvLevel(x):
    #  `f s1 S`(r, x)
    #elif mvLevel(r) < mvLevel(x):
    #  `f S s2`(r, x)
    #else:
      #`f s1 s2`(r, x)
#      `f12`(r, x)
macro makeLevel2(f,s1,t1,s2,t2,s3,t3:untyped):auto =
  #let f1SS = ident($f & $s1 & "SS")
  #let fSS3 = ident($f & "SS" & $s3)
  #let f1S3 = ident($f & $s1 & "S" & $s3)
  #let fS2S = ident($f & "S" & $s2 & "S")
  #let fS23 = ident($f & "S" & $s2 & $s3)
  #let f12S = ident($f & $s1 & $s2 & "S")
  let f123 = ident($f & $s1 & $s2 & $s3)
  result = quote do:
    #mixin mvLevel
    #template `f`*(r:`t1`; x:`t2`; y:`t3`):untyped =
    proc `f`*(r:`t1`; x:`t2`; y:`t3`) {.inline.} =
      #echoTyped(x.mvLevel)
      #when r.mvLevel > x.mvLevel:
      #  when r.mvLevel > y.mvLevel:
      #    `f1SS`(r, x, y)
      #  elif mvLevel(r) < mvLevel(y):
      #    `fSS3`(r, x, y)
      #  else:
      #    `f1S3`(r, x, y)
      #elif mvLevel(r) < mvLevel(x):
      #  when mvLevel(x) > mvLevel(y):
      #    `fS2S`(r, x, y)
      #  elif mvLevel(x) < mvLevel(y):
      #    `fSS3`(r, x, y)
      #  else:
      #    `fS23`(r, x, y)
      #else:
      #  when mvLevel(r) > mvLevel(y):
      #    `f12S`(r, x, y)
      #  elif mvLevel(r) < mvLevel(y):
      #    `fSS3`(r, x, y)
      #  else:
          `f123`(r, x, y)
  #echo result.repr

template assignIadd(x,y:untyped):untyped = iadd(x,y)
template negIadd(x,y:untyped):untyped = isub(x,y)
template iaddIadd(x,y:untyped):untyped = iadd(x,y)
template isubIadd(x,y:untyped):untyped = isub(x,y)
template makeMap1(op:untyped):untyped =
  #mixin op, `op Iadd`
  template `op SV`*(rr:typed; xx:typed):untyped =
    subst(r,rr,x,xx,i,_):
      #tmpvar(tr, r)
      op(r, x[0])
      forO i, 1, <x.len:
        `op Iadd`(r, x[i])
  template `op SM`*(rr:typed; xx:typed):untyped =
    subst(r,rr,x,xx,i,_):
      assert(x.nrows == x.ncols)
      op(r, x[0,0])
      forO i, 1, <x.nrows:
        `op Iadd`(r, x[i,i])
  #proc `op VS`*(r:Vany; x:any) =
  template `op VS`*(rr:typed; xx:typed):untyped =
    subst(r,rr,x,xx,tx,_,i,_):
      load(tx, deref(x))
      forO i, 0, <r.len:
        op(r[i], tx)
  #proc `op VV`*(r:Vany; x:any) {.inline.} =
  template `op VV`*(rr:typed; xx:typed):untyped =
    subst(r,rr,x,xx,i,_):
      assert(r.len == x.len)
      forO i, 0, <r.len:
        op(r[i], x[i])
  #proc `op MS`*(r:Vany; x:any) {.inline.} =
  template `op MS`*(rr:typed; xx:typed):untyped =
    subst(r,rr,x,xx,tx,_,i,_,j,_):
      #assert(r.nrows == r.ncols)
      load(tx, x)
      forO i, 0, <r.nrows:
        forO j, 0, <r.ncols:
          if i == j:
            op(r[i,j], tx)
          else:
            op(r[i,j], 0)
  template `op MV`*(rr:typed; xx:typed):untyped =
    subst(r,rr,x,xx,i,_,j,_):
      assert(r.nrows == x.len)
      assert(r.ncols == x.len)
      forO i, 0, <r.nrows:
        forO j, 0, <r.ncols:
          if i == j:
            op(r[i,j], x[i])
          else:
            op(r[i,j], 0)
  template `op MM`*(rr:typed; xx:typed):untyped =
    subst(r,rr,x,xx,i,_,j,_):
      assert(r.nrows == x.nrows)
      assert(r.ncols == x.ncols)
      forO i, 0, <r.nrows:
        forO j, 0, <r.ncols:
          op(r[i,j], x[i,j])
  makeLevel1(op, S, VarSca1, V, Vec2)
  makeLevel1(op, S, VarSca1, M, Mat2)
  makeLevel1(op, V, VarVec1, S, Sca2)
  makeLevel1(op, V, VarVec1, V, Vec2)
  makeLevel1(op, V, AsVarVector, S, Sca2)
  makeLevel1(op, V, AsVarVector, V, Vec2)
  makeLevel1(op, M, VarMat1, S, Sca2)
  makeLevel1(op, M, VarMat1, V, Vec2)
  makeLevel1(op, M, VarMat1, M, Mat2)
  makeLevel1(op, M, AsVarMatrix, S, Sca2)
  makeLevel1(op, M, AsVarMatrix, V, Vec2)
  makeLevel1(op, M, AsVarMatrix, M, Mat2)

makeMap1(assign)
makeMap1(neg)
makeMap1(iadd)
makeMap1(isub)
template `:=`*(x:VarVec1; y:SomeNumber) = assign(x, y)
template `:=`*(x:AsVarVector; y:SomeNumber) = assign(x, y)
template `:=`*(x:VarVec1; y:Vec2) = assign(x, y)
template `:=`*(x:AsVarVector; y:Vec2) = assign(x, y)
template `:=`*(x:VarMat1; y:SomeNumber) = assign(x, y)
template `:=`*(x:AsVarMatrix; y:SomeNumber) = assign(x, y)
template `:=`*(x:VarMat1; y:Vec2) = assign(x, y)
template `:=`*(x:AsVarMatrix; y:Vec2) = assign(x, y)
template `:=`*(x:VarMat1; y:Mat2) = assign(x, y)
template `:=`*(x:AsVarMatrix; y:Mat2) = assign(x, y)

template `+=`*(x:VarVec1; y:Vec2) = iadd(x, y)

template imulMS*(r:typed; x:typed):untyped =
  load(tx, x)
  forO i, 0, <r.nrows:
    forO j, 0, <r.ncols:
      imul(r[i,j], tx)
proc imul*(r:VarMat1; x:Sca2) {.inline.} = imulMS(r, x)
proc imul*(r:AsVarMatrix; x:Sca2) {.inline.} = imulMS(r, x)
proc `*=`*(r:var MV1; x:Sca2) {.inline.} = imul(r, x)
proc `*=`*(r:VarMV1; x:Sca2) {.inline.} = imul(r, x)

template makeMap2(op:untyped):untyped =
  template `op VVS`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,ty,_,i,_):
      assert(r.len == x.len)
      load(ty, y)
      forO i, 0, <r.len:
        op(r[i], x[i], ty)
  template `op VSV`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,tx,_,i,_):
      assert(r.len == y.len)
      load(tx, x)
      forO i, 0, <r.len:
        op(r[i], tx, y[i])
  template `op VVV`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,i,_):
      assert(r.len == y.len)
      assert(r.len == x.len)
      forO i, 0, <r.len:
        op(r[i], x[i], y[i])
  template `op MSS`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,tx,_,ty,_,i,_,j,_):
      assert(r.nrows == r.ncols)
      load(tx, x)
      load(ty, y)
      forO i, 0, <r.nrows:
        forO j, 0, <r.ncols:
          if i == j:
            op(r[i,j], tx, ty)
          else:
            op(r[i,j], 0, 0)
  #proc `op MVS`*(r:Vany; x,y:any) {.inline.} =
  template `op MVS`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,ty,_,i,_,j,_):
      assert(r.nrows == x.len)
      assert(r.ncols == x.len)
      load(ty, y)
      forO i, 0, <r.nrows:
        forO j, 0, <r.ncols:
          if i == j:
            op(r[i,j], x[i], ty)
          else:
            op(r[i,j], 0, 0)
  template `op MSV`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,tx,_,i,_,j,_):
      assert(r.nrows == y.len)
      assert(r.ncols == y.len)
      load(tx, x)
      forO i, 0, <r.nrows:
        forO j, 0, <r.ncols:
          if i == j:
            op(r[i,j], tx, y[i])
          else:
            op(r[i,j], 0, 0)
  template `op MVV`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,i,_,j,_):
      assert(r.nrows == x.len)
      assert(r.ncols == x.len)
      assert(x.len == y.len)
      forO i, 0, <r.nrows:
        forO j, 0, <r.ncols:
          if i == j:
            op(r[i,j], x[i], y[i])
          else:
            op(r[i,j], 0, 0)
  template `op MMS`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,ty,_,i,_,j,_):
      assert(r.nrows == r.ncols)
      assert(r.nrows == x.nrows)
      assert(r.ncols == x.ncols)
      load(ty, y)
      forO i, 0, <r.nrows:
        forO j, 0, <r.ncols:
          if i == j:
            op(r[i,j], x[i,j], ty)
          else:
            op(r[i,j], x[i,j], 0)
  template `op MSM`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,tx,_,i,_,j,_):
      assert(r.nrows == r.ncols)
      assert(r.nrows == y.nrows)
      assert(r.ncols == y.ncols)
      load(tx, x)
      forO i, 0, <r.nrows:
        forO j, 0, <r.ncols:
          if i == j:
            op(r[i,j], tx, y[i,j])
          else:
            op(r[i,j], 0, y[i,j])
  template `op MMV`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,i,_,j,_):
      assert(r.nrows == r.ncols)
      assert(r.nrows == x.nrows)
      assert(r.ncols == x.ncols)
      assert(r.nrows == y.len)
      forO i, 0, <r.nrows:
        forO j, 0, <r.ncols:
          if i == j:
            op(r[i,j], x[i,j], y[i])
          else:
            op(r[i,j], x[i,j], 0)
  template `op MVM`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,i,_,j,_):
      assert(r.nrows == r.ncols)
      assert(r.nrows == x.len)
      assert(r.nrows == y.nrows)
      assert(r.ncols == y.ncols)
      forO i, 0, <r.nrows:
        forO j, 0, <r.ncols:
          if i == j:
            op(r[i,j], x[i], y[i,j])
          else:
            op(r[i,j], 0, y[i,j])
  template `op MMM`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,i,_,j,_):
      assert(r.nrows == x.nrows)
      assert(r.ncols == x.ncols)
      assert(r.nrows == y.nrows)
      assert(r.ncols == y.ncols)
      forO i, 0, <r.nrows:
        forO j, 0, <r.ncols:
          op(r[i,j], x[i,j], y[i,j])

  makeLevel2(op, V, VarVec1, V, Vec2, S, Sca3)
  makeLevel2(op, V, VarVec1, S, Sca2, V, Vec3)
  makeLevel2(op, V, VarVec1, V, Vec2, V, Vec3)
  makeLevel2(op, M, VarMat1, S, Sca2, S, Sca3)
  makeLevel2(op, M, VarMat1, V, Vec2, S, Sca3)
  makeLevel2(op, M, VarMat1, S, Sca2, V, Vec3)
  makeLevel2(op, M, VarMat1, V, Vec2, V, Vec3)
  makeLevel2(op, M, VarMat1, M, Mat2, S, Sca3)
  makeLevel2(op, M, VarMat1, S, Sca2, M, Mat3)
  makeLevel2(op, M, VarMat1, M, Mat2, V, Vec3)
  makeLevel2(op, M, VarMat1, V, Vec2, M, Mat3)
  makeLevel2(op, M, VarMat1, M, Mat2, M, Mat3)

makeMap2(add)
makeMap2(sub)

proc `+`*(x:Vec1; y:Vec2):auto {.inline.} =
  assert(x.len==y.len)
  const n = x.len
  var r{.noInit.}:VectorArray[n,type(x[0]+y[0])]
  add(r, x, y)
  r
proc `+`*(x:Sca1; y:Mat2):auto {.inline.} =
  const nr = y.nrows
  const nc = y.ncols
  var r{.noInit.}:MatrixArray[nr,nc,type(x+y[0,0])]
  add(r, x, y)
  r
proc `+`*(x:Mat1; y:Mat2):auto {.inline.} =
  assert(x.nrows==y.nrows)
  assert(x.ncols==y.ncols)
  const nr = x.nrows
  const nc = x.ncols
  var r{.noInit.}:MatrixArray[nr,nc,type(x[0,0]+y[0,0])]
  add(r, x, y)
  r
proc `-`*(x:Vec1; y:Vec2):auto {.inline.} =
  assert(x.len==y.len)
  const n = x.len
  var r{.noInit.}:VectorArray[n,type(x[0]-y[0])]
  sub(r, x, y)
  r
proc `-`*(x:Sca1; y:Mat2):auto {.inline.} =
  const nr = y.nrows
  const nc = y.ncols
  var r{.noInit.}:MatrixArray[nr,nc,type(x-y[0,0])]
  sub(r, x, y)
  r

template mulSVV*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tr,_,i,_):
    mixin mul, imadd
    assert(x.len == y.len)
    tmpvar(tr, r)
    mul(tr, x[0], y[0])
    forO i, 1, <x.len:
      imadd(tr, x[i], y[i])
    assign(r, tr)
template mulVVS*(r:typed; x,y:typed):untyped =
  mixin mul
  assert(r.len == x.len)
  load(ty, y)
  forO i, 0, <r.len:
    mul(r[i], x[i], ty)
template mulVSV*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tx,_,i,_):
    mixin mul
    assert(r.len == y.len)
    load(tx, x)
    forO i, 0, <r.len:
      mul(r[i], tx, y[i])
proc mul*(r:var Vec1; x:Sca2; y:Vec3) {.inline.} = mulVSV(r, x, y)
proc mulMMS*(r:any; x,y:any) {.inline.} =
  #mixin mul
  assert(r.nrows == x.nrows)
  assert(r.ncols == x.ncols)
  #load(ty, y)
  forO i, 0, <r.nrows:
    forO j, 0, <r.ncols:
      #echo isComplex(r[i,j])
      mul(r[i,j], x[i,j], y)
template mulMSM*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tx,_,i,_,j,_):
    assert(r.nrows == y.nrows)
    assert(r.ncols == y.ncols)
    load(tx, x)
    forO i, 0, <r.nrows:
      forO j, 0, <r.ncols:
        mul(r[i,j], tx, y[i,j])
#proc mulVMV*(r:Vany; x,y:any) {.inline.} =
template mulVMV*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tr,_,ty,_,i,_,j,_,ty0r,_,ty0i,_,tyjr,_,tyji,_):
    mixin nrows, ncols, mul, imadd, assign, load1
    assert(x.nrows == r.len)
    assert(x.ncols == y.len)
    when true:
    #when false:
      tmpvar(tr, r)
      #var tr{.noInit.}:array[r.len,type(load1(r[0]))]
      load(ty0r, y[0].re)
      forO i, 0, <x.nrows:
        mulCCR(tr[i], x[i,0], ty0r)
      load(ty0i, y[0].im)
      forO i, 0, <x.nrows:
        imaddCCI(tr[i], x[i,0], ty0i)
      forO j, 1, <x.ncols:
        load(tyjr, y[j].re)
        forO i, 0, <x.nrows:
          imaddCCR(tr[i], x[i,j], tyjr)
        load(tyji, y[j].im)
        forO i, 0, <x.nrows:
          imaddCCI(tr[i], x[i,j], tyji)
      assign(r, tr)
      #forO i, 0, <r.len: assign(r[i], tr[i])
    else:
      tmpvar(tr, r)
      block:
        load(ty, y[0])
        forO i, 0, <x.nrows:
          mul(tr[i], x[i,0], ty)
      forO j, 1, <x.ncols:
        block:
          load(ty, y[j])
          forO i, 0, <x.nrows:
            imadd(tr[i], x[i,j], ty)
      assign(r, tr)
template mulMMM*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tr,_,i,_,j,_,k,_,txi0r,_,txi0i,_,txikr,_,txiki,_):
    mixin nrows, ncols, mul, imadd, assign, load1
    assert(r.nrows == x.nrows)
    assert(r.ncols == y.ncols)
    assert(x.ncols == y.nrows)
    #tmpvar(tr, r)
    var tr{.noInit.}:VectorArray[getConst(r.ncols),type(x[0,0]*y[0,0])]
    forO i, 0, <r.nrows:
      load(txi0r, x[i,0].re)
      forO j, 0, <r.ncols:
        mulCRC(tr[j], txi0r, y[0,j])
      load(txi0i, x[i,0].im)
      forO j, 0, <r.ncols:
        imaddCIC(tr[j], txi0i, y[0,j])
      forO k, 1, <x.ncols:
        load(txikr, x[i,k].re)
        forO j, 0, <r.ncols:
          imaddCRC(tr[j], txikr, y[k,j])
        load(txiki, x[i,k].im)
        forO j, 0, <r.ncols:
          imaddCIC(tr[j], txiki, y[k,j])
    #assign(r, tr)
      forO j, 0, <r.ncols:
        assign(r[i,j], tr[j])
#[
  assert(x.nrows == r.nrows)
  assert(x.ncols == y.nrows)
  assert(r.ncols == y.ncols)
  mixin mul, imadd
  forO i, 0, <r.nrows:
    var tr{.noInit.}:VectorArray[getConst(r.ncols),type(x[0,0]*y[0,0])]
    load(txi0, x[i,0])
    forO j, 0, <r.ncols:
      mul(tr[j], txi0, y[0,j])
    for k in 1..<x.ncols:
      load(txik, x[i,k])
      forO j, 0, <r.ncols:
        imadd(tr[j], txik, y[k,j])
    forO j, 0, <r.ncols:
      assign(r[i,j], tr[j])
]#
makeLevel2(mul, V, VarVec1, V, Vec2, S, Sca3)
#makeLevel2(mul, V, VarVec1, S, Sca2, V, Vec3)
#makeLevel2(op, S, Sca1, V, Vec2, V, Vec3)
makeLevel2(mul, M, VarMat1, M, Mat2, S, Sca3)
makeLevel2(mul, M, VarMat1, S, Sca2, M, Mat3)
makeLevel2(mul, V, VarVec1, M, Mat2, V, Vec3)
#makeLevel2(op, V, Vec1, V, Vec2, M, Mat3)
makeLevel2(mul, M, VarMat1, M, Mat2, M, Mat3)
#makeLevel2(op, M, Mat1, S, Sca2, S, Sca3)
#makeLevel2(op, M, Mat1, V, Vec2, S, Sca3)
#makeLevel2(op, M, Mat1, S, Sca2, V, Vec3)
#makeLevel2(op, M, Mat1, V, Vec2, V, Vec3)
#makeLevel2(op, M, Mat1, V, Vec2, M, Mat3)

template `*`*(x:Sca1; y:Vec2):expr =
  const n = y.len
  var r{.noInit.}:VectorArray[n,type(x*y[0])]
  mulVSV(r, x, y)
  r
#proc `*`*(x:Sca1; y:Vec2):auto {.inline.} =
#  const n = y.len
#  var r{.noInit.}:VectorArray[n,type(x*y[0])]
#  mulVSV(r, x, y)
#  r
proc `*`*(x:Vec1; y:Sca2):auto {.inline.} =
  const n = x.len
  var r{.noInit.}:VectorArray[n,type(x[0]*y)]
  mul(r, x, y)
  r
proc `*`*(x:Sca1; y:Mat2):auto {.inline.} =
  const nr = y.nrows
  const nc = y.ncols
  var r{.noInit.}:MatrixArray[nr,nc,type(x*y[0,0])]
  mul(r, x, y)
  r
proc `*`*(x:Mat1; y:Sca2):auto {.inline.} =
  const nr = x.nrows
  const nc = x.ncols
  var r{.noInit.}:MatrixArray[nr,nc,type(x[0,0]*y)]
  mul(r, x, y)
  r
template `*`*(x:Mat1; y:Vec2):expr =
  assert(x.ncols == y.len)
  const n = getConst(x.nrows)
  var r{.noInit.}:VectorArray[n,type(x[0,0]*y[0])]
  mul(r, x, y)
  r
#proc `*`*(x:Mat1; y:Vec2):auto {.inline.} =
#  assert(x.ncols == y.len)
#  const n = x.nrows
#  var r{.noInit.}:VectorArray[n,type(x[0,0]*y[0])]
#  mul(r, x, y)
#  r
proc `*`*(x:Mat1; y:Mat2):auto {.inline.} =
  assert(x.ncols == y.nrows)
  const nr = x.nrows
  const nc = y.ncols
  var r{.noInit.}:MatrixArray[nr,nc,type(x[0,0]*y[0,0])]
  mul(r, x, y)
  r

template imaddSVV*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tr,_,i,_):
    mixin imadd, assign
    assert(x.len == y.len)
    load(tr, r)
    forO i, 0, <x.len:
      imadd(tr, x[i], y[i])
    assign(r, tr)
template imaddVMV*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tr,_,ty,_,i,_,j,_,tyjr,_,tyji,_):
    mixin nrows, ncols, mul, imadd, assign, load1
    assert(x.nrows == r.len)
    assert(x.ncols == y.len)
    when true:
    #when false:
      load(tr, r)
      #var tr{.noInit.}:array[r.len,type(load1(r[0]))]
      #forO i, 0, <r.len: assign(tr[i], r[i])
      forO j, 0, <x.ncols:
        load(tyjr, y[j].re)
        forO i, 0, <x.nrows:
          imaddCCR(tr[i], x[i,j], tyjr)
        load(tyji, y[j].im)
        forO i, 0, <x.nrows:
          imaddCCI(tr[i], x[i,j], tyji)
      assign(r, tr)
      #forO i, 0, <r.len: assign(r[i], tr[i])
    else:
      load(tr, r)
      forO j, 0, <x.ncols:
        load(tyr, asReal(y[j].re))
        forO i, 0, <x.nrows:
          imadd(tr[i], x[i,j], tyr)
        load(tyi, asImag(y[j].im))
        forO i, 0, <x.nrows:
          imadd(tr[i], x[i,j], tyi)
      assign(r, tr)
proc imadd*(r:VarVec1; x:Mat2; y:Vec3) {.inline.} = imaddVMV(r, x, y)
proc imadd*(r:AsVarVector; x:Mat2; y:Vec3) {.inline.} = imaddVMV(r, x, y)

template imaddMMM*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tr,_,ty,_,i,_,j,_,k,_,tyjr,_,tyji,_):
    mixin nrows, ncols, mul, imadd, assign, load1
    assert(r.nrows == x.nrows)
    assert(r.ncols == y.ncols)
    assert(x.ncols == y.nrows)
    #when true:
    when false:
      load(tr, r)
      forO i, 0, <r.nrows:
        forO k, 0, <x.ncols:
          load(txikr, x[i,k].re)
          forO j, 0, <r.ncols:
            imaddCRC(tr[i,j], txikr, y[k,j])
          load(txiki, x[i,k].im)
          forO j, 0, <r.ncols:
            imaddCIC(tr[i,j], txiki, y[k,j])
      assign(r, tr)
    else:
      forO i, 0, <r.nrows:
        var tr{.noInit.}:VectorArray[getConst(r.ncols),type(x[0,0]*y[0,0])]
        forO j, 0, <r.ncols:
          assign(tr[j], r[i,j])
        forO k, 0, <x.ncols:
          load(txikr, x[i,k].re)
          forO j, 0, <r.ncols:
            imaddCRC(tr[j], txikr, y[k,j])
          load(txiki, x[i,k].im)
          forO j, 0, <r.ncols:
            imaddCIC(tr[j], txiki, y[k,j])
        forO j, 0, <r.ncols:
          assign(r[i,j], tr[j])
proc imadd*(r:VarMat1; x:Mat2; y:Mat3) {.inline.} = imaddMMM(r, x, y)

template imsubVSV*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tx,_,i,_):
    mixin imsub
    assert(r.len == y.len)
    load(tx, x)
    forO i, 0, <r.len:
      imsub(r[i], x, y[i])
proc imsub*(r:var Vec1; x:Sca2; y:Vec3) {.inline.} = imsubVSV(r, x, y)

template imsubVMV*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tr,_,ty,_,i,_,j,_):
    mixin imsub
    assert(x.nrows == r.len)
    assert(x.ncols == y.len)
    when true:
      load(tr, r)
      forO j, 0, <x.ncols:
        load(ty, y[j])
        forO i, 0, <x.nrows:
          imsub(tr[i], x[i,j], ty)
      assign(r, tr)
    else:
      load(tr, r)
      forO j, 0, <x.ncols:
        load(tyr, asReal(y[j].re))
        forO i, 0, <x.nrows:
          imsub(tr[i], x[i,j], tyr)
        load(tyi, asImag(y[j].im))
        forO i, 0, <x.nrows:
          imsub(tr[i], x[i,j], tyi)
      assign(r, tr)
proc imsub*(r:VarVec1; x:Mat2; y:Vec3) {.inline.} = imsubVMV(r, x, y)
proc imsub*(r:AsVarVector; x:Mat2; y:Vec3) {.inline.} = imsubVMV(r, x, y)

template imsubMMM*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tr,_,ty,_,i,_,j,_,k,_,txikr,_,txiki,_):
    mixin nrows, ncols, imsubCRC, imsubCIC, assign, load1
    assert(r.nrows == x.nrows)
    assert(r.ncols == y.ncols)
    assert(x.ncols == y.nrows)
    load(tr, r)
    forO i, 0, <r.nrows:
      forO k, 0, <x.ncols:
        load(txikr, x[i,k].re)
        forO j, 0, <r.ncols:
          imsubCRC(tr[i,j], txikr, y[k,j])
        load(txiki, x[i,k].im)
        forO j, 0, <r.ncols:
          imsubCIC(tr[i,j], txiki, y[k,j])
    assign(r, tr)
proc imsub*(r:VarMat1; x:Mat2; y:Mat3) {.inline.} = imsubMMM(r, x, y)

template msubVSVV*(rr:typed; xx,yy,zz:typed):untyped =
  subst(r,rr,x,xx,y,yy,z,zz,i,_):
    mixin msub
    assert(r.len == y.len)
    assert(r.len == z.len)
    forO i, 0, <r.len:
      msub(r[i], x, y[i], z[i])
proc msub*(r:VarVec1; x:any; y:Vec2; z:Vec3) {.inline.} = msubVSVV(r,x,y,z)

proc trace*(r:VarSca1; x:Mat2) {.inline.} =
  mixin nrows, ncols, trace, iadd
  let n = min(x.nrows, x.ncols)
  assign(r, 0)
  for i in 0..<n:
    let t = trace(x[i,i])
    iadd(r, t)
proc trace*(x:Mat1):auto {.inline,noInit.} =
  var t:type(trace(x[0,0]))
  trace(t, x)
  t

proc inorm2*(r:VarAny; x:Vec2) {.inline.} =
  mixin inorm2
  for i in 0..<x.len:
    #echo r
    inorm2(r, x[i])
proc inorm2*(r:VarAny; x:Mat2) {.inline.} =
  mixin nrows, ncols, inorm2
  for i in 0..<x.nrows:
    for j in 0..<x.ncols:
      inorm2(r, x[i,j])
proc norm2*(r:VarAny; x:Vec2) {.inline.} =
  mixin norm2, iadd
  assign(r, 0)
  for i in 0..<x.len:
    var t{.noInit.}:type(r)
    norm2(t, x[i])
    iadd(r, t)
proc norm2*(r:VarAny; x:Mat2) {.inline.} =
  mixin nrows, ncols, norm2, iadd
  assign(r, 0)
  for i in 0..<x.nrows:
    for j in 0..<x.ncols:
      var t{.noInit.}:type(r)
      norm2(t, x[i,j])
      iadd(r, t)
proc norm2*(x:Vec1):auto {.inline,noInit.} =
  var t{.noInit.}:type(norm2(x[0]))
  norm2(t, x)
  t
proc norm2*(x:Mat1):auto {.inline,noInit.} =
  var t{.noInit.}:type(norm2(x[0,0]))
  norm2(t, x)
  t
proc norm2X*(x:Vec1):auto {.inline,noInit.} =
  var t{.noInit.}:type(norm2X(x[0]))
  norm2(t, x)
  t
proc norm2X*(x:Mat1):auto {.inline,noInit.} =
  var t{.noInit.}:type(norm2X(x[0,0]))
  norm2(t, x)
  t

proc dot*(r:Sca1; x:Vec2; y:Vec3) {.inline.} =
  mulSVV(r, x.adj, y)
proc idot*(r:var Sca1; x:Vec2; y:Vec3) {.inline.} =
  imaddSVV(r, x.adj, y)
proc iredot*(r:var Sca1; x:Vec2; y:Vec3) {.inline.} =
  imaddSVV(r, x.adj, y)

  
when isMainModule:
  const nc = 3
  const ns = 4
  type
    CVec = VectorArray[nc,float]
    CMat = MatrixArray[nc,nc,float]
    SCVec = VectorArray[ns,CVec]
    SCMat = MatrixArray[ns,ns,CMat]

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
    assign(scv1, asScalar(cv1))
    echo "scv1: ", scv1[3][0]

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
