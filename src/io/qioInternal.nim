import physics/qcdTypes
export qcdTypes
import rng
export rng
import typetraits

# For IO with matching types; see read[T]/write[T]
#template IOtype*(x:typedesc[SVec0]):untyped = float32
#template IOtype*(x:typedesc[DVec0]):untyped = float64
#template IOtype*(x:typedesc[SComplexV]):untyped = SComplex
#template IOtype*(x:typedesc[DComplexV]):untyped = DComplex
#template IOtype*(x:typedesc[SColorMatrixV]):untyped = SColorMatrix
#template IOtype*(x:typedesc[DColorMatrixV]):untyped = DColorMatrix
template IOtype*[T](x:typedesc[T]):typedesc =
  mixin has, index
  when T.has Simd:
    eval(T.index(asSimd(int)))
  else:
    T

# For IO with mis-matching types; see read[T]/write[T]
#template IOtypeP*(x:typedesc[SVec0]):untyped = float64
#template IOtypeP*(x:typedesc[DVec0]):untyped = float32
#template IOtypeP*(x:typedesc[SComplexV]):untyped = DComplex
#template IOtypeP*(x:typedesc[DComplexV]):untyped = SComplex
#template IOtypeP*(x:typedesc[SColorMatrixV]):untyped = DColorMatrix
#template IOtypeP*(x:typedesc[DColorMatrixV]):untyped = SColorMatrix
template IOtypeP*[T](x:typedesc[T]):typedesc =
  mixin numberType, toDouble, toSingle
  when T.numberType is float32:
    eval(toDouble(IOtype(type T)))
  elif T.numberType is float64:
    eval(toSingle(IOtype(type T)))
  else:
    T

template IOnameDefault*[T](x:typedesc[T]):string =
  "QDP_" & T.name
template IOname*[T](x:typedesc[T]):string =
  T.IOnameDefault
#template IOname*[N:static int](x:typedesc[Color[MatrixArray[N,N,float]]]):string =
#  "QDP_F" & $N & "_ColorMatrix"
#template IOname*[N:static int](x:typedesc[Color[MatrixArray[N,N,DComplex]]]):string =
#  "QDP_D" & $N & "_ColorMatrix"
template IOname*[N:static int,T](x:typedesc[Color[MatrixArray[N,N,T]]]):string =
  when T is SComplex:
    "QDP_F" & $N & "_ColorMatrix"
  elif T is DComplex:
    "QDP_D" & $N & "_ColorMatrix"
  else:
    IOnameDefault T
