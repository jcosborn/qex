import physics/qcdTypes
export qcdTypes
import rng
export rng

# For IO with matching types; see read[T]/write[T]
<<<<<<< HEAD
template IOtype*(x:typedesc[SVec0]):untyped = float32
template IOtype*(x:typedesc[DVec0]):untyped = float64
template IOtype*(x:typedesc[SComplexV]):untyped = SComplex
template IOtype*(x:typedesc[DComplexV]):untyped = DComplex
template IOtype*(x:typedesc[SColorMatrixV]):untyped = SColorMatrix
template IOtype*(x:typedesc[DColorMatrixV]):untyped = DColorMatrix
template IOtype*(x:typedesc[RngMilc6]):untyped = RngMilc6
template IOtype*(x:typedesc[MRG32k3a]):untyped = MRG32k3a

# For IO with mis-matching types; see read[T]/write[T]
template IOtypeP*(x:typedesc[SVec0]):untyped = float64
template IOtypeP*(x:typedesc[DVec0]):untyped = float32
template IOtypeP*(x:typedesc[SComplexV]):untyped = DComplex
template IOtypeP*(x:typedesc[DComplexV]):untyped = SComplex
template IOtypeP*(x:typedesc[SColorMatrixV]):untyped = DColorMatrix
template IOtypeP*(x:typedesc[DColorMatrixV]):untyped = SColorMatrix
template IOtypeP*(x:typedesc[RngMilc6]):untyped = RngMilc6
template IOtypeP*(x:typedesc[MRG32k3a]):untyped = MRG32k3a
=======
#template IOtype*(x:typedesc[SVec0]):untyped = float32
#template IOtype*(x:typedesc[DVec0]):untyped = float64
#template IOtype*(x:typedesc[SComplexV]):untyped = SComplex
#template IOtype*(x:typedesc[DComplexV]):untyped = DComplex
#template IOtype*(x:typedesc[SColorMatrixV]):untyped = SColorMatrix
#template IOtype*(x:typedesc[DColorMatrixV]):untyped = DColorMatrix
#template IOtype*(x:typedesc[RngMilc6]):untyped = RngMilc6
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
#template IOtypeP*(x:typedesc[RngMilc6]):untyped = RngMilc6
template IOtypeP*[T](x:typedesc[T]):typedesc =
  mixin numberType, toDouble, toSingle
  when T.numberType is float32:
    eval(toDouble(IOtype(type T)))
  elif T.numberType is float64:
    eval(toSingle(IOtype(type T)))
  else:
    T
>>>>>>> upstream/devel
