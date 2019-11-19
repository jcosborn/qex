import physics/qcdTypes
export qcdTypes

# For IO with matching types; see read[T]/write[T]
template IOtype*(x:typedesc[SVec0]):untyped = float32
template IOtype*(x:typedesc[DVec0]):untyped = float64
template IOtype*(x:typedesc[SComplexV]):untyped = SComplex
template IOtype*(x:typedesc[DComplexV]):untyped = DComplex
template IOtype*(x:typedesc[SColorMatrixV]):untyped = SColorMatrix
template IOtype*(x:typedesc[DColorMatrixV]):untyped = DColorMatrix

# For IO with mis-matching types; see read[T]/write[T]
template IOtypeP*(x:typedesc[SVec0]):untyped = float64
template IOtypeP*(x:typedesc[DVec0]):untyped = float32
template IOtypeP*(x:typedesc[SComplexV]):untyped = DComplex
template IOtypeP*(x:typedesc[DComplexV]):untyped = SComplex
template IOtypeP*(x:typedesc[SColorMatrixV]):untyped = DColorMatrix
template IOtypeP*(x:typedesc[DColorMatrixV]):untyped = SColorMatrix
