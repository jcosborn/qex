import physics/qcdTypes
export qcdTypes

# For IO with matching types; see read[T]/write[T]
template IOtype*(x:typedesc[SComplexV]):typedesc = SComplex
template IOtype*(x:typedesc[DComplexV]):typedesc = DComplex
template IOtype*(x:typedesc[SColorMatrixV]):typedesc = SColorMatrix
template IOtype*(x:typedesc[DColorMatrixV]):typedesc = DColorMatrix

# For IO with mis-matching types; see read[T]/write[T]
template IOtypeP*(x:typedesc[SComplexV]):typedesc = DComplex
template IOtypeP*(x:typedesc[DComplexV]):typedesc = SComplex
template IOtypeP*(x:typedesc[SColorMatrixV]):typedesc = DColorMatrix
template IOtypeP*(x:typedesc[DColorMatrixV]):typedesc = SColorMatrix
