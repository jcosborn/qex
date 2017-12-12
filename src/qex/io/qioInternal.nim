import qex/physics/qcdTypes
export qcdTypes

# For IO with matching types; see read[T]/write[T]
template IOtype*(x:typedesc[SColorMatrixV]):typedesc = SColorMatrix
template IOtype*(x:typedesc[DColorMatrixV]):typedesc = DColorMatrix
# For IO with mis-matching types; see read[T]/write[T]
template IOtypeP*(x:typedesc[SColorMatrixV]):typedesc = DColorMatrix
template IOtypeP*(x:typedesc[DColorMatrixV]):typedesc = SColorMatrix
