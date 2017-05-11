import maths/complexConcept
export complexConcept
import maths/matrixConcept
export matrixConcept
import maths/matrixFunctions
export matrixFunctions

import simd
template `-`*(x: SimdD8, y: AsComplex): untyped =
  var r: DComplexV
  sub(r, x, y)
  r
template `*`*(x: SimdD8, y: AsComplex): untyped =
  var r: DComplexV
  mul(r, x, y)
  r
