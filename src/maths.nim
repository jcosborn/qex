import maths/complexNumbers
export complexNumbers
import maths/matrixConcept
export matrixConcept
import maths/matrixFunctions
export matrixFunctions


#[
import simd
when declared SimdD4:
  template `-`*(x: SimdD4, y: AsComplex): untyped =
    var r: DComplexV
    sub(r, x, y)
    r
  template `*`*(x: SimdD4, y: AsComplex): untyped =
    var r: DComplexV
    mul(r, x, y)
    r
  template sub*(r: var AsComplex, x: SimdD4, y: AsComplex2): untyped =
    r := asReal(x) - y

when declared SimdD8:
  template `-`*(x: SimdD8, y: AsComplex): untyped =
    #var r: DComplexV
    #sub(r, x, y)
    #r
    asReal(x) - y
  #template `*`*(x: SimdD8, y: AsComplex): untyped =
  #  var r: DComplexV
  #  mul(r, x, y)
  #  r
  template sub*(r: var AsComplex, x: SimdD8, y: AsComplex2): untyped =
    r := asReal(x) - y

when declared SimdD16:
  template `-`*(x: SimdD16, y: AsComplex): untyped =
    var r: DComplexV
    sub(r, x, y)
    r
  template `*`*(x: SimdD16, y: AsComplex): untyped =
    var r: DComplexV
    mul(r, x, y)
    r
]#
