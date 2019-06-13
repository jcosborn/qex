import eigens/lapack

template toPtrInt32(x: int): ptr int32 =
  var t = x.int32
  addr t
template toPtrScomplex(x: int): ptr scomplex =
  var t = scomplex(re: x.float32, im: 0'f32)
  addr t

template `&`(x: int): untyped = toPtrInt32(x)
template `&&`(x: int): untyped = toPtrScomplex(x)
template `&<`(x: ptr float32): untyped = cast[ptr scomplex](x)

# C^T = A^* * B^T
proc cmatmul(c,a,b: ptr scomplex, cr,cc,bc: int) =
  cgemm("C","N", &cr,&cc,&bc, &&1, a,&bc, b,&bc, &&0, c,&cc)
template cmatmul*(c,a,b: ptr float32, cr,cc,bc: int) =
  cmatmul(&<c, &<a, &<b, cr, cc, bc)

when isMainModule:
  import maths/complexType
  type Cmplx = ComplexType[float32]
  converter toPtrScomplex(x: ptr UncheckedArray[Cmplx]): ptr scomplex =
    cast[ptr scomplex](x)
  var m = 10
  var n = 20
  var k = 3
  var c = cast[ptr UncheckedArray[Cmplx]](alloc(m*n*sizeof(Cmplx)))
  var a = cast[ptr UncheckedArray[Cmplx]](alloc(m*k*sizeof(Cmplx)))
  var b = cast[ptr UncheckedArray[Cmplx]](alloc(n*k*sizeof(Cmplx)))

  cmatmul(c, a, b, m, n, k)
