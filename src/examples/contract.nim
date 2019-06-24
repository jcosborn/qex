import base, maths
import eigens/lapack

template toPtrInt32(x: int): ptr int32 =
  var t = x.int32
  addr t
template toPtrScomplex(x: int): ptr scomplex =
  var t = scomplex(re: x.float32, im: 0'f32)
  addr t
template toPtrDcomplex(x: int): ptr dcomplex =
  var t = dcomplex(re: x.float, im: 0'f64)
  addr t

template `&`(x: int): untyped = toPtrInt32(x)
template `&&`(x: int): untyped = toPtrScomplex(x)
template `&&&`(x: int): untyped = toPtrDcomplex(x)
template `&<`(x: ptr float32): untyped = cast[ptr scomplex](x)

type F32A = ptr UncheckedArray[float32]
type S32A = ptr UncheckedArray[ComplexType[float32]]
type F64A = ptr UncheckedArray[float64]
type S64A = ptr UncheckedArray[ComplexType[float64]]

# C = A x B'     ->  C^T = B^* x A^T
# C^T = A^* * B^T
proc cmatmul*(c,a,b: ptr scomplex, cr,cc,bc: int) =
  cgemm("C","N", &cc,&cr,&bc, &&1, b,&bc, a,&bc, &&0, c,&cc)
proc cmatmul*(c,a,b: ptr dcomplex, cr,cc,bc: int) =
  zgemm("C","N", &cc,&cr,&bc, &&&1, b,&bc, a,&bc, &&&0, c,&cc)

template cmatmul*(c,a,b: S32A, cr,cc,bc: int) =
  template `/`(x: S32A): untyped = cast[ptr scomplex](x)
  cmatmul(/c, /a, /b, cr, cc, bc)

template cmatmul*(c,a,b: S64A, cr,cc,bc: int) =
  template `/`(x: S64A): untyped = cast[ptr dcomplex](x)
  cmatmul(/c, /a, /b, cr, cc, bc)

#template cmatmul*(c,a,b: ptr float32, cr,cc,bc: int) =
#  cmatmul(&<c, &<a, &<b, cr, cc, bc)
#template cmatmul*(c,a,b: F32A, cr,cc,bc: int) =
#  cmatmul(&<c, &<a, &<b, cr, cc, bc)

#template adj(x: scomplex): untyped = scomplex(re: x.re, im: -x.im)
# C[cr,cc] = A[cr,bc] x B[cc,bc]'
proc cmatmulX*(c,a,b: S32A, cr,cc,bc: int) =
  for i in 0..<cr:
    for j in 0..<cc:
      var t = a[i*bc] * adj(b[j*bc])
      for k in 1..<bc:
        t += a[i*bc+k] * adj(b[j*bc+k])
      c[i*cc+j] = t
#template cmatmulX*(c,a,b: F32A, cr,cc,bc: int) =
#  cmatmulX(&<c, &<a, &<b, cr, cc, bc)

when isMainModule:
  import maths/complexType
  type Cmplx = ComplexType[float32]
  converter toPtrScomplex(x: ptr UncheckedArray[Cmplx]): ptr scomplex =
    cast[ptr scomplex](x)
  var m = 10
  var n = 20
  var k = 3
  var c = cast[ptr UncheckedArray[Cmplx]](alloc(m*n*sizeof(Cmplx)))
  var d = cast[ptr UncheckedArray[Cmplx]](alloc(m*n*sizeof(Cmplx)))
  var a = cast[ptr UncheckedArray[Cmplx]](alloc(m*k*sizeof(Cmplx)))
  var b = cast[ptr UncheckedArray[Cmplx]](alloc(n*k*sizeof(Cmplx)))

  for i in 0..<(m*k):
    a[i] := newComplex(i+1,i+2)
  for i in 0..<(n*k):
    b[i] := newComplex(i+2,i+1)

  cmatmulX(c, a, b, m, n, k)
  cmatmul(d, a, b, m, n, k)
  var s = 0.0
  for i in 0..<(m*n):
    let t = d[i] - c[i]
    s += t.norm2
  echo s
