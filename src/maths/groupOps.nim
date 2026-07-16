import base
import complexNumbers
import matrixConcept
import types
import matrixFunctions
import algorithms/numdiff

#[
anti-Hermitian generators, with normalization
	tr{T^a T^a} = -1/2

Basic relations with anti-symmetric f^{abc} and symmetric d^{abc}
	T^a T^b = 1/2 [ (f^{abc} + i d^{abc}) T^c - 1/n δ^{ab} ]

Convention:
	SU...: SU(N) group, special unitary matrix
	su...: su(n) algebra, traceless anti-Hermitian matrix
]#

const
  sqrt_1_3 = 0.57735026918962576451    # sqrt(1/3)

proc trysqrt(n2, n: static int): auto =
  when n*n == n2:
    n
  elif n*n < n2:
    trysqrt(n2, n+1)
  else:
    static: error("trysqrt n2=" & $n2 & " n=" & $n)
proc ncFromDim(d: static int): auto =
  trysqrt(d+1, 1)

proc matchGroupVec(label: static string, g: Mat1, v: Vec1) =
  const nc = g.nrows
  const dim = v.len
  when dim != nc*nc-1:
    static: error(label & " v[" & $v.len & "] g[" & $g.nrows & "," & $g.ncols & "]")
  when evalType(v[0]) is not evalType(g[0,0].re):
    static: error(label & " wrong type g:" & g.getType & " v:" & v.getType)

proc matchGroupAd(label: string, g: Mat1, a: Mat2) =
  const nc = g.nrows
  const dim = a.nrows
  when dim != nc*nc-1:
    static: error(label & " a[" & $a.nrows & "," & $a.ncols & "] g[" & $g.nrows & "," & $g.ncols & "]")
  when evalType(a[0,0]) is not evalType(g[0,0].re):
    static: error(label & " wrong type g:" & g.getType & " a:" & a.getType)

proc matchMatVec(label: string, m: Mat1, v: Vec1) =
  when m.nrows != m.ncols or m.nrows != v.len:
    static: error(label & " v[" & $v.len & "] m[" & $m.nrows & "," & $m.ncols & "]")
  when evalType(v[0]) is not evalType(m[0,0]):
    static: error(label & " wrong type m:" & m.getType & " v:" & v.getType)

proc suToVec*(r: var Vec1, m: Mat1) =
  ## Only for anti-Hermitian m, or in su(N) algebra.  Return real numbers X^a, such that X^a T^a = X - 1/N tr(X).
  ## Convention: tr{T^a T^a} = -1/2
  ## X^a = - 2 tr[T^a X]
  matchGroupVec("suToVec", m, r)
  const nc = m.nrows
  const c = -2.0
  when nc==3:
    # assuming anti-Hermitian m, uses the upper triangle
    let
      m00i = m[0,0].im
      m01r = m[0,1].re
      m01i = m[0,1].im
      m02r = m[0,2].re
      m02i = m[0,2].im
      m11i = m[1,1].im
      m12r = m[1,2].re
      m12i = m[1,2].im
      m22i = m[2,2].im
    r[0] := c*m01i
    r[1] := c*m01r
    r[2] := m11i-m00i
    r[3] := c*m02i
    r[4] := c*m02r
    r[5] := c*m12i
    r[6] := c*m12r
    r[7] := sqrt_1_3*(2.0*m22i-m11i-m00i)
  else:
    static: error"suToVec unimplemented for n!=3"

template suToVec*(m: Mat1): auto =
  const dim = m.nrows*m.nrows-1
  type R = evalType(m[0,0].re)
  type V = VectorArray[dim,R]
  var r {.noinit.}: V
  r.suToVec(m)
  r

func suToVec_mat*(m: Mat1): auto {.noinit.} =
  ## Only for anti-Hermitian m, or in su(N) algebra.  Return real numbers X^a, such that X^a T^a = X - 1/N tr(X).
  ## Convention: tr{T^a T^a} = -1/2
  ## X^a = - 2 tr[T^a X]
  ## Implemented with direct matrix ops
  const c = -2.0
  const dim = m.nrows*m.nrows-1
  var r {.noinit.}: VectorArray[dim, evalType(m[0,0].re)]
  when m.nrows==3:
    for a in 0..<dim:
      r[a] := c * trace(su3gen[a] * m).re
    r
  else:
    static: error"suToVec unimplemented for n!=3"

proc suFromVec*(r: var Mat1, v: Vec1) =
  ## Return su(N), X = X^a T^a
  matchGroupVec("suFromVec", r, v)
  const dim = v.len
  const c = -0.5
  when dim==8:
    let
      r01i = c*v[0]
      r01r = c*v[1]
      v2 = v[2]
      r02i = c*v[3]
      r02r = c*v[4]
      r12i = c*v[5]
      r12r = c*v[6]
      r2i = sqrt_1_3*v[7]
    r[0,0].re := 0
    r[0,0].im := c*(r2i+v2)
    r[0,1].re := r01r
    r[0,1].im := r01i
    r[0,2].re := r02r
    r[0,2].im := r02i
    r[1,0].re := -r01r
    r[1,0].im := r01i
    r[1,1].re := 0
    r[1,1].im := c*(r2i-v2)
    r[1,2].re := r12r
    r[1,2].im := r12i
    r[2,0].re := -r02r
    r[2,0].im := r02i
    r[2,1].re := -r12r
    r[2,1].im := r12i
    r[2,2].re := 0
    r[2,2].im := r2i
  else:
    static: error"suFromVec unimplemented for d!=8"

template suFromVec*(v: Vec1): auto =
  const nc = ncFromDim(v.len)
  type C = ComplexType[evalType(v[0])]
  type SU = MatrixArray[nc, nc, C]
  var r {.noinit.}: SU
  r.suFromVec(v)
  r

func suFromVec_mat*(v: Vec1): auto {.noinit.} =
  ## Return su(N), X = X^a T^a
  ## Implemented with direct matrix ops
  when v.len==8:
    var r : MatrixArray[3, 3, ComplexType[evalType(v[0])]]
    for a in 0..<v.len:
      r += v[a] * su3gen[a]
    r
  else:
    static: error"suFromVec unimplemented for d!=8"

proc sufabc*(r: var Mat1, v: Vec1) =
  ## returns f^{abc} v[c]
  ## [T^a, T^b] = f^abc T^c
  matchMatVec("sufabc", r, v)
  const d = v.len
  when d==8:
    const
      f012 = 1.0
      f036 = 0.5
      f045 = -f036
      f135 = f036
      f146 = f036
      f234 = f036
      f256 = f045
      f347 = 0.86602540378443864676    # sqrt(3/4)
      f567 = f347
    let
      v0 = v[0]
      v1 = v[1]
      v2 = v[2]
      v3 = v[3]
      v4 = v[4]
      v5 = v[5]
      v6 = v[6]
      v7 = v[7]
    r[0,0] := 0
    r[0,1] :=   f012  * v2
    r[0,2] := (-f012) * v1
    r[0,3] :=   f036  * v6
    r[0,4] :=   f045  * v5
    r[0,5] := (-f045) * v4
    r[0,6] := (-f036) * v3
    r[0,7] := 0
    r[1,0] := -r[0,1]
    r[1,1] := 0
    r[1,2] :=   f012  * v0
    r[1,3] :=   f135  * v5
    r[1,4] :=   f146  * v6
    r[1,5] := (-f135) * v3
    r[1,6] := (-f146) * v4
    r[1,7] := 0
    r[2,0] := -r[0,2]
    r[2,1] := -r[1,2]
    r[2,2] := 0
    r[2,3] :=   f234  * v4
    r[2,4] := (-f234) * v3
    r[2,5] :=   f256  * v6
    r[2,6] := (-f256) * v5
    r[2,7] := 0
    r[3,0] := -r[0,3]
    r[3,1] := -r[1,3]
    r[3,2] := -r[2,3]
    r[3,3] := 0
    r[3,4] :=   f347  * v7 + f234 * v2
    r[3,5] :=   f135  * v1
    r[3,6] :=   f036  * v0
    r[3,7] := (-f347) * v4
    r[4,0] := -r[0,4]
    r[4,1] := -r[1,4]
    r[4,2] := -r[2,4]
    r[4,3] := -r[3,4]
    r[4,4] := 0
    r[4,5] :=   f045  * v0
    r[4,6] :=   f146  * v1
    r[4,7] :=   f347  * v3
    r[5,0] := -r[0,5]
    r[5,1] := -r[1,5]
    r[5,2] := -r[2,5]
    r[5,3] := -r[3,5]
    r[5,4] := -r[4,5]
    r[5,5] := 0
    r[5,6] :=   f567  * v7 + f256 * v2
    r[5,7] := (-f567) * v6
    r[6,0] := -r[0,6]
    r[6,1] := -r[1,6]
    r[6,2] := -r[2,6]
    r[6,3] := -r[3,6]
    r[6,4] := -r[4,6]
    r[6,5] := -r[5,6]
    r[6,6] := 0
    r[6,7] :=   f567 * v5
    r[7,0] := 0
    r[7,1] := 0
    r[7,2] := 0
    r[7,3] := -r[3,7]
    r[7,4] := -r[4,7]
    r[7,5] := -r[5,7]
    r[7,6] := -r[6,7]
    r[7,7] := 0
  else:
    static: error"sufabc unimplemented for d!=8"

template sufabc*(v: Vec1): auto =
  const d = v.len
  var r {.noinit.}: MatrixArray[d, d, evalType(v[0])]
  r.sufabc(v)
  r

proc sudabc*(r: var Mat1, v: Vec1) =
  ## returns returns d^abc v[c]
  ## {T^a,T^b} = -1/n δ^ab + i d^abc T^c
  # NOTE: negative sign of what's on wikipedia, because of anti-Hermitian T
  matchMatVec("sufabc", r, v)
  const d = v.len
  when d==8:
    const
      d007 = -sqrt_1_3
      d035 = -0.5
      d046 = d035
      d117 = d007
      d136 = -d035
      d145 = d035
      d227 = d007
      d233 = d035
      d244 = d035
      d255 = d136
      d266 = d136
      d337 = 0.5*sqrt_1_3
      d447 = d337
      d557 = d337
      d667 = d337
      d777 = -d007
    let
      v0 = v[0]
      v1 = v[1]
      v2 = v[2]
      v3 = v[3]
      v4 = v[4]
      v5 = v[5]
      v6 = v[6]
      v7 = v[7]
    r[0,0] := d007*v7
    r[0,1] := 0
    r[0,2] := 0
    r[0,3] := d035*v5
    r[0,4] := d046*v6
    r[0,5] := d035*v3
    r[0,6] := d046*v4
    r[0,7] := d007*v0
    r[1,0] := 0
    r[1,1] := d117*v7
    r[1,2] := 0
    r[1,3] := d136*v6
    r[1,4] := d145*v5
    r[1,5] := d145*v4
    r[1,6] := d136*v3
    r[1,7] := d117*v1
    r[2,0] := 0
    r[2,1] := 0
    r[2,2] := d227*v7
    r[2,3] := d233*v3
    r[2,4] := d244*v4
    r[2,5] := d255*v5
    r[2,6] := d266*v6
    r[2,7] := d227*v2
    r[3,0] := r[0,3]
    r[3,1] := r[1,3]
    r[3,2] := r[2,3]
    r[3,3] := d337*v7+d233*v2
    r[3,4] := 0
    r[3,5] := d035*v0
    r[3,6] := d136*v1
    r[3,7] := d337*v3
    r[4,0] := r[0,4]
    r[4,1] := r[1,4]
    r[4,2] := r[2,4]
    r[4,3] := 0
    r[4,4] := d447*v7+d244*v2
    r[4,5] := d145*v1
    r[4,6] := d046*v0
    r[4,7] := d447*v4
    r[5,0] := r[0,5]
    r[5,1] := r[1,5]
    r[5,2] := r[2,5]
    r[5,3] := r[3,5]
    r[5,4] := r[4,5]
    r[5,5] := d557*v7+d255*v2
    r[5,6] := 0
    r[5,7] := d557*v5
    r[6,0] := r[0,6]
    r[6,1] := r[1,6]
    r[6,2] := r[2,6]
    r[6,3] := r[3,6]
    r[6,4] := r[4,6]
    r[6,5] := 0
    r[6,6] := d667*v7+d266*v2
    r[6,7] := d667*v6
    r[7,0] := r[0,7]
    r[7,1] := r[1,7]
    r[7,2] := r[2,7]
    r[7,3] := r[3,7]
    r[7,4] := r[4,7]
    r[7,5] := r[5,7]
    r[7,6] := r[6,7]
    r[7,7] := d777*v7
  else:
    static: error"sudabc unimplemented for d!=8"

template sudabc*(v: Vec1): auto =
  const d = v.len
  var r {.noinit.}: MatrixArray[d, d, evalType(v[0])]
  r.sudabc(v)
  r

proc suad*(r: var Mat1, v: var Vec1, x: Mat2) =
  # adX^{ab} = - f^{abc} X^c = f^{abc} 2 tr(X T^c) = 2 tr(X [T^a, T^b])
  # Input x must be in su(n) algebra.
  v.suToVec x
  v := -v
  r.sufabc v

proc suad*(r: var Mat1, x: Mat2) =
  var v = suToVec x
  v := -v
  r.sufabc v

template suad*(x: Mat1): auto =
  sufabc(-suToVec(x))

proc suadApply*(r:var Mat1, v: var Vec1, adx: Mat2, y: Mat3) =
  # adX(Y) = [X, Y]
  # adX(T^b) = T^a adX^{ab} = - T^a f^{abc} X^c = X^c f^{cba} T^a = X^c [T^c, T^b] = [X, T^b]
  # adX(Y) = T^a adX^{ab} Y^b = T^a adX^{ab} (-2) tr{T^b Y}
  # Input y must be in su(n) algebra.
  v.suToVec y
  v := adx*v
  r.suFromVec v

proc suadApply*(r:var Mat1, adx: Mat2, y: Mat3) =
  var v = suToVec y
  v := adx*v
  r.suFromVec v

template suadApply*(adx: Mat1, y: Mat2): auto =
  let v = adx*suToVec(y)
  suFromVec(v)

proc SUAd*(r: var Mat1, x: Mat2) =
  # X T^c X† = AdX T^c = T^b AdX^bc
  # Input x must be in SU(N) group.
  # AdX^bc = - 2 tr[T^b X T^c X†] = - 2 tr[T^c X† T^b X]
  #        = 2 tr[(T^b X) (X T^c)†]
  # Vij_{kl} = 1 if (i,j)=(k,l) else 0
  # tr[Vij X Vkl X†] = xjk conj(xil)
  # tr[Vji X Vlk X†] = xil conj(xjk) = conj(tr[Vij X Vlk X†])
  # T0 = -0.5i (V01 + V10)
  # T1 =  0.5 (-V01 + V10)
  # T2 = -0.5i (V00 - V11)
  # T3 = -0.5i (V02 + V20)
  # T4 =  0.5 (-V02 + V20)
  # T5 = -0.5i (V12 + V21)
  # T6 =  0.5 (-V12 + V21)
  # T7 = -0.5 sqrt(1/3) i (V00 + V11 - 2 V22)
  # -2 tr[T0 X T0 X†] = 0.5 ( x10 conj(x01) + x11 conj(x00) + complex conjugate ) = re(x10 conj(x01) + x11 conj(x00))
  # -2 tr[T0 X T1 X†] = 0.5i ( -x10 conj(x01) + x11 conj(x00) - complex conjugate ) = im(x10 conj(x01) - x11 conj(x00))
  # -2 tr[T0 X T2 X†] = 0.5 ( x10 conj(x00) - x11 conj(x01) + c.c. ) = re(x10 conj(x00) - x11 conj(x01))
  matchGroupAd("SUAd", x, r)
  when x.nrows==3:
    let
      x00 = x[0,0]
      x01 = x[0,1]
      x02 = x[0,2]
      x10 = x[1,0]
      x11 = x[1,1]
      x12 = x[1,2]
      x20 = x[2,0]
      x21 = x[2,1]
      x22 = x[2,2]
    let
      x00s = redot(x00,x00)
      x01s = redot(x01,x01)
      x02s = redot(x02,x02)
      x10s = redot(x10,x10)
      x11s = redot(x11,x11)
      x12s = redot(x12,x12)
      x20s = redot(x20,x20)
      x21s = redot(x21,x21)
      x22s = redot(x22,x22)
    let
      x00c01 = x00*adj(x01)
      x00c02 = x00*adj(x02)
      x00c10 = x00*adj(x10)
      x00c11 = x00*adj(x11)
      x00c12 = x00*adj(x12)
      x00c20 = x00*adj(x20)
      x00c21 = x00*adj(x21)
      x00c22 = x00*adj(x22)
    let
      x01c02 = x01*adj(x02)
      x01c10 = x01*adj(x10)
      x01c11 = x01*adj(x11)
      x01c12 = x01*adj(x12)
      x01c20 = x01*adj(x20)
      x01c21 = x01*adj(x21)
      x01c22 = x01*adj(x22)
    let
      x02c10 = x02*adj(x10)
      x02c11 = x02*adj(x11)
      x02c12 = x02*adj(x12)
      x02c20 = x02*adj(x20)
      x02c21 = x02*adj(x21)
      x02c22 = x02*adj(x22)
    let
      x10c11 = x10*adj(x11)
      x10c12 = x10*adj(x12)
      x10c20 = x10*adj(x20)
      x10c21 = x10*adj(x21)
      x10c22 = x10*adj(x22)
    let
      x11c12 = x11*adj(x12)
      x11c20 = x11*adj(x20)
      x11c21 = x11*adj(x21)
      x11c22 = x11*adj(x22)
    let
      x12c20 = x12*adj(x20)
      x12c21 = x12*adj(x21)
      x12c22 = x12*adj(x22)
    let
      x20c21 = x20*adj(x21)
      x20c22 = x20*adj(x22)
    let
      x21c22 = x21*adj(x22)
    r[0,0] :=   x00c11.re + x01c10.re
    r[0,1] :=   x00c11.im - x01c10.im
    r[0,2] :=   x00c10.re - x01c11.re
    r[0,3] :=   x00c12.re + x02c10.re
    r[0,4] :=   x00c12.im - x02c10.im
    r[0,5] :=   x01c12.re + x02c11.re
    r[0,6] :=   x01c12.im - x02c11.im
    r[0,7] :=  (x00c10.re + x01c11.re - 2.0*x02c12.re)*sqrt_1_3
    r[1,0] := -(x00c11.im + x01c10.im)
    r[1,1] :=   x00c11.re - x01c10.re
    r[1,2] :=   x01c11.im - x00c10.im
    r[1,3] := -(x00c12.im + x02c10.im)
    r[1,4] :=   x00c12.re - x02c10.re
    r[1,5] := -(x01c12.im + x02c11.im)
    r[1,6] :=   x01c12.re - x02c11.re
    r[1,7] :=  (2.0*x02c12.im - (x00c10.im + x01c11.im))*sqrt_1_3
    r[2,0] :=   x00c01.re - x10c11.re
    r[2,1] :=   x00c01.im - x10c11.im
    r[2,2] :=  (x00s - x01s - x10s + x11s)*0.5
    r[2,3] :=   x00c02.re - x10c12.re
    r[2,4] :=   x00c02.im - x10c12.im
    r[2,5] :=   x01c02.re - x11c12.re
    r[2,6] :=   x01c02.im - x11c12.im
    r[2,7] :=  (0.5*(x00s + x01s - (x10s + x11s)) - x02s + x12s)*sqrt_1_3
    r[3,0] :=   x00c21.re + x01c20.re
    r[3,1] :=   x00c21.im - x01c20.im
    r[3,2] :=   x00c20.re - x01c21.re
    r[3,3] :=   x00c22.re + x02c20.re
    r[3,4] :=   x00c22.im - x02c20.im
    r[3,5] :=   x01c22.re + x02c21.re
    r[3,6] :=   x01c22.im - x02c21.im
    r[3,7] :=  (x00c20.re + x01c21.re - 2.0*x02c22.re)*sqrt_1_3
    r[4,0] := -(x00c21.im + x01c20.im)
    r[4,1] :=   x00c21.re - x01c20.re
    r[4,2] :=   x01c21.im - x00c20.im
    r[4,3] := -(x00c22.im + x02c20.im)
    r[4,4] :=   x00c22.re - x02c20.re
    r[4,5] := -(x01c22.im + x02c21.im)
    r[4,6] :=   x01c22.re - x02c21.re
    r[4,7] :=  (2.0*x02c22.im - (x00c20.im + x01c21.im))*sqrt_1_3
    r[5,0] :=   x10c21.re + x11c20.re
    r[5,1] :=   x10c21.im - x11c20.im
    r[5,2] :=   x10c20.re - x11c21.re
    r[5,3] :=   x10c22.re + x12c20.re
    r[5,4] :=   x10c22.im - x12c20.im
    r[5,5] :=   x11c22.re + x12c21.re
    r[5,6] :=   x11c22.im - x12c21.im
    r[5,7] :=  (x10c20.re + x11c21.re - 2.0*x12c22.re)*sqrt_1_3
    r[6,0] := -(x10c21.im + x11c20.im)
    r[6,1] :=   x10c21.re - x11c20.re
    r[6,2] :=   x11c21.im - x10c20.im
    r[6,3] := -(x10c22.im + x12c20.im)
    r[6,4] :=   x10c22.re - x12c20.re
    r[6,5] := -(x11c22.im + x12c21.im)
    r[6,6] :=   x11c22.re - x12c21.re
    r[6,7] :=  (2.0*x12c22.im - (x10c20.im + x11c21.im))*sqrt_1_3
    r[7,0] :=  (x00c01.re + x10c11.re - 2.0*x20c21.re)*sqrt_1_3
    r[7,1] :=  (x00c01.im + x10c11.im - 2.0*x20c21.im)*sqrt_1_3
    r[7,2] :=  (0.5*(x00s + x10s - (x01s + x11s)) - x20s + x21s)*sqrt_1_3
    r[7,3] :=  (x00c02.re + x10c12.re - 2.0*x20c22.re)*sqrt_1_3
    r[7,4] :=  (x00c02.im + x10c12.im - 2.0*x20c22.im)*sqrt_1_3
    r[7,5] :=  (x01c02.re + x11c12.re - 2.0*x21c22.re)*sqrt_1_3
    r[7,6] :=  (x01c02.im + x11c12.im - 2.0*x21c22.im)*sqrt_1_3
    r[7,7] :=  (0.5*(x00s + x11s + x01s + x10s) + 2.0*x22s - (x02s + x20s + x12s + x21s))*0.33333333333333333333
  else:
    static: error"SUAd unimplemented for n!=3"

template SUAd*(x: Mat1): auto =
  const dim = x.nrows*x.nrows-1
  var r {.noinit.}: MatrixArray[dim, dim, evalType(x[0,0].re)]
  r.SUAd x
  r

func SUAd_mat*(x: Mat1): auto {.noinit.} =
  # X T^c X† = AdX T^c = T^b AdX^bc
  # Input x must be in SU(N) group.
  # AdX^bc = - 2 tr[T^b X T^c X†] = - 2 tr[T^c X† T^b X]
  ## Implemented with direct matrix ops
  const dim = x.nrows*x.nrows-1
  var r {.noinit.}: MatrixArray[dim, dim, evalType(x[0,0].re)]
  when x.nrows==3:
    for b in 0..<dim:
      for c in 0..<dim:
        r[b,c] := (-2.0) * trace(su3gen[b] * x * su3gen[c] * x.adj).re
    r
  else:
    static: error"SUAd_mat unimplemented for n!=3"

func mkGellMann*[T]: VectorArray[8, MatrixArray[3,3,ComplexType[T]]] =
  result[0][0,1].re = 1.0
  result[0][1,0].re = 1.0
  result[1][0,1].im = -1.0
  result[1][1,0].im = 1.0
  result[2][0,0].re = 1.0
  result[2][1,1].re = -1.0
  result[3][0,2].re = 1.0
  result[3][2,0].re = 1.0
  result[4][0,2].im = -1.0
  result[4][2,0].im = 1.0
  result[5][1,2].re = 1.0
  result[5][2,1].re = 1.0
  result[6][1,2].im = -1.0
  result[6][2,1].im = 1.0
  result[7][0,0].re = sqrt_1_3
  result[7][1,1].re = sqrt_1_3
  result[7][2,2].re = -2.0*sqrt_1_3

const gellMann* = mkGellMann[float]()

func mkSU3Gen*[T]: VectorArray[8, MatrixArray[3,3,ComplexType[T]]] =
  # let m_i_2 = newImag(-0.5)
  let m_i_2 = newComplex(0.0, -0.5)
  for i in 0..<result.len:
    result[i] = m_i_2 * gellMann[i]

const su3gen* = mkSU3Gen[float]()

func sugen*(nc:static int):auto {.noinit,inline.} =
  when nc==3:
    su3gen
  else:
    static: error"sugen unimplemented for n!=3"

proc diffProjectTAH*(r:var Mat1, m: Mat2, p: Mat3) =
  ## r_{ac} = ∂_c p^a = ∂_c projectTAH(m)^a = - tr[T^a (T^c M + M† T^c)] = -2 ReTr[T^a T^c M]
  #[
    P^a = -2 tr[T^a {- T^d tr[T^d (M - M†)]}]
        = - tr[T^a (M - M†)]
        = - ∂_a tr[M + M†]
    ∂_c P^a = - tr[T^a (T^c M + M† T^c)]
            = - 1/2 tr[{T^a,T^c} (M+M†) + [T^a,T^c] (M-M†)]
            = - 1/2 tr[d^acb T^b i (M+M†) - 1/N δ^ac (M+M†) + f^acb T^b (M-M†)]
            = - 1/2 { d^acb tr[T^b i(M+M†)] - 1/N δ^ac tr(M+M†) - f^acb F^b }
            = - 1/2 { d^acb tr[T^b i(M+M†)] - 1/N δ^ac tr(M+M†) + adF^ac }
    Note:
        T^a T^b = 1/2 {(f^abc + i d^abc) T^c - 1/N δ^ab}
  ]#
  const nc = m.nrows
  const ii = newComplex(0.0, 0.25)
  var t = m + m.adj
  let trMs = trace(t).re/(2.0*nc)
  t *= ii
  let v = suToVec(t)
  t := (-0.5) * p
  r.sudabc(v)
  # r += suad(t) + trMs    # FIXME: error in matrixOps.nim:/op MMS/+/op.* 0/
  r += trMs
  r += suad(t)

proc diffCrossProjectTAH*(r: var Mat1, Adx: Mat2, dp: Mat3) =
  ## R^ac = ∇_c p^a = ∇_c projectTAH(X Y)^a = - ∇_c ∂_a tr[X Y + Y† X†], where M = X Y
  ## The derivatives ∂ is on X and ∇ is on Y.
  ## Note that for M = X Y†
  ## we would have - ∇_c ∂_a tr[X Y† + Y X†] = tr[T^c T^a X Y† + T^a T^c Y X†],
  ## or tr[T^a (M T^c + T^c M†)] = tr[T^c (T^a M + M† T^a)] = 2 ReTr[T^c T^a M]
  ## which would be the same as diffProjectTAH with a negative sign and swap M for M†,
  ## or transposed diffProjectTAH with a negative sign.
  ## Adx = SU3Ad(x)
  ## dp = diffProjectTAH(m, p)
  #[
    ∇_c P^a = - 2 ReTr[T^a X T^c Y]
            = - tr[T^a (X T^c X† X Y + Y† X† X T^c X†)]
            = - tr[T^a (T^b M + M† T^b)] AdX^bc
  ]#
  r := dp * Adx

const diffExpC = [
  1.0,
  1.0/2.0,
  1.0/6.0,
  1.0/24.0,
  1.0/120.0,
  1.0/720.0,
  1.0/5040.0,
  1.0/40320.0,
  1.0/362880.0,
  1.0/3628800.0,
  1.0/39916800.0,
  1.0/479001600.0,
  1.0/6227020800.0,
  1.0/87178291200.0]

proc oddHalfOrder(order: int): int {.inline.} =
  ## order = 2*result - 1.
  if order < 1 or (order and 1) == 0:
    raise newException(ValueError, "projected-exponential gradients require a positive odd order")
  (order + 1) div 2

#[
  P_n(A) = sum(k=0..n) A^k/(k+1)!
  D P_n(A)[E]
    = sum(k=1..n) sum(j=0..k-1) A^j E A^(k-1-j)/(k+1)!

  Pair powers in A^2; for SU(3), reduce them by Cayley-Hamilton.
]#

proc squareSym(r: var Mat1, x: Mat2) =
  ## r = x*x; compute one symmetric triangle.
  mixin mul, imadd
  const n = r.nrows
  when r.ncols != n or x.nrows != n or x.ncols != n:
    static: error("squareSym requires equally sized square matrices")
  forStatic i, 0, n-1:
    forStatic j, i, n-1:
      var t {.noinit.}: evalType(r[i, j])
      mul(t, x[i, 0], x[0, j])
      forStatic k, 1, n-1:
        imadd(t, x[i, k], x[k, j])
      r[i, j] := t
  forStatic i, 1, n-1:
    forStatic j, 0, i-1:
      r[i, j] := r[j, i]

proc mulSymSkew(r: var Mat1, s, x: Mat2) =
  ## r = s*x for commuting symmetric s and skew x.
  mixin mul, imadd
  const n = r.nrows
  when r.ncols != n or s.nrows != n or s.ncols != n or x.nrows != n or x.ncols != n:
    static: error("mulSymSkew requires equally sized square matrices")
  forStatic i, 0, n-1:
    r[i, i] := 0
  forStatic i, 0, n-2:
    forStatic j, i+1, n-1:
      var t {.noinit.}: evalType(r[i, j])
      mul(t, s[i, 0], x[0, j])
      forStatic k, 1, n-1:
        imadd(t, s[i, k], x[k, j])
      r[i, j] := t
      r[j, i] := -t

proc mulCommSym(r: var Mat1, a, b: Mat2) =
  ## r = a*b for commuting symmetric a and b.
  mixin mul, imadd
  const n = r.nrows
  when r.ncols != n or a.nrows != n or a.ncols != n or b.nrows != n or b.ncols != n:
    static: error("mulCommSym requires equally sized square matrices")
  forStatic i, 0, n-1:
    forStatic j, i, n-1:
      var t {.noinit.}: evalType(r[i, j])
      mul(t, a[i, 0], b[0, j])
      forStatic k, 1, n-1:
        imadd(t, a[i, k], b[k, j])
      r[i, j] := t
  forStatic i, 1, n-1:
    forStatic j, 0, i-1:
      r[i, j] := r[j, i]

proc mulAddI(r: var Mat1, a, b: Mat2) =
  ## r = I + a*b.
  mixin mul, imadd
  const n = r.nrows
  when r.ncols != n or a.nrows != n or a.ncols != n or b.nrows != n or b.ncols != n:
    static: error("mulAddI requires equally sized square matrices")
  forStatic i, 0, n-1:
    forStatic j, 0, n-1:
      var t {.noinit.}: evalType(r[i, j])
      mul(t, a[i, 0], b[0, j])
      forStatic k, 1, n-1:
        imadd(t, a[i, k], b[k, j])
      when i == j:
        t += 1.0
      r[i, j] := t

# (x*a)[i,j], skipping structural zeros of the SU(3) adjoint x.
template su3AdMulAt(x, a: untyped, i, j: static int): untyped =
  when i == 0:
    x[0, 1]*a[1, j] + x[0, 2]*a[2, j] + x[0, 3]*a[3, j] + x[0, 4]*a[4, j] + x[0, 5]*a[5, j] + x[0, 6]*a[6, j]
  elif i == 1:
    x[1, 0]*a[0, j] + x[1, 2]*a[2, j] + x[1, 3]*a[3, j] + x[1, 4]*a[4, j] + x[1, 5]*a[5, j] + x[1, 6]*a[6, j]
  elif i == 2:
    x[2, 0]*a[0, j] + x[2, 1]*a[1, j] + x[2, 3]*a[3, j] + x[2, 4]*a[4, j] + x[2, 5]*a[5, j] + x[2, 6]*a[6, j]
  elif i == 3:
    x[3, 0]*a[0, j] + x[3, 1]*a[1, j] + x[3, 2]*a[2, j] + x[3, 4]*a[4, j] + x[3, 5]*a[5, j] + x[3, 6]*a[6, j] + x[3, 7]*a[7, j]
  elif i == 4:
    x[4, 0]*a[0, j] + x[4, 1]*a[1, j] + x[4, 2]*a[2, j] + x[4, 3]*a[3, j] + x[4, 5]*a[5, j] + x[4, 6]*a[6, j] + x[4, 7]*a[7, j]
  elif i == 5:
    x[5, 0]*a[0, j] + x[5, 1]*a[1, j] + x[5, 2]*a[2, j] + x[5, 3]*a[3, j] + x[5, 4]*a[4, j] + x[5, 6]*a[6, j] + x[5, 7]*a[7, j]
  elif i == 6:
    x[6, 0]*a[0, j] + x[6, 1]*a[1, j] + x[6, 2]*a[2, j] + x[6, 3]*a[3, j] + x[6, 4]*a[4, j] + x[6, 5]*a[5, j] + x[6, 7]*a[7, j]
  else:
    x[7, 3]*a[3, j] + x[7, 4]*a[4, j] + x[7, 5]*a[5, j] + x[7, 6]*a[6, j]

proc diffExp13X4(r: var Mat1, x, x2, x3, x4: Mat2) =
  ## P_13(x), blocked in x^4.
  var f {.noinit.}: evalType(x)
  f := diffExpC[12]
  f += diffExpC[13]*x
  r := x4*f
  r += diffExpC[8]
  r += diffExpC[9]*x
  r += diffExpC[10]*x2
  r += diffExpC[11]*x3
  f := x4*r
  f += diffExpC[4]
  f += diffExpC[5]*x
  f += diffExpC[6]*x2
  f += diffExpC[7]*x3
  r := x4*f
  r += diffExpC[0]
  r += diffExpC[1]*x
  r += diffExpC[2]*x2
  r += diffExpC[3]*x3

proc reduceSu3AdPoly[T](c: var array[7, T], s1, s3: T) =
  ## z^4 = -s1*z^3 - s1^2*z^2/4 - s3*z, z = suad(X)^2.
  let s2 = 0.25*s1*s1
  for k in countdown(6, 4):
    c[k-1] -= s1*c[k]
    c[k-2] -= s2*c[k]
    c[k-3] -= s3*c[k]

proc diffExp13Su3Ad(r: var Mat1, x, x2, x4, x6: Mat2) =
  ## P_13(x), reduced by the SU(3) adjoint identity above.
  type T = evalType(x[0, 0])
  var
    even {.noinit.}: array[7, T]
    odd {.noinit.}: array[7, T]
  for i in 0..<7:
    even[i] := diffExpC[2*i]
    odd[i] := diffExpC[2*i + 1]
  let
    s1 = -0.5*trace(x2)
    s3 = -(s1*s1*s1)/12.0 - trace(x6)/6.0
  even.reduceSu3AdPoly(s1, s3)
  odd.reduceSu3AdPoly(s1, s3)
  var
    o {.noinit.}: evalType(x)
    xo {.noinit.}: evalType(x)
  forStatic i, 0, 7:
    forStatic j, i, 7:
      var
        re {.noinit.}: T
        ro {.noinit.}: T
      when i == j:
        re := even[0]
        ro := odd[0]
      else:
        re := 0
        ro := 0
      re += even[1]*x2[i, j]
      re += even[2]*x4[i, j]
      re += even[3]*x6[i, j]
      ro += odd[1]*x2[i, j]
      ro += odd[2]*x4[i, j]
      ro += odd[3]*x6[i, j]
      r[i, j] := re
      o[i, j] := ro
  forStatic i, 1, 7:
    forStatic j, 0, i-1:
      r[i, j] := r[j, i]
      o[i, j] := o[j, i]
  xo.mulSymSkew(o, x)
  r += xo

proc diffExpX2(r: var Mat1, adX, x2: Mat2, order: int) =
  ## P_order(adX), given x2 = adX^2.
  if order <= 0:
    r := 1.0
    return
  if order == 1:
    r := 1.0 + 0.5*adX
    return
  var n = order
  if (n and 1) != 0:
    dec n
  var c = 1.0
  for k in 2..(n+1):
    c /= float(k)
  r := c
  if n < order:
    r += (c/float(n+2))*adX
  while n > 0:
    r := x2*r
    n -= 2
    c *= float(n+2)*float(n+3)
    r += c
    r += (c/float(n+2))*adX

proc diffExp*(r: var Mat1, adX: Mat2, order=13) =
  ## r = P_order(adX) = sum(k=0..order) adX^k/(k+1)!.
  if order <= 0:
    r := 1.0
    return
  if order == 1:
    r := 1.0 + 0.5*adX
    return
  let x2 = adX*adX
  if order == 13:
    let
      x3 = x2*adX
      x4 = x2*x2
    r.diffExp13X4(adX, x2, x3, x4)
  else:
    r.diffExpX2(adX, x2, order)

proc diffExpApply*(r: var Vec1, adX: Mat1, x: Vec2, order=13) {.inline.} =
  ## r = P_order(adX)*x.
  matchMatVec("diffExpApply", adX, x)
  when r.len != x.len:
    static: error("diffExpApply result and input lengths differ")
  if order <= 0:
    r := x
    return
  var n = order
  var c = 1.0
  for k in 2..(n+1):
    c /= float(k)
  r := c*x
  for k in countdown(n, 1):
    c *= float(k+1)
    r := adX*r + c*x

proc diffExpSuApply*(r: var Vec1, a: Mat1, x: Vec2, order=13) {.inline.} =
  ## r = P_order(suad(a))*x, using nonzero SU(3) f constants.
  matchGroupVec("diffExpSuApply", a, r)
  matchGroupVec("diffExpSuApply", a, x)
  when r.len != x.len:
    static: error("diffExpSuApply result and input lengths differ")
  if order <= 0:
    r := x
    return
  when r.len == 8:
    var v {.noinit.}: evalType(r)
    v.suToVec(a)
    v := -v
    const f = 0.86602540378443864676 # sqrt(3/4)
    let
      v0 = v[0]
      v1 = v[1]
      v2 = v[2]
      v3 = v[3]
      v4 = v[4]
      v5 = v[5]
      v6 = v[6]
      v7 = v[7]
      h0 = 0.5*v0
      h1 = 0.5*v1
      h2 = 0.5*v2
      h3 = 0.5*v3
      h4 = 0.5*v4
      h5 = 0.5*v5
      h6 = 0.5*v6
      f3 = f*v3
      f4 = f*v4
      f5 = f*v5
      f6 = f*v6
      f7 = f*v7
    template adApply(y, z: untyped) =
      block:
        let
          z0 = z[0]
          z1 = z[1]
          z2 = z[2]
          z3 = z[3]
          z4 = z[4]
          z5 = z[5]
          z6 = z[6]
          z7 = z[7]
        y[0] := v2*z1 - v1*z2 + h6*z3 - h5*z4 + h4*z5 - h3*z6
        y[1] := -v2*z0 + v0*z2 + h5*z3 + h6*z4 - h3*z5 - h4*z6
        y[2] := v1*z0 - v0*z1 + h4*z3 - h3*z4 - h6*z5 + h5*z6
        y[3] := -h6*z0 - h5*z1 - h4*z2 + h2*z4 + h1*z5 + h0*z6 + f7*z4 - f4*z7
        y[4] := h5*z0 - h6*z1 + h3*z2 - h2*z3 - h0*z5 + h1*z6 - f7*z3 + f3*z7
        y[5] := -h4*z0 + h3*z1 + h6*z2 - h1*z3 + h0*z4 - h2*z6 + f7*z6 - f6*z7
        y[6] := h3*z0 + h4*z1 - h5*z2 - h0*z3 - h1*z4 + h2*z5 - f7*z5 + f5*z7
        y[7] := f4*z3 - f3*z4 + f6*z5 - f5*z6
    var t {.noinit.}: evalType(r)
    if order == 13:
      type T = evalType(r[0])
      var even, odd: array[7, T]
      for i in 0..<7:
        even[i] := diffExpC[2*i]
        odd[i] := diffExpC[2*i + 1]
      let
        n2 = a.norm2
        di = determinant(a).im
        s1 = 3.0*n2
        s3 = 0.5*n2*n2*n2 - 27.0*di*di
      even.reduceSu3AdPoly(s1, s3)
      odd.reduceSu3AdPoly(s1, s3)
      r := odd[3]*x
      for k in countdown(6, 0):
        adApply(t, r)
        case k
        of 0: r := t + even[0]*x
        of 1: r := t + odd[0]*x
        of 2: r := t + even[1]*x
        of 3: r := t + odd[1]*x
        of 4: r := t + even[2]*x
        of 5: r := t + odd[2]*x
        else: r := t + even[3]*x
    else:
      var n = order
      var c = 1.0
      for k in 2..(n+1):
        c /= float(k)
      r := c*x
      for k in countdown(n, 1):
        c *= float(k+1)
        adApply(t, r)
        r := t + c*x
  else:
    static: error("diffExpSuApply is implemented only for SU(3)")

proc expProjectTAHPullback*(r: var Mat1, m: Mat2, x: Mat3, order=13) {.inline.} =
  ## A = projectTAH(m), E = exp(A), x = E† C.
  ## redot(r, dm) = redot(C, dE).
  when r.nrows != m.nrows or r.nrows != x.nrows:
    {.error: "expProjectTAHPullback requires matrices of the same size".}
  when r.nrows == 1:
    r.projectTAH(x)
  elif r.nrows == 3:
    var a, q {.noinit.}: evalType(r)
    a.projectTAH(m)
    q.projectTAH(x)
    var v, p {.noinit.}: evalType(suToVec(a))
    v.suToVec(q)
    p.diffExpSuApply(a, v, order)
    r.suFromVec(p)
  else:
    {.error: "expProjectTAHPullback supports only 1x1 and 3x3 matrices".}

#[
  F = projectTAH(X (Y V W†)†)
    = projectTAH(X W V† Y†)
    = projectTAH(M)
    = - T^b tr[T^b (M - M†)]
  M = X W V† Y†
  Z = exp(F) X
  exp(adF) J(F) = J(-F)
  ∂_X^c Z = exp(F) T^c X + exp(F) J(F)[∂_X^c F] X
          = exp(F) (T^c + J(F)[∂_X^c F]) exp(-F) Z
          = exp(adF)[T^c + J(F)[∂_X^c F]] Z
          = exp(adF)[T^c + T^e J(F)^eb [∂_X^c F]^b] Z
          = T^a { exp(adF)^ac + J(-F)^ab [∂_X^c F]^b } Z
  ∂_Y^c Z = exp(F) J(F)[∂_Y^c F] X
          = exp(F) J(F)[∂_Y^c F] exp(-F) Z
          = exp(adF)[J(F)[∂_Y^c F]] Z
          = exp(adF)[T^e J(F)^eb [∂_Y^c F]^b] Z
          = T^a J(-F)^ab [∂_Y^c F]^b Z
  ∂_V^c Z = T^a J(-F)^ab [∂_V^c F]^b Z
  ∂_W^c Z = T^a J(-F)^ab [∂_W^c F]^b Z
  ∂_X^c F = - T^b tr[T^b (T^c M + M† T^c)]
  ∂_Y^c F =   T^b tr[T^b (T^c M† + M T^c)]
  ∂_V^c F =   T^b tr[T^b (X W V† T^c Y† + Y T^c V W† X†)]
          =   T^b tr[T^b (T^d M† + M T^d)] AdY^dc
  ∂_W^c F = - T^b tr[T^b (Y V W† T^c X† + X T^c W V† Y†)]
          = - T^b tr[T^b (T^d M + M† T^d)] AdX^dc
]#

proc diffExpProjectTAHMul*(J,JF,dF,adF: var Mat1, F:var Mat2, M: Mat3, order=13) =
  ## return F = projectTAH(X Y†) = - T^b ∂_b tr[X Y† + Y X†], and
  ## the simplified J = δ^ac + J(F)^ab [∂_c F]^b
  ## Note it only has the same determinant as the actual J = exp(adF)^ac + J(-F)^ab [∂_c F]^b.
  ## with M = X Y†
  ## assuming X and Y are independent.
  ## Only works with positive det(∂Z/∂X).
  ## Derivative target: base Jacobian is with respect to the updating link (∂ on X inside M).
  #[
    Z = exp(- T^b ∂_b tr[X Y† + Y X†]) X,  for X,Y in G, and ∂_b X = T_b X
    M = X Y†,  M in G
    Z = exp(F) X,  F in g
    F = - T^b ∂_b tr[M + M†]
      = - T^b tr[T^b M - M† T^b]
      = - T^b tr[T^b (M - M†)]
    ∂_c Z = exp(F) T^c X + exp(F) J(F)[∂_c F] X
          = exp(F) { T^c + J(F)[∂_c F] } exp(-F) Z
          = exp(adF)[T^c + J(F)[∂_c F]] Z
          = exp(adF)[T^c + T^e J(F)^eb [∂_c F]^b] Z
          = T^a { exp(adF)^ac + [exp(adF)J(F)]^ab [∂_c F]^b } Z
          = T^a { exp(adF)^ac + [(exp(adF)-1)/adF]^ab [∂_c F]^b } Z
          = T^a { exp(adF)^ac + J(-F)^ab [∂_c F]^b } Z
    [∂_c F]^b = - ∂_c ∂_b tr[M + M†]
              = - tr[T^b T^c M + M† T^c T^b]
              = - 1/2 tr[{T^b,T^c} (M + M†)] - 1/2 tr[[T^b,T^c] (M - M†)]
              = - 1/2 { d^bcd tr[T^d i (M + M†)] - 1/N δ^bc tr[M + M†] + f^bcd tr[T^d (M - M†)] }
              = - 1/2 { d^bcd tr[T^d i (M + M†)] - 1/N δ^bc tr[M + M†] - f^bcd F^d }
              = - 1/2 { d^bcd tr[T^d i (M + M†)] - 1/N δ^bc tr[M + M†] + adF^bc }
    -2 tr[(∂_c Z) Z† T^a]
        = exp(adF)^ac + J(-F)^ab [∂_c F]^b
        = exp(adF)^ac
          - 1/2 [(exp(adF)-1)/adF]^ab {d^bcd tr[T^d i (M + M†)] - 1/N δ^bc tr[M + M†]}
          - 1/2 [(exp(adF)-1)]^ac
        = 1/2 { [(exp(adF)+1)]^ac - [(exp(adF)-1)/adF]^ab {d^bcd tr[T^d i (M + M†)] - 1/N δ^bc tr[M + M†]} }
        = 1/2 { [(exp(adF)+1)]^ac - J(-F)^ab {d^bcd tr[T^d i (M + M†)] - 1/N δ^bc tr[M + M†]} }
    det(-2 tr[(∂_c Z) Z† T^a])
        = det(exp(adF)^ac + J(-F)^ab [∂_c F]^b)
        = det(δ^ac + J(F)^ab [∂_c F]^b)
  ]#
  F.projectTAH(M)
  dF.diffProjectTAH(M,F)
  adF.suad F
  if order == 13:
    var x, x2, x4, x6 {.noinit.}: evalType(adF)
    x := -adF
    x2.squareSym(x)
    x4.squareSym(x2)
    x6.mulCommSym(x4, x2)
    JF.diffExp13Su3Ad(x, x2, x4, x6)
  else:
    JF.diffExp(-adF, order=order)
  J.mulAddI(JF, dF)

proc buildJAndInvFromM(invJ, J, JF, dF, adF: var Mat1, F: var Mat2, M: Mat3, order=13) =
  ## J = I + P(-adF)*dF; also form J^-1.
  J.diffExpProjectTAHMul(JF, dF, adF, F, M, order=order)
  invJ.inverse J

proc buildDirectionalTerms(d2F: var Mat1, dFd: var Vec1, TM: Mat2) =
  ## d2F = D projectTAH(TM), dFd = vec(projectTAH(TM)).
  var pTM: evalType(TM)
  pTM.projectTAH TM
  d2F.diffProjectTAH(TM, pTM)
  let vVec = suToVec(pTM)
  dFd := vVec

proc accumulateGrad(invJ, JF, dFbase, adF: Mat1, d2F: Mat2, dFd: Vec1, halfOrder: int): auto {.noinit.} =
  ## tr(invJ * (D P(-adF)[dadF]*dFbase + JF*d2F)).
  var dadF, dJF: evalType(invJ)
  dadF.sufabc(dFd)
  dJF.diffDiffExp(-adF, dadF, halfOrder=halfOrder)
  result = trace(invJ * (dJF * dFbase + JF * d2F))

proc diffDiffExp11X4(r: var Mat1, x, dx, x2, x3, x4: Mat2) =
  ## D P_11(x)[dx], blocked in x^4.
  var dx2, dx3, dx4, f, t {.noinit.}: evalType(x)
  dx2 := x*dx
  dx2 += dx*x
  dx3 := dx2*x
  dx3 += x2*dx
  dx4 := dx2*x2
  dx4 += x2*dx2

  f := diffExpC[8]
  f += diffExpC[9]*x
  f += diffExpC[10]*x2
  f += diffExpC[11]*x3
  r := diffExpC[9]*dx
  r += diffExpC[10]*dx2
  r += diffExpC[11]*dx3

  t := dx4*f
  t += x4*r
  t += diffExpC[5]*dx
  t += diffExpC[6]*dx2
  t += diffExpC[7]*dx3
  r := t
  t := x4*f
  t += diffExpC[4]
  t += diffExpC[5]*x
  t += diffExpC[6]*x2
  t += diffExpC[7]*x3
  f := t

  t := dx4*f
  t += x4*r
  t += diffExpC[1]*dx
  t += diffExpC[2]*dx2
  t += diffExpC[3]*dx3
  r := t

proc traceProductSym(a, s: Mat1): auto {.noinit.} =
  ## tr(a*s), using one triangle of symmetric s.
  const n = a.nrows
  when a.ncols != n or s.nrows != n or s.ncols != n:
    static: error("traceProductSym requires equally sized square matrices")
  var r {.noinit.}: evalType(a[0, 0])
  r := a[0, 0]*s[0, 0]
  forStatic i, 1, n-1:
    r += a[i, i]*s[i, i]
  forStatic i, 0, n-2:
    forStatic j, i+1, n-1:
      r += (a[i, j] + a[j, i])*s[i, j]
  r

proc traceProductSkew(a, x: Mat1): auto {.noinit.} =
  ## tr(a*x), using one triangle of skew x.
  const n = a.nrows
  when a.ncols != n or x.nrows != n or x.ncols != n:
    static: error("traceProductSkew requires equally sized square matrices")
  var r {.noinit.}: evalType(a[0, 0])
  r := (a[1, 0] - a[0, 1])*x[0, 1]
  forStatic i, 0, n-2:
    forStatic j, i+1, n-1:
      when i != 0 or j != 1:
        r += (a[j, i] - a[i, j])*x[i, j]
  r

proc diffDiffExp13Su3AdRevF(r: var Vec1, x, b, x2, x3, x4: Mat2, x3n2: auto) =
  ## r[g] = tr(b * D P_13(x)[sufabc(e_g)]).
  #[
    z = x^2,  P_13(x) = E(z) + x O(z)
    z^4 = -s1*z^3 - s1^2*z^2/4 - s3*z
    s1 = -tr(x^2)/2,  s3 = -s1^3/12 - tr(x^6)/6
    Differentiate the reduction; contract the 25 nonzero f-pairs.
  ]#
  type T = evalType(x[0, 0])
  var
    even {.noinit.}: array[7, T]
    odd {.noinit.}: array[7, T]
    evenA {.noinit.}: array[7, T]
    oddA {.noinit.}: array[7, T]
    evenD {.noinit.}: array[7, T]
    oddD {.noinit.}: array[7, T]
  let
    s1 = -0.5*trace(x2)
    s3 = (x3n2 - 0.5*s1*s1*s1)/6.0
    s2 = 0.25*s1*s1
  template reduce13(c, ca, cd: untyped, p: static int) =
    # Reduce z^6..z^4 and their derivatives in s1 and s3.
    for i in 0..<7:
      c[i] := diffExpC[p+2*i]
      ca[i] := 0
      cd[i] := 0
    for k in countdown(6, 4):
      let
        q = c[k]
        qa = ca[k]
        qd = cd[k]
      c[k-1] -= s1*q
      ca[k-1] -= q + s1*qa
      cd[k-1] -= s1*qd
      c[k-2] -= s2*q
      ca[k-2] -= 0.5*s1*q + s2*qa
      cd[k-2] -= s2*qd
      c[k-3] -= s3*q
      ca[k-3] -= s3*qa
      cd[k-3] -= q + s3*qd
  reduce13(even, evenA, evenD, 0)
  reduce13(odd, oddA, oddD, 1)

  var
    hi {.noinit.}: evalType(x)
    u4 {.noinit.}: evalType(x)
    uhi {.noinit.}: evalType(x)
    u2 {.noinit.}: evalType(x)
    u3 {.noinit.}: evalType(x)
  forStatic i, 0, 7:
    forStatic j, 0, 7:
      when i == j:
        hi[i, j] := even[2]
      else:
        hi[i, j] := 0
      hi[i, j] += odd[2]*x[i, j]
      hi[i, j] += even[3]*x2[i, j]
      hi[i, j] += odd[3]*x3[i, j]
  mul(u4, hi, b)
  mul(uhi, b, x4)
  forStatic i, 0, 7:
    forStatic j, 0, 7:
      u2[i, j] := even[1]*b[i, j]
      u2[i, j] += even[3]*uhi[i, j]
      u3[i, j] := odd[1]*b[i, j]
      u3[i, j] += odd[3]*uhi[i, j]

  let
    ge1 = traceProductSym(b, x2)
    go1 = traceProductSkew(b, x3)
    ge2 = trace(uhi)
    go2 = traceProductSkew(uhi, x)
    ge3 = traceProductSym(uhi, x2)
    go3 = traceProductSkew(uhi, x3)
  var ga = ge1*evenA[1] + go1*oddA[1]
  var gd = ge1*evenD[1] + go1*oddD[1]
  ga += ge2*evenA[2] + go2*oddA[2]
  gd += ge2*evenD[2] + go2*oddD[2]
  ga += ge3*evenA[3] + go3*oddA[3]
  gd += ge3*evenD[3] + go3*oddD[3]
  u3 -= (gd/3.0)*x3
  ga -= 0.25*gd*s1*s1
  forStatic i, 0, 7:
    u2[i, i] -= 0.5*ga

  imadd(u2, x2, u4)
  imadd(u2, u4, x2)
  # Only the symmetric part of u2 + x*u3 contributes below.
  forStatic i, 0, 7:
    forStatic j, i, 7:
      var up {.noinit.}: T
      up := u2[i, j] + su3AdMulAt(x, u3, i, j)
      when i == j:
        u4[i, j] := up + up
      else:
        var lo {.noinit.}: T
        lo := u2[j, i] + su3AdMulAt(x, u3, j, i)
        u4[i, j] := up + lo
        u4[j, i] := u4[i, j]
  template mulAsymAt(a, b: untyped, i, j: static int): untyped =
    block:
      var lo, up {.noinit.}: T
      lo := a[j, 0]*b[0, i]
      up := a[i, 0]*b[0, j]
      forStatic l, 1, 7:
        lo += a[j, l]*b[l, i]
        up += a[i, l]*b[l, j]
      lo - up
  template k(i, j: static int): untyped =
    odd[0]*(b[j, i] - b[i, j]) + odd[2]*(uhi[j, i] - uhi[i, j]) + mulAsymAt(u3, x2, i, j) + su3AdMulAt(x, u4, j, i) - su3AdMulAt(x, u4, i, j)
  const
    h = 0.5
    r3h = 0.86602540378443864676
  let
    k34 = k(3, 4)
    k56 = k(5, 6)
  r[0] := k(1, 2) + h*(k(3, 6) - k(4, 5))
  r[1] := -k(0, 2) + h*(k(3, 5) + k(4, 6))
  r[2] := k(0, 1) + h*(k34 - k56)
  r[3] := -h*(k(0, 6) + k(1, 5) + k(2, 4)) + r3h*k(4, 7)
  r[4] := h*(k(0, 5) - k(1, 6) + k(2, 3)) - r3h*k(3, 7)
  r[5] := h*(-k(0, 4) + k(1, 3) + k(2, 6)) + r3h*k(6, 7)
  r[6] := h*(k(0, 3) + k(1, 4) - k(2, 5)) - r3h*k(5, 7)
  r[7] := r3h*(k34 + k56)

proc diffDiffExpX2(r: var Mat1, adX, dadX, x2: Mat2, halfOrder: int) =
  ## r = D P_(2*halfOrder-1)(adX)[dadX], given x2 = adX^2.
  if halfOrder <= 1:
    r := 0.5*dadX
    return
  let dx2 = adX*dadX + dadX*adX
  var
    f: evalType(adX)
    n = 2*halfOrder-2
    c = 1.0
  for k in 2..(n+1):
    c /= float(k)
  var co = c/float(n+2)
  f := c
  f += co*adX
  r := co*dadX
  while n > 0:
    n -= 2
    c *= float(n+2)*float(n+3)
    co = c/float(n+2)
    r := dx2*f + x2*r + co*dadX
    if n > 0:
      f := x2*f
      f += c
      f += co*adX

proc diffDiffExp*(r: var Mat1, adX: Mat2, dadX: Mat3, halfOrder=6) =
  ## r = sum(k=1..n) sum(j=0..k-1) adX^j*dadX*adX^(k-1-j)/(k+1)!, n = 2*halfOrder-1.
  if halfOrder <= 1:
    r := 0.5*dadX
    return
  let x2 = adX*adX
  if halfOrder == 6:
    let
      x3 = x2*adX
      x4 = x2*x2
    r.diffDiffExp11X4(adX, dadX, x2, x3, x4)
  else:
    r.diffDiffExpX2(adX, dadX, x2, halfOrder)

proc expProjMulLogJac*[T](M: MatrixArray[1, 1, T], order=13): auto {.inline.} =
  ## ln J = ln(1 + Re M).
  discard order
  ln(1.0 + M[0, 0].re)

proc expProjMulLogJac*(M: Mat1, order=13): auto {.noinit.} =
  ## F = projectTAH(M), D = diffProjectTAH(M,F), J = I + P(-suad(F))*D.
  ## Return ln det J.
  ## Requires det J > 0 and nonzero LU pivots.
  when M.nrows != 3:
    {.error: "expProjMulLogJac requires a 3x3 matrix".}
  type T = evalType(M[0, 0].re)
  var J, JF, dF, adF {.noinit.}: MatrixArray[8, 8, T]
  var F {.noinit.}: evalType(M)
  J.diffExpProjectTAHMul(JF, dF, adF, F, M, order=order)
  ln(detNoPivot(J))

proc addSu3FAsym(qf: var Vec1, A: Mat1, c: static float) {.inline.} =
  ## qf[g] += c*f[g,a,b]*(A[b,a] - A[a,b]), a < b.
  template sk(i, j: static int): untyped =
    c*(A[j, i] - A[i, j])
  template addF(a, b, c: static int, f: untyped) =
    qf[a] += f*sk(b, c)
    qf[b] -= f*sk(a, c)
    qf[c] += f*sk(a, b)
  const
    h = 0.5
    r3h = 0.86602540378443864676
  addF(0, 1, 2, 1.0)
  addF(0, 3, 6, h)
  addF(0, 4, 5, -h)
  addF(1, 3, 5, h)
  addF(1, 4, 6, h)
  addF(2, 3, 4, h)
  addF(2, 5, 6, -h)
  addF(3, 4, 7, r3h)
  addF(5, 6, 7, r3h)

proc projJacPullbackSu3(G: var Mat1, qfA: Vec1, C: Mat2) =
  #[
    T^a T^b = ((f^abc + i*d^abc) T^c - delta^ab/3)/2
    Expand f,d,delta and form G with dL = redot(G,dM).
  ]#
  type T = evalType(C[0, 0])
  var qf, qd {.noinit.}: VectorArray[8, T]
  qf := qfA
  qf.addSu3FAsym(C, 0.5)
  qd := 0

  template addD(a, b, c: static int, d: untyped) =
    qd[a] += d*(C[b, c] + C[c, b])
    qd[b] += d*(C[a, c] + C[c, a])
    qd[c] += d*(C[a, b] + C[b, a])
  template addD2(a, c: static int, d: untyped) =
    qd[a] += d*(C[a, c] + C[c, a])
    qd[c] += d*C[a, a]
  const
    h = 0.5
    r13 = sqrt_1_3
    r16 = 0.5*sqrt_1_3
  addD(0, 3, 5, -h)
  addD(0, 4, 6, -h)
  addD(1, 3, 6, h)
  addD(1, 4, 5, -h)
  addD2(0, 7, -r13)
  addD2(1, 7, -r13)
  addD2(2, 7, -r13)
  addD2(3, 2, -h)
  addD2(4, 2, -h)
  addD2(5, 2, h)
  addD2(6, 2, h)
  addD2(3, 7, r16)
  addD2(4, 7, r16)
  addD2(5, 7, r16)
  addD2(6, 7, r16)
  qd[7] += r13*C[7, 7]

  var F, D {.noinit.}: evalType(G)
  F.suFromVec(qf)
  D.suFromVec(qd)
  for i in 0..<3:
    for j in 0..<3:
      G[i, j].re := 2.0*F[i, j].re + D[i, j].im
      G[i, j].im := 2.0*F[i, j].im - D[i, j].re
  var trc = C[0, 0]
  for i in 1..<8:
    trc += C[i, i]
  for i in 0..<3:
    G[i, i].re += trc/3.0

proc projJacPullbackSu3(G: var Mat1, A, C: Mat2) =
  type T = evalType(A[0, 0])
  var qf {.noinit.}: VectorArray[8, T]
  qf := 0
  qf.addSu3FAsym(A, 1.0)
  G.projJacPullbackSu3(qf, C)

#[
  F = projectTAH(M),  x = -suad(F)
  P = diffExp(x),  D = diffProjectTAH(M,F)
  J = I + P*D
  d ln det J
    = tr(J^-1 * (dP*D + P*dD))
    = tr((D*J^-1)*dP) + tr((J^-1*P)*dD)
  Pull back both terms: d ln det J = redot(G,dM).
]#
proc expProjMulLogJacGradSu3(G: var Mat1, M: Mat2, p: var Vec1, v: Vec2, apply: static bool, order: int) =
  when M.nrows != 3 or G.nrows != 3:
    {.error: "expProjMulLogJacGradSu3 requires 3x3 matrices".}
  type T = evalType(M[0, 0].re)
  let halfOrder = oddHalfOrder(order)
  var
    J {.noinit.}: MatrixArray[8, 8, T]
    JF {.noinit.}: MatrixArray[8, 8, T]
    dF {.noinit.}: MatrixArray[8, 8, T]
    x {.noinit.}: MatrixArray[8, 8, T]
    x2 {.noinit.}: MatrixArray[8, 8, T]
    x3 {.noinit.}: MatrixArray[8, 8, T]
    x4 {.noinit.}: MatrixArray[8, 8, T]
    x6 {.noinit.}: MatrixArray[8, 8, T]
  var F {.noinit.}: evalType(M)
  F.projectTAH(M)
  dF.diffProjectTAH(M, F)
  x.suad(F)
  x := -x
  x2.squareSym(x)
  if order == 13:
    x4.squareSym(x2)
    x6.mulCommSym(x4, x2)
    JF.diffExp13Su3Ad(x, x2, x4, x6)
  else:
    JF.diffExpX2(x, x2, order)
  when apply:
    for i in 0..<8:
      p[i] := JF[0, i] * v[0]
      for j in 1..<8:
        p[i] += JF[j, i] * v[j]
  J.mulAddI(JF, dF)
  J.solveLRNoPivot(JF, dF)
  if order == 13:
    var qf {.noinit.}: VectorArray[8, T]
    x3.mulSymSkew(x2, x)
    # x3 is skew-symmetric, so norm2(x3) = -trace(x6).
    let x3n2 = -trace(x6)
    qf.diffDiffExp13Su3AdRevF(x, dF, x2, x3, x4, x3n2)
    G.projJacPullbackSu3(qf, JF)
  else:
    J.diffDiffExpX2(x, dF, x2, halfOrder)
    G.projJacPullbackSu3(J, JF)

proc expProjMulLogJacGrad*(G: var Mat1, M: Mat2, order=13) =
  ## d ln det J = redot(G,dM). Requires det J > 0 and nonzero LU pivots.
  type V = evalType(suToVec(M))
  var p, v {.noinit.}: V
  G.expProjMulLogJacGradSu3(M, p, v, false, order)

proc expProjMulLogJacGrad*[T](G: var MatrixArray[1, 1, T], M: MatrixArray[1, 1, T], order=13) {.inline.} =
  ## G[0,0].re = 1/(1 + Re M).
  discard oddHalfOrder(order)
  G := 0
  G[0, 0].re := 1.0 / (1.0 + M[0, 0].re)

proc expProjMulLogJacGrad*(G: var Mat1, p: var Vec1, M: Mat2, v: Vec2, order=13) =
  ## Also set p = P(suad(projectTAH(M)))*v.
  matchGroupVec("expProjMulLogJacGrad", M, p)
  matchGroupVec("expProjMulLogJacGrad", M, v)
  when p.len != v.len:
    static: error("expProjMulLogJacGrad vector lengths differ")
  G.expProjMulLogJacGradSu3(M, p, v, true, order)

proc expProjMulLogJacGrad*(G: var Mat1, P: var Mat2, M: Mat3, X: Mat4, order=13) {.inline.} =
  ## A = projectTAH(M), E = exp(A), X = E† C.
  ## Also set P so redot(P, dM) = redot(C, dE).
  when G.nrows != P.nrows or G.nrows != M.nrows or G.nrows != X.nrows:
    {.error: "expProjMulLogJacGrad requires matrices of the same size".}
  when G.nrows == 1:
    G.expProjMulLogJacGrad(M, order)
    P.projectTAH(X)
  elif G.nrows == 3:
    var q {.noinit.}: evalType(G)
    q.projectTAH(X)
    var v, p {.noinit.}: evalType(suToVec(q))
    v.suToVec(q)
    G.expProjMulLogJacGrad(p, M, v, order)
    P.suFromVec(p)
  else:
    {.error: "expProjMulLogJacGrad supports only 1x1 and 3x3 matrices".}

proc diffLnDetDiffExpProjectTAHMul*(r: var Vec1, M: Mat1, order=13) =
  #[
    ∇_d ln det {δ^ac + J(F)^ab [∂_c F^b]}    # ∇ can act on different links, ∂ only on the updating link
        = m^{-1}^ca {[∇_d J(F)^ab] [∂_c F^b] + J(F)^ab [∇_d ∂_c F^b]}
    where
        m^ac = δ^ac + J(F)^ab [∂_c F^b]
    This is for ∇ = ∂.
    ∇_d ∂_c F^b = ∇_d tr[T^b T^c X Y + Y† X† T^c T^b]
                = tr[T^b T^c T^d X Y - Y† X† T^d T^c T^b]
                = tr[T^b (T^c T^d X Y - Y† X† T^d T^c)]  # identify M^d = T^d X Y
                = tr[T^b (T^c M^d + M^d† T^c)]  # diffProjectTAH*(r:var Mat1, m: Mat2, p: Mat3) =
                                                # r_{ac} = - tr[T^a (T^c M + M† T^c)]
    ∇_d adF^ce = ∇_d (- f^ceg F^g) = - f^ceg [∇_d F^g]
  ]#
  ## Derivative target: same-link gradient (∇ = ∂) acting on M.
  type T = evalType(M[0,0].re)
  const nc = M.nrows
  const dim = nc*nc-1
  const t = sugen(nc)
  type A = MatrixArray[dim, dim, T]
  type V = VectorArray[dim, T]
  let halfOrder = oddHalfOrder(order)
  var F,TM: evalType(M)
  var J,invJ,adF,JF,dF,d2F {.noinit.}: A
  var dFd {.noinit.}: V
  buildJAndInvFromM(invJ, J, JF, dF, adF, F, M, order=order)
  for d in 0..<dim:
    TM := t[d] * M
    buildDirectionalTerms(d2F, dFd, TM)
    # For same-link, ∇F equals the d-th column of dF (override dFd)
    for g in 0..<dim:
      dFd[g] = dF[g,d]
    r[d] = accumulateGrad(invJ, JF, dF, adF, d2F, dFd, halfOrder)

proc diffCrossGeneralLnDetDiffExpProjectTAHMul*(r: var Vec1, L: Mat1, Y: Mat2, R: Mat3, adjoint: static bool, order=13) =
  ## Unified cross-link gradient for ln det {I + J(F)[∂F]} with M = L · Y^σ · R.
  ## adjoint=false: σ=+1,  M = L · Y · R,  δM = L · (T^d Y) · R
  ## adjoint=true:  σ=−1,  M = L · Y† · R, δM = − (L · Y†) · T^d · R
  ## Derivative target: gradient with respect to the Y factor (or Y† when adjoint=true).
  type T = evalType(L[0,0].re)
  const nc = L.nrows
  const dim = nc*nc-1
  const t = sugen(nc)
  type A = MatrixArray[dim, dim, T]
  type V = VectorArray[dim, T]
  let halfOrder = oddHalfOrder(order)
  var M,F,TM: evalType(L)
  when adjoint:
    M := L * Y.adj * R
  else:
    M := L * Y * R
  var J,invJ,adF,JF,dF,d2F {.noinit.}: A
  var dFd {.noinit.}: V
  buildJAndInvFromM(invJ, J, JF, dF, adF, F, M, order=order)
  let dFbase = dF
  for d in 0..<dim:
    when adjoint:
      TM := - (L * Y.adj * t[d] * R)
    else:
      TM := L * t[d] * Y * R
    buildDirectionalTerms(d2F, dFd, TM)
    r[d] = accumulateGrad(invJ, JF, dFbase, adF, d2F, dFd, halfOrder)

proc diffCrossLnDetDiffExpProjectTAHMul*(r: var Vec1, X: Mat1, Y: Mat2, order=13) =
  #[
    ∇_d ln det {δ^ac + J(F)^ab [∂_c F^b]}    # ∇ can act on different links, ∂ only on the updating link
        = m^{-1}^ca {[∇_d J(F)^ab] [∂_c F^b] + J(F)^ab [∇_d ∂_c F^b]}
    where
        m^ac = δ^ac + J(F)^ab [∂_c F^b]
    This is for ∇ ≠ ∂, assuming ∇_d = ∇_Y^d and ∂_c = ∂_X^c, M = X Y
    ∇_d adF^ce = ∇_d (- f^ceg F^g) = - f^ceg [∇_d F^g]
    ∇_d F^g = ∇_d tr[T^g X Y - Y† X† T^g]
            = tr[T^g X T^d Y + Y† T^d X† T^g]
            = tr[T^g X T^d X† X Y + Y† X† X T^d X† T^g]
            = tr[T^g (T^f X Y + Y† X† T^f)] Ad(X)^fd
    ∇_d ∂_c F^b = ∇_d tr[T^b T^c X Y + Y† X† T^c T^b]
                = tr[T^b T^c X T^d Y - Y† T^d X† T^c T^b]
  ]#
  ## Derivative target: cross-link gradient with respect to Y (second factor in M = X · Y);
  ## base Jacobian is with respect to X (updating link).
  var I: evalType(X)
  I := 1.0
  diffCrossGeneralLnDetDiffExpProjectTAHMul(r, X, Y, I, adjoint=false, order=order)

# Cross-link for M = L · Y† · R, variation wrt Y (adjoint on middle factor)
proc diffCrossAdjLnDetDiffExpProjectTAHMul*(r: var Vec1, LYadj: Mat1, R: Mat2, order=13) =
  var I: evalType(LYadj)
  I := 1.0
  diffCrossGeneralLnDetDiffExpProjectTAHMul(r, LYadj, I, R, adjoint=true, order=order)

proc ndiffSUtoReal*(r: var Vec1, err: var Vec2, f: proc, x: Mat2, dx:float=2.0, scale:float=5.0, ordMax:static int=4) =
  ## for a function f: SU(N) → Real
  ## return the derivative in the vector space of su(n) algebra
  ## r_a = ∂_{l_a} f(exp(l_a T_a) x)
  ## using numerical differentiation, algarithms.numdiff.ndiff
  type T = evalType(x[0,0].re)
  const nc = x.nrows
  var z,d,dr,er:T
  z := 0.0
  d := dx
  const t = sugen(nc)
  for a in 0..<t.len:
    ndiff(dr, er,
      #proc (l:T):T {.noinit.} = f(exp(l*t[a])*x),
      proc (l:T):T = f(exp(l*t[a])*x),  # work around Nim cpp bug
      z, d, scale=scale, ordMax=ordMax)
    r[a] = dr
    err[a] = er

proc ndiffSUtoAlg*(r: var Mat1, err: var Mat2, f: proc, x: Mat3, dx:float=2.0, scale:float=5.0, ordMax:static int=4) =
  ## for a function f: SU(N) → su(n)
  ## return the jacobian in the vector space of su(n) algebra
  ## r_{ba} = ∂_{l_a} f_b(exp(l_a T_a) x)
  ## using numerical differentiation, algarithms.numdiff.ndiff
  type T = evalType(x[0,0].re)
  type V = evalType(suToVec(x))
  const nc = x.nrows
  var z,d:T
  var dr,er:V
  z := 0.0
  d := dx
  const t = sugen(nc)
  for a in 0..<t.len:
    ndiff(dr, er,
      proc (l:T):V {.noinit.} = result.suToVec(f(exp(l*t[a])*x)),
      z, d, scale=scale, ordMax=ordMax)
    for b in 0..<t.len:
      r[b,a] = dr[b]
      err[b,a] = er[b]

proc ndiffAlgtoSU*(r: var Mat1, err: var Mat2, f: proc, x: Mat3, dx:float=0.1, scale:float=5.0, ordMax:static int=4) =
  ## for a function f: su(n) → SU(N)
  ## return the jacobian in the vector space of su(n) algebra
  ## r_{ba} = - 2 ∂_{x_a} Tr[T^b f(x_a T_a) f(x)†]
  ## using numerical differentiation, algarithms.numdiff.ndiff
  type T = evalType(x[0,0].re)
  type V = evalType(suToVec(x))
  const nc = x.nrows
  var z,d:T
  var dr,er:V
  z := 0.0
  d := dx
  const t = sugen(nc)
  let fx = f(x)
  for a in 0..<t.len:
    ndiff(dr, er,
      proc (l:T):V {.noinit.} = result.suToVec(f(l*t[a]+x) * fx.adj),
      z, d, scale=scale, ordMax=ordMax)
    for b in 0..<t.len:
      r[b,a] = dr[b]
      err[b,a] = er[b]

proc ndiffSUtoSU*(r: var Mat1, err: var Mat2, f: proc, x: Mat3, dx:float=0.1, scale:float=5.0, ordMax:static int=4) =
  ## for a function f: SU(N) → SU(N)
  ## return the jacobian in the vector space of su(n) algebra
  ## r_{ba} = - 2 ∂_{l_a} Tr[T^b f(exp(l_a T_a) x) f(x)†]
  ## using numerical differentiation, algarithms.numdiff.ndiff
  type T = evalType(x[0,0].re)
  type V = evalType(suToVec(x))
  const nc = x.nrows
  var z,d:T
  var dr,er:V
  z := 0.0
  d := dx
  const t = sugen(nc)
  let fx = f(x)
  for a in 0..<t.len:
    ndiff(dr, er,
      proc (l:T):V {.noinit.} = result.suToVec(f(exp(l*t[a])*x) * fx.adj),
      z, d, scale=scale, ordMax=ordMax)
    for b in 0..<t.len:
      r[b,a] = dr[b]
      err[b,a] = er[b]

when isMainModule:
  import simd
  template check(s:string, x:untyped, n:SomeNumber):untyped =
    let r0 = x
    let r = simdSum(r0)/simdLength(r0)
    echo s, " error/eps: ", r/epsilon(r)
    doAssert(abs(r)<n*epsilon(r))
    # if abs(r)>=n*epsilon(r):
    #   echo "ERROR"

  proc testfabc =
    let fabctc = sufabc(su3gen)
    var del: evalType(su3gen[0])
    for a in 0..<8:
      for b in 0..<8:
        let f = su3gen[a] * su3gen[b] - su3gen[b] * su3gen[a]
        let d = fabctc[a,b] - f
        del += d
        # check("fabcTc[" & $a & "," & $b & "]", sqrt(norm2(d))/3, 1)
    check("fabc", sqrt(norm2(del))/64, 1)
  proc testdabc =
    const ii = newComplex(0, 1.0)
    let dabctc = sudabc(su3gen)
    var del: evalType(su3gen[0])
    for a in 0..<8:
      for b in 0..<8:
        let f = su3gen[a] * su3gen[b] + su3gen[b] * su3gen[a]
        var d = ii*dabctc[a,b] - f
        if a==b:
          d -= 1/3
        del += d
        # check("dabctc[" & $a & "," & $b & "]", sqrt(norm2(d))/3, 1)
    check("dabc", sqrt(norm2(del))/64, 1)

  testfabc()
  testdabc()

  proc test(T: typedesc) =
    const N = 3
    const D = N*N-1
    const v:array[D, float] = [-0.3, 0.2, 0.11, 0.23, -0.31, 0.03, -0.07, 0.17]
    type
      M = MatrixArray[N,N,ComplexType[T]]
      V = VectorArray[D,T]
      A = MatrixArray[D,D,T]
    var s1,s2,m0,m1,m2,m3: M
    var v0,v1,v2: V
    var a1,a2,a3: A
    const sl = simdLength(v1[0])
    for i in 0..<N:
      if sl==1:
        v0[i] := v[i]
      else:
        var a: array[sl, float]
        for k in 0..<sl:
          a[k] = v[i] + 0.3/(k.float-1.5)
        v0[i] := a
      for j in 0..<N:
        let fi = i.float
        let fj = j.float
        if sl==1:
          m1[i,j].re := 0.5 + 0.7/(0.9+1.3*fi-fj)
          m1[i,j].im := 0.1 + 0.3/(0.4+fi-1.1*fj)
        else:
          var a,b: array[sl,float]
          for k in 0..<sl:
            let fk = k.float
            a[k] = 0.5 + 0.7/(0.9+1.3*fi-fj+0.011*fk)
            b[k] = 0.1 + 0.3/(0.4+fi-1.1*fj-0.007*fk)
          m1[i,j].re := a
          m1[i,j].im := b
    m1.projectTAH
    s1 = exp(m1)
    echo "test " & $N & " " & $T
    # echo "m1: ", m1
    v1 = suToVec(m1)
    # echo "v1: ", v1
    v2 = suToVec_mat(m1)
    # echo "v2: ", v2
    check("suToVec", sqrt(norm2(v1-v2))/D, 1)
    m2 = suFromVec(v1)
    # echo "m2: ", m2
    check("suFromVec", sqrt(norm2(m1-m2))/N, 1)
    m2 = suFromVec_mat(v1)
    # echo "m2: ", m2
    check("suFromVec_mat", sqrt(norm2(m1-m2))/N, 1)

    a1.suad(v1, m1)
    a2.suad m1
    a3 := suad(m1)
    check("suad(a,v,m)", sqrt(norm2(a1-a3))/N, 1)
    check("suad(a,m)", sqrt(norm2(a2-a3))/N, 1)

    m0.suFromVec(v0)

    m1.suadApply(v1, a3, m0)
    m2.suadApply(a3, m0)
    m3 := suadApply(a3, m0)
    check("suadApply(m,v,a,m)", sqrt(norm2(m1-m2))/N, 1)
    check("suadApply(m,a,m)", sqrt(norm2(m1-m3))/N, 1)

    a1.SUAd(s1)
    a2 = SUAd_mat(s1)
    a3 = SUAd(s1)
    #[
    for b in 0..<D:
      for c in 0..<D:
        echo "a1[", b, ",", c, "]: ", a1[b,c]
        echo "a2[", b, ",", c, "]: ", a2[b,c]
        check("SUAd[" & $b & "," & $c & "]", sqrt(norm2(a1-a2))/N, 1)
    ]#
    check("SUAd(a,s)", sqrt(norm2(a1-a2))/D, 1)
    check("SUAd(s)", sqrt(norm2(a3-a2))/D, 1)

  template doTest(t:untyped) =
    when declared(t):
      test(t)
  doTest(float32)
  doTest(float64)
  doTest(SimdS1)
  doTest(SimdD1)
  doTest(SimdS2)
  doTest(SimdD2)
  doTest(SimdS4)
  doTest(SimdD4)
  doTest(SimdS8)
  doTest(SimdD8)
  doTest(SimdS16)
  doTest(SimdD16)
