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
    static: error(label & " v[" & $r.len & "] m[" & $m.nrows & "," & $m.ncols & "]")
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
  ## r_{ac} = ∂_c p^a = ∂_c projectTAH(m)^a = - tr[T^a (T^c M + M† T^c)]
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
  r.sudabc(v)
  t := (-0.5) * p
  # r += suad(t) + trMs    # FIXME: error in matrixOps.nim:/op MMS/+/op.* 0/
  r += trMs
  r += suad(t)

proc diffCrossProjectTAH*(r: var Mat1, Adx: Mat2, dp: Mat3) =
  ## R^ac = ∇_c p^a = ∇_c projectTAH(X Y)^a = - ∇_c ∂_a tr[X Y + Y† X†], where M = X Y
  ## The derivatives ∂ is on X and ∇ is on Y.
  ## Adx = SU3Ad(x)
  ## dp = diffprojectTAH(m, p)
  #[
    ∇_c P^a = - 2 ReTr[T^a X T^c Y]
            = - tr[T^a (X T^c X† X Y + Y† X† X T^c X†)]
            = - tr[T^a (T^b M + M† T^b)] AdX^bc
  ]#
  r := dp * Adx

proc diffExp*(r: var Mat1, adX: Mat2, order=13) =
  ## return J(X) = (1-exp{-adX})/adX = Σ_{k=0}^\infty 1/(k+1)! (-adX)^k  upto k=order
  #[
    [exp{-X(t)} d/dt exp{X(t)}]_ij = [J(X) d/dt X(t)]_ij = T^a_ij J(X)^ab (-2) T^b_kl [d/dt X(t)]_lk
    J(X) = 1 + 1/2 (-adX) (1 + 1/3 (-adX) (1 + 1/4 (-adX) (1 + ...)))
    J(x) ∂_t x
        = T^a J(x)^ab (-2) tr[T^b ∂_t x]
        = exp(-x) ∂_t exp(x)
    J(s x) ∂_t x = exp(-s x) ∂_t exp(s x)
    ∂_s J(s x) ∂_t x
        = - exp(-s x) x ∂_t exp(s x) + exp(-s x) ∂_t x exp(s x)
        = - exp(-s x) x ∂_t exp(s x) + exp(-s x) [∂_t x] exp(s x) + exp(-s x) x ∂_t exp(s x)
        = exp(-s x) [∂_t x] exp(s x)
        = exp(-s adx) ∂_t x
        = Σ_k 1/k! (-1)^k s^k (adx)^k ∂_t x
    J(0) = 0
    J(x) ∂_t x
        = ∫_0^1 ds Σ_{k=0} 1/k! (-1)^k s^k (adx)^k ∂_t x
        = Σ_{k=0} 1/(k+1)! (-1)^k (adx)^k ∂_t x
  ]#
  r := 1.0 + (-1.0)/(order+1.0) * adX
  for i in countdown(order, 2):
    r := 1.0 + (-1.0)/float(i) * (adX * r)

proc smearIndepLogDetJacobian*(F:var Mat1, X: Mat1, Y: Mat2, diffExpOrder=13): auto =
  ## return T^b ∂_b tr[X Y† + Y X†], and log det(∂Z/∂X)
  ## assuming X and Y are independent.
  ## Z = exp(T^b ∂_b tr[X Y† + Y X†]) X,  for X in G, and ∂_b X = T_b X
  let M = X * Y.adj
  F.projectTAH(M)
  var dadF = suad(F)
  dadF.diffExp(dadF, order=diffExpOrder)
  var j = diffprojectTAH(M,F)
  j := 1.0 + dadF * j
  let D = ln(determinant(j))
  D

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
      proc (l:T):T {.noinit.} = f(exp(l*t[a])*x),
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
  ## r_{ba} = - 2 ∂_{x_a} Tr[T^b f(x)† f(x_a T_a)]
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
      proc (l:T):V {.noinit.} = result.suToVec(fx.adj * f(l*t[a]+x)),
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
