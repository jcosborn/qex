import base
#import basicOps
import complexNumbers
#import complexType
import matrixConcept
import types
import matinv
export matinv
import projUderiv

proc determinantN*(a: auto): auto =
  const nc = a.nrows
  var c {.noInit.}: type(a)
  var row: array[nc,int]
  var nswaps = 0
  var r: type(a[0,0])
  r := 1

  for i in 0..<nc:
    for j in 0..<nc:
      c[i,j] = a[i,j]
    row[i] = i

  for j in 0..<nc:
    if j>0:
      for i in j..<nc:
        var t2 = c[i,j]
        for k in 0..<j:
          t2 -= c[i,k] * c[k,j]
        c[i,j] := t2

    var rmax = c[j,j].norm2
    #[
    var kmax = j
    for k in (j+1)..<nc:
      var rn = c[k,j].norm2
      if rn>rmax:
        rmax = rn
        kmax = k
    if rmax==0: # matrix is singular
      r := 0
      return r
    if kmax != j:
      swap(row[j], row[kmax])
      inc nswaps
    ]#

    r *= c[j,j]

    let ri = 1.0/rmax
    var Cjji = ri * c[j,j]
    for i in (j+1)..<nc:
      var t2 = c[j,i]
      for k in 0..<j:
        t2 -= c[j,k] * c[k,i]
      c[j,i] := t2 * Cjji

  if (nswaps and 1) != 0:
    r := -r

  r

proc determinant*(x: auto): auto =
  assert(x.nrows == x.ncols)
  when x.nrows==1:
    result = x[0,0]
  elif x.nrows==2:
    optimizeAst:
      result = x[0,0]*x[1,1] - x[0,1]*x[1,0]
  elif x.nrows==3:
    result = (x[0,0]*x[1,1]-x[0,1]*x[1,0])*x[2,2] +
             (x[0,2]*x[1,0]-x[0,0]*x[1,2])*x[2,1] +
             (x[0,1]*x[1,2]-x[0,2]*x[1,1])*x[2,0]
  else:
    result = determinantN(x)

proc eigs3(e0,e1,e2: var auto; tr,p2,det: auto) =
  mixin sin,cos,acos
  let tr3 = (1.0/3.0)*tr
  let p23 = (1.0/3.0)*p2
  let tr32 = tr3*tr3
  let q = abs(0.5*(p23-tr32))
  let r = 0.25*tr3*(5*tr32-p2) - 0.5*det
  let sq = sqrt(q)
  let sq3 = q*sq
  #let rsq3 = r/sq3
  #var minv,maxv {.noinit.}:type(rsq3)
  #minv := -1.0
  #maxv := 1.0
  #let rsq3r = min(maxv, max(minv,rsq3))
  let isq3 = 1.0/sq3
  var minv,maxv {.noinit.}: type(isq3)
  maxv := 3e38
  minv := -3e38
  let isq3c = min(maxv, max(minv,isq3))
  let rsq3c = r * isq3c
  maxv := 1
  minv := -1
  let rsq3 = min(maxv, max(minv,rsq3c))
  let t = (1.0/3.0)*acos(rsq3)
  let st = sin(t)
  let ct = cos(t)
  let sqc = sq*ct
  let sqs = 1.73205080756887729352*sq*st  # sqrt(3)
  let ll = tr3 + sqc
  e0 = tr3 - 2*sqc
  e1 = ll + sqs
  e2 = ll - sqs

template rsqrtPHM2(r:typed; x:typed):untyped =
  let x00 = x[0,0].re
  let x11 = x[1,1].re
  let x01r = 0.5*(x[0,1].re+x[1,0].re)
  let x01i = 0.5*(x[0,1].im-x[1,0].im)
  let det = abs(x00*x11 - x01r*x01r - x01i*x01i)
  let tr = x00 + x11
  let sdet = sqrt(det)
  let trsdet = tr + sdet
  let c1 = 1/(sdet*sqrt(trsdet+sdet))
  let c0 = trsdet*c1
  r := c0 - c1*x

proc rsqrtPHM3f(c0,c1,c2:var auto; tr,p2,det:auto) =
  #[
  mixin sin,cos,acos
  let tr3 = (1.0/3.0)*tr
  let p23 = (1.0/3.0)*p2
  let tr32 = tr3*tr3
  let q = abs(0.5*(p23-tr32))
  let r = 0.25*tr3*(5*tr32-p2) - 0.5*det
  let sq = sqrt(q)
  let sq3 = q*sq
  #let rsq3 = r/sq3
  #var minv,maxv {.noinit.}:type(rsq3)
  #minv := -1.0
  #maxv := 1.0
  #let rsq3r = min(maxv, max(minv,rsq3))
  let isq3 = 1.0/sq3
  var minv,maxv {.noinit.}: type(isq3)
  maxv := 3e38
  minv := -3e38
  let isq3c = min(maxv, max(minv,isq3))
  let rsq3c = r * isq3c
  maxv := 1
  minv := -1
  let rsq3 = min(maxv, max(minv,rsq3c))
  let t = (1.0/3.0)*acos(rsq3)
  let st = sin(t)
  let ct = cos(t)
  let sqc = sq*ct
  let sqs = 1.73205080756887729352*sq*st  # sqrt(3)
  let l0 = tr3 - 2*sqc
  let ll = tr3 + sqc
  let l1 = ll + sqs
  let l2 = ll - sqs
  ]#
  var l0,l1,l2 {.noInit.}: type(tr)
  eigs3(l0,l1,l2, tr,p2,det)
  let sl0 = sqrt(abs(l0))
  let sl1 = sqrt(abs(l1))
  let sl2 = sqrt(abs(l2))
  let u = sl0 + sl1 + sl2
  let w = sl0 * sl1 * sl2
  let d = w*(sl0+sl1)*(sl0+sl2)*(sl1+sl2)
  let di = 1/d
  c0 = (w*u*u+l0*sl0*(l1+l2)+l1*sl1*(l0+l2)+l2*sl2*(l0+l1))*di
  c1 = -(tr*u+w)*di
  c2 = u*di

template rsqrtPHM3(r:typed; x:typed):untyped =
  let tr = trace(x).re
  let x2 = x*x
  let p2 = trace(x2).re
  let det = determinant(x).re
  var c0,c1,c2:type(tr)
  rsqrtPHM3f(c0, c1, c2, tr, p2, det)
  r := c0 + c1*x + c2*x2

template rsqrtPHMN(r:typed; x:typed):untyped =
  let xi = 1/x
  let xi2 = xi.norm2
  let xit = trace(xi).re
  let ds = xit/xi2
  #var ds = x.norm2.simdMax
  #ds = 0.5*sqrt(ds)
  #echo "ds: ", ds

  var e = (0.5*ds)*xi - 0.5
  var s = 1 + e
  #echo "e: ", e
  #echo "s: ", s

  let estop = epsilon(ds.simdMax)^2
  let maxit = 20
  var nit = 0
  while true:
    inc nit
    #let t = (e/s) * e
    #let t = e * (s \ e)
    let si = 1/s
    let t = e * (si * e)
    e := -0.5 * t
    s += e
    let enorm = e.norm2.simdMax
    #echo nit, " enorm: ", enorm
    if nit>=maxit or enorm<estop: break
  let sds = 1/sqrt(ds)
  r := sds*s

# Bini (https://arxiv.org/pdf/1703.02456.pdf)
template rsqrtPHMN2(r:typed; x:typed):untyped =
  let xn = x.norm2
  let ds = sqrt(xn)
  let dsi = 3/ds
  #echo "ds: ", ds

  var a = dsi * x
  var b {.noInit.} :type(r)
  b := 1
  #var b = (1.4/xn)*x

  let estop = epsilon(ds.simdMax)^2
  let maxit = 20
  var nit = 0
  while true:
    let e = 1 - b*a*b
    let enorm = e.norm2.simdMax
    echo nit, " enorm: ", enorm
    if nit>=maxit or enorm<estop: break
    inc nit
    let t = b*e
    #let t2 = t.norm2
    #let c = 0.5/sqrt(t2)
    let c = 0.5
    b += c*t
  let sds = sqrt(dsi)
  r := sds*b
  #r := b

# -0.5[B+((1-BBA)^3-1)/(BA)] = -0.5B[1-3+3BBA-BBBBAA] = 0.5B[2-3BBA+BBBBAA]
# = 0.5B[1-BBA][2-BBA]
template rsqrtPHMN3(r:typed; x:typed):untyped =
  let xn = x.norm2
  let ds = sqrt(xn)
  let dsi = 1/ds
  #echo "ds: ", ds

  var a = dsi * x
  var b {.noInit.} :type(r)
  b := 1
  #var b = (1.4/xn)*x

  let estop = epsilon(ds.simdMax)^2
  let maxit = 20
  var nit = 0
  while true:
    let e = 1 - b*a*b
    let enorm = e.norm2.simdMax
    echo nit, " enorm: ", enorm
    if nit>=maxit or enorm<estop: break
    inc nit
    let t = b*e*(1+e)
    #let t2 = t.norm2
    #let c = 0.5/sqrt(t2)
    let c = 0.5
    b += c*t
  let sds = sqrt(dsi)
  r := sds*b
  #r := b


template rsqrtPHM(r:typed; x:typed):untyped =
  mixin rsqrt, nrows
  assert(r.nrows == x.nrows)
  assert(r.ncols == x.ncols)
  assert(r.nrows == r.ncols)
  when r.nrows==1:
    let t = rsqrt(x[0,0].re)
    r := t
  elif r.nrows==2:
    rsqrtPHM2(r, x)
  elif r.nrows==3:
    rsqrtPHM3(r, x)
  else:
    rsqrtPHMN(r, x)
    #rsqrtPHMN2(r, x)
    #rsqrtPHMN3(r, x)
proc rsqrtPH*(r: var Mat1; x: Mat2) = rsqrtPHM(r, x)
template rsqrtPH*[T:Mat1](x: T): T =
  var r {.noInit.}: T
  rsqrtPH(r, x)
  r

# x (x'x)^{-1/2}
proc projectU*(r: var Mat1; x: Mat2) =
  let t = x.adj * x
  var t2{.noInit.}: type(t)
  rsqrtPH(t2, t)
  mul(r, x, t2)

# (d/dX') Tr(U'C+C'U) / 2 = (d/dX') Tr(X'CZ+C'XZ) / 2
# = CZ - (1/2) < Z (X'C + C'X) Z (dY/dX') >
# (dY/dX') Y + Y (dY/dX') = 2X
# Z(dY/dX') Y + Y Z(dY/dX') = 2ZX
# S Y + Y S = U' C Z = Z X' C Z
# cz-xz^3(x'c+c'x)/2
# CH: 4528 flops
proc projectUderiv*(r: var Mat1, u: Mat2, x: Mat3, chain: Mat4) =
  # U = X (X'X)^{-1/2} = (XX')^{-1/2} X
  # Y = sqrt(X'X)
  # Z = (X'X)^{-1/2}
  # F = C Z - z (Cd U + Ud C) z (dY/dX)
  var y, z, t1, t2: Mat1
  y := x.adj * u
  inverse(z, y)
  #echo "inverse: ", z
  #QLA_M_eq_M_times_M(d, c, &z);
  r := chain * z
  #QLA_M_eq_Ma_times_M(&t1, p, d);
  t1 := u.adj * r
  #sylsolve_site(NCARG &t2, &y, &y, &t1);
  sylsolve(t2, y, t1)
  #QLA_M_eq_M(&t1, &t2);
  #QLA_M_peq_Ma(&t1, &t2);
  #QLA_M_eq_M_times_M(&t2, m, &t1);
  #QLA_M_meq_M(d, &t2);
  t1 := t2 + t2.adj
  r -= x * t1

proc projectUderiv*(r: var Mat1, x: Mat2, c: Mat3) =
  var u {.noInit.}: type(r)
  projectU(u, x)
  #echo u, x, c
  projectUderiv(r, u, x, c)

proc projectSU*(r: var Mat1; x: Mat2) =
  const nc = r.nrows
  var m{.noinit.}: type(r)
  m.projectU x
  var d = m.determinant    # already unitary: 1=|d
  let p = (1.0/float(-nc)) * atan2(d.im, d.re)
  d.re = cos p
  d.im = sin p
  r := d * m

proc projectTAH*(r: var Mat1; x: Mat2) =
  r := 0.5*(x-x.adj)
  const nc = x.nrows
  when nc > 1:
    let d = r.trace / nc.float
    r -= d

proc checkU*(x: Mat1): auto {.inline, noinit.} =
  ## Returns the sum of deviations of x^dag x and det(x) from unitarity.
  var d = norm2(-1.0 + x.adj * x)
  return d

proc checkSU*(x: Mat1): auto {.inline, noinit.} =
  ## Returns the sum of deviations of x^dag x and det(x) from unitarity.
  var d = norm2(-1.0 + x.adj * x)
  d += norm2(-1.0 + x.determinant)
  return d

discard """
template rsqrtM2(r:typed; x:typed):untyped =
  load(x00, x[0,0].re)
  load(x01, x[0,1])
  #load(x10, x[1,0])
  load(x11, x[1,1].re)
  let det := a00*a11 -
  QLA_r_eq_Re_c_times_c (det, a00, a11);
  QLA_r_meq_Re_c_times_c(det, a01, a10);
  tr = QLA_real(a00) + QLA_real(a11);
  sdet = sqrtP(fabsP(det));
  // c0 = (l2/sl1-l1/sl2)/(l2-l1) = (l2+sl1*sl2+l1)/(sl1*sl2*(sl1+sl2))
  // c1 = (1/sl2-1/sl1)/(l2-l1) = -1/(sl1*sl2*(sl1+sl2))
  c1 = 1/(sdet*sqrtP(fabsP(tr+2*sdet)));
  c0 = (tr+sdet)*c1;
  c1 = -c1;
  // c0 + c1*a
  QLA_c_eq_c_times_r_plus_r(QLA_elem_M(*r,0,0), a00, c1, c0);
  QLA_c_eq_c_times_r(QLA_elem_M(*r,0,1), a01, c1);
  QLA_c_eq_c_times_r(QLA_elem_M(*r,1,0), a10, c1);
  QLA_c_eq_c_times_r_plus_r(QLA_elem_M(*r,1,1), a11, c1, c0);

template rsqrtM(r:typed; x:typed):untyped =
  assert(r.nrows == x.nrows)
  assert(r.ncols == x.ncols)
  assert(r.nrows == r.ncols)
  if r.nrows==1:
    rsqrt(r[0,0], x[0,0])
  elif r.nrows==2:
    rsqrtM2(r, x)
  elif r.nrows==3:
    rsqrtM3(r, x)
  else:
    echo "unimplemented"
    quit(1)
proc rsqrt(r:var Mat1; x:Mat2) = rsqrt(r, x)
"""

proc exp*(m: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  when m.nrows == 1:
    r := exp(m[0,0])
  else:
    type ft = numberType(m)
    template term(n,x: typed): untyped =
      when x.type is nil.type: 1 + ft(n)*m
      else: 1 + ft(n)*m*x
    #template r3:untyped = nil
    let r12 = term(1.0/12.0, nil)
    let r11 = term(1.0/11.0, r12)
    let r10 = term(1.0/10.0, r11)
    let r9 = term(1.0/9.0, r10)
    let r8 = term(1.0/8.0, r9)
    let r7 = term(1.0/7.0, r8)
    let r6 = term(1.0/6.0, r7)
    let r5 = term(1.0/5.0, r6)
    let r4 = term(1.0/4.0, r5)
    let r3 = term(1.0/3.0, r4)
    let r2 = term(1.0/2.0, r3)
    r := 1 + m*r2
  r
proc ln*(m: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  when m.nrows == 1:
    r := ln(m[0,0])
  else:
    static: error("ln of matrix not implimented.")
  r

proc re*(m: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  for i in 0..<m.nrows:
    for j in 0..<m.ncols:
      r[i,j] := re(m[i,j])
  r
proc im*(m: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  for i in 0..<m.nrows:
    for j in 0..<m.ncols:
      r[i,j] := im(m[i,j])
  r

when isMainModule:
  import macros
  import simd
  template `+`(x: SimdS4): untyped = x
  template `+`(x: SimdS4, y: ComplexType): untyped =
    asReal(x) + y
  template `-`(x: SimdS4, y: ComplexType): untyped =
    asReal(x) - y
  template `*`(x: SimdS4, y: ComplexType): untyped =
    asReal(x) * y
  template `*`(x: ComplexType, y: SimdS4): untyped =
    x * asReal(y)
  template `/`(x: SomeFloat, y: SomeInteger): untyped = x/(type(x))(y)
  template add(r: ComplexType, x: SimdS4, y: ComplexType): untyped =
    add(r, asReal(x), y)
  template sub(r: ComplexType, x: SimdS4, y: ComplexType): untyped =
    sub(r, asReal(x), y)
  template mul(r: ComplexType, x: SimdS4, y: ComplexType): untyped =
    mul(r, asReal(x), y)
  template check(x:untyped, n:SomeNumber):untyped =
    let r0 = x
    let r = simdSum(r0)/simdLength(r0)
    echo "error/eps: ", r/epsilon(r)
    doAssert(abs(r)<n*epsilon(r))
  proc test(T: typedesc) =
    var m1,m2,m3,m4: T
    let N = m1.nrows
    for i in 0..<N:
      for j in 0..<N:
        let fi = i.float
        let fj = j.float
        m1[i,j].re := 0.5 + 0.7/(0.9+1.3*fi-fj)
        m1[i,j].im := 0.1 + 0.3/(0.4+fi-1.1*fj)
    echo "test " & $N & " " & $T
    #echo "m1: ", m1
    m2 := m1.adj * m1
    #echo m2
    rsqrtPH(m3, m2)
    #echo m3
    m4 := m3*m2*m3
    let err = sqrt((1-m4).norm2)/N
    #echo m4
    echo " rsqrtPH err: ", err
    check(err, 50)

    projectU(m2, m1)
    m3 := m2.adj*m2
    var err2 = sqrt((1-m3).norm2)/N
    echo " projectU err: ", err2
    check(err, 50)

    #m2 := 0.1*(m2 - (trace(m2)/N))
    let m2n = 1/(10*sqrt(m2.norm2))
    m3 := exp(m2n*m2)
    m4 := exp(-m2n*m2)
    m2 := m3*m4
    #echo "exp ",m2,"\n\t= ",m3
    err2 = sqrt((1-m2).norm2/(N*N))
    echo " exp err: ", err2
    check(err2, 5)

    inverse(m4, m1)
    m3 := m1*m4
    #echo m3
    err2 = sqrt((1-m3).norm2/(N*N))
    echo " inverse err: ", err2
    check(err2, 5)

    if N<5:
      projectU(m2, m1)
      let r1 = trace(m4.adj*m2).re
      let seps = sqrt(epsilon(simdSum(r1)))
      m3 := m1 + 3*seps*m1*m1
      #m3 := m1 + 1e-3'f32*m1*m1
      projectU(m2, m3)
      let r2 = trace(m4.adj*m2).re
      projectUderiv(m2, m1, m4)
      let dr = r2 - r1
      let dm = trace((m3-m1).adj * m2).re
      #echo " r1: ", r1
      #echo " r2: ", r2
      #echo " dr: ", dr
      #echo " dm: ", dm
      #echo "m1: ", m1
      #echo "m3: ", m3
      let dd = abs(dr - dm)
      echo " projectUderiv err: ", dd
      #doAssert(simdSum(dd)<simdLength(dd)*N*eps*40)
      check(dd, 20*N)


  type
    Cmplx[T] = ComplexType[T]
    CM[N:static[int],T] = MatrixArray[N,N,Cmplx[T]]
  template doTest(t:untyped) =
    when declared(t):
      test(CM[1,t])
      test(CM[2,t])
      test(CM[3,t])
      test(CM[4,t])
  doTest(float32)
  doTest(float64)
  doTest(SimdS4)
  doTest(SimdD4)
  doTest(SimdS8)
  doTest(SimdD8)
