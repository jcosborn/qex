import base
#import basicOps
import complexConcept
#import complexType
import matrixConcept
import types

proc determinant*(x:any):auto =
  assert(x.nrows == x.ncols)
  if x.nrows==1:
    result = x[0,0]
  elif x.nrows==2:
    result = x[0,0]*x[1,1] - x[0,1]*x[1,0]
  elif x.nrows==3:
    result = (x[0,0]*x[1,1]-x[0,1]*x[1,0])*x[2,2] +
             (x[0,2]*x[1,0]-x[0,0]*x[1,2])*x[2,1] +
             (x[0,1]*x[1,2]-x[0,2]*x[1,1])*x[2,0]
  else:
    echo "unimplemented"
    quit(1)

template rsqrtPHM2(r:typed; x:typed):untyped =
  let x00 = x[0,0].re
  let x11 = x[1,1].re
  let x01 = 0.5*(x[0,1]+adj(x[1,0]))
  let det = abs(x00*x11 - x01.norm2)
  let tr = x00 + x11
  let sdet = sqrt(det)
  let trsdet = tr + sdet
  let c1 = 1/(sdet*sqrt(trsdet+sdet))
  let c0 = trsdet*c1
  r := c0 - c1*x

proc rsqrtPHM3f(c0,c1,c2:var any; tr,p2,det:any) =
  mixin sin,cos,acos
  let tr3 = (1.0/3.0)*tr
  let p23 = (1.0/3.0)*p2
  let tr32 = tr3*tr3
  let q = abs(0.5*(p23-tr32))
  let r = 0.25*tr3*(5*tr32-p2) - 0.5*det
  let sq = sqrt(q)
  let sq3 = q*sq
  let rsq3 = r/sq3
  let t = (1.0/3.0)*acos(rsq3)
  let st = sin(t)
  let ct = cos(t)
  let sqc = sq*ct
  let sqs = 1.73205080756887729352*sq*st  # sqrt(3)
  let l0 = tr3 - 2*sqc
  let ll = tr3 + sqc
  let l1 = ll + sqs
  let l2 = ll - sqs
  #echo l0
  #echo l1
  #echo l2
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
  var ds = x.norm2
  #if ds == 0.0
  #    M_eq_d(r, 1./0.);
  #    return;
  #  }
  ds = sqrt(ds)

  var e = (0.5*ds)/x - 0.5
  var s = 1 + e

  let estop = epsilon(ds)
  let maxit = 20
  var nit = 0
  while true:
    inc nit
    #let t = (e/s) * e
    #e = -0.5 * t
    let t = e * (s \ e)
    e := -0.5 * t
    s += e
    let enorm = e.norm2
    #//printf("%i enorm = %g\n", nit, enorm);
    if nit>=maxit or enorm<estop: break
  r := x/sqrt(ds)

template rsqrtPHM(r:typed; x:typed):untyped =
  mixin rsqrt
  assert(r.nrows == x.nrows)
  assert(r.ncols == x.ncols)
  assert(r.nrows == r.ncols)
  if r.nrows==1:
    rsqrt(r[0,0].re, x[0,0].re)
    assign(r[0,0].im, 0)
  elif r.nrows==2:
    rsqrtPHM2(r, x)
  elif r.nrows==3:
    rsqrtPHM3(r, x)
  else:
    echo "unimplemented"
    quit(1)
proc rsqrtPH(r:var Mat1; x:Mat2) = rsqrtPHM(r, x)

proc projectU*(r:var Mat1; x:Mat2) =
  let t = x.adj * x
  var t2{.noInit.}:type(t)
  rsqrtPH(t2, t)
  mul(r, x, t2)

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
  type ft = numberType(m)
  template term(n,x: typed): untyped =
    when x.type is nil.type: 1 + ft(n)*m
    else: 1 + ft(n)*m*x
  #template r3:untyped = nil
  let r8 = term(1.0/8.0, nil)
  let r7 = term(1.0/7.0, r8)
  let r6 = term(1.0/6.0, r7)
  let r5 = term(1.0/5.0, r6)
  let r4 = term(1.0/4.0, r5)
  let r3 = term(1.0/3.0, r4)
  let r2 = term(1.0/2.0, r3)
  r := 1 + m*r2
  r

when isMainModule:
  import macros
  import simd
  template makeTest2(n,f:untyped):untyped =
    proc f[T]:auto =
      const N = n
      type
        Cmplx = ComplexType[T]
        M2 = MatrixArray[N,N,Cmplx]
      var m1,m2,m3,m4:M2
      for i in 0..<N:
        for j in 0..<N:
          #if i==j:
            m1[i,j] = cast[Cmplx](((0.5+i+j-i*j).to(T),(i-j-i*i+j*j).to(T)))
      m2 := m1.adj * m1
      #echo m2
      rsqrtPH(m3, m2)
      #echo m3
      m4 := m3*m2*m3
      let err = sqrt((1-m4).norm2/(N*N))
      echo "test " & $N & " " & $T
      #echo m4
      echo err
      result = err
      projectU(m2, m1)
      m3 := m2.adj*m2
      let err2 = sqrt((1-m3).norm2/(N*N))
      echo "err2: ", err2
      m3 := exp(m2)
  macro makeTest(n:untyped):auto =
    let f = ident("test" & n.repr)
    result = quote do: makeTest2(`n`,`f`)
  makeTest(1)
  makeTest(2)
  makeTest(3)
  block:
    template check(x:untyped):untyped =
      let r = x
      echo "error/eps: ", r/epsilon(r)
      doAssert(abs(r)<128*epsilon(r))
    check(test1[float32]())
    #check(test2[float32]())
    #check(test3[float32]())
    #check(test1[float64]())
    #check(test2[float64]())
    #check(test3[float64]())
  block:
    template check(x:untyped):untyped =
      let r0 = x
      let r = simdReduce(r0)/simdLength(r0)
      echo "error/eps: ", r/epsilon(r)
      doAssert(abs(r)<64*epsilon(r))
    template doTest(t:untyped) =
      when declared(t):
        check(test1[t]())
        check(test2[t]())
        check(test3[t]())
    #doTest(SimdS4)
    #doTest(SimdD4)
    #doTest(SimdS8)
    #doTest(SimdD8)
