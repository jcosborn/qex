import base
import complexNumbers
import matrixConcept
import types
#import strformat

# inverse(M,S,M) -> M=S/M
# inverse(M,M,M) -> M=M/M
# solve(V,M,V) -> V=M\V
# solve(M,M,M) -> M=M\M

proc inverseN*(r: var Mat1, c: SomeNumber, x: Mat2) =
  echo "inverseN not implemented"
  doAssert(false)

proc inverse*(r: var Mat1, c: SomeNumber, x: Mat2) =
  const nc = r.nrows
  when nc==1:
    r := c / x[0,0]
  elif nc==2:
    let x00 = x[0,0]
    let x01 = x[0,1]
    let x10 = x[1,0]
    let x11 = x[1,1]
    let det = x00*x11 - x01*x10
    let idet = c / det
    r[0,0] :=  idet * x11
    r[0,1] := -idet * x01
    r[1,0] := -idet * x10
    r[1,1] :=  idet * x00
  elif nc==3:
    let x00 = x[0,0]
    let x01 = x[0,1]
    let x02 = x[0,2]
    let x10 = x[1,0]
    let x11 = x[1,1]
    let x12 = x[1,2]
    let x20 = x[2,0]
    let x21 = x[2,1]
    let x22 = x[2,2]
    let det0 = x00 * x11 - x01 * x10
    let det1 = x02 * x10 - x00 * x12
    let det2 = x01 * x12 - x02 * x11
    let det = det0*x22 + det1*x21 + det2*x20
    let idet = c / det
    r[0,0] := idet*(x11*x22-x12*x21)
    r[0,1] := idet*(x21*x02-x22*x01)
    r[0,2] := idet*det2
    r[1,0] := idet*(x12*x20-x10*x22)
    r[1,1] := idet*(x22*x00-x20*x02)
    r[1,2] := idet*det1
    r[2,0] := idet*(x10*x21-x11*x20)
    r[2,1] := idet*(x20*x01-x21*x00)
    r[2,2] := idet*det0
  elif nc==4:
    let x00 = x[0,0]
    let x01 = x[0,1]
    let x02 = x[0,2]
    let x03 = x[0,3]
    let x10 = x[1,0]
    let x11 = x[1,1]
    let x12 = x[1,2]
    let x13 = x[1,3]
    let x20 = x[2,0]
    let x21 = x[2,1]
    let x22 = x[2,2]
    let x23 = x[2,3]
    let x30 = x[3,0]
    let x31 = x[3,1]
    let x32 = x[3,2]
    let x33 = x[3,3]
    let s0 = x00 * x11 - x10 * x01
    let s1 = x00 * x12 - x10 * x02
    let s2 = x00 * x13 - x10 * x03
    let s3 = x01 * x12 - x11 * x02
    let s4 = x01 * x13 - x11 * x03
    let s5 = x02 * x13 - x12 * x03
    let c5 = x22 * x33 - x32 * x23
    let c4 = x21 * x33 - x31 * x23
    let c3 = x21 * x32 - x31 * x22
    let c2 = x20 * x33 - x30 * x23
    let c1 = x20 * x32 - x30 * x22
    let c0 = x20 * x31 - x30 * x21
    let det = s0*c5 - s1*c4 + s2*c3 + s3*c2 - s4*c1 + s5*c0
    #echo det
    let idet = c / det
    r[0,0] := ( x11 * c5 - x12 * c4 + x13 * c3) * idet
    r[0,1] := (-x01 * c5 + x02 * c4 - x03 * c3) * idet
    r[0,2] := ( x31 * s5 - x32 * s4 + x33 * s3) * idet
    r[0,3] := (-x21 * s5 + x22 * s4 - x23 * s3) * idet
    r[1,0] := (-x10 * c5 + x12 * c2 - x13 * c1) * idet
    r[1,1] := ( x00 * c5 - x02 * c2 + x03 * c1) * idet
    r[1,2] := (-x30 * s5 + x32 * s2 - x33 * s1) * idet
    r[1,3] := ( x20 * s5 - x22 * s2 + x23 * s1) * idet
    r[2,0] := ( x10 * c4 - x11 * c2 + x13 * c0) * idet
    r[2,1] := (-x00 * c4 + x01 * c2 - x03 * c0) * idet
    r[2,2] := ( x30 * s4 - x31 * s2 + x33 * s0) * idet
    r[2,3] := (-x20 * s4 + x21 * s2 - x23 * s0) * idet
    r[3,0] := (-x10 * c3 + x11 * c1 - x12 * c0) * idet
    r[3,1] := ( x00 * c3 - x01 * c1 + x02 * c0) * idet
    r[3,2] := (-x30 * s3 + x31 * s1 - x32 * s0) * idet
    r[3,3] := ( x20 * s3 - x21 * s1 + x22 * s0) * idet
  else:
    inverseN(r, c, x)

template inverse*(r: var Mat1, x: Mat2) = inverse(r, 1, x)

template `/`*(x: SomeNumber, y: Mat1): auto =
  var r {.noInit.}: MatrixArray[y.nrows,y.ncols,type(x/y[0,0])]
  inverse(r, x, y)
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
    inverse(m2, 1.5, m1)
    #echo "m2: ", m2
    m3 := (1.0/1.5)*(m1*m2)
    #echo "m3: ", m3
    let err = sqrt((1-m3).norm2)/N
    echo "err: ", err
    check(err, 5)

    #[
    projectU(m2, m1)
    m3 := m2.adj*m2
    var err2 = sqrt((1-m3).norm2/(N*N))
    echo " projectU err: ", err2
    doAssert(simdSum(err2)<simdLength(err2)*N*eps*30)
    #m2 := 0.1*(m2 - (trace(m2)/N))
    let m2n = 1/(10*sqrt(m2.norm2))
    m3 := exp(m2n*m2)
    m4 := exp(-m2n*m2)
    m2 := m3*m4
    #echo "exp ",m2,"\n\t= ",m3
    err2 = sqrt((1-m2).norm2/(N*N))
    echo " exp err: ", err2
    doAssert(simdSum(err2)<simdLength(err2)*N*eps*2)
    inverse(m4, m1)
    m3 := m1*m4
    echo m3
    err2 = sqrt((1-m3).norm2/(N*N))
    echo " inverse err: ", err2
    doAssert(simdSum(err2)<simdLength(err2)*N*eps*2)
    projectU(m2, m1)
    let r1 = trace(m4.adj*m2).re
    m3 := m1 + 3*seps*m1*m1
    #m3 := m1 + 1e-3'f32*m1*m1
    projectU(m2, m3)
    let r2 = trace(m4.adj*m2).re
    projectUderiv(m2, m1, m4)
    let dr = r2 - r1
    let dm = trace((m3-m1).adj * m2).re
    echo " r1: ", r1
    echo " r2: ", r2
    echo " dr: ", dr
    echo " dm: ", dm
    #echo "m1: ", m1
    #echo "m3: ", m3
    let dd = abs(dr - dm)
    echo " dd: ", dd
    doAssert(simdSum(dd)<simdLength(dd)*N*eps*40)
    ]#


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
