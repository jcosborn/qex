import qex
import physics/qcdTypes
import testutils

qexInit()
type
  Cmplx[T] = ComplexType[T]
  RM[N:static[int],T] = MatrixArray[N,N,T]
  CM[N:static[int],T] = MatrixArray[N,N,Cmplx[T]]

template chkzero(x: SomeFloat, n: SomeNumber): untyped =
  let e = epsilon(x)
  check(x < n*e)

proc chkeq(x,y: auto): auto =
  let z = x - y
  var mx,md:type(x[0,0].re)
  for i in 0..<x.nrows:
    for j in 0..<x.ncols:
      mx = max(mx, abs(x[i,j].re))
      mx = max(mx, abs(x[i,j].im))
      md = max(md, abs(z[i,j].re))
      md = max(md, abs(z[i,j].im))
  let r = md / mx
  let rs = r.simdSum
  let s = rs / type(rs)(r.simdLength * x.nrows * x.ncols)
  #echo md
  #echo mx
  #echo s
  chkzero(s, 384*x.nrows)

proc rsqrtPH_test(x: auto): auto =
  #echo x
  #let y = x.adj * x
  #echo y
  var xa: type x
  xa := x.adj
  let y = xa * x
  #let y = x.adj * x
  #let xa = x.adj
  #let y2 = xa * x
  #echo y2
  #let z = y * y
  #let r = rsqrtPH(z)
  #let t = r * y
  let r = rsqrtPH(y)
  let t = r * y * r
  var o: type(t)
  o := 1
  #echo "x: ", x.norm2
  #echo "y: ", y.norm2
  #echo "t: ", t.norm2
  chkeq(t, o)

var rs: RngMilc6
rs.seed(13, 987654321)

suite "Test matrix rsqrtPH":
  #template trsqrtPH(T: typedesc) =
  proc trsqrtPH(T: typedesc) =
    var m: T
    for i in 0..<simdLength(m):
      when type(m[0,0].re) is SomeFloat:
        gaussian( m, rs )
      else:
        #gaussian( masked(m,1 shl i), rs )
        gaussian( m[asSimd(i)], rs )
      #m := 1
    test("rsqrtPH " & $m.type):
      subtest rsqrtPH_test(m)
  template doTest(t:untyped) =
    when declared(t):
      trsqrtPH(RM[1,t])
      trsqrtPH(RM[2,t])
      trsqrtPH(RM[3,t])
      trsqrtPH(RM[4,t])
      trsqrtPH(CM[1,t])
      trsqrtPH(CM[2,t])
      trsqrtPH(CM[3,t])
      trsqrtPH(CM[4,t])
  doTest(float32)
  doTest(float64)
  #doTest(SimdS1)
  #doTest(SimdD1)
  doTest(SimdS2)
  doTest(SimdD2)
  doTest(SimdS4)
  doTest(SimdD4)
  doTest(SimdS8)
  doTest(SimdD8)

proc solveLRTest(T: typedesc) =
  type A = MatrixArray[8,8,T]
  test("paired left/right solve " & $T):
    var a, a0, l, l0, r, r0:A
    for i in 0..<8:
      for j in 0..<8:
        when T is SomeFloat:
          a[i, j] := (if i == j: 3.0 + 0.1*i.float else: 0.01*(1 + 3*i - 2*j).float)
          l[i, j] := 0.02*(1 + 2*i + 3*j).float
          r[i, j] := 0.03*(2 - 3*i + j).float
        else:
          for k in 0..<simdLength(a[i, j]):
            let s = 1.0 + 0.01*k.float
            a[i, j][asSimd(k)] = s*(if i == j: 3.0 + 0.1*i.float else: 0.01*(1 + 3*i - 2*j).float)
            l[i, j][asSimd(k)] = s*0.02*(1 + 2*i + 3*j).float
            r[i, j][asSimd(k)] = s*0.03*(2 - 3*i + j).float
    a0 := a
    l0 := l
    r0 := r
    a.solveLRNoPivot(l, r)
    let
      el = sqrt((a0*l - l0).norm2.simdMax)
      er = sqrt((r*a0 - r0).norm2.simdMax)
    check el < 2e-12
    check er < 2e-12

suite "Test matrix paired solves":
  solveLRTest(float64)
  when declared(SimdD4):
    solveLRTest(SimdD4)

proc detNoPivotTest(T: typedesc) =
  type A = MatrixArray[8,8,T]
  test("unpivoted determinant " & $T):
    var a:A
    for i in 0..<8:
      for j in 0..<8:
        when T is SomeFloat:
          a[i, j] := (if i == j: 3.0 + 0.1*i.float else: 0.01*(1 + 3*i - 2*j).float)
        else:
          for k in 0..<simdLength(a[i, j]):
            let s = 1.0 + 0.01*k.float
            a[i, j][asSimd(k)] = s*(if i == j: 3.0 + 0.1*i.float else: 0.01*(1 + 3*i - 2*j).float)
    let
      x = detNoPivot(a)
      y = determinant(a)
      e = sqrt(norm2(x - y).simdMax)
      s = max(1.0, sqrt(norm2(y).simdMax))
    check e/s < 2e-12

suite "Test unpivoted determinant":
  detNoPivotTest(float64)
  when declared(SimdD4):
    detNoPivotTest(SimdD4)

proc expAHTest(T: typedesc) =
  type M = CM[3,T]
  test("adaptive AH exponential " & $T):
    for a in [0.0, 0.14, 0.16, 0.5, 2.0, 8.0]:
      var s:T
      when T is SomeFloat:
        s = a
      else:
        var v: array[simdLength(s), float]
        for k in 0..<v.len:
          v[k] = a*(1.0 + 0.1*k.float)
        s := v
      var m:M
      m := 0
      m[0, 0].im := 0.3*s
      m[1, 1].im := -0.1*s
      m[2, 2].im := -0.2*s
      m[0, 1].re := s
      m[1, 0].re := -s
      m[1, 2].im := 0.5*s
      m[2, 1].im := 0.5*s
      let
        r = expAH(m)
        er = exp(m)
        de = sqrt((r - er).norm2.simdMax)
        du = sqrt((r.adj*r - 1.0).norm2.simdMax)
        dd = sqrt(norm2(determinant(r) - 1.0).simdMax)
      check de < 2e-12
      check du < 2e-12
      check dd < 2e-12

proc expAH1Test(T: typedesc) =
  type M = CM[1,T]
  test("adaptive U(1) exponential " & $T):
    for a in [0.0, 0.14, 0.5, 2.0, 8.0]:
      var s:T
      when T is SomeFloat:
        s = a
      else:
        var v: array[simdLength(s), float]
        for k in 0..<v.len:
          v[k] = a*(1.0 + 0.1*k.float)
        s := v
      var m:M
      m := 0
      m[0, 0].im := s
      let
        r = expAH(m)
        er = exp(m)
        de = sqrt((r - er).norm2.simdMax)
        du = sqrt((r.adj*r - 1.0).norm2.simdMax)
      check de < 2e-12
      check du < 2e-12

suite "Test adaptive AH exponential":
  expAH1Test(float64)
  expAHTest(float64)
  when declared(SimdD4):
    expAH1Test(SimdD4)
    expAHTest(SimdD4)

qexFinalize()
