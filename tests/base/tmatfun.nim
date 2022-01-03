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
  echo s
  chkzero(s, 256*x.nrows)

proc rsqrtPH_test(x: auto): auto =
  #echo x
  #let y = x.adj * x
  #echo y
  var xa: type x
  xa := x.adj
  let y = xa * x
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
  echo "x: ", x.norm2
  echo "y: ", y.norm2
  echo "t: ", t.norm2
  chkeq(t, o)

var rs: RngMilc6
rs.seed(13, 987654321)

suite "Test matrix rsqrtPH":
  template trsqrtPH(T: typedesc) =
    var m: T
    for i in 0..<simdLength(m):
      gaussian( masked(m,1 shl i), rs )
      #gaussian( m[asSimd(i)], rs )
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
  doTest(SimdS4)
  doTest(SimdD4)
  doTest(SimdS8)
  doTest(SimdD8)

qexFinalize()
