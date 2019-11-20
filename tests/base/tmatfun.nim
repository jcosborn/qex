import qex
import physics/qcdTypes
import testutils

qexInit()
type
  MS[N:static[int]] = MatrixArray[N,N,SComplexV]
  MD[N:static[int]] = MatrixArray[N,N,DComplexV]

template chkzero(x: SomeFloat) =
  let e = epsilon(x)
  check(x < 20*e)

proc chkeq(x,y: any) =
  var mx = 0.0 * x[0,0].re
  var md = 0.0 * x[0,0].re
  for i in 0..<x.nrows:
    for j in 0..<x.ncols:
      let z = x - y
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
  chkzero(s)

proc rsqrtPH_test(x: any) =
  let y = x.adj * x
  #echo y
  #let z = y * y
  #let r = rsqrtPH(z)
  #let t = r * y
  let r = rsqrtPH(y)
  let t = r * y * r
  var o: type(t)
  o := 1
  chkeq(t, o)


var rs: RngMilc6
rs.seed(987654321, 1)

suite "Test matrix rsqrtPH":
  template trsqrtPH(n: typed) =
    var ms: MS[n]
    var md: MD[n]
    for i in 0..<simdLength(ms):
      gaussian( masked(ms,1 shl i), rs )
    test("n: " & $n & " single"):
      rsqrtPH_test(ms)
    for i in 0..<simdLength(md):
      gaussian( masked(md,1 shl i), rs )
    test("n: " & $n & " double"):
      rsqrtPH_test(md)
  trsqrtPH(1)
  trsqrtPH(2)
  trsqrtPH(3)
  #trsqrtPH(4)

qexFinalize()
