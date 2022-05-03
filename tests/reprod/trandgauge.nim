import qex

qexInit()
let
  lat = @[8,8,8,8]
  lo = lat.newLayout
var
  g = lo.newGauge
  #r = lo.newRNGField RngMilc6

g.random
echo g[0].norm2

when defined(FUELCompat):
  const P = @[0.0009367637705385087, -0.0008371956818981114, -0.0004590285518363084, 0.0005837830782273779, 0.0003457537836297183, 0.0005113416413136034]
else:
  const P = @[0.0006005738094166639, 0.0007744149733359666, 0.000491692592364555, -0.0002244585371871249, -0.000700363878755635, -4.121898341926528e-05]

let pl = g.plaq
var d = 0.0
for i in 0..<P.len:
  let s = P[i]-pl[i]
  d += s*s
echo "diff2: ",d
if d > 1e-30:
  echo "Failed"
  echo "Expecting ",P
  echo "Actual    ",pl
  qexExit 1
else:
  echo "Passed"
  qexFinalize()
