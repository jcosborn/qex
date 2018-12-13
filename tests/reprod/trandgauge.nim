import qex

qexInit()
let
  lat = @[8,8,8,8]
  lo = lat.newLayout
var
  g = lo.newGauge
  r = lo.newRNGField RngMilc6

g.random

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
