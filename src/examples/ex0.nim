import qex
qexInit()

var lat = [4,4,4,4]
var lo = newLayout(lat)
echoAll "Hello from rank ", lo.myRank, " of ", lo.nRanks

var v1 = lo.ColorVector()
var v2 = lo.ColorVector()
var m1 = lo.ColorMatrix()

threads:
  m1 := 1
  v1 := 2
  v2.even := m1 * v1
  v2.odd := 3
  echo "v2 even: ", v2.even.norm2
  echo "v2 odd: ", v2.odd.norm2

threads:
  shift(v1, dir=0, len=1, v2)
  shift(v2, dir=3, len=2, v1)
  echo "v2 even: ", v2.even.norm2
  echo "v2 odd: ", v2.odd.norm2

qexFinalize()
