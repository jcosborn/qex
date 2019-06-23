import qex

proc plaq*(g: seq[Field]): float =
  let lo = g[0].l
  let nd = lo.nDim
  var m = lo.ColorMatrix()
  var tr: type(trace(m))
  let t = newTransporters(g, g[0], 1)
  threads:
    m := 0
    for mu in 1..<nd:
      for nu in 0..<mu:
        m += (t[mu]^*g[nu]) * (t[nu]^*g[mu]).adj
    tr = trace(m)
  result = tr.re/(0.5*float(nd*(nd-1)*lo.physVol))

qexInit()

var lat = intSeqParam("lat", @[4,4,4,4])
var lo = newLayout(lat)
var g = lo.newGauge()
threads:
  g.unit
  #g.random

echo "average plaq: ", plaq(g)

qexFinalize()
