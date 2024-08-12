import qex, gauge/stoutsmear

qexInit()

let
  seed = 4321u
  lat = @[8,8,8,16]
  nd = lat.len
  lo = lat.newLayout
  vol = lo.physVol
  gc = GaugeActionCoeffs(plaq: 6.0)

var r = lo.newRNGField(RngMilc6, seed)

var gaugeu1 = newseq[Field[VLEN,Color[MatrixArray[1,1,DComplexV]]]](nd)
var gaugesu2 = newseq[Field[VLEN,Color[MatrixArray[2,2,DComplexV]]]](nd)
var gaugesu3= newseq[Field[VLEN,Color[MatrixArray[3,3,DComplexV]]]](nd)
for i in 0..<nd:
  gaugeu1[i].new lo
  gaugesu2[i].new lo
  gaugesu3[i].new lo

var fu1 = gaugeu1.newOneOf
var fsu2 = gaugesu2.newOneOf
var fsu3 = gaugesu3.newOneOf

threads:
  gaugeu1.random r
  gaugesu2.random r
  gaugesu3.random r

var eps = 0.1
for n in 0..<3:
  gc.gaugeForce(gaugeu1, fu1)
  gc.gaugeForce(gaugesu2, fsu2)
  gc.gaugeForce(gaugesu3, fsu3)
  threads:
    for mu in 0..<nd:
      for s in gaugeu1[mu]:
        gaugeu1[mu][s] := exp((-eps)*fu1[mu][s])*gaugeu1[mu][s]
      for s in gaugesu2[mu]:
        gaugesu2[mu][s] := exp((-eps)*fsu2[mu][s])*gaugesu2[mu][s]
      for s in gaugesu3[mu]:
        gaugesu3[mu][s] := exp((-eps)*fsu3[mu][s])*gaugesu3[mu][s]

echo gc.gaugeAction1(gaugeu1)
echoPlaq gaugeu1
echo gc.gaugeAction1(gaugesu2)
echoPlaq gaugesu2
echo gc.gaugeAction1(gaugesu3)
echoPlaq gaugesu3


proc fnorm2(f: auto): float =
  var s:float
  threads:
    var t = 0.0
    for i in 0..<nd:
      t += f[i].norm2
    threadSingle:
      s = t
  s

echo gaugeu1.fnorm2-float(4*vol)
echo gaugesu2.fnorm2-float(2*4*vol)
echo gaugesu3.fnorm2-float(3*4*vol)

echo fu1.fnorm2
echo fsu2.fnorm2
echo fsu3.fnorm2

qexFinalize()
