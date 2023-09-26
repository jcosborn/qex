import qex

qexInit()

proc testgauge(lat:seq[int]):seq[DLatticeColorMatrixV] =
  qexGC()
  let lo = lat.newLayout
  let g = lo.newGauge
  g.random
  qexGC()
  g

proc testfp(g:seq[DLatticeColorMatrixV]):float =
  let gc = GaugeActionCoeffs(plaq:6.0)
  qexGC()
  result = gc.gaugeAction1(g)
  qexGC()

proc runfp(f:proc, arg:auto):auto =
  f(arg)

proc testfun(g:auto):auto =
  runfp(testfp, g)

proc test(lat:auto):auto =
  var fail = 0
  qexGC()
  let g = runfp(testgauge,lat)
  let act = testfun(g)
  qexGC()
  for n in 0..<1024:
    qexGC()
    let a = testfun(g)
    if abs(a-act)/abs(act) > 1e-8:
      echo n," WRONG result: ",a," expecting: ",act
      inc fail
    if (n+1) mod 128 == 0:
      echo (n+1)," / 1024 failed ",fail
    qexGC()
  fail

var fail = 0
fail += test(@[8,8,8,8])
#fail += test(@[8,8,8,12])
fail += test(@[8,8,8,16])

if fail > 0:
  qexAbort()
else:
  qexFinalize()
