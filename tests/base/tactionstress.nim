import qex
import os

qexInit()

let
  fn = if paramCount() > 0: paramStr 1 else: ""
  lat = if fn.len == 0: @[8,8,8,8] else: fn.getFileLattice
  lo = lat.newLayout
var g = lo.newGauge
if fn.len == 0:
  g.random
elif 0 != g.loadGauge fn:
  echo "ERROR: couldn't load gauge file: ",fn
  qexFinalize()
  quit(-1)

proc testfp(g:seq[DLatticeColorMatrixV]):float =
  let gc = GaugeActionCoeffs(plaq:6.0)
  qexGC()
  result = gc.gaugeAction1(g)
  qexGC()

proc runfp(f:proc, arg:auto):auto =
  f(arg)

proc testfun(g:auto):auto =
  runfp(testfp, g)

proc test(g:auto):auto =
  var fail = 0
  qexGC()
  let act = testfun(g)
  qexGC()
  for n in 0..<1024:
    qexGC()
    let a = testfun(g)
    if abs(a-act)/abs(act) > 1e-8:
      echo "Wrong result: ",a
      inc fail
    if n mod 128 == 0:
      echo n," / 1024 failed ",fail
    qexGC()
  fail

if test(g) > 0:
  qexAbort()
else:
  qexFinalize()
