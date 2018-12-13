import qex
import gauge/gaugefix
import ../base/testutils

qexInit()
let defaultLat = @[8,8,8,8]
defaultSetup()
var t = lo.ColorMatrix()
var t1 = lo.ColorMatrix()
var t2 = lo.ColorMatrix()
var g1 = lo.newGauge()
var g2 = lo.newGauge()
let lo1 = newLayout(lo.physGeom, 1)
var rnd = newRNGField(RngMilc6, lo1)

template getDiff(g,g2: any): float =
  var d = 0.0
  for i in 0..<g.len:
    d += norm2(g[i]-g2[i])
  d/g.len.float

suite "Test Gauge Fixing":
  test "gauge transform from unity":
    threads:
      g.unit
      g1.random rnd
      g2.random rnd
      t1.randomU rnd
      t2 := t1.adj
    let p0 = g.plaq
    g1.gaugeTransform(g, t1)
    let p1 = g1.plaq
    CT = 1e-13
    check(p0 ~ p1)
    g2.gaugeTransform(g1, t2)
    let d = getDiff(g, g2).sqrt
    echo "RMS diff: ", d
    check(d <= 1e-10)

  test "gauge transform from random":
    threads:
      g.random rnd
      g1.random rnd
      g2.random rnd
      t1.randomU rnd
      t2 := t1.adj
    let p0 = g.plaq
    g1.gaugeTransform(g, t1)
    let p1 = g1.plaq
    CT = 1e-10
    check(p0 ~ p1)
    g2.gaugeTransform(g1, t2)
    let d = getDiff(g, g2).sqrt
    echo "RMS diff: ", d
    check(d <= 1e-9)

  test "gauge fix from transform of unity":
    threads:
      g.unit
      g1.random rnd
      g2.random rnd
      t1.randomSU rnd
    let p0 = g.plaq
    g1.gaugeTransform(g, t1)
    var gstop = 1e-10
    var orf = 1.5
    let tdirs = @[ @[0,1,2,3], @[0,1,2], @[3], @[1,2,3], @[0,2,3] ]
    for dirs in tdirs:
      echo dirs
      let l0 = g.linkTrace dirs
      threads: t := 1
      getGaugeFixTransform(t, g1, dirs, gstop, orf)
      g2.gaugeTransform(g1, t)
      let p2 = g2.plaq
      let l2 = g2.linkTrace dirs
      CT = 1e-11
      check(p0 ~ p2)
      check(l0 ~ l2)

  test "gauge fix from random":
    threads:
      g.random rnd
      g1.random rnd
      g2.random rnd
      t1.randomSU rnd
      t1 += 50
      t1.projectSU
    let p0 = g.plaq
    var gstop = 1e-7
    var orf = 1.8
    let tdirs = @[ @[0,1,2,3], @[0,1,2], @[3], @[1,2,3], @[0,2,3] ]
    let sf = 1.0/(lo.physVol.float*t[0].nrows.float^2)
    for dirs in tdirs:
      echo dirs
      threads: t := 1
      getGaugeFixTransform(t, g, dirs, gstop, orf)
      g1.gaugeTransform(g, t)
      let p1 = g1.plaq
      let l1 = g1.linkTrace dirs
      g2.gaugeTransform(g1, t1)
      threads: t2 := 1
      getGaugeFixTransform(t2, g2, dirs, gstop, orf)
      g.gaugeTransform(g2, t2)
      let p2 = g.plaq
      let l2 = g.linkTrace dirs
      CT = 1e-6
      check(p1 ~ p2)
      check(l1 ~ l2)
      let d = sqrt(sf * getDiff(g, g1))
      echo "RMS diff: ", d
      check(d <= 4e-2/dirs.len.float)

qexFinalize()
