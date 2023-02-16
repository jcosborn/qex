import testutils
import qex
import sequtils, strformat

AT = 1e-12

proc linkTrace(g: auto):auto =
  let n = g[0][0].ncols * g[0].l.physVol * g.len
  var lt: evalType(g[0].trace)
  threads:
    var t = g[0].trace
    for i in 1..<g.len: t += g[i].trace
    threadSingle: lt := t/n.float
  return lt

const nd = 4
proc replicate(g2,g1: openarray[Field]) =
  let
    lov1 = g1[0].l
    lo1 = newLayout(lov1.physGeom, 1)
    lov2 = g2[0].l
    lo2 = newLayout(lov2.physGeom, 1)
    cm1 = lo1.ColorMatrix1()
    cm2 = lo2.ColorMatrix1()
  for mu in 0..<nd:
    #echo fmt"g1[{mu}]: {g1[mu].norm2}"
    #for i in lo1.sites:
    #  cm1[i] := g1[mu]{i}
    cm1.remapLocalFrom g1[mu]
    #echo fmt"{mu} cm1: {cm1.norm2}"
    cm2.replicateFrom(cm1)
    #echo fmt"{mu} cm2: {cm2.norm2}"
    #for i in lo2.sites:
    #  g2[mu]{i} := cm2[i]
    g2[mu].remapLocalFrom cm2
    #echo fmt"g2[{mu}]: {g2[mu].norm2}"
  #echo "g1[0][[0,0,0,0]]: ", g1[0][[0,0,0,0]]
  #echo "g2[0][[0,0,0,0]]: ", g2[0][[0,0,0,0]]

qexInit()

suite "Multi-Layout test":
  let
    lat1 = latticeFromLocalLattice([8,8,8,8], nRanks)
    lat2 = latticeFromLocalLattice([16,16,16,16], nRanks)
  var
    lo1 = lat1.newLayout
    g1 = lo1.newGauge
    lo2 = lat2.newLayout
    g2 = lo2.newGauge
    rs1 = newRNGField(RngMilc6, lo1, 987654321)
    #rs2 = newRNGField(RngMilc6, lo2, 987654321)
    #rsX: RngMilc6  # workaround Nim codegen bug
  if g1[0].numberType is float32: AT = 1e-8

  test "unit gauge":
    let
      l1 = g1.linkTrace
      l2 = g2.linkTrace
      p1 = g1.plaq
      p2 = g2.plaq
    const
      le = 1.0
      pe = mapit(@[1.0,1,1,1,1,1],it/6)
    check(l1.re~le)
    check(l1.im~0)
    check(l2.re~le)
    check(l2.im~0)
    check(p1~pe)
    check(p2~pe)

  test "random gauge, replicate and double size":
    g1.random rs1
    #echo fmt"g1[0]: {g1[0].norm2}"
    g2.replicate g1
    let
      l1 = g1.linkTrace
      l2 = g2.linkTrace
      p1 = g1.plaq
      p2 = g2.plaq
    withCT(1e-12):
      check(l1.re~l2.re)
      check(l1.im~l2.im)
      check(p1~p2)

qexFinalize()
