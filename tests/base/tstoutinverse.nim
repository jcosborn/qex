import qex, gauge/stoutsmear
import os, sequtils

proc printPlaq(g: auto) =
  let
    p = g.plaq
    sp = 2.0*(p[0]+p[1]+p[2])
    tp = 2.0*(p[3]+p[4]+p[5])
  echo "plaq ",p
  echo "plaq ss: ",sp," st: ",tp," tot: ",p.sum

qexInit()

let
  fn = if paramCount() > 0: paramStr 1 else: ""
  lat = if fn.len == 0: @[8,8,8,8] else: fn.getFileLattice
  lo = lat.newLayout
var
  g = lo.newGauge
  f = lo.newGauge
  u = lo.newGauge
var ss = lo.newStoutSmear(0.02)    # large values could lead to diverging inverse iterations
if fn.len == 0:
  g.random
  for n in 0..<10:
    ss.smear(g, g)
elif 0 != g.loadGauge fn:
  echo "ERROR: couldn't load gauge file: ",fn
  qexFinalize()
  quit(-1)
g.printPlaq

ss.smear(g, f)
f.printPlaq
let (iter, r2) = ss.inverse(u, f)
qexLog "inverse iter ",iter," r2 ",r2
u.printPlaq

var del2 = 0.0
const nc = u[0][0].nrows
let vol = u[0].l.physVol
let nd = u.len
threads:
  var d2 = 0.0
  for mu in 0..<nd:
    for e in u[0]:
      let d = norm2(-1.0 + u[mu][e] * g[mu][e].adj)
      d2 += d.simdSum
  threadRankSum d2
  threadMaster:
    del2 = d2/float(2*(nc*nc+1)*nd*vol)

echo del2

if del2 > 1e-24:
  qexAbort()
else:
  qexFinalize()
