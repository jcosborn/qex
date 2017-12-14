#[
  Replicate a 4D gauge configuration 16 times and stitch them together
  to get a double sized lattice.  It serves as an example of handling two
  lattices of different sizes.  It also shows indexing routines of the on-memory
  order and the lexicographical order, and the transformation between those.

  This program is designed for use with single rank to avoid dealing with
  data movement.
]#

import qex
import seqUtils

proc printPlaq(g: any) =
  let
    p = g.plaq
    sp = 2.0*(p[0]+p[1]+p[2])
    tp = 2.0*(p[3]+p[4]+p[5])
  echo "plaq ",p
  echo "plaq ss: ",sp," \tst: ",tp," \ttot: ",p.sum

proc linkTrace(g: any): auto =
  let n = g[0][0].ncols * g[0].l.physVol * g.len
  var lt: type(g[0].trace)
  threads:
    var t = g[0].trace
    for i in 1..<g.len: t += g[i].trace
    threadSingle: lt := t / n.float
  return lt

qexInit()
if nRanks > 1:
  echo "Error: `replicate' only designed for running with single rank."
  qexExit 1
threads: echo "thread ", threadNum, "/", numThreads

const nd = 4
var
  (lo, g, r) = setupLattice([4,4,4,4])
  lat = lo.physGeom
  lat2 = lat.mapIt(2*it)
  lo2 = lat2.newLayout
  g2 = lo2.newGauge

threads:
  let d = g.checkSU
  echo "unitary deviation avg: ",d.avg," max: ",d.max

threads: g.projectSU

threads:
  let d = g.checkSU
  echo "new unitary deviation avg: ",d.avg," max: ",d.max
g.printPlaq
echo "LinkTrace: ",g.linkTrace

tic()

threads:
  for mu in 0..<nd:
    for j in lo2.sites:
      var cv:array[nd,cint]
      lo2.coord(cv,(lo2.myRank,j))
      # echo cv,"  <-"
      for k in 0..<nd:
        if cv[k] >= lat[k]: cv[k] -= lat[k].cint
      # echo "    ",cv
      let i = lo.rankIndex(cv).index
      # echo "copy ", j," <- ",i
      forO a, 0, 2:
        forO b, 0, 2:
          g2[mu]{j}[a,b].re := g[mu]{i}[a,b].re
          g2[mu]{j}[a,b].im := g[mu]{i}[a,b].im

toc("copy lattice")

# Check the plaquette
echo "Generated:"
g2.printPlaq
echo "LinkTrace: ",g2.linkTrace

let outfn = strParam("out","outlat.lime")
if 0 != g2.saveGauge outfn:
  echo "Error: failed to save gauge to ",outfn
  qexExit 1

qexFinalize()
