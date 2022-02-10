import qex
import testutils

suite "Test ordered product":
  qexInit()
  echo "rank ", myRank, "/", nRanks
  threads: echo "thread ",threadNum," / ",numThreads
  var
    (l,g,_) = setupLattice([8,8,8,8])
    p = g.plaq

  test "plaquette":
    var pl = newseq[float]()
    for x in g.wilsonLines [@[1,2,-1,-2], @[1,3,-1,-3], @[2,3,-2,-3], @[1,4,-1,-4], @[2,4,-2,-4], @[3,4,-3,-4]]:
      pl.add x.re/6.0
    check(p~pl)

  test "gaugeProd":
    let ps = g.gaugeProd [
      @[1,2,-1,-2], @[2,1,-2,-1], @[1,-2,-1,2], @[2,-1,-2,1], @[-1,-2,1,2],
      @[1,2,-1,-2,1,-2,-1,2], @[-1,-2,1,2,-1,2,1,-2], @[-2,-1,2,1,1,-2,-1,2], @[-2,1,2,-1,-1,-2,1,2]]
    var t = newOneOf ps[0]
    var d2:float
    threads:
      t := ps[1].adj
    check(ps[0]~t)
    threads:
      t := ps[0]*ps[2]
    check(ps[5]~t)
    threads:
      t := ps[4]*ps[3].adj
    check(ps[6]~t)
    threads:
      t := ps[4].adj*ps[2]
    check(ps[7]~t)
    threads:
      t := ps[2].adj*ps[4]
    check(ps[8]~t)

  qexFinalize()
