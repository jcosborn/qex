import qex
import testutils
import sequtils

proc linkTrace(g: auto):auto =
  let n = g[0][0].ncols * g[0].l.physVol * g.len
  var lt: type(g[0].trace)
  threads:
    var t = g[0].trace
    for i in 1..<g.len: t += g[i].trace
    threadSingle: lt := t/n.float
  return lt

suite "Stress plaquette test":
  qexInit()
  const
    #nd = 4
    lat = [8,8,8,8]
  var
    lo = lat.newLayout
    g = lo.newGauge
    rs = newRNGField(RngMilc6, lo, 987654321)
    #rsX: RngMilc6  # workaround Nim codegen bug

  test "unit gauge":
    let
      l = g.linkTrace
      p = g.plaq
    const
      le = 1.0
      pe = mapit(@[1.0,1,1,1,1,1],it/6)
    check(l.re~le)
    check(l.im~0)
    check(p~pe)

  test "change single link":
    for s in 0..<lo.physVol:
      var cr,ci: float
      let ri = lo.rankIndex(s)
      let r = ri.rank
      if lo.myrank == r:
        let i = ri.index
        g[0]{i}.gaussian rs{i}
        #g[0]{i}.projectSU
        var t: float
        for a in 0..<g[0][0].ncols:
          t := g[0]{i}[a,a].re
          cr += t
        cr /= 3.0
        for a in 0..<g[0][0].ncols:
          t := g[0]{i}[a,a].im
          ci += t
        ci /= 3.0
        #echo "i: ",i
      rankSum(cr)
      rankSum(ci)
      let
        lr = 1 - (1-cr) / float(lo.physVol * g.len)
        li = ci / float(lo.physVol * g.len)
        l = g.linkTrace
        pr = 1 - 2*(1-cr) / float(lo.physVol)
        pe = mapit(@[pr,pr,1.0,pr,1.0,1.0],it/6)
        p = g.plaq
      check(l.re~lr)
      check(l.im~li)
      check(p~pe)
      if lo.myrank == r:
        let i = ri.index
        g[0]{i} := 1

  qexFinalize()
