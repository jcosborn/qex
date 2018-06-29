import qex
import testutils
import sequtils

proc linkTrace(g: any):auto =
  let n = g[0][0].ncols * g[0].l.physVol * g.len
  var lt: type(g[0].trace)
  threads:
    var t = g[0].trace
    for i in 1..<g.len: t += g[i].trace
    threadSingle: lt := t/n.float
  return lt

const nd = 4
proc replicate(g2,g1:openarray[Field]) =
  let
    lo1 = g1[0].l
    lat1 = lo1.physGeom
    lo2 = g2[0].l
    lat2 = lo2.physGeom
  if lat2 != lat1.mapit(2*it):
    echo "Error: the output lattice is not exactly twice as the input."
    qexAbort()
  # indexing with {} load/store the whole vector, unsafe in threads
  for mu in 0..<nd:
    for j in lo2.sites:
      var cv:array[nd,cint]
      lo2.coord(cv,(lo2.myRank,j))
      # echo cv,"  <-"
      for k in 0..<nd:
        if cv[k] >= lat1[k]: cv[k] -= lat1[k].cint
      # echo "    ",cv
      let i = lo1.rankIndex(cv).index
      # echo "copy ", j," <- ",i
      g2[mu]{j} := g1[mu]{i}
      #for a in 0..2:
      #  for b in 0..2:
      #    g2[mu]{j}[a,b] := g1[mu]{i}[a,b]
      #    #g2[mu]{j}[a,b].re := g1[mu]{i}[a,b].re
      #    #g2[mu]{j}[a,b].im := g1[mu]{i}[a,b].im
      #    #echotyped:
      #    #  g2[mu]{j}[a,b] := g1[mu]{i}[a,b]
      #[
      block:
        var x,y:float
        x := g1[mu]{i}[0,0].re
        y := g2[mu]{j}[0,0].re
        var cv1,cv2:array[nd,cint]
        lo2.coord(cv2,(lo2.myRank,j))
        lo1.coord(cv1,(lo1.myRank,i))
        if x != 1:
          echo mu,":copy ",cv2," (",j,") = ",y," <- ",cv1," (",i,") = ",x
      ]#


suite "Multi-Layout test":
  qexInit()
  const
    lat1 = [4,4,4,4]
    lat2 = [8,8,8,8]
  var
    lo1 = lat1.newLayout
    g1 = lo1.newGauge
    lo2 = lat2.newLayout
    g2 = lo2.newGauge
    rs1 = newRNGField(RngMilc6, lo1, 987654321)
    rs2 = newRNGField(RngMilc6, lo2, 987654321)

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

  test "change single link":
    #block:
    #  let i = 128
    for i in lo1.sites:
    #for i in 0..128:
    #for i in {0,128}:
      g1[0]{i}.gaussian rs1{i}
      #g1[0]{i}.projectSU
      #g1[0]{i}[a,a].re := c
      #g1[0]{i}[a,a].im := c
      var t:float
      var cr1,ci1:float
      for a in 0..<g1[0][0].ncols:
        t := g1[0]{i}[a,a].re
        cr1 += t
      cr1 /= 3.0
      for a in 0..<g1[0][0].ncols:
        t := g1[0]{i}[a,a].im
        ci1 += t
      ci1 /= 3.0
      #echo "i: ",i
      #for j in lo2.sites:
      block:
        let j = 0
        #echo "j: ",j
        g2[0]{j}.gaussian rs2{i}
        #g2[0]{j}.projectSU
        #g2[0]{j}[a,a].re := c
        #g2[0]{j}[a,a].im := c
        var cr2,ci2:float
        for a in 0..<g2[0][0].ncols:
          t := g2[0]{j}[a,a].re
          cr2 += t
        cr2 /= 3.0
        for a in 0..<g2[0][0].ncols:
          t := g2[0]{j}[a,a].im
          ci2 += t
        ci2 /= 3.0
        let
          lr1 = 1 - (1-cr1) / float(lo1.physVol * g1.len)
          li1 = ci1 / float(lo1.physVol * g1.len)
          lr2 = 1 - (1-cr2) / float(lo2.physVol * g2.len)
          li2 = ci2 / float(lo2.physVol * g2.len)
          l1 = g1.linkTrace
          l2 = g2.linkTrace
          pr1 = 1 - 2*(1-cr1) / float(lo1.physVol)
          pr2 = 1 - 2*(1-cr2) / float(lo2.physVol)
          pe1 = mapit(@[pr1,pr1,1.0,pr1,1.0,1.0],it/6)
          pe2 = mapit(@[pr2,pr2,1.0,pr2,1.0,1.0],it/6)
          p1 = g1.plaq
          p2 = g2.plaq
        check(l1.re~lr1)
        check(l1.im~li1)
        check(l2.re~lr2)
        check(l2.im~li2)
        check(p1~pe1)
        check(p2~pe2)
        g2[0]{j} := 1
      g1[0]{i} := 1

  test "change single link, replicate and double size":
    #for i in lo1.sites:
    block:
      let i = 0
      #echo "i: ",i
      #for mu in 0..<nd:
      block:
        let mu = 1
        g1[mu]{i}.gaussian rs1{i}
        #g1[mu]{i}[a,a].re := c
        #g1[mu]{i}[a,a].im := c
        var t:float
        var cr,ci:float
        for a in 0..<g1[mu][0].ncols:
          t := g1[mu]{i}[a,a].re
          cr += t
        cr /= 3.0
        for a in 0..<g1[mu][0].ncols:
          t := g1[mu]{i}[a,a].im
          ci += t
        ci /= 3.0
        g2.replicate g1
        let
          lr = 1 - (1-cr) / float(lo1.physVol * g1.len)
          li = ci / float(lo1.physVol * g1.len)
          l1 = g1.linkTrace
          l2 = g2.linkTrace
          pr = 1 - 2*(1-cr) / float(lo1.physVol)
          p1 = g1.plaq
          p2 = g2.plaq
        var pe = @[1.0,1,1,1,1,1]
        if mu == 0:
          pe[0] = pr
          pe[1] = pr
          pe[3] = pr
        elif mu == 1:
          pe[0] = pr
          pe[2] = pr
          pe[4] = pr
        elif mu == 2:
          pe[1] = pr
          pe[2] = pr
          pe[5] = pr
        elif mu == 3:
          pe[3] = pr
          pe[4] = pr
          pe[5] = pr
        pe.applyit(it/6.0)
        check(l1.re~lr)
        check(l1.im~li)
        check(l2.re~lr)
        check(l2.im~li)
        check(p1~pe)
        check(p2~pe)
        g1[mu]{i} := 1
        #[
        for i in lo2.sites:
          var x:float
          x := g2[mu]{i}[0,0].re
          if x == 1: continue
          var cv:array[nd,cint]
          lo2.coord(cv,(lo2.myRank,i))
          echo cv," (",i,") = ",x
        ]#

  test "change two links, replicate and double size":
    let mu = 0
    proc chg(n:array[nd,int]) =
      let j = lo1.rankIndex(n).index
      #echo j," : ",n
      g1[mu]{j} := 100*(n[0]+10*(n[1]+10*(n[2]+10*n[3])))
    chg([0,0,0,0])
    chg([0,1,0,0])
    g2.replicate g1
    let
      l1 = g1.linkTrace
      l2 = g2.linkTrace
      p1 = g1.plaq
      p2 = g2.plaq
    check(l1.re~l2.re)
    check(l1.im~l2.im)
    check(p1~p2)
    #echo "P1d: ",p1.mapit((1.0-6.0*it)*float(lo1.physVol))
    #echo "P2d: ",p2.mapit((1.0-6.0*it)*float(lo2.physVol))
    # look at each link in g2 and print out the location of non unit links
    #[
    for i in lo2.sites:
      var x:float
      x := g2[mu]{i}[0,0].re
      if x == 1: continue
      var cv:array[nd,cint]
      lo2.coord(cv,(lo2.myRank,i))
      echo cv," (",i,") = ",x
    ]#

  qexFinalize()
