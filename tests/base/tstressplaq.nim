import qex
import unittest, sequtils

const CT = 1e-14
proc `~`(x,y:float):bool =
  if x==0 or y==0: result = x==y
  else: result = abs(x-y)/max(abs(x),abs(y)) < CT
  if not result:
    echo "LHS: ",x
    echo "RHS: ",y
proc `~`[T](x,y:openarray[T]):bool =
  result = x.len == y.len
  for i in 0..<x.len:
    result = result and (x[i] ~ y[i])
  if not result:
    proc show[T](x:openArray[T]):string =
      result = $x.len & " #[ "
      for c in x:
        result &= $c & " "
      result &= "]"
    echo "LHS: ",x.show
    echo "RHS: ",y.show
proc linkTrace(g: any):auto =
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
    nd = 4
    lat = [4,4,4,4]
  var
    lo = lat.newLayout
    g = lo.newGauge
    rs: RngMilc6
  rs.seed(7,11)

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
    for i in lo.sites:
      g[0]{i}.gaussian rs
      #g[0]{i}.projectSU
      var t:float
      var cr,ci:float
      for a in 0..<g[0][0].ncols:
        t := g[0]{i}[a,a].re
        cr += t
      cr /= 3.0
      for a in 0..<g[0][0].ncols:
        t := g[0]{i}[a,a].im
        ci += t
      ci /= 3.0
      #echo "i: ",i
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
      g[0]{i} := 1

  qexFinalize()
