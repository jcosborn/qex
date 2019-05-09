import base
import physics/qcdTypes

proc relResid*(r,x: any, a: float): float =
  mixin norm2, simdReduce
  var t: type(r[0].norm2)
  for e in r:
    let r2 = r[e].norm2
    let x2 = x[e].norm2
    t += r2/(x2+a)
  result = simdReduce(t)
  threadRankSum(result)

proc relR*(p,r,x: any, a: float) =
  mixin norm2
  for e in p:
    #let w = asReal(1.0/(a+x[e].norm2))
    let w = asReal(sqrt(1.0/(a+x[e].norm2)))
    #let w = 1.0
    #let w = asReal(sqrt((a+x[e].norm2)))
    #let w = asReal(a+x[e].norm2)
    p[e] := w * r[e]

proc relResidUpdate*(x: any, Ax1,Ax2: any, b: any, a: float):
  array[2,DComplex] =
  var t: array[6,type(dot(Ax1[0],Ax1[0]))]
  for e in x:
    let w = asReal(1.0/(a+x[e].norm2))
    let wAx1 = w * Ax1[e]
    let wAx2 = w * Ax2[e]
    t[0] += dot(Ax1[e], wAx1)
    t[1] += dot(Ax1[e], wAx2)
    t[2] += dot(Ax2[e], wAx1)
    t[3] += dot(Ax2[e], wAx2)
    t[4] += dot(wAx1, b[e])
    t[5] += dot(wAx2, b[e])
  var t2: array[6,DComplex]
  for i in 0..5:
    t2[i].re = simdReduce(t[i].re)
    t2[i].im = simdReduce(t[i].im)
  threadRankSum(cast[ptr array[12,float]](t2[0].addr)[])
  let d = t2[0]*t2[3] - t2[2]*t2[1]
  let di = 1.0/d
  result[0] = di*(t2[4]*t2[3] - t2[5]*t2[1])
  result[1] = di*(t2[0]*t2[5] - t2[2]*t2[4])


when isMainModule:
  import qex
  import physics/qcdTypes
  qexInit()
  echo "rank ", myRank, "/", nRanks
  #var lat = [8,8,8,8]
  var lat = [4,4,4,4]
  var lo = newLayout(lat)
  var m = lo.ColorMatrix()
  var v1 = lo.ColorVector()
  var v2 = lo.ColorVector()
  var v3 = lo.ColorVector()
  var r = lo.ColorVector()
  var b = lo.ColorVector()
  type opArgs = object
    m: type(m)
  var oa = opArgs(m: m)
  proc apply*(oa: opArgs; r: type(v1); x: type(v1)) =
    r := oa.m*x
    #mul(r, m, x)
  threads:
    m.even := 1
    m.odd := 10
    threadBarrier()
    tfor i, 0..<lo.nSites:
      m{i} := i+1
    threadBarrier()
    b.even := 1
    b.odd := 2
    v1 := 0
    echo b.norm2
    echo m.norm2
  template resid(r,b,x,oa: untyped) =
    oa.apply(r, x)
    r := b - r

  var ra = 1e-10

  r.resid(b,v1,oa)
  echo "rsq: ", r.norm2
  echo "rel: ", relResid(r,v1,ra)
  for e in v1:
    v1[e] := asReal(1.0/m[e][0,0].re)*b[e]

  r.resid(b,v1,oa)
  echo "rsq: ", r.norm2
  echo "rel: ", relResid(r,v1,ra)

  oa.apply(v2, v1)
  v3 := 0.9 * v2
  var rru = relResidUpdate(v1, v3,v1, b, ra)
  echo rru

  var Av2 = lo.ColorVector()
  var Ap = lo.ColorVector()
  var p = lo.ColorVector()
  v2 := b
  oa.apply(Av2, v2)
  r := b - Av2
  var nits = intParam("nits", 20)
  for c in 1..nits:

    #p := r
    for e in p:
      p[e] := asReal(1.0/m[e][0,0].re)*r[e]
    p.relR(p, v2, ra)
    #oa.apply(p, r)
    #p.relR(r, v2, ra)

    oa.apply(Ap, p)
    var rru = relResidUpdate(v2, Av2,Ap, b, ra)
    #echo rru
    v2 := rru[0]*v2 + rru[1]*p
    Av2 := rru[0]*Av2 + rru[1]*Ap
    r := (1-rru[0])*b + rru[0]*r - rru[1]*Ap
    let r2 = r.norm2
    #v3.resid(b,v2,oa)
    #let tr2 = v3.norm2
    oa.apply(Av2, v2)
    r := b - Av2
    echo "r2: ", r2
    #echo "r2: ", r2, "   tr2: ", tr2
    echo "rel: ", relResid(r,v2,ra)
