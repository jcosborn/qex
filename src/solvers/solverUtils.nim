import base
import physics/qcdTypes

# av = a V xp2 / x2
# Sum_s (1+av)rs2/(1+av*xs2/x2)
# a=0 -> r2
# a->infty -> x2 Sum_s rs2/xs2
proc relResid*(r,x: auto, a: float): float =
  mixin norm2, simdReduce
  var t: type(r[0].norm2)
  #var u: type(r[0].norm2)
  let av = a * x.l.physVol.float
  let ax = av / (x.norm2 + float.epsilon)
  for e in r:
    let r2 = r[e].norm2
    let x2 = x[e].norm2
    let s = (1.0+av)/(1.0+ax*x2)
    t += s*r2
    #u += s
  var ts = simdSum(t)
  #var us = simdSum(u)
  #var rs = [ts,us]
  #threadRankSum(rs)
  #result = x.l.physVol.float * rs[0]/rs[1]
  #result = rs[0]
  threadRankSum(ts)
  result = ts

proc relR*(p,r,x: auto, a: float) =
  mixin norm2
  for e in p:
    #let w = asReal(1.0/(a+x[e].norm2))
    let w = asReal(sqrt(1.0/(a+x[e].norm2)))
    #let w = 1.0
    #let w = asReal(sqrt((a+x[e].norm2)))
    #let w = asReal(a+x[e].norm2)
    p[e] := w * r[e]

proc relResidUpdate*(x: auto, Ax1,Ax2: auto, b: auto, a: float):
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
    t2[i].re = simdSum(t[i].re)
    t2[i].im = simdSum(t[i].im)
  threadRankSum(cast[ptr array[12,float]](t2[0].addr)[])
  let d = t2[0]*t2[3] - t2[2]*t2[1]
  let di = 1.0/d
  result[0] = di*(t2[4]*t2[3] - t2[5]*t2[1])
  result[1] = di*(t2[0]*t2[5] - t2[2]*t2[4])


#[  old test
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
]#

when isMainModule:
  import qex, physics/stagSolve, observables/sources
  qexInit()
  var defaultLat = @[8,8,8,8]
  #var defaultLat = @[12,12,12,12]
  defaultSetup()
  var v1 = lo.ColorVector()
  var v2 = lo.ColorVector()
  var v3 = lo.ColorVector()
  var r = lo.ColorVector()
  var rs = newRNGField(RngMilc6, lo, intParam("seed", 987654321).uint64)
  threads:
    g.random rs
    g.setBC
    g.stagPhase
    v1.gaussian rs
  echo v1.norm2
  var s = newStag(g)
  var m = floatParam("m", 0.01)
  var sp = newSolverParams()
  sp.verbosity = intParam("verb", 2)
  sp.subset.layoutSubset(lo, "all")
  sp.maxits = int(1e9/lo.physVol.float)
  sp.r2req = floatParam("rsq", 1e-12)

  proc test =
    v2 := 0
    s.solve(v2, v1, m, sp)
    threads:
      s.D(v3, v2, m)
      v1 := 0
    resetTimers()
    s.solve(v1, v3, m, sp)
    echo "x2:      ", norm2slice(v1, 3)
    echo "x2 even: ", norm2slice(v1.even, 3)
    echo "x2 odd:  ", norm2slice(v1.odd, 3)
    threads:
      r := v1 - v2
      let e2 = r.norm2
      s.D(r, v1, m)
      threadBarrier()
      r := v3 - r
      let r2 = r.norm2
      echo "r2:    ", r2
      let rr0 = relResid(r, v1, 0)
      echo "rr(0): ", rr0
      var a = 1
      for i in 0..16:
        let rra = relResid(r, v1, float a)
        echo "rr(10^",i,"): ", rra
        a *= 10
      echo "err2:  ", e2
    echo "r2:      ", norm2slice(r, 3)
    echo "r2 even: ", norm2slice(r.even, 3)
    echo "r2 odd:  ", norm2slice(r.odd, 3)
    s.Ddag(v3, r, m)
    r := v3
    echo "Dr2:      ", norm2slice(r, 3)
    echo "Dr2 even: ", norm2slice(r.even, 3)
    echo "Dr2 odd:  ", norm2slice(r.odd, 3)
    r := v1 - v2
    echo "e2:      ", norm2slice(r, 3)
    echo "e2 even: ", norm2slice(r.even, 3)
    echo "e2 odd:  ", norm2slice(r.odd, 3)


  block:
    v1 := 0
    let p = lo.rankIndex([0,0,0,0])
    if myRank==p.rank:
      v1{p.index}[0] := 1
    echo "even point"
    test()
    echo sp.getStats()

  qexFinalize()
