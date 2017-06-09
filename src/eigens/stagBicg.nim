import base
import eigens/svdLanczos
import eigens/linalgFuncs
import strUtils

proc makeTri(t0,t1: var any; d,e: any; m: float) =
  let n = t0.len
  let m2 = m*m
  t0[0] = d[0]*d[0] + m2
  for i in 0..(n-2):
    t1[i] = d[i]*e[i]
    t0[i+1] = e[i]*e[i] + d[i+1]*d[i+1] + m2

proc cholTri(c0,c1: var any; t0,t1: any) =
  let n = t0.len
  var x = t0[0]
  var xs = sqrt(x)
  c0[0] = xs
  for i in 0..(n-2):
    var l = t1[i]/xs
    c1[i] = l
    x = t0[i+1] - l*l
    xs = sqrt(x)
    c0[i+1] = xs

proc solveTriU(x: any; c0,c1: any; b: any) =
  let n = x.nrows
  let nc = x.ncols
  for j in 0..<nc:
    x[n-1,j] = b[n-1,j]/c0[n-1]
  for i in countdown(n-2, 0):
    for j in 0..<nc:
      x[i,j] = (b[i,j]-c1[i]*x[i+1,j])/c0[i]

proc solveTriL(x: any; c0,c1: any; b: any) =
  let n = x.nrows
  let nc = x.ncols
  for j in 0..<nc:
    x[0,j] = b[0,j]/c0[0]
  for i in countup(1, n-1):
    for j in 0..<nc:
      x[i,j] = (b[i,j]-c1[i-1]*x[i-1,j])/c0[i]


# A Vn = Un B
# A' Un = Vn1 B'
# (A'A +s) x = b
# (A'A +s) Vn y = b
# Vn' (A'A+s) Vn y = Vn' b
# (B'B+s) y = Vn' b

proc solve(linop: any, dests: any, srcs: any, m: float, maxit: int) =
  var nsrc = srcs.len
  var src = srcs[0].newOneOf
  var d = newSeq[float](maxit)
  var e = newSeq[float](maxit)
  var qv = newSeq[type(src)](nsrc)
  var dv = newZmat(maxit, nsrc)
  var qu = newSeq[type(src)](0)
  var du = newZmat(maxit, 0)
  template `&`(s: seq[float]): untyped = cast[ptr carray[float]](addr s[0])

  src := 0
  for i in 0..<nsrc:
    qv[i] = srcs[i].newOneOf
    qv[i] := srcs[i]
    src += srcs[i]
  #src := srcs[nsrc-1]
  getBidiagLanczos(linop, src.even, d, e, qv, dv, qu, du, maxit, 1)
  #for i in 0..<maxit:
  #  echo "dv[$#]: $#"%[$i,$dv[i,0]]

  var sv = newSeq[float](maxit)
  svdbi(&sv, &d, &e, maxit)
  echo sv[0]

  var vd = newDmat(maxit, nsrc)
  var ud = newDmat(maxit, nsrc)
  var nvout = svdBi4(sv, vd, ud, d, e, maxit, nsrc, nsrc, 0.0, 9e99)

  for i in 0..<nsrc:
    dests[i] := 0
    #qv[i] = dests[i]
  runBidiagLanczos(linop, src.even, d, e, vd, dests, du, qu, maxit, 1)
  src := 0
  for i in 0..<1:
    src += dests[0]

  #src := 0
  #srcs[0] := src
  for i in 0..<nsrc:
    qv[i] := srcs[i]
    src += srcs[i]
  srcs[0] := src
  qv[0] := src
  #src := srcs[0]
  getBidiagLanczos(linop, src.even, d, e, qv, dv, qu, du, maxit, 1)
  for i in 0..<nsrc:
    echo "berr[", i, "]: ", sqrt(qv[i].norm2/srcs[i].norm2)

  svdbi(&sv, &d, &e, maxit)
  echo sv[0]

  var t0 = newSeq[float](maxit)
  var t1 = newSeq[float](maxit)
  makeTri(t0, t1, d, e, m)
  var c0 = newSeq[float](maxit)
  var c1 = newSeq[float](maxit)
  cholTri(c0, c1, t0, t1)

  var s = 0.0
  for i in 0..<maxit:
    let df = d[i] - c0[i]
    s += df*df
  for i in 0..(maxit-2):
    let df = e[i] - c1[i]
    s += df*df
  echo "err: ", s

  var dv0 = newZmat(maxit, nsrc)
  var dv1 = newZmat(maxit, nsrc)
  solveTriL(dv0, c0, c1, dv)
  solveTriU(dv1, c0, c1, dv0)

  for i in 0..<nsrc:
    dests[i] := 0
    qv[i] = dests[i]
  runBidiagLanczos(linop, src.even, d, e, dv1, qv, du, qu, maxit, 1)

when isMainModule:
  import qex
  import physics/qcdTypes
  import gauge
  import physics/stagD
  import physics/hisqLinks
  import rng/rng

  qexInit()
  var defaultGaugeFile = "l88.scidac"
  #var defaultLat = [8,8,8,8]
  #var defaultLat = [16,16,16,16]
  defaultSetup()
  threads:
    #g.random
    g.setBC
    g.stagPhase
  #var s = newStag(g)
  var hc: HisqCoefs
  hc.init()
  echo hc
  var fl = lo.newGauge()
  var ll = lo.newGauge()
  hc.smear(g, fl, ll)
  var s = newStag3(fl, ll)

  var lo1 = newLayout(lat, 1)
  var rs: Field[1,RngMilc6]
  rs.new(lo1)
  var seed = 987654321
  threads:
    for j in lo.sites:
      var l = lo.coords[lo.nDim-1][j].int
      for i in countdown(lo.nDim-2, 0):
        l = l * lo.physGeom[i].int + lo.coords[i][j].int
      rs[j].seed(seed, l)

  type MyOp = object
    s: type(s)
    r: type(rs)
    lo: type(lo)
  var op = MyOp(r:rs,s:s,lo:lo)
  proc gaussian(x: AsVarComplex, r: var any) =
    x.re = gaussian(r)
    x.im = gaussian(r)
  proc gaussian(x: AsVarVector, r: var any) =
    for i in 0..<x.len:
      gaussian(x[i], r)
  proc gaussian(v: Field, r: Field2) =
    for i in v.l.sites:
      gaussian(v{i}, r[i])
  template rand(op: MyOp, v: any) =
    gaussian(v, op.r)
  template newVector(op: MyOp): untyped =
    op.lo.ColorVectorD()
  template apply(op: MyOp, r,v: typed) =
    threadBarrier()
    stagD(op.s.so, r.field, op.s.g, v.field, 0.0)
  template applyAdj(op: MyOp, r,v: typed) =
    threadBarrier()
    stagD(op.s.se, r.field, op.s.g, v.field, 0.0, -1)
  template newRightVec(op: MyOp): untyped = newVector(op).even
  template newLeftVec(op: MyOp): untyped = newVector(op).odd

  var m = floatParam("mass", 0.0)
  var maxit = intParam("maxit", 100)
  var r = lo.ColorVectorD()
  var t = lo.ColorVectorD()
  var m2 = m*m
  proc applyD2(op: any, dest: any, src: any, m: float) =
    apply(op, t.odd, src.even)
    applyAdj(op, dest.even, t.odd)
    let m2 = m*m
    dest.even += m2*src
  proc getResid(op: any, r: any, d,s: any) =
    #echo "r.odd: ", r.odd.norm2
    #echo "d.even: ", d.even.norm2
    apply(op, r.odd, d.even)
    #echo "r.odd: ", r.odd.norm2
    applyAdj(op, t.even, r.odd)
    r.even := s - t - m2*d

  var nsrc = intParam("nsrc", 1)
  var srcs = newSeq[type(r)](nsrc)
  var dests = newSeq[type(r)](nsrc)
  for i in 0..<nsrc:
    srcs[i] = lo.ColorVectorD()
    dests[i] = lo.ColorVectorD()
  for i in 0..<nsrc:
    srcs[i] := 0
    dests[i] := 0
    let c = i mod 3
    let s = i div 3
    srcs[i]{s}[c] := 1

  solve(op, dests, srcs, m, maxit)

  for i in 0..<nsrc:
    getResid(op, r, dests[i], srcs[i])
    #echo "src2[", i, "]: ", srcs[i].even.norm2
    #echo "dest2[", i, "]: ", dests[i].even.norm2
    echo "r2[", i, "]/src2: ", r.even.norm2/srcs[i].even.norm2

  var sp = initSolverParams()
  sp.maxits = 10_000
  sp.r2req = 1e-16
  for i in 0..<nsrc:
    srcs[i].odd := 0
    dests[i] := 0
    s.solve(dests[i], srcs[i], m, sp)
    srcs[i].odd := 0
    dests[i].odd := 0
    getResid(op, r, dests[i], srcs[i])
    echo "src2[", i, "]: ", srcs[i].even.norm2
    echo "dest2[", i, "]: ", dests[i].even.norm2
    echo "r2[", i, "]: ", r.even.norm2

  qexFinalize()
