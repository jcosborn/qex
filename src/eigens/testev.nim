import base
import eigens/svdLanczos
import eigens/linalgFuncs
import strUtils

template `&`(s: seq[float]): untyped = cast[ptr carray[float]](addr s[0])
template nothreads(x: untyped): untyped =
  block:
    x
template nothreads(v,x: untyped): untyped =
  block:
    x

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

proc rayleighRitz(op: any, v: any) =
  tic()
  let n = v.len
  let n2 = n*n
  var m1 = newSeq[type(dot(v[0].odd,v[0]))](n2)
  var m2 = newSeq[type(dot(v[0].even,v[0]))](n2)
  toc("rr setup")
  nothreads:
    for i in 0..<n:
      let n2 = v[i].even.norm2
      if threadNum==0: m2[i*n+i] := n2
      for j in countdown(i-1,0):
        let m2ji = dot(v[i].even,v[j])
        if threadNum==0: m2[j*n+i] = m2ji
  toc("rr m2")
  nothreads(t):
    for i in 0..<n:
      op.apply(v[i].odd, v[i].even)
    for i in 0..<n:
      for j in 0..i:
        let m1t = dot(v[i].odd,v[j])
        if threadNum==0: m1[j*n+i] = m1t
  toc("rr m1")
  for i in 0..<n:
    echo i, "  ", sqrt(m1[i*(n+1)].re/m2[i*(n+1)].re)

  var ev = newSeq[float64](n)
  var t0 = getElapsedTime()
  zeigsgv(cast[ptr float64](addr m1[0]), cast[ptr float64](addr m2[0]),
          addr ev[0], n)
  toc("rr zeigsgv")
  var t1 = getElapsedTime()
  #for i in 0..<n:
  #  echo i, "  ", sqrt(ev[i])
  #  for j in 0..<n:
  #    echo m[i*n+j]
  #echo sqrt(ev[0]), "\t", t[0].sv
  var r = newAlignedMem[type(v[0][0])](n)
  for e in v[0]:
    for i in 0..<n:
      mixin `:=`
      r[i] := 0
      for j in 0..<n:
        r[i] += m1[i*n+j] * v[j][e]
    for i in 0..<n:
      v[i][e] := r[i]
  #sort(t, a, b)
  for i in 0..<n:
    op.apply(v[i].odd, v[i].even)
    echo i, "  ", sqrt(v[i].odd.norm2/v[i].even.norm2)
  toc("rr end")
  var t2 = getElapsedTime()
  echo "rr dots: ", t0|-6, "  zeigsgv: ", (t1-t0)|-6, "  vecs: ", (t2-t1)|-6

proc cdot(m: any, c0,c1: int): float =
  let nr = m.nrows
  var s = 0.0
  for i in 0..<nr:
    s += m[i,c0]*m[i,c1]
  s

proc orthocols(m: any, c0: int) =
  let nr = m.nrows
  let nc = m.ncols
  var si = 1.0/sqrt(cdot(m, c0, c0))
  for i in 0..<nr:
    m[i,c0] *= si
  for c in (c0+1)..<nc:
    let di = cdot(m, c0, c)
    for i in 0..<nr:
      m[i,c] -= di * m[i,c0]
    let ei = 1.0/sqrt(cdot(m, c, c))
    for i in 0..<nr:
      m[i,c] *= ei

proc mulTriU(y: any, d: any, e: any, x: any) =
  let n1 = d.len - 1
  for i in 0..<n1:
    y[i,0] = d[i]*x[i,0] + e[i]*x[i+1,0]
  y[n1,0] = d[n1]*x[n1,0]

proc mulTriL(y: any, d: any, e: any, x: any) =
  y[0,0] = d[0]*x[0,0]
  for i in 1..<d.len:
    y[i,0] = d[i]*x[i,0] + e[i-1]*x[i-1,0]

proc getSv(d: any, e: any, x: any, y: any, r: any): auto =
  mulTriU(y, d, e, x)
  mulTriL(r, d, e, y)
  var x2,s2,r2 = 0.0
  for i in 0..<d.len:
    x2 += x[i,0]*x[i,0]
    s2 += y[i,0]*y[i,0]
    r2 += r[i,0]*r[i,0]
  let sv = sqrt(s2/x2)
  var er = 0.5*sqrt(r2/s2-s2/x2)
  (sv, er)

proc lowsv(linop: any, dest: any, src: any, maxit: int) =
  let nrr = dest.len
  var d = newSeq[float](maxit)
  var e = newSeq[float](maxit)
  var qv = newSeq[type(src.even)](nrr)
  var dv = newDmat(maxit, nrr)
  var qu = newSeq[type(src)](0)
  var du = newZmat(maxit, 0)
  getBidiagLanczos(linop, src, d, e, qu, du, qu, du, maxit, 1)
  e[maxit-1] = 0
  #echo "min d: ", d.min
  #echo "min e: ", e.min

  #var sv = newSeq[float](maxit)
  #svdbi(&sv, &d, &e, maxit)
  #echo "svdbi: ", sv[0], "  ", sv[1]
  #[
  var ur = newDmat(maxit, nrr)
  var nvout = svdBi4(sv, dv, ur, d, e, maxit, nrr, nrr, 0.0, 9e99)
  ]#

  var x = newDmat(maxit,1)
  var y = newDmat(maxit,1)
  var r = newDmat(maxit,1)
  for i in 0..<maxit:
    x[i,0] = 1.0/sqrt(maxit.float)
    #x[i,0] = 0.0
  #x[0,0] = 1.0

  for c in 0..<nrr:
    solveTriL(y, d, e, x)
    solveTriU(x, d, e, y)
    let sv = getSv(d, e, x, y, r)
    echo sv
    let k = nrr-1-c
    for i in 0..<maxit:
      dv[i,k] = x[i,0]
    if c>0:
      let d00 = cdot(dv, k, k)
      let d01 = cdot(dv, k, k+1)
      let d11 = cdot(dv, k+1, k+1)
      let c01 = d01/sqrt(d00*d11)
      echo c," c01: ", c01
    orthocols(dv, k)
  # ]#

  for i in 0..<dest.len:
    qv[i] = dest[i].even
    qv[i] := 0
  runBidiagLanczos(linop, src, d, e, dv, qv, du, qu, maxit, 1)

  rayleighRitz(linop, dest)
  rayleighRitz(linop, dest)

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
  echo sv[0], "  ", sv[1]

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
  import rng/milcrng

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
  proc getRay(op: any, v: any): auto =
    apply(op, r.odd, v.even)
    applyAdj(op, t.even, r.odd)
    let v2 = v.even.norm2
    let s2 = r.odd.norm2
    let sv = sqrt(s2/v2)
    t.even -= (s2/v2)*v
    (sv, 0.5*sqrt(t.even.norm2/s2))

  var nsrc = intParam("nsrc", 10)
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

  while true:
    var re, re0 = 0.0
    template pr(s: string, x: typed) =
      re = x[1]/x[0]
      echo s, x, "  ", re
    var s0 = getRay(op, srcs[0])
    pr "d: ", s0
    lowsv(op, dests, srcs[0], maxit)
    for i in 0..<dests.len:
      let si = getRay(op, dests[i])
      pr("sv[" & $i & "]: ", si)
      if i==0:
        s0 = si
        re0 = re
    if re0<1e-4: break
    srcs[0] := dests[0]

  #[
  #solve(op, dests, srcs, m, maxit)
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
  ]#

  qexFinalize()
