import base
import layout
import field
import gaugeUtils
import maths
import simd

proc gaugeTransform*(gt: any, g: any, t: any) =
  tic()
  let n = g.len
  var sf = newSeq[Shifter[type(t),type(t[0])]](n)
  for mu in 0..<n:
    sf[mu] = newShifter(t, mu, 1)
  toc("gt newShifters")
  threads:
    for mu in 0..<n:
      discard sf[mu] ^* t
    for mu in 0..<n:
      gt[mu] := t * (g[mu] * sf[mu].field.adj)
  toc("gt done")

proc gtGradient(grad: Field, g: any, t: any, dirs: array|seq) =
  mixin imadd
  tic()
  let n = dirs.len
  var sf = newSeq[shiftBType(t)](n)
  var sb = newSeq[shiftBType(t)](n)
  for i,mu in dirs.pairs:
    sf[i].initShiftB(t, mu,  1, "all")
    sb[i].initShiftB(t, mu, -1, "all")
  toc("grad2 initShiftB")
  threads:
    for i in 0..<n:
      startSB(sf[i], t[ix].adj)
    for i in 0..<n:
      let mu = dirs[i]
      startSB(sb[i], g[mu][ix].adj*t[ix].adj)
    for ir in grad:
      var m: type(grad[0])
      for i in 0..<n:
        let mu = dirs[i]
        localSB(sf[i], ir, imadd(m, g[mu][ir], it), t[ix].adj)
        localSB(sb[i], ir, iadd(m, it), g[mu][ix].adj*t[ix].adj)
      assign(grad[ir], m)
    for i in 0..<n:
      let mu = dirs[i]
      boundarySB(sf[i], imadd(grad[ir], g[mu][ir], it))
    for i in 0..<n:
      boundarySB(sb[i], iadd(grad[ir], it))
  #grad[ir][].projectTAH(m[])
  toc("grad2 done")

proc gtGradientX(grad: Field, gt: any, dirs: array|seq) =
  tic()
  let n = dirs.len
  var sb = newSeq[Shifter[type(gt[0]),type(gt[0][0])]](n)
  for i,mu in dirs.pairs:
    sb[i] = newShifter(gt[mu], mu, -1)
  toc("gradX newShifters")
  threads:
    for i,mu in dirs.pairs:
      discard sb[i] ^* gt[mu]
    for e in grad:
      var m = gt[dirs[0]][e] + sb[0].field[e].adj
      for k in 1..<n:
        m += gt[dirs[k]][e] + sb[k].field[e].adj
      grad[e] := m
  toc("gradX done")

proc gtGradient(grad: Field, gt: any, dirs: array|seq) =
  tic()
  let n = dirs.len
  var sb = newSeq[Shifter[type(gt[0]),type(gt[0][0])]](n)
  for i,mu in dirs.pairs:
    sb[i] = newShifter(gt[mu], mu, -1)
  toc("grad newShifters")
  threads:
    for i,mu in dirs.pairs:
      discard sb[i] ^* gt[mu]
    for e in grad:
      var m = gt[dirs[0]][e] + sb[0].field[e].adj
      for k in 1..<n:
        m += gt[dirs[k]][e] + sb[k].field[e].adj
      grad[e][].projectTAH(m[])
  toc("grad done")

proc gtUpdate(t: any, a: any, s: float) =
  tic()
  threads:
    for e in t:
      let m = t[e]
      t[e] := exp(-s*a[e]) * m
  toc("gtUpdate done")

proc gtUpdate(t: any, a: any, s: Field) =
  tic()
  threads:
    for e in t:
      let m = t[e]
      t[e] := exp(-s[e]*a[e]) * m
  toc("gtUpdateF done")

proc gtDistance(t0: any, t: any): float =
  let sf = 1.0/(t.l.physVol.float*t[0].nrows)
  let d = trace(t*t0.adj).re
  1.0 - sf*d

proc gfLocalMetric(m: Field, g: any, dirs: array|seq) =
  tic()
  let n = g.len
  var sb = newSeq[Shifter[type(g[0]),type(g[0][0])]](n)
  for mu in dirs:
    sb[mu] = newShifter(g[0], mu, -1)
  toc("LocalMetric newShifters")
  threads:
    for mu in dirs:
      discard sb[mu] ^* g[mu]
    for e in m:
      var t = trace(g[dirs[0]][e]).re + trace(sb[0].field[e]).re
      for k in 1..<dirs.len:
        t += trace(g[dirs[k]][e]).re + trace(sb[k].field[e]).re
      m[e] := t
  toc("LocalMetric done")

proc gfMetric(g: any, dirs: array|seq): float =
  tic()
  let sf = 1.0/(dirs.len.float*g[0].l.physVol.float*g[0][0].nrows)
  result = 0.0
  for d in dirs:
    result += trace(g[d]).re
  toc("metric done")
  result *= sf

proc gfMetrics(gd: any, t: any, nd: int): auto =
  tic()
  let sf = 0.5/(nd.float*gd.l.physVol.float*gd[0].nrows)
  let sfg = 2.0 * sf * nd.float
  var met, gre, gro = 0.0
  threads:
    var mt, ge, go: type(t[0][0,0].re)
    for e in gd.even:
      var tgd = t[e][]*gd[e][]
      mt += trace(tgd).re
      tgd.projectTAH
      ge += tgd.norm2
    for e in gd.odd:
      var tgd = t[e]*gd[e]
      mt += trace(tgd).re
      tgd.projectTAH
      go += tgd.norm2
    #var mt2 = simdSum(mt)
    #var ge2 = simdSum(ge)
    #var go2 = simdSum(go)
    #mt2.threadRankSum
    #ge2.threadRankSum
    #go2.threadRankSum
    var r = [simdSum(mt), simdSum(ge), simdSum(go)]
    r.threadRankSum
    met = r[0]
    gre = r[1]
    gro = r[2]
  result = (met: sf*met, gre: sfg*gre, gro: sfg*gro)
  toc("metric")

proc gfMetric(gd: any, t: any, nd: int): float =
  tic()
  let sf = 0.5/(nd.float*gd.l.physVol.float*gd[0].nrows)
  var met = 0.0
  threads:
    var mt: type(t[0][0,0].re)
    for e in gd:
      mt += trace(t[e]*gd[e]).re
    var mt2 = simdSum(mt)
    mt2.threadRankSum
    met = mt2
  result = sf * met
  toc("metric")

proc gfMetric(g: any, t: any, dirs: array|seq): float =
  tic()
  var gd = t.newOneOf
  gd.gtGradient(g, t, dirs)
  result = gfMetric(gd, t, dirs.len)
  toc("metric done")

proc gfLineMin(g: any, gd: any, t,t0: any,
               dirs: array|seq, eps: var float, m0: float) =
  tic()
  t0 := t
  threads:
    for e in gd:
      var tgd = t[e][]*gd[e][]
      gd[e][].projectTAH(tgd)
  echo "eps: ", 0.0, "  ", m0
  #t := t0
  t.gtUpdate(gd, eps)
  #t.projectSU
  var m1 = gfMetric(g, t, dirs)
  echo "eps: ", eps, "  ", m1
  t.gtUpdate(gd, eps)
  #t.projectSU
  var m2 = gfMetric(g, t, dirs)
  echo "eps: ", 2*eps, "  ", m2
  # m = m0 (x-e)(x-2e)/(2e^2) - m1 (x)(x-2e)/(e^2) + m2 (x)(x-e)/(2e^2)
  # m = [m0(x-e)(x-2e) - 2m1 x(x-2e) + m2 x(x-e) ]/(2e^2)
  # m' = [m0(2x-3e) - 4m1 (x-e) + m2 (2x-e) ]/(2e^2)
  var x = eps*(3*m0-4*m1+m2)/(2*m0-4*m1+2*m2)
  x = max(0,x)
  x = min(2*eps,x)
  eps = x
  t := t0
  t.gtUpdate(gd, eps)
  t.projectSU
  var m3 = gfMetric(g, t, dirs)
  echo "eps: ", eps, "  ", m3
  toc("gfLineMin")

proc moveFromZero(x: float, eps: float): float =
  result = x
  if abs(x)<eps:
    if x<0: result = -eps
    else: result = eps

proc moveFromZero(x: SimdD4, e: float): SimdD4 =
  result[0] = moveFromZero(x[0], e)
  result[1] = moveFromZero(x[1], e)
  result[2] = moveFromZero(x[2], e)
  result[3] = moveFromZero(x[3], e)

proc overRelaxSu2(r: var any, x: any, i,j: int, o: float) =
  mixin rsqrt
  var r0 =  x[i,i].re + x[j,j].re
  var r1 = -x[j,i].im - x[i,j].im
  var r2 =  x[j,i].re - x[i,j].re
  var r3 =  x[j,j].im - x[i,i].im
  #r0 = moveFromZero(r0, 1e-12)
  #var n1 = r1*r1 + r2*r2 + r3*r3
  #var n2 = r0*r0 + n1
  #var n = sqrt(n2)
  var n = sqrt(r0*r0+r1*r1+r2*r2+r3*r3)
  r0 += n*(1-o)/o
  r0 = moveFromZero(r0, 1e-12)
  #let r0n = r0/n
  #let t2 = 1 - r0n
  #let t2o = o*t2
  #let r0o = 1 - t2o
  #let o2 = o*o
  #let r0o = max(-1, 1 - o2 + o2*r0n)
  #r0 = n*r0o
  #let d = max(1e-24, 1-r0n*r0n)
  #let s = sqrt((1-r0o*r0o)/d)/(-n)
  #let r0o = r0n
  #let s = -n
  #r0 = r0o
  #r1 = s*r1
  #r2 = s*r2
  #r3 = s*r3
  #let u00 = newComplex(r0o, s*r3)
  #let u01 = newComplex(s*r2, s*r1)
  #r0 *= r0*r0
  #r1 *= r1*r1
  #r2 *= r2*r2
  #r3 *= r3*r3
  var nn = rsqrt(r0*r0 + r1*r1 + r2*r2 + r3*r3)
  let u00 = newComplex(nn*r0, nn*r3)
  let u01 = newComplex(nn*r2, nn*r1)
  for l in 0..2:
    let ti = u00*r[i,l] + u01*r[j,l]
    let tj = u00.adj*r[j,l] - u01.adj*r[i,l]
    r[i,l] := ti
    r[j,l] := tj

proc relaxE(t: any, gd: Field2, g: any, dirs: array|seq, orf: float) =
  tic()
  threads:
    for e in t.even:
      var m0 = t[e][] * gd[e][]
      overRelaxSu2(t[e][], m0, 0,1, orf)
      var m1 = t[e][] * gd[e][]
      overRelaxSu2(t[e][], m1, 1,2, orf)
      var m2 = t[e][] * gd[e][]
      overRelaxSu2(t[e][], m2, 0,2, orf)
      #t[e][].projectSU
  toc("relaxE")

proc relaxO(t: any, gd: Field2, g: any, dirs: array|seq, orf: float) =
  tic()
  threads:
    for e in t.odd:
      var m0 = t[e][] * gd[e][]
      overRelaxSu2(t[e][], m0, 0,1, orf)
      var m1 = t[e][] * gd[e][]
      overRelaxSu2(t[e][], m1, 1,2, orf)
      var m2 = t[e][] * gd[e][]
      overRelaxSu2(t[e][], m2, 0,2, orf)
      #t[e][].projectSU
  toc("relaxO")

proc getGaugeFixTransform(t: Field, g: any, dirs: seq[int],
                          gstop=1e-5, orf=1.8) =
  tic()
  var gd = t.newOneOf
  var t0 = t.newOneOf
  var eps = 0.1
  var polish = 0
  var its = 0
  while true:
    tic()
    gd.gtGradient(g, t, dirs)
    toc("grad")
    let mets = gfMetrics(gd, t, dirs.len)
    let met = mets.met
    let gdsq = mets.gre+mets.gro
    echo its, "  tr: ", met
    echo "  gradE: ", 2.0*mets.gre
    echo "  gradO: ", 2.0*mets.gro
    echo "  grad:  ", gdsq
    toc("metrics")
    if gdsq <= gstop:
      inc polish
    else:
      polish = 0
    if polish>10: break
    inc its
    var updateType = its mod 2
    if polish > 0:
      updateType = 2
    case updateType
    of 0:
      relaxE(t, gd, g, dirs, orf)
    of 1:
      relaxO(t, gd, g, dirs, orf)
    of 2:
      gfLineMin(g, gd, t, t0, dirs, eps, met)
    else:
      t.projectSU
    toc("loop end")
  toc("main2")

when isMainModule:
  import qex
  qexInit()
  let defaultLat = @[8,8,8,8]
  defaultSetup()
  var t = lo.ColorMatrix()
  var t0 = lo.ColorMatrix()
  var t1 = lo.ColorMatrix()
  var ml = lo.Real()
  var m1 = lo.Real()
  var m2 = lo.Real()
  var gt = lo.newGauge()
  var gd = lo.ColorMatrix()
  var gd2 = lo.ColorMatrix()
  if fn == "":
    g.random
    t := g[0]
    gt.unit
    g.gaugeTransform(gt, t)
  t := 1
  var eps = floatParam("eps", 0.1)
  var nh = intParam("nh", 100)
  var gstop = floatParam("gstop", 1e-6)
  var orf = floatParam("orf", 1.5)
  var outfn = stringParam("o", "")

  echo "gradient^2 stopping condition (gstop): ", gstop
  echo "overrelaxation factor (orf): ", orf
  if outfn == "":
    echo "not saving result, no output file specified (o)"
  else:
    echo "output file (o): ", outfn

  template disp(g: typed, dirs: typed) =
    echo g.gfMetric(dirs)

  template pdisp(g: typed, dirs: typed) =
    let p = g.plaq
    let sp = 2.0*(p[0]+p[1]+p[2])
    let tp = 2.0*(p[3]+p[4]+p[5])
    echo "plaqs: ", p
    echo sp
    echo tp
    disp(g, dirs)

  let sf = 1.0/(g[0].l.physVol.float*g[0][0].nrows.float)
  #var dirs = newSeq[int]()
  #for dr in 0..2:
  #  dirs.add dr
  let dirs = @[0,1,2]
  var its = 0
  var polish = 0
  tic()
  block:
    var gdsq = 2*gstop
    pdisp(g, dirs)
    gt.gaugeTransform(g, t)
    pdisp(gt, dirs)
    var met = gfMetric(gt, dirs)
    toc("main1")
    getGaugeFixTransform(t, g, dirs, gstop, orf)
    toc("main2")

    gt.gaugeTransform(g, t)
    pdisp(gt, dirs)

  toc("main3")
  echoTimers()
  qexFinalize()
