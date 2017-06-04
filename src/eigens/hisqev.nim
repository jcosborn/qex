import base
import layout
import field
import svdLanczos
import linalgFuncs
import math
import strUtils

type EigTable*[T] = object
  v*: T
  vn0*: float
  vn*: float
  Dvn*: float
  sv*: float
  err*: float
  updates*: int

proc sort[T](t: var seq[EigTable[T]], a=0, bb= -1) =
  var b = bb
  if b<0: b += t.len
  for i in a..b:
    var jmin = i
    for j in (i+1)..b:
      if t[j].sv<t[jmin].sv:
        jmin = j
    if jmin!=i: swap(t[i], t[jmin])

proc sett[T](t: var EigTable[T], op: any) =
  mixin `:=`
  t.vn = t.v.even.norm2
  if t.vn==0:
    #echo "sett: norm is zero"
    #quit 1
    t.vn0 = 0
    t.sv = 1e99
    return
  t.v.even *= 1/sqrt(t.vn)
  t.vn0 = t.vn
  t.vn = 1
  #let Dv = op.newVector
  let Dv = t.v
  op.apply(Dv.odd, t.v.even)
  t.Dvn = Dv.odd.norm2
  t.sv = sqrt(t.Dvn/t.vn)
  #echo t.sv

proc maketable[T,U](t: var seq[EigTable[T]], v: var seq[U], op: any) =
  for i in 0..<v.len:
    let tv = t[i].v
    t[i].v = v[i]
    v[i] = tv
    sett(t[i], op)
    #echo i, "  ", t[i].sv
  sort(t)

proc setterr(t: var EigTable, op: any) =
  mixin `:=`
  let Dv = op.newVector
  op.apply(Dv.odd, t.v.even)
  let DDv = op.newVector
  op.applyAdj(DDv.even, Dv.odd)
  let e = t.Dvn/t.vn
  Dv.even := DDv - e*t.v
  let r2 = Dv.even.norm2/t.vn
  t.err = sqrt(t.sv*t.sv+sqrt(r2)) - t.sv

proc geterr(v: var seq[EigTable], op: any) =
  for i in 0..<v.len:
    setterr(v[i], op)

# x = x - (<y.x>/<y.y>) y
proc projectOut(x,y: EigTable) =
  let c = y.v.even.dot(x.v)/y.vn
  x.v.even -= c*y.v

proc ortho1(t: var any; imn,imx: int, op: any) =
  for j in 0..<imn:
    for i in imn..imx:
      t[i].projectOut(t[j])
  for i in imn..imx:
    sett(t[i], op)
    if t[i].vn0<1e-10:
      #echo "ortho1 i: ", i, "  vn0: ", t[i].vn0
      if t[i].vn0==0: op.rand(t[i].v)
      for j in 0..<imn:
        t[i].projectOut(t[j])
      sett(t[i], op)

proc ortho2(t: var any; imn,imx: int, op: any) =
  for i in imn..imx:
    var i0 = i
    for j in (i+1)..imx:
      if t[j].sv < t[i0].sv: i0 = j
    if i0!=i: swap(t[i],t[i0])
    for j in (i+1)..imx:
      t[j].projectOut(t[i])
      sett(t[j], op)
      if t[j].vn0<1e-10:
        echo "j: ", j, "  vn0: ", t[j].vn0
        if t[j].vn0==0: op.rand(t[j].v)
        sett(t[j], op)
        for k in 0..<j:
          t[j].projectOut(t[k])
        sett(t[j], op)

# proc orthonormalize(v, vg, "even")
proc rayleighRitz2(t: var any, a,b: int) =
  mixin dot
  let n = b-a+1
  let n2 = n*n
  var m = newSeq[type(dot(t[a].v.odd,t[b].v))](n2)
  for i in 0..<n:
    for j in i..<n:
      m[i*n+j] = dot(t[a+i].v.odd,t[a+j].v)
  var ev = newSeq[float64](n)
  zeigs(cast[ptr float64](addr m[0]), addr ev[0], n)
  #for i in 0..<n:
  #  echo ev[i]
  #  for j in 0..<n:
  #    echo m[i*n+j]
  #echo sqrt(ev[0]), "\t", t[0].sv
  var r = newAlignedMem[type(t[a].v[0])](n)
  for e in t[a].v:
    for i in 0..<n:
      mixin `:=`
      r[i] := 0
      for j in 0..<n:
        r[i] += m[i*n+j].adj * t[a+j].v[e]
    for i in 0..<n:
      t[a+i].v[e] := r[i]
  for i in a..<n:
    t[a+i].Dvn = ev[i]
    #echo t[i].vn
    t[a+i].sv = sqrt(ev[i])

proc rayleighRitz(t: var any, a,b: int, op: any) =
  mixin adj
  tic()
  mixin dot, simdSum, rankSum, `:=`
  let n = b-a+1
  let n2 = n*n
  var m1 = newSeq[type(dot(t[a].v.odd,t[b].v))](n2)
  var m2 = newSeq[type(dot(t[a].v.even,t[b].v))](n2)
  #for i in 0..<n:
  #  for j in i..<n:
  #    m1[i*n+j] = dot(t[a+i].v.odd,t[a+j].v)
  #    m2[i*n+j] = dot(t[a+i].v.even,t[a+j].v)
  #for i in 0..<n:
  #  for j in 0..i:
  #    m1[j*n+i] = dot(t[a+i].v.odd,t[a+j].v)
  #    m2[j*n+i] = dot(t[a+i].v.even,t[a+j].v)
  toc("rr setup")
  threads(t):
    for i in 0..<n:
      let n2 = t[a+i].v.even.norm2
      if threadNum==0: m2[i*n+i] := n2
      for j in countdown(i-1,0):
        let m2ji = dot(t[a+i].v.even,t[a+j].v)
        if threadNum==0: m2[j*n+i] = m2ji
  toc("rr m2")
  threads(t):
    for i in 0..<n:
      for j in 0..i:
        let m1t = dot(t[a+i].v.odd,t[a+j].v)
        if threadNum==0: m1[j*n+i] = m1t
  toc("rr m1")

  #[
  let nu = n*(n+1) div 2
  var dte = newAlignedMem[type(dot(t[a].v[0],t[b].v[0]))](nu)
  for e in t[a].v.even:
    for i in 0..<n:
      let k = i*n - ((i*(i-1)) div 2) - i
      for j in i..<n:
        dte[k+j] += dot(t[a+i].v[e],t[a+j].v[e])
  var dtes = newSeq[type(simdSum(dte[0]))](nu)
  for k in 0..<nu:
    dtes[k] = simdSum(dte[k])
  rankSum(dtes)
  var dto = newAlignedMem[type(dot(t[a].v[0],t[b].v[0]))](nu)
  for e in t[a].v.odd:
    for i in 0..<n:
      let k = i*n - ((i*(i-1)) div 2) - i
      for j in i..<n:
        dto[k+j] += dot(t[a+i].v[e],t[a+j].v[e])
  var dtos = newSeq[type(simdSum(dto[0]))](nu)
  for k in 0..<nu:
    dtos[k] = simdSum(dto[k])
  rankSum(dtos)
  for i in 0..<n:
    let k = i*n - ((i*(i-1)) div 2) - i
    for j in i..<n:
      m1[i*n+j] = dtos[k+j]
      m2[i*n+j] = dtes[k+j]
  ]#

  #toc("rr mat")
  var ev = newSeq[float64](n)
  var t0 = getElapsedTime()
  zeigsgv(cast[ptr float64](addr m1[0]), cast[ptr float64](addr m2[0]),
          addr ev[0], n)
  toc("rr zeigsgv")
  var t1 = getElapsedTime()
  #for i in 0..<n:
  #  echo i, "  ", ev[i]
  #  for j in 0..<n:
  #    echo m[i*n+j]
  #echo sqrt(ev[0]), "\t", t[0].sv
  var r = newAlignedMem[type(t[a].v[0])](n)
  for e in t[a].v:
    for i in 0..<n:
      mixin `:=`
      r[i] := 0
      for j in 0..<n:
        r[i] += m1[i*n+j] * t[a+j].v[e]
    for i in 0..<n:
      t[a+i].v[e] := r[i]
  for i in 0..<n:
    #t[a+i].Dvn = ev[i]
    #echo t[i].vn
    #t[a+i].sv = sqrt(ev[i])
    sett(t[a+i], op)
  sort(t, a, b)
  toc("rr end")
  var t2 = getElapsedTime()
  echo "rr dots: ", t0|-6, "  zeigsgv: ", (t1-t0)|-6, "  vecs: ", (t2-t1)|-6

proc sortdown(t: var seq[EigTable], n0,n: int): int =
  var imin = n
  for i in countdown(n,n0+1):
    if t[i].sv < t[i-1].sv:
      #echo "swap: ", i, ": ", t[i].sv, "  ", t[i-1].sv
      swap(t[i], t[i-1])
      imin = min(imin, i-1)
  imin
template sortdown(t: var seq[EigTable], n: int): int = sortdown(t, 0, n)

proc merge2(t1: var seq[EigTable], t2: seq[EigTable], rrbs: int, op: any) =
  for i in range(t2.len):
    t1.add t2[i]

  let rrBack = 2
  var rrMax = 0
  var rrLen = rrbs

  while rrMax < t1.len-1:

    var rrMin = sortdown(t1, t1.len-1)
    rrMin = max(rrMin-rrBack, 0)
    rrMax = min(rrMin+rrLen, t1.len) - 1
    while true:
      let k = sortdown(t1, t1.len-1)
      if k >= rrMax: break

    echo "> rrMin: ", rrMin, "  rrMax: ", rrMax
    rayleighRitz(t1, rrMin, rrMax, op)
    rayleighRitz(t1, rrMin, rrMax, op)

proc merget(vt1: var seq[EigTable]; vt2: var seq[EigTable];
            ng,nmx,rrbs: int, op: any) =
  tic()
  var nt = vt1.len + vt2.len
  var ngd = min(ng,nt)
  var nmax = min(nmx,nt)
  var t = newSeq[type(vt1[0])](nt)
  var i1,i2: int
  for i in 0..<t.len:
    t[i].updates = 0
    if i2>=vt2.len or (i1<vt1.len and vt1[i1].sv<=vt2[i2].sv):
      t[i] = vt1[i1]
      inc i1
    else:
      t[i] = vt2[i2]
      inc i2
    var k = i - 1
    while k>0 and t[k].sv>t[k+1].sv:
      swap(t[k], t[k+1])
      dec k

  var rrLen = rrbs
  var rrMin = 0
  var rrMax = 0
  var to1 = 0.0
  var tr1 = 0.0

  toc("merge setup")
  #while rrMax < nt-1:
  while rrMin < ng:
    tic()
    rrMax = rrMin + rrLen - 1
    var rrMinNext = (rrMin + rrMax + 1) div 2
    if rrMax >= nt: rrMax = nt-1

    echo "rrMin: ", rrMin, "  rrMax: ", rrMax

    while true:
      tic()
      ortho1(t, rrMin, rrMax, op)
      var t0 = getElapsedTime()
      rayleighRitz(t, rrMin, rrMax, op)
      rayleighRitz(t, rrMin, rrMax, op)
      var t1 = getElapsedTime()
      #ortho1(t, rrMin, rrMax, op)
      var t2 = getElapsedTime()
      to1 += t0 + (t2-t1)
      tr1 += t1 - t0

      for i in rrMin..rrMax: inc t[i].updates
      for i in (rrMax+1)..<nt: t[i].updates = 0
      #if rrMax>=nt-1: break
      var s = sortdown(t, rrMin, nt-1)
      #if s>rrMax: break
      while s<min(nt-1,rrMax+rrLen):
        s = sortdown(t, rrMin, nt-1)
      break

    rrMin = 0
    while rrMin<nt and t[rrMin].updates>=1: inc rrMin
    rrMin = max(0,rrMin-5)
    #rrMin += rrLen div 2
    toc("merge loop end")

  toc("merge loop")

  var tm = getElapsedTime()
  echo "merget time = $1 secs (o1 $2 r1 $3)"%[tm|-6,to1|-6,tr1|-6]

  for i in 0..<vt2.len:
    if nmax+i<vt1.len: vt2[i].v = vt1[nmax+i].v
    else: vt2[i].v = nil
  vt1.setLen(nmax)
  for i in 0..<vt1.len: vt1[i] = t[i]
  toc("merge end")

proc svd(op: any, src: any, v: var any, sits: int, emin,emax: float) =
  let n = v.len
  var sv = newSeq[float](n)
  var qv = newSeq[type(v[0].even)](n)
  var qva = newSeq[type(v[0].odd)](0)
  var rsq = 0.0
  var lverb = 3
  for i in 0..<n:
    if v[i].isNil: v[i] = op.newVector()
    qv[i] = v[i].even
    #qva[i] = v[i].odd
  let nevs = svdLanczos(op, src.even, sv, qv, qva, rsq, sits,
                        emin, emax, lverb)
  echo "svd nevs: ", nevs

type EigOpts* = object
  nev*: int
  nvecs*: int
  relerr*: float
  abserr*: float
  svdits*: int
  maxup*: int
  rrbs*: int

proc initOpts*(eo: var EigOpts) =
  eo.nev = 10
  eo.nvecs = 20
  eo.relerr = 1e-4
  eo.abserr = 1e-6
  eo.svdits = 100
  eo.maxup = 10
  eo.rrbs = 20

## op: operator object
## opts: options onject
## vv: seq of vectors
proc hisqev*(op: any, opts: any, vv: any): auto =
  let ng = opts.nev
  let nvt = opts.nvecs
  let relerr = opts.relerr
  let abserr = opts.abserr
  let svdits = opts.svdits
  let maxup = opts.maxup
  let rrbs = opts.rrbs

  echo "starting hisqev"
  echo "ng = ", ng
  echo "nvt = ", nvt
  echo "relerr = ", relerr
  echo "abserr = ", abserr
  echo "svdits = ", svdits
  echo "maxup = ", maxup
  echo "rrbs = ", rrbs

  tic()
  let tt0 = getTics()

  var emin = 0.0
  var emax = 1e99
  var vt1 = newSeq[EigTable[type(op.newVector)]](ng)
  var vt2 = newSeq[EigTable[type(op.newVector)]](ng)
  var src = op.newVector
  var v: seq[type(op.newVector)]
  if not vv.isNil and vv.len>0:
    shallowCopy(v, vv)
    vt1.maketable(v, op)
  else:
    op.rand(src)
    #mixin `:=`
    #src := 0
    #src[0][0] := 1
    #let sits = 5*ng
    let sits = svdits
    v.newSeq(ng)
    svd(op, src, v, sits, emin, emax)
    vt1.maketable(v, op)

  geterr(vt1, op)
  for i in 0..<vt1.len:
    #cprintf("  %i\t%-18.12g%-18.12g%-18.12g\n", i, vt1[i].sv, vt1[i].err,
    #        vt1[i].err/vt1[i].sv)
    echo "  $1\t$2$3$4"%[$i, vt1[i].sv|(-18,12), vt1[i].err|(-18,12),
                         (vt1[i].err/vt1[i].sv)|(-18,12)]

  var ngcount = 0
  var re = relerr
  var ae = abserr
  for iter in 1..maxup:
    tic()
    var iv = 0
    while(iv<vt1.len-1 and (vt1[iv].err<ae or vt1[iv].err<re*vt1[iv].sv)):
      inc iv
    #cprintf("iv = %i  sv[iv] = %g\t%g\t%g\n", iv, vt1[iv].sv,
    #        vt1[iv].err, vt1[iv].err/vt1[iv].sv)
    echo "iv = $1  sv[iv] = $2\t$3\t$4"%[$iv, $vt1[iv].sv, $vt1[iv].err,
                                         $(vt1[iv].err/vt1[iv].sv)]
    if iv>ng:
      ngcount += 1
    else:
      ngcount = 0
    if ngcount>1: break
    let vin = vt1[iv].v
    emin = 0
    emax = 1e99
    if vt1.len>=ng: emax = vt1[min(vt1.len,nvt)-1].sv
    #cprintf("emin %g  emax %g\n", emin, emax)
    echo "emin $1  emax $2\n"%[$emin, $emax]
    op.rand(src)
    let srcn2 = src.norm2
    vin += (vt1[iv].sv*vt1[iv].sv/sqrt(srcn2)) * src
    let sits = svdits
    svd(op, vin, v, sits, emin, emax)
    toc("svd")
    GC_fullCollect()
    toc("GC")
    vt2.maketable(v, op)
    toc("maketable")
    merget(vt1, vt2, ng, nvt, rrbs, op)
    toc("merge")
    geterr(vt1, op)
    #cprintf("pass %i\n", iter)
    echo "pass ", iter
    for i in 0..<vt1.len:
      #cprintf("  %i\t%-18.12g%-18.12g%-18.12g\n", i, vt1[i].sv, vt1[i].err,
      #        vt1[i].err/vt1[i].sv)
      echo "  $1\t$2$3$4"%[$i, vt1[i].sv|(-18,12), vt1[i].err|(-18,12),
                           (vt1[i].err/vt1[i].sv)|(-18,12)]
    toc("end iter")
    let tt1 = getTics()
    #cprintf("iteration %i time = %.2f seconds\n", iter, toSeconds(tt1-tt0))
    echo "iteration $1 time = $2 seconds"%[$iter, toSeconds(tt1-tt0)|-6]
  toc()
  let t1 = getElapsedTime()
  #cprintf("total time = %.2f seconds\n", t1)
  echo "total time = $1 seconds"%[t1|-6]
  vt1

proc hisqev*(op: any, opts: any): auto =
  var v = newSeq[type(op.newVector)](0)
  result = hisqev(op, opts, v)

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
  var r: Field[1,RngMilc6]
  r.new(lo1)
  var seed = 987654321
  threads:
    for j in lo.sites:
      var l = lo.coords[lo.nDim-1][j].int
      for i in countdown(lo.nDim-2, 0):
        l = l * lo.physGeom[i].int + lo.coords[i][j].int
      r[j].seed(seed, l)

  type MyOp = object
    s: type(s)
    r: type(r)
    lo: type(lo)
  var op = MyOp(r:r,s:s,lo:lo)
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
    op.lo.ColorVector()
  template apply(op: MyOp, r,v: typed) =
    threadBarrier()
    stagD(op.s.so, r.field, op.s.g, v.field, 0.0)
  template applyAdj(op: MyOp, r,v: typed) =
    threadBarrier()
    stagD(op.s.se, r.field, op.s.g, v.field, 0.0, -1)
  template newRightVec(op: MyOp): untyped = newVector(op).even
  template newLeftVec(op: MyOp): untyped = newVector(op).odd

  var opts: EigOpts
  opts.initOpts
  opts.nev = intParam("nev", 200)
  opts.nvecs = intParam("nvecs", (opts.nev*11) div 10)
  opts.rrbs = intParam("rrbs", opts.nvecs)
  opts.relerr = 1e-4
  opts.abserr = 1e-6
  #opts.relerr = 1e-6
  #opts.abserr = 1e-8
  opts.svdits = intParam("svdits", opts.nev*2)
  opts.maxup = 100

  var evals = hisqev(op, opts)

  qexFinalize()

# nev: number of converged eigenvectors requested
# nvecs: number of vectors to keep in between passes
