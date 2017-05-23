import base
import math
import lapack
import linalgFuncs
import times
import strUtils

template QMP_time() = epochTime()
var verb: int

proc getSvals*(e: dvec; a: dvec; b: dvec; n: int) =
  #for i in 0..<n:
  #  echo i, ": ", a[i], "  ", b[i]
  svdbi(e.dat, a.dat, b.dat, n)

proc svd_bi3*(ev: dvec; m: dmat; ma: dmat; a: dvec; b: dvec) =
  var n = mat_nrows(m)
  var k = mat_ncols(m)
  var d = cast[ptr carray[float64]](alloc(n*sizeof(float64)))
  var e = cast[ptr carray[float64]](alloc(n*sizeof(float64)))
  for i in 0..<n:
    d[i] = a[i]
    if i<n-1: e[i] = b[i]
  svdBidiag(d, e, cast[ptr carray[float64]](addr m[0,0]),
            cast[ptr carray[float64]](addr ma[0,0]), n, k)
  for i in 0..<n:
    ev[i] = d[i]
  dealloc(d)
  dealloc(e)
  if verb>0:
    #template pv(i: int) =
    #  cprintf("sv[%i] %14.12g ", i, ev[i])
    #pv(0)
    #pv(1)
    #pv(k-1)
    #pv(n-1)
    #echo()
    template pv(i: int): untyped = "sv[$1] $2 "%[$i, ev[i]|(14,12)]
    echo pv(0), pv(1), pv(k-1), pv(n-1)

template nothreads(x: untyped): untyped =
  block: x

# A V = U B, Ad U = V Bd, B = [[a0,b0,0,...][0,a1,b1,0,...]...]
proc svdLanczos*(linop: any; src: any; sv: var any; qv: any; qva: any;
                 rsq: float; kmx: int; emin,emax: float, lverb: int): int =
  mixin `:=`
  verb = lverb
  var nv = qv.len
  var nva = qva.len
  var kmax = kmx
  var
    a: dvec
    b: dvec
    ev: dvec
  dvec_alloc(a, kmax)
  dvec_alloc(b, kmax)
  dvec_alloc(ev, kmax)
  var r = linop.newLeftVec
  var u = linop.newLeftVec
  var p = linop.newRightVec
  var v = linop.newRightVec

  tic()
  nothreads:
    let sn = sqrt(src.norm2)
    #echo "sn: ", sn
    v := src / sn
    p := v
    u := 0
  var alpha = 0.0
  var beta = 1.0
  var kcheck = 1
  var k = 0

  toc("setup")
  while true:
    tic()
    nothreads:
      tic()
      v := p / beta
      toc("loop1 thread1 exp1")
      #threadBarrier()
      linop.apply(r, v)
      toc("loop1 thread1 op")
      #threadBarrier()
      #echo "v: ", v.norm2
      #echo "r: ", r.norm2
      r -= beta*u
      let alph = sqrt(r.norm2)
      if threadNum==0:
        alpha = alph
        a[k] = alph
        inc k
      toc("loop1 thread1 end")
    toc("loop1 thread1")

    #check singular values
    if k >= kcheck or k >= kmax:
      getSvals(ev, a, b, k)
      if verb>0:
        #cprintf("%-5i", k)
        #cprintf(" sv%i %-16.12g", 0, dvec_get(ev, 0))
        #cprintf(" sv%i %-16.12g", 1, dvec_get(ev, 1))
        #cprintf(" sv%i %-16.12g", nv-1, dvec_get(ev, nv-1))
        #cprintf(" sv%i %-16.12g\n", k-1, dvec_get(ev, k-1))
        template sv(n): untyped = " sv$1 $2"%[$n,ev[n]|(-16,12)]
        echo k|-5, sv(0), sv(1), sv(nv-1), sv(k-1)
      kcheck = 1 + (1.5 * kcheck.float).int
      if k >= kmax: break
    toc("loop1 out")

    nothreads:
      tic()
      u := r/alpha
      toc("loop1 thread2 exp1")
      linop.applyAdj(p, u)
      toc("loop1 thread2 op")
      p -= alpha*v
      let bet = sqrt(p.norm2)
      if threadNum==0:
        beta = bet
        b[k-1] = bet
    #cprintf("%i\t%16.12g %16.12g\n", k-1, alpha, beta);
      toc("loop1 thread2 end")
    toc("loop1 thread2")

  toc("done iterations 1")
  var dtime1 = getElapsedTime()
  if verb>0:
    #cprintf("svd_lanczos %g secs\n", dtime1)
    echo "svd_lanczos $1 secs"%[dtime1|-6]

  kmax = k
  if nv > kmax:
    for i in kmax..<nv:
      sv[i] = 1e99
    nv = kmax
  for i in 0..<nv:
    sv[i] = ev[i]
  var
    vr: dmat
    ur: dmat
    #vr2: dmat
    #ur2: dmat
  dmat_alloc(vr, kmax, nv)
  dmat_alloc(ur, kmax, nva)
  #dmat_alloc(vr2, kmax, nv)
  #dmat_alloc(ur2, kmax, nva)
  #svd_bi3(ev, vr, ur, a, b)
  #var nvout = nv
  var nvout = svdBi4(ev, vr, ur, a, b, kmax, nv, nva, emin, emax)
  nv = min(nv,nvout)
  nva = min(nva,nvout)
  #var s2 = 0.0
  #for i in 0..<vr.nrows:
  #  #for j in 0..<vr.ncols:
  #  for j in 0..<1:
  #    let d = vr[i,j]*vr2[0,j] - vr2[i,j]*vr[0,j]
  #    s2 += d*d
  #echo s2

  toc("svd")
  var dtime2 = getElapsedTime()
  if verb>0:
    #cprintf("svd_lanczos %g secs %g\x0A", dtime2-dtime1, dtime2)
    echo "svd_lanczos $1 secs $2"%[(dtime2-dtime1)|-6, dtime2|-6]

  for i in 0..<qv.len:
    qv[i] := 0
  for i in 0..<qva.len:
    qva[i] := 0
  beta = sqrt(src.norm2)
  v := src / beta
  p := v
  beta = 1.0
  u := 0
  k = 0
  while true:
    tic()
    v := p/beta
    toc("loop2 eq1")
    #cprintf("v[%i,0]: %.12g\n", k, vr[k,0]/vr[0,0])
    for i in 0..<nv:
      qv[i] += vr[k,i] * v
    toc("loop2 qv")
    linop.apply(r, v)
    toc("loop2 linop1")
    r -= beta*u
    let alpha = sqrt(r.norm2)
    u := r/alpha
    toc("loop2 eq2")
    for i in 0..<nva:
      qva[i] += ur[k,i] * u
    toc("loop2 qva")
    inc k
    if k >= kmax: break
    linop.applyAdj(p, u)
    toc("loop2 linop2")
    p -= alpha * v
    beta = sqrt(p.norm2)
    toc("loop 2 end")

  toc("done")
  var dtime3 = getElapsedTime()
  if verb>0:
    #cprintf("svd_lanczos %g secs %g\x0A", dtime3-dtime2, dtime3)
    echo "svd_lanczos $1 secs $2"%[(dtime3-dtime2)|-6, dtime3|-6]
  result = nvout


when isMainModule:
  import qex
  import qcdTypes
  import gaugeUtils
  import stagD
  import rng

  qexInit()
  #var defaultGaugeFile = "l88.scidac"
  var defaultLat = [8,8,8,8]
  #var defaultLat = [16,16,16,16]
  defaultSetup()
  threads:
    g.setBC
    g.stagPhase
  var s = newStag(g)

  var lo1 = newLayout(lat, 1)
  var r: Field[1,RngMilc6]
  r.new(lo1)
  var seed = 987654321
  threads:
    for s in lo1.sites:
      var l = lo1.coords[lo1.nDim-1][s].int
      for i in countdown(lo1.nDim-2, 0):
        l = l * lo1.physGeom[i].int + lo1.coords[i][s].int
      r[s].seed(seed, l)

  type MyOp = object
    s: type(s)
    r: type(r)
    lo: type(lo)
  var op = MyOp(s:s,r:r,lo:lo)
  proc gaussian(x: C1, r: var any) =
    x.re = gaussian(r)
    x.im = gaussian(r)
  proc gaussian(x: Vec1, r: var any) =
    for i in 0..<x.len:
      gaussian(x[i], r)
  proc gaussian(v: SomeField, r: SomeField2) =
    for i in v.sites:
      gaussian(v{i}, r[i])
  template rand(op: MyOp, v: any) =
    gaussian(v, op.r)
  template newVector(op: MyOp): untyped =
    op.lo.ColorVector()
  template apply(op: MyOp, r,v: typed) =
    stagD(op.s.so, r.field, op.s.g, v.field, 0.0)
  template applyAdj(op: MyOp, r,v: typed) =
    stagD(op.s.se, r.field, op.s.g, v.field, 0.0, -1)
  template newRightVec(op: MyOp): untyped = newVector(op).even
  template newLeftVec(op: MyOp): untyped = newVector(op).odd

  var src = op.newRightVec
  op.rand(src)
  var nv = 20
  var nva = nv
  var sv = newDvec(nv)
  var qv = newSeq[type(op.newRightVec)](nv)
  var qva = newSeq[type(op.newLeftVec)](nv)
  for i in 0..<nv:
    qv[i] = op.newRightVec()
    qva[i] = op.newLeftVec()
  var rsq = 0.0
  var kmax = 100
  var lverb = 3
  var k = svdLanczos(op, src, sv, qv, qva, rsq, kmax, lverb)

  op.apply(src, qv[0])
  echo "sv: ", sv[0], "  ave: ", qva[0].dot(src.field)
  op.applyAdj(src, src)
  let vn = qv[0].norm2
  let Dvn = src.field.odd.norm2
  echo "l: ", sqrt(Dvn/vn)

  qexFinalize()
