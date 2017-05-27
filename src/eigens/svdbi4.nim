import math
import algorithm
import strUtils
import base
import eigens/linalgFuncs

# A V = U B, Ad U = V Bd, B = [[a0,b0,0,...][0,a1,b1,0,...]...]
# B  v = u ev
# Bd u = v ev

proc checksv(v: dvec, u: dvec, a: dvec, b: dvec, ev: float): float =
  let n = v.len
  for i in 0..<n:
    var x0 = ev*v[i]
    var y0 = a[i]*u[i]
    if i>0: y0 += b[i-1]*u[i-1]
    var d0 = x0 - y0
    result += d0*d0
    var x1 = ev*u[i]
    var y1 = a[i]*v[i]
    if i<n-1: y1 += b[i]*v[i+1]
    var d1 = x1 - y1
    result += d1*d1

proc dvsort(tg: dvec, itg: var seq[int]) =
  let n = tg.len
  for i in 0..<n: itg[i] = i
  proc sortdv(x,y: int): int = cmp(tg[x], tg[y])
  itg.sort(sortdv)

proc getTwistFac*(a: dvec; b: dvec; ev: float64;
                  tl: dvec; tu: dvec; tg: dvec): int =
  var
    ai: float64
    bi: float64
    lp: float64
    dp: float64
    g: float64
    gmin = 1e99
    e2 = ev * ev
    sp = -e2
    n = tg.len
    imin = -1
  for i in 0..(n-2):
    tg[i] = sp + e2
    ai = a[i]
    bi = b[i]
    dp = ai * ai + sp
    lp = ai * bi / dp
    tl[i] = lp
    sp = lp * sp * bi / ai - e2
  tg[n-1] = sp + e2
  ai = a[n-1]
  # dp =  ai*ai + sp;
  #  sp=pp, dp=rp, lp=up
  sp = ai * ai - e2
  g = sp + tg[n-1]
  if abs(g) < gmin:
    gmin = abs(g)
    imin = n - 1
  tg[n-1] = g
  for i in countdown(n-2,0):
    ai = a[i]
    bi = b[i]
    dp = bi * bi + sp
    lp = ai * bi / dp
    tu[i] = lp
    sp = sp * ai * ai / dp - e2
    g = sp + tg[i]
    if abs(g) < gmin:
      gmin = abs(g)
      imin = i
    tg[i] = g
  # dp = sp;
  return imin

proc twistSolve*(tl: dvec; tu: dvec; v: dvec; k: int) =
  var s = 1.0
  var q = 1.0
  var n = v.len
  v[k] = q
  for i in countdown(k-1,0):
    q = - (tl[i] * q)
    v[i] = q
    s += q * q
  q = 1.0
  for i in (k+1)..<n:
    q = - (tu[i-1] * q)
    v[i] = q
    s += q * q
  s = 1.0 / sqrt(s)
  for i in 0..<n:
    v[i] = s * v[i]

proc getu*(a: dvec; b: dvec; v: dvec; u: dvec) =
  var
    t = 0.0
    s = 0.0
    k = v.len
  for i in 0..(k-2):
    t = a[i]*v[i] + b[i]*v[i+1]
    u[i] = t
    s += t * t
  t = a[k-1] * v[k-1]
  u[k-1] = t
  s += t * t
  s = 1.0/sqrt(s)
  for i in 0..<k:
    u[i] = s*u[i]

## svd of bidiagonal matrix
## e: singular values
## m: A ma = e m
## k: number of singular values wanted
proc svdBi4*(e: dvec; m: dmat; ma: dmat; a: dvec; b: dvec; k: int;
             nx: int; nax: int; emin: float64; emax: float64): int =
  #svdbi(&e, &a, &b, k)

  var n = nx
  var na = nax
  var nn = 0
  if n > 0 or na > 0:
    var imin = -1
    var ne = 0
    for i in 0..<k:
      var ee = e[i]
      if emin <= ee and ee <= emax:
        if imin < 0: imin = i
        inc ne
    if n > ne: n = ne
    if na > ne: na = ne
    nn = max(n, na)
    # printf0("nn %i n %i na %i ne %i\n", nn, n, na, ne);

    if nn > 0:
      var
        tl = newDvec(k-1)
        tu = newDvec(k-1)
        tg = newDvec(k)
        v = newSeq[dvec](nn)
        va = newSeq[dvec](nn)
      var dotstop = 1e-4
      var ib = 0
      var itg = newSeq[int](k)
      for i in 0..<nn:
        var ie = imin + i
        var ev = e[ie]
        var rsqstop = 1e-6 * ev * ev
        v[i].colvec(m, i)
        va[i].colvec(ma, i)
        var i0 = getTwistFac(a, b, ev, tl, tu, tg)
        var ig = 0
        twistSolve(tl, tu, v[i], i0)
        while true:
          var x = 0.0
          for j in ib..<i:
            var xx = dot(v[j], v[i])
            daxpy(-xx, v[j], v[i])
            x += abs(xx)
          normalize(v[i])
          getu(a, b, v[i], va[i])
          var rsq = checksv(v[i], va[i], a, b, ev)
          #echo i, ": ", e[i], "  ", rsq
          if ig > 0:
            #printf0("svdvec %3i %3i %5i %13g %13g %13g %13g\x0A",
            #        i, ib, ig, ev, x, rsq, dvec_get(tg, i0))
            echo "svdvec $# $# $# $# $# $# $#" %
              [$i, $ib, $ig, $ev, $x, $rsq, $tg[i0]]
          if ig > 10 or rsq < rsqstop:
            if ig == 0 and x < dotstop: ib = i
            break
          if ig == 0:
            dvsort(tg, itg)
          inc ig
          i0 = itg[ig]
          twistSolve(tl, tu, v[i], i0)
          # getu(a, b, v, va);
  return nn


when isMainModule:
  import base/basicOps
  import maths/complexType
  template adj*(x: SomeNumber): untyped = x
  import times

  proc echoCol(m: dmat, i: int) =
    var s = ""
    for j in 0..<m.nrows:
      s &= $m[j,i] & " "
    echo s

  proc testSvdbd(n: int) =
    var v = newSeq[float](n)
    var d = newSeq[float](n)
    var e = newSeq[float](n)
    template `&`(x: seq): untyped = cast[ptr carray[float]](addr x[0])
    for i in 0..<n:
      d[i] = (i+1).float
      e[i] = -(i+1).float
    svdbi(&v, &d, &e, n)
    svdbi(&v, &d, &e, n)
    let t0 = epochTime()
    svdbi(&v, &d, &e, n)
    let t1 = epochTime()
    #for i in 0..<n:
    #  echo i, ": ", v[i]
    echo v[0], "  ", v[n-1]
    echo n, " time: ", t1-t0
  #testSvdbd(100)
  #testSvdbd(200)
  #testSvdbd(400)
  #testSvdbd(800)
  #testSvdbd(1600)
  #testSvdbd(3200)
  #testSvdbd(6400)
  #testSvdbd(12800)

  proc testSvdbdv(n: int): auto =
    var s = newDvec(n)
    var d = newDvec(n)
    var e = newDvec(n-1)
    var v = newDmat(n, n)
    var u = newDmat(n, n)
    template `&`(x: seq): untyped = cast[ptr carray[float]](addr x[0])
    for i in 0..(n-2):
      d[i] = (i+1).float
      e[i] = -(i+1).float
    d[n-1] = (n).float

    let t0 = epochTime()
    var nn = svdBi4(s, v, u, d, e, n, n, n, 0.0, 9999999.0)
    let t1 = epochTime()
    #for j in 0..<nr:
    #    v[i+nr*j]
    echo s[0], "  ", s[n-1]
    echoCol(v,0)
    echoCol(v,1)
    echoCol(v,2)
    echoCol(v,3)
    echo nn, " time: ", t1-t0
    v
  #testSvdbdv(100)

  proc testSvdbdv2(n: int): auto =
    var s = newDvec(n)
    var d = newDvec(n)
    var e = newDvec(n-1)
    var v = newDmat(n, n)
    var u = newDmat(n, n)
    template `&`(x: dvec): untyped = cast[ptr carray[float]](addr x[0])
    template `&`(x: dmat): untyped = cast[ptr carray[float]](addr x[0,0])
    for i in 0..(n-2):
      d[i] = (i+1).float
      e[i] = -(i+1).float
    d[n-1] = (n).float

    let t0 = epochTime()
    svdBidiag(&d, &e, &v, &u, n, n)
    let t1 = epochTime()
    #for j in 0..<nr:
    #    v[i+nr*j]
    echo d[0], "  ", d[n-1]
    echoCol(v,0)
    echoCol(v,1)
    echoCol(v,2)
    echoCol(v,3)
    echo n, " time: ", t1-t0
    v

  template test2(n: int) =
    #testSvdbd(n)
    #testSvdbdv(n)
    #testSvdbdv2(n)
    let v1 = testSvdbdv(n)
    let v2 = testSvdbdv2(n)
    var s2 = 0.0
    for i in 0..<v1.nrows:
      for j in 0..<v1.ncols:
        let d = v1[i,j]*v2[0,j] - v2[i,j]*v1[0,j]
        s2 += d*d
    echo s2

  test2(4)
  #test2(10)
  #test2(100)
  #test2(1000)
  #test2(10000)
