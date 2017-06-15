import math
import linalgFuncs

proc netri(d: any, e: any, x: float): int =
  let n = d.len
  let eps = 1e-20
  var s0 = -x
  var t0 = 1.0
  if s0<=eps:
    s0 = min(s0, -eps)
    inc result
  var s = -s0*(x*s0 + d[0]*d[0])
  var t = s0*s0
  if s<=eps*t:
    s = min(s, -eps*t)
    inc result
  for i in 1..<n:
    #s0 = -x - e[i-1]*e[i-1]/s
    #if s0<=eps:
    #  s0 = min(s0, -eps)
    #  inc result
    #s = -x - d[i]*d[i]/s0
    s0 = -s*(x*s + e[i-1]*e[i-1]*t)
    t0 = s*s
    if s0<=eps*t0:
      s0 = min(s0, -eps*t0)
      inc result
    s = -s0*(x*s0 + d[i]*d[i]*t0)
    t = s0*s0
    if s<=eps*t:
      s = min(s, -eps*t)
      inc result
    #echo i, ": ", s, " : ", t
    if abs(s)<1e-16 or abs(s)*1e-16>1.0:
      s = s/t
      t = 1.0

proc fz*(d: any, e: any, x0,x1: float, n: int): auto =
  template f(x: float): int = netri(d, e, x)
  var xm = x0
  var ym = f(xm)
  #echo "xm: ", xm, " : ", ym
  while ym > n:
    xm -= abs(x1-xm)
    ym = f(xm)
  var xp = x1
  var yp = f(xp)
  #echo "xp: ", xp, " : ", yp
  while yp <= n:
    xp += abs(xp-x0)
    yp = f(xp)
  #echo "xm: ", xm, "  xp: ", xp
  var xp0 = xp
  var yp0 = yp
  while true:
    var xn = 0.5*(xm+xp)
    if xn==xm or xn==xp: break
    var yn = f(xn)
    if yn <= n:
      xm = xn
    else:
      xp = xn
      if yp0>n+1:
        xp0 = xn
        yp0 = yn
    #echo "xm: ", xm, "  xp: ", xp, "  yn: ", yn
  (0.5*(xm+xp), xp0)

proc svbinX*(sv: var any, d: any, e: any, ioff,nmn,nmx: int,
             x0: float, y0: int, x1: float, y1: int): int =
  if nmn>=y1: return
  if nmn<y0:
    echo "error: nmn<y0: ", nmn, " < ", y0
    quit(-1)
  var xn = 0.5*(x0+x1)
  if xn==x0 or xn==x1:
    let nn = y1 - y0
    for i in 0..<nn:
      sv[ioff+i] = xn
    return nn
  var yn = netri(d, e, xn)
  var nd0,nd1 = 0
  if nmn<yn:
    nd0 = svbinX(sv, d, e, ioff, nmn, nmx, x0, y0, xn, yn)
  if nmn+nd0<y1 and yn<nmx:
    nd1 = svbinX(sv, d, e, ioff+nd0, nmn+nd0, nmx, xn, yn, x1, y1)
  result = nd0 + nd1

proc svbin*(sv: var any, d: any, e: any, off,imin,imax: int) =
  ## find the smallest singular values of a bidiagonal matrix
  let n = d.len
  let ns = imax - imin + 1
  var nd = 0
  var x0 = 0.0
  var x1 = 1.0
  var n0 = n + imin
  var n1 = n + imax + 1
  var ioff = off
  while true:
    var y0 = netri(d, e, x0)
    var y1 = netri(d, e, x1)
    var nd0 = svbinX(sv, d, e, ioff, n0, n1, x0, y0, x1, y1)
    nd += nd0
    if nd >= ns: break
    ioff += nd0
    n0 += nd0
    x0 = x1
    x1 = x0 + 1.0

template `&`(x: seq[float]): untyped = cast[ptr carray[float]](addr x[0])

when isMainModule:
  import times
  var n = 10000
  var ns = min(n, 100)
  var d = newSeq[float](n)
  var e = newSeq[float](n-1)
  var sv = newSeq[float](ns)
  var sv2 = newSeq[float](n)
  for i in 0..<n:
    d[i] = (i+1).float / (n+1).float
  for i in 0..(n-2):
    e[i] = 0.0

  #var nt = 15
  #for i in 0..nt:
  #  let x = i.float / nt.float
  #  let y = svbidet(d, e, x)
  #  echo x, ": ", y

  template chk =
    for i in 0..<ns:
      let d = sv[i] - sv2[i]
      if abs(d) > 1e-15:
        echo "error: ", sv2[i], ": ", sv[i], " : ", d

  template test =
    let t0 = epochTime()
    svbin(sv, d, e, 0, 0, ns-1)
    let t1 = epochTime()
    svdbi(&sv2, &d, &e, n)
    let t2 = epochTime()
    chk()
    echo "svbiz: ", t1-t0, "  svdbi: ", t2-t1

  test()
  #for i in 0..n:
  #  let t = i.float/n.float
  #  echo t, " : ", netri(d, e, t)

  for i in 0..(n-2):
    e[i] = (i+0.5).float / n.float

  #for i in 0..n:
  #  let t = i.float/n.float
  #  echo t, " : ", netri(d, e, t)

  test()
