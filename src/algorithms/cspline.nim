## cubic spline and related routines

import std/algorithm
import base

proc estimateDerivative*[N:static[int],T](dx,dy:array[N,T]):T =
  ## estimate the derivative given a list of dx and dy, using Taylor series
  ## dx and dy must be the difference away from a single point.
  when N==1:
    return dy[0]/dx[0]
  elif N==2:
    let h10 = dx[1]-dx[0]
    return dy[0]*dx[1]/(dx[0]*h10) - dx[0]*dy[1]/(h10*dx[1])
  elif N==3:
    let
      h10 = dx[1]-dx[0]
      h02 = dx[0]-dx[2]
      h21 = dx[2]-dx[1]
      a = dy[0]*dx[1]*dx[2]/(dx[0]*h10*h02)
      b = dx[0]*dy[1]*dx[2]/(dx[1]*h10*h21)
      c = dx[0]*dx[1]*dy[2]/(dx[2]*h02*h21)
    return -(a+b+c)
  else:
    error("estimateDerivative: Unimplemented for N = " & $N)

type CSpline*[T] = object
  x: seq[T]
  ys: seq[array[2,T]]  ## y and computed second derivatives of y

type
  CSplineBoundDyKind = enum
    CSBEstimateDy, CSBZeroD2y, CSBSetDy
  CSplineBoundDy = object
    case kind: CSplineBoundDyKind
    of CSBEstimateDy, CSBZeroD2y:
      discard
    of CSBSetDy:
      dy: float
  CSplineBounds = object
    lo,hi: CSplineBoundDy

const CSplineBoundEstimateDy* = CSplineBoundDy(kind:CSBEstimateDy)
const CSplineBoundZeroD2y* = CSplineBoundDy(kind:CSBZeroD2y)
proc dyBound*(dy:float):CSplineBoundDy = CSplineBoundDy(kind:CSBSetDy,dy:dy)  ## Fixed endpoint first derivative.

proc csplineBounds*(lo=CSplineBoundEstimateDy, hi=CSplineBoundEstimateDy):CSplineBounds =
  CSplineBounds(lo:lo, hi:hi)

func triSolve[T](a,b,c,r:openarray[T]):seq[T] =
  ## Solve a tridiagonal system. a[0] and c[^1] are unused.
  let n = b.len
  result = newseq[T](n)
  var gam = newseq[T](n)
  var bet = b[0]
  result[0] = r[0]/bet
  for j in 1..<n:
    gam[j] = c[j-1]/bet
    bet = b[j] - a[j]*gam[j]
    result[j] = (r[j] - a[j]*result[j-1])/bet
  for j in countdown(n-2,0):
    result[j] -= gam[j+1]*result[j+1]

proc newCSpline*[T](x,y:openarray[T], bounds=csplineBounds()):CSpline[T] =
  let n = x.len
  if y.len != n:
    qexError "different length in x and y: ",n," != ",y.len
  var r = CSpline[T](x:newseq[T](n), ys:newseq[array[2,T]](n))
  for i in 0..<n:
    r.ys[i] = [y[i], x[i]]
  r.ys.sort do (a,b:array[2,T]) -> int:
    cmp(a[1], b[1])  # sort by x
  for i in 0..<n:
    r.x[i] = r.ys[i][1]  # copy sorted x

  template x(i:int):auto = r.x[i]
  template y(i:int):auto = r.ys[i][0]
  template d2y(i:int):auto = r.ys[i][1]
  var
    a = newseq[T](n)
    b = newseq[T](n)
    c = newseq[T](n)
    rhs = newseq[T](n)

  if bounds.lo.kind==CSBZeroD2y:
    b[0] = T(1.0)
    rhs[0] = T(0.0)
  else:
    let
      dy =
        if bounds.lo.kind==CSBEstimateDy:
          if n>3: estimateDerivative([x(1)-x(0), x(2)-x(0), x(3)-x(0)], [y(1)-y(0), y(2)-y(0), y(3)-y(0)])
          elif n==3: estimateDerivative([x(1)-x(0), x(2)-x(0)], [y(1)-y(0), y(2)-y(0)])
          elif n==2: estimateDerivative([x(1)-x(0)], [y(1)-y(0)])
          else: T(0.0)
        else:
          bounds.lo.dy
      d = y(1)-y(0)
      h = x(1)-x(0)
    b[0] = T(1.0)
    c[0] = T(0.5)
    rhs[0] = T(3.0)*(d/h-dy)/h

  for j in 1..<n-1:
    let
      xm = x(j-1)
      xj = x(j)
      xp = x(j+1)
      ym = y(j-1)
      yj = y(j)
      yp = y(j+1)
      hm = xj-xm
      hj = xp-xj
      hjm = hj/hm
      dhm = (yj-ym)/hm
      dhj = (yp-yj)/hj
      bj = T(2.0)*(T(1.0)+hjm)
    a[j] = T(1.0)
    b[j] = bj
    c[j] = hjm
    rhs[j] = T(6.0)*(dhj-dhm)/hm

  if bounds.hi.kind==CSBZeroD2y:
    b[n-1] = T(1.0)
    rhs[n-1] = T(0.0)
  else:
    let
      dy =
        if bounds.hi.kind==CSBEstimateDy:
          if n>3: estimateDerivative([x(n-2)-x(n-1), x(n-3)-x(n-1), x(n-4)-x(n-1)], [y(n-2)-y(n-1), y(n-3)-y(n-1), y(n-4)-y(n-1)])
          elif n==3: estimateDerivative([x(n-2)-x(n-1), x(n-3)-x(n-1)], [y(n-2)-y(n-1), y(n-3)-y(n-1)])
          elif n==2: estimateDerivative([x(n-2)-x(n-1)], [y(n-2)-y(n-1)])
          else: T(0.0)
        else:
          bounds.hi.dy
      d = y(n-1)-y(n-2)
      h = x(n-1)-x(n-2)
    a[n-1] = T(1.0)
    b[n-1] = T(2.0)
    rhs[n-1] = T(6.0)*(dy-d/h)/h

  let m = triSolve(a,b,c,rhs)
  for j in 0..<n:
    d2y(j) = m[j]

  return r

func bisect*[T](xs:openarray[T], x:T): int =
  ## Segment i with xs[i] <= x, clamped so xs[i+1] is valid. Assumes sorted xs.
  var
    lo = 0
    hi = xs.len - 1
  while hi - lo > 1:
    let mid = (lo + hi) div 2
    if xs[mid] <= x: lo = mid
    else: hi = mid
  lo

func interpolate*[T](csp:CSpline[T], x:T):T =
  let
    i = csp.x.bisect x
    x0 = csp.x[i]
    x1 = csp.x[i+1]
    y0 = csp.ys[i][0]
    d2y0 = csp.ys[i][1]
    y1 = csp.ys[i+1][0]
    d2y1 = csp.ys[i+1][1]
    h = x1-x0
    a = (x1-x)/h
    b = (x-x0)/h
    c = (a*a*a-a)*h*h/T(6.0)
    d = (b*b*b-b)*h*h/T(6.0)
  a*y0 + b*y1 + c*d2y0 + d*d2y1

func interpolateDy*[T](csp:CSpline[T], x:T):T =
  let
    i = csp.x.bisect x
    x0 = csp.x[i]
    x1 = csp.x[i+1]
    y0 = csp.ys[i][0]
    d2y0 = csp.ys[i][1]
    y1 = csp.ys[i+1][0]
    d2y1 = csp.ys[i+1][1]
    h = x1-x0
    a = (x1-x)/h
    b = (x-x0)/h
    c = (a*a*a-a)*h*h/T(6.0)
    d = (b*b*b-b)*h*h/T(6.0)
  (y1-y0)/h - (T(3.0)*a*a-T(1.0))*h*d2y0/T(6.0) + (T(3.0)*b*b-T(1.0))*h*d2y1/T(6.0)

func interpolateD2y*[T](csp:CSpline[T], x:T):T =
  let
    i = csp.x.bisect x
    x0 = csp.x[i]
    x1 = csp.x[i+1]
    y0 = csp.ys[i][0]
    d2y0 = csp.ys[i][1]
    y1 = csp.ys[i+1][0]
    d2y1 = csp.ys[i+1][1]
    h = x1-x0
    a = (x1-x)/h
    b = (x-x0)/h
    c = (a*a*a-a)*h*h/T(6.0)
    d = (b*b*b-b)*h*h/T(6.0)
  a*d2y0 + b*d2y1

# --- periodic uniform cubic spline ---------------------------------------------

func cyclicSolve[T](a,b,c:openarray[T]; alpha,beta:T; r:openarray[T]):seq[T] =
  ## Solve a tridiagonal system with corner entries A[0][^1] and A[^1][0].
  let n = b.len
  let gamma = -b[0]
  var bb = newseq[T](n)
  bb[0] = b[0] - gamma
  bb[n-1] = b[n-1] - alpha*beta/gamma
  for i in 1..<n-1: bb[i] = b[i]
  let x = triSolve(a, bb, c, r)
  var u = newseq[T](n)
  u[0] = gamma
  u[n-1] = alpha
  let z = triSolve(a, bb, c, u)
  let fact = (x[0] + beta*x[n-1]/gamma) / (T(1.0) + z[0] + beta*z[n-1]/gamma)
  result = newseq[T](n)
  for i in 0..<n: result[i] = x[i] - fact*z[i]

proc newCSplinePeriodic*[T](y:openarray[T]; x0,h:T):CSpline[T] =
  ## Uniform periodic cubic spline. Values wrap, so y[n] == y[0].
  let n = y.len
  if n < 3:
    qexError "periodic cubic spline needs at least 3 knots, got ", n
  var r = newseq[T](n)
  for j in 0..<n:
    let jm = (j + n - 1) mod n
    let jp = (j + 1) mod n
    r[j] = T(6.0)*(y[jp] - T(2.0)*y[j] + y[jm])/(h*h)
  var a = newseq[T](n)
  var b = newseq[T](n)
  var c = newseq[T](n)
  for j in 0..<n:
    a[j] = T(1.0); b[j] = T(4.0); c[j] = T(1.0)
  let m = cyclicSolve(a, b, c, T(1.0), T(1.0), r)
  result = CSpline[T](x: newseq[T](n+1), ys: newseq[array[2,T]](n+1))
  for j in 0..<n:
    result.x[j] = x0 + T(j)*h
    result.ys[j] = [y[j], m[j]]
  result.x[n] = x0 + T(n)*h
  result.ys[n] = [y[0], m[0]]

func segmentCoeffs*[T](csp:CSpline[T]; i:int):tuple[a,b,c,d:T] =
  ## Coefficients for S(x) = a + b*u + c*u^2 + d*u^3, u = x - csp.x[i].
  let
    h = csp.x[i+1] - csp.x[i]
    y0 = csp.ys[i][0]
    y1 = csp.ys[i+1][0]
    m0 = csp.ys[i][1]
    m1 = csp.ys[i+1][1]
  result.a = y0
  result.b = (y1 - y0)/h - h*(T(2.0)*m0 + m1)/T(6.0)
  result.c = T(0.5)*m0
  result.d = (m1 - m0)/(T(6.0)*h)

when isMainModule:
  import qex
  import utils/test

  proc fun0(x:float):auto =
    return (1.0+x, 1.0, 0.0, 0.0)
  proc fun1(x:float):auto =
    return ((1.0+x)*(2.0-x), 1.0-2.0*x, -2.0, 0.0)
  proc fun2(x:float):auto =
    return ((1.0+x)*(2.0-x)*(1.0-x), (3.0*x-4.0)*x-1.0, 6.0*x-4.0, 6.0)

  proc testEstD(test:QEXTest, ord:int, dx,dy:array[3,float], actual:float) =
    let d = [estimateDerivative([dx[0]],[dy[0]]),
             estimateDerivative([dx[0],dx[1]],[dy[0],dy[1]]),
             estimateDerivative(dx,dy)]
    let test = test.newTest("estimate derivative")
    for o in ord..3:
      test.assertAlmostEqual(d[o-1], actual)

  proc testCSp(test:QEXTest, spline:CSPline[float], ord:int, f:proc, checkValues=true) =
    let n = spline.x.len
    for i in 0..<n:
      test.logInfo i," x: ",spline.x[i]," y: ",spline.ys[i][0]," y'': ",spline.ys[i][1]
      let fx = f(spline.x[i])
      test.logInfo "  exact y': ",fx[1],"  y'': ",fx[2],"  y''': ",fx[3]
      #if i<n-1:
      #  let h = spline.x[i+1]-spline.x[i]
      #  let d = spline.ys[i+1][0]-spline.ys[i][0]
      #  echo "  y'",i," ",d/h+h*spline.ys[i][1]/(-3.0)+h*spline.ys[i+1][1]/(-6.0)
      #  echo "  y'",i+1," ",d/h+h*spline.ys[i][1]/6.0+h*spline.ys[i+1][1]/3.0
    let testcontdy = test.newTest("Continuous Derivatives", hidden=1)
    for i in 1..<n-1:
      let
        hm = spline.x[i]-spline.x[i-1]
        hp = spline.x[i+1]-spline.x[i]
        dym = (spline.ys[i][0]-spline.ys[i-1][0])/hm + hm*spline.ys[i-1][1]/6.0 + hm*spline.ys[i][1]/3.0
        dyp = (spline.ys[i+1][0]-spline.ys[i][0])/hp + hp*spline.ys[i][1]/(-3.0) + hp*spline.ys[i+1][1]/(-6.0)
      testcontdy.assertAlmostEqual(dyp, dym)
    if checkValues:
      for x in [spline.x[0], spline.x[n-1], spline.x[0]+0.05, 0.0, spline.x[n-1]-0.05]:
        let
          testp = test.newTest("x=" & $x, hidden=1)
          yi = spline.interpolate(x)
          dyi = spline.interpolateDy(x)
          d2yi = spline.interpolateD2y(x)
          (y, dy, d2y, d3y) = f(x)
        if ord<4:
          testp.newTest("y", hidden=1).assertAlmostEqual(yi, y)
          testp.newTest("dy", hidden=1).assertAlmostEqual(dyi, dy)
          testp.newTest("d2y", hidden=1).assertAlmostEqual(d2yi, d2y)

  proc run(test:QEXTest, ord:int, f:proc) =
    let test = test.newTest("polynomial degree " & $ord)
    let
      n = 7
      m = 4
    var
      xs = newseq[float](n+m)
      ys = newseq[float](n+m)
      dys = newseq[array[3,float]](n+m)
    for i in 0..<n:
      let x = float(i)*5.0/float(n-1) - 2.0
      let fx = f(x)
      xs[i] = x
      ys[i] = fx[0]
      dys[i][0] = fx[1]
      dys[i][1] = fx[2]
      dys[i][2] = fx[3]
    for i in 0..<m:
      let x = float(i)*5.0/float(m-1) - 1.9
      let fx = f(x)
      xs[n+i] = x
      ys[n+i] = fx[0]
      dys[n+i][0] = fx[1]
      dys[n+i][1] = fx[2]
      dys[n+i][2] = fx[3]
    testEstD(test, ord, [xs[3]-xs[2],xs[4]-xs[2],xs[5]-xs[2]], [ys[3]-ys[2],ys[4]-ys[2],ys[5]-ys[2]], dys[2][0])
    testCSp(
      test.newTest("cspline default (est. 1st deriv.)"),
      newCSpline(xs,ys),
      ord, f)
    testCSp(
      test.newTest("cspline set 1st deriv. bounds"),
      newCSpline(xs,ys,csplineBounds(dyBound(dys[0][0]),dyBound(dys[^1][0]))),
      ord, f)
    testCSp(
      test.newTest("cspline natural (zero 2nd deriv.)"),
      newCSpline(xs,ys,csplineBounds(CSplineBoundZeroD2y,CSplineBoundZeroD2y)),
      ord, f, checkValues=false)

  proc runPeriodic(test:QEXTest) =
    ## Periodic spline for a smooth periodic function.
    let test = test.newTest("periodic cubic spline")
    const n = 64
    let x0 = -PI
    let h = 2.0*PI/float(n)
    proc f(p:float):tuple[v,d1:float] =
      let v = exp(0.5*cos(p) + 0.3*sin(p))
      (v, (-0.5*sin(p) + 0.3*cos(p))*v)
    var ys = newseq[float](n)
    for j in 0..<n: ys[j] = f(x0 + float(j)*h).v
    let csp = newCSplinePeriodic(ys, x0, h)

    let tknot = test.newTest("knot values")
    for j in 0..<n:
      tknot.assertAlmostEqual(csp.interpolate(x0 + float(j)*h), ys[j])

    let tseg = test.newTest("segment coeffs vs interpolate")
    for j in 0..<n:
      let (a,b,c,d) = csp.segmentCoeffs(j)
      for t in [0.17, 0.5, 0.83]:
        let u = t*h
        let x = x0 + float(j)*h + u
        tseg.assertAlmostEqual(a + b*u + c*u*u + d*u*u*u, csp.interpolate(x))
        tseg.assertAlmostEqual(b + 2.0*c*u + 3.0*d*u*u, csp.interpolateDy(x))

    # Check derivative continuity across the wrap.
    test.newTest("seam derivative").assertAlmostEqual(
      csp.interpolateDy(x0 + 1e-7), csp.interpolateDy(x0 + 2.0*PI - 1e-7),
      absTol=1e-6, relTol=1e-6)

    # Simpson must match the exact cubic integral.
    let tint = test.newTest("segment integral")
    var total = 0.0
    for j in 0..<n:
      let (a,b,c,d) = csp.segmentCoeffs(j)
      let segInt = a*h + 0.5*b*h*h + (c/3.0)*h*h*h + 0.25*d*h*h*h*h
      let hm = 0.5*h
      let yL = a
      let yMid = a + b*hm + c*hm*hm + d*hm*hm*hm
      let yR = a + b*h + c*h*h + d*h*h*h
      tint.assertAlmostEqual(segInt, (h/6.0)*(yL + 4.0*yMid + yR))
      total += segInt
    var refInt = 0.0
    const nf = 200000
    for i in 0..<nf: refInt += f(x0 + (float(i)+0.5)*2.0*PI/float(nf)).v
    refInt *= 2.0*PI/float(nf)
    test.newTest("period integral").assertAlmostEqual(total, refInt, absTol=1e-4, relTol=1e-4)

    let tapp = test.newTest("approximation")
    for t in [0.1, 1.3, 2.7, 4.5, 6.0]:
      let x = x0 + t
      tapp.assertAlmostEqual(csp.interpolate(x), f(x).v, absTol=1e-5, relTol=1e-5)
      tapp.assertAlmostEqual(csp.interpolateDy(x), f(x).d1, absTol=1e-3, relTol=1e-3)

  qexInit()
  let thetest = newQEXTest("CSpline")
  thetest.run(1,fun0)
  thetest.run(2,fun1)
  thetest.run(3,fun2)
  thetest.runPeriodic
  thetest.qexFinalize
