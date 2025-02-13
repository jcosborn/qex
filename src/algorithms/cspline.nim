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
converter toCSplineBoundDy*(dy:float):CSplineBoundDy = CSplineBoundDy(kind:CSBSetDy,dy:dy)

proc csplineBounds*(lo=CSplineBoundEstimateDy, hi=CSplineBoundEstimateDy):CSplineBounds =
  CSplineBounds(lo:lo, hi:hi)

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
  var g = newseq[T](n-1)

  var beta:T

  if bounds.lo.kind==CSBZeroD2y:
    d2y(0) = T(0)
    g[0] = T(0)
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
    d2y(0) = T(3.0)*(d/h-dy)/h
    g[0] = T(0.5)

  # solve the tridiagonal system
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
    beta = bj - g[j-1]
    d2y(j) = (T(6.0)*(dhj-dhm)/hm-d2y(j-1))/beta
    g[j] = hjm/beta
    #echo "# iter j = ",j
    #echo "(h",j-1,"+h",j,")/3 ",(xp-xm)/3.0
    #echo "d",j-1,"/h",j-1," ",dhm,"  d",j,"/h",j," ",dhj
    #echo "2(1+h",j,"/h",j-1,") ",bj
    #echo "h",j,"/h",j-1," ",hjm
    #echo "6(d",j,"/h",j,"-d",j-1,"/h",j-1,")/h",j-1," ",T(6.0)*(dhj-dhm)/hm
    #echo "beta ",beta,"  g",j," ",g[j]
    #echo "forward: d2y(",j,") ",d2y(j)

  # last row
  if bounds.hi.kind==CSBZeroD2y:
    d2y(n-1) = T(0)
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
    d2y(n-1) = (T(6.0)*(dy-d/h)/h-d2y(n-2))/(T(2.0)-g[n-2])

  # back substitute
  for j in countdown(n-2,0):
    d2y(j) -= g[j]*d2y(j+1)
    #echo "back: d2y(",j,") ",d2y(j)

  return r

func bisect*[T](xs:openarray[T], x:T): int =
  ## assuming ascending xs
  ## no boundary check
  var
    tot = xs.len
    mid = tot div 2
  while tot>1:
    if x<xs[mid]:
      tot = mid
      mid -= tot div 2
    else:
      tot -= mid
      mid += tot div 2
  mid

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
      newCSpline(xs,ys,csplineBounds(dys[0][0],dys[^1][0])),
      ord, f)
    testCSp(
      test.newTest("cspline natural (zero 2nd deriv.)"),
      newCSpline(xs,ys,csplineBounds(CSplineBoundZeroD2y,CSplineBoundZeroD2y)),
      ord, f, checkValues=false)

  qexInit()
  let thetest = newQEXTest("CSpline")
  thetest.run(1,fun0)
  thetest.run(2,fun1)
  thetest.run(3,fun2)
  thetest.qexFinalize
