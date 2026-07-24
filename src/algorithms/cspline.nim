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

type
  CSpline*[T] = object
    x: seq[T]
    ys: seq[array[2,T]]  ## y and computed second derivatives of y

  CubicCoeffs*[T] = tuple[a,b,c,d:T]

  BSplineWeights*[T] = object
    w*, dw*, iw*: array[4,T]

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
    qexError("different length in x and y: ", n, " != ", y.len)
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

func bSplineCoeffs*[T](ctl:array[4,T]; h:T):CubicCoeffs[T] =
  ## ctl = [c[i-1],c[i],c[i+1],c[i+2]], h > 0.
  ## p(u) = a+b*u+c*u^2+d*u^3 on 0 <= u <= h.
  let
    hi = T(1.0)/h
    hi2 = hi*hi
    hi3 = hi2*hi
  result.a = (ctl[0]+T(4.0)*ctl[1]+ctl[2])/T(6.0)
  result.b = (-ctl[0]+ctl[2])*hi/T(2.0)
  result.c = (ctl[0]-T(2.0)*ctl[1]+ctl[2])*hi2/T(2.0)
  result.d = (-ctl[0]+T(3.0)*ctl[1]-T(3.0)*ctl[2]+ctl[3])*hi3/T(6.0)

func bSplineWeights*[T](u,h:T):BSplineWeights[T] =
  ## h > 0 and 0 <= u <= h.
  ## w[k] = B_k(u/h), dw[k] = dB_k(u/h)/du,
  ## iw[k] = integral_0^u B_k(v/h) dv.
  let
    f = u/h
    f2 = f*f
    f3 = f2*f
    fm = T(1.0)-f
    fm2 = fm*fm
    hi = T(1.0)/h
    hs = h/T(24.0)
  result.w = [fm2*fm/T(6.0),
    (T(3.0)*f3-T(6.0)*f2+T(4.0))/T(6.0),
    (-T(3.0)*f3+T(3.0)*f2+T(3.0)*f+T(1.0))/T(6.0),
    f3/T(6.0)]
  result.dw = [-T(0.5)*fm2*hi,
    (T(1.5)*f2-T(2.0)*f)*hi,
    (-T(1.5)*f2+f+T(0.5))*hi,
    T(0.5)*f2*hi]
  result.iw = [
    hs*f*(T(4.0)+f*(-T(6.0)+f*(T(4.0)-f))),
    hs*f*(T(16.0)+f2*(-T(8.0)+T(3.0)*f)),
    hs*f*(T(4.0)+f*(T(6.0)+f*(T(4.0)-T(3.0)*f))),
    hs*f*f3]

proc periodicBSplineControls*[T](y:openarray[T]):seq[T] =
  ## Controls c_j of the uniform periodic cubic B-spline interpolating y_j:
  ## (c_{j-1} + 4 c_j + c_{j+1}) / 6 = y_j.
  let n = y.len
  if n < 3:
    qexError("periodic cubic B-spline needs at least 3 values, got ", n)
  var
    a = newseq[T](n)
    b = newseq[T](n)
    c = newseq[T](n)
    r = newseq[T](n)
  for j in 0..<n:
    a[j] = T(1.0)
    b[j] = T(4.0)
    c[j] = T(1.0)
    r[j] = T(6.0)*y[j]
  cyclicSolve(a, b, c, T(1.0), T(1.0), r)

proc newCSplinePeriodic*[T](y:openarray[T]; x0,h:T):CSpline[T] =
  ## Uniform periodic cubic spline. Values wrap, so y[n] == y[0].
  let n = y.len
  if n < 3:
    qexError("periodic cubic spline needs at least 3 knots, got ", n)
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

func segmentCoeffs*[T](csp:CSpline[T]; i:int):CubicCoeffs[T] =
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
