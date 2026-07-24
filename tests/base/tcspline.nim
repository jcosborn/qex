#RUNCMD $RUN1

import std/[algorithm, math]
from qex import qexInit
import algorithms/cspline
import utils/test

proc fun0(x:float):auto =
  (1.0+x, 1.0, 0.0, 0.0)

proc fun1(x:float):auto =
  ((1.0+x)*(2.0-x), 1.0-2.0*x, -2.0, 0.0)

proc fun2(x:float):auto =
  ((1.0+x)*(2.0-x)*(1.0-x), (3.0*x-4.0)*x-1.0, 6.0*x-4.0, 6.0)

proc testEstD(test:QEXTest, ord:int, dx,dy:array[3,float], actual:float) =
  let
    d = [estimateDerivative([dx[0]],[dy[0]]),
         estimateDerivative([dx[0],dx[1]],[dy[0],dy[1]]),
         estimateDerivative(dx,dy)]
    test = test.newTest("estimate derivative")
  for o in ord..3:
    test.assertAlmostEqual(d[o-1], actual)

proc testCSp(test:QEXTest, spline:CSpline[float], x:openArray[float],
             f:proc, checkValues=true) =
  var knots = @x
  knots.sort

  let tknot = test.newTest("knot values")
  for x in knots:
    tknot.assertAlmostEqual(spline.interpolate(x), f(x)[0])

  let tcont = test.newTest("continuous derivatives")
  for i in 1..<knots.len-1:
    let
      h = knots[i]-knots[i-1]
      lo = spline.segmentCoeffs(i-1)
      hi = spline.segmentCoeffs(i)
    tcont.assertAlmostEqual(lo.b+2.0*lo.c*h+3.0*lo.d*h*h, hi.b)
    tcont.assertAlmostEqual(2.0*lo.c+6.0*lo.d*h, 2.0*hi.c)

  if checkValues:
    for x in [knots[0], knots[^1], knots[0]+0.05, 0.0, knots[^1]-0.05]:
      let
        testp = test.newTest("x=" & $x, hidden=1)
        fx = f(x)
      testp.newTest("y", hidden=1).assertAlmostEqual(spline.interpolate(x), fx[0])
      testp.newTest("dy", hidden=1).assertAlmostEqual(spline.interpolateDy(x), fx[1])
      testp.newTest("d2y", hidden=1).assertAlmostEqual(spline.interpolateD2y(x), fx[2])

proc run(test:QEXTest, ord:int, f:proc) =
  let
    test = test.newTest("polynomial degree " & $ord)
    n = 7
    m = 4
  var
    xs = newseq[float](n+m)
    ys = newseq[float](n+m)
    dys = newseq[array[3,float]](n+m)
  for i in 0..<n:
    let
      x = float(i)*5.0/float(n-1)-2.0
      fx = f(x)
    xs[i] = x
    ys[i] = fx[0]
    dys[i] = [fx[1],fx[2],fx[3]]
  for i in 0..<m:
    let
      x = float(i)*5.0/float(m-1)-1.9
      fx = f(x)
    xs[n+i] = x
    ys[n+i] = fx[0]
    dys[n+i] = [fx[1],fx[2],fx[3]]
  testEstD(test, ord,
    [xs[3]-xs[2],xs[4]-xs[2],xs[5]-xs[2]],
    [ys[3]-ys[2],ys[4]-ys[2],ys[5]-ys[2]], dys[2][0])
  testCSp(
    test.newTest("cspline default (est. 1st deriv.)"),
    newCSpline(xs,ys), xs, f)
  testCSp(
    test.newTest("cspline set 1st deriv. bounds"),
    newCSpline(xs,ys,csplineBounds(dyBound(dys[0][0]),dyBound(dys[^1][0]))),
    xs, f)
  testCSp(
    test.newTest("cspline natural (zero 2nd deriv.)"),
    newCSpline(xs,ys,csplineBounds(CSplineBoundZeroD2y,CSplineBoundZeroD2y)),
    xs, f, checkValues=false)

proc runPeriodic(test:QEXTest) =
  let test = test.newTest("periodic cubic spline")
  const n = 64
  let
    x0 = -PI
    h = 2.0*PI/float(n)
  proc f(p:float):tuple[v,d1:float] =
    let v = exp(0.5*cos(p)+0.3*sin(p))
    (v, (-0.5*sin(p)+0.3*cos(p))*v)
  var ys = newseq[float](n)
  for j in 0..<n:
    ys[j] = f(x0+float(j)*h).v
  let csp = newCSplinePeriodic(ys, x0, h)

  let tknot = test.newTest("knot values")
  for j in 0..<n:
    tknot.assertAlmostEqual(csp.interpolate(x0+float(j)*h), ys[j])

  let tseg = test.newTest("segment coeffs vs interpolate")
  for j in 0..<n:
    let p = csp.segmentCoeffs(j)
    for t in [0.17, 0.5, 0.83]:
      let
        u = t*h
        x = x0+float(j)*h+u
      tseg.assertAlmostEqual(p.a+u*(p.b+u*(p.c+u*p.d)), csp.interpolate(x))
      tseg.assertAlmostEqual(p.b+u*(2.0*p.c+3.0*u*p.d), csp.interpolateDy(x))

  test.newTest("seam derivative").assertAlmostEqual(
    csp.interpolateDy(x0+1e-7), csp.interpolateDy(x0+2.0*PI-1e-7),
    absTol=1e-6, relTol=1e-6)

  let tint = test.newTest("segment integral")
  var total = 0.0
  for j in 0..<n:
    let
      p = csp.segmentCoeffs(j)
      segInt = h*(p.a+h*(0.5*p.b+h*(p.c/3.0+0.25*h*p.d)))
      hm = 0.5*h
      ymid = p.a+hm*(p.b+hm*(p.c+hm*p.d))
      yr = p.a+h*(p.b+h*(p.c+h*p.d))
    tint.assertAlmostEqual(segInt, (h/6.0)*(p.a+4.0*ymid+yr))
    total += segInt
  var refInt = 0.0
  const nf = 200000
  for i in 0..<nf:
    refInt += f(x0+(float(i)+0.5)*2.0*PI/float(nf)).v
  refInt *= 2.0*PI/float(nf)
  test.newTest("period integral").assertAlmostEqual(
    total, refInt, absTol=1e-4, relTol=1e-4)

  let tapp = test.newTest("approximation")
  for t in [0.1, 1.3, 2.7, 4.5, 6.0]:
    let x = x0+t
    tapp.assertAlmostEqual(csp.interpolate(x), f(x).v, absTol=1e-5, relTol=1e-5)
    tapp.assertAlmostEqual(csp.interpolateDy(x), f(x).d1, absTol=1e-3, relTol=1e-3)

  let
    tb = test.newTest("periodic B-spline controls")
    ctl = periodicBSplineControls(ys)
  for j in 0..<n:
    let v = (ctl[(j+n-1) mod n]+4.0*ctl[j]+ctl[(j+1) mod n])/6.0
    tb.assertAlmostEqual(v, ys[j])

  let
    y3 = [0.7, -0.2, 1.1]
    ctl3 = periodicBSplineControls(y3)
  for j in 0..<3:
    let v = (ctl3[(j+2) mod 3]+4.0*ctl3[j]+ctl3[(j+1) mod 3])/6.0
    tb.assertAlmostEqual(v, y3[j])

  let
    tweight = test.newTest("physical B-spline weights")
    bc = [0.7, -0.2, 1.3, 0.4]
    bh = 1.7
    p = bSplineCoeffs(bc, bh)
  for u in [0.0, bh*1e-12, 0.17*bh, 0.5*bh, bh]:
    let w = bSplineWeights(u, bh)
    var v, dv, iv, sw, sdw, siw = 0.0
    for k in 0..3:
      v += w.w[k]*bc[k]
      dv += w.dw[k]*bc[k]
      iv += w.iw[k]*bc[k]
      sw += w.w[k]
      sdw += w.dw[k]
      siw += w.iw[k]
    tweight.assertAlmostEqual(sw, 1.0)
    tweight.assertAlmostEqual(sdw, 0.0)
    tweight.assertAlmostEqual(siw, u)
    tweight.assertAlmostEqual(v, p.a+u*(p.b+u*(p.c+u*p.d)))
    tweight.assertAlmostEqual(dv, p.b+u*(2.0*p.c+3.0*u*p.d))
    tweight.assertAlmostEqual(iv,
      u*(p.a+u*(0.5*p.b+u*(p.c/3.0+0.25*u*p.d))))
  let
    endw = bSplineWeights(bh, bh)
    full = [1.0, 11.0, 11.0, 1.0]
  for k in 0..3:
    tweight.assertAlmostEqual(endw.iw[k], bh*full[k]/24.0)
  let small = bSplineWeights(bh*1e-12, bh)
  tweight.assertAlmostEqual(if small.iw[0] > 0.0: 1.0 else: 0.0, 1.0)
  tweight.assertAlmostEqual(small.iw[0]/(bh*1e-12/6.0), 1.0,
    absTol=2e-11, relTol=2e-11)

qexInit()
let thetest = newQEXTest("CSpline")
thetest.run(1,fun0)
thetest.run(2,fun1)
thetest.run(3,fun2)
thetest.runPeriodic
thetest.qexFinalize
