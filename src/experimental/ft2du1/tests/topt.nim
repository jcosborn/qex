const optTest = true
include ../opt

import maths/special

proc analyticCdf(x: float): float =
  (x+PI+0.3*sin(x)-(0.2/3.0)*(cos(3.0*x)+1.0))/TAU

proc cdfTests() =
  let smooth: LogProc = proc(x: float): float =
    ln(1.0+0.3*cos(x)+0.2*sin(3.0*x))
  for np in [1, 2, 4, 8, 16]:
    let c = buildCdf(smooth, np, 1e-15, 20)
    var cdfErr, invErr = 0.0
    for i in 0..256:
      let
        x = -PI+TAU*float(i)/256.0
        q = cdfPos(c, x)
        p = if q.upper: 1.0-q.p else: q.p
      cdfErr = max(cdfErr, abs(p-analyticCdf(x)))
      invErr = max(invErr, abs(inverseCdf(c, q)-x))
    let zErr = abs(c.logz-ln(TAU))
    doAssert zErr < 2e-15
    doAssert cdfErr < 3e-15
    doAssert invErr < 8e-15
    echo "cdf analytic panels=", np, " logZ=", zErr,
      " cdf=", cdfErr, " inverse=", invErr

  for k in [1.0, 6.0, 12.0, 24.0]:
    let
      kk = k
      peaked: LogProc = proc(x: float): float = kk*cos(x)
      c = buildCdf(peaked, 4, 1e-15, 20)
    var invErr, symErr = 0.0
    for i in 0..256:
      let
        x = -PI+TAU*float(i)/256.0
        q = cdfPos(c, x)
        p = if q.upper: 1.0-q.p else: q.p
      invErr = max(invErr, abs(inverseCdf(c, q)-x))
      symErr = max(symErr, abs(cdfAt(c, -x)+p-1.0))
    let zErr = abs(c.logz-ln(TAU*besselI0(k)))
    echo "cdf peaked k=", k, " logZ=", zErr,
      " symmetry=", symErr, " inverse=", invErr
    doAssert zErr < 8e-15
    doAssert invErr < 2e-13
    doAssert symErr < 5e-15

  let
    peaked: LogProc = proc(x: float): float = 96.0*cos(x)
    cref = buildCdf(peaked, 64, 1e-15, 24)
  var lastCdfErr: float
  for tol in [1e-6, 1e-9, 1e-12, 1e-15]:
    let c = buildCdf(peaked, 1, tol, 24)
    var cdfErr: float
    for i in 0..256:
      let x = -PI+TAU*float(i)/256.0
      cdfErr = max(cdfErr, abs(cdfAt(c, x)-cdfAt(cref, x)))
    let zErr = abs(c.logz-ln(TAU*besselI0(96.0)))
    echo "cdf tolerance tol=", tol, " panels=", c.panels.len,
      " logZ=", zErr, " cdf=", cdfErr
    lastCdfErr = cdfErr
  doAssert lastCdfErr < 3e-15

proc linkCdfConvergence() =
  let o = Opt(target: "boundary", beta: 6.0, bias: 3.0, kappa: 2.0)
  var zErr, yErr, gpErr: float
  var at: array[4, int]
  for stage in 0..<3:
    let
      t0 = 1.0-float(stage)/3.0
      t1 = 1.0-float(stage+1)/3.0
    for im in 0..<128:
      let
        m = -PI+TAU*float(im)/128.0
        src: LogProc = proc(x: float): float = linkLog(x, m, t0, o)
        dst: LogProc = proc(x: float): float = linkLog(x, m, t1, o)
        a1 = buildCdf(src, 1, 1e-15, 20, true)
        b1 = buildCdf(dst, 1, 1e-15, 20, true)
        a64 = buildCdf(src, 64, 1e-15, 20, true)
        b64 = buildCdf(dst, 64, 1e-15, 20, true)
      zErr = max(zErr, abs(a1.logz-a64.logz))
      zErr = max(zErr, abs(b1.logz-b64.logz))
      for j in 0..<512:
        let
          x = -PI+TAU*float(j)/512.0
          y1 = inverseCdf(b1, cdfPos(a1, x))
          y64 = inverseCdf(b64, cdfPos(a64, x))
          g1 = exp(src(x)-a1.logz-dst(y1)+b1.logz)
          g64 = exp(src(x)-a64.logz-dst(y64)+b64.logz)
        if abs(g1-g64) > gpErr:
          gpErr = abs(g1-g64)
          at = [stage, im, j, 0]
        yErr = max(yErr, abs(y1-y64))
  echo "link CDF panels 1/64 logZ=", zErr, " y=", yErr,
    " slope=", gpErr, " at=", at
  doAssert zErr < 3e-15
  doAssert yErr < 2e-13
  doAssert gpErr < 1e-12

proc blockCdfConvergence() =
  var o = Opt(target: "boundary", ctxBasis: "bspline",
    beta: 6.0, bias: 1.0, kappa: 0.5, floor: 1e-6,
    depth: 2, knots: 512, ctxKnots: 128,
    cdfPanels: 1, cdfTol: 1e-15, cdfDepth: 20)
  let a = blockFit(o)
  o.cdfPanels = 64
  let b = blockFit(o)
  doAssert a.controls.len == b.controls.len
  var rms, mx: float
  for i in 0..<a.controls.len:
    let d = a.controls[i]-b.controls[i]
    rms += d*d
    mx = max(mx, abs(d))
  rms = sqrt(rms/float(a.controls.len))
  echo "block CDF panels 1/64 controls=", a.controls.len,
    " rms=", rms, " max=", mx
  doAssert rms < 2e-13
  doAssert mx < 2e-12

proc peakConvergence() =
  var o = Opt(geometry: "link2", construction: "scalar",
    target: "boundary", ctxBasis: "bspline",
    beta: 6.0, bias: 3.0, kappa: 2.0, floor: 1e-6,
    scanStep: 2e-5, depth: 3, knots: 512, ctxKnots: 128,
    cdfPanels: 16, cdfTol: 1e-15, cdfDepth: 20, scan: 8)
  let
    fit = linkFit(o)
    spec = buildMapSpec(mapParams(o, fit.controls), o.beta)
  let a = landscapeMetrics(o, spec)
  o.scan = 96
  let b = landscapeMetrics(o, spec)
  echo "peak scan 8/96 stiffness=", a.stiffnessMax, "/", b.stiffnessMax
  doAssert abs(a.stiffnessMax-b.stiffnessMax) < 2e-4

proc validationTests() =
  var o = Opt(geometry: "link2", ctxBasis: "fourier",
    ctxOrder: 4, contextGrid: 8)
  var rejected = false
  try:
    o.validateContextGrid
  except ValueError:
    rejected = true
  doAssert rejected
  o.contextGrid = 9
  o.validateContextGrid
  o.ctxBasis = "bspline"
  o.contextGrid = 2
  o.validateContextGrid

cdfTests()
linkCdfConvergence()
blockCdfConvergence()
peakConvergence()
validationTests()
