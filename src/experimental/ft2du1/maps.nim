## Analytic circle and contextual maps used by the finite-transformation HMC.

import math
import algorithms/cspline

const
  mapDerivativeFloor* = 1.5e-8
  mapDerivativeCeiling = 1.0 / mapDerivativeFloor

# Scalar circle map: y(x+2π)=y(x)+2π, y'(x)>0.
type
  CircleEval* = object
    y*, dy*, ddy*, dddy*: float

  CircleStageKind = enum
    cskFourier, cskSpline, cskBSpline

  CircleStage = object
    kind: CircleStageKind
    sa, ca: seq[float]
    a, b, c, d, prefix: seq[float]
    h: float

  CircleMap* = object
    stages: seq[CircleStage]

proc circleEval*(m: CircleMap; x: float): CircleEval

func principalAngle*(x: float): float {.inline.} =
  x - TAU*floor((x + PI)/TAU)

proc mapFail*(msg: string) {.noreturn.} =
  raise newException(ValueError, msg)

proc segmentMin(a, b, c, d, h: float): float =
  result = min(a, a + h*(b + h*(c + h*d)))
  let
    qa = 3.0*d
    qb = 2.0*c
  var roots: array[2, float]
  var n = 0
  if abs(qa) < 1e-15:
    if abs(qb) > 1e-15:
      roots[0] = -b/qb
      n = 1
  else:
    let disc = qb*qb - 4.0*qa*b
    if disc >= 0.0:
      let s = sqrt(disc)
      roots = [(-qb-s)/(2.0*qa), (-qb+s)/(2.0*qa)]
      n = 2
  for k in 0..<n:
    let u = roots[k]
    if u > 0.0 and u < h:
      result = min(result, a + u*(b + u*(c + u*d)))

proc polyStage(kind: CircleStageKind; a, b, c, d: seq[float]; h: float): CircleStage =
  let n = a.len
  if n < 1 or b.len != n or c.len != n or d.len != n:
    mapFail("bad circle-map polynomial arrays")
  result = CircleStage(kind: kind, a: a, b: b, c: c, d: d,
    prefix: newSeq[float](n), h: h)
  var total = 0.0
  for j in 0..<n:
    result.prefix[j] = total
    total += h*(a[j] + h*(0.5*b[j] + h*(c[j]/3.0 + 0.25*h*d[j])))
    if segmentMin(a[j], b[j], c[j], d[j], h) < mapDerivativeFloor:
      mapFail("circle-map derivative is below the safe floor on segment " & $j)
  if abs(total-TAU) > 5e-12:
    mapFail("circle-map derivative is not normalized")

proc identityCircle*(): CircleMap =
  CircleMap()

proc fourierCircle(sa, ca: openArray[float]; derivativeFloor = mapDerivativeFloor;
                   certify = true): CircleMap =
  if derivativeFloor < mapDerivativeFloor:
    mapFail("Fourier derivative floor is below the safe floor")
  let n = max(sa.len, ca.len)
  var
    s = newSeq[float](n)
    c = newSeq[float](n)
    bound = 0.0
    slopeBound = 0.0
  for j in 0..<n:
    if j < sa.len: s[j] = sa[j]
    if j < ca.len: c[j] = ca[j]
    let a = abs(s[j])+abs(c[j])
    bound += float(j+1)*a
    slopeBound += float((j+1)*(j+1))*a
  if certify and 1.0-bound < derivativeFloor:
    var
      points = max(64, 8*n)
      certified = false
    while points <= 1048576:
      var lo = 1.0
      for j in 0..<points:
        let x = -PI+TAU*(float(j)+0.5)/float(points)
        var gp = 1.0
        for k in 0..<n:
          let nk = float(k+1)
          gp += nk*(s[k]*cos(nk*x)-c[k]*sin(nk*x))
        lo = min(lo, gp)
      if lo-PI*slopeBound/float(points) >= derivativeFloor:
        certified = true
        break
      points *= 2
    if not certified:
      mapFail("Fourier displacement does not certify the requested derivative floor")
  if n == 0: return identityCircle()
  CircleMap(stages: @[CircleStage(kind: cskFourier, sa: s, ca: c)])

proc sineCircle*(coeffs: openArray[float]; strength = 1.0;
                 derivativeFloor = mapDerivativeFloor): CircleMap =
  if strength < 0.0 or strength > 1.0:
    mapFail("sine strength must lie in [0,1]")
  var s = newSeq[float](coeffs.len)
  for j, v in coeffs: s[j] = strength*v
  fourierCircle(s, [], derivativeFloor)

proc cosineProduct(a, b: openArray[float]): seq[float] =
  result = newSeq[float](a.len+b.len-1)
  for i, x in a:
    for j, y in b:
      let v = x*y
      if i == 0 and j == 0: result[0] += v
      elif i == 0: result[j] += v
      elif j == 0: result[i] += v
      else:
        result[abs(i-j)] += 0.5*v
        result[i+j] += 0.5*v

proc sqFourierCircle*(coeffs: openArray[float]; epsilon = 0.05;
                      strength = 1.0; derivativeFloor = mapDerivativeFloor): CircleMap =
  if epsilon <= 0.0:
    mapFail("squared-Fourier epsilon must be positive")
  if strength < 0.0 or strength > 1.0:
    mapFail("squared-Fourier strength must lie in [0,1]")
  var q = newSeq[float](coeffs.len+1)
  q[0] = 1.0
  for j, v in coeffs: q[j+1] = v
  var rho = cosineProduct(q, q)
  rho[0] += epsilon
  if rho[0] <= 0.0:
    mapFail("squared-Fourier normalization is nonpositive")
  let lower = 1.0-strength+strength*epsilon/rho[0]
  var sa = newSeq[float](rho.len-1)
  for n in 1..<rho.len:
    sa[n-1] = strength*rho[n]/(rho[0]*float(n))
  fourierCircle(sa, [], derivativeFloor, certify = lower < derivativeFloor)

proc cSplineCircle*(logRho: openArray[float]; strength = 1.0;
                    derivativeFloor = mapDerivativeFloor): CircleMap =
  let n = logRho.len
  if n < 6 or (n and 1) != 0:
    mapFail("cubic spline needs an even number of knots of at least six")
  if strength < 0.0 or strength > 1.0:
    mapFail("cubic-spline strength must lie in [0,1]")
  var mx = logRho[0]
  for v in logRho: mx = max(mx, v)
  var y = newSeq[float](n)
  for j in 0..<n: y[j] = exp(logRho[j]-mx)
  let
    h = TAU/float(n)
    s = newCSplinePeriodic(y, -PI, h)
  var
    a = newSeq[float](n)
    b = newSeq[float](n)
    c = newSeq[float](n)
    d = newSeq[float](n)
    total = 0.0
  for j in 0..<n:
    let p = s.segmentCoeffs(j)
    a[j] = p.a
    b[j] = p.b
    c[j] = p.c
    d[j] = p.d
    total += h*(p.a + h*(0.5*p.b + h*(p.c/3.0 + 0.25*h*p.d)))
  if total <= 0.0: mapFail("bad cubic-spline normalization")
  let z = TAU/total
  for j in 0..<n:
    a[j] = 1.0-strength+strength*z*a[j]
    b[j] = strength*z*b[j]
    c[j] = strength*z*c[j]
    d[j] = strength*z*d[j]
  let st = polyStage(cskSpline, a, b, c, d, h)
  for j in 0..<n:
    if segmentMin(st.a[j], st.b[j], st.c[j], st.d[j], h) < derivativeFloor:
      mapFail("cubic-spline derivative is below the requested floor")
  CircleMap(stages: @[st])

proc bSplineCircle*(logControls: openArray[float]; identityFloor = 0.05;
                    strength = 1.0; derivativeFloor = mapDerivativeFloor): CircleMap =
  let n = logControls.len
  if n < 4: mapFail("cubic B-spline needs at least four controls")
  if identityFloor <= 0.0 or identityFloor >= 1.0:
    mapFail("B-spline identity floor must lie in (0,1)")
  if strength < 0.0 or strength > 1.0:
    mapFail("B-spline strength must lie in [0,1]")
  var mx = logControls[0]
  for v in logControls: mx = max(mx, v)
  var ctl = newSeq[float](n)
  var mean = 0.0
  for j in 0..<n:
    ctl[j] = exp(logControls[j]-mx)
    mean += ctl[j]
  mean /= float(n)
  let
    h = TAU/float(n)
    w = strength*(1.0-identityFloor)
  var
    a = newSeq[float](n)
    b = newSeq[float](n)
    c = newSeq[float](n)
    d = newSeq[float](n)
    lo = ctl[0]
    hi = ctl[0]
  for j in 0..<n:
    let
      cm = ctl[(j+n-1) mod n]
      c0 = ctl[j]
      cp = ctl[(j+1) mod n]
      c2 = ctl[(j+2) mod n]
      p = bSplineCoeffs([cm,c0,cp,c2], h)
    a[j] = 1.0-w+w*p.a/mean
    b[j] = w*p.b/mean
    c[j] = w*p.c/mean
    d[j] = w*p.d/mean
    lo = min(lo, c0)
    hi = max(hi, c0)
  let
    lower = 1.0-w+w*lo/mean
    upper = 1.0-w+w*hi/mean
  if lower < derivativeFloor or upper > mapDerivativeCeiling:
    mapFail("B-spline derivative bounds are not numerically safe")
  CircleMap(stages: @[polyStage(cskBSpline, a, b, c, d, h)])

proc fejerValidate(orders: openArray[int]; centers, weights: openArray[float]):
    tuple[ncomp: int; total: float] =
  ## Validate broadcastable Fejer inputs; return the component count and weight sum.
  if orders.len < 1: mapFail("Fejer basis needs at least one kernel")
  let ncomp = max(orders.len, max(centers.len, weights.len))
  if orders.len notin [1, ncomp] or centers.len notin [0, 1, ncomp] or
      weights.len notin [1, ncomp]:
    mapFail("Fejer orders and weights must broadcast or match; centers may also be empty for a uniform grid")
  var total = 0.0
  for j in 0..<ncomp:
    let order = orders[if orders.len == 1: 0 else: j]
    if order < 2: mapFail("Fejer orders must be at least two")
    let w = weights[if weights.len == 1: 0 else: j]
    if w < 0.0: mapFail("Fejer weights must be nonnegative")
    total += w
  if total <= 0.0: mapFail("Fejer weights must have positive sum")
  (ncomp, total)

proc fejerCircle*(orders: openArray[int]; centers, weights: openArray[float];
                  strength = 0.8; derivativeFloor = mapDerivativeFloor): CircleMap =
  if strength < 0.0 or strength >= 1.0 or 1.0-strength < derivativeFloor:
    mapFail("Fejer strength must leave the requested identity floor")
  let (ncomp, total) = fejerValidate(orders, centers, weights)
  var
    sa: seq[float]
    ca: seq[float]
  for j in 0..<ncomp:
    let order = orders[if orders.len == 1: 0 else: j]
    if sa.len < order-1:
      sa.setLen(order-1)
      ca.setLen(order-1)
    let
      mu = if centers.len == 0:
          (if ncomp == 1: 0.0 else: -PI+TAU*float(j)/float(ncomp))
        else: centers[if centers.len == 1: 0 else: j]
      w = strength*weights[if weights.len == 1: 0 else: j]/total
    for n in 1..<order:
      let a = 2.0*w*(1.0-float(n)/float(order))/float(n)
      sa[n-1] += a*cos(float(n)*mu)
      ca[n-1] -= a*sin(float(n)*mu)
  fourierCircle(sa, ca, derivativeFloor, certify = false)

proc composeCircle*(maps: openArray[CircleMap]): CircleMap =
  if maps.len < 1: mapFail("circle-map composition is empty")
  for m in maps: result.stages.add m.stages

proc stageEval(s: CircleStage; x: float): CircleEval =
  case s.kind
  of cskFourier:
    result = CircleEval(y: x, dy: 1.0)
    for j in 0..<max(s.sa.len, s.ca.len):
      let
        n = float(j+1)
        a = if j < s.sa.len: s.sa[j] else: 0.0
        b = if j < s.ca.len: s.ca[j] else: 0.0
        nx = n*x
        sn = sin(nx)
        cn = cos(nx)
      result.y += a*sn+b*cn
      result.dy += n*(a*cn-b*sn)
      result.ddy -= n*n*(a*sn+b*cn)
      result.dddy += n*n*n*(-a*cn+b*sn)
  of cskSpline, cskBSpline:
    let
      xw = principalAngle(x)
      lift = x-xw
      n = s.a.len
      jf = floor((xw+PI)/s.h)
      j = max(0, min(n-1, int(jf)))
      u = xw-(-PI+float(j)*s.h)
      a = s.a[j]
      b = s.b[j]
      c = s.c[j]
      d = s.d[j]
    result.y = -PI+s.prefix[j]+u*(a+u*(0.5*b+u*(c/3.0+0.25*u*d)))+lift
    result.dy = a+u*(b+u*(c+u*d))
    result.ddy = b+u*(2.0*c+3.0*u*d)
    result.dddy = 2.0*c+6.0*u*d

proc circleEval*(m: CircleMap; x: float): CircleEval =
  result = CircleEval(y: x, dy: 1.0)
  for s in m.stages:
    let e = stageEval(s, result.y)
    result.dddy = e.dddy*result.dy*result.dy*result.dy+
      3.0*e.ddy*result.dy*result.ddy+e.dy*result.dddy
    result.ddy = e.ddy*result.dy*result.dy+e.dy*result.ddy
    result.dy = e.dy*result.dy
    result.y = e.y
  if result.dy < mapDerivativeFloor:
    mapFail("circle-map derivative is below the safe floor")

proc circleInv*(m: CircleMap; y: float; tol = 2e-14; maxIter = 80): float =
  if tol <= 0.0 or maxIter < 1: mapFail("invalid circle-map inverse controls")
  let
    y0 = circleEval(m, -PI).y
    winding = floor((y-y0)/TAU)
    target = y-winding*TAU
  var
    lo = -PI
    hi = PI
    x = lo+(hi-lo)*(target-y0)/TAU
  x = max(lo, min(hi, x))
  for _ in 0..<maxIter:
    let
      e = circleEval(m, x)
      r = e.y-target
    if abs(r) <= tol*max(1.0, abs(target)):
      return x+winding*TAU
    if r < 0.0: lo = x else: hi = x
    let trial = x-r/e.dy
    if trial > lo and trial < hi: x = trial
    else: x = 0.5*(lo+hi)
  mapFail("circle-map inverse did not converge")

const maxMapContext* = 4

# Conditional scalar map y=g(x;c), including ∂x, ∂xx, ∂c, and ∂x∂c.
type
  ContextKind* = enum
    ckIdentity, ckSine, ckSqFourier, ckBSpline, ckTensorBSpline, ckFejer

  ContextEval* = object
    y*, dx*, dxx*: float
    dc*, dxc*: array[maxMapContext, float]

  ContextMap* = object
    kind*: ContextKind
    nctx*, order*, ctxOrder*, knots*, ctxKnots*: int
    strength*, floor*, epsilon*: float
    coeffs*, controls*, prefix*, means*, center*: seq[float]
    fejerOrders*: seq[int]
    fejerCenters*, fejerWeights*: seq[float]
    twisted*, phase*: bool

  PairCoreEval* = object
    y*, du*, dm*, duu*, dum*: float

  PairMap* = object
    stages*: seq[ContextMap]
    invTol*: float
    invIter*: int

  PairEval* = object
    pplus*, pminus*, u*, q*, delta*: float
    deltaPlus*, deltaMinus*: float
    physicalPlus*, physicalMinus*: float
    jac*: array[2, array[2, float]]
    det*, detPlus*, detMinus*: float
    logdet*, logdetPlus*, logdetMinus*: float

proc featureCount*(nctx, order: int; phase = false): int {.inline.} =
  if phase: 1+2*order else: 1+nctx*order

# a(c)=a0+Σ_k Σ_r a[k,r] cos(r c_k), stored row-major.
proc featureValue(data: openArray[float]; row, nctx, order: int;
                  ctx: array[maxMapContext, float]; phase = false): tuple[v: float; d: array[maxMapContext, float]] =
  let nf = featureCount(nctx, order, phase)
  result.v = data[row*nf]
  if phase:
    var
      p = ctx[nctx-1]
      q = row*nf+1
    for k in 0..<nctx-1: p -= ctx[k]
    let
      cp = cos(p)
      sp = sin(p)
    var
      cr = 1.0
      sr = 0.0
    for r in 1..order:
      let
        rf = float(r)
        ac = data[q]
        asi = data[q+1]
        nr = cr*cp-sr*sp
        ni = sr*cp+cr*sp
        dp = rf*(-ac*ni+asi*nr)
      cr = nr
      sr = ni
      result.v += ac*cr+asi*sr
      for k in 0..<nctx-1: result.d[k] -= dp
      result.d[nctx-1] += dp
      q += 2
  else:
    var q = row*nf+1
    for k in 0..<nctx:
      let
        cp = cos(ctx[k])
        sp = sin(ctx[k])
      var
        cr = 1.0
        sr = 0.0
      for r in 1..order:
        let
          a = data[q]
          nr = cr*cp-sr*sp
          ni = sr*cp+cr*sp
        cr = nr
        sr = ni
        result.v += a*cr
        result.d[k] -= float(r)*a*sr
        inc q

proc validateContextCommon(nctx, ctxOrder: int; strength, floor: float) =
  if nctx < 0 or nctx > maxMapContext:
    mapFail("context dimension is outside the supported range")
  if ctxOrder < 0:
    mapFail("context order must be nonnegative")
  if strength < 0.0 or strength > 1.0:
    mapFail("map strength must lie in [0,1]")
  if floor < mapDerivativeFloor or floor >= 1.0:
    mapFail("map derivative floor is outside the safe range")

proc identityContextMap*(nctx: int): ContextMap =
  validateContextCommon(nctx, 0, 0.0, 0.5)
  ContextMap(kind: ckIdentity, nctx: nctx, floor: 1.0)

proc sineContextMap*(nctx, order, ctxOrder: int; coeffs: openArray[float];
                     strength = 1.0; derivativeFloor = 0.05;
                     twisted = false): ContextMap =
  validateContextCommon(nctx, ctxOrder, strength, derivativeFloor)
  if order < 1: mapFail("conditional sine basis needs a positive order")
  let nf = featureCount(nctx, ctxOrder)
  if coeffs.len != order*nf:
    mapFail("wrong conditional sine coefficient count: expected " & $(order*nf) &
      ", got " & $coeffs.len)
  result = ContextMap(kind: ckSine, nctx: nctx, order: order,
    ctxOrder: ctxOrder, strength: strength, floor: derivativeFloor,
    coeffs: newSeq[float](coeffs.len), twisted: twisted)
  for i, v in coeffs: result.coeffs[i] = v
  if twisted:
    if nctx != 1: mapFail("twisted conditional maps need exactly one context")
    for n in 1..order:
      if (n and 1) != 0 and abs(result.coeffs[(n-1)*nf]) > 1e-14:
        mapFail("twisted sine constant coefficients require an even active mode")
      for r in 1..ctxOrder:
        if ((n+r) and 1) != 0 and abs(result.coeffs[(n-1)*nf+r]) > 1e-14:
          mapFail("twisted sine coefficients require active+context mode parity to be even")
  var bound = 0.0
  for n in 1..order:
    var a = abs(result.coeffs[(n-1)*nf])
    for j in 1..<nf: a += abs(result.coeffs[(n-1)*nf+j])
    bound += float(n)*a
  if 1.0-strength*bound < derivativeFloor:
    mapFail("conditional sine coefficients do not certify the requested derivative floor")

proc sqFourierContextMap*(nctx, order, ctxOrder: int; coeffs: openArray[float];
                          epsilon = 0.05; strength = 1.0;
                          derivativeFloor = 0.05; twisted = false): ContextMap =
  validateContextCommon(nctx, ctxOrder, strength, derivativeFloor)
  if order < 1: mapFail("conditional squared-Fourier basis needs a positive order")
  if epsilon <= 0.0: mapFail("conditional squared-Fourier epsilon must be positive")
  let nf = featureCount(nctx, ctxOrder)
  if coeffs.len != order*nf:
    mapFail("wrong conditional squared-Fourier coefficient count: expected " &
      $(order*nf) & ", got " & $coeffs.len)
  result = ContextMap(kind: ckSqFourier, nctx: nctx, order: order,
    ctxOrder: ctxOrder, strength: strength, floor: derivativeFloor,
    epsilon: epsilon, coeffs: newSeq[float](coeffs.len), twisted: twisted)
  for i, v in coeffs: result.coeffs[i] = v
  if twisted:
    if nctx != 1: mapFail("twisted conditional maps need exactly one context")
    for n in 1..order:
      if (n and 1) != 0 and abs(result.coeffs[(n-1)*nf]) > 1e-14:
        mapFail("twisted squared-Fourier constants require an even active mode")
      for r in 1..ctxOrder:
        if ((n+r) and 1) != 0 and abs(result.coeffs[(n-1)*nf+r]) > 1e-14:
          mapFail("twisted squared-Fourier coefficients require active+context mode parity to be even")

func wrapIndex(i, n: int): int {.inline.} =
  ((i mod n)+n) mod n

proc bSplineContextMap*(nctx, ctxOrder, knots: int; controls: openArray[float];
                        strength = 1.0; derivativeFloor = 0.05;
                        twisted = false; phase = false): ContextMap =
  validateContextCommon(nctx, ctxOrder, strength, derivativeFloor)
  if phase and nctx < 2:
    mapFail("phase B-splines need at least two contexts")
  if phase and twisted:
    mapFail("phase and twisted B-spline symmetries cannot be combined")
  if knots < 4 or twisted and (knots and 1) != 0:
    mapFail("conditional B-spline needs at least four controls and an even count when twisted")
  let nf = featureCount(nctx, ctxOrder, phase)
  result = ContextMap(kind: ckBSpline, nctx: nctx, ctxOrder: ctxOrder,
    knots: knots, strength: strength, floor: derivativeFloor,
    controls: newSeq[float](knots*nf), twisted: twisted, phase: phase)
  if twisted:
    if nctx != 1: mapFail("twisted B-splines need exactly one context")
    let half = knots div 2
    if controls.len != half*nf:
      mapFail("wrong independent twisted B-spline control count: expected " &
        $(half*nf) & ", got " & $controls.len)
    for j in 0..<half:
      for f in 0..<nf:
        let v = controls[j*nf+f]
        result.controls[j*nf+f] = v
        let r = if f == 0: 0 else: f
        result.controls[(j+half)*nf+f] = (if (r and 1) == 0: v else: -v)
  else:
    if controls.len != knots*nf:
      mapFail("wrong conditional B-spline control count: expected " &
        $(knots*nf) & ", got " & $controls.len)
    for i, v in controls: result.controls[i] = v
  result.means = newSeq[float](nf)
  result.prefix = newSeq[float]((knots+1)*nf)
  let
    h = TAU/float(knots)
    iw = bSplineWeights(h, h).iw
  for j in 0..<knots:
    for f in 0..<nf:
      result.means[f] += result.controls[j*nf+f]/float(knots)
      result.prefix[(j+1)*nf+f] = result.prefix[j*nf+f]
      for k in 0..3:
        result.prefix[(j+1)*nf+f] +=
          iw[k]*result.controls[wrapIndex(j+k-1, knots)*nf+f]
  for j in 0..<knots:
    var
      base = result.controls[j*nf]
      bound = 0.0
      slope = 0.0
    if phase:
      for r in 1..ctxOrder:
        let
          ac = result.controls[j*nf+2*r-1]
          asi = result.controls[j*nf+2*r]
        bound += abs(ac)+abs(asi)
        slope += float(r)*(abs(ac)+abs(asi))
    else:
      for k in 0..<nctx:
        for r in 1..ctxOrder:
          let a = result.controls[j*nf+1+k*ctxOrder+r-1]
          bound += abs(a)
          slope += float(r)*abs(a)
    if base-bound <= 0.0:
      if not phase and nctx != 1:
        mapFail("conditional B-spline controls do not certify positivity")
      let points = max(1024, 64*ctxOrder)
      var lo = base
      for i in 0..<points:
        var c: array[maxMapContext, float]
        let p = -PI+TAU*(float(i)+0.5)/float(points)
        if phase: c[nctx-1] = p else: c[0] = p
        lo = min(lo, featureValue(result.controls, j, nctx, ctxOrder,
          c, phase).v)
      if lo-PI*slope/float(points) <= 0.0:
        mapFail("conditional B-spline density control is not positive")

proc fejerContextMap*(orders: openArray[int]; centers, weights: openArray[float];
                      strength = 0.8; derivativeFloor = 0.05): ContextMap =
  validateContextCommon(1, 0, strength, derivativeFloor)
  let ncomp = fejerValidate(orders, centers, weights).ncomp
  if 1.0-strength < derivativeFloor:
    mapFail("conditional Fejer strength does not leave the requested derivative floor")
  result = ContextMap(kind: ckFejer, nctx: 1, strength: strength,
    floor: derivativeFloor, twisted: true,
    fejerOrders: newSeq[int](ncomp),
    fejerCenters: newSeq[float](ncomp),
    fejerWeights: newSeq[float](ncomp))
  for j in 0..<ncomp:
    result.fejerOrders[j] = orders[if orders.len == 1: 0 else: j]
    result.fejerCenters[j] = if centers.len == 0:
        (if ncomp == 1: 0.0 else: -PI+TAU*float(j)/float(ncomp))
      else: centers[if centers.len == 1: 0 else: j]
    result.fejerWeights[j] = weights[if weights.len == 1: 0 else: j]

type
  SplineContextEval = object
    v, dv, integ: float
    dc, ic: array[maxMapContext, float]

  SplineFeatureEval = object
    v, dv, integ: float

  SplinePoint = object
    s: int
    w: BSplineWeights[float]

func splinePoint(x: float; n: int): SplinePoint =
  let
    h = TAU/float(n)
    t = (x+PI)/h
    s = min(n-1, max(0, int(floor(t))))
    u = h*min(1.0, max(0.0, t-float(s)))
  result.s = s
  result.w = bSplineWeights(u, h)

func sampleSplineFeature(m: ContextMap; s, f: int;
                         w: BSplineWeights[float]): SplineFeatureEval =
  let nf = featureCount(m.nctx, m.ctxOrder, m.phase)
  result.integ = m.prefix[s*nf+f]
  for k in 0..3:
    let a = m.controls[wrapIndex(s+k-1, m.knots)*nf+f]
    result.v += w.w[k]*a
    result.dv += w.dw[k]*a
    result.integ += w.iw[k]*a

proc sampleContextSpline(m: ContextMap; x: float;
                         ctx: array[maxMapContext, float]): SplineContextEval =
  let
    q = splinePoint(x, m.knots)
    a0 = sampleSplineFeature(m, q.s, 0, q.w)
  result.v = a0.v
  result.dv = a0.dv
  result.integ = a0.integ
  if m.phase:
    var p = ctx[m.nctx-1]
    for k in 0..<m.nctx-1: p -= ctx[k]
    let
      cp = cos(p)
      sp = sin(p)
    var
      cr = 1.0
      sr = 0.0
    for r in 1..m.ctxOrder:
      let
        ac = sampleSplineFeature(m, q.s, 2*r-1, q.w)
        asi = sampleSplineFeature(m, q.s, 2*r, q.w)
        nr = cr*cp-sr*sp
        ni = sr*cp+cr*sp
        rf = float(r)
        dc = rf*(-ac.v*ni+asi.v*nr)
        di = rf*(-ac.integ*ni+asi.integ*nr)
      cr = nr
      sr = ni
      result.v += ac.v*cr+asi.v*sr
      result.dv += ac.dv*cr+asi.dv*sr
      result.integ += ac.integ*cr+asi.integ*sr
      for k in 0..<m.nctx-1:
        result.dc[k] -= dc
        result.ic[k] -= di
      result.dc[m.nctx-1] += dc
      result.ic[m.nctx-1] += di
  else:
    for k in 0..<m.nctx:
      let
        cp = cos(ctx[k])
        sp = sin(ctx[k])
      var
        cr = 1.0
        sr = 0.0
      for r in 1..m.ctxOrder:
        let
          a = sampleSplineFeature(m, q.s, 1+k*m.ctxOrder+r-1, q.w)
          nr = cr*cp-sr*sp
          ni = sr*cp+cr*sp
          rf = float(r)
        cr = nr
        sr = ni
        result.v += a.v*cr
        result.dv += a.dv*cr
        result.integ += a.integ*cr
        result.dc[k] -= rf*a.v*sr
        result.ic[k] -= rf*a.integ*sr

func contextAngle(m: ContextMap;
                  ctx: array[maxMapContext, float]): float =
  if m.phase:
    result = ctx[m.nctx-1]
    for k in 0..<m.nctx-1: result -= ctx[k]
  else:
    result = ctx[0]

func samplePeriodicRow(data: openArray[float]; n: int;
                       q: SplinePoint): tuple[v, d: float] =
  for k in 0..3:
    let a = data[wrapIndex(q.s+k-1, n)]
    result.v += q.w.w[k]*a
    result.d += q.w.dw[k]*a

proc setContextDerivative(m: ContextMap; d: float;
                          dst: var array[maxMapContext, float]) =
  if m.phase:
    for k in 0..<m.nctx-1: dst[k] = -d
    dst[m.nctx-1] = d
  else:
    dst[0] = d

proc tensorBSplineContextMap*(nctx, knots, ctxKnots: int;
                              controls: openArray[float];
                              strength = 1.0; derivativeFloor = 0.05;
                              twisted = false; phase = false): ContextMap =
  validateContextCommon(nctx, 0, strength, derivativeFloor)
  if nctx != 1 and not phase:
    mapFail("tensor B-splines need one context or one phase combination")
  if phase and nctx < 2:
    mapFail("phase tensor B-splines need at least two contexts")
  if twisted and (nctx != 1 or phase):
    mapFail("twisted tensor B-splines need exactly one context")
  if knots < 4 or ctxKnots < 4:
    mapFail("tensor B-splines need at least four knots in each direction")
  if twisted and ((knots and 1) != 0 or (ctxKnots and 1) != 0):
    mapFail("twisted tensor B-splines need even active and context knot counts")
  if controls.len != knots*ctxKnots:
    mapFail("wrong tensor B-spline control count: expected " &
      $(knots*ctxKnots) & ", got " & $controls.len)
  result = ContextMap(kind: ckTensorBSpline, nctx: nctx, knots: knots,
    ctxKnots: ctxKnots, strength: strength, floor: derivativeFloor,
    controls: newSeq[float](controls.len), twisted: twisted, phase: phase,
    prefix: newSeq[float]((knots+1)*ctxKnots),
    means: newSeq[float](ctxKnots), center: newSeq[float](ctxKnots))
  for i, v in controls:
    if v <= 0.0:
      mapFail("tensor B-spline density control is not positive")
    result.controls[i] = v
  if twisted:
    let
      hx = knots div 2
      hc = ctxKnots div 2
    for i in 0..<hx:
      for j in 0..<ctxKnots:
        let
          a = result.controls[(i+hx)*ctxKnots+j]
          b = result.controls[i*ctxKnots+wrapIndex(j+hc, ctxKnots)]
        if abs(a-b) > 2e-12*(1.0+abs(a)+abs(b)):
          mapFail("tensor B-spline controls violate twisted periodicity")
  let
    hx = TAU/float(knots)
    iw = bSplineWeights(hx, hx).iw
  for i in 0..<knots:
    for j in 0..<ctxKnots:
      result.means[j] += result.controls[i*ctxKnots+j]/float(knots)
      result.prefix[(i+1)*ctxKnots+j] =
        result.prefix[i*ctxKnots+j]
      for k in 0..3:
        result.prefix[(i+1)*ctxKnots+j] +=
          iw[k]*result.controls[wrapIndex(i+k-1, knots)*ctxKnots+j]
  if not phase:
    let q = splinePoint(0.0, knots)
    for j in 0..<ctxKnots:
      var z = result.prefix[q.s*ctxKnots+j]
      for k in 0..3:
        z += q.w.iw[k]*
          result.controls[wrapIndex(q.s+k-1, knots)*ctxKnots+j]
      result.center[j] = 0.5*z

proc sampleTensorSpline(m: ContextMap; x: float; cq: SplinePoint): SplineContextEval =
  let
    xq = splinePoint(x, m.knots)
    nc = m.ctxKnots
  var
    vc, ic: float
  for k in 0..3:
    let j = wrapIndex(cq.s+k-1, nc)
    let a = m.prefix[xq.s*nc+j]
    result.integ += cq.w.w[k]*a
    ic += cq.w.dw[k]*a
  for i in 0..3:
    let ix = wrapIndex(xq.s+i-1, m.knots)
    var a, ac: float
    for j in 0..3:
      let v = m.controls[ix*nc+wrapIndex(cq.s+j-1, nc)]
      a += cq.w.w[j]*v
      ac += cq.w.dw[j]*v
    result.v += xq.w.w[i]*a
    result.dv += xq.w.dw[i]*a
    result.integ += xq.w.iw[i]*a
    vc += xq.w.w[i]*ac
    ic += xq.w.iw[i]*ac
  setContextDerivative(m, vc, result.dc)
  setContextDerivative(m, ic, result.ic)

proc evalContext*(m: ContextMap; xRaw: float;
                  ctx: array[maxMapContext, float]): ContextEval =
  case m.kind
  of ckIdentity:
    result.y = xRaw
    result.dx = 1.0
  of ckSine:
    result.y = xRaw
    result.dx = 1.0
    let x = principalAngle(xRaw)
    for n in 1..m.order:
      let
        a = featureValue(m.coeffs, n-1, m.nctx, m.ctxOrder, ctx)
        nx = float(n)*x
        sn = sin(nx)
        cn = cos(nx)
        nf = float(n)
      result.y += m.strength*a.v*sn
      result.dx += m.strength*nf*a.v*cn
      result.dxx -= m.strength*nf*nf*a.v*sn
      for k in 0..<m.nctx:
        result.dc[k] += m.strength*a.d[k]*sn
        result.dxc[k] += m.strength*nf*a.d[k]*cn
  of ckSqFourier:
    let
      x = principalAngle(xRaw)
      nc = m.order+1
    var
      q = newSeq[float](nc)
      qc = newSeq[array[maxMapContext, float]](nc)
    q[0] = 1.0
    for n in 1..m.order:
      let a = featureValue(m.coeffs, n-1, m.nctx, m.ctxOrder, ctx)
      q[n] = a.v
      qc[n] = a.d
    var b = cosineProduct(q, q)
    b[0] += m.epsilon
    var bc: array[maxMapContext, seq[float]]
    for k in 0..<m.nctx:
      var qk = newSeq[float](nc)
      for n in 0..<nc: qk[n] = qc[n][k]
      bc[k] = cosineProduct(q, qk)
      for j in 0..<bc[k].len: bc[k][j] *= 2.0
    let
      b0 = b[0]
      blend = m.strength*(1.0-m.floor)
    if b0 <= 0.0: mapFail("conditional squared-Fourier normalization is nonpositive")
    result.y = xRaw
    result.dx = 1.0
    for n in 1..<b.len:
      let
        nf = float(n)
        sn = sin(nf*x)
        cn = cos(nf*x)
        a = b[n]/b0
      result.y += blend*a*sn/nf
      result.dx += blend*a*cn
      result.dxx -= blend*nf*a*sn
      for k in 0..<m.nctx:
        let da = (bc[k][n]*b0-b[n]*bc[k][0])/(b0*b0)
        result.dc[k] += blend*da*sn/nf
        result.dxc[k] += blend*da*cn
  of ckBSpline:
    let
      x = principalAngle(xRaw)
      lift = xRaw-x
      blend = m.strength*(1.0-m.floor)
    let
      sx = sampleContextSpline(m, x, ctx)
      mu = featureValue(m.means, 0, m.nctx, m.ctxOrder, ctx, m.phase)
      sz = if m.phase: SplineContextEval()
        else: sampleContextSpline(m, 0.0, ctx)
      integ = if m.phase: sx.integ else: sx.integ-0.5*sz.integ
      gx = if m.phase: -PI+integ/mu.v+lift
        else: integ/mu.v-0.5*PI+lift
      rhox = sx.v/mu.v
    if mu.v <= 0.0:
      mapFail("conditional B-spline normalization is nonpositive")
    result.y = xRaw+blend*(gx-xRaw)
    result.dx = 1.0+blend*(rhox-1.0)
    result.dxx = blend*sx.dv/mu.v
    for k in 0..<m.nctx:
      let
        ic = if m.phase: sx.ic[k] else: sx.ic[k]-0.5*sz.ic[k]
      result.dc[k] = blend*(ic*mu.v-integ*mu.d[k])/(mu.v*mu.v)
      result.dxc[k] = blend*(sx.dc[k]*mu.v-sx.v*mu.d[k])/(mu.v*mu.v)
  of ckTensorBSpline:
    let
      x = principalAngle(xRaw)
      lift = xRaw-x
      blend = m.strength*(1.0-m.floor)
      cq = splinePoint(principalAngle(contextAngle(m, ctx)), m.ctxKnots)
      sx = sampleTensorSpline(m, x, cq)
      mu = samplePeriodicRow(m.means, m.ctxKnots, cq)
      z = samplePeriodicRow(m.center, m.ctxKnots, cq)
      integ = if m.phase: sx.integ else: sx.integ-z.v
      gx = if m.phase: -PI+integ/mu.v+lift
        else: integ/mu.v-0.5*PI+lift
      rhox = sx.v/mu.v
    if mu.v <= 0.0:
      mapFail("tensor B-spline normalization is nonpositive")
    result.y = xRaw+blend*(gx-xRaw)
    result.dx = 1.0+blend*(rhox-1.0)
    result.dxx = blend*sx.dv/mu.v
    var md, zd: array[maxMapContext, float]
    setContextDerivative(m, mu.d, md)
    setContextDerivative(m, z.d, zd)
    for k in 0..<m.nctx:
      let ic = if m.phase: sx.ic[k] else: sx.ic[k]-zd[k]
      result.dc[k] = blend*(ic*mu.v-integ*md[k])/(mu.v*mu.v)
      result.dxc[k] = blend*(sx.dc[k]*mu.v-sx.v*md[k])/(mu.v*mu.v)
  of ckFejer:
    var total = 0.0
    for j in 0..<m.fejerOrders.len:
      total += m.fejerWeights[if m.fejerWeights.len == 1: 0 else: j]
    var h, f, fp = 0.0
    for j, order in m.fejerOrders:
      let
        w = m.fejerWeights[if m.fejerWeights.len == 1: 0 else: j]/total
        mu = m.fejerCenters[if m.fejerCenters.len == 1: 0 else: j]
        v = xRaw-ctx[0]-mu
      f += w
      for n in 1..<order:
        let
          nf = float(n)
          a = 2.0*w*(1.0-nf/float(order))
        f += a*cos(nf*v)
        fp -= nf*a*sin(nf*v)
        h += a*sin(nf*v)/nf
    result.y = xRaw+m.strength*h
    result.dx = 1.0+m.strength*(f-1.0)
    result.dxx = m.strength*fp
    result.dc[0] = -m.strength*(f-1.0)
    result.dxc[0] = -m.strength*fp
  if result.dx < m.floor and m.kind != ckIdentity:
    mapFail("conditional map derivative is below its requested floor")

proc evalContextStages*(stages: openArray[ContextMap]; x: float;
                        ctx: array[maxMapContext, float]): ContextEval =
  result = ContextEval(y: x, dx: 1.0)
  for s in stages:
    let e = evalContext(s, result.y, ctx)
    for k in 0..<s.nctx:
      result.dxc[k] = e.dxx*result.dx*result.dc[k]+e.dxc[k]*result.dx+
        e.dx*result.dxc[k]
      result.dc[k] = e.dx*result.dc[k]+e.dc[k]
    result.dxx = e.dxx*result.dx*result.dx+e.dx*result.dxx
    result.dx = e.dx*result.dx
    result.y = e.y

proc invertContext*(m: ContextMap; y: float; ctx: array[maxMapContext, float];
                    tol = 2e-14; maxIter = 80): float =
  if tol <= 0.0 or maxIter < 1: mapFail("invalid conditional inverse controls")
  let
    y0 = evalContext(m, -PI, ctx).y
    winding = floor((y-y0)/TAU)
    target = y-winding*TAU
  var
    lo = -PI
    hi = PI
    x = min(PI, max(-PI, target))
  for _ in 0..<maxIter:
    let
      e = evalContext(m, x, ctx)
      r = e.y-target
    if abs(r) <= tol*max(1.0, abs(target)):
      return x+winding*TAU
    if r < 0.0: lo = x else: hi = x
    let xn = x-r/e.dx
    if xn > lo and xn < hi: x = xn
    else: x = 0.5*(lo+hi)
  mapFail("conditional map inverse did not converge")

proc evalPairCore*(m: PairMap; u, context: float): PairCoreEval =
  var c: array[maxMapContext, float]
  c[0] = context
  let e = evalContextStages(m.stages, u, c)
  result = PairCoreEval(y: e.y, du: e.dx, dm: e.dc[0],
    duu: e.dxx, dum: e.dxc[0])

proc invertPairCore(m: PairMap; q, context: float): float =
  result = q
  var c: array[maxMapContext, float]
  c[0] = context
  for j in countdown(m.stages.len-1, 0):
    result = invertContext(m.stages[j], result, c, m.invTol, m.invIter)

# m=(p+ + p-)/2, u=(p+ - p-)/2, (P+,P-)=(m+g(u;m),m-g(u;m)).
proc evalPair*(m: PairMap; pplus, pminus: float): PairEval =
  let
    mean = 0.5*(pplus+pminus)
    u = 0.5*(pplus-pminus)
    e = evalPairCore(m, u, mean)
    ym = e.dm
  result.pplus = pplus
  result.pminus = pminus
  result.u = u
  result.q = e.y
  result.delta = e.y-u
  result.deltaPlus = 0.5*(e.du+ym-1.0)
  result.deltaMinus = 0.5*(1.0+ym-e.du)
  result.physicalPlus = mean+e.y
  result.physicalMinus = mean-e.y
  result.jac = [[0.5*(1.0+ym+e.du), 0.5*(1.0+ym-e.du)],
                [0.5*(1.0-ym-e.du), 0.5*(1.0-ym+e.du)]]
  result.det = e.du
  result.detPlus = 0.5*(e.dum+e.duu)
  result.detMinus = 0.5*(e.dum-e.duu)
  result.logdet = ln(result.det)
  result.logdetPlus = result.detPlus/result.det
  result.logdetMinus = result.detMinus/result.det

proc invertPair*(m: PairMap; physicalPlus, physicalMinus: float): array[2, float] =
  let
    mean = 0.5*(physicalPlus+physicalMinus)
    q = 0.5*(physicalPlus-physicalMinus)
    u = invertPairCore(m, q, mean)
  [mean+u, mean-u]

# Plaq4: transform one plaquette and distribute its change over four links.
