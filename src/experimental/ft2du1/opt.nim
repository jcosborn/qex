## Dependency-free transport fitter for the finite-transformation HMC maps.

import base/globals
setDefaultNc(1)
setVLENmax(4)

import qex
import math
import std/[deques, heapqueue, os, strformat, strutils]
import algorithms/cspline
import scan

type
  Opt = object
    geometry, construction, target, output, scanOutput, ctxBasis: string
    beta, alpha, bias, kappa, power, gateStart, floor, scanStep, cdfTol: float
    depth, knots, ctxOrder, ctxKnots, cdfPanels, cdfDepth: int
    contextGrid, validation, scan: int

  LogProc = proc(x: float): float {.closure.}

  Quad = object
    val, err, round, resabs: float

  CdfPanel = object
    a, b, mass, err, round: float

  Cdf = object
    f: LogProc
    panels: seq[CdfPanel]
    cdf, surv: seq[float]
    lo, hi, shift, mass, logz, tol, err, round: float
    maxDepth: int
    centered: bool

  CdfPos = object
    p: float
    upper, neg: bool

  ShapeEval = object
    v, d1, d2: float

  FitResult = object
    controls: seq[float]
    residualRms, residualMax, minControl, minSlope, maxSlope: float
    cdfError, cdfRound, cdfResidual: float

  StageControls = tuple[controls: seq[float];
    minCtl, minSlope, maxSlope, cdfError, cdfRound, cdfResidual: float]

  Landscape = object
    potentialMin, potentialMax, minX, minY, maxX, maxY: float
    barrier, barrierFromMin, boundaryFree, boundaryFreeFromMin: float
    boundaryEffective, area, exitFraction, stiffnessRoughness: float
    forceMax, forceX, forceY, stiffnessMax, stiffnessX, stiffnessY: float

proc validateContextGrid(o: Opt) =
  if o.geometry != "plaq4" and o.ctxBasis == "fourier" and
      o.contextGrid < 2*o.ctxOrder+1:
    raise newException(ValueError,
      "Fourier context needs contextGrid >= 2*ctxOrder+1")

proc readOpt(): Opt =
  letParam:
    # plaq4 | link2 | block5
    geometry = "plaq4"
    # scalar | chain | coupling; inferred when empty
    construction = ""
    # reduced | exp | power | smooth | smoother | gate | boundary
    target = "power"
    # fitted parameter file
    output = "cvmap.txt"
    # optional local scan table
    scanOutput = ""
    # fourier | bspline
    ctxBasis = "fourier"
    beta = 6.0
    # reduced-beta fraction
    alpha = 0.75
    # high-action log-density amplitude D
    bias = 2.0
    # high-action localization
    kappa = 0.5
    power = 2.0
    gateStart = 0.35
    # identity floor removed before fitting controls
    floor = 0.000001
    # auxiliary finite-difference step for scan Hessians
    scanStep = 0.00002
    # internal transport stages
    depth = 1
    # active-coordinate B-spline knots
    knots = 512
    # Fourier context modes
    ctxOrder = 16
    # context B-spline knots
    ctxKnots = 128
    # Adaptive Gauss-Kronrod CDF integration
    cdfPanels = 16
    cdfTol = 1e-15
    cdfDepth = 16
    # Fourier projection needs contextGrid >= 2*ctxOrder+1.
    contextGrid = 128
    # target-fit validation points
    validation = 8192
    # local scan points per physical axis
    scan = 64
  result = Opt(
    geometry: geometry.toLowerAscii,
    construction: construction.toLowerAscii,
    target: target.toLowerAscii,
    output: output,
    scanOutput: scanOutput,
    ctxBasis: ctxBasis.toLowerAscii,
    beta: beta,
    alpha: alpha,
    bias: bias,
    kappa: kappa,
    power: power,
    gateStart: gateStart,
    floor: floor,
    scanStep: scanStep,
    depth: depth,
    knots: knots,
    ctxOrder: ctxOrder,
    ctxKnots: ctxKnots,
    cdfPanels: cdfPanels,
    cdfTol: cdfTol,
    cdfDepth: cdfDepth,
    contextGrid: contextGrid,
    validation: validation,
    scan: scan)
  if result.construction.len == 0:
    result.construction = if result.geometry == "block5": "coupling" else: "scalar"
  if result.geometry notin ["plaq4", "link2", "block5"]:
    raise newException(ValueError, "geometry must be plaq4, link2, or block5")
  if result.geometry == "block5" and result.construction notin ["chain", "coupling"]:
    raise newException(ValueError, "block5 construction must be chain or coupling")
  if result.geometry != "block5" and result.construction != "scalar":
    raise newException(ValueError, "plaq4 and link2 use construction=scalar")
  if result.target notin ["reduced", "exp", "power", "smooth", "smoother", "gate", "boundary"]:
    raise newException(ValueError,
      "target must be reduced, exp, power, smooth, smoother, gate, or boundary")
  if result.ctxBasis notin ["fourier", "bspline"]:
    raise newException(ValueError, "ctxBasis must be fourier or bspline")
  if result.beta <= 0.0 or result.alpha <= 0.0 or result.alpha > 1.0:
    raise newException(ValueError, "need beta>0 and 0<alpha<=1")
  if result.bias < 0.0 or result.kappa < 0.0 or result.power < 2.0:
    raise newException(ValueError, "need bias,kappa>=0 and power>=2")
  if result.gateStart < 0.0 or result.gateStart >= 1.0:
    raise newException(ValueError, "gateStart must lie in [0,1)")
  if result.floor <= 0.0 or result.floor >= 1.0:
    raise newException(ValueError, "floor must lie in (0,1)")
  if result.scanStep <= 0.0:
    raise newException(ValueError, "scanStep must be positive")
  if result.depth < 1 or result.knots < 4:
    raise newException(ValueError,
      "need depth>=1 and knots>=4")
  if result.geometry != "plaq4":
    if result.ctxBasis == "fourier" and result.ctxOrder < 0:
      raise newException(ValueError, "Fourier context needs ctxOrder>=0")
    if result.ctxBasis == "bspline" and result.ctxKnots < 4:
      raise newException(ValueError, "B-spline context needs ctxKnots>=4")
  if result.geometry in ["link2"] or
      result.geometry == "block5" and result.construction == "chain":
    if (result.knots and 1) != 0:
      raise newException(ValueError, "link2 and block5/chain need an even knot count")
    if result.ctxBasis == "bspline" and (result.ctxKnots and 1) != 0:
      raise newException(ValueError,
        "link2 and block5/chain need an even context knot count")
  if result.cdfPanels < 1 or result.cdfTol <= 0.0 or result.cdfDepth < 1:
    raise newException(ValueError, "need cdfPanels>=1, cdfTol>0, and cdfDepth>=1")
  if result.contextGrid < 2 or result.validation < 256:
    raise newException(ValueError, "increase contextGrid and validation")
  result.validateContextGrid

proc shape(z, dz, d2z: float; o: Opt): ShapeEval =
  var hz, hzz: float
  case o.target
  of "exp":
    result.v = exp(-2.0*o.kappa*(1.0-z))
    hz = 2.0*o.kappa*result.v
    hzz = 4.0*o.kappa*o.kappa*result.v
  of "power":
    result.v = pow(z, o.power)
    hz = o.power*pow(z, o.power-1.0)
    hzz = o.power*(o.power-1.0)*pow(z, o.power-2.0)
  of "smooth":
    result.v = z*z*(3.0-2.0*z)
    hz = 6.0*z*(1.0-z)
    hzz = 6.0-12.0*z
  of "smoother":
    result.v = z*z*z*(10.0+z*(-15.0+6.0*z))
    hz = 30.0*z*z*(1.0-z)*(1.0-z)
    hzz = 60.0*z*(1.0-z)*(1.0-2.0*z)
  of "gate":
    if z <= o.gateStart:
      return
    let
      w = 1.0-o.gateStart
      t = (z-o.gateStart)/w
    result.v = t*t*t*(10.0+t*(-15.0+6.0*t))
    hz = 30.0*t*t*(1.0-t)*(1.0-t)/w
    hzz = 60.0*t*(1.0-t)*(1.0-2.0*t)/(w*w)
  else:
    discard
  result.d1 = hz*dz
  result.d2 = hzz*dz*dz+hz*d2z

proc scalarLog(x, t: float; o: Opt): ShapeEval =
  if o.target == "reduced":
    let b = o.beta*(1.0-t*(1.0-o.alpha))
    return ShapeEval(v: b*cos(x), d1: -b*sin(x), d2: -b*cos(x))
  let
    z = 0.5*(1.0-cos(x))
    h = if o.target == "boundary":
        let
          v = exp(-o.kappa*(1.0+cos(x)))
          d1 = o.kappa*sin(x)*v
          d2 = (o.kappa*cos(x)+o.kappa*o.kappa*sin(x)*sin(x))*v
        ShapeEval(v: v, d1: d1, d2: d2)
      else: shape(z, 0.5*sin(x), 0.5*cos(x), o)
  ShapeEval(v: o.beta*cos(x)+t*o.bias*h.v,
    d1: -o.beta*sin(x)+t*o.bias*h.d1,
    d2: -o.beta*cos(x)+t*o.bias*h.d2)

proc linkLogEval(u, m, t: float; o: Opt): ShapeEval =
  # ell=2 beta cos(m) cos(u)+D h; V=-ell+C, F=-ell_u, K=-ell_uu.
  let
    cm = cos(m)
    a = 2.0*o.beta*cm*cos(u)
    au = -2.0*o.beta*cm*sin(u)
    auu = -2.0*o.beta*cm*cos(u)
  if o.target == "reduced":
    let s = 1.0-t*(1.0-o.alpha)
    return ShapeEval(v: s*a, d1: s*au, d2: s*auu)
  if o.target == "boundary":
    let
      pp = m+u
      pm = m-u
      r1 = exp(-o.kappa*(1.0+cos(pp)))
      r2 = exp(-o.kappa*(1.0+cos(pm)))
      h = r1+r2-r1*r2
      r1u = o.kappa*sin(pp)*r1
      r2u = -o.kappa*sin(pm)*r2
      r1uu = (o.kappa*cos(pp)+o.kappa*o.kappa*sin(pp)*sin(pp))*r1
      r2uu = (o.kappa*cos(pm)+o.kappa*o.kappa*sin(pm)*sin(pm))*r2
      hu = r1u*(1.0-r2)+r2u*(1.0-r1)
      huu = r1uu*(1.0-r2)+r2uu*(1.0-r1)-2.0*r1u*r2u
    return ShapeEval(v: a+t*o.bias*h, d1: au+t*o.bias*hu,
      d2: auu+t*o.bias*huu)
  let
    z = 0.5*(1.0-cm*cos(u))
    h = shape(z, 0.5*cm*sin(u), 0.5*cm*cos(u), o)
  ShapeEval(v: a+t*o.bias*h.v, d1: au+t*o.bias*h.d1,
    d2: auu+t*o.bias*h.d2)

proc linkLog(u, m, t: float; o: Opt): float =
  linkLogEval(u, m, t, o).v

proc blockLog(x, phase, t: float; o: Opt): float =
  let a = o.beta*(cos(x)+cos(phase-x))
  if o.target == "reduced":
    return (1.0-t*(1.0-o.alpha))*a
  if o.target == "boundary":
    let
      r1 = exp(-o.kappa*(1.0+cos(x)))
      r2 = exp(-o.kappa*(1.0+cos(phase-x)))
      h = r1+r2-r1*r2
    return a+t*o.bias*h
  let
    z = 0.25*(2.0-cos(x)-cos(phase-x))
    h = shape(z, 0.0, 0.0, o).v
  a+t*o.bias*h

const
  fpEps = 2.2204460492503130808e-16
  xgk: array[8, float] = [
    0.9914553711208126392, 0.9491079123427585245,
    0.8648644233597690728, 0.7415311855993944399,
    0.5860872354676911303, 0.4058451513773971669,
    0.2077849550078984676, 0.0]
  wgk: array[8, float] = [
    0.02293532201052922496, 0.06309209262997855329,
    0.1047900103222501838, 0.1406532597155259187,
    0.1690047266392679028, 0.1903505780647854099,
    0.2044329400752988924, 0.2094821410847278280]
  wg: array[4, float] = [
    0.1294849661688696933, 0.2797053914892766679,
    0.3818300505051189450, 0.4179591836734693878]

proc kahanAdd(s, c: var float; x: float) =
  let
    y = x-c
    t = s+y
  c = (t-s)-y
  s = t

proc gk15(f: LogProc; a, b, shift: float): Quad =
  let
    mid = 0.5*(a+b)
    half = 0.5*(b-a)
    fc = exp(f(mid)-shift)
  var
    fv1, fv2: array[7, float]
    rg = wg[3]*fc
    rk = wgk[7]*fc
    ra = wgk[7]*abs(fc)
  for j in 0..6:
    let d = half*xgk[j]
    fv1[j] = exp(f(mid-d)-shift)
    fv2[j] = exp(f(mid+d)-shift)
    let fs = fv1[j]+fv2[j]
    rk += wgk[j]*fs
    ra += wgk[j]*(abs(fv1[j])+abs(fv2[j]))
  for j in 0..2:
    rg += wg[j]*(fv1[2*j+1]+fv2[2*j+1])
  let mean = 0.5*rk
  var rc = wgk[7]*abs(fc-mean)
  for j in 0..6:
    rc += wgk[j]*(abs(fv1[j]-mean)+abs(fv2[j]-mean))
  result.val = rk*half
  result.resabs = ra*abs(half)
  rc *= abs(half)
  result.err = abs((rk-rg)*half)
  if rc != 0.0 and result.err != 0.0:
    result.err = rc*min(1.0, pow(200.0*result.err/rc, 1.5))
  result.round = 50.0*fpEps*result.resabs

proc quadSum(a, b: Quad): Quad =
  result.val = a.val+b.val
  result.err = a.err+b.err
  result.round = a.round+b.round
  result.resabs = a.resabs+b.resabs

func quadDone(q: Quad; rtol: float): bool =
  q.err <= max(rtol*q.resabs, q.round)

proc integrate(f: LogProc; a, b, shift, rtol: float; depth: int;
               q0: Quad): Quad =
  if q0.quadDone(rtol):
    return q0
  if depth <= 0:
    raise newException(ValueError,
      &"CDF quadrature did not converge on [{a:.8g},{b:.8g}]; increase -cdfDepth")
  let
    m = 0.5*(a+b)
    ql = gk15(f, a, m, shift)
    qr = gk15(f, m, b, shift)
  quadSum(
    integrate(f, a, m, shift, rtol, depth-1, ql),
    integrate(f, m, b, shift, rtol, depth-1, qr))

proc appendPanels(c: var Cdf; a, b, rtol: float; depth: int; q: Quad) =
  if q.quadDone(rtol):
    c.panels.add CdfPanel(a: a, b: b, mass: q.val,
      err: q.err, round: q.round)
    return
  if depth <= 0:
    raise newException(ValueError,
      &"CDF quadrature did not converge on [{a:.8g},{b:.8g}]; increase -cdfDepth")
  let
    m = 0.5*(a+b)
    ql = gk15(c.f, a, m, c.shift)
    qr = gk15(c.f, m, b, c.shift)
  c.appendPanels(a, m, rtol, depth-1, ql)
  c.appendPanels(m, b, rtol, depth-1, qr)

proc buildCdf(f: LogProc; panels: int; tol: float; depth: int;
              centered = false): Cdf =
  result.f = f
  result.tol = tol
  result.maxDepth = depth
  result.centered = centered
  result.lo = if centered: 0.0 else: -PI
  result.hi = PI
  let ns = max(256, 16*panels)
  result.shift = f(result.lo)
  for i in 1..ns:
    result.shift = max(result.shift,
      f(result.lo+(result.hi-result.lo)*float(i)/float(ns)))
  let h = (result.hi-result.lo)/float(panels)
  var seeds = newSeq[Quad](panels)
  for i in 0..<panels:
    seeds[i] = gk15(f, result.lo+float(i)*h,
      result.lo+float(i+1)*h,
      result.shift)
  for i in 0..<panels:
    result.appendPanels(result.lo+float(i)*h,
      result.lo+float(i+1)*h,
      tol, depth, seeds[i])
  result.cdf = newSeq[float](result.panels.len+1)
  result.surv = newSeq[float](result.panels.len+1)
  var cz: float
  for i, p in result.panels:
    kahanAdd(result.mass, cz, p.mass)
    result.cdf[i+1] = result.mass
    result.err += p.err
    result.round += p.round
  var
    sz, csz: float
    i = result.panels.len-1
  while i >= 0:
    kahanAdd(sz, csz, result.panels[i].mass)
    result.surv[i] = sz
    dec i
  result.logz = result.shift+ln(result.mass*(if centered: 2.0 else: 1.0))
  for i in 1..<result.cdf.len:
    result.cdf[i] /= result.mass
  for i in 0..<result.surv.len-1:
    result.surv[i] /= result.mass
  result.cdf[^1] = 1.0
  result.surv[0] = 1.0

proc lowerMass(c: Cdf; x: float): Quad =
  if x <= c.lo: return
  if x >= c.hi:
    return Quad(val: c.mass, err: c.err, round: c.round)
  let q = gk15(c.f, c.lo, x, c.shift)
  integrate(c.f, c.lo, x, c.shift, c.tol, c.maxDepth, q)

proc upperMass(c: Cdf; x: float): Quad =
  if x <= c.lo:
    return Quad(val: c.mass, err: c.err, round: c.round)
  if x >= c.hi: return
  let q = gk15(c.f, x, c.hi, c.shift)
  integrate(c.f, x, c.hi, c.shift, c.tol, c.maxDepth, q)

proc cdfPos(c: Cdf; x: float): CdfPos =
  let z =
    if c.centered:
      result.neg = x < 0.0
      min(PI, abs(x))
    else:
      if x <= -PI: return
      if x >= PI: return CdfPos(upper: true)
      x
  let lo = c.lowerMass(z).val
  if lo <= 0.5*c.mass:
    result.p = lo/c.mass
  else:
    result.p = c.upperMass(z).val/c.mass
    result.upper = true

proc cdfAt(c: Cdf; x: float): float =
  let q = c.cdfPos(x)
  let r = if q.upper: 1.0-q.p else: q.p
  if c.centered:
    if q.neg: 0.5*(1.0-r) else: 0.5*(1.0+r)
  else:
    r

proc inverseCdf(c: Cdf; q: CdfPos): float =
  if q.p <= 0.0:
    if c.centered:
      let x = if q.upper: PI else: 0.0
      return (if q.neg: -x else: x)
    return (if q.upper: PI else: -PI)
  var
    ilo = 0
    ihi = c.cdf.len-1
  if q.upper:
    while ihi-ilo > 1:
      let m = (ilo+ihi) div 2
      if c.surv[m] >= q.p: ilo = m
      else: ihi = m
  else:
    while ihi-ilo > 1:
      let m = (ilo+ihi) div 2
      if c.cdf[m] <= q.p: ilo = m
      else: ihi = m
  let panel = c.panels[ilo]
  let target = q.p*c.mass
  var
    lo = c.lo
    hi = c.hi
    lastR, lastTol: float
    x =
      if q.upper:
        panel.b-(panel.b-panel.a)*
          (q.p-c.surv[ilo+1])/(c.surv[ilo]-c.surv[ilo+1])
      else:
        panel.a+(panel.b-panel.a)*
          (q.p-c.cdf[ilo])/(c.cdf[ilo+1]-c.cdf[ilo])
  for _ in 0..<64:
    let
      m = if q.upper: c.upperMass(x) else: c.lowerMass(x)
      r = m.val-target
      d = exp(c.f(x)-c.shift)
      rtol = max(c.tol*target,
        4.0*fpEps*max(1.0, abs(x))*d)
      etol = max(rtol, m.err+m.round+q.p*(c.err+c.round))
    lastR = r
    lastTol = etol
    if q.upper:
      if r > 0.0: lo = x else: hi = x
    else:
      if r < 0.0: lo = x else: hi = x
    if abs(r) <= rtol or
        hi-lo <= 8.0*fpEps*max(1.0, abs(x)):
      return (if c.centered and q.neg: -x else: x)
    let xn = if q.upper: x+r/d else: x-r/d
    let next = if xn > lo and xn < hi: xn else: 0.5*(lo+hi)
    if next == x and abs(r) <= etol:
      return (if c.centered and q.neg: -x else: x)
    x = next
  if abs(lastR) <= lastTol:
    return (if c.centered and q.neg: -x else: x)
  raise newException(ValueError,
    &"CDF inverse did not converge: residual={lastR:.8g} tolerance={lastTol:.8g}")

proc mergeStage(acc: var FitResult; q: StageControls; firstSample: bool) =
  acc.cdfError = max(acc.cdfError, q.cdfError)
  acc.cdfRound = max(acc.cdfRound, q.cdfRound)
  acc.cdfResidual = max(acc.cdfResidual, q.cdfResidual)
  if firstSample:
    acc.minControl = q.minCtl
    acc.minSlope = q.minSlope
    acc.maxSlope = q.maxSlope
  else:
    acc.minControl = min(acc.minControl, q.minCtl)
    acc.minSlope = min(acc.minSlope, q.minSlope)
    acc.maxSlope = max(acc.maxSlope, q.maxSlope)

proc stageControls(src, dst: LogProc; o: Opt; centered: bool): StageControls =
  let
    a = buildCdf(src, o.cdfPanels, o.cdfTol, o.cdfDepth, centered)
    b = buildCdf(dst, o.cdfPanels, o.cdfTol, o.cdfDepth, centered)
  result.cdfError = max(a.err/a.mass, b.err/b.mass)
  result.cdfRound = max(a.round/a.mass, b.round/b.mass)
  var rho = newSeq[float](o.knots)
  for j in 0..<o.knots:
    let
      x = -PI+TAU*float(j)/float(o.knots)
      p = cdfPos(a, x)
      y = inverseCdf(b, p)
      gp = exp(src(x)-a.logz-dst(y)+b.logz)
      q = cdfPos(b, y)
      r =
        if q.upper == p.upper and q.neg == p.neg: abs(q.p-p.p)
        else:
          let
            qr = if q.upper: 1.0-q.p else: q.p
            pr = if p.upper: 1.0-p.p else: p.p
            qf = if centered and q.neg: -qr else: qr
            pf = if centered and p.neg: -pr else: pr
          abs(qf-pf)
    result.cdfResidual = max(result.cdfResidual, r)
    if j == 0:
      result.minSlope = gp
      result.maxSlope = gp
    else:
      result.minSlope = min(result.minSlope, gp)
      result.maxSlope = max(result.maxSlope, gp)
    rho[j] = (gp-o.floor)/(1.0-o.floor)
    if rho[j] <= 0.0:
      raise newException(ValueError,
        &"target slope {gp:.8g} is below floor {o.floor:.8g}; lower -floor")
  let ctl = periodicBSplineControls(rho)
  result.minCtl = ctl[0]
  for x in ctl: result.minCtl = min(result.minCtl, x)
  if result.minCtl <= 0.0:
    raise newException(ValueError,
      &"B-spline control {result.minCtl:.8g} is nonpositive; increase -knots or use more -depth")
  result.controls = ctl

proc scalarFit(o: Opt): FitResult =
  for stage in 0..<o.depth:
    let
      t0 = 1.0-float(stage)/float(o.depth)
      t1 = 1.0-float(stage+1)/float(o.depth)
      src: LogProc = proc(x: float): float = scalarLog(x, t0, o).v
      dst: LogProc = proc(x: float): float = scalarLog(x, t1, o).v
      q = stageControls(src, dst, o, true)
    for x in q.controls: result.controls.add ln(x)
    mergeStage(result, q, stage == 0)

func tensorNode(c: openArray[float]; row, col, nc: int): float =
  (c[row*nc+(col+nc-1) mod nc]+4.0*c[row*nc+col]+
    c[row*nc+(col+1) mod nc])/6.0

proc contextSplineControls(y: openArray[float]; what: string): seq[float] =
  result = periodicBSplineControls(y)
  var lo = result[0]
  for x in result: lo = min(lo, x)
  if lo <= 0.0:
    raise newException(ValueError,
      &"{what} context B-spline control {lo:.8g} is nonpositive; increase -ctxKnots")

proc linkFit(o: Opt): FitResult =
  o.validateContextGrid
  let
    nm = if o.ctxBasis == "bspline": o.ctxKnots else: o.contextGrid
    half = o.knots div 2
    nf = 1+o.ctxOrder
  var nerr = 0
  for stage in 0..<o.depth:
    let
      t0 = 1.0-float(stage)/float(o.depth)
      t1 = 1.0-float(stage+1)/float(o.depth)
    var rows = newSeq[seq[float]](nm)
    for im in 0..<nm:
      let m = -PI+TAU*float(im)/float(nm)
      let
        src: LogProc = proc(x: float): float = linkLog(x, m, t0, o)
        dst: LogProc = proc(x: float): float = linkLog(x, m, t1, o)
        q = stageControls(src, dst, o, true)
      rows[im] = q.controls
      mergeStage(result, q, stage == 0 and im == 0)
    var c: seq[float]
    if o.ctxBasis == "bspline":
      c = newSeq[float](o.knots*nm)
      for j in 0..<half:
        var y = newSeq[float](nm)
        for im in 0..<nm:
          y[im] = 0.5*(rows[im][j]+
            rows[(im+nm div 2) mod nm][j+half])
        let ctl = contextSplineControls(y, &"link stage {stage} active knot {j}")
        for im in 0..<nm:
          c[j*nm+im] = ctl[im]
          c[(j+half)*nm+im] = ctl[(im+nm div 2) mod nm]
          result.minControl = min(result.minControl, ctl[im])
      for im in 0..<nm:
        for j in 0..<o.knots:
          let e = tensorNode(c, j, im, nm)-rows[im][j]
          result.residualRms += e*e
          result.residualMax = max(result.residualMax, abs(e))
          inc nerr
    else:
      c = newSeq[float](half*nf)
      for j in 0..<half:
        for im in 0..<nm: c[j*nf] += rows[im][j]/float(nm)
        for r in 1..o.ctxOrder:
          for im in 0..<nm:
            let m = -PI+TAU*float(im)/float(nm)
            c[j*nf+r] += 2.0*rows[im][j]*cos(float(r)*m)/float(nm)
      for im in 0..<nm:
        let m = -PI+TAU*float(im)/float(nm)
        for j in 0..<o.knots:
          let
            jj = j mod half
            shift = if j < half: 1.0 else: -1.0
          var v = c[jj*nf]
          for r in 1..o.ctxOrder:
            v += c[jj*nf+r]*cos(float(r)*m)*
              (if j < half or (r and 1) == 0: 1.0 else: shift)
          let e = v-rows[im][j]
          result.residualRms += e*e
          result.residualMax = max(result.residualMax, abs(e))
          inc nerr
    result.controls.add c
  result.residualRms = sqrt(result.residualRms/float(nerr))

proc blockFit(o: Opt): FitResult =
  o.validateContextGrid
  let
    nc = if o.ctxBasis == "bspline": o.ctxKnots else: o.contextGrid
    nf = 1+2*o.ctxOrder
  var nerr = 0
  for stage in 0..<o.depth:
    let
      t0 = 1.0-float(stage)/float(o.depth)
      t1 = 1.0-float(stage+1)/float(o.depth)
    var rows = newSeq[seq[float]](nc)
    for ic in 0..<nc:
      let p = -PI+TAU*float(ic)/float(nc)
      let
        src: LogProc = proc(x: float): float = blockLog(x, p, t0, o)
        dst: LogProc = proc(x: float): float = blockLog(x, p, t1, o)
        q = stageControls(src, dst, o, false)
      rows[ic] = q.controls
      mergeStage(result, q, stage == 0 and ic == 0)
    var c: seq[float]
    if o.ctxBasis == "bspline":
      c = newSeq[float](o.knots*nc)
      for j in 0..<o.knots:
        var y = newSeq[float](nc)
        for ic in 0..<nc: y[ic] = rows[ic][j]
        let ctl = contextSplineControls(y, &"block stage {stage} active knot {j}")
        for ic in 0..<nc:
          c[j*nc+ic] = ctl[ic]
          result.minControl = min(result.minControl, ctl[ic])
      for ic in 0..<nc:
        for j in 0..<o.knots:
          let e = tensorNode(c, j, ic, nc)-rows[ic][j]
          result.residualRms += e*e
          result.residualMax = max(result.residualMax, abs(e))
          inc nerr
    else:
      c = newSeq[float](o.knots*nf)
      for j in 0..<o.knots:
        for ic in 0..<nc: c[j*nf] += rows[ic][j]/float(nc)
        for r in 1..o.ctxOrder:
          for ic in 0..<nc:
            let
              p = -PI+TAU*float(ic)/float(nc)
              w = 2.0*rows[ic][j]/float(nc)
            c[j*nf+2*r-1] += w*cos(float(r)*p)
            c[j*nf+2*r] += w*sin(float(r)*p)
      for ic in 0..<nc:
        let p = -PI+TAU*float(ic)/float(nc)
        for j in 0..<o.knots:
          var v = c[j*nf]
          for r in 1..o.ctxOrder:
            v += c[j*nf+2*r-1]*cos(float(r)*p)+
              c[j*nf+2*r]*sin(float(r)*p)
          let e = v-rows[ic][j]
          result.residualRms += e*e
          result.residualMax = max(result.residualMax, abs(e))
          inc nerr
    result.controls.add c
  result.residualRms = sqrt(result.residualRms/float(nerr))

proc mapParams(o: Opt; c: seq[float]): MapParams =
  MapParams(geometry: o.geometry, construction: o.construction,
    basis: "bspline", ctxBasis: o.ctxBasis,
    mapDepth: o.depth, flowDepth: 1,
    mapStrengths: @[1.0], mapFloor: o.floor, mapEpsilon: 0.05,
    mapOrder: 1, ctxOrder: o.ctxOrder, mapKnots: o.knots,
    ctxKnots: o.ctxKnots,
    mapControls: c, mapStageOrder: "0123", mapDirs: @[0],
    mapParities: @[0], mapOffsets: @[0, 0],
    mapStride: (if o.geometry == "block5": 4 else: 0),
    mapInvTol: 2e-14, mapScanStep: o.scanStep, mapInvIter: 100,
    mapScan: o.scan)

proc scalarMetrics(o: Opt; spec: MapSpec) =
  let
    e0 = circleEval(spec.circle, 0.0)
    t0 = scalarLog(0.0, 1.0, o)
    v0 = -o.beta*cos(e0.y)-ln(e0.dy)+t0.v
  var
    vr, fr, kr = 0.0
    vmax, fmax, kmax = 0.0
    mingp, maxgp: float
  for i in 0..<o.validation:
    let
      x = -PI+TAU*(float(i)+0.5)/float(o.validation)
      e = circleEval(spec.circle, x)
      t = scalarLog(x, 1.0, o)
      v = -o.beta*cos(e.y)-ln(e.dy)
      f = o.beta*sin(e.y)*e.dy-e.ddy/e.dy
      k = o.beta*(cos(e.y)*e.dy*e.dy+sin(e.y)*e.ddy)-
        e.dddy/e.dy+(e.ddy/e.dy)*(e.ddy/e.dy)
    if i == 0:
      mingp = e.dy
      maxgp = e.dy
    else:
      mingp = min(mingp, e.dy)
      maxgp = max(maxgp, e.dy)
    let
      ev = v+t.v-v0
      ef = f+t.d1
      ek = k+t.d2
    vr += ev*ev
    fr += ef*ef
    kr += ek*ek
    vmax = max(vmax, abs(ev))
    fmax = max(fmax, abs(ef))
    kmax = max(kmax, abs(ek))
  echo &"target fit V rms={sqrt(vr/float(o.validation)):.8g} max={vmax:.8g}"
  echo &"target fit F rms={sqrt(fr/float(o.validation)):.8g} max={fmax:.8g}"
  echo &"target fit K rms={sqrt(kr/float(o.validation)):.8g} max={kmax:.8g}"
  echo &"composed slope=[{mingp:.8g}, {maxgp:.8g}]"

proc linkMetrics(o: Opt; spec: MapSpec) =
  let
    nm = min(o.contextGrid, 64)
    nu = min(o.validation, 2048)
    h = o.scanStep
  var
    vr, fr, kr = 0.0
    vmax, fmax, kmax = 0.0
    kmaxU, kmaxM = 0.0
    n = 0
  for im in 0..<nm:
    let m = -PI+TAU*(float(im)+0.5)/float(nm)
    var p0: Vec5
    p0[0] = m
    p0[1] = m
    let
      e0 = evalLocal(spec, p0)
      t0 = linkLogEval(0.0, m, 1.0, o)
    for iu in 0..<nu:
      let u = -PI+TAU*(float(iu)+0.5)/float(nu)
      var
        p: Vec5
        pp: Vec5
        pm: Vec5
      p[0] = m+u
      p[1] = m-u
      pp = p
      pm = p
      pp[0] += h
      pp[1] -= h
      pm[0] -= h
      pm[1] += h
      let
        e = evalLocal(spec, p)
        t = linkLogEval(u, m, 1.0, o)
        v = e.seff-e0.seff
        f = e.force[0]
        k = (evalLocal(spec, pp).force[0]-evalLocal(spec, pm).force[0])/(2.0*h)
        ev = v+t.v-t0.v
        ef = f+t.d1
        ek = k+t.d2
      vr += ev*ev
      fr += ef*ef
      kr += ek*ek
      vmax = max(vmax, abs(ev))
      fmax = max(fmax, abs(ef))
      if abs(ek) > kmax:
        kmax = abs(ek)
        kmaxU = u
        kmaxM = m
      inc n
  echo &"target fit V rms={sqrt(vr/float(n)):.8g} max={vmax:.8g}"
  echo &"target fit F rms={sqrt(fr/float(n)):.8g} max={fmax:.8g}"
  echo &"target fit K rms={sqrt(kr/float(n)):.8g} max={kmax:.8g} " &
    &"at (u={kmaxU:.8g}, m={kmaxM:.8g})"

proc wrapAngle(x: float): float =
  result = x
  while result < -PI: result += TAU
  while result >= PI: result -= TAU

proc refinePeak(o: Opt; spec: MapSpec; x0, y0, h: float;
                stiffness: bool): tuple[v, x, y: float] =
  result.x = x0
  result.y = y0
  let p0 = scanPoint(spec, scanPhysical(spec, x0, y0))
  result.v = if stiffness: p0.stiffMax else: sqrt(p0.force2)
  var step = h
  for iteration in 0..<96:
    var
      bx = result.x
      by = result.y
      bv = result.v
    for dx in -1..1:
      for dy in -1..1:
        if dx == 0 and dy == 0: continue
        let
          x = wrapAngle(result.x+float(dx)*step)
          y = wrapAngle(result.y+float(dy)*step)
          p = scanPoint(spec, scanPhysical(spec, x, y))
          v = if stiffness: p.stiffMax else: sqrt(p.force2)
        if v > bv:
          bx = x
          by = y
          bv = v
    if bv > result.v:
      result = (bv, bx, by)
    else:
      step *= 0.5
      if step < 1e-7: break

proc landscapeMetrics(o: Opt; spec: MapSpec): Landscape =
  let
    n = o.scan
    h = TAU/float(n)
    src = (n div 2)*n+n div 2
    bins = min(6, n)
  var
    v = newSeq[float](n*n)
    k = newSeq[float](n*n)
    forceBins = newSeq[float](bins*bins)
    stiffBins = newSeq[float](bins*bins)
    forcePoints = newSeq[int](bins*bins)
    stiffPoints = newSeq[int](bins*bins)
    forceSeen = newSeq[bool](bins*bins)
    stiffSeen = newSeq[bool](bins*bins)
  for ix in 0..<n:
    let x = -PI+float(ix)*h
    for iy in 0..<n:
      let
        y = -PI+float(iy)*h
        p = scanPoint(spec, scanPhysical(spec, x, y))
        j = ix*n+iy
        b = ((ix*bins) div n)*bins+(iy*bins) div n
      v[j] = p.dseff
      k[j] = p.stiffMax
      if not forceSeen[b] or p.force2 > forceBins[b]:
        forceSeen[b] = true
        forceBins[b] = p.force2
        forcePoints[b] = j
      if not stiffSeen[b] or p.stiffMax > stiffBins[b]:
        stiffSeen[b] = true
        stiffBins[b] = p.stiffMax
        stiffPoints[b] = j
  for b in 0..<bins*bins:
    if forceSeen[b]:
      let
        j = forcePoints[b]
        x = -PI+float(j div n)*h
        y = -PI+float(j mod n)*h
      for step in [h, TAU/float(bins)]:
        let fp = refinePeak(o, spec, x, y, step, false)
        if fp.v > result.forceMax:
          result.forceMax = fp.v
          result.forceX = fp.x
          result.forceY = fp.y
    if stiffSeen[b]:
      let
        j = stiffPoints[b]
        x = -PI+float(j div n)*h
        y = -PI+float(j mod n)*h
      for step in [h, TAU/float(bins)]:
        let kp = refinePeak(o, spec, x, y, step, true)
        if kp.v > result.stiffnessMax:
          result.stiffnessMax = kp.v
          result.stiffnessX = kp.x
          result.stiffnessY = kp.y
  var
    vmin = v[0]
    vmax = v[0]
    imin = 0
    imax = 0
  for i in 1..<v.len:
    if v[i] < vmin:
      vmin = v[i]
      imin = i
    if v[i] > vmax:
      vmax = v[i]
      imax = i
  result.potentialMin = vmin
  result.potentialMax = vmax
  result.minX = -PI+float(imin div n)*h
  result.minY = -PI+float(imin mod n)*h
  result.maxX = -PI+float(imax div n)*h
  result.maxY = -PI+float(imax mod n)*h

  let far = vmax+abs(vmax-vmin)+1.0
  proc minimax(start: int): float =
    var dist = newSeq[float](n*n)
    for x in dist.mitems: x = far
    dist[start] = v[start]
    var q: HeapQueue[(float, int)]
    q.push((v[start], start))
    while q.len > 0:
      let (d, j) = q.pop()
      if d != dist[j]: continue
      let
        ix = j div n
        iy = j mod n
      if ix == 0 or o.geometry == "link2" and iy == 0:
        return d
      for dx in -1..1:
        for dy in -1..1:
          if dx == 0 and dy == 0: continue
          let
            x = ix+dx
            y = iy+dy
          if x >= 0 and x < n and y >= 0 and y < n:
            let
              z = x*n+y
              nd = max(d, v[z])
            if nd < dist[z]:
              dist[z] = nd
              q.push((nd, z))
    far

  let
    vacuumLevel = minimax(src)
    minLevel = minimax(imin)
  result.barrier = vacuumLevel-v[src]
  result.barrierFromMin = minLevel-vmin

  var boundary: seq[float]
  for ix in 0..<n:
    for iy in 0..<n:
      if ix == 0 or o.geometry == "link2" and iy == 0:
        boundary.add v[ix*n+iy]-v[src]
  var bm = boundary[0]
  for x in boundary: bm = min(bm, x)
  var sw, sw2 = 0.0
  for x in boundary:
    let w = exp(-(x-bm))
    sw += w
    sw2 += w*w
  result.boundaryFree = bm-ln(sw/float(boundary.len))
  result.boundaryFreeFromMin = result.boundaryFree+v[src]-vmin
  result.boundaryEffective = sw*sw/(float(boundary.len)*sw2)

  let level = vacuumLevel+0.5
  var
    seen = newSeq[bool](n*n)
    bfs = initDeque[int]()
    exits = 0
    boundaryCount = 0
  if v[src] <= level:
    seen[src] = true
    bfs.addLast(src)
  while bfs.len > 0:
    let
      j = bfs.popFirst()
      ix = j div n
      iy = j mod n
    for dx in -1..1:
      for dy in -1..1:
        if dx == 0 and dy == 0: continue
        let
          x = ix+dx
          y = iy+dy
        if x >= 0 and x < n and y >= 0 and y < n:
          let z = x*n+y
          if not seen[z] and v[z] <= level:
            seen[z] = true
            bfs.addLast(z)
  for ix in 0..<n:
    for iy in 0..<n:
      let boundaryPoint = ix == 0 or o.geometry == "link2" and iy == 0
      if seen[ix*n+iy]: result.area += 1.0
      if boundaryPoint:
        inc boundaryCount
        if seen[ix*n+iy]: inc exits
  result.area /= float(n*n)
  result.exitFraction = float(exits)/float(boundaryCount)

  var ss = 0.0
  var nk = 0
  for ix in 0..<n:
    for iy in 0..<n:
      let j = ix*n+iy
      if ix+1 < n:
        let d = k[(ix+1)*n+iy]-k[j]
        ss += d*d
        inc nk
      if iy+1 < n:
        let d = k[ix*n+iy+1]-k[j]
        ss += d*d
        inc nk
  result.stiffnessRoughness = sqrt(ss/float(nk))/h

proc writeControls(path: string; c: seq[float]) =
  let dir = path.parentDir
  if dir.len > 0: createDir(dir)
  var f = open(path, fmWrite)
  defer: f.close()
  for x in c: f.writeLine(&"{x:.17g}")

proc main() =
  let o = readOpt()
  installStandardParams()
  echoParams()
  processHelpParam()
  let fit =
    if o.geometry == "plaq4": scalarFit(o)
    elif o.geometry == "block5" and o.construction == "coupling": blockFit(o)
    else: linkFit(o)
  writeControls(o.output, fit.controls)
  echo &"wrote {fit.controls.len} controls to {o.output}"
  echo &"control fit rms={fit.residualRms:.8g} max={fit.residualMax:.8g}"
  echo &"CDF estimated relative error={fit.cdfError:.8g} " &
    &"roundoff bound={fit.cdfRound:.8g} " &
    &"transport residual={fit.cdfResidual:.8g}"
  echo &"sampled stage slope=[{fit.minSlope:.8g}, {fit.maxSlope:.8g}] minControl={fit.minControl:.8g}"
  let spec = buildMapSpec(mapParams(o, fit.controls), o.beta)
  if o.geometry == "plaq4": scalarMetrics(o, spec)
  elif o.geometry == "link2": linkMetrics(o, spec)
  let scan = mapScan(spec, o.scan, o.scanOutput)
  echoMapSummary(spec, scan, o.scan)
  let land = landscapeMetrics(o, spec)
  echo &"landscape barrier={land.barrier:.8g} boundaryFree={land.boundaryFree:.8g} " &
    &"boundaryEffective={land.boundaryEffective:.8g}"
  echo &"landscape dSeff=[{land.potentialMin:.8g}, {land.potentialMax:.8g}] " &
    &"minAt=({land.minX:.8g}, {land.minY:.8g}) maxAt=({land.maxX:.8g}, {land.maxY:.8g})"
  echo &"landscape barrierFromMin={land.barrierFromMin:.8g} " &
    &"boundaryFreeFromMin={land.boundaryFreeFromMin:.8g}"
  echo &"landscape area(+0.5)={land.area:.8g} exitFraction(+0.5)={land.exitFraction:.8g} " &
    &"stiffnessRoughness={land.stiffnessRoughness:.8g}"
  echo &"refined forceMax={land.forceMax:.8g} at ({land.forceX:.8g}, {land.forceY:.8g}) " &
    &"stiffnessMax={land.stiffnessMax:.8g} at ({land.stiffnessX:.8g}, {land.stiffnessY:.8g})"
  echo "HMC map arguments:"
  var args = &"-geometry:{o.geometry} -construction:{o.construction} -basis:bspline " &
    &"-mapDepth:{o.depth} -mapKnots:{o.knots} -ctxBasis:{o.ctxBasis} "
  if o.ctxBasis == "bspline": args.add &"-ctxKnots:{o.ctxKnots} "
  else: args.add &"-ctxOrder:{o.ctxOrder} "
  args.add &"-mapStrengths:1 -mapFloor:{o.floor:.17g} -mapParamFile:{o.output}"
  echo args
  processSaveParams()
  writeParamFile()

when isMainModule and not declared(optTest):
  qexInit()
  main()
  qexFinalize()
