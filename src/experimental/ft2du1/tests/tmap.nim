import base/globals
setDefaultNc(1)
setVLENmax(4)

import qex, math, algorithms/[cspline, numdiff]
import std/[os, strutils]
import gauge/gaugefix
import ../../graph/[core, scalar, gauge]
import ../../graph/hmcgauge/ftstout
import ../[flow, scan]

proc baseParams(geometry, construction, basis: string): MapParams =
  MapParams(geometry: geometry, construction: construction, basis: basis,
    ctxBasis: "fourier",
    mapDepth: 1, flowDepth: 1, mapStrengths: @[0.45],
    mapFloor: 0.15, mapEpsilon: 0.08, mapRho: 0.17, mapOrder: 3,
    ctxOrder: 1, mapKnots: 8, ctxKnots: 16,
    mapCoeffs: @[], mapControls: @[],
    mapParamFile: "", fejerOrders: @[8], fejerCenters: @[0.23],
    fejerWeights: @[1.0], mapStageOrder: "2031",
    mapDirs: @[0], mapParities: @[0], mapOffsets: @[0, 0],
    mapStride: (if geometry == "block5": 4 else: 0),
    mapInvTol: 2e-14, mapScanStep: 2e-5, mapInvIter: 100, mapScan: 24,
    mapDump: "", mapFunctionScan: 64, mapFunctionDump: "",
    mapFunctionContexts: @[], monitorEvery: 0, checkMap: false,
    startWidth: -1.0)

proc scalarTests() =
  type D3 = VectorArray[3, float]
  let maps = [
    ("identity", identityCircle()),
    ("sine", sineCircle([0.24, -0.05, 0.02], 0.8, 0.1)),
    ("sine-certified", sineCircle([0.6, 0.2], 1.0, 0.1)),
    ("sqfourier", sqFourierCircle([-0.60, 0.20, -0.10], 0.08, 0.8, 0.1)),
    ("cspline", cSplineCircle([-0.30, -0.18, 0.02, 0.21, 0.35, 0.16,
      -0.04, -0.25], 0.8, 0.1)),
    ("bspline", bSplineCircle([-0.30, -0.18, 0.02, 0.21, 0.35, 0.16,
      -0.04, -0.25], 0.1, 0.8, 0.1)),
    ("fejer", fejerCircle([7, 10], [0.23, -0.71], [1.0, 0.4], 0.72, 0.1)),
    ("fejer-grid", fejerCircle([12], [],
      [0.4, 0.7, 1.1, 1.6, 1.2, 0.8, 0.5, 0.6], 0.72, 0.1)),
    ("composed", composeCircle([
      sineCircle([0.15, -0.03], 0.7, 0.1),
      bSplineCircle([-0.20, -0.08, 0.14, 0.31, 0.22, -0.05],
        0.1, 0.6, 0.1)]))]
  let xs = [-7.1, -2.83, -1.71, -0.93, 0.17, 1.22, 2.74, 8.3]
  const h = 1e-3
  for (name, map0) in maps:
    var m = map0
    var dyErr, ddyErr, dddyErr, invErr, seamErr = 0.0
    for x in xs:
      proc vals(t: float): D3 =
        let e = circleEval(m, t)
        result[0] = e.y
        result[1] = e.dy
        result[2] = e.ddy
      var d, de: D3
      ndiff(d, de, vals, x, h, ordMax = 5)
      let
        e = circleEval(m, x)
        back = circleInv(m, e.y)
        seam = circleEval(m, x+TAU)
      dyErr = max(dyErr, abs(d[0]-e.dy))
      ddyErr = max(ddyErr, abs(d[1]-e.ddy))
      dddyErr = max(dddyErr, abs(d[2]-e.dddy))
      invErr = max(invErr, abs(back-x))
      seamErr = max(seamErr, abs(seam.y-e.y-TAU))
      seamErr = max(seamErr, abs(seam.dy-e.dy))
      seamErr = max(seamErr, abs(seam.ddy-e.ddy))
      seamErr = max(seamErr, abs(seam.dddy-e.dddy))
      doAssert e.dy >= 0.1-2e-13
    doAssert dyErr < 3e-9
    doAssert ddyErr < 3e-8
    doAssert dddyErr < 5e-7
    doAssert invErr < 2e-11
    doAssert seamErr < 2e-11
    echo "scalar[", name, "]: dy=", dyErr, " ddy=", ddyErr,
      " dddy=", dddyErr, " inv=", invErr, " seam=", seamErr

proc pairTests() =
  type D6 = VectorArray[6, float]
  for name in ["identity", "sine", "sqfourier", "bspline",
               "tensor-bspline", "fejer", "stout"]:
    let basis = if name == "tensor-bspline": "bspline" else: name
    var p = baseParams("link2", "scalar", basis)
    if name == "tensor-bspline": p.ctxBasis = "bspline"
    if basis == "sqfourier":
      p.mapDepth = 2
      p.mapStrengths = @[0.25, 0.40]
    let
      spec = buildMapSpec(p, 3.5)
      m = spec.pair
    var jacErr, detErr, logjErr, forceErr, invErr, seamErr, mixedErr = 0.0
    for (pp0, pm0) in [(0.77, -0.29), (0.82, -0.46),
                       (0.20, 0.10), (2.7, -2.1)]:
      var
        pp = pp0
        pm = pm0
      proc vals(xp, xm: float): D6 =
        let e = evalPair(m, xp, xm)
        result[0] = e.physicalPlus
        result[1] = e.physicalMinus
        result[2] = e.delta
        result[3] = e.det
        result[4] = e.logdet
        result[5] = e.action
      proc fx(x: float): D6 = vals(x, pm)
      proc fm0(x: float): D6 = vals(pp, x)
      var dx, dy, ex, ey: D6
      ndiff(dx, ex, fx, pp, 1e-3, ordMax = 5)
      ndiff(dy, ey, fm0, pm, 1e-3, ordMax = 5)
      let
        e = evalPair(m, pp, pm)
        back = invertPair(m, e.physicalPlus, e.physicalMinus)
        det = e.jac[0][0]*e.jac[1][1]-e.jac[0][1]*e.jac[1][0]
      for i in 0..1:
        jacErr = max(jacErr, abs(dx[i]-e.jac[i][0]))
        jacErr = max(jacErr, abs(dy[i]-e.jac[i][1]))
      jacErr = max(jacErr, abs(dx[2]-e.deltaPlus))
      jacErr = max(jacErr, abs(dy[2]-e.deltaMinus))
      detErr = max(detErr, abs(dx[3]-e.detPlus))
      detErr = max(detErr, abs(dy[3]-e.detMinus))
      logjErr = max(logjErr, abs(dx[4]-e.logdetPlus))
      logjErr = max(logjErr, abs(dy[4]-e.logdetMinus))
      forceErr = max(forceErr, abs(dx[5]-e.forcePlus))
      forceErr = max(forceErr, abs(dy[5]-e.forceMinus))
      invErr = max(invErr, abs(back[0]-pp))
      invErr = max(invErr, abs(back[1]-pm))
      doAssert det > 0.0 and abs(det-e.det) < 2e-13
      doAssert abs(ln(det)-e.logdet) < 2e-13
      if name == "stout":
        let
          mean = 0.5*(pp+pm)
          u = 0.5*(pp-pm)
          q = u-2.0*p.mapRho*cos(mean)*sin(u)
          exactDet = 1.0-p.mapRho*(cos(pp)+cos(pm))
        doAssert abs(e.physicalPlus-(mean+q)) < 2e-14
        doAssert abs(e.physicalMinus-(mean-q)) < 2e-14
        doAssert abs(e.det-exactDet) < 2e-14
      for k in 0..1:
        let s = if k == 0: evalPair(m, pp+TAU, pm)
          else: evalPair(m, pp, pm+TAU)
        seamErr = max(seamErr, abs(s.physicalPlus-e.physicalPlus-
          (if k == 0: TAU else: 0.0)))
        seamErr = max(seamErr, abs(s.physicalMinus-e.physicalMinus-
          (if k == 1: TAU else: 0.0)))
        seamErr = max(seamErr, abs(s.logdet-e.logdet))
    for (u, mean) in [(-2.3, -1.1), (0.41, 0.7), (1.7, -2.0)]:
      let
        e = evalPairCore(m, u, mean)
        t = evalPairCore(m, u+PI, mean+PI)
        hh = 2e-5
        num = (evalPairCore(m, u+hh, mean+hh).du-
          evalPairCore(m, u-hh, mean-hh).du)/(2.0*hh)
      seamErr = max(seamErr, abs(t.y-e.y-PI))
      seamErr = max(seamErr, abs(t.du-e.du))
      mixedErr = max(mixedErr, abs(num-(e.duu+e.dum)))
    doAssert jacErr < 2e-8
    doAssert detErr < 3e-8
    doAssert logjErr < 3e-8
    doAssert forceErr < 5e-8
    doAssert invErr < 3e-11
    doAssert seamErr < 3e-10
    doAssert mixedErr < 3e-7
    echo "pair[", name, "]: jac=", jacErr, " det=", detErr, " logj=", logjErr,
      " force=", forceErr, " inv=", invErr, " seam=", seamErr

  let
    c = 0.4
    eps = 0.08
    s = 0.3
    f = 0.2
    x = 0.37
    b0 = 1.0+0.5*c*c+eps
    mode = (2.0*c*cos(x)+0.5*c*c*cos(2.0*x))/b0
    a = circleEval(sqFourierCircle([c], eps, s, f), x)
    b = evalContext(sqFourierContextMap(1, 1, 0, [c], eps, s, f),
      x, [0.0, 0.0, 0.0, 0.0])
  doAssert abs(a.dy-(1.0+s*mode)) < 2e-14
  doAssert abs(b.dx-(1.0+s*(1.0-f)*mode)) < 2e-14

proc contextSplineTests() =
  const
    nk = 8
    nc = 128
    nr = 2
  let half = nk div 2
  var f = newSeq[float](half*(1+nr))
  for i in 0..<half:
    let x = -PI+TAU*float(i)/float(nk)
    f[i*3] = 1.0+0.08*cos(2.0*x)
    f[i*3+1] = 0.16*cos(x)
    f[i*3+2] = 0.04*cos(2.0*x)
  let fm = bSplineContextMap(1, nr, nk, f, 0.85, 0.1,
    twisted = true)
  var t = newSeq[float](nk*nc)
  for i in 0..<half:
    var y = newSeq[float](nc)
    for j in 0..<nc:
      let c = -PI+TAU*float(j)/float(nc)
      y[j] = f[i*3]+f[i*3+1]*cos(c)+f[i*3+2]*cos(2.0*c)
    let ctl = periodicBSplineControls(y)
    for j in 0..<nc:
      t[i*nc+j] = ctl[j]
      t[(i+half)*nc+j] = ctl[(j+nc div 2) mod nc]
  let tm = tensorBSplineContextMap(1, nk, nc, t, 0.85, 0.1,
    twisted = true)
  var yerr, dxerr, dxxerr, dcerr, dxcerr = 0.0
  for ix in 0..<37:
    let x = -PI+TAU*(float(ix)+0.37)/37.0
    for ic in 0..<41:
      let
        c = -PI+TAU*(float(ic)+0.23)/41.0
        ctx = [c, 0.0, 0.0, 0.0]
        a = evalContext(fm, x, ctx)
        b = evalContext(tm, x, ctx)
      yerr = max(yerr, abs(a.y-b.y))
      dxerr = max(dxerr, abs(a.dx-b.dx))
      dxxerr = max(dxxerr, abs(a.dxx-b.dxx))
      dcerr = max(dcerr, abs(a.dc[0]-b.dc[0]))
      dxcerr = max(dxcerr, abs(a.dxc[0]-b.dxc[0]))
  doAssert yerr < 2e-7
  doAssert dxerr < 3e-7
  doAssert dxxerr < 2e-6
  doAssert dcerr < 4e-6
  doAssert dxcerr < 5e-6
  echo "context[Fourier/tensor-bspline]: y=", yerr, " dx=", dxerr,
    " dxx=", dxxerr, " dc=", dcerr, " dxc=", dxcerr

proc stagedParameterTests() =
  var p = baseParams("link2", "scalar", "sine")
  p.mapDepth = 2
  p.mapStrengths = @[1.0]
  p.mapOrder = 2
  p.ctxOrder = 1
  p.mapCoeffs = @[0.0, 0.025, 0.012, 0.0,
                  0.0, -0.031, -0.009, 0.0]
  let
    m = buildMapSpec(p, 3.0).pair
    want = PairMap(beta: 3.0, invTol: p.mapInvTol, invIter: p.mapInvIter,
      stages: @[
        sineContextMap(1, 2, 1, p.mapCoeffs[0..3], 1.0, p.mapFloor,
          twisted = true),
        sineContextMap(1, 2, 1, p.mapCoeffs[4..7], 1.0, p.mapFloor,
          twisted = true)])
  doAssert m.stages.len == 2
  for (pp, pm) in [(0.7, -0.4), (2.1, 0.3), (-1.6, 2.4)]:
    let
      a = evalPair(m, pp, pm)
      b = evalPair(want, pp, pm)
    doAssert abs(a.physicalPlus-b.physicalPlus) < 2e-14
    doAssert abs(a.physicalMinus-b.physicalMinus) < 2e-14
    doAssert abs(a.logdet-b.logdet) < 2e-14
    doAssert abs(a.forcePlus-b.forcePlus) < 2e-13
    doAssert abs(a.forceMinus-b.forceMinus) < 2e-13

  var q = baseParams("plaq4", "scalar", "sine")
  q.mapDepth = 2
  q.mapStrengths = @[1.0]
  q.mapOrder = 2
  q.mapCoeffs = @[0.03, -0.01, -0.02, 0.015]
  let
    c = buildMapSpec(q, 3.0).circle
    cw = composeCircle([sineCircle(q.mapCoeffs[0..1], 1.0, q.mapFloor),
                        sineCircle(q.mapCoeffs[2..3], 1.0, q.mapFloor)])
  for x in [-2.3, -0.4, 0.8, 2.7]:
    let
      a = circleEval(c, x)
      b = circleEval(cw, x)
    doAssert abs(a.y-b.y) < 2e-14
    doAssert abs(a.dy-b.dy) < 2e-14
    doAssert abs(a.ddy-b.ddy) < 2e-14

  var r = baseParams("block5", "coupling", "bspline")
  r.mapDepth = 2
  r.mapStrengths = @[1.0]
  r.mapKnots = 4
  r.ctxOrder = 1
  let unit = r.mapKnots*(1+2*r.ctxOrder)
  r.mapControls = newSeq[float](4*r.mapDepth*unit)
  for stage in 0..<r.mapDepth:
    for coord in 0..3:
      let off = (4*stage+coord)*unit
      for knot in 0..<r.mapKnots:
        r.mapControls[off+knot*3] = 1.0+0.05*float(4*stage+coord)
        r.mapControls[off+knot*3+1] = 0.01
        r.mapControls[off+knot*3+2] = -0.005
  let bm = buildMapSpec(r, 3.0).blockMap
  for stage in 0..<r.mapDepth:
    for coord in 0..3:
      doAssert bm.stages[coord].len == r.mapDepth
      doAssert abs(bm.stages[coord][stage].controls[0]-
        (1.0+0.05*float(4*stage+coord))) < 2e-14

proc blockTests() =
  type D7 = VectorArray[7, float]
  type D8 = VectorArray[8, float]
  let
    z: Vec4 = [0.31, -0.52, 0.77, -0.18]
    flux = 0.43
  for name in ["identity", "sine", "sqfourier", "bspline", "tensor-bspline"]:
    let basis = if name == "tensor-bspline": "bspline" else: name
    var p = baseParams("block5", "coupling", basis)
    if name == "tensor-bspline":
      p.ctxBasis = "bspline"
      p.mapDepth = 2
      p.mapStrengths = @[0.20, 0.35]
    elif basis == "bspline":
      p.mapDepth = 2
      p.mapStrengths = @[0.20, 0.35]
    let
      spec = buildMapSpec(p, 3.0)
      m = spec.blockMap
      e = evalBlockMap(m, z, flux)
    var jacErr, detErr, logjErr, forceErr, fluxErr, invErr = 0.0
    for jj in 0..3:
      var j = jj
      proc values(x: float): D7 =
        var q = z
        q[j] = x
        let a = evalBlockMap(m, q, flux)
        for i in 0..3: result[i] = a.physicalZ[i]
        result[4] = a.det
        result[5] = a.logdet
        result[6] = a.action
      var d, de: D7
      ndiff(d, de, values, z[j], 1e-3, ordMax = 5)
      for i in 0..3: jacErr = max(jacErr, abs(d[i]-e.jacobian[i][j]))
      detErr = max(detErr, abs(d[4]-e.detGradient[j]))
      logjErr = max(logjErr, abs(d[5]-e.logdetGradient[j]))
      forceErr = max(forceErr, abs(d[6]-e.gradient[j]))
    proc fluxValues(x: float): D8 =
      let a = evalBlockMap(m, z, x)
      for i in 0..4: result[i] = a.physical[i]
      result[5] = a.det
      result[6] = a.logdet
      result[7] = a.action
    var df, dfe: D8
    ndiff(df, dfe, fluxValues, flux, 1e-3, ordMax = 5)
    for i in 0..4: fluxErr = max(fluxErr, abs(df[i]-e.derivatives[i][4]))
    detErr = max(detErr, abs(df[5]-e.detFluxGradient))
    fluxErr = max(fluxErr, abs(df[6]-e.logdetFluxGradient))
    fluxErr = max(fluxErr, abs(df[7]-e.fluxForce))
    let
      d = det4(e.jacobian)
      back = invertBlockMap(m, e.physicalZ, flux)
    doAssert d > 0.0 and abs(d-e.det) < 3e-12
    doAssert abs(ln(d)-e.logdet) < 3e-12
    for i in 0..3: invErr = max(invErr, abs(back[i]-z[i]))
    doAssert jacErr < 3e-8
    doAssert detErr < 8e-8
    doAssert logjErr < 5e-8
    doAssert forceErr < 8e-8
    doAssert fluxErr < 8e-8
    doAssert invErr < 5e-11
    echo "block[", name, "]: jac=", jacErr, " det=", detErr, " logj=", logjErr,
      " force=", forceErr, " flux=", fluxErr, " inv=", invErr

proc moveLocal(p: Vec5; np, j: int; h: float): Vec5 =
  result = p
  if np == 2:
    result[0] += h
    result[1] -= h
  else:
    result[0] += h
    result[j+1] -= h

proc numericalLocalDet(spec: MapSpec; p: Vec5; h: float): float =
  let e = evalLocal(spec, p)
  if e.nv == 1:
    let
      ep = evalLocal(spec, moveLocal(p, e.np, 0, h))
      em = evalLocal(spec, moveLocal(p, e.np, 0, -h))
    return (ep.physical[0]-em.physical[0])/(2.0*h)
  var a: flow.Mat4
  for j in 0..<e.nv:
    let
      ep = evalLocal(spec, moveLocal(p, e.np, j, h))
      em = evalLocal(spec, moveLocal(p, e.np, j, -h))
    for i in 0..<e.nv:
      a[i][j] = -(ep.physical[i+1]-em.physical[i+1])/(2.0*h)
  det4(a)

proc localTests() =
  let cases = [
    ("plaq4", "scalar", "identity"),
    ("plaq4", "scalar", "sine"),
    ("plaq4", "scalar", "sqfourier"),
    ("plaq4", "scalar", "cspline"),
    ("plaq4", "scalar", "bspline"),
    ("plaq4", "scalar", "fejer"),
    ("link2", "scalar", "identity"),
    ("link2", "scalar", "sine"),
    ("link2", "scalar", "sqfourier"),
    ("link2", "scalar", "bspline"),
    ("link2", "scalar", "tensor-bspline"),
    ("link2", "scalar", "fejer"),
    ("link2", "scalar", "stout"),
    ("block5", "chain", "identity"),
    ("block5", "chain", "sine"),
    ("block5", "chain", "sqfourier"),
    ("block5", "chain", "bspline"),
    ("block5", "chain", "tensor-bspline"),
    ("block5", "chain", "fejer"),
    ("block5", "coupling", "identity"),
    ("block5", "coupling", "sine"),
    ("block5", "coupling", "sqfourier"),
    ("block5", "coupling", "bspline"),
    ("block5", "coupling", "tensor-bspline")]
  const h = 2e-5
  for (geometry, construction, name) in cases:
    let basis = if name == "tensor-bspline": "bspline" else: name
    var p0 = baseParams(geometry, construction, basis)
    if name == "tensor-bspline": p0.ctxBasis = "bspline"
    let
      spec = buildMapSpec(p0, 3.0)
      p: Vec5 = if geometry == "link2":
        [0.41, -0.73, 0.0, 0.0, 0.0]
        else: [0.31, -0.52, 0.77, -0.18, 0.43]
      e = evalLocal(spec, p)
      back = invertLocal(spec, e.physical)
      det = numericalLocalDet(spec, p, h)
      point = scanPoint(spec, e.physical)
    var invErr, forceErr = 0.0
    for i in 0..<e.np:
      invErr = max(invErr, abs(principalAngle(back[i]-p[i])))
    for j in 0..<e.nv:
      let
        ep = evalLocal(spec, moveLocal(p, e.np, j, h))
        em = evalLocal(spec, moveLocal(p, e.np, j, -h))
        g = (ep.seff-em.seff)/(2.0*h)
      forceErr = max(forceErr, abs(g-e.force[j]))
    doAssert invErr < 8e-11
    doAssert det > 0.0 and abs(ln(det)-e.logj) < 2e-7
    doAssert forceErr < 2e-7
    doAssert point.invErr < 8e-11
    doAssert point.hessAsym < 2e-5
    if geometry == "block5" and construction == "chain":
      doAssert spec.chainOrder == [2, 0, 3, 1]
    echo "local[", geometry, "/", construction, "/", name,
      "]: det=", abs(ln(det)-e.logj), " force=", forceErr,
      " inv=", invErr, " hessAsym=", point.hessAsym

proc identityScanTests() =
  const beta = 2.5
  for (geometry, construction) in [
      ("link2", "scalar"), ("plaq4", "scalar"),
      ("block5", "chain"), ("block5", "coupling")]:
    let
      spec = buildMapSpec(baseParams(geometry, construction, "identity"), beta)
      point = scanPoint(spec, [0.0, 0.0, 0.0, 0.0, 0.0])
    doAssert abs(point.logj) < 1e-14
    doAssert abs(point.dseff) < 1e-14
    doAssert point.force2 < 1e-24
    doAssert point.invErr < 1e-14
    if geometry == "link2":
      doAssert abs(point.hessian[0][0]-2.0*beta) < 2e-8
      doAssert abs(point.stiffMax-2.0*beta) < 2e-8
      doAssert abs(point.stiffRms-2.0*beta) < 2e-8
      doAssert abs(point.curvMin-2.0*beta) < 2e-8
      doAssert abs(point.curvMax-2.0*beta) < 2e-8
    else:
      for i in 0..3:
        for j in 0..3:
          let exact = beta*(if i == j: 2.0 else: 1.0)
          doAssert abs(point.hessian[i][j]-exact) < 2e-8
      doAssert abs(point.stiffMax-5.0*beta) < 3e-8
      doAssert abs(point.stiffRms-beta*sqrt(7.0)) < 3e-8
      doAssert abs(point.curvMin-beta) < 3e-8
      doAssert abs(point.curvMax-5.0*beta) < 3e-8

proc hessianStepTests() =
  for (geometry, construction) in [
      ("link2", "scalar"), ("plaq4", "scalar"),
      ("block5", "chain"), ("block5", "coupling")]:
    var coarse = buildMapSpec(baseParams(geometry, construction, "sine"), 3.0)
    let
      p: Vec5 = if geometry == "link2":
        [0.36, -0.61, 0.0, 0.0, 0.0]
        else: [0.27, -0.44, 0.63, -0.16, 0.38]
      physical = evalLocal(coarse, p).physical
      a = scanPoint(coarse, physical)
    var fine = coarse
    fine.scanStep = 0.5*coarse.scanStep
    let b = scanPoint(fine, physical)
    var err, scale = 0.0
    let nv = if geometry == "link2": 1 else: 4
    for i in 0..<nv:
      for j in 0..<nv:
        err = max(err, abs(a.hessian[i][j]-b.hessian[i][j]))
        scale = max(scale, abs(b.hessian[i][j]))
    doAssert err < 2e-5*(1.0+scale)

proc scanOutputTests() =
  let
    spec = buildMapSpec(baseParams("link2", "scalar", "identity"), 3.0)
    path = "/tmp/qex_ft2du1_scan_" & $myRank & ".dat"
  if fileExists(path): removeFile(path)
  let s = mapScan(spec, 8, path)
  let lines = readFile(path).splitLines
  var rows, gaps = 0
  var hasMass, hasColumns = false
  for line in lines:
    if line.len == 0:
      inc gaps
    elif line[0] == '#':
      hasMass = hasMass or line.contains("mass=1 localRounds=1")
      hasColumns = hasColumns or
        line.contains("force2 stiffMax stiffRms curvMin") and
        line.contains("curvMax")
    else:
      inc rows
  removeFile(path)
  doAssert rows == 64 and gaps >= 8
  doAssert hasMass and hasColumns
  doAssert abs(s.minLogj) < 1e-14 and abs(s.maxLogj) < 1e-14
  doAssert s.maxInverseError < 1e-13 and s.maxHessAsym < 1e-8

  let
    stoutSpec = buildMapSpec(baseParams("link2", "scalar", "stout"), 3.0)
    stoutPath = "/tmp/qex_ft2du1_stout_scan_" & $myRank & ".dat"
  if fileExists(stoutPath): removeFile(stoutPath)
  let stoutScan = mapScan(stoutSpec, 8, stoutPath)
  let stoutHeader = readFile(stoutPath)
  removeFile(stoutPath)
  doAssert stoutHeader.contains("basis=stout") and
    stoutHeader.contains("mapRho=0.17")
  doAssert abs(stoutScan.minLogj-ln(1.0-2.0*stoutSpec.mapRho)) < 2e-14
  doAssert abs(stoutScan.maxLogj-ln(1.0+2.0*stoutSpec.mapRho)) < 2e-14
  doAssert stoutScan.maxInverseError < 3e-13

proc functionOutputTests() =
  let cases = [
    ("plaq4", "scalar", "cspline", newSeq[float](), 1),
    ("link2", "scalar", "stout", @[-0.7, 0.0, 0.7], 3),
    ("block5", "coupling", "bspline",
      @[-0.4, 0.2, 0.7, 0.1, 0.3, -0.5, 0.8, -0.2], 8)]
  const n = 16
  for (geometry, construction, basis, contexts, sets) in cases:
    let
      spec = buildMapSpec(baseParams(geometry, construction, basis), 3.0)
      path = "/tmp/qex_ft2du1_function_" & geometry & "_" & $myRank & ".dat"
    if fileExists(path): removeFile(path)
    dumpMapFunction(spec, n, contexts, path)
    var rows, datasets = 0
    for line in readFile(path).splitLines:
      if line.startsWith("# dataset="):
        inc datasets
      elif line.len > 0 and line[0] != '#':
        let v = line.splitWhitespace
        doAssert v.len == 22
        let
          gp = parseFloat(v[6])
          logGp = parseFloat(v[8])
        doAssert gp > 0.0
        doAssert abs(logGp-ln(gp)) < 2e-13
        inc rows
    removeFile(path)
    doAssert datasets == sets
    doAssert rows == sets*(n+1)

proc validationTests() =
  var rejected = 0
  for p in [baseParams("plaq4", "chain", "sine"),
            baseParams("link2", "coupling", "bspline"),
            baseParams("block5", "scalar", "sqfourier"),
            baseParams("block5", "coupling", "fejer")]:
    try:
      discard buildMapSpec(p, 3.0)
    except ValueError:
      inc rejected
  doAssert rejected == 4
  var bad = baseParams("plaq4", "scalar", "sine")
  bad.mapCoeffs = @[3.0]
  try:
    discard buildMapSpec(bad, 3.0)
  except ValueError:
    inc rejected
  doAssert rejected == 5
  bad = baseParams("link2", "scalar", "bspline")
  bad.mapKnots = 7
  try:
    discard buildMapSpec(bad, 3.0)
  except ValueError:
    inc rejected
  doAssert rejected == 6

  bad = baseParams("link2", "scalar", "sine")
  bad.mapScanStep = 0.0
  try:
    discard buildMapSpec(bad, 3.0)
  except ValueError:
    inc rejected
  doAssert rejected == 7

  let stoutGeometry = baseParams("plaq4", "scalar", "stout")
  var stoutDepth = baseParams("link2", "scalar", "stout")
  stoutDepth.mapDepth = 2
  var stoutRho = baseParams("link2", "scalar", "stout")
  stoutRho.mapRho = 0.5
  var stoutFloor = baseParams("link2", "scalar", "stout")
  stoutFloor.mapRho = 0.43
  var stoutCoeffs = baseParams("link2", "scalar", "stout")
  stoutCoeffs.mapCoeffs = @[0.1]
  for p in [stoutGeometry, stoutDepth, stoutRho, stoutFloor, stoutCoeffs]:
    try:
      discard buildMapSpec(p, 3.0)
    except ValueError:
      inc rejected
  doAssert rejected == 12

  var badContext = baseParams("link2", "scalar", "bspline")
  badContext.ctxBasis = "polynomial"
  var oddContext = baseParams("link2", "scalar", "bspline")
  oddContext.ctxBasis = "bspline"
  oddContext.ctxKnots = 15
  for p in [badContext, oddContext]:
    try:
      discard buildMapSpec(p, 3.0)
    except ValueError:
      inc rejected
  doAssert rejected == 14

  var noContext = baseParams("link2", "scalar", "sqfourier")
  noContext.ctxOrder = 0
  let e = evalPair(buildMapSpec(noContext, 3.0).pair, 0.4, -0.2)
  doAssert e.jac[0][0]*e.jac[1][1]-e.jac[0][1]*e.jac[1][0] > 0.0

  var badFourierOrder = baseParams("link2", "scalar", "bspline")
  badFourierOrder.ctxOrder = -1
  try:
    discard buildMapSpec(badFourierOrder, 3.0)
  except ValueError:
    inc rejected
  doAssert rejected == 15

  var tensorUnusedOrder = baseParams("link2", "scalar", "bspline")
  tensorUnusedOrder.ctxBasis = "bspline"
  tensorUnusedOrder.ctxOrder = -1
  discard buildMapSpec(tensorUnusedOrder, 3.0)

  var fourierUnusedKnots = baseParams("link2", "scalar", "bspline")
  fourierUnusedKnots.ctxKnots = 1
  discard buildMapSpec(fourierUnusedKnots, 3.0)

  var scalarUnusedContext = baseParams("plaq4", "scalar", "bspline")
  scalarUnusedContext.ctxBasis = "bspline"
  scalarUnusedContext.ctxOrder = -1
  scalarUnusedContext.ctxKnots = 1
  discard buildMapSpec(scalarUnusedContext, 3.0)

  for basis in ["sine", "sqfourier", "cspline", "bspline", "fejer"]:
    for floor in [0.5*mapDerivativeFloor, 1.0]:
      var p = baseParams("plaq4", "scalar", basis)
      p.mapFloor = floor
      try:
        discard buildMapSpec(p, 3.0)
      except ValueError:
        inc rejected
  doAssert rejected == 25

  for basis in ["sine", "sqfourier", "cspline", "bspline", "fejer"]:
    var p = baseParams("plaq4", "scalar", basis)
    p.mapFloor = mapDerivativeFloor
    discard buildMapSpec(p, 3.0)

  for floor in [0.0, 1.0]:
    var p = baseParams("plaq4", "scalar", "identity")
    p.mapFloor = floor
    discard buildMapSpec(p, 3.0)

  var stoutFloorRange = baseParams("link2", "scalar", "stout")
  stoutFloorRange.mapFloor = 1.0
  try:
    discard buildMapSpec(stoutFloorRange, 3.0)
  except ValueError:
    inc rejected
  doAssert rejected == 26

proc runPureMapTests*() =
  scalarTests()
  pairTests()
  contextSplineTests()
  stagedParameterTests()
  blockTests()
  localTests()
  identityScanTests()
  hessianStepTests()
  scanOutputTests()
  functionOutputTests()
  validationTests()

proc maxLinkDiff2(a, b: seq): float =
  var m = 0.0
  threads:
    var mm = 0.0
    for mu in 0..<a.len:
      for x in a[mu]: mm = max(mm, (a[mu][x]-b[mu][x]).norm2.simdMax)
    mm.threadRankMax
    threadSingle: m = mm
  m

proc copyGauge(g: gauge.Gauge): gauge.Gauge =
  let r = g[0].l.newgauge
  threads:
    for mu in 0..<g.len: r[mu] := g[mu]
  r

proc rotateLink(g: gauge.Gauge; mu, x, y: int; a: float) =
  let
    lo = g[0].l
    xx = ((x mod lo[0])+lo[0]) mod lo[0]
    yy = ((y mod lo[1])+lo[1]) mod lo[1]
    ri = lo.rankIndex([xx, yy])
  if ri.rank == lo.myRank:
    let
      z = g[mu]{ri.index}[0, 0]
      zr = z.re[][]
      zi = z.im[][]
      c = cos(a)
      s = sin(a)
    g[mu]{ri.index}[0, 0].re = c*zr-s*zi
    g[mu]{ri.index}[0, 0].im = s*zr+c*zi

proc globalJacobianTests() =
  const
    l = 4
    n = 2*l*l
    h = 1e-2
    beta = 3.0
  type Col = VectorArray[n, float]
  let lo = @[l, l].newLayout
  var
    rng = lo.newRNGField(MRG32k3a, 246813579'u)
    v = lo.newgauge
  v.warm(0.18, rng)
  let cases = [
    ("plaq4", "scalar", "sine"),
    ("plaq4", "scalar", "sqfourier"),
    ("plaq4", "scalar", "cspline"),
    ("plaq4", "scalar", "bspline"),
    ("plaq4", "scalar", "fejer"),
    ("link2", "scalar", "sqfourier"),
    ("link2", "scalar", "bspline"),
    ("link2", "scalar", "tensor-bspline"),
    ("link2", "scalar", "fejer"),
    ("link2", "scalar", "stout"),
    ("block5", "chain", "sqfourier"),
    ("block5", "chain", "tensor-bspline"),
    ("block5", "coupling", "sqfourier"),
    ("block5", "coupling", "bspline"),
    ("block5", "coupling", "tensor-bspline")]
  for (geometry, construction, name) in cases:
    let basis = if name == "tensor-bspline": "bspline" else: name
    var p = baseParams(geometry, construction, basis)
    if name == "tensor-bspline": p.ctxBasis = "bspline"
    if geometry == "link2" and basis == "sqfourier":
      p.mapDepth = 2
      p.mapStrengths = @[0.20, 0.30]
      p.flowDepth = 2
    if basis == "stout":
      p.mapStride = 4
      p.mapOffsets = @[1, 2]
      p.mapDirs = @[0, 1]
    let
      spec = buildMapSpec(p, beta)
      layout = buildMapLayout(v[0], p, spec)
      sm = mapHost(v, spec, layout)
    var jac: MatrixArray[n, n, float]
    for r in 0..<n:
      for c in 0..<n: jac[r, c] = 0.0
    for c in 0..<n:
      let
        lex = c div 2
        mu = c mod 2
      proc column(a: float): Col =
        var vt = copyGauge(v)
        rotateLink(vt, mu, lex mod l, lex div l, a)
        let u = mapHost(vt, spec, layout).u
        for olex in 0..<l*l:
          let ri = lo.rankIndex(olex)
          if ri.rank == lo.myRank:
            for nu in 0..1:
              let
                z = u[nu]{ri.index}[0, 0]
                z0 = sm.u[nu]{ri.index}[0, 0]
                re = z.re[][]*z0.re[][]+z.im[][]*z0.im[][]
                im = z.im[][]*z0.re[][]-z.re[][]*z0.im[][]
              result[2*olex+nu] = atan2(im, re)
      var col, ce: Col
      ndiff(col, ce, column, 0.0, h, ordMax = 4)
      for r in 0..<n: jac[r, c] = col[r]
    rankSum(jac)
    let d = determinant(jac)
    doAssert d > 0.0
    let err = abs(ln(d)-sm.lndet)
    echo "global[", geometry, "/", construction, "/", name,
      "]: det=", d, " logdet error=", err
    doAssert err < 2e-7

proc latticeTests() =
  const beta = 3.0
  let
    lo = @[8, 8].newLayout
    grt = initGraphRuntime()
    gc = actWilson(scalar.toGvalue(grt, beta))
  var
    rng = lo.newRNGField(MRG32k3a, 1234567891'u)
    v0 = lo.newgauge
    u0 = lo.newgauge
    vgt = lo.newgauge
    gt = lo.ColorMatrix()
  v0.warm(0.22, rng)
  u0.warm(0.17, rng)
  threads: gt.randomU rng
  vgt.gaugeTransform(v0, gt)

  proc checkGrad(name: string; action: proc(V: Ggauge): Gscalar; tol: float) =
    var r = lo.newgauge
    r.randomTAH rng
    let
      vg = gauge.toGvalue(grt, v0)
      rg = gauge.toGvalue(grt, r)
      t = scalar.toGvalue(grt, 0.0)
      vt = axexpmuly(t, rg, vg)
      f = action(vt)
      analytic = f.grad(t).eval.sval
    proc value(x: float): float =
      t.update(x)
      f.eval.sval
    var numeric, err: float
    ndiff(numeric, err, value, 0.0, 1e-3, ordMax = 8)
    t.update(0.0)
    let rel = abs(numeric-analytic)/(abs(numeric)+abs(analytic)+1e-30)
    echo "graph[", name, "]: analytic=", analytic, " numeric=", numeric,
      " rel=", rel
    doAssert rel < tol

  let cases = [
    ("plaq4", "scalar", "identity"),
    ("plaq4", "scalar", "sine"),
    ("plaq4", "scalar", "sqfourier"),
    ("plaq4", "scalar", "cspline"),
    ("plaq4", "scalar", "bspline"),
    ("plaq4", "scalar", "fejer"),
    ("link2", "scalar", "sine"),
    ("link2", "scalar", "sqfourier"),
    ("link2", "scalar", "bspline"),
    ("link2", "scalar", "tensor-bspline"),
    ("link2", "scalar", "fejer"),
    ("link2", "scalar", "stout"),
    ("block5", "chain", "sqfourier"),
    ("block5", "chain", "bspline"),
    ("block5", "chain", "tensor-bspline"),
    ("block5", "coupling", "sine"),
    ("block5", "coupling", "sqfourier"),
    ("block5", "coupling", "bspline"),
    ("block5", "coupling", "tensor-bspline")]

  for (geometry, construction, profile) in cases:
    let basis = if profile == "tensor-bspline": "bspline" else: profile
    var p = baseParams(geometry, construction, basis)
    if profile == "tensor-bspline": p.ctxBasis = "bspline"
    if geometry == "block5" and construction == "chain" and basis == "sqfourier":
      p.mapDepth = 2
      p.mapStrengths = @[0.20, 0.30]
      p.flowDepth = 2
    if basis == "stout":
      p.mapStride = 4
      p.mapOffsets = @[1, 2]
      p.mapDirs = @[0, 1]
    let
      spec = buildMapSpec(p, beta)
      layout = buildMapLayout(v0[0], p, spec)
      action = mapAction(grt, gc, spec, layout)
      hs = mapHost(v0, spec, layout)
      ghs = mapHost(vgt, spec, layout)
      vg = gauge.toGvalue(grt, v0)
      vgg = gauge.toGvalue(grt, vgt)
      se = action(vg)
      seg = action(vgg)
      seff = se.eval.sval
      gseff = seg.eval.sval
      literal = gaugeAction(gc, gauge.toGvalue(grt, hs.u)).eval.sval-hs.lndet
      force = contractProjTAH(grad(se, vg), vg)
      forceg = contractProjTAH(grad(seg, vgg), vgg)
      name = geometry & "/" & construction & "/" & profile
    let actionErr = abs(seff-literal)/(1.0+abs(literal))
    echo "literal[", name, "]: err=", actionErr, " logdet=", hs.lndet
    doAssert actionErr < 2e-10

    var transformed = lo.newgauge
    transformed.gaugeTransform(hs.u, gt)
    let
      mapGaugeErr = maxLinkDiff2(ghs.u, transformed)
      actionGaugeErr = abs(gseff-seff)/(1.0+abs(seff))
      logdetGaugeErr = abs(ghs.lndet-hs.lndet)
      forceGaugeErr = (forceg-force).norm2.eval.sval/(1.0+force.norm2.eval.sval)
    doAssert mapGaugeErr < 2e-19
    doAssert actionGaugeErr < 3e-10
    doAssert logdetGaugeErr < 3e-10
    doAssert forceGaugeErr < 2e-18

    var back = copyGauge(hs.u)
    invertMapFlow(back, spec, layout)
    let invErr = maxLinkDiff2(back, v0)
    var aux = copyGauge(u0)
    invertMapFlow(aux, spec, layout)
    let redo = mapHost(aux, spec, layout).u
    let physInvErr = maxLinkDiff2(redo, u0)
    echo "inverse[", name, "]: inverse-forward2=", invErr,
      " forward-inverse2=", physInvErr
    doAssert invErr < 2e-16
    doAssert physInvErr < 2e-16
    checkGrad(name, action, 8e-6)

  block:
    var dense = baseParams("link2", "scalar", "stout")
    dense.flowDepth = 2
    dense.mapDirs = @[0, 1, 0, 1]
    dense.mapParities = @[0, 0, 1, 1]
    let
      spec = buildMapSpec(dense, beta)
      layout = buildMapLayout(v0[0], dense, spec)
      pair = mapHost(v0, spec, layout)
      qex = smearedField(gauge.toGvalue(grt, v0), spec.mapRho,
        spec.flowDepth)
    discard qex.smeared.eval
    let
      qexGauge = qex.smeared.gaugeSnapshot
      qexLogj = qex.lndet.eval.sval
      vg = gauge.toGvalue(grt, v0)
      pairSeff = mapAction(grt, gc, spec, layout)(vg)
      qexSeff = stoutAction(gc, spec.mapRho, spec.flowDepth).action(vg)
      pairForce = contractProjTAH(grad(pairSeff, vg), vg)
      qexForce = contractProjTAH(grad(qexSeff, vg), vg)
      forceErr = (pairForce-qexForce).norm2.eval.sval/
        (1.0+qexForce.norm2.eval.sval)
    echo "stout analytic/QEX dense: gauge2=", maxLinkDiff2(pair.u, qexGauge),
      " logdet=", abs(pair.lndet-qexLogj),
      " seff=", abs(pairSeff.eval.sval-qexSeff.eval.sval),
      " force2=", forceErr
    doAssert maxLinkDiff2(pair.u, qexGauge) < 2e-28
    doAssert abs(pair.lndet-qexLogj) < 2e-12
    doAssert abs(pairSeff.eval.sval-qexSeff.eval.sval) < 2e-11
    doAssert forceErr < 2e-20

  block:
    let
      p = baseParams("plaq4", "scalar", "identity")
      spec = buildMapSpec(p, beta)
      layout = buildMapLayout(v0[0], p, spec)
      vg = gauge.toGvalue(grt, v0)
      sf = contractProjTAH(grad(mapAction(grt, gc, spec, layout)(vg), vg), vg)
      wf = gaugeForce(gc, vg)
    doAssert (sf-wf).norm2.eval.sval < 1e-20*(1.0+wf.norm2.eval.sval)

  globalJacobianTests()
  echo "ft2du1 analytic maps OK"

when isMainModule:
  qexInit()
  runPureMapTests()
  latticeTests()
  qexFinalize()
