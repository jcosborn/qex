## Local effective-action, force, Hessian, and map-function scans.

import math

import config
export config

type LocalJac = array[5, Vec4]

const
  linkInc: array[2, float] = [1.0, -1.0]
  blockInc: LocalJac = [
    [1.0, 1.0, 1.0, 1.0],
    [-1.0, 0.0, 0.0, 0.0],
    [0.0, -1.0, 0.0, 0.0],
    [0.0, 0.0, -1.0, 0.0],
    [0.0, 0.0, 0.0, -1.0]]

proc finishLocal(e: var LocalEval; d: LocalJac; dl: Vec4; beta: float) =
  # g_j=Σ_i (∂P_i/∂v_j) β sin(P_i)-∂logJ/∂v_j.
  for i in 0..<e.np: e.sw -= beta*cos(e.physical[i])
  e.seff = e.sw-e.logj
  for j in 0..<e.nv:
    e.force[j] = -dl[j]
    for i in 0..<e.np:
      e.force[j] += d[i][j]*beta*sin(e.physical[i])

proc evalLocal*(spec: MapSpec; p: Vec5): LocalEval =
  case spec.geometry
  of mgPlaq4:
    let
      c = circleEval(spec.circle, p[0])
      delta = c.y-p[0]
    var
      d: LocalJac
      dl: Vec4
    result = LocalEval(np: 5, nv: 4, auxiliary: p, logj: ln(c.dy))
    result.physical[0] = c.y
    for i in 0..3: result.physical[i+1] = p[i+1]-0.25*delta
    for j in 0..3:
      d[0][j] = c.dy
      dl[j] = c.ddy/c.dy
      for i in 0..3:
        d[i+1][j] = (if i == j: -1.0 else: 0.0)-0.25*(c.dy-1.0)
    finishLocal(result, d, dl, spec.beta)
  of mgLink2:
    let a = evalPair(spec.pair, p[0], p[1])
    var
      d: LocalJac
      dl: Vec4
    result = LocalEval(np: 2, nv: 1, auxiliary: p, logj: a.logdet)
    result.physical[0] = a.physicalPlus
    result.physical[1] = a.physicalMinus
    d[0][0] = a.jac[0][0]-a.jac[0][1]
    d[1][0] = a.jac[1][0]-a.jac[1][1]
    dl[0] = a.logdetPlus-a.logdetMinus
    finishLocal(result, d, dl, spec.beta)
  of mgBlock5:
    case spec.construction
    of mcChain:
      var
        q = p
        d = blockInc
        dl: Vec4
        logj = 0.0
      for oi in 0..3:
        let
          k = spec.chainOrder[oi]
          plus = if k < 2: 0 else: k+1
          minus = if k < 2: k+1 else: 0
          a = evalPair(spec.pair, q[plus], q[minus])
        var dp, dm: Vec4
        for j in 0..3:
          dl[j] += a.logdetPlus*d[plus][j]+a.logdetMinus*d[minus][j]
          dp[j] = a.jac[0][0]*d[plus][j]+a.jac[0][1]*d[minus][j]
          dm[j] = a.jac[1][0]*d[plus][j]+a.jac[1][1]*d[minus][j]
        d[plus] = dp
        d[minus] = dm
        q[plus] = a.physicalPlus
        q[minus] = a.physicalMinus
        logj += a.logdet
      result = LocalEval(np: 5, nv: 4, auxiliary: p, physical: q, logj: logj)
      finishLocal(result, d, dl, spec.beta)
    of mcCoupling:
      let
        z: Vec4 = [p[0], p[1], p[2], p[3]]
        flux = p[0]+p[1]+p[2]+p[3]+p[4]
        a = evalBlockMap(spec.blockMap, z, flux)
      var
        d: LocalJac
        dl: Vec4
      result = LocalEval(np: 5, nv: 4, auxiliary: p,
        physical: a.physical, logj: a.logdet)
      for j in 0..3:
        for k in 0..3:
          dl[j] += a.logdetGradient[k]*blockInc[k][j]
          for i in 0..4:
            d[i][j] += a.derivatives[i][k]*blockInc[k][j]
      finishLocal(result, d, dl, spec.beta)
    of mcScalar:
      mapFail("invalid block5 scalar construction")

proc invertLocal*(spec: MapSpec; physical: Vec5): Vec5 =
  result = physical
  case spec.geometry
  of mgPlaq4:
    let
      c = circleInv(spec.circle, physical[0], spec.invTol, spec.invIter)
      delta = physical[0]-c
    result[0] = c
    for i in 1..4: result[i] = physical[i]+0.25*delta
  of mgLink2:
    let p = invertPair(spec.pair, physical[0], physical[1])
    result[0] = p[0]
    result[1] = p[1]
  of mgBlock5:
    case spec.construction
    of mcChain:
      for oi in countdown(3, 0):
        let
          k = spec.chainOrder[oi]
          plus = if k < 2: 0 else: k+1
          minus = if k < 2: k+1 else: 0
          p = invertPair(spec.pair, result[plus], result[minus])
        result[plus] = p[0]
        result[minus] = p[1]
    of mcCoupling:
      let
        flux = physical[0]+physical[1]+physical[2]+physical[3]+physical[4]
        z = invertBlockMap(spec.blockMap,
          [physical[0], physical[1], physical[2], physical[3]], flux)
      for i in 0..3: result[i] = z[i]
      result[4] = flux-z[0]-z[1]-z[2]-z[3]
    of mcScalar:
      mapFail("invalid block5 scalar construction")

func localInc(np, i, j: int): float {.inline.} =
  if np == 2: linkInc[i]
  else: blockInc[i][j]

proc symmetricEigen(a0: Mat4; n: int): Vec4 =
  var a = a0
  for _ in 0..<24:
    for p in 0..<n:
      for q in p+1..<n:
        let apq = a[p][q]
        if abs(apq) > 1e-14*(1.0+abs(a[p][p])+abs(a[q][q])):
          let
            tau = (a[q][q]-a[p][p])/(2.0*apq)
            t = if tau >= 0.0: 1.0/(tau+sqrt(1.0+tau*tau))
                else: -1.0/(-tau+sqrt(1.0+tau*tau))
            c = 1.0/sqrt(1.0+t*t)
            s = t*c
            app = a[p][p]
            aqq = a[q][q]
          for k in 0..<n:
            if k != p and k != q:
              let
                akp = a[k][p]
                akq = a[k][q]
              a[k][p] = c*akp-s*akq
              a[p][k] = a[k][p]
              a[k][q] = s*akp+c*akq
              a[q][k] = a[k][q]
          a[p][p] = c*c*app-2.0*s*c*apq+s*s*aqq
          a[q][q] = s*s*app+2.0*s*c*apq+c*c*aqq
          a[p][q] = 0.0
          a[q][p] = 0.0
  for i in 0..<n: result[i] = a[i][i]

proc vacuumLocal(spec: MapSpec): LocalEval =
  var p: Vec5
  evalLocal(spec, invertLocal(spec, p))

proc makeScanPoint(spec: MapSpec; physical: Vec5; base: LocalEval;
                   x, y: float): ScanPoint =
  let
    aux = invertLocal(spec, physical)
    e = evalLocal(spec, aux)
    h = spec.scanStep
  var raw: Mat4
  for j in 0..<e.nv:
    var pp = aux
    var pm = aux
    for i in 0..<e.np:
      let b = localInc(e.np, i, j)
      pp[i] += h*b
      pm[i] -= h*b
    let
      gp = evalLocal(spec, pp).force
      gm = evalLocal(spec, pm).force
    for i in 0..<e.nv: raw[i][j] = (gp[i]-gm[i])/(2.0*h)

  result = ScanPoint(x: x, y: y, sw: e.sw, dsw: e.sw-base.sw,
    logj: e.logj, logjPerVar: e.logj/float(e.nv), seff: e.seff,
    dseff: e.seff-base.seff, auxiliary: aux, physical: physical)
  result.dseffPerPlaq = result.dseff/float(e.np)
  for i in 0..<e.nv:
    result.force2 += e.force[i]*e.force[i]/float(e.nv)
    result.hessian[i][i] = raw[i][i]
    for j in i+1..<e.nv:
      result.hessAsym = max(result.hessAsym, abs(raw[i][j]-raw[j][i]))
      let v = 0.5*(raw[i][j]+raw[j][i])
      result.hessian[i][j] = v
      result.hessian[j][i] = v
  for i in 0..<e.np:
    result.invErr = max(result.invErr,
      abs(principalAngle(e.physical[i]-physical[i])))
  let lam = symmetricEigen(result.hessian, e.nv)
  result.stiffMax = abs(lam[0])
  result.curvMin = lam[0]
  result.curvMax = lam[0]
  for i in 0..<e.nv:
    result.stiffMax = max(result.stiffMax, abs(lam[i]))
    result.curvMin = min(result.curvMin, lam[i])
    result.curvMax = max(result.curvMax, lam[i])
    result.stiffRms += lam[i]*lam[i]/float(e.nv)
  result.stiffRms = sqrt(result.stiffRms)

proc scanPoint*(spec: MapSpec; physical: Vec5): ScanPoint =
  makeScanPoint(spec, physical, vacuumLocal(spec), physical[0], physical[1])

proc dumpMapFunction*(spec: MapSpec; n: int; contexts: openArray[float];
                      path: string) =
  if n < 8: mapFail("mapFunctionScan must be at least eight")
  if path.len == 0: mapFail("mapFunctionDump path is empty")
  let
    nctx = case spec.geometry
      of mgPlaq4: 0
      of mgLink2: 1
      of mgBlock5: (if spec.construction == mcCoupling: 4 else: 1)
    ncoord = if spec.geometry == mgBlock5 and
      spec.construction == mcCoupling: 4 else: 1
  if nctx == 0 and contexts.len != 0:
    mapFail("plaq4 map functions do not take contexts")
  if nctx > 0 and contexts.len mod nctx != 0:
    mapFail("mapFunctionContexts must contain complete context vectors of length " & $nctx)
  let nslice = if contexts.len == 0: 1 else: contexts.len div nctx
  var f = open(path, fmWrite)
  defer: f.close()
  f.writeLine("# geometry=" & geometryName(spec.geometry) &
    " construction=" & constructionName(spec.construction) &
    " basis=" & basisName(spec.basis) & " n=" & $n &
    " intervals=" & $n & " endpoints=included")
  f.writeLine("# Each data set is a complete scalar map at fixed context; block5/coupling emits its conditional stages.")
  f.writeLine("# stage=-1 denotes the complete composed scalar map.")
  f.writeLine(case spec.geometry
    of mgPlaq4: "# context: none"
    of mgLink2: "# context: c0=mean=(Pplus+Pminus)/2"
    of mgBlock5:
      if spec.construction == mcChain:
        "# context: c0=pair mean"
      else:
        "# context: c0,c1,c2=other active coordinates in ascending order; c3=flux")
  f.writeLine("# slice coordinate stage x g delta gp gpp logGp dLogGp dc0 dc1 dc2 dc3 dxc0 dxc1 dxc2 dxc3 c0 c1 c2 c3")

  var dataset = 0
  for slice in 0..<nslice:
    var ctx: array[maxMapContext, float]
    if contexts.len > 0:
      for k in 0..<nctx: ctx[k] = contexts[slice*nctx+k]
    for coord in 0..<ncoord:
      let nstage = case spec.geometry
        of mgPlaq4, mgLink2: 1
        of mgBlock5:
          if spec.construction == mcCoupling:
            spec.blockMap.stages[coord].len
          else:
            1
      for part in 0..<nstage:
        let stage = if spec.geometry == mgBlock5 and
          spec.construction == mcCoupling: part else: -1
        if dataset > 0:
          f.writeLine("")
          f.writeLine("")
        f.writeLine("# dataset=" & $dataset & " slice=" & $slice &
          " coordinate=" & $coord & " stage=" & $stage)
        for i in 0..n:
          let x = -PI+TAU*float(i)/float(n)
          var e: ContextEval
          case spec.geometry
          of mgPlaq4:
            let q = circleEval(spec.circle, x)
            e = ContextEval(y: q.y, dx: q.dy, dxx: q.ddy)
          of mgLink2:
            e = evalContextStages(spec.pair.stages, x, ctx)
          of mgBlock5:
            if spec.construction == mcCoupling:
              e = evalContext(spec.blockMap.stages[coord][stage], x, ctx)
            else:
              e = evalContextStages(spec.pair.stages, x, ctx)
          var line = $slice & " " & $coord & " " & $stage & " " & $x &
            " " & $e.y & " " & $(e.y-x) & " " & $e.dx & " " & $e.dxx &
            " " & $ln(e.dx) & " " & $(e.dxx/e.dx)
          for k in 0..<maxMapContext: line.add " " & $e.dc[k]
          for k in 0..<maxMapContext: line.add " " & $e.dxc[k]
          for k in 0..<maxMapContext: line.add " " & $ctx[k]
          f.writeLine(line)
        inc dataset

proc scanPhysical*(spec: MapSpec; x, y: float): Vec5 =
  result[0] = x
  if spec.geometry == mgLink2:
    result[1] = y
  else:
    for i in 1..4: result[i] = y

proc scanLine(p: ScanPoint; np: int): string =
  result = $p.x & " " & $p.y & " " & $p.sw & " " & $p.dsw & " " &
    $p.logj & " " & $p.logjPerVar & " " & $p.seff & " " & $p.dseff &
    " " & $p.dseffPerPlaq & " " & $p.force2 & " " & $p.stiffMax &
    " " & $p.stiffRms & " " & $p.curvMin & " " & $p.invErr
  for i in 0..<np: result.add " " & $p.physical[i]
  for i in 0..<np: result.add " " & $p.auxiliary[i]
  result.add " " & $p.curvMax

proc mapScan*(spec: MapSpec; n = 96; path = ""): MapScan =
  if n < 8: mapFail("local map scan needs at least eight points per axis")
  let base = vacuumLocal(spec)
  var f: File
  if path.len > 0:
    f = open(path, fmWrite)
    f.writeLine("# geometry=" & geometryName(spec.geometry) &
      " construction=" & constructionName(spec.construction) &
      " basis=" & basisName(spec.basis) & " beta=" & $spec.beta &
      (if spec.basis == mbStout: " mapRho=" & $spec.mapRho else: "") &
      " n=" & $n & " mapScanStep=" & $spec.scanStep &
      " mass=1 localRounds=1 configuredFlowDepth=" & $spec.flowDepth &
      " flowDepthApplied=false")
    f.writeLine(if spec.geometry == mgLink2:
      "# slice: x=Pplus y=Pminus"
      else: "# slice: x=Pc y=P0=P1=P2=P3")
    f.writeLine(if spec.geometry == mgLink2:
      "# x y sw dsw logj logjPerVar seff dseff dseffPerPlaq force2 stiffMax stiffRms curvMin invErr Pplus Pminus pplus pminus curvMax"
      else: "# x y sw dsw logj logjPerVar seff dseff dseffPerPlaq force2 stiffMax stiffRms curvMin invErr Pc P0 P1 P2 P3 pc p0 p1 p2 p3 curvMax")
  defer:
    if path.len > 0: f.close()

  var first = true
  for ix in 0..<n:
    let x = -PI+TAU*float(ix)/float(n)
    for iy in 0..<n:
      let
        y = -PI+TAU*float(iy)/float(n)
        physical = scanPhysical(spec, x, y)
        p = makeScanPoint(spec, physical, base, x, y)
      if first:
        result = MapScan(minLogj: p.logj, maxLogj: p.logj,
          minDseff: p.dseff, maxDseff: p.dseff,
          minDseffPerPlaq: p.dseffPerPlaq,
          maxDseffPerPlaq: p.dseffPerPlaq, maxForce2: p.force2,
          maxStiffness: p.stiffMax, minCurvature: p.curvMin,
          maxCurvature: p.curvMax,
          maxInverseError: p.invErr, maxHessAsym: p.hessAsym)
        first = false
      else:
        result.minLogj = min(result.minLogj, p.logj)
        result.maxLogj = max(result.maxLogj, p.logj)
        result.minDseff = min(result.minDseff, p.dseff)
        result.maxDseff = max(result.maxDseff, p.dseff)
        result.minDseffPerPlaq = min(result.minDseffPerPlaq, p.dseffPerPlaq)
        result.maxDseffPerPlaq = max(result.maxDseffPerPlaq, p.dseffPerPlaq)
        result.maxForce2 = max(result.maxForce2, p.force2)
        result.maxStiffness = max(result.maxStiffness, p.stiffMax)
        result.minCurvature = min(result.minCurvature, p.curvMin)
        result.maxCurvature = max(result.maxCurvature, p.curvMax)
        result.maxInverseError = max(result.maxInverseError, p.invErr)
        result.maxHessAsym = max(result.maxHessAsym, p.hessAsym)
      if path.len > 0: f.writeLine(scanLine(p, base.np))
    if path.len > 0: f.writeLine("")

proc echoMapSummary*(spec: MapSpec; s: MapScan; n: int) =
  echo "map geometry=", geometryName(spec.geometry),
    " construction=", constructionName(spec.construction),
    " basis=", basisName(spec.basis),
    (if spec.basis == mbStout: " mapRho=" & $spec.mapRho else: ""),
    " flowDepth=", spec.flowDepth
  echo "local scan grid=", n, "x", n,
    " mass=1 localRounds=1 scanStep=", spec.scanStep
  echo "local logJ=[", s.minLogj, ", ", s.maxLogj,
    "] dSeff=[", s.minDseff, ", ", s.maxDseff,
    "] dSeff/plaq=[", s.minDseffPerPlaq, ", ", s.maxDseffPerPlaq, "]"
  echo "local maxForce2=", s.maxForce2,
    " maxStiffness=", s.maxStiffness,
    " curvature=[", s.minCurvature, ", ", s.maxCurvature, "]",
    " maxInverseError=", s.maxInverseError,
    " maxHessianAsymmetry=", s.maxHessAsym

# HMC application.
