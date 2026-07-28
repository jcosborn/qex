## Map selection, parameter validation, construction, and lattice layout.

import qex
import math
import std/strutils

import flow
export flow

import ../graph/core
import ../graph/scalar
import ../graph/gauge
import ../graph/gauge/shared

type
  MapGeometry* = enum
    mgPlaq4, mgLink2, mgBlock5

  MapConstruction* = enum
    mcScalar, mcChain, mcCoupling

  MapBasis* = enum
    mbIdentity, mbSine, mbSqFourier, mbCSpline, mbBSpline, mbFejer, mbStout

  MapParams* = object
    geometry*, construction*, basis*: string
    ctxBasis*: string
    mapDepth*, flowDepth*: int
    # mapStrengths broadcasts or has mapDepth entries. Scalar/link parameter blocks
    # may likewise be shared or concatenated in stage order.
    mapStrengths*: seq[float]
    mapFloor*, mapEpsilon*: float
    mapRho*: float
    mapOrder*, ctxOrder*, mapKnots*, ctxKnots*: int
    # Conditional Fourier rows store [a0, a(k=0,r=1..R), ...].
    # Scalar spline controls are log densities. Conditional B-spline controls
    # use Fourier context rows or a direct positive KxctxKnots tensor.
    mapCoeffs*, mapControls*: seq[float]
    mapParamFile*: string
    fejerOrders*: seq[int]
    fejerCenters*, fejerWeights*: seq[float]
    mapStageOrder*: string
    mapDirs*, mapParities*, mapOffsets*: seq[int]
    mapStride*: int
    mapInvTol*, mapScanStep*: float
    mapInvIter*, mapScan*: int
    mapDump*: string
    mapFunctionScan*: int
    mapFunctionDump*: string
    mapFunctionContexts*: seq[float]
    monitorEvery*: int
    checkMap*: bool
    startWidth*: float

  MapSpec* = object
    geometry*: MapGeometry
    construction*: MapConstruction
    basis*: MapBasis
    beta*: float
    mapRho*: float
    flowDepth*: int
    invTol*, scanStep*: float
    invIter*: int
    chainOrder*: array[4, int]
    circle*: CircleMap
    pair*: PairMap
    blockMap*: BlockMap

  LocalEval* = object
    np*, nv*: int
    auxiliary*, physical*: Vec5
    logj*, sw*, seff*: float
    force*: Vec4

  ScanPoint* = object
    x*, y*: float
    sw*, dsw*, logj*, logjPerVar*: float
    seff*, dseff*, dseffPerPlaq*: float
    force2*, stiffMax*, stiffRms*, curvMin*, curvMax*: float
    invErr*, hessAsym*: float
    auxiliary*, physical*: Vec5
    hessian*: flow.Mat4

  MapScan* = object
    minLogj*, maxLogj*: float
    minDseff*, maxDseff*: float
    minDseffPerPlaq*, maxDseffPerPlaq*: float
    maxForce2*, maxStiffness*, minCurvature*, maxCurvature*: float
    maxInverseError*, maxHessAsym*: float

  MapLayout* = object
    circleMasks*: seq[PortalMask]
    pairLayers*: seq[PairLayer]
    blockMask*: PortalMask

proc geometryName*(x: MapGeometry): string =
  case x
  of mgPlaq4: "plaq4"
  of mgLink2: "link2"
  of mgBlock5: "block5"

proc constructionName*(x: MapConstruction): string =
  case x
  of mcScalar: "scalar"
  of mcChain: "chain"
  of mcCoupling: "coupling"

proc basisName*(x: MapBasis): string =
  case x
  of mbIdentity: "identity"
  of mbSine: "sine"
  of mbSqFourier: "sqfourier"
  of mbCSpline: "cspline"
  of mbBSpline: "bspline"
  of mbFejer: "fejer"
  of mbStout: "stout"

proc parseGeometry(s: string): MapGeometry =
  case s.toLowerAscii
  of "plaq4": mgPlaq4
  of "link2": mgLink2
  of "block5": mgBlock5
  else: mapFail("unknown map geometry: " & s)

proc parseConstruction(s: string): MapConstruction =
  case s.toLowerAscii
  of "scalar": mcScalar
  of "chain": mcChain
  of "coupling": mcCoupling
  else: mapFail("unknown map construction: " & s)

proc parseBasis(s: string): MapBasis =
  case s.toLowerAscii
  of "identity": mbIdentity
  of "sine": mbSine
  of "sqfourier": mbSqFourier
  of "cspline": mbCSpline
  of "bspline": mbBSpline
  of "fejer": mbFejer
  of "stout": mbStout
  else: mapFail("unknown map basis: " & s)

proc splineContext(p: MapParams): bool =
  case p.ctxBasis.toLowerAscii
  of "", "fourier": false
  of "bspline": true
  else: mapFail("ctxBasis must be fourier or bspline")

proc readNumberFile(path: string): seq[float] =
  for word in readFile(path).splitWhitespace:
    result.add parseFloat(word)

proc readMapParams*(d: MapParams): MapParams =
  letParam:
    geometry = d.geometry
    construction = d.construction
    basis = d.basis
    ctxBasis = d.ctxBasis
    mapDepth = d.mapDepth
    flowDepth = d.flowDepth
    mapStrengths = d.mapStrengths
    mapFloor = d.mapFloor
    mapEpsilon = d.mapEpsilon
    mapRho = d.mapRho
    mapOrder = d.mapOrder
    ctxOrder = d.ctxOrder
    mapKnots = d.mapKnots
    ctxKnots = d.ctxKnots
    mapCoeffs = d.mapCoeffs
    mapControls = d.mapControls
    mapParamFile = d.mapParamFile
    fejerOrders = d.fejerOrders
    fejerCenters = d.fejerCenters
    fejerWeights = d.fejerWeights
    mapStageOrder = d.mapStageOrder
    mapDirs = d.mapDirs
    mapParities = d.mapParities
    mapOffsets = d.mapOffsets
    mapStride = d.mapStride
    mapInvTol = d.mapInvTol
    mapScanStep = d.mapScanStep
    mapInvIter = d.mapInvIter
    mapScan = d.mapScan
    mapDump = d.mapDump
    mapFunctionScan = d.mapFunctionScan
    mapFunctionDump = d.mapFunctionDump
    mapFunctionContexts = d.mapFunctionContexts
    monitorEvery = d.monitorEvery
    checkMap = d.checkMap
    startWidth = d.startWidth
  result = MapParams(geometry: geometry, construction: construction, basis: basis,
    ctxBasis: ctxBasis,
    mapDepth: mapDepth, flowDepth: flowDepth, mapStrengths: mapStrengths,
    mapFloor: mapFloor, mapEpsilon: mapEpsilon, mapRho: mapRho, mapOrder: mapOrder,
    ctxOrder: ctxOrder, mapKnots: mapKnots, ctxKnots: ctxKnots,
    mapCoeffs: mapCoeffs,
    mapControls: mapControls, mapParamFile: mapParamFile,
    fejerOrders: fejerOrders, fejerCenters: fejerCenters,
    fejerWeights: fejerWeights, mapStageOrder: mapStageOrder,
    mapDirs: mapDirs, mapParities: mapParities, mapOffsets: mapOffsets,
    mapStride: mapStride, mapInvTol: mapInvTol, mapScanStep: mapScanStep,
    mapInvIter: mapInvIter,
    mapScan: mapScan, mapDump: mapDump, mapFunctionScan: mapFunctionScan,
    mapFunctionDump: mapFunctionDump, mapFunctionContexts: mapFunctionContexts,
    monitorEvery: monitorEvery,
    checkMap: checkMap, startWidth: startWidth)
  if result.mapParamFile.len > 0:
    let v = readNumberFile(result.mapParamFile)
    case parseBasis(result.basis)
    of mbSine, mbSqFourier:
      if result.mapCoeffs.len > 0: mapFail("use either mapCoeffs or mapParamFile, not both")
      result.mapCoeffs = v
    of mbCSpline, mbBSpline:
      if result.mapControls.len > 0: mapFail("use either mapControls or mapParamFile, not both")
      result.mapControls = v
    of mbFejer:
      if result.fejerWeights.len > 0: mapFail("use either fejerWeights or mapParamFile, not both")
      result.fejerWeights = v
    of mbIdentity, mbStout:
      mapFail(basisName(parseBasis(result.basis)) & " basis does not accept a parameter file")

proc stageStrengths(p: MapParams): seq[float] =
  if p.mapDepth < 1: mapFail("mapDepth must be positive")
  let s = if p.mapStrengths.len == 0: @[1.0] else: p.mapStrengths
  if s.len notin [1, p.mapDepth]:
    mapFail("mapStrengths must have length one or mapDepth")
  result = newSeq[float](p.mapDepth)
  for i in 0..<p.mapDepth:
    result[i] = s[if s.len == 1: 0 else: i]
    if result[i] < 0.0 or result[i] > 1.0:
      mapFail("map strengths must lie in [0,1]")

proc stageValues(v: seq[float]; unit, depth, stage: int; what: string): seq[float] =
  if v.len == unit:
    return v
  if v.len != depth*unit:
    mapFail(what & " must contain one shared block or mapDepth stage blocks")
  result = newSeq[float](unit)
  for i in 0..<unit:
    result[i] = v[stage*unit+i]

proc defaultCircleValues(p: MapParams; beta: float): seq[float] =
  let n = p.mapOrder
  case parseBasis(p.basis)
  of mbSine:
    result = newSeq[float](n)
    result[0] = -min(0.35, 0.08*beta)
  of mbSqFourier:
    result = newSeq[float](n)
    result[0] = -0.6
    if n > 1: result[1] = 0.2
    if n > 2: result[2] = -0.1
  of mbCSpline, mbBSpline:
    let knots = p.mapKnots
    result = newSeq[float](knots)
    for j in 0..<knots:
      result[j] = -min(beta, 3.0)*cos(-PI+TAU*float(j)/float(knots))
  else: discard

proc buildCircle(p: MapParams; beta: float): CircleMap =
  let
    basis = parseBasis(p.basis)
    ss = stageStrengths(p)
  case basis
  of mbSine, mbSqFourier:
    if p.mapOrder < 1: mapFail("mapOrder must be positive for Fourier bases")
    if p.mapCoeffs.len > 0 and p.mapCoeffs.len notin [p.mapOrder, p.mapDepth*p.mapOrder]:
      mapFail("scalar Fourier mapCoeffs must contain one shared block or mapDepth stage blocks")
  of mbCSpline:
    if p.mapKnots < 6 or (p.mapKnots and 1) != 0:
      mapFail("cspline mapKnots must be even and at least six")
    if p.mapControls.len > 0 and p.mapControls.len notin [p.mapKnots, p.mapDepth*p.mapKnots]:
      mapFail("scalar cspline mapControls must contain one shared block or mapDepth stage blocks")
  of mbBSpline:
    if p.mapKnots < 4: mapFail("bspline mapKnots must be at least four")
    if p.mapControls.len > 0 and p.mapControls.len notin [p.mapKnots, p.mapDepth*p.mapKnots]:
      mapFail("scalar bspline mapControls must contain one shared block or mapDepth stage blocks")
  of mbFejer:
    if p.fejerOrders.len == 0 and p.mapOrder < 2:
      mapFail("mapOrder must be at least two for the default Fejer kernel")
  of mbIdentity: discard
  of mbStout:
    mapFail("stout basis requires geometry=link2 construction=scalar")
  var stages: seq[CircleMap]
  for stage, strength in ss:
    case basis
    of mbIdentity:
      stages.add identityCircle()
    of mbSine:
      let all = if p.mapCoeffs.len > 0: p.mapCoeffs else: defaultCircleValues(p, beta)
      let c = stageValues(all, p.mapOrder, p.mapDepth, stage, "scalar Fourier mapCoeffs")
      stages.add sineCircle(c, strength, p.mapFloor)
    of mbSqFourier:
      let all = if p.mapCoeffs.len > 0: p.mapCoeffs else: defaultCircleValues(p, beta)
      let c = stageValues(all, p.mapOrder, p.mapDepth, stage, "scalar Fourier mapCoeffs")
      stages.add sqFourierCircle(c, p.mapEpsilon, strength, p.mapFloor)
    of mbCSpline:
      let all = if p.mapControls.len > 0: p.mapControls else: defaultCircleValues(p, beta)
      let c = stageValues(all, p.mapKnots, p.mapDepth, stage, "scalar cspline mapControls")
      stages.add cSplineCircle(c, strength, p.mapFloor)
    of mbBSpline:
      let all = if p.mapControls.len > 0: p.mapControls else: defaultCircleValues(p, beta)
      let c = stageValues(all, p.mapKnots, p.mapDepth, stage, "scalar bspline mapControls")
      stages.add bSplineCircle(c, p.mapFloor, strength, p.mapFloor)
    of mbFejer:
      let
        orders = if p.fejerOrders.len > 0: p.fejerOrders else: @[p.mapOrder]
        centers = p.fejerCenters
        weights = if p.fejerWeights.len > 0: p.fejerWeights else: @[1.0]
      stages.add fejerCircle(orders, centers, weights,
        strength*(1.0-p.mapFloor), p.mapFloor)
    of mbStout:
      discard
  composeCircle(stages)

proc defaultPairCoeffs(p: MapParams): seq[float] =
  let
    n = p.mapOrder
    r = p.ctxOrder
    nf = featureCount(1, r)
  result = newSeq[float](n*nf)
  let a = if parseBasis(p.basis) == mbSine: -0.24 else: -0.65
  if r > 0: result[1] = a
  elif n > 1: result[nf] = a
  if parseBasis(p.basis) == mbSqFourier and n > 1 and r > 0:
    result[nf] = 0.12

proc defaultPairControls(p: MapParams; beta: float): seq[float] =
  let
    knots = p.mapKnots
    half = knots div 2
    r = p.ctxOrder
    nf = featureCount(1, r)
  if (knots and 1) != 0: mapFail("twisted pair B-spline needs an even mapKnots")
  result = newSeq[float](half*nf)
  for j in 0..<half:
    let
      u = -PI+TAU*float(j)/float(knots)
      a = min(0.4, 0.06*beta)
    result[j*nf] = 1.0
    if r > 0: result[j*nf+1] = a*cos(u)
    else: result[j*nf] += a*cos(2.0*u)

proc defaultPairTensorControls(p: MapParams; beta: float): seq[float] =
  let
    knots = p.mapKnots
    ctxKnots = p.ctxKnots
    a = min(0.4, 0.06*beta)
  if (knots and 1) != 0 or (ctxKnots and 1) != 0:
    mapFail("twisted pair tensor B-spline needs even active and context knot counts")
  result = newSeq[float](knots*ctxKnots)
  for i in 0..<knots:
    let u = -PI+TAU*float(i)/float(knots)
    for j in 0..<ctxKnots:
      let m = -PI+TAU*float(j)/float(ctxKnots)
      result[i*ctxKnots+j] = 1.0+a*cos(u)*cos(m)

proc buildPair(p: MapParams; beta: float): PairMap =
  let
    basis = parseBasis(p.basis)
    order = p.mapOrder
    ctxOrder = p.ctxOrder
    tensor = basis == mbBSpline and splineContext(p)
  if basis == mbStout:
    if p.mapDepth != 1:
      mapFail("stout has one fixed local stage; use flowDepth for full stout sweeps")
    if p.mapCoeffs.len > 0 or p.mapControls.len > 0 or p.mapParamFile.len > 0:
      mapFail("stout accepts mapRho, not map coefficients or controls")
    if abs(p.mapRho) >= 0.5:
      mapFail("stout requires abs(mapRho) < 1/2 for a positive link Jacobian")
    if 1.0-2.0*abs(p.mapRho) < p.mapFloor:
      mapFail("stout link Jacobian is below mapFloor")
    result = PairMap(beta: beta, invTol: p.mapInvTol, invIter: p.mapInvIter)
    # QEX: g(u;m)=u-2 rho cos(m) sin(u), J=1-rho(cos p+ + cos p-).
    result.stages.add sineContextMap(1, 1, 1, [0.0, -2.0*p.mapRho],
      1.0, p.mapFloor, twisted = true)
    return
  let ss = stageStrengths(p)
  case basis
  of mbSine, mbSqFourier:
    if ctxOrder < 0: mapFail("ctxOrder must be nonnegative")
    if order < 1: mapFail("mapOrder must be positive for Fourier bases")
    let unit = order*featureCount(1, ctxOrder)
    if p.mapCoeffs.len > 0 and p.mapCoeffs.len notin [unit, p.mapDepth*unit]:
      mapFail("link Fourier mapCoeffs must contain one shared block or mapDepth stage blocks")
  of mbBSpline:
    if p.mapKnots < 4 or (p.mapKnots and 1) != 0:
      mapFail("link2 bspline mapKnots must be even and at least four")
    if tensor and (p.ctxKnots < 4 or (p.ctxKnots and 1) != 0):
      mapFail("link2 tensor B-spline ctxKnots must be even and at least four")
    if not tensor and ctxOrder < 0:
      mapFail("ctxOrder must be nonnegative for Fourier context")
    let unit = if tensor: p.mapKnots*p.ctxKnots
      else: (p.mapKnots div 2)*featureCount(1, ctxOrder)
    if p.mapControls.len > 0 and p.mapControls.len notin [unit, p.mapDepth*unit]:
      mapFail("link B-spline mapControls must contain one shared block or mapDepth stage blocks")
  of mbFejer:
    if p.fejerOrders.len == 0 and order < 2:
      mapFail("mapOrder must be at least two for the default Fejer kernel")
  of mbIdentity, mbCSpline: discard
  of mbStout: discard
  result = PairMap(beta: beta, invTol: p.mapInvTol, invIter: p.mapInvIter)
  for stage, strength in ss:
    case basis
    of mbIdentity:
      result.stages.add identityContextMap(1)
    of mbSine:
      let all = if p.mapCoeffs.len > 0: p.mapCoeffs else: defaultPairCoeffs(p)
      let c = stageValues(all, order*featureCount(1, ctxOrder), p.mapDepth,
        stage, "link Fourier mapCoeffs")
      result.stages.add sineContextMap(1, order, ctxOrder, c,
        strength, p.mapFloor, twisted = true)
    of mbSqFourier:
      let all = if p.mapCoeffs.len > 0: p.mapCoeffs else: defaultPairCoeffs(p)
      let c = stageValues(all, order*featureCount(1, ctxOrder), p.mapDepth,
        stage, "link Fourier mapCoeffs")
      result.stages.add sqFourierContextMap(1, order, ctxOrder, c,
        p.mapEpsilon, strength, p.mapFloor, twisted = true)
    of mbBSpline:
      if tensor:
        let
          all = if p.mapControls.len > 0: p.mapControls
            else: defaultPairTensorControls(p, beta)
          c = stageValues(all, p.mapKnots*p.ctxKnots, p.mapDepth, stage,
            "link tensor B-spline mapControls")
        result.stages.add tensorBSplineContextMap(1, p.mapKnots, p.ctxKnots,
          c, strength, p.mapFloor, twisted = true)
      else:
        let
          all = if p.mapControls.len > 0: p.mapControls else: defaultPairControls(p, beta)
          c = stageValues(all, (p.mapKnots div 2)*featureCount(1, ctxOrder),
            p.mapDepth, stage, "link B-spline mapControls")
        result.stages.add bSplineContextMap(1, ctxOrder, p.mapKnots, c,
          strength, p.mapFloor, twisted = true)
    of mbFejer:
      let
        orders = if p.fejerOrders.len > 0: p.fejerOrders else: @[p.mapOrder]
        centers = p.fejerCenters
        weights = if p.fejerWeights.len > 0: p.fejerWeights else: @[1.0]
      result.stages.add fejerContextMap(orders, centers, weights,
        strength*(1.0-p.mapFloor), p.mapFloor)
    of mbCSpline:
      mapFail("cspline is context-free; use bspline for link2 and block chains")
    of mbStout:
      discard

proc parseBlockOrder(s: string): array[4, int] =
  let t = if s.len == 0: "0123" else: s
  if t.len != 4: mapFail("mapStageOrder must be a permutation of 0123")
  var seen: array[4, bool]
  for i in 0..3:
    let j = ord(t[i])-ord('0')
    if j notin 0..3 or seen[j]: mapFail("mapStageOrder must be a permutation of 0123")
    seen[j] = true
    result[i] = j

proc blockValueSlice(v: seq[float]; unit, depth, stage, coord: int): seq[float] =
  if v.len == unit: return v
  var off = 0
  if v.len == depth*unit:
    off = stage*unit
  elif v.len == 4*depth*unit:
    off = (4*stage+coord)*unit
  else:
    mapFail("block coupling parameters must contain one shared block, one block per depth, or four coordinate blocks per depth")
  result = newSeq[float](unit)
  for i in 0..<unit: result[i] = v[off+i]

proc defaultBlockCoeffs(p: MapParams): seq[float] =
  let
    order = p.mapOrder
    ctxOrder = p.ctxOrder
    nf = featureCount(4, ctxOrder)
  result = newSeq[float](order*nf)
  result[0] = if parseBasis(p.basis) == mbSine: -0.16 else: -0.5
  if ctxOrder > 0:
    for k in 0..<4: result[1+k*ctxOrder] = 0.06

proc defaultBlockControls(p: MapParams): seq[float] =
  let
    knots = p.mapKnots
    ctxOrder = p.ctxOrder
    nf = featureCount(4, ctxOrder, phase = true)
  result = newSeq[float](knots*nf)
  for j in 0..<knots:
    let u = -PI+TAU*float(j)/float(knots)
    result[j*nf] = 1.0-0.2*cos(u)
    if ctxOrder > 0:
      result[j*nf+1] = -0.08*cos(u)

proc defaultBlockTensorControls(p: MapParams): seq[float] =
  result = newSeq[float](p.mapKnots*p.ctxKnots)
  for i in 0..<p.mapKnots:
    let u = -PI+TAU*float(i)/float(p.mapKnots)
    for j in 0..<p.ctxKnots:
      let c = -PI+TAU*float(j)/float(p.ctxKnots)
      result[i*p.ctxKnots+j] = 1.0-0.2*cos(u)-0.08*cos(u)*cos(c)

proc buildBlock(p: MapParams; beta: float): BlockMap =
  let
    basis = parseBasis(p.basis)
    order = p.mapOrder
    ctxOrder = p.ctxOrder
    ss = stageStrengths(p)
    tensor = basis == mbBSpline and splineContext(p)
  case basis
  of mbSine, mbSqFourier:
    if ctxOrder < 0: mapFail("ctxOrder must be nonnegative")
    if order < 1: mapFail("mapOrder must be positive for Fourier bases")
  of mbBSpline:
    if p.mapKnots < 4: mapFail("block5 bspline mapKnots must be at least four")
    if tensor and p.ctxKnots < 4:
      mapFail("block5 tensor B-spline ctxKnots must be at least four")
    if not tensor and ctxOrder < 0:
      mapFail("ctxOrder must be nonnegative for Fourier context")
  of mbIdentity, mbCSpline, mbFejer: discard
  of mbStout:
    mapFail("stout basis requires geometry=link2 construction=scalar")
  result = BlockMap(beta: beta, order: parseBlockOrder(p.mapStageOrder),
    invTol: p.mapInvTol, invIter: p.mapInvIter)
  for coord in 0..3:
    for stage, strength in ss:
      case basis
      of mbIdentity:
        result.stages[coord].add identityContextMap(4)
      of mbSine:
        let nf = featureCount(4, ctxOrder)
        let all = if p.mapCoeffs.len > 0: p.mapCoeffs else: defaultBlockCoeffs(p)
        result.stages[coord].add sineContextMap(4, order, ctxOrder,
          blockValueSlice(all, order*nf, p.mapDepth, stage, coord),
          strength, p.mapFloor)
      of mbSqFourier:
        let nf = featureCount(4, ctxOrder)
        let all = if p.mapCoeffs.len > 0: p.mapCoeffs else: defaultBlockCoeffs(p)
        result.stages[coord].add sqFourierContextMap(4, order, ctxOrder,
          blockValueSlice(all, order*nf, p.mapDepth, stage, coord),
          p.mapEpsilon, strength, p.mapFloor)
      of mbBSpline:
        if tensor:
          let
            unit = p.mapKnots*p.ctxKnots
            all = if p.mapControls.len > 0: p.mapControls
              else: defaultBlockTensorControls(p)
          result.stages[coord].add tensorBSplineContextMap(4, p.mapKnots,
            p.ctxKnots, blockValueSlice(all, unit, p.mapDepth, stage, coord),
            strength, p.mapFloor, phase = true)
        else:
          let
            knots = p.mapKnots
            nf = featureCount(4, ctxOrder, phase = true)
            all = if p.mapControls.len > 0: p.mapControls else: defaultBlockControls(p)
          result.stages[coord].add bSplineContextMap(4, ctxOrder, knots,
            blockValueSlice(all, knots*nf, p.mapDepth, stage, coord),
            strength, p.mapFloor, phase = true)
      of mbCSpline, mbFejer, mbStout:
        mapFail("block5 coupling supports identity, sine, sqfourier, and bspline")

proc buildMapSpec*(p: MapParams; beta: float): MapSpec =
  if beta <= 0.0: mapFail("beta must be positive")
  if p.flowDepth < 1: mapFail("flowDepth must be positive")
  if p.mapInvTol <= 0.0 or p.mapInvIter < 1:
    mapFail("invalid inverse solver controls")
  if p.mapScan < 8: mapFail("mapScan must be at least eight")
  if p.mapScanStep <= 0.0: mapFail("mapScanStep must be positive")
  if p.mapStride < 0: mapFail("mapStride must be nonnegative")
  if p.monitorEvery < 0: mapFail("monitorEvery must be nonnegative")
  result = MapSpec(geometry: parseGeometry(p.geometry),
    construction: parseConstruction(p.construction), basis: parseBasis(p.basis),
    beta: beta, mapRho: p.mapRho, flowDepth: p.flowDepth, invTol: p.mapInvTol,
    scanStep: p.mapScanStep, invIter: p.mapInvIter)
  if result.basis != mbIdentity and
      (p.mapFloor < mapDerivativeFloor or p.mapFloor >= 1.0):
    mapFail("map derivative floor is outside the safe range")
  case result.geometry
  of mgPlaq4:
    if result.construction != mcScalar:
      mapFail("plaq4 supports only construction=scalar")
    result.circle = buildCircle(p, beta)
  of mgLink2:
    if result.construction != mcScalar:
      mapFail("link2 supports only construction=scalar")
    result.pair = buildPair(p, beta)
  of mgBlock5:
    if result.basis == mbStout:
      mapFail("stout basis requires geometry=link2 construction=scalar")
    case result.construction
    of mcChain:
      result.chainOrder = parseBlockOrder(p.mapStageOrder)
      result.pair = buildPair(p, beta)
    of mcCoupling:
      if p.flowDepth != 1:
        mapFail("block5 coupling is one sparse lattice layer; mapDepth controls its internal coupling depth")
      result.blockMap = buildBlock(p, beta)
    of mcScalar:
      mapFail("block5 supports construction=chain or construction=coupling")

proc mapOffset(p: MapParams): tuple[x, y: int] =
  if p.mapOffsets.len == 0: return (0, 0)
  if p.mapOffsets.len != 2: mapFail("mapOffsets must contain x,y")
  (p.mapOffsets[0], p.mapOffsets[1])

proc parityMasks(proto: PortalMask; parities: seq[int]): seq[PortalMask] =
  let both = portalMasksEvenOdd(proto)
  let ps = if parities.len == 0: @[0, 1] else: parities
  for p in ps:
    if p notin 0..1: mapFail("map parities must be zero or one")
    result.add both[p]

proc buildMapLayout*(proto: PortalMask; p: MapParams; spec: MapSpec): MapLayout =
  let off = mapOffset(p)
  case spec.geometry
  of mgPlaq4:
    if p.mapStride > 0:
      result.circleMasks = portalMasksStride(proto, p.mapStride, off.x, off.y)
    else:
      result.circleMasks = parityMasks(proto, p.mapParities)
  of mgLink2:
    let dirs = if p.mapDirs.len == 0: @[0] else: p.mapDirs
    if p.mapStride > 0:
      if p.mapParities.len > 1:
        mapFail("strided link2 layers do not use multiple parity entries")
      let mask = portalMasksStride(proto, p.mapStride, off.x, off.y)[0]
      for dir in dirs:
        if dir notin 0..1: mapFail("map directions must be zero or one")
        result.pairLayers.add PairLayer(mask: mask, dir: dir)
    else:
      let
        ps = if p.mapParities.len == 0: @[0, 1] else: p.mapParities
        n = max(dirs.len, ps.len)
      if dirs.len notin [1, n] or ps.len notin [1, n]:
        mapFail("mapDirs and mapParities must broadcast or have equal length")
      let both = portalMasksEvenOdd(proto)
      for i in 0..<n:
        let
          dir = dirs[if dirs.len == 1: 0 else: i]
          parity = ps[if ps.len == 1: 0 else: i]
        if dir notin 0..1 or parity notin 0..1:
          mapFail("map directions and parities must be zero or one")
        result.pairLayers.add PairLayer(mask: both[parity], dir: dir)
  of mgBlock5:
    let stride = if p.mapStride > 0: p.mapStride else: 4
    if stride < 4: mapFail("block5 needs mapStride >= 4")
    case spec.construction
    of mcChain:
      const
        dirs = [0, 1, 0, 1]
        dx = [0, 1, 0, 0]
        dy = [0, 0, 1, 0]
      let order = parseBlockOrder(p.mapStageOrder)
      for j in order:
        let mask = portalMasksStride(proto, stride, off.x+dx[j], off.y+dy[j])[0]
        result.pairLayers.add PairLayer(mask: mask, dir: dirs[j])
    of mcCoupling:
      result.blockMask = portalMasksStride(proto, stride, off.x, off.y)[0]
    of mcScalar: discard

proc mapAction*(grt: GraphRuntime; gc: Gactcoeff; spec: MapSpec;
                layout: MapLayout): proc(V: Ggauge): Gscalar =
  case spec.geometry
  of mgPlaq4:
    onePlaqAction(grt, gc, spec.circle, layout.circleMasks, spec.flowDepth)
  of mgLink2:
    pairAction(gc, spec.pair, layout.pairLayers, spec.flowDepth)
  of mgBlock5:
    case spec.construction
    of mcChain: pairAction(gc, spec.pair, layout.pairLayers, spec.flowDepth)
    of mcCoupling: blockAction(gc, layout.blockMask, spec.blockMap)
    of mcScalar: mapFail("invalid block5 scalar construction")

proc mapHost*(V: shared.Gauge; spec: MapSpec; layout: MapLayout): tuple[u: shared.Gauge, lndet: float] =
  case spec.geometry
  of mgPlaq4:
    smearFlowHost(V, spec.circle, layout.circleMasks, spec.flowDepth)
  of mgLink2:
    pairFlowHost(V, spec.pair, layout.pairLayers, spec.flowDepth)
  of mgBlock5:
    case spec.construction
    of mcChain: pairFlowHost(V, spec.pair, layout.pairLayers, spec.flowDepth)
    of mcCoupling: blockSmearHost(V, layout.blockMask, spec.blockMap)
    of mcScalar: mapFail("invalid block5 scalar construction")

proc invertMapFlow*(U: auto; spec: MapSpec; layout: MapLayout) =
  case spec.geometry
  of mgPlaq4:
    invertPortalFlow(U, spec.circle, layout.circleMasks, spec.flowDepth,
      spec.invTol, spec.invIter)
  of mgLink2:
    invertPairFlow(U, spec.pair, layout.pairLayers, spec.flowDepth)
  of mgBlock5:
    case spec.construction
    of mcChain: invertPairFlow(U, spec.pair, layout.pairLayers, spec.flowDepth)
    of mcCoupling: invertBlockSmear(U, layout.blockMask, spec.blockMap)
    of mcScalar: mapFail("invalid block5 scalar construction")
