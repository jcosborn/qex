## Native traces of S = -k sum cos(theta), with angular force k sin(theta).
import qex
import ../../[core, scalar, gauge, plan]
import ../../hmcgauge/[integrator, trajectory, config]
import std/[math, sequtils]

type
  IntegratorSample* = object
    coefficient*: float
    theta*, gaugeRe*, gaugeIm*, momentum*, momentumRe*, force*: seq[float]
  IntegratorCase* = object
    dt*, lambda*, kappa*: float
    theta0*, momentum0*: seq[float]
    events*: seq[IntegratorSample]
    initialHamiltonian*, finalHamiltonian*: float
    graph*: TrajectoryGraph

proc values(x: Ggauge; imaginary: bool): seq[float] =
  discard x.eval
  let fs = x.gval
  let lo = fs[0].l
  result = newSeq[float](fs.len*lo.physVol)
  for d in 0..<fs.len:
    for s in 0..<lo.nSites:
      let r = lo.coords[0][s].int
      let c = lo.coords[1][s].int
      let i = (d*lo.physGeom[0]+r)*lo.physGeom[1]+c
      if imaginary: result[i] := fs[d]{s}[0,0].im
      else: result[i] := fs[d]{s}[0,0].re
  lo.comm.rankSum(result)

proc sampleIntegrator*(n: int; planned = true; kind = ik2MNp): IntegratorCase =
  let lo = newLayout(@[8,8])
  let g = lo.newGauge
  let p = lo.newGauge
  result.dt = 0.35
  result.kappa = 0.47
  result.theta0 = newSeq[float](2*lo.physVol)
  result.momentum0 = newSeq[float](2*lo.physVol)
  for d in 0..1:
    for s in 0..<lo.nSites:
      let r = lo.coords[0][s].int
      let c = lo.coords[1][s].int
      let i = (d*lo.physGeom[0]+r)*lo.physGeom[1]+c
      let a = 0.11+0.20*float(r)-0.17*float(c)+0.31*float(d)
      let b = 1.9+0.09*float(r)-0.11*float(c)+0.4*sin(0.29*float(d)+0.13*float(r))
      result.theta0[i] = a
      result.momentum0[i] = b
      g[d]{s}[0,0].re := cos(a)
      g[d]{s}[0,0].im := sin(a)
      p[d]{s}[0,0].re := 0.0
      p[d]{s}[0,0].im := b
  lo.comm.rankSum(result.theta0)
  lo.comm.rankSum(result.momentum0)
  let rt = initGraphRuntime()
  let kappa = result.kappa
  let coeffs = parseIntegratorCoeffs(kind,[])
  result.lambda = coeffs.lambda
  let action: GaugeAction = proc(x: Ggauge): Gscalar =
    scalar.toGvalue(x.runtime,-kappa)*retr(x)
  let conf = RunConfig(dt:result.dt,gsteps:n,integratorCoeffs:coeffs)
  let graph = buildTrajectoryGraph(rt,g,p,action,conf,buildTraining=false,trace=true)
  result.graph = graph
  var roots = @[Gvalue(graph.initialState.hamiltonian), Gvalue(graph.finalState.hamiltonian)]
  for event in graph.trace:
    roots.add event.gauge
    roots.add event.momentum
    if event.force != nil: roots.add event.force
    if event.coefficient != nil: roots.add event.coefficient
  var evals: GraphPlan
  defer:
    if evals != nil: evals.clear
  var outputs: seq[Gvalue]
  if planned:
    evals = plan(roots)
    outputs = evals.eval
  else:
    for root in roots: outputs.add root.eval
  var k = 2
  var theta = result.theta0.mapIt(it)
  for event in graph.trace:
    var record: IntegratorSample
    let gauge = Ggauge(outputs[k])
    let momentum = Ggauge(outputs[k+1])
    k += 2
    if event.force != nil:
      record.force = values(Ggauge(outputs[k]),true)
      inc k
    if event.coefficient != nil:
      record.coefficient = Gscalar(outputs[k]).sval
      inc k
    record.momentum = values(momentum,true)
    if event.kind == ieDrift:
      for i in 0..<theta.len: theta[i] += record.coefficient*record.momentum[i]
    record.theta = theta.mapIt(it)
    record.gaugeRe = values(gauge,false)
    record.gaugeIm = values(gauge,true)
    record.momentumRe = values(momentum,false)
    result.events.add record
  result.initialHamiltonian = Gscalar(outputs[0]).sval
  result.finalHamiltonian = Gscalar(outputs[1]).sval
