## Traced integration of the separable U(1) action against its analytic dynamics.
import base/globals
setDefaultNc(1)
setVLENmax(4)
import qex
import ../[core, gauge]
import ../hmcgauge/integrator
import hmcgauge/integrator_case
import std/[math, sequtils, unittest]

addOutputFormatter(newConsoleOutputFormatter(colorOutput=false))

template checkArray(actual, expected: seq[float]; label: string) =
  # A template so `check` reports into the enclosing case; bind the arguments once.
  block:
    let a = actual
    let b = expected
    check a.len == b.len
    if a.len == b.len:
      for i in 0..<a.len:
        if abs(a[i]-b[i]) >= 2e-12:
          checkpoint(label & " differs at " & $i)
          check abs(a[i]-b[i]) < 2e-12

proc energy(theta, momentum: seq[float]; kappa: float): float =
  for i in 0..<theta.len:
    result += -kappa*cos(theta[i])+0.5*momentum[i]*momentum[i]

type Step = tuple[kind: IntegrationEventKind, coefficient: float, step: int]

proc schedule2MN(kind: IntegratorKind; n: int; lam, dt: float): seq[Step] =
  ## Both minimal-norm orderings in the integrator's coefficient arithmetic;
  ## kicks apply the negated coefficient.
  let first = lam*dt
  let between = 2.0*first
  let middle = dt-between
  let half = 0.5*dt
  let mf = kind == ik2MNp
  let outer = if mf: ieKick else: ieDrift
  let inner = if mf: ieDrift else: ieKick
  proc add(r: var seq[Step]; k: IntegrationEventKind; c: float; s: int) =
    if k == ieKick:
      r.add (ieForce, 0.0, s)
      r.add (k, -c, s)
    else:
      r.add (k, c, s)
  result.add(outer, first, 0)
  for i in 0..<n:
    result.add(inner, half, i)
    result.add(outer, middle, i)
    result.add(inner, half, i)
    if i+1 < n: result.add(outer, between, i)
  result.add(outer, first, n)

proc scheduleFG(kind: IntegratorKind; n: int): seq[Step] =
  ## Event kinds of the force-gradient schedules; the coefficients are not pinned.
  let kick = @[ieForce, ieKick]
  let gradKick = @[ieForce, ieShift, ieForce, ieKick]
  let step =
    if kind == ik4MN3F1GP: @[ieDrift] & kick & @[ieDrift] & gradKick & @[ieDrift] & kick
    else: @[ieDrift] & kick & @[ieDrift] & gradKick & @[ieDrift] & kick & @[ieDrift] & gradKick & @[ieDrift] & kick
  for i in 0..<n:
    for k in step: result.add (k, 0.0, i)
  result.add (ieDrift, 0.0, n)

qexInit()
suite "traced U(1) integration":
  for kind in [ik2MN, ik2MNp, ik4MN3F1GP, ik4MN5F2GP]:
    for planned in [false,true]:
      for n in [1,2,10]:
        test $kind & (if planned: " planned " else: " direct ") & $n & " steps follow the analytic dynamics":
          let data = sampleIntegrator(n, planned, kind)
          let graph = data.graph
          let dt = data.dt
          let lam = data.lambda
          let kappa = data.kappa
          let minimalNorm = kind in {ik2MN, ik2MNp}
          if minimalNorm: check lam == 0.1931833275037836
          let expected = if minimalNorm: schedule2MN(kind, n, lam, dt) else: scheduleFG(kind, n)
          check graph.trace.len == expected.len
          let nforces = expected.countIt(it.kind == ieForce)
          check graph.mdForces.len == nforces
          check graph.trace[^1].gauge.nodeKey == graph.finalState.gauge.nodeKey
          check graph.trace[^1].momentum.nodeKey == graph.finalState.momentum.nodeKey
          var theta = data.theta0.mapIt(it)
          var mom = data.momentum0.mapIt(it)
          var fcur, shifted: seq[float]
          var cur = graph.initialState.gauge      # gauge node of the state
          var curp = graph.initialState.momentum  # momentum node of the state
          var shiftGauge: Ggauge                  # gauge displaced along a force
          var driftTotal, kickTotal = 0.0
          var nf = 0
          for j,event in graph.trace:
            let record = data.events[j]
            let want = expected[j]
            check event.kind == want.kind and event.step == want.step
            if minimalNorm and event.kind != ieForce:
              check record.coefficient == want.coefficient
            let displaced = shiftGauge != nil and event.gauge.nodeKey == shiftGauge.nodeKey
            let at = if displaced: shifted else: theta
            case event.kind
            of ieForce:
              check event.coefficient == nil
              check displaced or event.gauge.nodeKey == cur.nodeKey
              check event.momentum.nodeKey == curp.nodeKey
              check event.force.nodeKey == graph.mdForces[nf].nodeKey
              inc nf
              fcur = at.mapIt(kappa*sin(it))
              checkArray(record.force, fcur, "force")
            of ieShift:
              check event.momentum.nodeKey == curp.nodeKey
              check event.force.nodeKey == graph.trace[j-1].force.nodeKey
              shifted = newSeq[float](theta.len)
              for i in 0..<theta.len: shifted[i] = theta[i]+record.coefficient*fcur[i]
              shiftGauge = event.gauge
            of ieKick:
              check event.gauge.nodeKey == cur.nodeKey
              check event.force.nodeKey == graph.trace[j-1].force.nodeKey
              for i in 0..<theta.len: mom[i] += record.coefficient*fcur[i]
              kickTotal += record.coefficient
              curp = event.momentum
            of ieDrift:
              check event.force == nil
              check event.momentum.nodeKey == curp.nodeKey
              for i in 0..<theta.len: theta[i] += record.coefficient*mom[i]
              driftTotal += record.coefficient
              cur = event.gauge
            # Shift events and the force events they precede see the displaced gauge.
            let angles = if event.kind == ieShift or displaced: shifted else: theta
            checkArray(record.gaugeRe,angles.mapIt(cos(it)),"link real part")
            checkArray(record.gaugeIm,angles.mapIt(sin(it)),"link imaginary part")
            checkArray(record.momentum,mom,"momentum")
            checkArray(record.theta,theta,"unwrapped angle diagnostic")
            checkArray(record.momentumRe,newSeq[float](mom.len),"momentum real part")
          check nf == nforces
          check cur.nodeKey == graph.finalState.gauge.nodeKey
          check curp.nodeKey == graph.finalState.momentum.nodeKey
          check abs(driftTotal-float(n)*dt) < 1e-12
          check abs(kickTotal+float(n)*dt) < 1e-12
          check abs(data.initialHamiltonian-energy(data.theta0,data.momentum0,kappa)) < 1e-10
          check abs(data.finalHamiltonian-energy(theta,mom,kappa)) < 1e-10
          if n == 10:
            check theta.max > PI
qexFinalize()
