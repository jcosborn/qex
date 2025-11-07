import base, maths
import gaugeUtils, gaugeAction

#[

d/dt Vt = Z(Vt) Vt

Runge-Kutta:

W0 <- Vt
W1 <- exp(1/4 Z0) W0
W2 <- exp(8/9 Z1 - 17/36 Z0) W1
V(t+eps) <- exp(3/4 Z2 - 8/9 Z1 + 17/36 Z0) W2

where

Zi = eps Z(Wi)

]#

template gaugeFlow*(g: array|seq, steps: int, eps: float, measure: untyped) =
  ## Wilson flow.
  ## The input gauge field will be modified.
  ## `wflowT` is injected for `measure`.
  proc flowProc {.gensym.} =
    tic("flowProc")
    const nc = g[0][0].nrows.float
    var
      p = g[0].l.newGauge  # mom
      f = g[0].l.newGauge  # force
      n = 1
    while true:
      let t = n * eps
      let epsnc = eps * nc  # compensate force normalization
      f.gaugeForce g
      threads:
        for mu in 0..<f.len:
          for e in g[mu]:
            var v {.noinit.}:type(load1(f[0][0]))
            v := (-1.0/4.0)*epsnc*f[mu][e]
            let t = exp(v)*g[mu][e]
            p[mu][e] := v
            g[mu][e] := t
      f.gaugeForce g
      threads:
        for mu in 0..<f.len:
          for e in g[mu]:
            var v {.noinit.}:type(load1(f[0][0]))
            v := (-8.0/9.0)*epsnc*f[mu][e] + (-17.0/9.0)*p[mu][e]
            let t = exp(v)*g[mu][e]
            p[mu][e] := v
            g[mu][e] := t
      f.gaugeForce g
      threads:
        for mu in 0..<f.len:
          for e in g[mu]:
            var v {.noinit.}:type(load1(f[0][0]))
            v := (-3.0/4.0)*epsnc*f[mu][e] - p[mu][e]
            let t = exp(v)*g[mu][e]
            g[mu][e] := t
      let wflowT {.inject,used.} = t
      measure
      inc n
      if steps>0 and n>steps:
        break
    toc("end")
  flowProc()
template gaugeFlow*(g: array|seq, eps: float, measure: untyped) =
  gaugeFlow(g, 0, eps):
    measure

import ../algorithms/rk

type
  GaugeFlowRKOp*[G] = object of RK2NOp[G,G]
    f: G
  GaugeFlowRKAdaptiveOp*[G] = object of RK2NAdaptiveOp[G,G]
    f: G

proc mkGaugeFlowRKOp*[G](g: G): GaugeFlowRKOp[G] =
  ## the returned op will hold a reference to the input g
  var op = GaugeFlowRKOp[G](
    y: g,
    delta: g[0].l.newGauge,
    f: g[0].l.newGauge)
  proc adv(base: var type(g), d: var type(g), a: float, b: float) =
    const nc = op.y[0][0].nrows.float
    op.f.gaugeForce(base)
    let u = cast[ptr cArray[type(base[0])]](unsafeAddr(base[0]))
    let f = cast[ptr cArray[type(d[0])]](unsafeAddr(d[0]))
    let ff = cast[ptr cArray[type(d[0])]](unsafeAddr(op.f[0]))
    let n = base.len
    if a == 0.0:
      threads:
        for mu in 0..<n:
          for e in u[mu]:
            var v {.noinit.}:type(load1(u[0][0]))
            v := (-b * nc) * ff[mu][e]
            let t = exp(v) * u[mu][e]
            f[mu][e] := v
            u[mu][e] := t
    else:
      threads:
        for mu in 0..<n:
          for e in u[mu]:
            var v {.noinit.}:type(load1(u[0][0]))
            v := a * f[mu][e] + (-b * nc) * ff[mu][e]
            let t = exp(v) * u[mu][e]
            f[mu][e] := v
            u[mu][e] := t
  op.advance = adv
  op

proc mkGaugeFlowRKAdaptiveOp*[G](g: G): GaugeFlowRKAdaptiveOp[G] =
  ## the returned op will hold a reference to the input g
  mixin simdMax
  var rk = g.mkGaugeFlowRKOp
  var op = GaugeFlowRKAdaptiveOp[G](
    y: rk.y,
    delta: rk.delta,
    f: rk.f,
    advance: rk.advance,
    y0: g[0].l.newGauge,
  )
  proc assignG(dst: var type(g), src: type(g)) =
    let u = cast[ptr cArray[type(dst[0])]](unsafeAddr(dst[0]))
    let s = cast[ptr cArray[type(src[0])]](unsafeAddr(src[0]))
    let n = dst.len
    threads:
      for mu in 0..<n:
        for e in u[mu]:
          u[mu][e] := s[mu][e]
  op.assign = assignG
  # op.errDelta = proc(d: type(g)): float =
  #   var r = 0.0
  #   threads:
  #     var p2t = 0.0
  #     for i in 0..<d.len:
  #       p2t += d[i].norm2
  #     threadMaster: r = p2t
  #   const nc = d[0][0].nrows
  #   sqrt(r / float((nc*nc-1)*d.len*d[0].l.physVol))
  op.errDelta = proc(d: type(g)): float =
    var res = 0.0
    threads:
      var r = 0.0
      for i in 0..<d.len:
        for s in d[i]:
          let e = d[i][s].norm2.simdMax
          if r<e: r = e
      threadRankMax r
      threadMaster: res = r
    const nc = d[0][0].nrows
    sqrt(res/float(nc*nc-1))
  op

proc gaugeFlowRK*[G](g: var G, steps: int, eps: float, coeffs: auto, measure: proc(wflowT: float) {.closure.} = nil) =
  ## Wilson flow using encapsulated 2N RK operator with arbitrary RK2N coefficients.
  tic("flowProcRK")
  var op = g.mkGaugeFlowRKOp
  rk2n(op, coeffs, steps, eps, measureCb=measure)
  toc("endRK")

proc gaugeFlowRKAdaptive*[G](g: var G, t: float, h0: float, coeffs: auto, tol: float, safety: float = 0.95, maxSteps:int = 100000, controllerExp: float = 1.0/3.0, measure: proc(wflowT: float) {.closure.} = nil): RKAdaptiveStats =
  ## Adaptive Wilson flow to total time t using encapsulated 2N RK operator and arbitrary coefficients.
  ## See comments under rk2nAdaptive for more details.
  tic("flowProcAdaptive")
  var op = g.mkGaugeFlowRKAdaptiveOp
  result = rk2nAdaptive(op, coeffs, t, h0, deltaTol=tol, safety=safety, maxSteps=maxSteps, controllerExp=controllerExp, measureCb=measure)
  toc("endAdaptive")

when isMainModule:
  import qex, gauge, physics/qcdTypes
  import os, sequtils

  proc printPlaq(g: auto) =
    let
      p = g.plaq
      sp = 2.0*(p[0]+p[1]+p[2])
      tp = 2.0*(p[3]+p[4]+p[5])
    echo "plaq ",p
    echo "plaq ss: ",sp," st: ",tp," tot: ",p.sum

  qexInit()

  let
    fn = if paramCount() > 0: paramStr 1 else: ""
    lat = if fn.len == 0: @[8,8,8,8] else: fn.getFileLattice
    lo = lat.newLayout
  var g = lo.newGauge
  if fn.len == 0:
    g.random
  elif 0 != g.loadGauge fn:
    echo "ERROR: couldn't load gauge file: ",fn
    qexFinalize()
    quit(-1)
  g.printPlaq

  # Compare plaquettes
  proc relDiff(a, b: auto): float =
    let pa = a.plaq
    let pb = b.plaq
    zip(pa, pb).foldl(a + abs(b[0] - b[1]), 0.0) / max(1e-30, pa.sum)

  var gRef = g.newGauge
  gRef.gaugeFlowRK(120, 0.0005, RK64_2N, (proc(t:float) =
    var nt {.global.} = 0.01
    if t < (1.0-1.0e-12)*nt:
      return
    nt += 0.01
    echo "WFLOW RK64 ref ", t
    gRef.printPlaq))

  template runFlow(coeff:auto) =
    echo "gauge flow ",astToStr(coeff)
    var d = newseq[float](0)
    for (steps, eps) in [(6,0.01), (12,0.005)]:
      var gRK = g.newGauge
      gRK.gaugeFlowRK(steps, eps, coeff)
      d.add relDiff(gRef, gRK)
      echo "    steps = ",steps," eps = ",eps,"  rel plaq diff = ",d[^1]
    echo "  error scaling: ",(ln(d[0]/d[1])/ln(2.0)),"  per stage coeff: ",d[^1]/pow(0.005/coeff.beta.len,ln(d[0]/d[1])/ln(2.0))
  runFlow(RK3W6_2N)
  runFlow(RK3W7_2N)
  runFlow(RK43_1_2N)
  runFlow(RK43_2_2N)
  runFlow(RK43_3_2N)
  runFlow(RK43_4_2N)
  runFlow(RK53_1_2N)
  runFlow(RK53_2_2N)
  runFlow(RK53_3_2N)
  runFlow(RK53_4_2N)
  runFlow(RK4CK_2N)
  runFlow(RK4BBB_2N)
  runFlow(RK64_2N)

  template runFlowAd(tol) =
    block:
      var gRK = g.newGauge
      let res = gRK.gaugeFlowRKAdaptive(6*0.01, 0.01, RK53_4_2N, tol)
      echo "gauge flow adaptive RK53_4 rel plaq diff = ",relDiff(gRef, gRK), "  steps=", res.steps, " acc=", res.accepts, " rej=", res.rejects, "  h∈[", res.minH, ", ", res.maxH, "]"
  runFlowAd(1e-5)
  runFlowAd(1e-7)
  runFlowAd(1e-9)
  runFlowAd(1e-11)

  # Forward/backward gauge flow tests (fixed-step)
  template runFlowRev(coeff:auto) =
    echo "Fwd/Back Fixed ",astToStr(coeff)
    var d = newseq[float](0)
    for (steps, eps) in [(10,0.01), (20,0.005)]:
      var gr = g.newGauge
      var op = gr.mkGaugeFlowRKOp
      op.rk2n(coeff, steps, eps)
      op.rk2n(coeff, steps, -eps)
      d.add relDiff(g, gr)
      echo "    steps = ",steps," eps = ",eps,"  rel plaq diff = ",d[^1]
    echo "  error scaling: ",(ln(d[0]/d[1])/ln(2.0)),"  per stage coeff: ",d[^1]/pow(0.005/coeff.beta.len,ln(d[0]/d[1])/ln(2.0))
  runFlowRev(RK3W6_2N)
  runFlowRev(RK3W7_2N)
  runFlowRev(RK43_1_2N)
  runFlowRev(RK43_2_2N)
  runFlowRev(RK43_3_2N)
  runFlowRev(RK43_4_2N)
  runFlowRev(RK53_1_2N)
  runFlowRev(RK53_2_2N)
  runFlowRev(RK53_3_2N)
  runFlowRev(RK53_4_2N)
  runFlowRev(RK4CK_2N)
  runFlowRev(RK4BBB_2N)
  runFlowRev(RK64_2N)

  template runFlowAdRev(tol) =
    block:
      var gRK = g.newGauge
      var op = gRK.mkGaugeFlowRKAdaptiveOp
      let resF = rk2nAdaptive(op, RK53_4_2N, 0.1, 0.01, tol)
      let resB = rk2nAdaptive(op, RK53_4_2N, -0.1, 0.01, tol)
      echo "Fwd/Back adaptive RK53_4 rel plaq diff = ",relDiff(g, gRK)
      echo "    Fwd steps=", resF.steps, " acc=", resF.accepts, " rej=", resF.rejects, "  h∈[", resF.minH, ", ", resF.maxH, "]"
      echo "    Bwd steps=", resB.steps, " acc=", resB.accepts, " rej=", resB.rejects, "  h∈[", resB.minH, ", ", resB.maxH, "]"
  runFlowAdRev(1e-5)
  runFlowAdRev(1e-7)
  runFlowAdRev(1e-9)
  runFlowAdRev(1e-11)

  # Default flow (original 3-stage wrapper)
  g.gaugeFlow(6, 0.01):
    echo "WFLOW default ", wflowT
    g.printPlaq

  when g[0][0].nrows == 1 or g[0][0].nrows == 3:
    if fn.len == 0:
      when g[0][0].nrows == 1:
        when defined(FUELCompat):
          const
            p0 = @[
              0.01554978485185948,
              0.01864915539989097,
              0.02185632226074754,
              0.01631465138128996,
              0.01599186401283072,
              0.02143528647233249]
        else:
          const
            p0 = @[
              0.02071819961326493,
              0.01982444094659434,
              0.01998727952235725,
              0.01837579824880835,
              0.01528029895697336,
              0.01782410235480616]
      when g[0][0].nrows == 3:
        when defined(FUELCompat):
          const
            p0 = @[
              0.01984173412117658,
              0.01775826434617258,
              0.01814446406048498,
              0.01960829801483532,
              0.01913142046361973,
              0.01953858380137068]
        else:
          const
            p0 = @[
              0.01960725848281519,
              0.01982378149813489,
              0.01938877647467847,
              0.0185899778070918,
              0.0180821938831715,
              0.01876842496122964]
      let p = g.plaq
      let d = zip(p,p0).foldl(a + abs(b[0]-b[1]), 0.0) / sum(p0)
      if d > 2e-14:
        echo "Test failed, relative diff: ", d
        echo "Expected:\t",p0
        echo "Actual:\t",p
        qexExit -1

  qexFinalize()
