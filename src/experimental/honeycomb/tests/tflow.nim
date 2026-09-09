## `hcflow.nim` (gradient flow + its normalisation).
##
## The defining test (FORMULATION.md 4.4): for a weak Abelian plane wave
## A_mu(x) = eps_mu cos(p.x), embedded via T = diag(1,-1,0) with link phases
## from EXACT straight-line integrals, the flow must act as the heat kernel:
## the mode amplitude decays as exp(-t p^2), i.e. the (gauge-invariant,
## quadratic) action decays as exp(-2 lambda t) with lambda = p^2 (1+O(a^2p^2)).
##
## 1. The identical measurement on a cubic lattice with QEX's own `gaugeFlow`
##    validates the harness: there the exact linearised rate is
##    lambda = phat^2 = 4 sin^2(p/2), so rate/phat^2 = 1 to integrator
##    precision and rate/p^2 = 1 + O(p^2) (reported).
## 2. The honeycomb rate with a provisional cflow = 1 gives c_HC = lambda/p^2
##    extrapolated to p^2 -> 0 (three momenta, O(p^2) artifact removed);
##    cflow = 1/c_HC.  Measured: c_HC = 1/6, cflow = 6 (pinned exact; the
##    Rayleigh-quotient estimate 2/9 of hcflow.nim's docs is only an upper
##    bound and is measurably wrong -- see there).
## 3. Verification runs with the calibrated cflow, plus an integrator
##    step-size independence check.
## 4. On a warm configuration: E (hexagon clover) decreases monotonically and
##    SU(3) is preserved along a long flow.

import std/[math, strformat, unittest]
import qex except epsilon
import physics/qcdTypes
import gauge, gauge/wflow
import ../hcgeom
import ../hcgauge
import ../hcaction
import ../hcflow
import ../hctopo
import helpers

qexInit()

const flowEps = 0.05

proc cubicRate(ns, k: int): tuple[rate, drift, p2: float] =
  ## heat-kernel decay rate of the plane-wave mode under QEX's gaugeFlow
  var pv, eps: array[4, float]
  pv[0] = 2.0*PI*float(k)/float(ns)
  for mu in 0..<4: eps[mu] = epsAmp*ehat[mu]
  let lo = newLayout(@[ns, ns, ns, ns])
  var g = lo.newGauge
  setAbelianCubic(g, lo, eps, pv)
  let gact = GaugeActionCoeffs(plaq: 1.0)
  var ss = @[wilsonAction(gact, g)]
  let nsteps = int(2.5/flowEps + 0.5)
  g.gaugeFlow(nsteps, flowEps):
    ss.add wilsonAction(gact, g)
    discard wflowT
  let (r, d) = plateauRate(ss, flowEps)
  (r, d, pv[0]*pv[0])

proc hcRate(ns, k: int, cflow: float, eps = flowEps, tmax = 4.0):
    tuple[rate, drift, p2: float] =
  ## same measurement on the honeycomb with flow at given cflow.
  ## The rate window sits late (last quarter of [0,tmax]) so that the fast
  ## "optical" transients (relative axis/diagonal phase mismatches, rate O(1)
  ## at cflow = 1) have died out; `drift` monitors the residual contamination.
  var pv, epsv: array[4, float]
  pv[0] = 2.0*PI*float(k)/float(ns)
  for mu in 0..<4: epsv[mu] = epsAmp*ehat[mu]
  let hl = newLayout([ns, ns, ns, ns])
  var g = newHcGauge(hl)
  setAbelianHc(g, epsv, pv)
  var w = newActionWork(g)
  var ss = @[action(w, 1.0, g)]
  let nsteps = int(tmax/eps + 0.5)
  g.flow(nsteps, eps, cflow):
    ss.add action(w, 1.0, g)
    discard wflowT
  let (r, d) = plateauRate(ss, eps)
  (r, d, pv[0]*pv[0])

suite "hcflow":

  test "1. cubic harness: QEX gaugeFlow is the heat kernel (rate = phat^2)":
    # This pins the measurement harness AND the flow-time convention `t` to
    # QEX's correctly normalised cubic flow before touching the honeycomb.
    echo "  Ns k   p^2       rate        rate/p^2   rate/phat^2  plateau-drift"
    var xs, ys: array[3, float]
    let cases = [(12, 1), (8, 1), (12, 2)]
    var allSharp = true
    for i in 0..<3:
      let (ns, k) = cases[i]
      let (r, d, p2) = cubicRate(ns, k)
      let p = sqrt(p2)
      let ph2 = 4.0*sin(0.5*p)^2
      xs[i] = p2
      ys[i] = r/p2
      echo &"  {ns:2d} {k:2d}  {p2:.5f}  {r:.8f}  {r/p2:.6f}  {r/ph2:.8f}  {d:.2e}"
      if abs(r/ph2 - 1.0) > 1e-3: allSharp = false
    ok("rate/phat^2 = 1 to 1e-3 at all momenta (exact lattice prediction)",
       allSharp)
    let c = fit3(xs, ys)
    echo &"  rate/p^2 extrapolated to p^2 -> 0: {c[0]:.6f}",
         &"   (O(p^2) coeff {c[1]:.4f}, cf. exact -1/12 = {-1.0/12.0:.4f})"
    ok(&"cubic rate/p^2 -> 1 as p -> 0 (extrap {c[0]:.6f})",
       abs(c[0] - 1.0) < 3e-3)
    ok(&"cubic O(p^2) artifact coeff = -1/12 ({c[1]:.4f})",
       abs(c[1] + 1.0/12.0) < 0.02)

  test "2. honeycomb flow calibration: measure cflow from the heat kernel":
    echo "  provisional cflow = 1:"
    echo "  Ns k   p^2       rate         rate/p^2    plateau-drift"
    var xs, ys: array[3, float]
    let cases = [(12, 1), (8, 1), (12, 2)]
    for i in 0..<3:
      let (ns, k) = cases[i]
      let (r, d, p2) = hcRate(ns, k, 1.0)
      xs[i] = p2
      ys[i] = r/p2
      echo &"  {ns:2d} {k:2d}  {p2:.5f}  {r:.8f}   {r/p2:.7f}   {d:.2e}"
    let cq = fit3(xs, ys)
    let clin = ys[0] - xs[0]*(ys[1]-ys[0])/(xs[1]-xs[0])
    let cHC = cq[0]
    let cErr = abs(cq[0] - clin)
    echo &"  c_HC = lambda/p^2 (p->0):  quad-3pt {cq[0]:.7f}   lin-2pt {clin:.7f}"
    echo &"  c_HC = {cHC:.7f} +- {cErr:.1e}   (O(p^2) coeff {cq[1]:.4f})"
    let cflow0 = 1.0/cHC
    echo &"  ==> cflow = 1/c_HC = {cflow0:.7f} +- {cErr/(cHC*cHC):.1e}"
    # simple-rational scan: smallest denominator within 1e-3 wins
    var bestN, bestD = 0
    var bestErr = 1.0
    for den in 1..12:
      let num = int(round(cflow0*den.float))
      if num < 1: continue
      let e = abs(cflow0 - num.float/den.float)
      if e < 1e-3:
        bestErr = e
        bestN = num
        bestD = den
        break
      if e < bestErr:
        bestErr = e
        bestN = num
        bestD = den
    if bestErr < 1e-3:
      echo &"  *** cflow is the SIMPLE RATIONAL {bestN}/{bestD} = ",
           &"{bestN.float/bestD.float:.6f} (|dev| = {bestErr:.1e} < 1e-3):",
           " pinning cflow = 6 exactly ***"
    else:
      echo &"  cflow is NOT within 1e-3 of a simple rational (closest ",
           &"{bestN}/{bestD}, dev {bestErr:.1e})"
    ok(&"cflow rational pin: 6 (dev {abs(cflow0-6.0):.1e})",
       bestErr < 1e-3 and bestN == 6 and bestD == 1)
    ok(&"shipped cflow = {cflow} agrees with the measurement " &
       &"(rel {abs(cflow0-cflow)/cflow:.1e})",
       abs(cflow0 - cflow) < 1e-2*cflow)

  test "3. calibrated flow is the heat kernel; step-size independence":
    echo &"  cflow = {cflow}:"
    echo "  Ns k   p^2       rate         rate/p^2    plateau-drift"
    # the calibrated flow decays 6x faster, so the measurement window must
    # end earlier or S(t) reaches the cancellation-noise floor of the action
    var xs2, ys2: array[2, float]
    let cases = [(12, 1, 3.0), (8, 1, 2.0)]
    for i in 0..<2:
      let (ns, k, tm) = cases[i]
      let (r, d, p2) = hcRate(ns, k, cflow, flowEps, tm)
      xs2[i] = p2
      ys2[i] = r/p2
      echo &"  {ns:2d} {k:2d}  {p2:.5f}  {r:.8f}   {r/p2:.7f}   {d:.2e}"
      echo &"        O(a^2 p^2) artifact at this p: {ys2[i]-1.0: .5f}"
    let cver = ys2[0] - xs2[0]*(ys2[1]-ys2[0])/(xs2[1]-xs2[0])
    echo &"  rate/p^2 extrapolated to p^2 -> 0: {cver:.6f}"
    ok(&"calibrated flow: |rate/p^2 - 1| = {abs(cver-1.0):.2e} < 1e-2 " &
       "after artifact extrapolation", abs(cver - 1.0) < 1e-2)
    # integrator step-size independence (8^4, k=1)
    let (rA, _, _) = hcRate(8, 1, cflow, flowEps, 2.0)
    let (rB, _, _) = hcRate(8, 1, cflow, 0.5*flowEps, 2.0)
    let dEps = abs(rB/rA - 1.0)
    echo &"  8^4 k=1: rate(eps={flowEps}) = {rA:.9f},  ",
         &"rate(eps={0.5*flowEps}) = {rB:.9f},  rel diff {dEps:.2e}"
    ok(&"RK3 step-size independence ({dEps:.2e} < 1e-3)", dEps < 1e-3)

  test "4. warm config: E monotone along flow, SU(3) preserved":
    # NB: [6,6,6,6] cells are rejected by the V=4 SIMD layout ("can't lay out
    # inner geom", cf. task C's note); 8^4 works.
    let hl = newLayout([8, 8, 8, 8])
    var r = hl.newRNGField(RngMilc6, 246813579'u64)
    var g = newHcGauge(hl)
    threads:
      g.warm(0.35, r)
    var wt = newTopoWork(g)
    var es = @[EQ(wt, g).e]
    var sus = @[g.checkSU.max]
    g.flow(60, flowEps, cflow):
      es.add EQ(wt, g).e
      if (int(wflowT/flowEps + 0.5) mod 10) == 0:
        sus.add g.checkSU.max
    var mono = true
    for i in 1..<es.len:
      if es[i] >= es[i-1]:
        mono = false
        echo &"  NOT monotone at step {i}: E {es[i-1]:.6e} -> {es[i]:.6e}"
    echo &"  E(t=0) = {es[0]:.6f}  E(t=3) = {es[^1]:.6e}  ({es.len-1} steps)"
    ok("E decreases monotonically along the flow (60 steps to t = 3)", mono)
    let sumax = max(sus)
    echo &"  checkSU max over the flow: {sumax:.2e}  (no reunitarisation done)"
    ok(&"flow preserves SU(3) to 1e-10 without reunit ({sumax:.2e})",
       sumax < 1e-10)

qexFinalize()
