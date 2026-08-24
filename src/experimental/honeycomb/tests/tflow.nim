## Task W test, part 1: `hcflow.nim` (gradient flow + its normalisation).
##
##   export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
##   cd build_mac && make src/experimental/honeycomb/tests/tflow.nim && \
##     OMP_NUM_THREADS=4 ./bin/tflow
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
## 3. Verification runs with the calibrated hcFlowCflow, plus an integrator
##    step-size independence check.
## 4. On a warm configuration: E (hexagon clover) decreases monotonically and
##    SU(3) is preserved along a long flow.

import std/[math, strformat, unittest]
import qex except epsilon
import physics/qcdTypes
import gauge, gauge/wflow
import ../hcgeom
import ../hclayout
import ../hcgauge
import ../hcaction
import ../hcflow
import ../hctopo

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

template ok(msg: string, cond: bool) =
  let c = cond
  if c: echo "PASS: ", msg
  else: echo "FAIL: ", msg
  check c

qexInit()

const nc = getDefaultNc()

# ---------------------------------------------------------------------------
# weak Abelian plane wave from exact line integrals (recipe of tests/taction.nim)
# ---------------------------------------------------------------------------

proc phaseIntegral(a, b: float): float =
  ## int_0^1 cos(a + t b) dt
  if abs(b) < 1e-12: cos(a)
  else: (sin(a+b) - sin(a))/b

proc linkPhase(x, n, eps, pv: array[4, float]): float =
  ## int_x^{x+n} A.dl  for  A_mu(x) = eps_mu cos(p.x)
  var a = 0.0
  var b = 0.0
  var en = 0.0
  for mu in 0..<4:
    a += pv[mu]*x[mu]
    b += pv[mu]*n[mu]
    en += eps[mu]*n[mu]
  en*phaseIntegral(a, b)

template setPhase(m: untyped, th: float) =
  ## m := exp(i th T), T = diag(1,-1,0); assumes m starts as the identity
  m[0, 0].re := cos(th)
  m[0, 0].im := sin(th)
  m[1, 1].re := cos(th)
  m[1, 1].im := -sin(th)

proc setAbelianHc(hg: auto, eps, pv: array[4, float]) =
  ## honeycomb links U = exp(i phi T) from exact line integrals of A
  hg.unit
  let lo = hg.lo
  for i in lo.sites:
    var y, yb: array[4, float]
    for mu in 0..<4:
      y[mu] = lo.coords[mu][i].float    # A(y) at integer coords
      yb[mu] = y[mu] + 0.5              # B(y) at half-integer coords
    for mu in 0..<4:
      var n: array[4, float]
      n[mu] = 1.0
      setPhase(hg.uA[mu]{i}, linkPhase(y, n, eps, pv))
      setPhase(hg.uB[mu]{i}, linkPhase(yb, n, eps, pv))
    for d in 0..<nDiag:
      var n: array[4, float]
      for mu in 0..<4:
        n[mu] = float((d shr mu) and 1) - 0.5
      setPhase(hg.uD[d]{i}, linkPhase(yb, n, eps, pv))

proc setAbelianCubic(g: auto, lo: auto, eps, pv: array[4, float]) =
  ## cubic links from the same continuum field
  g.unit
  for i in lo.sites:
    var y: array[4, float]
    for mu in 0..<4:
      y[mu] = lo.coords[mu][i].float
    for mu in 0..<4:
      var n: array[4, float]
      n[mu] = 1.0
      setPhase(g[mu]{i}, linkPhase(y, n, eps, pv))

proc wilsonAction(gact: GaugeActionCoeffs, g: auto): float =
  ## S = beta sum_{x,mu<nu} (1 - Re Tr P/N); gaugeAction1 omits the constant
  gact.gaugeAction1(g) + gact.plaq*6.0*float(g[0].l.physVol)

# ---------------------------------------------------------------------------
# decay-rate extraction
# ---------------------------------------------------------------------------

proc plateauRate(ss: seq[float], dt: float, frac = 0.25):
    tuple[rate, drift: float] =
  ## rate_i = ln(S_{i-1}/S_i)/(2 dt); average over the last `frac` of the
  ## series, with the max deviation inside that window as `drift`.
  var r: seq[float]
  for i in 1..<ss.len:
    r.add ln(ss[i-1]/ss[i])/(2.0*dt)
  let m = max(2, int(frac*r.len.float))
  var s = 0.0
  for i in r.len-m..<r.len: s += r[i]
  let avg = s/m.float
  var d = 0.0
  for i in r.len-m..<r.len: d = max(d, abs(r[i]-avg))
  (avg, d)

proc fit3(x, y: array[3, float]): array[3, float] =
  ## exact solve of y = c0 + c1 x + c2 x^2 through 3 points
  let
    d = (x[1]-x[0])*(x[2]-x[0])*(x[2]-x[1])
    c2 = ((y[2]-y[0])*(x[1]-x[0]) - (y[1]-y[0])*(x[2]-x[0]))/d
    c1 = (y[1]-y[0])/(x[1]-x[0]) - c2*(x[1]+x[0])
    c0 = y[0] - c1*x[0] - c2*x[0]*x[0]
  [c0, c1, c2]

const
  epsAmp = 3.0e-3
    ## weak-field amplitude.  Large enough that S(t) in the late measurement
    ## windows stays >~1e5 above the cancellation-noise floor of the action
    ## sum (at 1e-3 the extracted rate moved by ~5e-5 rel between thread
    ## counts through that floor); small enough that nonlinear corrections
    ## O(epsAmp^2) ~ 1e-5 are far below the 1e-3 rational-pinning gate.
  ehat = [0.0, 1.0, -0.6, 0.3]        # transverse polarisation, p along dir 0
  flowEps = 0.05

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

proc hcRate(ns, k: int, cflow: float, eps = flowEps, tmax = 8.0):
    tuple[rate, drift, p2: float] =
  ## same measurement on the honeycomb with hcGaugeFlow at given cflow.
  ## The rate window sits late (last quarter of [0,tmax]) so that the fast
  ## "optical" transients (relative axis/diagonal phase mismatches, rate O(1))
  ## have died out; `drift` monitors the residual contamination.
  var pv, epsv: array[4, float]
  pv[0] = 2.0*PI*float(k)/float(ns)
  for mu in 0..<4: epsv[mu] = epsAmp*ehat[mu]
  let hl = newHcLayout([ns, ns, ns, ns])
  var g = newHcGauge(hl)
  setAbelianHc(g, epsv, pv)
  var w = newHcActionWork(g)
  var ss = @[hcAction(w, 1.0, g)]
  let nsteps = int(tmax/eps + 0.5)
  g.hcGaugeFlow(nsteps, eps, cflow):
    ss.add hcAction(w, 1.0, g)
    discard wflowT
  let (r, d) = plateauRate(ss, eps)
  (r, d, pv[0]*pv[0])

# ---------------------------------------------------------------------------

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
    ok(&"shipped hcFlowCflow = {hcFlowCflow} agrees with the measurement " &
       &"(rel {abs(cflow0-hcFlowCflow)/hcFlowCflow:.1e})",
       abs(cflow0 - hcFlowCflow) < 1e-2*hcFlowCflow)

  test "3. calibrated flow is the heat kernel; step-size independence":
    echo &"  cflow = hcFlowCflow = {hcFlowCflow}:"
    echo "  Ns k   p^2       rate         rate/p^2    plateau-drift"
    # the calibrated flow decays 6x faster, so the measurement window must
    # end earlier or S(t) reaches the cancellation-noise floor of the action
    var xs2, ys2: array[2, float]
    let cases = [(12, 1, 3.0), (8, 1, 2.0)]
    for i in 0..<2:
      let (ns, k, tm) = cases[i]
      let (r, d, p2) = hcRate(ns, k, hcFlowCflow, flowEps, tm)
      xs2[i] = p2
      ys2[i] = r/p2
      echo &"  {ns:2d} {k:2d}  {p2:.5f}  {r:.8f}   {r/p2:.7f}   {d:.2e}"
      echo &"        O(a^2 p^2) artifact at this p: {ys2[i]-1.0: .5f}"
    let cver = ys2[0] - xs2[0]*(ys2[1]-ys2[0])/(xs2[1]-xs2[0])
    echo &"  rate/p^2 extrapolated to p^2 -> 0: {cver:.6f}"
    ok(&"calibrated flow: |rate/p^2 - 1| = {abs(cver-1.0):.2e} < 1e-2 " &
       "after artifact extrapolation", abs(cver - 1.0) < 1e-2)
    # integrator step-size independence (8^4, k=1)
    let (rA, _, p2A) = hcRate(8, 1, hcFlowCflow, flowEps, 2.0)
    let (rB, _, _) = hcRate(8, 1, hcFlowCflow, 0.5*flowEps, 2.0)
    let dEps = abs(rB/rA - 1.0)
    echo &"  8^4 k=1: rate(eps={flowEps}) = {rA:.9f},  ",
         &"rate(eps={0.5*flowEps}) = {rB:.9f},  rel diff {dEps:.2e}"
    ok(&"RK3 step-size independence ({dEps:.2e} < 1e-3)", dEps < 1e-3)
    discard p2A

  test "4. warm config: E monotone along flow, SU(3) preserved":
    # NB: [6,6,6,6] cells are rejected by the V=4 SIMD layout ("can't lay out
    # inner geom", cf. task C's note); 8^4 works.
    let hl = newHcLayout([8, 8, 8, 8])
    var r = hl.lo.newRNGField(RngMilc6, 246813579'u64)
    var g = newHcGauge(hl)
    threads:
      g.warm(0.35, r)
    var wt = newHcTopoWork(g)
    var es: seq[float]
    var sus: seq[float]
    block:
      let (e0, q0) = hcEQ(wt, g)
      es.add e0
      discard q0
    let d0 = g.checkSU
    sus.add d0.max
    g.hcGaugeFlow(60, flowEps, hcFlowCflow):
      let (e, q) = hcEQ(wt, g)
      es.add e
      discard q
      if (int(wflowT/flowEps + 0.5) mod 10) == 0:
        let d = g.checkSU
        sus.add d.max
    var mono = true
    for i in 1..<es.len:
      if es[i] >= es[i-1]:
        mono = false
        echo &"  NOT monotone at step {i}: E {es[i-1]:.6e} -> {es[i]:.6e}"
    echo &"  E(t=0) = {es[0]:.6f}  E(t=3) = {es[^1]:.6e}  ({es.len-1} steps)"
    ok("E decreases monotonically along the flow (60 steps to t = 3)", mono)
    var sumax = 0.0
    for s in sus: sumax = max(sumax, s)
    echo &"  checkSU max over the flow: {sumax:.2e}  (no reunitarisation done)"
    ok(&"flow preserves SU(3) to 1e-10 without reunit ({sumax:.2e})",
       sumax < 1e-10)

qexFinalize()
