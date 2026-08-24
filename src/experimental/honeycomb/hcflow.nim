## Gradient flow for the triangle action on the 16-cell honeycomb (task **W1**).
##
## Flow equation (Luscher, JHEP 1008:071):  dV_l/dt = Z_l(V) V_l, integrated
## with the same 3-stage Runge-Kutta as `src/gauge/wflow.nim:4-19`:
##
##   W0 <- Vt
##   W1 <- exp(1/4 Z0) W0
##   W2 <- exp(8/9 Z1 - 17/36 Z0) W1
##   V(t+eps) <- exp(3/4 Z2 - 8/9 Z1 + 17/36 Z0) W2,   Zi = eps Z(Wi)
##
## `wflow.nim` is hardwired to the plaquette force, so this is a fork driving
## the RK with task A's `hcForce` (the triangle action).
##
## Flow-time normalisation and `cflow`
## -----------------------------------
## QEX's cubic `gaugeFlow` uses Z = -eps * nc * gaugeForce(plaq=1); with that
## normalisation a weak Abelian plane wave decays exactly as the heat kernel,
## A(p,t) = exp(-t p^2) A(p,0) (up to O(a^2 p^2) artifacts), i.e. `t` is the
## standard continuum flow time in lattice units a^2.  (Verified in
## tests/tflow.nim test 1, and independently by task C's t0 values.)
##
## Here we define, in exact structural parallel,
##
##   Z = -eps * cflow * nc * hcForce(beta=1)     [ = -eps * hcForce(beta=cflow*nc) ]
##
## with `hcForce` the redot-convention gradient of the triangle action
## S = (beta/2) sum_x sum_32 (1 - Re Tr P/N)  (task A, doc/STATUS.md).  Because
## the honeycomb has a different link density (12 links/site vs 4) the naive
## transcription cflow = 1 does NOT give the heat kernel; the induced continuum
## rate is  lambda(p) = c_HC * cflow * p^2 + O(p^4)  with a lattice-geometry
## constant c_HC that FORMULATION.md section 4.4 mandates be fixed numerically.
##
## Measured value (tests/tflow.nim test 2; full numbers in doc/STATUS.md
## task W):  c_HC = 0.1666791 +- 8e-6 after removing the (tiny) O(p^2)
## artifact over three momenta -- i.e. **c_HC = 1/6 and cflow = 6, exactly**
## (pinned at the 4.5e-4 level; FORMULATION 4.4's rough estimate "1/6" was
## right).  The rate scales exactly linearly in cflow (verified to 1.3e-4).
##
## Analytic remark, so nobody "re-derives" a wrong value: linearising the
## flow in the link phases phi_l = A.n_l gives per-link-class rates
##   axis:     dphi/dt = -(cflow/3) (dF)_mu,     (dF)_mu = sum_nu d_nu F_nu mu
##   diagonal: dphi/dt = -(cflow/6) d.(dF)
## and the Rayleigh quotient of the SAMPLED continuum mode is
## (2/9) cflow p^2.  That is only an UPPER bound on the slow rate: the true
## acoustic eigenvector differs from the sampled mode at O(p) (a relative
## axis/diagonal "optical" admixture is allowed within the 24 phases per
## cell), which lowers the eigenvalue at leading order in p^2.  The heat
## kernel measurement is the normative definition; it gives 1/6, not 2/9.
##
## `hcFlowCflow` below is the numerically pinned value; tests/tflow.nim
## re-derives it from scratch (heat-kernel decay at three momenta, O(p^2)
## artifact extrapolated away) and asserts agreement.  Bonus finding: the
## O(p^2) artifact of the honeycomb flow rate is ~35x smaller than the cubic
## one (coefficient -0.002 vs the exact cubic -1/12).
##
## Flow-scale setup (`flowScaleSetup`)
## -----------------------------------
## With `cflow = hcFlowCflow` the flow time `t` handed to `measure` (as
## `wflowT`) is CONTINUUM flow time in units of a^2, directly comparable to
## the cubic-lattice flow time of `gauge/wflow.nim`:  t0 defined by
## t^2<E> = 0.3 on the honeycomb can be compared with cubic t0/a^2 values
## 1:1, and the smoothing radius is sqrt(8t) in both cases.  <E> must be the
## intensive average energy density from `hctopo.hcEQ` (no a^4/2 factor).

import base, layout, field, maths, rng
import physics/qcdTypes
import gauge
import hcgeom, hclayout, hcgauge, hcaction

export hcaction

const hcFlowCflow* = 6.0
  ## Calibrated flow normalisation: with this value the weak-field flow is the
  ## heat kernel exp(-t p^2).  Pinned as the exact rational 6 -- the measured
  ## 1/c_HC = 5.99955 +- 2.8e-4 after O(p^2) extrapolation (tests/tflow.nim
  ## test 2; numbers in doc/STATUS.md task W).

template hcFlowStageU(u, pu, fu: untyped; ca, cb: float; first: static bool) =
  ## one RK stage on one link field: v = cb*f [+ ca*p];  p <- v;  u <- exp(v)*u
  for e in u:
    var v {.noinit.}: type(load1(fu[0]))
    when first:
      v := cb*fu[e]
    else:
      v := cb*fu[e] + ca*pu[e]
    let t = exp(v)*u[e]
    pu[e] := v
    u[e] := t

template hcFlowStage(gg, pp, ff: untyped; ca, cb: float; first: static bool) =
  threads:
    for mu in 0..<nDim:
      hcFlowStageU(gg.uA[mu], pp.uA[mu], ff.uA[mu], ca, cb, first)
      hcFlowStageU(gg.uB[mu], pp.uB[mu], ff.uB[mu], ca, cb, first)
    for d in 0..<nDiag:
      hcFlowStageU(gg.uD[d], pp.uD[d], ff.uD[d], ca, cb, first)

template hcGaugeFlow*(g: HcGauge; steps: int; eps: float; cflow: float;
                      measure: untyped) =
  ## Wilson (triangle-action) flow on the honeycomb.  The input gauge field is
  ## modified in place.  `wflowT` (the flow time after the current step) is
  ## injected for `measure`; `break` inside `measure` stops the flow.
  ## `steps <= 0` flows until `measure` breaks.
  ## `cflow` is the overall normalisation constant (see module docs); pass
  ## `hcFlowCflow` for continuum-normalised flow time.
  block:
    proc flowProc {.gensym.} =
      const nc = g.uA[0][0].nrows.float
      var
        fp = newOneOf(g)          # RK momentum accumulator
        ff = newOneOf(g)          # force
        fw = newHcActionWork(g)
        n = 1
      let betaFlow = cflow*nc
      while true:
        let t = n.float*eps
        hcForce(fw, betaFlow, g, ff)
        hcFlowStage(g, fp, ff, 0.0, -0.25*eps, true)
        hcForce(fw, betaFlow, g, ff)
        hcFlowStage(g, fp, ff, -17.0/9.0, (-8.0/9.0)*eps, false)
        hcForce(fw, betaFlow, g, ff)
        hcFlowStage(g, fp, ff, -1.0, -0.75*eps, false)
        let wflowT {.inject, used.} = t
        measure
        inc n
        if steps > 0 and n > steps:
          break
    flowProc()

template hcGaugeFlow*(g: HcGauge; eps: float; cflow: float;
                      measure: untyped) =
  hcGaugeFlow(g, 0, eps, cflow, measure)

template hcGaugeFlow*(g: HcGauge; eps: float; measure: untyped) =
  ## continuum-normalised flow (cflow = hcFlowCflow)
  hcGaugeFlow(g, 0, eps, hcFlowCflow, measure)

when isMainModule:
  import std/monotimes, std/times
  qexInit()
  let hl = newHcLayout([4, 4, 4, 4])
  var r = hl.lo.newRNGField(RngMilc6, 2468'u64)
  var g = newHcGauge(hl)
  threads:
    g.warm(0.35, r)
  var w = newHcActionWork(g)
  let s0 = hcAction(w, 1.0, g)
  echo "warm 4^4: S(beta=1) = ", s0
  let t0 = getMonoTime()
  var nsteps = 0
  g.hcGaugeFlow(10, 0.02, hcFlowCflow):
    inc nsteps
    if nsteps mod 5 == 0:
      echo "t = ", wflowT, "  S = ", hcAction(w, 1.0, g)
  let t1 = getMonoTime()
  let s1 = hcAction(w, 1.0, g)
  echo "after flow to t = 0.2: S = ", s1, "  (must be < initial ", s0, ")"
  doAssert s1 < s0, "flow must decrease the action"
  let d = g.checkSU
  echo "checkSU after flow: avg ", d.avg, " max ", d.max
  echo "10 RK3 steps on 4^4 cells: ",
       (t1-t0).inMicroseconds.float/1e6, " s total"
  qexFinalize()
