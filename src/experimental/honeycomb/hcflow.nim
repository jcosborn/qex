## Gradient flow of the triangle action and the stout smearing step.
##
## Flow (Luscher, JHEP 1008:071, App. C; the scheme of gauge/wflow.nim):
##   Z_i = -eps force(W_i; beta = cflow N)
##   W1 = exp(Z0/4) W0,  W2 = exp(8/9 Z1 - 17/36 Z0) W1,  V(t+eps) = exp(3/4 Z2 - 8/9 Z1 + 17/36 Z0) W2
## cflow = 6 exactly: linearising the 8 triangle fluxes of a link gives, for
## axis and diagonal links alike, dtheta_l/dt = (cflow/6) n_l . d_nu F_nu mu,
## so cflow = 6 is Luscher's normalisation, dB_mu/dt = D_nu G_nu mu, with t in
## units of a^2 comparable 1:1 with the cubic flow (24 links per unit volume
## instead of 4, and Luscher's S_W sums oriented plaquettes).  The sampled
## plane wave is the exact leading-order eigenmode of the linearised flow;
## tests/tflow.nim measures 1/c_HC = 5.9996(3).  t0 uses the intensive <E> of
## hctopo.EQ.
## RK3 steps: eps <= 0.05 tracks eps = 0.02 on rough 8^4 configurations; eps
## >= 0.1 stays stable but changes which dislocations annihilate (Q history).
##
## Stout (Morningstar-Peardon) with the raw 8-triangle staple sum Sigma
## (= actionDeriv at beta = 2N):
##   U' = exp(-rho projectTAH(U Sigma^dag)) U
## which is one Euler flow step of size rho*stoutKappa, stoutKappa = 2/cflow
## = 1/3 (cubic StoutSmear: 1).  Equal smearing radius sqrt(8 t) therefore
## needs rho_hc = 3 rho_cubic.

import base, layout, field, maths, rng
import physics/qcdTypes
import gauge
import hcgeom, hcgauge, hcaction

export hcaction

const
  cflow* = 6.0
    ## flow normalisation: with force(beta = cflow N) a weak Abelian mode
    ## decays as exp(-t p^2)
  stoutKappa* = 2.0/cflow
    ## heat-kernel time per honeycomb stout step, in units of rho

template flowStage(gg, pp, ff: untyped; ca, cb: float; first: static bool) =
  ## v = cb f [+ ca p];  p <- v;  u <- exp(v) u   on all 24 fields
  threads:
    for i in 0..<nDirs:
      let u = gg.links[i]
      let pu = pp.links[i]
      let fu = ff.links[i]
      for e in u:
        var v {.noinit.}: type(load1(fu[0]))
        when first:
          v := cb*fu[e]
        else:
          v := cb*fu[e] + ca*pu[e]
        let t = exp(v)*u[e]
        pu[e] := v
        u[e] := t

template flow*(g: HcGauge; steps: int; eps: float; c: float;
               measure: untyped) =
  ## RK3 flow in place; `wflowT` (time after the step) is injected for
  ## `measure`; `break` in `measure` stops; steps <= 0 runs until it breaks.
  ## Pass c = cflow for continuum-normalised flow time.
  proc flowProc {.gensym.} =
    const nc = g.uA[0][0].nrows.float
    var
      fp = newOneOf(g)
      ff = newOneOf(g)
      fw = newActionWork(g)
      n = 1
    let betaFlow = c*nc
    while true:
      let t = n.float*eps
      force(fw, betaFlow, g, ff)
      flowStage(g, fp, ff, 0.0, -0.25*eps, true)
      force(fw, betaFlow, g, ff)
      flowStage(g, fp, ff, -17.0/9.0, (-8.0/9.0)*eps, false)
      force(fw, betaFlow, g, ff)
      flowStage(g, fp, ff, -1.0, -0.75*eps, false)
      let wflowT {.inject, used.} = t
      measure
      inc n
      if steps > 0 and n > steps:
        break
  flowProc()

template flow*(g: HcGauge; eps: float; c: float; measure: untyped) =
  flow(g, 0, eps, c, measure)

template flow*(g: HcGauge; eps: float; measure: untyped) =
  ## continuum-normalised flow until `measure` breaks
  flow(g, 0, eps, cflow, measure)

proc stout*[F, SH, SS](w: ActionWork[F, SH, SS], g: HcGauge[F], rho: float,
                       f: var HcGauge[F], n = 1) =
  ## n stout steps in place; `f` is scratch for the staple sums.
  ## Outside `threads:`.
  const nc = g.uA[0][0].nrows
  let ff = f
  let mrho = -rho
  for step in 0..<n:
    actionDeriv(w, 2.0*float(nc), g, f)
    threads:
      for i in 0..<nDirs:
        let u = g.links[i]
        let d = ff.links[i]
        for e in u:
          let sm = u[e]*d[e].adj
          var q {.noinit.}: type(load1(u[0]))
          q.projectTAH sm
          q := exp(mrho*q)
          let t = q*u[e]
          u[e] := t
