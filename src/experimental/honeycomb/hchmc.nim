## HMC for the triangle action.
##
##   H = S(U) + (1/2) sum_l redot(p_l, p_l),   p ~ exp(-redot(p,p)/2)  (randomTAH)
##   T: U_l <- exp(dt p_l) U_l        V: p_l <- p_l - dt f_l   (hcaction.force)
##
## Integrators from mdevolve (leapfrog | 2MN | 4MN5FV, tau split into nsteps);
## accept/reject, reversibility hooks and statistics from hmc/metropolis
## (`update`; fields deltaH, expmDeltaH, accepted, acceptRatio, ...).  The
## accept random number comes from a lattice-global RNG, so the decision is
## identical on all ranks.  Nothing allocates inside `threads:`.

import std/math
import base, layout, field, maths, rng
import physics/qcdTypes
import gauge
import hmc/metropolis
import algorithms/integrator
import hcgeom, hcgauge, hcaction

export hcaction, metropolis

type
  HcHmc*[F, W, RF, GR] = ref object of MetropolisRoot
    beta*, tau*: float
    w*: W                     ## ActionWork
    g*: HcGauge[F]            ## the gauge field being updated
    p*, f*, g0*: HcGauge[F]   ## momenta, force scratch, saved gauge
    r*: RF                    ## per-site RNG field (momenta)
    R*: GR                    ## lattice-global RNG (accept)
    H*: Integrator
    nForce*: int

proc mdt*(h: HcHmc, dt: float) =
  ## g := exp(dt p) g
  expMul(h.g, h.g, h.p, dt)

proc mdv*(h: HcHmc, dt: float) =
  ## p := p - dt force(g)
  force(h.w, h.beta, h.g, h.f)
  inc h.nForce
  let p = h.p
  let f = h.f
  threads:
    for i in 0..<nDirs:
      p.links[i] -= dt*f.links[i]

proc newHcHmc*[F, W, RF, GR](g: HcGauge[F], w: W, beta, tau: float,
                             nsteps: int, r: RF, R: GR,
                             algo = "2MN"): auto =
  ## allocates momenta and two scratch gauges; outside `threads:`
  var h = HcHmc[F, W, RF, GR](beta: beta, tau: tau, w: w, g: g, r: r, R: R)
  var m = MetropolisRoot h
  init(m)
  h.p = newOneOf(g)
  h.f = newOneOf(g)
  h.g0 = newOneOf(g)
  proc vstep(dt: float) = h.mdv(dt)
  proc tstep(dt: float) = h.mdt(dt)
  let (V, T) = newIntegratorPair(vstep, tstep)
  h.H = case algo
    of "leapfrog": mkLeapFrog(T = T, V = V, steps = nsteps)
    of "2MN": mkOmelyan2MN(T = T, V = V, steps = nsteps)
    of "4MN5FV": mkOmelyan4MN5FV(T = T, V = V, steps = nsteps)
    else: toIntegratorProc(algo)(T, V, nsteps)
  h

proc integrate*(h: HcHmc) =
  ## the MD trajectory alone (no momentum refresh, no accept/reject)
  h.H.evolve h.tau
  h.H.finish

proc flipMomenta*(h: HcHmc) =
  let p = h.p
  threads:
    for u in p.links:
      u := -1*u

proc hamiltonian*(h: HcHmc): tuple[s, t, h: float] =
  let s = action(h.w, h.beta, h.g)
  let t = 0.5*redot(h.p, h.p)
  (s, t, s + t)

# MetropolisRoot routines

proc start*(h: HcHmc) =
  let p = h.p
  let g0 = h.g0
  let g = h.g
  threads:
    p.randomTAH h.r
    g0 := g

proc getH*(h: HcHmc): float = h.hamiltonian.h

proc generate*(h: HcHmc) = h.integrate

proc globalRand*(h: HcHmc): float = float(h.R.uniform)

proc accept*(h: HcHmc) =
  let g = h.g
  threads:
    g.reunit

proc reject*(h: HcHmc) =
  let g = h.g
  let g0 = h.g0
  threads:
    g := g0

proc revCheck*(h: HcHmc): tuple[dHf, dHb, sumdH, linkDiff: float] =
  ## fresh momenta, forward, flip, back; returns (dH_fwd, dH_bwd,
  ## dH_fwd + dH_bwd, sqrt(sum |U_back - U_0|^2 / nLinks)); g restored
  h.start
  let h0 = h.getH
  h.integrate
  let h1 = h.getH
  h.flipMomenta
  h.integrate
  let h2 = h.getH
  let g = h.g
  let g0 = h.g0
  var d2 = 0.0
  threads:
    let d = norm2diff(g, g0)
    threadMaster: d2 = d
    g := g0
  (h1 - h0, h2 - h1, h2 - h0, sqrt(d2/float(nDirs*g.lo.physVol)))
