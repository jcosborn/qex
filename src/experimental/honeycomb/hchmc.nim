## HMC for the triangle gauge action on the 16-cell honeycomb (task **M1**).
##
## Structure and conventions follow `src/examples/puregaugehmc.nim` /
## `refCubicGen.nim`, adapted to the 24 link fields of `HcGauge`:
##
##   momenta   one traceless anti-Hermitian field per link field
##             (an `HcGauge` holding TAH matrices, drawn by `randomTAH`, which
##             fills each matrix with N^2-1 unit gaussians in an orthonormal
##             TAH basis, i.e. P(p) ~ exp(-|p|^2/2) = exp(-T(p)))
##   T(p)    = (1/2) sum_l redot(p_l, p_l) = (1/2) sum_l |p_l|^2
##   H       = T(p) + S(U),   S = hcaction.hcAction
##   update    U_l <- exp(dt p_l) U_l          (mdt)
##             p_l <- p_l - dt f_l             (mdv, f = hcaction.hcForce)
##
## Consistency of T with the Task A force convention (verified in
## tests/thmc.nim): hcForce satisfies d/ds S(exp(s P)U)|_0 = sum_l redot(P_l,f_l)
## and randomTAH draws p from exp(-redot(p,p)/2), so the leapfrog above
## conserves H = S + (1/2) redot(p,p) and `<exp(-dH)> = 1` holds by the
## standard argument.  This is exactly the puregaugehmc.nim normalisation
## (there written `0.5*p2` with `p2 = sum norm2`).
##
## Integrators (per MD step of size eps = tau/nsteps; T = position update,
## V = momentum update; adjacent same-kind stages across step boundaries are
## merged, so the force counts below are per trajectory):
##
##   "leapfrog"  T(1/2) V(1) T(1/2)                        nsteps   force calls
##   "2MN"       Omelyan 2nd order minimum norm,           2*nsteps
##               T(lam) V(1/2) T(1-2lam) V(1/2) T(lam),
##               lam = 0.1931833275037836   (Omelyan et al 2003, eq 31;
##               same stage list as QEX's mdevolve mkOmelyan2MN)
##   "4MN5FV"    Omelyan 4th order, 5 force evals,         4*nsteps+1
##               velocity version (mdevolve mkOmelyan4MN5FV coefficients);
##               the most cost-effective choice found for the cubic reference
##               runs (STATUS.md task C)
##
## The Metropolis accept/reject uses a *lattice-global* random number from a
## caller-supplied global RNG (e.g. `MRG32k3a` seeded with `seed(seed,ix)`,
## which broadcasts the seed from rank 0), so the decision is identical on all
## ranks -- the puregaugehmc.nim pattern.
##
## Nothing here allocates inside a `threads:` block; all fields and shifters
## live in the `HcHmc` object, created once by `newHcHmc`.

import std/math
import base, layout, field, maths, rng
import physics/qcdTypes
import gauge
import hcgeom, hclayout, hcgauge, hcaction

export hcaction

const omelyan2MNLambda* = 0.1931833275037836

type
  MdStage* = tuple[isV: bool, c: float]

  HcHmc*[V: static[int], F, W] = object
    ## Persistent HMC state: parameters, MD schedule, the shared
    ## `HcActionWork` (shifters + scratch), momenta and scratch gauges.
    ## Create with `newHcHmc` (allocates); everything per trajectory is then
    ## allocation free.  Not bound to a particular gauge field: `trajectory`
    ## takes `g` explicitly and works for any same-shape `HcGauge`.
    beta*, tau*: float
    nsteps*: int
    algo*: string
    sched*: seq[MdStage]
    w*: W                     ## HcActionWork
    p*: HcGauge[V, F]         ## momenta (TAH)
    f*: HcGauge[V, F]         ## force scratch
    g0*: HcGauge[V, F]        ## saved gauge for reject / reversibility
    nForce*: int              ## force evaluations so far (diagnostics)

proc mdSchedule*(algo: string, nsteps: int): seq[MdStage] =
  ## Build the merged stage list for `nsteps` MD steps of the named
  ## integrator.  Coefficients are in units of eps = tau/nsteps; the T (=
  ## position) coefficients of one step sum to 1.
  doAssert nsteps > 0, "hchmc: nsteps must be positive"
  var stage: seq[MdStage]
  case algo
  of "leapfrog":
    stage = @[(false, 0.5), (true, 1.0), (false, 0.5)]
  of "2MN", "omelyan":
    const l = omelyan2MNLambda
    stage = @[(false, l), (true, 0.5), (false, 1.0 - 2.0*l), (true, 0.5),
              (false, l)]
  of "4MN5FV":
    const
      rho = 0.2539785108410595
      theta = -0.03230286765269967
      vartheta = 0.08398315262876693
      lam = 0.6822365335719091
    stage = @[(true, vartheta), (false, rho), (true, lam), (false, theta),
              (true, 0.5 - (lam + vartheta)), (false, 1.0 - 2.0*(theta + rho)),
              (true, 0.5 - (lam + vartheta)), (false, theta), (true, lam),
              (false, rho), (true, vartheta)]
  else:
    qexError "hchmc: unknown integrator '" & algo &
      "' (leapfrog, 2MN, 4MN5FV)"
  for k in 0..<nsteps:
    for s in stage:
      if result.len > 0 and result[^1].isV == s.isV:
        result[^1].c += s.c
      else:
        result.add s
  # sanity: the position coefficients must sum to nsteps (i.e. to tau)
  var st = 0.0
  for s in result:
    if not s.isV: st += s.c
  doAssert abs(st - nsteps.float) < 1e-12*nsteps.float

proc nForcePerTraj*(sched: seq[MdStage]): int =
  for s in sched:
    if s.isV: inc result

proc newHcHmc*[V: static[int], F](g: HcGauge[V, F], beta, tau: float,
                                  nsteps: int, algo = "2MN"): auto =
  ## Allocate the HMC state for gauge fields shaped like `g`.
  ## Call outside `threads:`.
  var w = newHcActionWork(g)
  var h = HcHmc[V, F, type(w)](beta: beta, tau: tau, nsteps: nsteps,
                               algo: algo, w: w)
  h.sched = mdSchedule(algo, nsteps)
  h.p = newOneOf(g)
  h.f = newOneOf(g)
  h.g0 = newOneOf(g)
  h

# ---------------------------------------------------------------------------
# MD updates
# ---------------------------------------------------------------------------

proc mdt*[V: static[int], F, W](h: var HcHmc[V, F, W], g: HcGauge[V, F],
                                dt: float) =
  ## g := exp(dt*p) * g on all 24 link fields
  let p = h.p
  threads:
    for mu in 0..<nDim:
      for e in g.uA[mu]:
        g.uA[mu][e] := exp(dt*p.uA[mu][e])*g.uA[mu][e]
      for e in g.uB[mu]:
        g.uB[mu][e] := exp(dt*p.uB[mu][e])*g.uB[mu][e]
    for d in 0..<nDiag:
      for e in g.uD[d]:
        g.uD[d][e] := exp(dt*p.uD[d][e])*g.uD[d][e]

proc mdv*[V: static[int], F, W](h: var HcHmc[V, F, W], g: HcGauge[V, F],
                                dt: float) =
  ## p := p - dt * hcForce(g)
  hcForce(h.w, h.beta, g, h.f)
  inc h.nForce
  let p = h.p
  let f = h.f
  threads:
    for mu in 0..<nDim:
      p.uA[mu] -= dt*f.uA[mu]
      p.uB[mu] -= dt*f.uB[mu]
    for d in 0..<nDiag:
      p.uD[d] -= dt*f.uD[d]

proc integrate*[V: static[int], F, W](h: var HcHmc[V, F, W],
                                      g: HcGauge[V, F]) =
  ## Run the MD schedule over trajectory length tau (no momentum refresh, no
  ## accept/reject) -- the piece the reversibility test drives directly.
  let eps = h.tau/h.nsteps.float
  for s in h.sched:
    if s.isV: h.mdv(g, s.c*eps)
    else: h.mdt(g, s.c*eps)

proc flipMomenta*[V: static[int], F, W](h: var HcHmc[V, F, W]) =
  let p = h.p
  threads:
    eachLink(p, u):
      u := -1*u

proc kinetic*[V: static[int], F, W](h: HcHmc[V, F, W]): float =
  ## T(p) = (1/2) sum_l |p_l|^2
  0.5*redot(h.p, h.p)

proc hamiltonian*[V: static[int], F, W](h: var HcHmc[V, F, W],
                                        g: HcGauge[V, F]):
    tuple[s, t, h: float] =
  let s = hcAction(h.w, h.beta, g)
  let t = kinetic(h)
  (s, t, s + t)

# ---------------------------------------------------------------------------
# one trajectory
# ---------------------------------------------------------------------------

proc trajectory*[V: static[int], F, W, RF, GR](h: var HcHmc[V, F, W],
                                               g: HcGauge[V, F],
                                               r: RF, R: var GR,
                                               alwaysAccept = false):
    tuple[dH, acc: float, accepted: bool] =
  ## One HMC trajectory: refresh momenta from `r` (per-site RNG field),
  ## integrate, and Metropolis-accept with a single *lattice-global* uniform
  ## from `R` (identical on all ranks: seed `R` with the broadcast `seed`
  ## proc, and call `trajectory` the same number of times everywhere).
  ## On accept the gauge is reunitarised (projectSU); on reject it is restored
  ## from the saved copy.  Returns (dH, exp(-dH), accepted).
  let p = h.p
  let g0 = h.g0
  threads:
    eachLink(p, u):
      u.randomTAH r      # gaugeUtils.randomTAH(Field, RNGField)
    g0 := g
  let (_, _, h0) = h.hamiltonian(g)
  h.integrate(g)
  let (_, _, h1) = h.hamiltonian(g)
  let dH = h1 - h0
  let acc = exp(-dH)
  let rnd = R.uniform    # lattice global: same on every rank
  let accepted = (rnd.float <= acc) or alwaysAccept
  if accepted:
    threads:
      g.reunit
  else:
    threads:
      g := g0
  (dH, acc, accepted)

proc revCheck*[V: static[int], F, W, RF](h: var HcHmc[V, F, W],
                                         g: HcGauge[V, F], r: RF):
    tuple[dHf, dHb, sumdH, linkDiff: float] =
  ## Reversibility diagnostic: fresh momenta, integrate forward, flip p,
  ## integrate back; returns (dH_fwd, dH_bwd, dH_fwd+dH_bwd,
  ## sqrt(sum|U_back-U_0|^2/nLinks)).  `g` is restored to its initial state
  ## (up to the roundtrip roundoff just measured).
  let p = h.p
  let g0 = h.g0
  threads:
    eachLink(p, u):
      u.randomTAH r
    g0 := g
  let (_, _, h0) = h.hamiltonian(g)
  h.integrate(g)
  let (_, _, h1) = h.hamiltonian(g)
  h.flipMomenta
  h.integrate(g)
  let (_, _, h2) = h.hamiltonian(g)
  # per-link deviation from the starting configuration
  var d2 = 0.0
  let f = h.f
  threads:
    for mu in 0..<nDim:
      f.uA[mu] := g.uA[mu] - g0.uA[mu]
      f.uB[mu] := g.uB[mu] - g0.uB[mu]
    for d in 0..<nDiag:
      f.uD[d] := g.uD[d] - g0.uD[d]
  d2 = redot(f, f)
  threads:
    g := g0
  (h1 - h0, h2 - h1, h2 - h0, sqrt(d2/float(nLinks(g.hl))))

when isMainModule:
  import std/[strformat, monotimes, times]
  qexInit()
  echo "hchmc smoke test + timing"
  let hl = newHcLayout([4, 4, 4, 4])
  var r = hl.lo.newRNGField(MRG32k3a, 987654321'u64)
  var R: MRG32k3a
  R.seed(987654321'u64, 17'u64)
  var g = newHcGauge(hl)
  threads:
    g.warm(0.35, r)
  var h = newHcHmc(g, 6.0, 1.0, 8, "2MN")
  echo &"schedule {h.algo} nsteps {h.nsteps}: {nForcePerTraj(h.sched)} force calls/traj"
  for n in 1..5:
    let t0 = getMonoTime()
    let (dH, acc, accepted) = h.trajectory(g, r, R)
    let t1 = getMonoTime()
    echo &"traj {n}: dH {dH: .6e}  exp(-dH) {acc:.6f}  acc {accepted}  ",
         &"S/vol {hcAction(h.w, h.beta, g)/float(nLinks(hl)):.6f}  ",
         &"{(t1-t0).inMicroseconds.float*1e-6:.3f} s"
  let rc = h.revCheck(g, r)
  echo &"revcheck: dHf {rc.dHf: .3e}  dHb {rc.dHb: .3e}  sum {rc.sumdH: .3e}  linkdiff {rc.linkDiff:.3e}"
  qexFinalize()
