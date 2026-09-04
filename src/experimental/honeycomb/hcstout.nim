## Stout smearing on the 16-cell honeycomb and repeated steps for
## QEX's cubic `gauge/stoutsmear` (task **D2**, part 1).
##
## Morningstar-Peardon (hep-lat/0311018) stout step, per link:
##
##   C_l     = rho * Sigma_l          Sigma_l = the plain staple sum
##   Omega_l = C_l U_l^dag
##   Q_l     = projectTAH(Omega_l)    (traceless anti-Hermitian, = "i Q_MP")
##   U'_l    = exp(Q_l) U_l
##
## Honeycomb staples.  `Sigma_l` is the sum of the 8 triangle staples of the
## link, in QEX convention (loop P_k = U_l V_k^dag, Sigma = sum_k V_k): exactly
## `hcActionDeriv` at `beta = 2 N` (its prefactor is beta/2N, see hcaction.nim
## and the task-M pin `staples == hcActionDeriv(2N)` to 1.7e-32).  Since
## projectTAH(C U^dag) = -rho * projectTAH(U Sigma^dag) = -rho * hcForce(2N),
## the step is
##
##   U'_l = exp( -rho * projectTAH( U_l Sigma_l^dag ) ) U_l ,
##
## the gradient-descent direction of the triangle action: smearing INCREASES
## the triangle sum (tests/tstout.nim test 4).
##
## Cubic smearing.  QEX's `StoutSmear.smear` computes
## U' = exp(-alpha*nc * projectTAH(U ds^dag)) U with
## ds = gaugeActionDeriv(plaq=1) = (1/nc) * (plaquette staple sum), i.e.
## U' = exp(projectTAH(C U^dag)) U with C = rho * (staple sum),
## the isotropic MP step with rho = alpha.  The `smearN` overload below
## repeats this step; `StoutSmear.alpha` controls its strength.
##
## Heat-kernel normalisation (the numbers task D4 needs)
## -----------------------------------------------------
## For a weak field, n stout steps act on every gauge mode as the Wilson-flow
## heat kernel exp(-t_eff p^2) with
##
##   t_eff = n * rho * kappa ,     kappa_cubic     = 1      (exact)
##                                 kappa_honeycomb = 1/3    (exact, given cflow=6)
##
## Cubic: one MP stout step is exactly one Euler step of QEX's `gaugeFlow`
## (Z = -eps*nc*gaugeForce(plaq=1)) with eps = rho, and that flow's t is
## continuum heat-kernel time (tests/tflow.nim test 1; Luscher 2010: the
## Wilson flow is the continuous-stout limit).  Hence kappa = 1.
##
## Honeycomb: one stout step is Z = -rho*hcForce(beta = 2 nc), while the
## continuum-normalised flow (hcflow.nim) uses Z = -eps*hcForce(beta = cflow*nc)
## with the pinned cflow = 6.  hcForce is linear in beta, so one stout step
## = one Euler flow step of size eps = rho * (2/cflow) = rho/3.  Hence
## kappa = 2/cflow = 1/3 exactly, inheriting the cflow = 6 pin (2.3e-5).
##
## Both kappas are re-measured from scratch in tests/tstout.nim test 5 (plane
## wave from exact line integrals, three momenta, O(p^2) artifact removed).
##
## **rho equivalence for task D4**: at the same rho the honeycomb smears 3x
## less than the cubic lattice.  Equal smearing radius sqrt(8 t_eff) needs
## rho_hc = 3 * rho_cubic (or 3x the steps).  The paper's "6 steps at
## rho = 0.05 on both lattices", read with the plain MP staple-sum convention
## above, therefore means t_eff = 0.30 a^2 (radius 1.55 a) on the cubic
## lattice but t_eff = 0.10 a^2 (radius 0.89 a) on the honeycomb.
##
## Usage: create once (allocates), then `smear`/`smearN` are allocation-free
## (honeycomb) and open their own `threads:` blocks -- do NOT call from inside
## one.  `gout` may be the same object as `gin` (in-place smearing).
## No force/chain rule (`smearDeriv`) on the honeycomb: quenched measurements
## only (`StoutSmear` provides QEX's cubic chain rule).

import base, layout, field, maths
import physics/qcdTypes
import gauge
import gauge/stoutsmear
import hcgeom, hclayout, hcgauge, hcaction

export hcaction
export stoutsmear

const hcStoutKappa* = 1.0/3.0
  ## t_eff per honeycomb stout step, in units of rho: exact 2/cflow = 1/3
  ## (hcflow.hcFlowCflow = 6); re-measured in tests/tstout.nim test 5.

const cubicStoutKappa* = 1.0
  ## t_eff per cubic (QEX StoutSmear) step in units of rho: exact 1.

# ---------------------------------------------------------------------------
# honeycomb
# ---------------------------------------------------------------------------

type
  HcStout*[V: static[int], F, W] = ref object
    ## Persistent stout-smearing work space: the staple machinery of task A
    ## (`HcActionWork`) plus one scratch gauge for the staple sums.
    ## Create once with `newHcStout`; each `smear` call is then allocation
    ## free.  Rebound to the gauge passed to each call (any same-shape gauge).
    rho*: float
    w*: W                    ## HcActionWork (staple shifters + scratch)
    d*: HcGauge[V, F]        ## staple sums Sigma_l (24 fields)

proc newHcStout*[V: static[int], F](g: HcGauge[V, F], rho: float): auto =
  ## Allocate the work space for gauges shaped like `g`.  Outside `threads:`.
  var w = newHcActionWork(g)
  var d = newOneOf(g)
  HcStout[V, F, typeof(w)](rho: rho, w: w, d: d)

proc smear*[V: static[int], F, W](s: HcStout[V, F, W], gin: HcGauge[V, F],
                                  gout: var HcGauge[V, F]) =
  ## One Morningstar-Peardon stout step with the 8-triangle staples:
  ## gout_l = exp(-rho * projectTAH(gin_l Sigma_l^dag)) gin_l.
  ## Opens its own `threads:`.  `gout` may alias `gin` (the staple sums are
  ## fully computed first; the update is site-local).  `gin` must not be
  ## `s.d` itself.
  const nc = gin.uA[0][0].nrows
  hcActionDeriv(s.w, 2.0*float(nc), gin, s.d)    # s.d = raw staple sums Sigma
  let sd = s.d
  let go = gout            # value copy: the 24 fields are refs (hcaction style)
  let mrho = -s.rho
  threads:
    template upd(uin, ud, uout: untyped) =
      for e in uin:
        let sm = uin[e] * ud[e].adj
        var q {.noinit.}: type(load1(uin[0]))
        q.projectTAH sm
        q := exp(mrho*q)
        let t = q*uin[e]
        uout[e] := t
    for mu in 0..<nDim:
      upd(gin.uA[mu], sd.uA[mu], go.uA[mu])
      upd(gin.uB[mu], sd.uB[mu], go.uB[mu])
    for d in 0..<nDiag:
      upd(gin.uD[d], sd.uD[d], go.uD[d])

proc smearN*[V: static[int], F, W](s: HcStout[V, F, W], g: HcGauge[V, F],
                                   gout: var HcGauge[V, F], n: int) =
  ## n stout steps; n = 0 copies g to gout.  Opens its own `threads:` blocks.
  if n <= 0:
    let go = gout            # value copy: threads cannot capture var params
    threads:
      go := g
  else:
    smear(s, g, gout)
    for i in 1..<n:
      smear(s, gout, gout)

# ---------------------------------------------------------------------------
# cubic
# ---------------------------------------------------------------------------

proc smearN*[G](s: var StoutSmear[G], g: G, gout: G, n: int) =
  ## n stout steps; n = 0 copies g to gout.
  if n <= 0:
    threads:
      for mu in 0..<gout.len:
        gout[mu] := g[mu]
  else:
    s.smear(g, gout)
    for i in 1..<n:
      s.smear(gout, gout)

when isMainModule:
  import std/monotimes, std/times
  import rng
  qexInit()
  echo "hcstout: kappa_hc = ", hcStoutKappa, "  kappa_cubic = ", cubicStoutKappa
  block:                        # honeycomb smoke + timing, 8^4 cells
    let hl = newHcLayout([8, 8, 8, 8])
    var r = hl.lo.newRNGField(RngMilc6, 31415'u64)
    var g = newHcGauge(hl)
    threads:
      g.warm(0.35, r)
    var aw = newHcActionWork(g)
    var st = newHcStout(g, 0.05)
    echo "hc 8^4 warm: triangleSum(0) = ", g.triangleSum
    st.smear(g, g)              # warm up buffers
    let n = 10
    let t0 = getMonoTime()
    for i in 0..<n:
      st.smear(g, g)
    let t1 = getMonoTime()
    echo "hc 8^4: triangleSum(11 steps) = ", g.triangleSum
    let d = g.checkSU
    echo "hc checkSU: avg ", d.avg, " max ", d.max
    echo "hc 8^4 cells: smear ", (t1-t0).inMicroseconds.float/(1e3*n.float),
         " ms/step"
    discard hcAction(aw, 1.0, g)
  block:                        # cubic smoke + timing, 8^4 sites
    let lo = newLayout(@[8, 8, 8, 8])
    var r = lo.newRNGField(RngMilc6, 27182'u64)
    var g = lo.newGauge
    threads:
      g.warm(0.35, r)
    var st = newStoutSmear(lo, 0.05)
    echo "cubic 8^4: plaq(0) = ", g.plaq.sum
    st.smear(g, g)
    let n = 10
    let t0 = getMonoTime()
    for i in 0..<n:
      st.smear(g, g)
    let t1 = getMonoTime()
    echo "cubic 8^4: plaq(11 steps) = ", g.plaq.sum
    echo "cubic 8^4: smear ", (t1-t0).inMicroseconds.float/(1e3*n.float),
         " ms/step"
  qexFinalize()
