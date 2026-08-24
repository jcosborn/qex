## Tree-level clover (Sheikholeslami-Wohlert) improvement of the Wilson-Dirac
## operator, on the 16-cell honeycomb and on the cubic lattice (task **D2**,
## part 2).
##
## The improved operator
## ---------------------
##   D_c = D  -  (c_SW r / 4) sigma_munu F_munu           (a = 1, sum over ALL
##                                                          ordered pairs mu,nu)
## with sigma_munu = (i/2)[gamma_mu, gamma_nu] (Hermitian) and F_munu the
## HERMITIAN field strength of the convention U = exp(+i int A.dl) that the
## whole honeycomb stack uses (hctopo pins Fhat = +i a^2 F T on Abelian
## backgrounds).  In terms of the code's traceless ANTI-Hermitian clover field
## Fhat_ab = +i a^2 F_ab (stored for a > b, `pairIndex` order) this collapses,
## using sigma_munu F_munu = gamma_mu gamma_nu Fhat_munu (mu != nu), to
##
##   D_c psi(x) = D psi(x)  -  (c_SW r / 2) sum_{a>b} gamma_a gamma_b
##                                             Fhat_ab(x) psi(x)      ... (*)
##
## Why this exact coefficient (derivation, both lattices at once).  The naive
## term of either operator is gamma.D + O(a^2) (the O(a) piece is odd in the
## neighbour vector n and cancels in the +-n pairs), and the Wilson term is
##   cubic:      -(r a/2) (1/2) sum_{8}  (n_i.D)^2 = -(r a/2) D^2 + O(a^2)
##   honeycomb:  -(r a/2) (1/6) sum_{24} (n_i.D)^2 = -(r a/2) D^2 + O(a^2)
## (the same "a r p^2/2" normalisation, FORMULATION 5.1, pinned numerically by
## tasks F/D1), with NO F term of its own: n_mu n_nu is symmetric, so (n.D)^2
## only contains {D_mu, D_nu}.  Since (gamma.D)^2 = D^2 + (1/2) sigma_munu
## F_munu ([D_mu,D_nu] = i F_munu in this U convention), adding
## -(c_SW r a/4) sigma.F with c_SW = 1 turns the lattice operator into
##   m + gamma.D - (r a/2)(gamma.D)^2 + O(a^2),
## whose O(a) piece is a pure on-shell mass shift: the standard tree-level
## O(a) improvement.  BOTH lattices share the -(r/4) coefficient because they
## share the Wilson-term normalisation -- there is no honeycomb-specific
## factor (unlike the flow constant cflow).  The prefactor of (*) is pinned
## numerically, convention-free, by tests/tclover.nim test 5: on the exact
## constant-flux Atiyah-Singer background the measured on-site matrix is
##   (D_c - D)(x) = (i/2)[f1 s1 gamma2gamma1 + f2 s2 gamma4gamma3] (x) T
## with s = the exact clover artifact factor, s_hc = 4 sin(f/4)/f and
## s_cubic = sin(f)/f -- measured to ~1e-12, artifact 1 - s = O(1/L^4).
##
## Structure of the added term: it is site-diagonal (on-site), gauge
## covariant, HERMITIAN ((gamma_a gamma_b)^dag = -gamma_a gamma_b and
## Fhat^dag = -Fhat), and commutes with gamma5 -- hence D_c stays
## gamma5-Hermitian, gamma5 D_c gamma5 = D_c^dag.  (For a self-dual
## background F_01 = F_23 the term vanishes on one chirality: sigma_01
## sigma_23 = -gamma5, so the gamma5 = +1 spin sector gets zero shift for
## f1 = f2 -- the chirality-splitting cross-check of tests/tclover.nim.)
##
## Honeycomb field strength: `hctopo.hcFmunu` (hexagon clover, normalisation
## Fhat = a^2 F verified by task W).  Cubic: QEX `fmunu(g, 1)` (1x1 clover);
## its storage fm[a][b] (a > b) was measured on the constant-flux background
## to equal +i sin(phi_ab) T = +i F_ab T + O(F^3), i.e. the SAME +i F_ab
## convention as hctopo pair storage -- no sign flip needed
## (`cubicFmunuSign = +1`, asserted in tests/tclover.nim test 5).
##
## The clover field (6 colour-matrix fields per sublattice) is precomputed
## ONCE per configuration by `gaugeRefresh`; each D_c application then costs
## the bare D plus 6 on-site colour matrix-vector products per site, all
## allocation free.  `r_w` (the Wilson parameter) multiplies the clover term
## at apply time (the tree-level improvement is proportional to r); the cubic
## wrapper wraps QEX `physics/wilsonD`, which hardwires r_w = 1.
##
## Boundary conditions: the (anti)periodic fermion BC sign flips (`setBC`)
## cancel in every closed clover loop (each loop crosses the time cut an even
## number of times), so F computed before or after `setBC` is identical; call
## `gaugeRefresh` after any change of the links all the same (smearing DOES
## change F).

import base, layout, field, maths
import physics/qcdTypes
import physics/wilsonD
import gauge
import hcgeom, hclayout, hcgauge, hcwilson, hctopo

export hcwilson, hctopo
export wilsonD

const cubicFmunuSign* = 1.0
  ## QEX fmunu(g,1)[a][b] = +i F_ab T (measured, constant-flux background;
  ## asserted in tests/tclover.nim test 5).  Kept as a named constant so a
  ## core-QEX convention change fails loudly there.

# ---------------------------------------------------------------------------
# the six spin matrices gamma_a gamma_b (a > b, hctopo.pairIndex order)
# ---------------------------------------------------------------------------

proc buildCloverGammas(): array[6, SpinMat] =
  var t: SpinMat
  template setp(a, b: int, ga, gb: typed) =
    t := ga*gb
    result[pairIndex(a, b)] = t
  setp(1, 0, gamma2, gamma1)
  setp(2, 0, gamma3, gamma1)
  setp(2, 1, gamma3, gamma2)
  setp(3, 0, gamma4, gamma1)
  setp(3, 1, gamma4, gamma2)
  setp(3, 2, gamma4, gamma3)

let hcCloverGam* = buildCloverGammas()
  ## hcCloverGam[pairIndex(a,b)] = gamma_a gamma_b (0-based directions,
  ## DeGrand-Rossi gammas; anti-Hermitian for a != b)

# ---------------------------------------------------------------------------
# honeycomb
# ---------------------------------------------------------------------------

type
  HcCloverWilson*[W, TW] = ref object
    ## Clover-improved honeycomb Wilson operator: wraps a bare `HcWilson`
    ## (which owns the gauge ref) plus an `HcTopoWork` whose f[sub][pair]
    ## fields hold the precomputed hexagon-clover Fhat.  Construct once with
    ## `newHcCloverWilson`; call `gaugeRefresh` whenever the links change
    ## in place (smearing, setBC).  ref object: `threads:` can capture it.
    w*: W          ## the bare HcWilson operator
    tw*: TW        ## HcTopoWork; tw.f[sub][pairIndex(a,b)] = Fhat_ab
    cSW*: float

proc newHcCloverWilson*[V: static[int], MF](g: HcGauge[V, MF],
                                            cSW: float): auto =
  ## Build the improved operator on `g` (allocates; outside `threads:`).
  ## The clover field is computed from the current links.
  var w = newHcWilson(g)
  var tw = newHcTopoWork(g)
  var c = HcCloverWilson[typeof(w), typeof(tw)](w: w, tw: tw, cSW: cSW)
  hcFmunu(tw, g)
  c

proc gaugeRefresh*(c: HcCloverWilson) =
  ## Recompute the pre-shifted links of the bare operator AND the clover
  ## field, from the current contents of the gauge field the operator was
  ## built on.  Call after smearing into it, after `setBC`, etc.
  ## Outside `threads:` (opens its own).
  c.w.gaugeRefresh
  hcFmunu(c.tw, c.w.g)

proc applyClover(c: HcCloverWilson, r: HcFermion, x: HcFermion, cf: float) =
  ## r += cf * sum_p gam_p (Fhat_p x) on both sublattices; opens `threads:`.
  var gm {.noinit.}: array[6, SpinMat]
  for p in 0..<6:
    gm[p] := cf*hcCloverGam[p]
  threads:
    for p in 0..<6:
      for e in r.a:
        r.a[e] += gm[p] * (c.tw.f[0][p][e] * x.a[e])
      for e in r.b:
        r.b[e] += gm[p] * (c.tw.f[1][p][e] * x.b[e])

proc D*(c: HcCloverWilson, r: var HcFermion, x: HcFermion,
        m: float, rw: float = 1.0) =
  ## r = D x - (cSW rw/2) sum_{a>b} gamma_a gamma_b Fhat_ab x, the
  ## tree-level-improved operator (*) (cSW = 1 = standard normalisation).
  ## Opens its own `threads:` blocks; r must not alias x.
  c.w.D(r, x, m, rw)
  if c.cSW != 0.0:
    applyClover(c, r, x, -0.5*c.cSW*rw)

proc Ddag*(c: HcCloverWilson, r: var HcFermion, x: HcFermion,
           m: float, rw: float = 1.0) =
  ## Hermitian conjugate of `D`: the clover term is Hermitian and commutes
  ## with gamma5, so it enters D^dag unchanged (tests/tclover.nim test 3).
  c.w.Ddag(r, x, m, rw)
  if c.cSW != 0.0:
    applyClover(c, r, x, -0.5*c.cSW*rw)

# ---------------------------------------------------------------------------
# cubic
# ---------------------------------------------------------------------------

type
  CubicCloverWilson*[W, MF] = ref object
    ## Same shape on the cubic lattice: wraps QEX `physics/wilsonD.Wilson`
    ## (which owns g and hardwires r_w = 1) plus the 6 precomputed clover
    ## fields from QEX `fmunu(g, 1)` in hctopo pair order/convention.
    s*: W               ## physics/wilsonD Wilson operator (holds s.g)
    f*: array[6, MF]    ## Fhat_ab (a > b, pairIndex order), = +i F_ab T conv.
    cSW*: float

proc gaugeRefresh*(c: CubicCloverWilson) =
  ## Recompute the clover field from the current links c.s.g.  ALLOCATES
  ## (QEX fmunu builds work fields); once per configuration, outside
  ## `threads:`.
  let fm = fmunu(c.s.g, 1)
  threads:
    for a in 1..<4:
      for b in 0..<a:
        c.f[pairIndex(a, b)] := cubicFmunuSign*fm[a][b]

proc newCubicCloverWilson*[G](g: seq[G], cSW: float): auto =
  ## VLEN-layout constructor (lo.DiracFermion fermions, like newWilson(g)).
  var s = newWilson(g)
  var c = CubicCloverWilson[typeof(s), G](s: s, cSW: cSW)
  for p in 0..<6:
    c.f[p] = g[0].newOneOf
  c.gaugeRefresh
  c

proc newCubicCloverWilson*[G, T](g: seq[G], v: T, cSW: float): auto =
  ## generic-V constructor; `v` = prototype fermion field (the tarnoldi
  ## makeWilson pattern for V = 1 layouts).
  var s = newWilson(g, v)
  var c = CubicCloverWilson[typeof(s), G](s: s, cSW: cSW)
  for p in 0..<6:
    c.f[p] = g[0].newOneOf
  c.gaugeRefresh
  c

proc applyCloverCubic(c: CubicCloverWilson, r: auto, x: auto, cf: float) =
  var gm {.noinit.}: array[6, SpinMat]
  for p in 0..<6:
    gm[p] := cf*hcCloverGam[p]
  threads:
    for p in 0..<6:
      for e in r:
        r[e] += gm[p] * (c.f[p][e] * x[e])

proc D*(c: CubicCloverWilson, r: var auto, x: auto, m: SomeNumber) =
  ## r = (4+m) x - (1/2) sum_{+-mu} (r_w -+ gamma_mu) U x(x+-mu)
  ##     - (cSW/2) sum_{a>b} gamma_a gamma_b Fhat_ab x        (r_w = 1)
  ## Opens its own `threads:` blocks; r must not alias x.
  let rr = r                  # Field is a ref: value copy for threads capture
  threads:
    c.s.D(rr, x, m)
  if c.cSW != 0.0:
    applyCloverCubic(c, rr, x, -0.5*c.cSW)

proc Ddag*(c: CubicCloverWilson, r: var auto, x: auto, m: SomeNumber) =
  let rr = r
  threads:
    c.s.Ddag(rr, x, m)
  if c.cSW != 0.0:
    applyCloverCubic(c, rr, x, -0.5*c.cSW)

when isMainModule:
  qexInit()
  echo "hcclover: pair gammas built; cubicFmunuSign = ", cubicFmunuSign
  block:                        # honeycomb smoke: unit gauge, D_c == D form
    let hl = newHcLayout([4, 4, 4, 4])
    var g = newHcGauge(hl)
    var c = newHcCloverWilson(g, 1.0)
    var x = newHcFermion(hl)
    var r = newHcFermion(hl)
    x.a := 1
    c.D(r, x, 0.1)
    echo "hc unit gauge, constant psiA: |D_c psi|^2 = ", r.norm2
  block:                        # cubic smoke
    let lo = newLayout(@[4, 4, 4, 4])
    var g = lo.newGauge
    var c = newCubicCloverWilson(g, 1.0)
    var x = lo.DiracFermion()
    var r = lo.DiracFermion()
    threads:
      x := 1
    c.D(r, x, 0.1)
    echo "cubic unit gauge, constant psi: |D_c psi|^2 = ", r.norm2
  qexFinalize()
