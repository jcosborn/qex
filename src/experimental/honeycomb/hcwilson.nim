## Interacting Wilson-Dirac operator on the 16-cell honeycomb (task **D1**).
##
## The operator (doc/FORMULATION.md section 5.1, `a = 1`):
##
##   D psi(x) = (m + 4 r) psi(x) + (1/6) sum_{i=1}^{24} (gamma.n_i - r) U_i(x) psi(x+n_i)
##
## with the 24 unit neighbour vectors of `hcgeom` (8 axis, 16 diagonal), the
## link fields of `hcgauge` (FORMULATION 1.4) and QEX's DeGrand-Rossi gamma
## matrices (`physics/spinOld`, `gamma1..gamma4`; direction mu = 0..3 maps to
## gamma(mu+1)).  Correctness first (PLAN 1.2): every hop applies the full 4x4
## `(gamma.n - r)` spin matrix; no half-spinor projection.
##
## A honeycomb fermion is one Dirac fermion field per sublattice, both on the
## CELL layout:  `psiA(y)` at the A site `y`, `psiB(y)` at the B site `y+1/2`.
##
## Hop structure per cell `y` (`tests/tgeom.nim` checks the link endpoints;
## `tests/twilson.nim` checks the operator's locality and free spectrum):
##
##   A row, axis +mu :  uA[mu](y)          psiA(y+e_mu)     gamma factor (+gamma_mu - r)
##   A row, axis -mu :  uA[mu](y-e_mu)^dag psiA(y-e_mu)                   (-gamma_mu - r)
##   A row, diag d(delta): the A(y) -> B(y-dbar) link is uD[dbar](y-dbar)^dag,
##          dbar = delta xor 15:
##                      uD[dbar](y-dbar)^dag psiB(y-dbar)    (gamma.d(delta) - r)
##   B row, axis     :  same with uB
##   B row, diag d(delta): uD[delta](y)    psiA(y+delta)     (gamma.d(delta) - r)
##
## Machinery: QEX `Transporter`s for the 8 axis hops (fused U*shift, both
## signs), the task-L `HcShift16` binary trees for the 16-way fermion shifts
## (`psiA(y+delta)` forward, `psiB(y-delta)` backward), and 16 pre-shifted link
## copies `uDsh[delta](y) = uD[delta](y-delta)` rebuilt from the gauge field by
## `gaugeRefresh` (call it after any change to the links, e.g. `setBC`).
##
## D^dag is implemented directly: since the gammas are Hermitian and the
## reverse of hop `n_i` is hop `n_{opp(i)}` through the same link,
## `D^dag` is the SAME hopping sum with `(gamma.n_i - r)` replaced by
## `(-gamma.n_i - r)`, i.e. the spin factor of the opposite direction.
## gamma5-hermiticity `gamma5 D gamma5 = D^dag` is then a real test
## (tests/twilson.nim), not an identity built into the code.
##
## Plane-wave / momentum convention (used by tests and by tasks D3/D4)
## -------------------------------------------------------------------
## The free operator at cell momentum `k` is `hcfree.freeD8(k)`, an 8x8 matrix
## in (spin x sublattice), index `2*spin + sub`.  Its eigenvectors correspond
## to fields
##
##   psiA(y) = a_s e^{i k.y},   psiB(y) = b_s e^{i k.y}
##
## **both phased with the integer CELL coordinate `y`** -- the B sublattice
## carries NO extra half-site phase `e^{i k.(1/2,1/2,1/2,1/2)}`.  Antiperiodic
## time (`setBC`) shifts the allowed `k_3` by `pi/N_t` in this convention.

import base, layout, field, maths
import physics/qcdTypes
import rng
import hcgeom, hclayout, hcgauge

export hcgeom, hclayout, hcgauge

type
  HcFermion*[F] = object
    ## one Dirac fermion field per sublattice, on the cell layout
    a*, b*: F

  SpinMat* = typeof(gamma0)
    ## a constant 4x4 complex spin matrix (color-diagonal)

# ---------------------------------------------------------------------------
# fermion helpers
# ---------------------------------------------------------------------------

proc newDiracField*[V: static[int]](hl: HcLayout[V],
                                    nc: static[int] = getDefaultNc()): auto =
  ## one Dirac fermion field on the cell layout.  Built from the layout's own
  ## SIMD complex type (like `ColorMatrix(l, n)`), so it works for any `V`,
  ## including V = 1; for V = VLEN it is exactly `lo.DiracFermion()`'s type.
  type C = typeof(hl.lo.newDComplexV)
  type DF = Spin[VectorArray[4, Color[VectorArray[nc, C]]]]
  hl.lo.newField(DF)

proc newHcFermion*[V: static[int]](hl: HcLayout[V]): auto =
  ## a zeroed honeycomb fermion; **allocates**, call outside `threads:`
  type F = typeof(newDiracField(hl))
  var r: HcFermion[F]
  r.a = newDiracField(hl)
  r.b = newDiracField(hl)
  r.a := 0
  r.b := 0
  r

proc newOneOf*[F](x: HcFermion[F]): HcFermion[F] =
  ## same shape, zeroed; **allocates**, call outside `threads:`
  result.a = x.a.newOneOf
  result.b = x.b.newOneOf
  result.a := 0
  result.b := 0

# The following are allocation free: safe inside (or outside) `threads:`.

proc `:=`*[F](r: HcFermion[F], x: HcFermion[F]) =
  r.a := x.a
  r.b := x.b

proc `:=`*[F](r: HcFermion[F], v: SomeNumber) =
  r.a := v
  r.b := v

proc gaussian*(x: HcFermion, r: var RNGField) =
  x.a.gaussian r
  x.b.gaussian r

proc norm2*(x: HcFermion): float =
  ## |psiA|^2 + |psiB|^2 (global)
  norm2(x.a) + norm2(x.b)

proc dot*(x, y: HcFermion): auto =
  ## QEX field dot summed over both sublattices (global).
  ## QEX convention: the FIRST argument is conjugated (verified in twilson).
  dot(x.a, y.a) + dot(x.b, y.b)

proc norm2diff*(x, y: HcFermion): float =
  ## |x - y|^2 over both sublattices (global)
  norm2diff(x.a, y.a) + norm2diff(x.b, y.b)

proc applyGamma5*(r: HcFermion, x: HcFermion) =
  ## r = (gamma5 (x) 1_sublattice) x;  allocation free
  for e in r.a:
    r.a[e] := gamma5 * x.a[e]
  for e in r.b:
    r.b[e] := gamma5 * x.b[e]

# ---------------------------------------------------------------------------
# missing QEX overload
# ---------------------------------------------------------------------------

template mul*(r: var Spin, x: Color, y: Spin2) =
  ## Fused color-matrix times Dirac-fermion product, `r = x*y`.  QEX's
  ## `spinOld.mul` only handles all-`Spin` arguments; `Transporter`'s forward
  ## apply needs the mixed form (the color matrix acts as a scalar on the spin
  ## vector, `Sca2` in `maths/matrixConcept`).
  mixin mul
  mul(r[], x, y[])

# ---------------------------------------------------------------------------
# spin matrices
# ---------------------------------------------------------------------------

proc gammaDotDir*(dir: int): SpinMat =
  ## `gamma.n` for one of the 24 direction indices of `hcgeom`
  let n = toFloat(dirVec(dir))
  var t: SpinMat
  t := n[0]*gamma1 + n[1]*gamma2 + n[2]*gamma3 + n[3]*gamma4
  t

proc hopMat*(dir: int, rw: float): SpinMat =
  ## `(gamma.n_dir - rw)/6`, the hop factor of direction `dir`
  var t: SpinMat
  t := gammaDotDir(dir) - rw*gamma0
  var s: SpinMat
  s := (1.0/6.0)*t
  s

# ---------------------------------------------------------------------------
# the operator
# ---------------------------------------------------------------------------

type
  HcWilson*[V: static[int], MF, FF, TR, SHF, SHM] = ref object
    ## Persistent work object: gauge refs, pre-shifted diagonal links, all
    ## shifters and hop matrices.  Construct once with `newHcWilson`; call
    ## `gaugeRefresh` whenever the links change in place (e.g. after `setBC`).
    ## A `ref object` so that it can be captured by `threads:` blocks.
    g*: HcGauge[V, MF]
    uDsh*: array[nDiag, MF]      ## uDsh[d](y) = uD[d](y-d)
    tfA, tbA, tfB, tbB: array[nDim, TR]  ## axis transporters, forward/backward
    shAf: HcShift16[FF, SHF]     ## f[d](y) = psiA(y+d)
    shBb: HcShift16[FF, SHF]     ## f[d](y) = psiB(y-d)
    shDb: array[nDim, SHM]       ## single-axis backward link shifters
    gm: array[nDirs, SpinMat]    ## (gamma.n_dir - rw)/6
    rwCur: float
    rwValid: bool

proc setHopMats(w: HcWilson, rw: float) =
  ## rebuild the 24 hop matrices if the Wilson parameter changed;
  ## call outside `threads:`
  if w.rwValid and w.rwCur == rw: return
  for dir in 0..<nDirs:
    w.gm[dir] = hopMat(dir, rw)
  w.rwCur = rw
  w.rwValid = true

proc gaugeRefresh*(w: HcWilson) =
  ## Rebuild `uDsh[d](y) = uD[d](y-d)` from `w.g`.  Chained single-axis
  ## backward shifts along the set bits of `d` (sum_d popcount(d) = 32 shifts).
  ## Call **outside** `threads:` (opens its own).
  threads:
    w.uDsh[0] := w.g.uD[0]
    for d in 1..<nDiag:
      var cur = w.g.uD[d]
      for mu in 0..<nDim:
        if ((d shr mu) and 1) == 1:
          cur = w.shDb[mu] ^* cur
      w.uDsh[d] := cur

proc newHcWilson*[V: static[int], MF](g: HcGauge[V, MF]): auto =
  ## Build the persistent work object for the gauge field `g`.
  ## **Allocates** (shifter buffers, 16 shifted link fields, one prototype
  ## fermion); call outside `threads:`.
  var proto = newDiracField(g.hl)
  type FF = typeof(proto)
  type TR = typeof(newTransporter(g.uA[0], proto, 0, 1))
  type SHF = typeof(newShifter(proto, 0, 1))
  type SHM = typeof(newShifter(g.uA[0], 0, 1))
  var w = HcWilson[V, MF, FF, TR, SHF, SHM]()
  w.g = g
  for d in 0..<nDiag:
    w.uDsh[d] = g.uD[d].newOneOf
  for mu in 0..<nDim:
    w.tfA[mu] = newTransporter(g.uA[mu], proto, mu, 1)
    w.tbA[mu] = newTransporter(g.uA[mu], proto, mu, -1)
    w.tfB[mu] = newTransporter(g.uB[mu], proto, mu, 1)
    w.tbB[mu] = newTransporter(g.uB[mu], proto, mu, -1)
    w.shDb[mu] = newShifter(g.uA[0], mu, -1)
  w.shAf = newHcShift16(proto, 1)
  w.shBb = newHcShift16(proto, -1)
  w.gaugeRefresh
  w

proc applyDirac(w: HcWilson, r: HcFermion, x: HcFermion,
                m, rw: float, dag: bool) =
  ## r = D x (dag = false) or r = D^dag x (dag = true).
  ## Opens its own `threads:` block; r must not alias x.
  doAssert not (r.a == x.a or r.b == x.b), "hcwilson: r must not alias x"
  w.setHopMats(rw)
  # D^dag uses the spin factor of the opposite direction (see module docs)
  var cf: array[nDirs, int]
  for i in 0..<nDirs:
    cf[i] = if dag: opposite(i) else: i
  w.shAf.setSrc x.a
  w.shBb.setSrc x.b
  let mass = m + 4.0*rw
  threads:
    w.shAf.run
    w.shBb.run
    r.a := mass*x.a
    r.b := mass*x.b
    for mu in 0..<nDim:
      block:                                  # A row, axis +mu
        let gmm = w.gm[cf[axisIndex(mu, false)]]
        let h = w.tfA[mu] ^* x.a
        for e in r.a:
          r.a[e] += gmm * h[e]
      block:                                  # A row, axis -mu
        let gmm = w.gm[cf[axisIndex(mu, true)]]
        let h = w.tbA[mu] ^* x.a
        for e in r.a:
          r.a[e] += gmm * h[e]
      block:                                  # B row, axis +mu
        let gmm = w.gm[cf[axisIndex(mu, false)]]
        let h = w.tfB[mu] ^* x.b
        for e in r.b:
          r.b[e] += gmm * h[e]
      block:                                  # B row, axis -mu
        let gmm = w.gm[cf[axisIndex(mu, true)]]
        let h = w.tbB[mu] ^* x.b
        for e in r.b:
          r.b[e] += gmm * h[e]
    for d in 0..<nDiag:
      let gmm = w.gm[cf[diagIndex(d)]]        # gamma.d(delta) for both rows
      block:                                  # B row: uD[d](y) psiA(y+d)
        for e in r.b:
          r.b[e] += gmm * (w.g.uD[d][e] * w.shAf.f[d][e])
      block:               # A row: uD[dbar](y-dbar)^dag psiB(y-dbar)
        let dl = d xor 15
        for e in r.a:
          r.a[e] += gmm * (w.uDsh[dl][e].adj * w.shBb.f[dl][e])

proc D*(w: HcWilson, r: var HcFermion, x: HcFermion,
        m: float, rw: float = 1.0) =
  ## `r = (m + 4 rw) x + (1/6) sum_i (gamma.n_i - rw) U_i x(.+n_i)`.
  ## Opens its own `threads:` block.
  applyDirac(w, r, x, m, rw, false)

proc Ddag*(w: HcWilson, r: var HcFermion, x: HcFermion,
           m: float, rw: float = 1.0) =
  ## the Hermitian conjugate of `D` (direct implementation, not gamma5 D gamma5)
  applyDirac(w, r, x, m, rw, true)

# ---------------------------------------------------------------------------
# boundary conditions
# ---------------------------------------------------------------------------

proc setBC*(g: HcGauge) =
  ## Antiperiodic boundary conditions in time (mu = 3): flip the sign of every
  ## link crossing the time boundary, i.e. `uA[3]`, `uB[3]` and the 8
  ## `uD[delta]` with bit 3 of delta set, all at cell time `N_t - 1`.
  ## Mirrors QEX `setBC` (gauge/gaugeUtils); call inside a `threads:` block.
  ## Any `HcWilson` built on `g` needs `gaugeRefresh` afterwards.
  template flip(u: untyped) =
    let nt1 = u.l.physGeom[3] - 1
    for i in u.l.sites:
      if u.l.coords[3][i] == nt1:
        u{i} *= -1
  flip g.uA[3]
  flip g.uB[3]
  for d in 0..<nDiag:
    if ((d shr 3) and 1) == 1:
      flip g.uD[d]

when isMainModule:
  qexInit()
  let hl = newHcLayout([4, 4, 4, 4])
  var g = newHcGauge(hl)
  var w = newHcWilson(g)
  var x = newHcFermion(hl)
  var r = newHcFermion(hl)
  x.a := 1
  w.D(r, x, 0.1)
  echo "unit gauge, constant psiA: |D psi|^2 = ", r.norm2
  qexFinalize()
