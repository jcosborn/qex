# Module layout and interface contracts

**This file is normative.** Every work package codes against the signatures here. If you need to
change one, edit this file in the same commit and say so in [`06-status.md`](06-status.md) — do not
diverge silently, other packages are being written against it in parallel.

---

## 0. Global decisions (settled — do not relitigate)

| decision | choice | why |
|---|---|---|
| lattice framework | **plain Nim `seq`, no QEX `Layout`/`Field`** | the S² lattice is simplicial, not hypercubic; QEX Field buys nothing and costs the ET protocol |
| complex type | **`std/complex.Complex64`** | plain value semantics; bit-identical layout to LAPACK `dcomplex`, so `cast[ptr dcomplex]` is legal |
| spinor | **`array[2, Complex64]`** | 2-component, as in the paper |
| gauge field | **two `seq[float]`** (spatial links, temporal links) | non-compact; avoids any "inactive link" phantom dof |
| threading | **serial**, hot kernels isolated in single procs | target volumes are 10³–10⁴ sites; the prior attempt lost 60× to OpenMP barriers on this scale. LAPACK/BLAS supplies the parallelism where it matters. Revisit only with a measurement. |
| timers | call **`freezeTimers()`** in every app after `qexInit` | QEX `tic/toc` metadata allocation caused a 900× RSS blow-up on small volumes in the prior attempt |
| solver initial guess | **always zero** | MD reversibility |
| Zolotarev poles | stored in **σ² units** (`spec(X†X) ⊂ [σmin², σmax²]`) | the single easiest thing to get wrong |
| imports | `import base` (never `import qex`) plus explicit `rng/threefry4x64`, `eigens/lapack`, `eigens/linalgFuncs`, `algorithms/rk`, `algorithms/integrator`, `utils/resample`, `hmc/metropolis` | keeps the hypercubic machinery out |

Every app starts:
```nim
import base
...
qexInit()
freezeTimers()          # mandatory, see above
let p = readParams(...)
installStandardParams(); echoParams(); processHelpParam()
...
processSaveParams(); writeParamFile()
qexFinalize()
```

Build (from `<worktree>/build_mac`, which is already configured):
```bash
make run experimental/radial/tests/tgeom
```

---

## 1. Directory layout

```
src/experimental/radial/
  doc/            01-slides 02-formulation 03-targets 04-interfaces 05-plan 06-status
  core/           types.nim  geom.nim  lattice.nim  spinor.nim  analytic.nim
  ops/            wilson.nim  zolotarev.nim  solve.nim  overlap.nim  gaugeact.nim  flow.nim
  hmc/            pseudofermion.nim  trajectory.nim
  meas/           harmonics.nim  observables.nim  gevp.nim  fit.nim  dataio.nim
  tests/          tgeom.nim tspinor.nim twilson.nim tzolo.nim tsolve.nim toverlap.nim
                  tgauge.nim tflow.nim thmc.nim tmeas.nim tanalytic.nim tfit.nim
  campaign/       t2.sh (WP-L interacting driver)  free/ (WP-K figure scripts)
  rgeom.nim  rfree.nim  rspec.nim  rgauge.nim  rflow.nim  rhmc.nim  rmeas.nim  ranalyze.nim
  README.md
```
Library modules live in subdirectories; **apps live at the subtree root** (same convention as
`src/experimental/graph`). `make experimental/radial` then builds exactly the apps.

---

## 2. `core/types.nim` — shared value types

```nim
import std/complex
export complex

type
  Vec3* = array[3, float]
  Spinor* = array[2, Complex64]          ## 2-component Dirac spinor
  Spin* = seq[Spinor]                    ## a spinor field, flat index n = v + nv*t
  Mat2* = array[2, array[2, Complex64]]  ## 2x2 matrix (Pauli algebra, spin connection)

const
  sig1* = ...   ## Pauli matrices as Mat2
  sig2* = ...
  sig3* = ...
  id2*  = ...

proc dot*(a, b: Vec3): float
proc cross*(a, b: Vec3): Vec3
proc norm*(a: Vec3): float
proc unit*(a: Vec3): Vec3
proc `*`*(s: float, a: Vec3): Vec3
proc `+`*(a, b: Vec3): Vec3
proc `-`*(a, b: Vec3): Vec3

proc `*`*(m: Mat2, x: Spinor): Spinor
proc `*`*(m, n: Mat2): Mat2
proc adj*(m: Mat2): Mat2
proc expIsig3*(w: float): Mat2   ## exp(i*sigma3*w/2) -- NOTE the built-in factor 1/2
```
`expIsig3(w)` returns `diag(exp(i w/2), exp(-i w/2))`, i.e. exactly \(\Omega\) given \(\omega\).

---

## 3. `core/geom.nim` — the S² icosahedral lattice  *(WP-A)*

`geom.nim` re-exports `types.nim`, and `lattice.nim` re-exports `geom.nim`, so
`import core/lattice` is the single entry point for everything below.

```nim
const defaultTilt* = 0.7    ## chart tilt; maximises `poleGap` over L = 1..16

type
  Edge* = object
    a*, b*: int              ## endpoints; canonical orientation a -> b, a < b
    len*: float              ## geodesic length l_e
    dl*: array[2, float]     ## signed dual lengths l*, one per incident face (order matches `f`)
    dual*: float             ## dl[0] + dl[1]
    area*: float             ## diamond area A_e, EXACT spherical kite sum
                             ## kiteArea(len,dl[0])+kiteArea(len,dl[1]); the flat
                             ## form 0.5*len*dual is O(abar^2) apart (doc/06
                             ## "THE COUPLING CONVENTION" -- both published
                             ## Delta_0 values pin the exact form)
    ea*: array[2, float]     ## e^a_{ab}(a) = (cos alpha_a, sin alpha_a)
    eb*: array[2, float]     ## e^a_{ab}(b), the tangent at b continuing away from a,
                             ## so the reverse hop needs e^a_{ba}(b) = -eb
    omega*: float            ## spin-connection angle, principal value in (-pi, pi];
                             ## Omega_ab = sgn * expIsig3(omega), Omega_ba = its adjoint
    sgn*: float              ## +1.0 or -1.0, the Z2 spin-structure sign
    f*: array[2, int]        ## incident faces: f[0] traverses a->b, f[1] traverses b->a

  Face* = object
    v*: array[3, int]        ## vertices, counterclockwise seen from outside the sphere
    e*: array[3, int]        ## e[i] is the edge joining v[i] and v[(i+1) mod 3]
    s*: array[3, int]        ## +1 if e[i] runs v[i]->v[i+1], else -1
    area*: float             ## spherical area A_tri
    cc*: Vec3                ## dual point: circumcenter projected onto S^2
    sub*: array[3, float]    ## tilde A_i for v[0], v[1], v[2]; sum = area

  Sphere* = ref object
    lev*: int                ## refinement level L
    nv*, ne*, nf*: int       ## 10L^2+2, 30L^2, 20L^2
    pos*: seq[Vec3]          ## unit position vectors
    area*: seq[float]        ## A_y, dual polygon area per site
    nbr*: seq[seq[int]]      ## neighbour site indices (5 or 6), counterclockwise
    nbe*: seq[seq[int]]      ## matching edge indices, nbe[y][k] joins y and nbr[y][k]
    nbf*: seq[seq[int]]      ## matching faces, nbf[y][k] = (y, nbr[y][k], nbr[y][k+1])
    edges*: seq[Edge]
    faces*: seq[Face]
    abar*: float             ## mean edge length
    tilt*: float
    chart*: array[3, Vec3]   ## chart axes (x', y', z'); z' is the polar axis

proc newSphere*(lev: int, tilt = defaultTilt): Sphere
  ## Build the level-`lev` refined icosahedron projected onto the unit sphere.
  ## `tilt` selects the local-Lorentz gauge (see doc/02 section 2.3); observables must not depend on it.

proc omegaChart*(s: Sphere, e: int, nq = 32): float
  ## Independent oracle: -int cos(theta) dphi along the geodesic, by continuous tracking
  ## in the tilted polar chart. Must equal `s.edges[e].omega` modulo the Z2 sign.

proc holonomy*(s: Sphere, f: int): Mat2
  ## Product of Omega around face f in the orientation of `Face.s`.

proc checkClosure*(s: Sphere): tuple[offDiag, diag: float]
  ## Simplicial closure relation (IV.6): returns max |off-diagonal| and max |diag - A_tri|.

proc poleGap*(s: Sphere): tuple[site, link: float]
  ## Angular distance from the chart poles to the nearest site and to the nearest point of
  ## any link geodesic. Quality of the gauge choice; both must stay well clear of zero.

# spherical primitives, also used directly by the tests and by WP-G/WP-I
func sphArea*(a, b, c: Vec3): float        ## signed spherical excess, CCW positive
func circum*(a, b, c: Vec3): Vec3          ## dual point, sign fixed by a+b+c
func dualLen*(cc, p, q: Vec3): float       ## signed l* = arcsin(cc . unit(p x q))
func ptrans*(p, q, t: Vec3): Vec3          ## parallel transport of t along the geodesic p->q
func chartFrame*(tilt: float): array[3, Vec3]
func pframe*(ch: array[3, Vec3], p: Vec3): array[2, Vec3]   ## (e_theta, e_phi)
```

**Spherical vs flat dual identities.** The complex is intrinsic, so the flat relations
\(\sum_i\ell_i\ell^*_i/2=A_\triangle\) and \(A_e=\tfrac12\ell(\ell^*_1+\ell^*_2)\) hold only up to
\(O(\bar a_s^2)\). `Edge.area` uses the flat form **by definition** (it is what enters
\(\kappa\) and \(\beta_\ell\)); the exact spherical statement, which does hold to \(10^{-16}\), is
\[
A_\triangle=\sum_i 4\arctan\!\big(\tan(\ell_i/4)\,\tan(\ell^*_i/2)\big).
\]
See the WP-A entry in [`06-status.md`](06-status.md).

Construction notes (normative):
* Icosahedron vertices: the 12 cyclic permutations of \((0,\pm1,\pm\phi)/\sqrt{1+\phi^2}\), \(\phi=(1+\sqrt5)/2\).
* Refinement: for each of the 20 faces \((A,B,C)\), vertices \(v_{ijk}=(iA+jB+kC)/L\), \(i+j+k=L\),
  normalized; deduplicate shared vertices by rounding coordinates to \(10^{-9}\).
* Circumcenter: \(\hat c=\mathrm{unit}\big((b-a)\times(c-a)\big)\), sign fixed by \(\hat c\cdot(a{+}b{+}c)>0\).
* \(\ell^*\) is **signed**: with \(\hat n=\mathrm{unit}(y_i\times y_{i+1})\) (so that \(\hat n\) points to the
  side of the edge containing the third vertex), \(\ell^*=\arcsin(\hat c\cdot\hat n)\).
* Spin connection, **primary construction** (coordinate-free): deterministic vertex frame → SO(3)
  great-circle transport → SO(2) angle mod \(2\pi\) → Spin(2) lift with one \(\mathbb Z_2\) sign per
  edge → **solve the \(N_F\) face constraints over \(\mathbb F_2\)** for the signs (rank \(N_F-1\);
  one redundant constraint) → verify every face holonomy.
  `omegaChart` is the independent oracle; they must agree.

---

## 4. `core/lattice.nim` — S² × ℝ  *(WP-A)*

```nim
type
  Lat* = ref object
    sph*: Sphere
    nt*: int
    at*: float
    nsite*: int              ## sph.nv * nt
    kap*: seq[float]         ## kappa_e = 2 A_e / (abar * l_e), per edge (t-independent)
    kapT*: seq[float]        ## kappa'_y = A_y / (abar * at), per site
    volw*: seq[float]        ## A_y * at, the local volume element
    volbar*: float           ## mean of volw

proc newLat*(sph: Sphere, nt: int, at: float): Lat
proc sIdx*(l: Lat, v, t: int): int   ## site index v + nv*t, t taken mod nt
proc eIdx*(l: Lat, e, t: int): int   ## spatial-link index e + ne*t
proc tIdx*(l: Lat, v, t: int): int   ## temporal-link index v + nv*t
proc betaFace*(l: Lat, f: int, g2: float): float  ## (1/g2) * at / A_tri
proc betaEdge*(l: Lat, e: int, g2: float): float  ## (1/g2) * 2 A_e / (l_e^2 * at)
proc asOverAt*(l: Lat): float                     ## abar / at, must be >= 4/3
proc maxM*(l: Lat): float                         ## 0.9 * min(4/sqrt(3), sqrt(3)*abar/at)
```

Boundary conditions: gauge fields **periodic** in t; fermion temporal hops carry an explicit
\(-1\) across the \(t=n_t-1 \to 0\) seam.

---

## 5. `core/spinor.nim` — spinor-field linear algebra  *(WP-D)*

```nim
proc newSpin*(n: int): Spin                         ## n sites, zero initialized
proc zero*(x: var Spin)
proc `:=`*(x: var Spin, y: Spin)
proc axpy*(x: var Spin, a: float, y: Spin)          ## x += a*y
proc axpy*(x: var Spin, a: Complex64, y: Spin)
proc axpby*(x: var Spin, a: float, y: Spin, b: float)  ## x = a*y + b*x
proc scale*(x: var Spin, a: float)
proc dot*(x, y: Spin): Complex64                    ## sum conj(x)*y
proc redot*(x, y: Spin): float
proc norm2*(x: Spin): float
proc gaussian*(x: var Spin, r: var Threefry4x64)    ## <|x_i|^2> = 1 per complex component
proc pointSource*(n, site, comp: int): Spin
```
`spinor.nim` depends only on `core/types.nim` and `rng/threefry4x64` — it does **not** import
`lattice.nim`, so fields are sized by a plain `int` and `pointSource` takes a flat site index.
`lattice.nim` owns the convenience overloads `newSpin(l: Lat)` and `pointSource(l: Lat, v, t, c)`.
`axpby` is the CG search-direction update; it is exported because `ops/solve.nim` needs it.

## 6. `core/analytic.nim` — continuum formulas  *(WP-C, done)*

```nim
func jacobiP*(j: int, alpha, beta, z: float): float
  ## P_j^(alpha,beta)(z), (C.11), generalized real parameters, j >= 0. P_l(z) = jacobiP(l,0,0,z).
func fermionG*(t: float, nmax = -1): float
  ## G^(1,1)(t) = sign(t)/(16 pi sinh^2(|t|/2)); if nmax >= 0 use the truncated sum (V.3).
func fermionGPeriodic*(t, T: float, nmax = -1): float
  ## image-summed, antiperiodic in T
func gaugeG*(t: float, nmax = 200): float           ## (V.14)
func gaugeGPeriodic*(t, T: float, nmax = 200): float
func flatSpectrum*(kap, kapT: float, n1, n2, nt: int): seq[Complex64]   ## (IV.8)
func s2FermionProp*(theta: float, nmax: int): float ## (C.55) vs sigma1/(4 pi sin(theta/2)) (C.54)
func s2CurrentCorr*(z: float, nmax: int): float     ## (C.57)
func effDim*(g: openArray[float], at, T: float): seq[float]  ## (V.4)-(V.5)
```

Conventions (all settled in WP-C, see [`06-status.md`](06-status.md)):

* **Truncation.** `nmax >= 0` keeps exactly the modes up to `nmax` — that is what the \(n_{\max}\)
  fits need. `nmax < 0` means the untruncated correlator: the closed form for the fermion, and for
  the gauge tower (irrational exponents, no closed form) a sum carried until the terms underflow.
  `gaugeG`/`gaugeGPeriodic` keep the interface default `nmax = 200`.
* **Momentum grid.** `flatSpectrum` uses \(k=2\pi i/n\), \(i=0..n-1\), and emits both branches
  consecutively, so the result has `2*n1*n2*nt` entries and is closed under conjugation.
* **`s2FermionProp`** returns the scalar coefficient of \(\sigma_1\) (the paper plots the (1,2)
  spinor component and the matrix is \(\sigma_1\) times this scalar).
* **`effDim`** indexes `g[i]` at \(t=i\,a_t\) and returns `g.len-1` values, `result[i]` =
  \(\Delta_{\rm eff}(i\,a_t)\); \(G(T/2)\) is read at index `round(T/(2 at))`. Past \(T/2\) the sign
  flips, as it must.
* No `legendreP` was added: `jacobiP(l, 0, 0, x)` is the Legendre polynomial and nothing in WP-C
  needs the associated form. WP-I should add `legendreP` in `meas/harmonics.nim` where `ylm` lives.

## 7. `ops/wilson.nim` — Wilson-Dirac  *(WP-E, done — signatures below are the implemented ones)*

```nim
type
  Gauge* = object
    s*: seq[float]           ## spatial link angles, index eIdx(e,t), length ne*nt
    t*: seq[float]           ## temporal link angles, index tIdx(v,t), length nv*nt

  DwPart* = enum
    dwSpatial                ## the kappa sum over spatial neighbours of (IV.1)
    dwTemporal               ## the kappa' terms of (IV.1)
    dwNaive                  ## the e^a sigma_a / sigma_3 pieces: the antihermitian C of (IV.4)
    dwWilson                 ## the 1 and diagonal pieces: the hermitian B of (IV.4)
  DwParts* = set[DwPart]

const
  dwAll*   = {dwSpatial, dwTemporal, dwNaive, dwWilson}
  dwC*     = {dwSpatial, dwTemporal, dwNaive}    ## antihermitian part
  dwB*     = {dwSpatial, dwTemporal, dwWilson}   ## hermitian part
  dwSpace* = {dwSpatial, dwNaive, dwWilson}      ## the whole spatial operator
  dwTime*  = {dwTemporal, dwNaive, dwWilson}     ## the whole temporal operator

proc newGauge*(l: Lat): Gauge
proc zero*(u: var Gauge)
proc nlinkS*(l: Lat): int                        ## ne*nt
proc nlinkT*(l: Lat): int                        ## nv*nt

proc applyDw*(l: Lat, dst: var Spin, src: Spin, u: Gauge, m = 0.0, parts = dwAll)
proc applyDwAdj*(l: Lat, dst: var Spin, src: Spin, u: Gauge, m = 0.0, parts = dwAll)
  ## dst = (D_W - m) src   /   (D_W - m)^dag src, exactly Eq. (IV.1).
  ## `dst` must already have length l.nsite and must not alias `src`.  Allocation free.
proc applyDwDeriv*(l: Lat, dst: var Spin, src: Spin, u: Gauge, du: Gauge, parts = dwAll)
  ## dst = (delta D_W)[du] src
proc dwPullback*(l: Lat, f: var Gauge, left, right: Spin, u: Gauge,
                 scale = 1.0, add = false, parts = dwAll)
  ## f_link += scale * d[2 Re <left, D_W right>] / d theta_link
proc hatScale*(l: Lat): seq[float]               ## abar/sqrt(A_y), per sphere site
proc applyDwHat*(l: Lat, dst: var Spin, src: Spin, u: Gauge, work: var Spin,
                 m = 0.0, parts = dwAll)
proc applyDwHatAdj*(l: Lat, dst: var Spin, src: Spin, u: Gauge, work: var Spin,
                    m = 0.0, parts = dwAll)
  ## dst = (Dhat_W - m) src, Dhat_W = abar^2 A^{-1/2} D_W A^{-1/2}, A = diag(A_y).
  ## `work` is a caller-owned scratch field of length l.nsite; nothing is allocated.
proc denseDw*(l: Lat, u: Gauge, m = 0.0, parts = dwAll): seq[Complex64]
  ## column-major dense matrix, dimension 2*nsite, row index 2*sIdx + spin;
  ## for tests and small-lattice spectra only.  Assembled block by block from Mat2
  ## algebra, i.e. an independent second coding of (IV.1) — that is what makes
  ## T1.2d ("matrix-free == dense") a real test.
proc gaugeTransform*(l: Lat, u: var Gauge, alpha: openArray[float])
proc spinGaugeTransform*(l: Lat, x: var Spin, alpha: openArray[float])
```

**Link conventions (settled by WP-E, used by WP-F/G/H).**
`u.s[eIdx(e,t)]` is theta on the edge in its canonical orientation a → b, so the hop a → b
(row b, column a) carries \(U=e^{+i\theta}\) and the reverse hop its conjugate;
`u.t[tIdx(y,t)]` is theta on the temporal link (y,t) → (y,t+1). Both are periodic in t; the
fermion temporal hop carries the explicit −1 across the t = nt−1 → 0 seam, in **both**
directions of that link (which is what makes the C/B split exactly (anti)hermitian).
`gaugeTransform` implements doc/02 §5, \(\theta_e\to\theta_e+\alpha_b-\alpha_a\) and
\(\theta^t_{y,t}\to\theta^t+\alpha_{y,t+1}-\alpha_{y,t}\), so
\(D(\theta+d\alpha)=P D(\theta)P^\dagger\) with \(P={\rm diag}(e^{i\alpha})\) and the Ward
identity of §15 item 6 reads \(\delta_{d\alpha}D=i(\alpha D-D\alpha)\) — the **opposite sign**
to the one written there, which corresponds to the opposite link convention.
**Normalization used by the overlap** (see `overlap.nim`): the volume-normalized kernel is
\(\hat D_W = \bar a_s^2\,A^{-1/2}D_WA^{-1/2}\), \(A={\rm diag}(A_y)\) — a real diagonal
similarity, so the solver inner product stays plain, the adjoint is the same routine with
`applyDwAdj` inside, gauge covariance is manifest and a physical mass is site independent.
Its spectrum is \({\rm eig}(\hat D_W)=\bar a_s\,{\rm eig}(D_{\rm cont})\) (measured, see WP-E in
[`06-status.md`](06-status.md)). The raw `applyDw` is exactly (IV.1) and is untouched.
**Caution for WP-F:** the domain-wall height \(M\) of (IV.10) is quoted in **raw** \(D_{\rm lat}\)
units — \(M_0=\min(4/\sqrt3,\sqrt3\bar a_s/a_t)\) is the flat-limit \(\min(4\kappa,2\kappa')\) of
the raw operator — and \(\hat D_W\) differs from \(D_{\rm lat}\) by the site-dependent factor
\(\bar a_s^2/A_y\in[1.07,1.54]\) at L = 4 (\(\to 2/\sqrt3=1.1547\) in the flat limit). WP-E's
measurement of the T2.1 free column (below) confirms the paper's \(D_W\) is the raw one.

## 8. `ops/zolotarev.nim` — rational approximation  *(WP-B)*

```nim
proc ellipticK*(k: float): float                    ## AGM
proc jacobiSn*(u, k: float): float                  ## descending Landen / AGM

type Rat* = object
  order*, npole*: int
  smin*, smax*: float        ## sigma bounds (NOT squared)
  cst*: float                ## constant term
  pole*, res*, zero*: seq[float]   ## in sigma^2 units
  maxRelErr*, maxAbsErr*: float
  hash*: uint64

proc newRat*(smin, smax: float, order: int, nsample = 20001): Rat
  ## order must be odd and >= 3; npole = (order-1) div 2
proc ratValue*(r: Rat, x: float): float             ## cst + sum res_j/(x + pole_j) ~ 1/sqrt(x)
```
**Settled (WP-B).** `order = 31` → **15 shifts**, `order = 11` → **5 shifts**. The talk's
"6 poles (n=11)" counts partial-fraction *terms* (5 poles + the constant, which costs no solve);
its "15 poles (n=31)" counts poles. Multishift CG therefore sees 15 systems for the action
operator and 5 for the force operator. Measured `maxRelErr` at a 1/200 window: 4.6e-10 (order 31),
1.2e-3 (order 11).

## 9. `ops/solve.nim` — CG and multishift CG  *(WP-D)*

```nim
type
  CgInfo* = object
    iters*: int
    r2*: float               ## |b - A x|^2 / |b|^2, RELATIVE and always recomputed
    converged*: bool
  MultiCgInfo* = object
    iters*: int              ## multishift recurrence iterations only
    r2*: seq[float]          ## per shift, |b - (A+s_j) x_j|^2 / |b|^2, recomputed
    r2pre*: seq[float]       ## same, measured before refinement
    converged*: bool
    refined*: int            ## how many shifts needed a single-shift refinement pass
    refits*: int             ## CG iterations spent in those refinement passes

proc cgSolve*(x: var Spin, b: Spin, r2req: float, maxits: int,
              op: proc(dst: var Spin, src: Spin)): CgInfo
proc cgmSolve*(xs: var seq[Spin], b: Spin, shifts: openArray[float],
               r2req: float, maxits: int,
               op: proc(dst: var Spin, src: Spin)): MultiCgInfo
  ## op applies the UNSHIFTED positive operator. `shifts` must be ascending and > 0.
  ## After the recurrence converges, the TRUE residual is recomputed for every shift and any
  ## shift that misses tolerance is re-solved with a plain shifted CG. Always starts from x = 0.
```
`r2` is **relative** (divided by `|b|^2`) in both info records, and is always the recomputed
residual, never the recursive one — so the guard below reads literally. `cgSolve` therefore costs
one extra operator apply per solve.

Tolerance guard: accept when `r2 <= 1.001*r2req` (recursive and recomputed residuals differ at the
roundoff floor).

Two implementation points that callers depend on:
* **Converged shifts are frozen.** Once `z_j^2 r2 <= r2req |b|^2` shift `j` stops being updated.
  This is required, not an optimization: `z_j` decays like `prod 1/(1 + sg alpha_i)` and underflows
  to 0, after which the `z` recurrence is `0/0` and NaN-poisons `x_j`.
* **The refinement pass skips `j = 0`.** The recurrence *is* plain CG on `A + shifts[0]`, step for
  step, so re-solving the seed system reproduces `xs[0]` bit for bit.

## 10. `ops/overlap.nim` — overlap  *(WP-F, done — signatures below are the implemented ones)*

```nim
type
  SolveStats* = object
    nx*: int          ## D_W-level applies: X and X^dag count 1 each, H counts 2
    nmulti*: int      ## multishift solves
    miters*: int      ## multishift recurrence iterations
    mrefits*: int     ## multishift refinement iterations
    ncg*: int         ## plain CG solves (solveNormal outer, kernelWindow inner)
    cgiters*: int     ## their iterations
    ok*: bool         ## AND of every inner-solve `converged` since the last clear

  Ov* = ref object
    l*: Lat
    m*: float                  ## domain-wall height M (paper: rho), RAW units; default 1.0
    rat*: Rat                  ## frozen rational; poles in sigma^2 units
    work*: seq[Spin]           ## persistent scratch, 6 fixed slots -- see overlap.nim header
    xs*, xt*: seq[Spin]        ## multishift solution banks, npole fields each, preallocated
    r2inner*, r2outer*: float
    maxits*: int
    stats*: SolveStats

proc clearStats*(o: Ov)                                       ## stats = 0, ok = true
proc newOv*(l: Lat, m: float, rat: Rat, r2inner, r2outer: float, maxits: int): Ov
proc applyX*(o: Ov, dst: var Spin, src: Spin, u: Gauge)       ## X = D_W - M, RAW operator
proc applyXAdj*(o: Ov, dst: var Spin, src: Spin, u: Gauge)
proc applyH*(o: Ov, dst: var Spin, src: Spin, u: Gauge)       ## H = X^dag X; uses work[0]
proc applyOv*(o: Ov, dst: var Spin, src: Spin, u: Gauge, mass = 0.0)
proc applyOvAdj*(o: Ov, dst: var Spin, src: Spin, u: Gauge, mass = 0.0)
proc applyNormal*(o: Ov, dst: var Spin, src: Spin, u: Gauge, mass = 0.0)  ## D^dag D
proc solveNormal*(o: Ov, x: var Spin, b: Spin, u: Gauge, mass = 0.0): CgInfo
proc ovGradient*(o: Ov, f: var Gauge, left, right: Spin, u: Gauge,
                 scale = 1.0, add = false, mass = 0.0)
  ## THE single pullback: f (+)= scale * d[2 Re <left, D(mass) right>]/d theta.  Used by the
  ## HMC force, the Hasenbusch frames, the conserved current, and the Ward test.  There
  ## must not be a second, separately derived force kernel.
proc kernelWindow*(o: Ov, u: Gauge, iters = 32): tuple[smin, smax, lo, hi: float; inside: bool]
  ## inverse/power iteration; `inside` false => the frozen rational window is violated => STOP.
proc denseOv*(o: Ov, u: Gauge): seq[Complex64]
  ## exact dense D_ov = 1 + X (X^dag X)^{-1/2} via denseDw + zheev, tests only
```
Mass convention: **standard overlap at \(\rho=1\)**,
\(D(m)=(1-m/2)D_{\rm ov}+m\), with always-on validation \(0\le m<2\).
It is recorded as `massConvention=standard-overlap-rho1` and `overlapRho=1` in every
output manifest. Legacy additive mass \(\mu\) maps as
\(\mu=m/(1-m/2)\), or \(m=\mu/(1+\mu/2)\); equal numeric values are not equivalent.
The rational window is **frozen** at construction. Never rebuild it mid-ensemble.

**Implemented contract (WP-F; measured numbers in [`06-status.md`](06-status.md)):**
* Solve counts: `applyOv`/`applyOvAdj` = 1 multishift solve each; `applyNormal` and
  `ovGradient` = 2; `solveNormal` = `2*(iters+1)` (cgSolve recomputes the true residual).
  All counts and D_W-level applies accumulate in `stats`; `stats.ok` goes false if any
  inner solve missed its tolerance — check it, nothing raises in the apply path.
* `dst` must not alias `src` in any apply (same rule as `applyDw`).
* `newOv` preallocates `work`, `xs`, and `xt`. Solver callbacks capture gauge and mass
  within each call; their environments and `cgmSolve`/`cgSolve` scratch are temporary.
  Repeated calls retain no growing state. Each `Ov` requires exclusive use of its scratch.
* `kernelWindow`: `smin`/`smax` are Rayleigh quotients of \(X^\dagger X\) (σ units), so
  `smin` ≥ σ_min and `smax` ≤ σ_max unconditionally; `lo`/`hi` expand them by the eigenpair
  residual \(|Hv-\rho v|/|v|\), so [lo, hi] brackets the truth once the iteration is
  anywhere near converged. Both loops start from a fixed seeded vector, run at most
  `iters` steps and stop early below 1e-8 relative residual; inverse iteration costs one
  CG solve of \(Hx=v\) at `r2inner` per step. `inside` = `[lo,hi] ⊂ [rat.smin, rat.smax]`.
  It is a monitor and allocates its own vectors.
* `ovGradient` pullback, with \(R(H)=c_0+\sum_j r_jG_j\), \(G_j=(H+q_j)^{-1}\),
  \(z=R(H)\,\text{right}\), \(s_j=G_j\,\text{right}\), \(t_j=G_jX^\dagger\text{left}\):
  \(dF=2\,\mathrm{Re}[\langle\text{left},\delta X\,z\rangle-\sum_jr_j(\langle Xs_j,\delta X\,t_j\rangle+\langle Xt_j,\delta X\,s_j\rangle)]\),
  each term one `dwPullback` (\(\delta X=\delta D_W\), M constant); cost 2 multishift solves
  + \(2n_{\rm pole}{+}1\) X applies and pullbacks. The result is multiplied by
  \(1-m/2\), because \(\delta D(m)=(1-m/2)\delta D_{\rm ov}\).

**Kernel convention (settled 2026-08-21, after WP-E).** \(X=D_{\rm lat}-M\) uses the **raw**
operator of (IV.1) with **plain matrix adjoints**, \(M=1\) default — exactly the paper's (IV.9).
Rationale: the paper's entire construction is plain-matrix (the volume weights are *inside*
\(D_{\rm lat}\); the measure is flat; \(\det\mathcal D_{\rm lat}\) is a plain determinant; the GW
relation (IV.22) and the parity identity (IV.14) are stated with plain adjoints), and WP-E showed
the slide-8 legend \(\min|D_W-1|\) reproduces only with the raw operator, consistent with (IV.10)
being the raw flat doubler position. The hat kernel \(\hat D_W\) (a *congruence*, not a
similarity — it changes the spectrum whenever \(A_y\) varies) remains available for spectra and
diagnostics only. Do NOT build \(X\) from \(\hat D_W\).

## 11. `ops/gaugeact.nim` — non-compact U(1)  *(WP-G, done)*

Every signature below is implemented. `Gauge`, `newGauge`, `zero` and `gaugeTransform` are
imported from `ops/wilson.nim` (§7) and re-exported, so `import ops/gaugeact` is enough for
gauge-only code and there is exactly one gauge-field type in the tree.
Note `wilson.gaugeTransform(l, u, alpha)` is `u += gradient(alpha)`.

```nim
type
  GeomConv* = enum           ## which O(abar^2) transcription of (IV.26) to use
    gcGeodesic               ## geodesic l, spherical A_tri, A_e = l(l*_1+l*_2)/2
    gcExactArea              ## as gcGeodesic but A_e = exact spherical diamond area
    gcFlat                   ## chord l, planar A_tri, in-plane l*
  Beta* = object             ## precomputed plaquette couplings (IV.26)
    face*, edge*: seq[float]
    afac*: seq[float]        ## the A_tri that went into beta_tri (J^t normalization)
    g2*: float
    conv*: GeomConv

proc newBeta*(l: Lat, g2: float, conv = gcExactArea): Beta   ## default = the paper's convention
func nlink*(l: Lat): int                 ## (ne + nv)*nt
func slink*(l: Lat): int                 ## ne*nt, where the temporal links start

# vector space (the CG needs it; same names as core/spinor.nim's Spin ops)
proc `:=`*(x: var Gauge, y: Gauge)
proc axpy*(x: var Gauge, a: float, y: Gauge)
proc axpby*(x: var Gauge, a: float, y: Gauge, b: float)
proc scale*(x: var Gauge, a: float)
func dot*(x, y: Gauge): float
func norm2*(x: Gauge): float
func toSeq*(u: Gauge): seq[float]        ## flat link vector, spatial then temporal
proc fromSeq*(l: Lat, v: openArray[float]): Gauge
proc unitSource*(l: Lat, b: var Gauge, link: int)

func plaqSpatial*(l: Lat, u: Gauge, f, t: int): float            ## Theta_face
func plaqTemporal*(l: Lat, u: Gauge, e, t: int): float
func jtop*(l: Lat, u: Gauge, f, t: int): float                   ## Theta_face/A_tri  (V.12)
func jtop*(l: Lat, u: Gauge, b: Beta, f, t: int): float          ## with b's area convention

proc gaugeActionParts*(l: Lat, u: Gauge, b: Beta): tuple[sp, tp: float]
proc gaugeAction*(l: Lat, u: Gauge, b: Beta): float
proc gaugeAction*(l: Lat, u: Gauge, g2: float): float
proc gaugeForce*(l: Lat, f: var Gauge, u: Gauge, b: Beta)        ## f = dS/dtheta = M u
proc gaugeForce*(l: Lat, f: var Gauge, u: Gauge, g2: float)
template applyM*(l: Lat, dst: var Gauge, src: Gauge, b: Beta)    ## = gaugeForce
proc mDiagonal*(l: Lat, d: var Gauge, b: Beta)                   ## diag(M), sizes the flow step
proc mDiagMax*(l: Lat, b: Beta): float

proc gradient*(l: Lat, p: var Gauge, alpha: openArray[float])    ## d alpha
proc divergence*(l: Lat, d: var seq[float], p: Gauge)            ## d^dagger p
proc laplace*(l: Lat, dst: var seq[float], src: openArray[float])  ## d^dagger d
proc projectGauge*(l: Lat, p: var Gauge, r2req = 1e-24, maxits = 10000): CgInfo
  ## remove the gauge-orbit component; applied to the committed field, the refreshed
  ## momentum, AND every MD force at every level.
proc projectFlat*(l: Lat, p: var Gauge)                          ## uniform Polyakov mode
proc projectKernel*(l: Lat, p: var Gauge, r2req = 1e-24, maxits = 10000): CgInfo

proc cgM*(l: Lat, x: var Gauge, b: Gauge, bt: Beta, r2req = 1e-20, maxits = 20000): CgInfo
proc pseudoSolve*(l: Lat, x: var Gauge, b: Gauge, bt: Beta,
                  r2req = 1e-20, maxits = 20000): tuple[proj, sol: CgInfo]
  ## the literal (V.16)-(V.17) double CG; use for sources that are NOT transverse
type RegOp* = object                     ## A = M + sig d d^dagger + tau P P^T/|P|^2
  bt*: Beta
  sig*, tau*: float
proc newRegOp*(l: Lat, bt: Beta, sig = 0.0, tau = 0.0): RegOp
proc applyReg*(l: Lat, o: var RegOp, dst: var Gauge, src: Gauge)
proc regSolve*(l: Lat, x: var Gauge, b: Gauge, o: var RegOp,
               r2req = 1e-24, maxits = 100000): CgInfo
  ## SPD, so stable to the roundoff floor; equals Mtilde^{-1} b for transverse b

proc heatbathSource*(l: Lat, b: var Gauge, bt: Beta, r: var Threefry4x64)
proc heatbath*(l: Lat, u: var Gauge, bt: Beta, r: var Threefry4x64,
               r2req = 1e-20, maxits = 20000): CgInfo
proc heatbath*(l: Lat, u: var Gauge, g2: float, r: var Threefry4x64)
  ## exact Gaussian sampling of the transverse modes (the action is quadratic)
proc gaugePropagator*(l: Lat, g2: float, srcLink: int,
                      r2req = 1e-20, maxits = 20000): seq[float]
  ## <theta_m theta_n> via the double-CG pseudo-inverse (V.16)-(V.17)
proc triSource*(l: Lat, b: var Gauge, f, t: int)
  ## incidence row of the spatial plaquette, so <Theta Theta> = b^T Mtilde^{-1} b'
```

**`gcExactArea` is the paper's convention for \(A_\ell\).** With it, Δ₀(L=1, L_t=120, T=16)
comes out 1.332430 against the published 1.33242 — every digit — where the flat-form
`Edge.area` gives 1.356697. See the WP-G entry in [`06-status.md`](06-status.md) T1.5b; the
same question is open for the fermion \(\kappa\).

**ker M is one dimension bigger than the gauge orbit.** The uniform temporal mode
\(\theta^t_v(t)=c\), \(\theta^s=0\) costs no action and is *not* a gauge mode (a gauge function
periodic in \(t\) can only make temporal shifts that sum to zero round the circle). It is
orthogonal to the gauge orbit, so `projectKernel` = `projectGauge` + `projectFlat`, and
\(\dim\ker M = n_V L_t\), \({\rm rank}\,M = n_E L_t\).

## 12. `ops/flow.nim` — gradient flow  *(WP-G, done)*
Built on **`algorithms/rk.nim`** (generic, closure-driven, zero QEX coupling).
```nim
proc newFlowOp*(l: Lat, u: Gauge, b: Beta): RK2NOp[Gauge, Gauge]
proc flowStep*[S: static[int]](l: Lat, u: var Gauge, b: Beta, eps: float,
                               coeffs: RK2NCoeffs[S])
proc flowStep*[S: static[int]](l: Lat, u: var Gauge, g2, eps: float,
                               coeffs: RK2NCoeffs[S])
proc flowRun*[S: static[int]](l: Lat, u: var Gauge, b: Beta, times: openArray[float],
                              step: float, coeffs: RK2NCoeffs[S],
                              measure: proc(t: float, u: Gauge))
proc flowRun*(l: Lat, u: var Gauge, g2: float, times: openArray[float], step: float,
              measure: proc(t: float, u: Gauge))                 ## RK4CK_2N
proc energyDensity*(l: Lat, u: Gauge, b: Beta): float            ## E_s = S_spatial/(4 pi T)
proc energyDensity*(l: Lat, u: Gauge, g2: float): float
proc energyDensityT*(l: Lat, u: Gauge, b: Beta): float           ## the electric partner
proc energyDensityT*(l: Lat, u: Gauge, g2: float): float
```
`coeffs: auto` became `coeffs: RK2NCoeffs[S]` (a static-int generic): `auto` cannot carry a
default and rk.nim's coefficient sets are distinct types. `measure` may be `nil`.
For the Gaussian action the flow is **linear** (\(\dot\theta=-M\theta\)) — the RK result is
checked against the exact matrix exponential on a small lattice.
**Normalization**: the flow is generated by `S(theta; b)`, so passing `newBeta(l, g2)` makes the
flow rate carry the \(1/g^2\) of (IV.26); `newBeta(l, 1.0)` is the standard (Lüscher) convention
in which the flow time is a length\(^2\). For the free theory the two differ only by \(s\to s/g^2\).
`energyDensity` is the **action density**, \(S_{\rm spatial}/(4\pi T)\), which is
\(\tfrac14\langle F_cF_c\rangle\) for the canonically normalized field and is \(g^2\)-independent
in the free theory.

## 13. `hmc/pseudofermion.nim`, `hmc/trajectory.nim`  *(WP-H, done — signatures below are the implemented ones)*
Structurally like `src/experimental/graph/hmcgauge/trajectory.nim`, built on
**`mdevolve`** (`newIntegratorPair` with a combined updater, one `mkOmelyan2MN` per level,
`newParallelEvolution` sharing the position update) and **`hmc/metropolis.nim`**
(`RadialHmc` subclasses `MetropolisRoot`; hooks `start/getH/generate/globalRand/accept/reject`).

```nim
# pseudofermion.nim -- the Hasenbusch ladder, standard mass D(m) = (1-m/2)D_ov + m,
# Q_i = D(m_i)^dag D(m_i); one pseudofermion pair per two flavors.
#   ratio frame i < K:  S_i = phi^dag D_{i+1} Q_i^{-1} D_{i+1}^dag phi,
#                       heatbath phi = D_{i+1} Q_{i+1}^{-1} D_i^dag xi  => S_i = |xi|^2
#   heaviest frame K:   S_K = phi^dag Q_K^{-1} phi, heatbath phi = D_K^dag xi
# Forces via ovGradient ONLY:
#   heaviest: dS = -2Re<y, dD_K eta>, eta = Q_K^{-1} phi, y = D_K eta
#   ratio:    dS = +2Re<phi, dD_{i+1} eta> - 2Re<y, dD_i eta>,
#             chi = D_{i+1}^dag phi, eta = Q_i^{-1} chi, y = D_i eta
type Pf* = ref object
  l*: Lat
  actOp*, frcOp*: Ov         ## SAME lattice and M, DIFFERENT rational order (31 vs 11)
  nf*: int                   ## even, >= 2; copies = nf div 2
  masses*: seq[float]        ## strictly increasing Hasenbusch ladder, masses[0] = physical
  phi*: seq[seq[Spin]]       ## [copy][frame]
  xi2*: seq[seq[float]]      ## |xi|^2 recorded by the last refresh (heatbath-identity test)

func ncopy*(p: Pf): int
func nframe*(p: Pf): int
proc newPf*(l: Lat, actOp, frcOp: Ov, nf: int, masses: seq[float]): Pf
proc refreshFrame*(p: Pf, u: Gauge, c, i: int, r: var Threefry4x64)
proc refresh*(p: Pf, u: Gauge, r: var Threefry4x64)   ## all frames from one stream
proc frameAction*(p: Pf, o: Ov, u: Gauge, c, i: int): float
proc pfAction*(p: Pf, u: Gauge): float                ## actOp, all copies and frames
proc frameForce*(p: Pf, o: Ov, f: var Gauge, u: Gauge, c, i: int, add = false)
proc pfForce*(p: Pf, f: var Gauge, u: Gauge, level: int)  ## frcOp, summed over copies

# trajectory.nim -- MD driver, projections, checkpointing
const rkMomentum* = 1; rkAccept* = 2; rkPseudo* = 3    ## RNG purpose keys
func mixKey*(seed: uint64, traj, purpose: int, copy = 0, frame = 0): uint64
proc keyedRng*(seed: uint64, traj, purpose: int, copy = 0, frame = 0): Threefry4x64
proc gaussian*(x: var Gauge, r: var Threefry4x64)      ## one unit normal per link

type RadialHmc* = ref object of MetropolisRoot
  l*: Lat
  bt*: Beta
  pf*: Pf                    ## nil = pure gauge
  u*, p*: Gauge              ## committed field, momentum
  tau*: float
  seed*: uint64
  traj*: int                 ## trajectory counter (incremented by `start`)
  fcount*: seq[int]          ## force evaluations per level
  projR2*: float
  forceAccept*, forceReject*: bool

proc newRadialHmc*(l: Lat, bt: Beta, pf: Pf, tau: float, steps: seq[int],
                   seed: uint64, projR2 = 1e-24): RadialHmc
  ## steps per level: [gauge, heaviest frame, ..., lightest frame]
proc setSteps*(m: RadialHmc, steps: openArray[int])
proc clearForceCounts*(m: RadialHmc)
proc mdEvolve*(m: RadialHmc)                           ## one tau of MD (evolve + finish)
proc hmcH*(m: RadialHmc): float                        ## |p|^2/2 + S_g + S_pf(actOp)
proc refreshMomentum*(m: RadialHmc, traj: int)
proc refreshPseudo*(m: RadialHmc, traj: int)           ## keyed per (copy, frame)
# metropolis hooks: start getH generate globalRand accept reject (exported)
proc transversality*(l: Lat, p: Gauge): tuple[divP, flatP: float]
type RevInfo* = object
  dh*, duRms*, duMax*, dpRms*, dpMax*, divP*, flatP*: float
proc reversibilityCheck*(m: RadialHmc): RevInfo       ## fwd, p -> -p, fwd; restores state
proc windowCheck*(m: RadialHmc): tuple[smin, smax, lo, hi: float, inside: bool]
proc saveCheckpoint*(m: RadialHmc, path: string)
proc loadCheckpoint*(m: RadialHmc, path: string)      ## raises ValueError on any mismatch
```

Force levels: level 0 = gauge (innermost, most steps), levels 1..nframe = Hasenbusch frames
heaviest→lightest, one `mkOmelyan2MN` each over the shared position update.
Never use the low-order rational in the Hamiltonian or the heatbath; the actOp/frcOp split
makes the accept/reject `dH` saturate at the dt-independent action mismatch
\(\sim 2\,\mathrm{maxRelErr}(11)\,S_{\rm pf}\times O(2\%)\) — measured, see WP-H in
[`06-status.md`](06-status.md).

**Gauge zero modes (non-negotiable).** `projectKernel` (WP-G) is applied to (a) the committed
field — at construction/load and at every commit, so `reject` stays a bitwise restore of the
trajectory's start field; (b) the refreshed momentum; (c) every MD force at every level.
All three are required: at fixed \(\phi\) the extended action has a longitudinal force even
though the integrated determinant is gauge invariant.  Momentum dof
\(=(n_E{+}n_V)L_t-\dim\ker M=n_EL_t\).

**Randomness is trajectory addressed**: every draw comes from a Threefry stream seeded by a
splitmix64 mix of (baseSeed, trajectory number, purpose, copy, frame); no generator state is
serialized and a restart depends only on the committed trajectory counter (restart draws are
bitwise identical — tested).

**Checkpoint** (versioned binary, trailing FNV-1a over the payload): magic/version, lev, nt,
at, g2, geometry convention, nf, M, both rational hashes (`Rat.hash`), masses, tau, per-level
steps, seed, trajectory counter, the gauge field. Version 2 also stores the explicit
`standard-overlap-rho1` convention identifier. `loadCheckpoint` validates every manifest
field against the live configuration and refuses mismatches; legacy version-1 additive-mass
checkpoints are deliberately rejected.

**Solver targets.** Keep `r2inner`/`r2outer` ~2 decades above the operators' roundoff floors
(WP-D): the WP-F suite values (1e-26, 1e-22) sit ON the floor once `solveNormal` outer solves
and MD-moved fields are in play — `thmc` measured multishift floors ~8e-27 and outer floors
~5e-23 at L=1, n_t=6, and uses (1e-24, 1e-20); `rhmc` defaults to (1e-22, 1e-18).

## 14. `meas/*`  *(WP-I, WP-J)*
```nim
# harmonics.nim  -- IMPLEMENTED (WP-I).  REAL harmonics (doc/07 1.3), so the
# earlier Complex64 sketch became float throughout.
func legendreP*(l: int, x: float): float    ## recurrence; the ylm oracle (addition theorem)
func ylm*(l, m: int, v: Vec3): float        ## real orthonormal Y_lm, l <= 4, hard-coded
func siteProject*(sph: Sphere, x: openArray[float], l, m: int): float
  ## sum_y A_y Y_lm(pos_y) x_y   (x per site, one time slice)
func faceProject*(sph: Sphere, x: openArray[float], l, m: int): float
  ## sum_f Y_lm(cc_f) x_f -- NO area weight: for x = Theta_f the area cancels (doc/07 4.1)
func siteGram*(sph: Sphere, l1, l2: int): seq[float]   ## (2l1+1)x(2l2+1) row-major
func faceGram*(sph: Sphere, l1, l2: int): seq[float]   ## same, sum_f A_tri Y Y'
type IcosaGroup* = object
  rot*: seq[array[3, Vec3]]          ## the 60 rotations of I (rows)
  site*, face*: seq[seq[int]]        ## induced site/face permutations
  siteDev*, faceDev*: float          ## measured equivariance residuals (~1e-16)
proc icosaGroup*(sph: Sphere, tol = 1e-6): IcosaGroup

# observables.nim  -- IMPLEMENTED (WP-I)
func fiveFoldSite*(sph: Sphere): int            ## first original icosahedron vertex
proc propSolve*(o: Ov, x: var Spin, b: Spin, u: Gauge, mass = 0.0): CgInfo
  ## x = S b = (D^dag D)^{-1} D^dag b  (adjoint FIRST; D = D(mass))
proc propSolveDag*(o: Ov, x: var Spin, b: Spin, u: Gauge, mass = 0.0): CgInfo
  ## x = S^dag b = D (D^dag D)^{-1} b
proc propagatorT*(o: Ov, u: Gauge, mass: float, src: int): seq[Spinor]
  ## G(t) at the coincident site: source (src, t=0, component 0)
proc condensatePS*(o: Ov, u: Gauge, mass: float, nnoise: int,
                   r: var Threefry4x64): tuple[v, e: float]
  ## Re tr[(1 - D_ov/2)D(m)^{-1}]/(2 nsite), Gaussian noise, stderr
proc condensateDense*(o: Ov, u: Gauge, mass: float): float   ## exact spectral oracle
proc denseS*(o: Ov, u: Gauge, mass = 0.0): seq[Complex64]    ## exact dense S (tests)
proc denseOvDeriv*(o: Ov, u: Gauge, du: Gauge, mass = 0.0): seq[Complex64]
  ## exact dense delta D(mass)[du] of the RATIONAL operator (tests); pinned
  ## against ovGradient -- the ONE current kernel (doc/07 1.1)
proc tsliceForm*(l: Lat, w: openArray[float], t: int): Gauge ## du on slice-t temporal links
proc linkCurrent*(o: Ov, u: Gauge, left, right: Spin, mass = 0.0): tuple[re, im: Gauge]
  ## <left, K_l right> per link via TWO ovGradient calls (left and i*left)
proc tsliceAmp*(l: Lat, g: tuple[re, im: Gauge], w: openArray[float]): seq[Complex64]
  ## A(t) = sum_y w_y <left, K_{(y,t)} right>; w = Y_lm site values (the
  ## dual area sits inside J_link, so l = 0 means w = 1: the total charge)
proc wardChargeScan*(o: Ov, u: Gauge, mass: float, v0, ta, tb: int):
    tuple[c: seq[Complex64], jump: Complex64]
  ## total charge inserted in a fermion line: exactly two plateaus with
  ## jumps -+jump = -+i S_ba at ta, tb (any mass)
type CurrentSample* = object
  a*, b*, d*: seq[Complex64]         ## flattened [iop*nt + t]
proc currentSample*(o: Ov, u: Gauge, mass: float, w: openArray[seq[float]],
                    r: var Threefry4x64): CurrentSample
  ## one (eta, xi) noise pair; E[a1 b2] = tr[K1 S K2 S], E[d] = tr[K S];
  ## 2 sequential solves + 3 linkCurrent, no tangent kernel anywhere
proc currentCorrConn*(samples: openArray[CurrentSample], k1, k2: int):
    tuple[v, e: Complex64]           ## mean/stderr of tr[K1 S K2 S]
proc currentTraceDisc*(samples: openArray[CurrentSample], k1, k2: int):
    tuple[v, e: float]               ## 2Re tr[K1 S] 2Re tr[K2 S], cross-noise products
proc scalarCorrDense*(o: Ov, u: Gauge, mass = 0.0): tuple[ps, fs: seq[float]]
proc scalarCorrPoint*(o: Ov, u: Gauge, mass: float, v0, t0: int): tuple[ps, fs: seq[float]]
  ## sigma_PS / sigma_FS connected timeslice correlators (doc/07 section 3);
  ## F=(1-D_ov^dag)S^dag=((1+m/2)S^dag-1)/(1-m/2);
  ## at mass 0 PS and FS are IDENTICAL at every dt, configuration by configuration
proc jtopProject*(l: Lat, u: Gauge, lh, mh: int): seq[float] ## sum_f Y Theta_f per t
proc f2Project*(l: Lat, u: Gauge, lh, mh: int): seq[float]
  ## sum_f Y Theta^2/A_f + sum_e Y(mid) Theta_e^2 2A_e/(l_e a_t)^2, RAW
  ## (vacuum-subtract the l = 0 channel at analysis time)
type LoopShape* = enum               ## GEVP basis; 5-7 are the I_h-covariant,
  lsTri, lsRhomb, lsStar, lsQuad,    ## time-reflection-EVEN realizations --
  lsTPlaq, lsTRect2, lsTRhomb        ## see doc/06 WP-I (change vs doc/07 4.2)
func loopCount*(l: Lat, sh: LoopShape): int
func loopCenter*(sph: Sphere, sh: LoopShape, i: int): Vec3
func loopFlux*(l: Lat, u: Gauge, sh: LoopShape, i, t: int): float
proc loopProject*(l: Lat, u: Gauge, sh: LoopShape, lh, mh: int): seq[float]
proc loopOps*(l: Lat, u: Gauge, shapes: openArray[LoopShape], lh, mh: int): seq[seq[float]]
proc loopSource*(l: Lat, b: var Gauge, sh: LoopShape, lh, mh, t: int)
  ## incidence vector of the projected shape operator; exactly transverse
proc jtopCorrExact*(l: Lat, bt: Beta, lh: int, r2req = 1e-26,
                    maxits = 200000): seq[seq[float]]
  ## exact free-theory <O_lm(dt) O_lm'(0)> matrices via regSolve; result[dt]
  ## is (2l+1)^2 row-major -- the machine-precision degeneracy object
proc loopCorrExact*(l: Lat, bt: Beta, shapes: openArray[LoopShape], lh, mh: int,
                    r2req = 1e-26, maxits = 200000): seq[seq[seq[float]]]

# gevp.nim  -- IMPLEMENTED (WP-I)
proc gevp*(c: seq[seq[seq[float]]], t0, t: int, cut = 1e-10): seq[float]
  ## lambda_n(t, t0) DESCENDING on the subspace where C(t0) > cut*evmax;
  ## full rank -> the committed eigens/linalgFuncs.zeigsgv, rank deficient ->
  ## truncate + whiten + zeigs.  Rank deficiency is the NORMAL case (fewer
  ## states than operators; the L=1 l=1 gauge channel is exactly rank 1).
proc gevpDims*(c: seq[seq[seq[float]]], t0: int, at: float, cut = 1e-10): seq[seq[float]]
  ## Delta_n(t) = -ln(lambda_n(t+at)/lambda_n(t))/at, result[t][n]
proc gevpCheck*(c: seq[seq[seq[float]]], t0: int): tuple[evmin, evmax, cond, asym: float]
  ## spectrum of the symmetrized C(t0) and the relative antisymmetry of the raw input

# fit.nim   -- IMPLEMENTED (WP-J); the signatures below supersede the earlier sketch
type
  FitStatus* = enum
    fitLimit, fitOk, fitShort, fitSingular, fitStalled
  PlateauFit* = object
    d0*, c*, dp*: float                ## Delta_0, c, Delta' of (V.6)
    ed0*, ec*, edp*: float
    cov*: array[3, array[3, float]]    ## order (d0, c, dp)
    chi2*: float
    dof*: int                          ## npoint - 3
    iters*: int
    status*: FitStatus                 ## only fitOk is usable for extrapolation
  LineFit* = object
    a*, b*, ea*, eb*, cab*: float      ## intercept, slope, errors, cov(a,b)
    chi2*: float
    dof*: int                          ## npoint - 2
  PlaneFit* = object
    a*, cs*, ct*: float                ## Delta_0^cont, c_s, c_t of (V.7)
    ea*, ecs*, ect*: float
    cov*: array[3, array[3, float]]    ## order (a, cs, ct)
    chi2*: float
    dof*: int                          ## npoint - 3
  NmaxFit* = object
    nmax*: int
    c*: float                          ## fitted overall normalization
    res*: float
    dof*: int                          ## npoint - 2 (C and nmax)
    resDof*: float
  SeriesStat* = object
    mean*, err*, bias*: float
    hasErr*: bool                      ## false when fewer than two groups are available
    tau2*, tau2p*: float               ## 2 tau_int: Wolff window / positive sequence
    neff*: float                       ## n / max(1, tau2)
    n*, stride*, blockSize*, nblock*: int
    provisional*: bool                 ## nblock < 4
  Estimate* = object
    v*: float
    ok*: bool
  Jk* = object
    full*: Estimate
    reps*: seq[Estimate]               ## failed replicas retain their group positions
    trajs*: seq[int]                   ## ordered identities within one ensemble
    bs*: int

func chi2dof*(f: PlateauFit|LineFit|PlaneFit): float

proc effMass*(c: openArray[float], at, T: float): seq[float]
  ## (V.4)-(V.5).  c sampled at t = i*at; result[i] is Delta_eff at the same t, len = c.len-1.
proc plateauFit*(m: openArray[float], t0, t1: int, at = 1.0,
                 e: openArray[float] = [], maxit = 200): PlateauFit
  ## (V.6) over the index window [t0, t1), t = at*i.  `e` is indexed like `m`.
  ## Empty `e` = unit weights, and then the covariance is scaled by chi2/dof.
proc contFit*(x, y, e: openArray[float]): LineFit                  ## weighted y = a + b x
proc contFit2*(as2, at2, y, e: openArray[float]): PlaneFit
  ## (V.7) y = a + c_s*as2 + c_t*at2.  `as2`/`at2` are ALREADY SQUARED.
  ## Raises ValueError if the points do not span a plane (e.g. one a_t only).
proc nmaxFit*(g: openArray[float], model: proc(nmax: int): seq[float],
              nmaxRange: Slice[int], e: openArray[float] = []): NmaxFit
  ## (V.9).  Empty `e` selects the relative residual w_t = 1/g_t^2.
proc jack*(s: openArray[float], stride = 1, bs = 0): SeriesStat    ## uses utils/resample
  ## blockSize = ceil(max(1, 2 tau_int)) unless `bs > 0` overrides it.
proc effLocal*(c: openArray[float], at: float, positive = false): seq[Estimate]
  ## Log ratios require adjacent values with the same nonzero sign; positive=true for GEVP.
proc deltaFit*(m: openArray[Estimate], t0, t1: int, at: float): Estimate
proc deltaFit*(c: openArray[float], t0, t1: int, at: float): Estimate
  ## Fit the longest available run in [t0,t1), retaining the first on a tie.
proc jkFrom*[V: float|Estimate](cs: seq[seq[float]], est: proc(c: seq[float]): V,
                              trajs: seq[int], bs = 1): Jk
proc jkStat*(j: Jk): tuple[v, e: Estimate]
proc jkMean*(js: openArray[Jk]): Jk
proc ratioVE*(a, b: Jk): tuple[v, e: Estimate]
  ## Paired groups preserve covariance; overlapping unpaired groups have no error estimate.
  ## Disjoint identities use independent-error propagation. All IDs share one ensemble domain.
  ## A failed replica makes uncertainty unavailable. Availability is serialized as text.

# dataio.nim  -- IMPLEMENTED (WP-J).  TSV, `#key=value` header lines, %.17g numbers.
proc writeTsv*[T: float|string = float](path: string, header: openArray[(string, string)],
                                      colNames: openArray[string], cols: openArray[seq[T]])
  ## Atomically replace a complete file. String cells are literal text without tabs or newlines.
proc requireTsvMeta*(path: string, expected: openArray[(string, string)])
proc outputsDone*(paths: openArray[string], expected: openArray[(string, string)]): bool
  ## All members must exist and match identifying metadata; incompatible files raise ValueError.
proc readTsv*(path: string): tuple[meta: Table[string, string], names: seq[string],
                                   cols: seq[seq[float]]]
  ## `columns` is the reserved key carrying the column names; it comes back as
  ## `names`, not in `meta`.  `#` lines with no `=` are comments.  Round trip is
  ## bit identical.
```

---

## 15. Test conventions

Tests live in `src/experimental/radial/tests/t*.nim`, use `std/unittest`, and are run with
`make run experimental/radial/tests/t<name>`. First line `#RUNCMD $RUN1` (single rank).

The **test ladder** every operator must pass, in this order:
1. explicit adjoint \(\langle x,Dy\rangle=\langle D^\dagger x,y\rangle\)
2. gauge covariance (including the \(\varphi\) seam and the temporal antiperiodic wrap)
3. matrix-free apply == dense basis assembly
4. analytic tangent vs centered finite differences over 4–5 step sizes, taking the **best**
5. pullback vs tangent contraction: \(\langle f,\delta u\rangle=2\,{\rm Re}\langle x,\delta D\,y\rangle\)
6. Ward identity \(\delta_{d\alpha}D=i(D\alpha-\alpha D)\)

and for HMC:
7. heatbath identity \(S=\|\xi\|^2\)
8. dense determinant / Hasenbusch telescoping oracle
9. per-frame force vs finite differences
10. reversibility; bitwise field restore on reject
11. exact integer force counts per trajectory
12. \(|\Delta H|\propto dt^2\)

Also mandatory: an **allocation regression** (`getRawMemAllocated()` unchanged over 64 applies).
