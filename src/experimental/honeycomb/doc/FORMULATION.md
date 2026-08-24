# The 16-cell honeycomb: formulation, conventions and verified numbers

This is the normative reference for `src/experimental/honeycomb`.  Everything here has
been derived and, where marked ✅, **verified numerically** (see
[VERIFIED.md](VERIFIED.md) for the reproduction commands).
If code disagrees with this document, one of the two is wrong — resolve it, don't paper
over it.

Companion documents: [SLIDES.md](SLIDES.md) (what the talk shows),
[PLAN.md](PLAN.md) (work breakdown), [STATUS.md](STATUS.md) (running log).

---

## 1. The lattice

### 1.1 Two equivalent descriptions

The 16-cell honeycomb {3,3,4,3} has as its vertex set the **D₄ lattice**.  The slides give
two descriptions; they are related by a similarity transformation (D₄ is similar to its
dual D₄* in 4D, with ratio 1/√2):

| name | site set | nearest-neighbour vectors | NN distance |
|---|---|---|---|
| **D₄** | `{n ∈ Z⁴ : Σ n_μ even}` | 24 vectors `(±1,±1,0,0)` & perms | √2 |
| **D₄\*** ("4D bcc") | `Z⁴ ∪ (Z+½)⁴` | 8 × `(±1,0,0,0)` + 16 × `(±½,±½,±½,±½)` | 1 |

**We use D₄\* throughout.**  Reasons: (a) it is the description on slides 7–8;
(b) the nearest-neighbour distance equals the cubic-sublattice spacing, so `a` has the
same meaning as on a cubic lattice and volumes/`N_t` are directly comparable;
(c) it maps onto a hypercubic QEX `Layout` of *cells* with a 2-site basis.

### 1.2 Cell + sublattice encoding (the implementation representation)

A hypercubic **cell** lattice of size `N_s × N_s × N_s × N_t` (QEX `Layout`), with two
sites per cell:

```
  sublattice A (s=0):  x = y                     y ∈ Z⁴ (the cell index)
  sublattice B (s=1):  x = y + (½,½,½,½)
```

* number of sites  = `2 · N_s³ · N_t`
* number of links  = `12 · N_sites` = `24 · N_cells`
* volume per site  = `a⁴/2`   ← **this factor of 2 appears in every extensive observable**
* physical extents `L = N_s a`, `T = N_t a` — directly comparable to a cubic `N_s³×N_t`.

Boundary conditions (paper): periodic with period `N_s` in the (1,0,0,0),(0,1,0,0),
(0,0,1,0) directions, (anti)periodic with period `N_t` in (0,0,0,1).  Since these
periods are integer vectors, they map A→A and B→B, so any `N_s`, `N_t` is allowed.

### 1.3 The 24 neighbour vectors

```
  axis      n = ± ê_μ                            μ = 0..3            (8 vectors)
  diagonal  n = ½(σ₀,σ₁,σ₂,σ₃),  σ_μ ∈ {±1}                        (16 vectors)
```
All have `|n| = 1`.  Axis vectors connect A↔A and B↔B; diagonal vectors connect A↔B.

✅ `Σ_{i=1}^{24} n_{μi} n_{νi} = 6 δ_{μν}`  (cubic: `Σ_{i=1}^{8} = 2 δ_{μν}`).
This single identity fixes every `1/6` normalisation below.

### 1.4 Link indexing (canonical, used by all modules)

24 gauge-link fields per cell.  Let `y` be the cell index, `δ ∈ {0,1}⁴`, `μ ∈ {0,1,2,3}`.

| field | index | link | displacement |
|---|---|---|---|
| `uA[μ]` | `μ = 0..3` | `A(y) → A(y+ê_μ)` | `+ê_μ` |
| `uB[μ]` | `μ = 0..3` | `B(y) → B(y+ê_μ)` | `+ê_μ` |
| `uD[δ]` | `δ = 0..15`, bits `δ_μ` | `B(y) → A(y+δ)` | `½(2δ₀−1, 2δ₁−1, 2δ₂−1, 2δ₃−1)` |

The B→A displacement for `uD[δ]` is `(y+δ) − (y + ½1) = δ − ½1 = ½(2δ−1)`, which runs
over all 16 diagonal vectors as `δ` runs over `{0,1}⁴`.  Every diagonal link is
represented exactly once, oriented B→A.  Conversely `A(y) → B(y−δ)` is `uD[δ]†` living
on cell `y−δ`.

Flat index convention: `δ = δ₀ + 2δ₁ + 4δ₂ + 8δ₃`; `uD[δ]` diagonal direction
`d(δ)_μ = δ_μ − ½`.

### 1.5 Point symmetry group

The stabiliser of a site — the set of 4×4 orthogonal matrices mapping the neighbour set
onto itself — has order

✅ **1152** for the 16-cell honeycomb (Weyl group of F₄)
✅ **384** for the cubic lattice (Weyl group of B₄ = signed permutations, `2⁴·4!`)

exactly the numbers on slide 4.  (The 24-cell honeycomb has the same 1152.)

---

## 2. Triangles, hexagons and the counting on slide 9

### 2.1 Triangles

A triangle is a closed 3-link loop.  Solving `n_i + n_j + n_k = 0` over the 24 vectors:

* ✅ every triangle has **exactly one axis edge and two diagonal edges**;
* ✅ **32** unordered zero-sum triples of neighbour vectors;
* ✅ **96** triangles pass through each site (192 counting orientation — slide 9's number);
* since a triangle has 3 corners, `96/3 = ` **32 triangles per site** — the `Σ_{i=1}^{32}`
  on slide 9;
* ✅ **8** triangles contain each link (`64 triangles/cell × 3 / 24 links/cell`).

All triangles are **equilateral with side `a`**, hence area `√3 a²/4`.

### 2.2 Canonical enumeration: the apex labelling

Each triangle has a unique **apex** — the corner at which both edges are diagonal (the
corner not on the axis edge).  This gives a bijection

```
  triangles  ⟷  (site x, unordered pair {d, d'} of diagonal directions
                 differing in exactly one component sign)
```

Per site: 16 diagonals × 4 sign-flips / 2 = **32**.  ✅

Label a triangle by `(x, d, μ)` where `d` is a diagonal direction and `μ` the flipped
component, `d' = d − 2 d_μ ê_μ`.  Its three vertices are `x`, `x+d`, `x+d'`; the axis
edge `x+d → x+d'` has displacement `−2 d_μ ê_μ = ∓ê_μ`.

*Implementation note.* Because the apex of every triangle is unique, iterating
`for x in sites: for (d,μ) in 32 apex labels:` covers each triangle exactly once with no
double-counting logic — this is the loop the action uses.

In cell/sublattice terms (32 per site × 2 sites per cell = **64 triangles per cell**):

* **apex on A**, cell `y`: edges `A(y)→B(y−δ)` and `A(y)→B(y−δ')` with `δ' = δ ⊕ 2^μ`;
  the axis edge joins the two B sites.  Loop:
  `uD[δ](y−δ) · uB[μ] or uB[μ]† · uD[δ'](y−δ')†`  (see §2.4).
* **apex on B**, cell `y`: edges `B(y)→A(y+δ)` and `B(y)→A(y+δ')`; the axis edge joins
  two A sites.

### 2.3 Hexagons (needed for the clover)

Six coplanar neighbour vectors `{±u, ±v, ±(u+v)}` with `|u|=|v|=|u+v|=1` (i.e.
`u·v = −½`) form a regular hexagon, and the 6 triangles spanned by consecutive pairs
tile it.

✅ **16 hexagons at each site**, exactly as the paper states.  Structure: each hexagon
contains exactly one axis pair `±ê_μ` and two diagonal pairs `±d⁺, ±d⁻` with
`d^± = ½(±1, σ₁, σ₂, σ₃)` (component `μ` flipped between them).  Labelling:
`μ ∈ {0..3}` (4) × `(σ_ν)_{ν≠μ}` up to overall sign (4) = **16**.  ✅

Consistency: 16 hexagons × 6 triangles = 96 triangles through the site ✅, and each
hexagon contributes exactly 2 apex-triangles at its centre, 16×2 = 32 ✅.

Cyclic order of a hexagon (60° steps): `d⁻, d⁺, ê_μ, −d⁻, −d⁺, −ê_μ`, using
`d⁺ = d⁻ + ê_μ`.

Unit area 2-form of hexagon `h = (μ, σ)`: plane spanned by `ê_μ` and
`f = (1/√3) Σ_{ν≠μ} σ_ν ê_ν`, so

```
  Ω^{(h)}_{μν} = σ_ν/√3  (ν ≠ μ),      Ω^{(h)}_{νμ} = −σ_ν/√3,   else 0.
```

---

## 3. Gauge action

```
  S = (β/2) Σ_x Σ_{i=1}^{32} ( 1 − (1/N) Re Tr P_i(x) ),     β = 2N/g₀²
```
with `x` over all `2 N_s³ N_t` sites and `i` over the 32 apex-triangles at `x`.

### 3.1 Normalisation check (why β = 2N/g² is right)

For a planar loop with area 2-form `S_{αβ}`, `P = exp(i g Φ)`, `Φ = ½ S_{αβ} F_{αβ}`,
and `1 − (1/N)Re Tr P ≈ (g²/2N) Tr Φ²`.  For apex-triangle `(d,μ₀)`,
`S_{μ₀ν} = d_{μ₀} d_ν` (ν≠μ₀), giving `Φ = d_{μ₀} Σ_{ν≠μ₀} d_ν F_{μ₀ν}`.  Summing over
the 32 triangles at a site with `Σ_{16 d} d_ν d_ρ = 4 δ_{νρ}`:

```
  Σ_{i=1}^{32} Tr Φ_i² = ½ Σ_{μ≠ν} Tr F_{μν}²
```

so `S = ¼ Σ_x Σ_{μν} Tr F² `.  With `Σ_x = (2/a⁴) ∫d⁴x`:

```
  S = ½ ∫ d⁴x Σ_{μν} Tr F_{μν} F_{μν}
```

— **identical** to the cubic Wilson action `S = β Σ_x Σ_{μ<ν}(1−Re Tr P/N)` with
`β = 2N/g²`.  So the `β/2` prefactor together with the 32-triangle sum is exactly right,
and `β` means the same thing on both lattices. ✅ (analytic; to be re-checked numerically
by the smooth-field test, task **G3**).

### 3.2 Loop ordering and the 8 staples per link

For apex `x`, pair `(d, d')`:
`P = U(x→x+d) · U(x+d→x+d') · U(x+d'→x)`.
The staple of a link is the product of the other two edges; each link sits in 8
triangles, so the force needs an 8-term staple sum per link.

* axis link `A(y)→A(y+ê_μ)`: apexes are the 8 B-sites adjacent to both, i.e. cells
  `y−δ` with `δ_μ = 0`.
* diagonal link `B(y)→A(y+δ)`: 4 triangles with apex `B(y)` (flip any of the 4
  components of `δ`) and 4 with apex `A(y+δ)`.

---

## 4. Field strength, energy density, topological charge

### 4.1 Hexagon clover

For hexagon `h` at site `x`, let `C_h(x) = Σ_{k=1}^{6} P_k(x)` be the sum of its 6
triangle loops, all traversed with the same orientation and all based at `x`.  Each has
area `√3a²/4`, so with `F_Ω ≡ ½ Ω_{αβ} F_{αβ}` (the field strength in the hexagon plane)

```
  TAH[ (1/6) C_h(x) ]  ≈  (√3/4) a² F_Ω^{(h)}(x)
```

where `TAH` is the traceless anti-Hermitian projection `(M − M†)/2 − tr/N`.  Hence

```
  F̂_Ω^{(h)} ≡ a² F_Ω^{(h)}  =  s · (4/√3) · TAH[ (1/6) C_h(x) ]     ... (4.1)
```

with an orientation sign `s = ±1` fixed by the handedness of the hexagon traversal
relative to `Ω^{(h)}`.  ✅ Task W found that `hcgeom.hexTriPaths`' ring order is
*clockwise* w.r.t. `omega`, so the code uses `hcCloverSign = −1`; pinned by the
site-wise weak-field test (a sign error gives ratio −1, magnitudes `4/√3`, `3/8`, `½`
all confirmed) and by the Atiyah–Singer constant-flux test.

### 4.2 Reconstructing F_μν

Summing the 16 hexagons gives `Σ_h Ω^{(h)}_{μν} F_Ω^{(h)} = (8/3) F_{μν}`, so

```
  F̂_μν(x) = (3/8) Σ_{h=1}^{16} Ω^{(h)}_{μν} F̂_Ω^{(h)}(x)             ... (4.2)
```

Because `Ω^{(h)}_{μν}` is nonzero only when one of `μ,ν` equals the hexagon's axis
direction, (4.2) collapses to: for each `(μ,ν)`, the 4 hexagons with axis `μ` plus the
4 with axis `ν` contribute, with weights `±σ/√3`.

### 4.3 E and Q

```
  E(x) = −½ Σ_{μν} Tr[ F̂_μν F̂_μν ] / a⁴          ⟨E⟩ = (1/N_sites) Σ_x E(x)·a⁴
  q(x) = (1/32π²) ε_{μνρσ} Tr[ F̂_μν F̂_ρσ ]
  Q    = ½ Σ_{x ∈ sites} q(x)                    ← the ½ is a⁴/2 per site!
```

**The `½` in Q is the single most likely source of a wrong answer.**  It comes from
`Σ_x (a⁴/2) = ∫d⁴x`.  `⟨E⟩` is an average so it needs no such factor.

QEX's `topoQ` (`src/gauge/gaugeUtils.nim:1274`) uses prefactor `−1/(4π²)` on
`(a−b+c)` with `a = ΣRe tr(F₁₀F₃₂)` etc.  ✅ **That normalisation is CORRECT** — an
earlier revision of this document suspected a factor 2; the suspicion came from
mis-pairing `1/32π²` with `ε F^a F^a` (the textbook identities are
`Q = (1/32π²)∫ ε_{μνρσ} Tr[F_{μν}F_{ρσ}] = (1/64π²)∫ ε F^a F^a`, since
`F^aF^a = −2Tr FF`).  Task C proved it exactly with an Atiyah–Singer constant-flux
Cartan configuration (`T = diag(1,−1,0)`, fluxes `n₁,n₂`, `Q = 2n₁n₂`): QEX returns
1.999598 for exact 2 on 16⁴, the residual being the pure `1−(φ₁²+φ₂²)/6` clover
artefact, scaling as `1/L⁴`.  **Use the same constant-flux test on the honeycomb** —
it is far sharper than "Q near an integer" (with an unimproved clover at
`a²/t₀ ≳ 0.5`, flowed `Q(t₀)` routinely misses the nearest integer by 0.2–0.45).

### 4.4 Gradient flow

Lüscher flow `V̇_ℓ = −g₀² (∂_ℓ S) V_ℓ`, integrated with the RK3 scheme already in
`src/gauge/wflow.nim`.

**Normalisation caveat.**  Both `S` and the flow equation are formally the same as on the
cubic lattice, but the *density of links* differs (12 links/site here vs 4, and half the
volume per site), so the induced continuum flow rate `∂_t B_μ = c · D_ν G_{νμ}` has a
different `c`.  A naive transcription is expected to be off by a constant factor
(estimated `1/6`, but the derivation has factor-order ambiguities and **must not be
trusted**).

The constant is therefore **fixed numerically, not analytically**, by the free-field
test:

> For a weak Abelian plane-wave gauge field, the flow must act as the heat kernel,
> `Â_μ(p,t) = e^{−t p²} Â_μ(p,0)`.

This determines the flow-time normalisation to machine precision and is cheap.  It is a
required deliverable (task **F2**).  The same test on the cubic lattice validates the
harness against `src/gauge/wflow.nim`.

Independently, the *observable* `E` needs the triangle-area factor `4/√3` of (4.1) — this
is what the paper means by "one only needs to take into account the non unit area of a
basic triangle".

`t₀` is defined by `t² ⟨E⟩ |_{t=t₀} = 0.3` (Lüscher).  QEX has **no t₀ finder**; we write
one (interpolate the `t²E` curve).  Scale setting in the paper: `√(8t₀) = 0.47 fm`.

---

## 5. Fermions

### 5.1 The operator

Naive term (slide 11), with `∇_i f(x) = (U_i(x) f(x+n_i) − f(x))/a`:

```
  D₀ = (1/6) Σ_{i=1}^{24} (γ·n_i) ∇_i
```
The `−f(x)` pieces cancel because `Σ_i n_i = 0`.  ✅ In momentum space
`D₀(p) = (i/6a) Σ_i (γ·n_i) sin(a p·n_i) → i γ·p` — this is what the `1/6` is for,
given `Σ_i n_μ n_ν = 6 δ_{μν}` (§1.3).

Wilson term.  Slide 11 writes `+ a(r/6) Σ_i ∇*_i ∇_i`, but also states that the
normalisation is chosen so the Wilson term is `a r p²/2` at small `p`, "as usual".
Those two statements differ by a factor 2, because `Σ_{i=1}^{24} ∇*_i∇_i` double-counts
the ±n pairs.  **The `a r p²/2` convention is the correct one** — it is fixed
independently by the maximal real eigenvalue (§5.3), which matches the slide-14 plot.
So:

```
  W(p) = (r/(6a)) Σ_{i=1}^{24} ( 1 − cos(a p·n_i) )   →   a r p²/2      ✅
       ≡ −(a r/12) Σ_{i=1}^{24} ∇*_i ∇_i
```

**Position space (the form to implement):**

```
  D ψ(x) = (m + 4r/a) ψ(x) + (1/(6a)) Σ_{i=1}^{24} ( γ·n_i − r ) U_i(x) ψ(x+n_i)
```

Cubic analogue for reference (same expression with the 8 axis vectors and `1/(2a)`):
`D ψ = (m+4r/a)ψ + (1/2a) Σ_{i=1}^{8} (γ·n_i − r) U_i ψ(x+n_i)`.  ✅ reduces to the
textbook Wilson operator.

### 5.2 Projectors

`(γ·n_i)² = |n_i|² = 1` for **all 24** directions, so for `r = 1`

```
  (γ·n_i − 1)/2 = −P_i ,    P_i = (1 − γ·n_i)/2   a rank-2 projector
```

exactly as in the cubic case.  The half-spinor trick therefore carries over verbatim,
with 24 projectors instead of 8.  QEX's `spproj*` templates are hardcoded per axis
direction (`src/physics/spinOld.nim:351-379`) and cannot express `γ·n` for a diagonal
`n`; we need our own.  *For the first correct implementation, apply the full 4×4
`(γ·n − r)` and skip projection; optimise later.*

For a diagonal direction `n = ½(σ₀,σ₁,σ₂,σ₃)`,
`γ·n = ½ Σ_μ σ_μ γ_μ` — a fixed 4×4 complex matrix per `σ`, computable once at
start-up from QEX's `gamma1..gamma4` (`src/physics/spinOld.nim:242-266`, DeGrand–Rossi
convention, 1-based names, `gamma5 = diag(1,1,−1,−1)`).

### 5.3 Free spectrum — the numbers to hit

Free (`U = 1`), massless, lattice units `a = 1`.  `D(p) = M(p) + i γ·K(p)`, so the
eigenvalues are `M(p) ± i|K(p)|`, each doubly degenerate:

```
  M(p) = (r/6) Σ_i ( 1 − cos(p·n_i) )
       = 4r − (r/6)[ 2 Σ_μ cos p_μ + 16 Π_μ cos(p_μ/2) ]
  K_μ(p) = (1/6)[ 2 sin p_μ + 8 sin(p_μ/2) Π_{ν≠μ} cos(p_ν/2) ]
```

| quantity | 16-cell (r=1) | cubic (r=1) |
|---|---|---|
| max Re λ | ✅ **16/3 = 5.33333** (at `p=(π,π,π,π)`) | ✅ **8** |
| max \|Im λ\| | ✅ **1.467890 = 3^{1/4}(1+√3)/√6** (`p ∝ (1,1,1,1)`; on-axis at `p=(2·arccos((√3−1)/2),0,0,0) = (2.392124,0,…)`, NOT at `2π/3` which gives only 1.44338) | ✅ **2** |
| small-p `\|K\|/\|p\|` | ✅ 1.000000 | ✅ 1.000000 |
| small-p `M/(r p²/2)` | ✅ 1.000000 | ✅ 1.000000 |

The slide-14 plot shows exactly this: the blue 16-cell cloud reaching `Re λ ≈ 5.3` and
`|Im λ| ≈ 1.47`, versus red cubic `[0,8] × [−2,2]`.

Brillouin zone: the site set `Z⁴ ∪ (Z+½)⁴` (covolume ½) has reciprocal lattice `2π·D₄`
(covolume 2), so the BZ has volume `2·(2π)⁴` — twice the cubic one, consistent with 2
sites per unit hypercube.  *(An earlier revision of this document said `(2π)⁴/2`; that
was wrong — caught and pinned by `tests/tfree.nim`.)*  On a finite `N_s³×N_t`
cell lattice there are `2 N_s³ N_t` momenta; with the cell encoding they are most easily
enumerated as the `N_s³N_t` cell momenta × the 2-dim sublattice space (the free operator
becomes an 8×8 matrix in (spin ⊗ sublattice) at each cell momentum).

### 5.4 Free pressure

```
  O = ( p(T) − p(0) ) / T⁴
```
normalised per fermionic degree of freedom, so the continuum value is
`7π²/720 = (7/8)(π²/90)`.  Targets from slide 12:

```
  O_cubic   = (7π²/720)[ 1 + (248/147)(π²/N_t²) + (635/147)(π⁴/N_t⁴) + … ]
  O_16-cell = (7π²/720)[ 1 + (127/980)(π⁴/N_t⁴) + (73/4158)(π⁶/N_t⁶) + … ]
```

Numeric coefficients: `248/147 = 1.68707`, `635/147 = 4.31973`,
`127/980 = 0.129592`, `73/4158 = 0.0175565`.

Slide 13 data (read off the plot, `p/p_cont`):

| N_t | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 12 | 20 |
|---|---|---|---|---|---|---|---|---|---|
| cubic | 3.83 | 2.80 | 2.11 | 1.70 | 1.47 | 1.32 | 1.23 | 1.14 | ≈1.04 |
| 16-cell | ≈1.06 | ≈1.02 | ≈1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |

(The paper quotes "at N_t = 4 the correction is 7 % and 283 %", consistent with 1.07 and
3.83.)

**How to compute it.**  Free lattice fermions, massless, antiperiodic in time with `N_t`
sites, `N_s → ∞`:

```
  p = (T/V) ln Z ,     ln Z = Σ_{p} ln det D(p)
  O = ( p(T) − p(0) ) / T⁴ ,   T = 1/(N_t a)
```
i.e. `O = N_t⁴ · [ (1/N_t) (1/V_s) Σ_{p} ln det D(p) ]_{T} − (same at T=0)`, with the
spatial momentum sum replaced by an integral (`N_s → ∞`) and the `T=0` piece being the
same integral with a continuous `p₄`.  The `det` is over spin (×sublattice for the
16-cell, ×colour trivially, divided out by the per-dof normalisation).
Practical recipe: evaluate the spatial integral by a large-`N_s` Riemann sum with
Richardson/`N_s→∞` extrapolation, and subtract the `T=0` term computed as a 4D integral
with the same spatial rule.  Cross-check against the analytic series above.

### 5.5 Improvement / smearing used in the paper for the interacting spectra

Tree-level clover improvement (`c_SW = 1`) and 6 steps of stout smearing with `ϱ = 0.05`.
1215 quenched 12⁴ configurations at `a ≈ 0.12 fm`; lowest 300 eigenvalues (100 on the
single 24⁴ configuration at `a ≈ 0.06 fm`).

`Q_Dirac = n₊ − n₋` counting real low-lying eigenvalues by the sign of their chirality
`⟨v|γ₅|v⟩/⟨v|v⟩`; `Q_flow = Q(t = t₀)`.

---

## 6. Cost accounting (slide 21)

| | cubic | 16-cell |
|---|---|---|
| sites per `a⁴` | 1 | 2 |
| neighbours / site | 8 | 24 |
| links / site | 4 | 12 |
| loops in the action / site | 6 plaquettes | 32 triangles |
| work ratio | 1 | **6** |

HMC cost `∼ V^{3/2}`; a factor 2 coarser `a` at fixed physics gives `2^{4·3/2} = 64`,
against the fixed 6× penalty → ≈ 10× net gain.

---

## 7. Conventions fixed for the code

* `a = 1` in all code; the cubic-sublattice spacing equals the NN distance.
* Directions: `μ = 0,1,2,3`; `μ = 3` is time (matches QEX's `setBC`, `densityE`).
* Sublattice: `0 = A` (integer coords), `1 = B` (half-integer).
* Diagonal index `δ ∈ 0..15`, bit `μ` = `δ_μ`; direction `d(δ)_μ = δ_μ − ½`.
* Gauge links: `uA[μ]`, `uB[μ]` (μ=0..3), `uD[δ]` (δ=0..15) — see §1.4.
* γ matrices: QEX `gamma1..gamma4` from `src/physics/spinOld.nim`
  (1-based names; `gamma0` is the identity, `gamma5 = diag(1,1,−1,−1)`).
* `TAH(M) = (M − M†)/2 − Tr(M − M†)/(2N)`, i.e. QEX `projectTAH`.
* Extensive gluonic observables carry the `a⁴/2 per site` measure; intensive ones
  (averages) do not.
