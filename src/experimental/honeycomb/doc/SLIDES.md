# Lattice QCD on the 16-cell honeycomb — slide-by-slide record

Source: <https://indico.global/event/16565/contributions/161563/attachments/74158/144712/16cell.pdf>

* **Talk:** "Lattice QCD on the 16-cell honeycomb"
* **Authors:** Sandor D. Katz and Daniel Nogradi (Eotvos University, Budapest)
* **Venue:** Lattice 2026 (43rd International Symposium on Lattice Field Theory),
  University of Maryland, College Park, 31 July 2026, 15:20,
  session *Algorithms and artificial intelligence*.
* **Companion paper:** arXiv:2512.10604 [hep-lat], "QCD on the 16-cell honeycomb",
  6 pages, 9 figures, submitted to Phys. Rev. D.
* **Deck:** 22 slides (LaTeX/beamer-like, black + crimson highlight, gnuplot figures).

Abstract of the talk (from Indico):

> The speakers propose discretizing lattice QCD on a 16-cell honeycomb lattice rather
> than the standard 4-dimensional hypercubic lattice.  This geometry offers "a larger
> symmetry" with three times more symmetry group elements.  Findings indicate "better
> scaling for gauge observables and better chiral properties for the Wilson-type Dirac
> operator."  Despite each lattice site having 24 neighbors instead of 8, they anticipate
> "significant" cost gains in dynamical simulations due to improved scaling behavior.

Abstract of the paper:

> We formulate QCD discretized on the four dimensional 16-cell honeycomb.  The advantage
> is a higher degree of rotational symmetry as compared to a traditional cubic lattice
> leading to much smaller cut-off effects.  We demonstrate in quenched QCD, through both
> gluonic and fermionic observables, that the scaling properties are indeed superior to
> the cubic lattice and much larger lattice spacings are sufficient for controlled
> continuum extrapolations.  Chiral and topological properties also show remarkable
> improvement.

---

## 1. Slide-by-slide transcription

### Slide 1 — Title
> **Lattice QCD on the 16-cell honeycomb**
> Sandor Katz and Daniel Nogradi
> Eotvos University, Budapest
> 2512.10604 [hep-lat]   Phys.Rev.D

### Slide 2 — Context and motivation
> * Symmetry important, lattice breaks it
> * Lorentz, Poincaré, chiral, …
> * Lots of effort: break as little as possible
> * Lattice: almost always 4D cubic
> * **Not most symmetric!**

### Slide 3 — Context and motivation (2D analogue)
> In 4D, 3 different regular space-filling tesselations
> 2D analogue
>
> *Figure:* three 2D tilings drawn side by side — square lattice (left), triangular
> lattice (middle), honeycomb lattice (right).  These are the three regular tilings of
> the plane {4,4}, {3,6}, {6,3}.

### Slide 4 — Context and motivation (symmetry orders)
> * 4D: also 3 regular space-filling tesselations
> * Cubic, 16-cell, 24-cell, **like in 2D, but not 3D**
> * Order of symmetry group:
> * cubic: **384**, 16-cell and 24-cell: **1152**
> * Very similar again to 2D

### Slide 5 — Context and motivation
> **Let's do lattice QCD on the 16-cell honeycomb**
> Larger symmetry → closer to continuum → larger lattice spacing → computational gain
> Old idea, 1980s,  *William Celmaster*

### Slide 6 — Context and motivation (history)
> Old idea
> * William Celmaster: gauge theory + naive fermions
> * Lorentz violating naive fermion doublers
> * Abandoned …
> * Bhanot/Heller/Neuberger: scalar field theory

### Slide 7 — The lattice
> Many ways to define, sites:
> * Cubic lattice ∪ (1/2, 1/2, 1/2, 1/2) shifted cubic lattice → body centered cubic lattice
> * (n₁, n₂, n₃, n₄) with n₁+n₂+n₃+n₄ even → D₄ lattice
>
> Links: nearest neighbors, **24 total** (8 for cubic)

### Slide 8 — The lattice (neighbour vectors)
> Nearest neighbor directions, 24 total
>
> (±1,0,0,0),  (0,±1,0,0),  (0,0,±1,0),  (0,0,0,±1)
> Usual cubic neighbors: 8
>
> (±1/2, ±1/2, ±1/2, ±1/2)
> 16 more

### Slide 9 — Gauge theory
> 12 links per site *x* (4 for cubic)
> Shortest closed loop: triangle, *P(x) = U₃U₂U₁*
>
> ```
> S = (β/2) Σ_x Σ_{i=1}^{32} ( 1 − (1/N) Re Tr P_i(x) )
> ```
>
> β = 2N/g² just as before
> Gradient flow straightforward
> Topological charge → clover-style hexagons from triangles

### Slide 10 — Topological susceptibility, SU(3)   *(= paper Fig. 2)*
> *Figure:* `10⁴ t₀² χ` (y, range 4 → 7) vs `a²/t₀` (x, range −0.1 → 2.1).
> * red squares "cubic [Cè et.al.]" for a²/t₀ ≈ 0.15 … 0.35, values ≈ 6.13 … 6.45,
>   with a straight red `O(a²)` fit extrapolating to ≈ 6.65 at a² = 0.
> * blue dots "16-cell" for a²/t₀ ≈ 0 … 2.0, values ≈ 6.78 down to ≈ 4.43,
>   with a curved blue `O(a⁴)` fit extrapolating to ≈ 6.78 at a² = 0.
> * The two extrapolations agree at a²=0 within errors (≈ 6.7).
>
> **O(a²) consistent with zero → probably very small**

### Slide 11 — Fermions
> Naive discretization, `n_i` the 24 neighbor vectors
>
> ```
> D₀ = (1/6) Σ_{i=1}^{24} γ_μ n_{μi} ∇_i
> ```
>
> Many doublers, some break Lorentz invariance  *(Celmaster)*
> Introduce Wilson term
>
> ```
> D = D₀ + a (r/6) Σ_{i=1}^{24} ∇*_i ∇_i
> ```
>
> Only 1 physical mode

### Slide 12 — Fermions, free pressure
> ```
> O = ( p(T) − p(0) ) / T⁴
>
> O_cubic   = (7π²/720) ( 1 + (248/147)(π²/N_t²) + (635/147)(π⁴/N_t⁴) … )
> O_16-cell = (7π²/720) ( 1 + (127/980)(π⁴/N_t⁴) + (73/4158)(π⁶/N_t⁶) + … )
> ```
> **16-cell: leading correction O(a⁴)**

### Slide 13 — Fermions, free pressure   *(= paper Fig. 4)*
> *Figure:* `p/p_cont` (y, 1 → 4) vs `N_t` (x, 3 → 20).
> * red squares "cubic": ≈ 3.83 (N_t=4), 2.80 (5), 2.11 (6), 1.70 (7), 1.47 (8),
>   1.32 (9), 1.23 (10), … slowly → 1 at N_t=20 (≈1.04).
> * blue dots "16-cell": ≈ 1.06 (N_t=4) and indistinguishable from 1 from N_t ≳ 6.
> * Solid curves are the truncated analytic series of slide 12.
>
> (Paper text: at N_t = 4 the correction to the continuum is 7 % and 283 % on the
> 16-cell and cubic lattices respectively.)

### Slide 14 — Fermions, free spectrum   *(= paper Fig. 3)*
> *Figure:* `Im(λ)` (y, −3 → 3) vs `Re(λ)` (x, 0 → 9) — free Wilson–Dirac eigenvalues
> on 16⁴ volumes.
> * red "cubic": the familiar three-lobed Wilson spectrum filling Re λ ∈ [0, 8],
>   |Im λ| ≤ 2.
> * blue "16-cell": a much more compact, nearly elliptical band, Re λ ∈ [0, ≈5.33],
>   |Im λ| ≲ 1.47, with two small empty "holes" near Re λ ≈ 4.3 and 4.7 on the real axis.
>
> **Close to a circle (ellipse) → Ginsparg-Wilson … ???**

### Slide 15 — Fermions, interacting
> Free improvement amazing, what about interacting case?
> * Dirac-spectrum, additive mass normalization (Wilson term)
> * Index theorem, topological charge, chirality

### Slide 16 — Additive mass normalization (Wilson term)   *(= paper Figs. 5 and 6)*
> *Left figure:* `Im(λ)` (−0.8 → 0.8) vs `Re(λ)` (0 → 0.2), quenched.
> * red squares "cubic a≈0.12 fm" — scattered points, left edge of the band around
>   Re λ ≈ 0.08.
> * magenta diamonds "cubic a≈0.06 fm" (a single 24⁴ configuration) — tight arc with
>   left edge ≈ 0.043.
> * blue dots "16-cell a≈0.12 fm" — smooth parabola-like arc with left edge ≈ 0.030.
>
> *Right figure:* histogram of `Re(λ₀)` (x from 0.02 to 0.12), 12⁴, a≈0.12 fm.
> * pink "cubic": broad distribution centred ≈ 0.09, width ≈ 0.01, peak height ≈ 0.14.
> * blue "16-cell": narrow distribution centred ≈ 0.031, peak height ≈ 0.38.
>
> **Fluctuates much less → real gain in dynamical simulations**

### Slide 17 — Index theorem, topology, chirality
> Wilson-term: no chiral symmetry
> `n_±`: number of real low-lying eigenvalues with non-zero chirality
>
> ```
> Q_Dirac = n₊ − n₋            Q_flow = Q_flow(t = t₀)
> ```
>
> Of course Q_Dirac ≠ Q_flow and chirality ≠ ±1

### Slide 18 — Index theorem, topology, chirality   *(= paper Fig. 7)*
> *Two histograms of `Q_Dirac − Q_flow`* (x from −3 to 3), grey band marking [−0.5, 0.5].
> * left, red: "cubic, 12⁴, a≈0.12 fm" — broad, peak height ≈ 0.085 at 0, visible
>   tails out to ±2.5.
> * right, blue: "16-cell, 12⁴, a≈0.12 fm" — much narrower, peak height ≈ 0.20 at 0,
>   tails essentially gone beyond ±1.
>
> **Distribution of Q_Dirac − Q_flow**

### Slide 19 — Index theorem, topology, chirality   *(= paper Fig. 8)*
> *Figure:* histogram of `chirality` (x from −1 to 1), 12⁴, a≈0.12 fm.
> * pink "cubic": broad, mass spread over |chirality| ∈ [0.5, 1], peaks ≈ 0.11–0.13
>   at ±0.9.
> * blue "16-cell": sharply peaked at ±1 (heights ≈ 0.44 and ≈ 0.40), almost nothing
>   in between.
>
> **Distribution of real low-mode chiralities**

### Slide 20 — Summary and outlook
> * Cubic → 16-cell lattice: larger space-time symmetry
> * Expect better scaling
> * Confirmed in both gauge and fermion sector
> * Currently: dynamical simulations — very promising

### Slide 21 — Summary and outlook (cost estimate)
> 3× more neighbors, 2× more sites: **6× cost**
> HMC ∼ V^{3/2}
> 2× gain in lattice spacing → 64× gain in HMC → fixed 6× cost penalty →
> **order of magnitude cheaper than cubic**
> Plus more from better chiral symmetry

### Slide 22 — Thank you for your attention!

---

## 2. What the authors actually did

Distilled from the deck plus the companion paper.

1. **Chose the 16-cell honeycomb** ({3,3,4,3}) as the space-time lattice instead of the
   hypercubic {4,3,3,4}.  Its vertex set is the D₄ lattice (equivalently the 4D
   body-centred cubic lattice `Z⁴ ∪ (Z+½)⁴`).  Each site has 24 nearest neighbours
   (the 4D kissing number) and 12 links.  The point symmetry group has 1152 elements
   (the Weyl group of F₄) versus 384 (the Weyl group of B₄) for the cubic lattice.

2. **Gauge action from triangles.**  The shortest closed loop is an equilateral triangle
   of side *a*.  Each site is the *apex* of 32 triangles (equivalently 96 triangles pass
   through each site, 192 counting orientation; each triangle has 3 corners).
   `S = (β/2) Σ_x Σ_{i=1}^{32} (1 − (1/N) Re Tr P_i(x))` with `β = 2N/g₀²`.

3. **Generated quenched SU(3) ensembles** with local heatbath + overrelaxation, O(50 000)
   independent configurations each, on volumes 6⁴, 8⁴, 9⁴, 10⁴, 11⁴, 12⁴, with
   periodic boundaries in the three (1,0,0,0)-type spatial directions (N_s sites) and
   (anti)periodicity in the (0,0,0,1) direction (N_t sites).

4. **Gradient flow and topology.**  Flow follows from the action; the only extra input is
   the non-unit area of a basic triangle.  A clover-style field strength is built from the
   16 *hexagons* at each site (6 coplanar triangles sharing a corner form a hexagon).
   From this, the topological charge density and Q.  Scale set with t₀
   (√(8t₀) = 0.47 fm convention quoted as "8t₀ = 0.47 fm").

5. **Measured the topological susceptibility** χ in t₀ units as a function of a²/t₀ and
   compared with the cubic-lattice continuum extrapolation of Cè et al.
   (arXiv:1506.06052).  The 16-cell data are flat out to a²/t₀ ≈ 1 and require an
   O(a⁴) fit; the O(a²) coefficient is consistent with zero.  Both extrapolate to the
   same continuum value ≈ 6.7 ×10⁻⁴/t₀², with χ²/dof 0.34 and 0.30.

6. **Wilson-type fermions.**  Naive operator `D₀ = (1/6) Σ_i γ·n_i ∇_i` plus a Wilson
   term over the same 24 directions; only one physical mode survives.

7. **Free-field analysis.**  Computed the free lattice pressure `O = (p(T)−p(0))/T⁴`
   exactly for N_t = 3…20 (N_s → ∞), and derived the asymptotic series in 1/N_t: the
   cubic lattice has an O(a²) leading correction, the 16-cell lattice's leading
   correction is O(a⁴).  Also plotted the full free Wilson–Dirac spectrum on 16⁴.

8. **Interacting Wilson–Dirac spectra.**  On 1215 quenched 12⁴ configurations at
   a ≈ 0.12 fm (plus one 24⁴ configuration at a ≈ 0.06 fm), with tree-level clover
   improvement and 6 steps of stout smearing (ϱ = 0.05), computed the lowest 300
   eigenvalues (100 on 24⁴).  Compared: the shape of the low-lying spectrum, the
   distribution of the real part of the lowest eigenvalue λ₀ (additive mass
   renormalisation), the distribution of Q_Dirac − Q_flow, and the distribution of the
   chirality of the real low modes.

9. **Cost argument.**  6× more work per site-volume, but HMC ∼ V^{3/2} and a factor 2
   coarser lattice spacing buys 2^6 = 64, so an order of magnitude net gain.

## 3. Cited works

| # | Reference |
|---|---|
| 1–7 | W. Celmaster, gauge theory + naive fermions on non-hypercubic lattices, 1982–1986 |
| 8–10 | G. Bhanot, U. Heller, H. Neuberger, scalar/Higgs field theory on such lattices, 1987–1992 |
| 11 | J. H. Conway, N. J. A. Sloane, *Sphere Packings, Lattices and Groups*, Springer 1999 |
| 12 | M. Lüscher, JHEP **1008** (2010) 071, arXiv:1006.4518 — gradient flow |
| 13 | M. Cè et al., Phys. Rev. D **92** (2015) 074502, arXiv:1506.06052 — cubic χ_top continuum limit |
| 14 | B. Sheikholeslami, R. Wohlert, Nucl. Phys. B **259** (1985) 572 — clover |
| 15 | C. Morningstar, M. Peardon, Phys. Rev. D **69** (2004) 054501, hep-lat/0311018 — stout smearing |
| 16 | P. Hernandez, Nucl. Phys. B **536** (1998) 345, hep-lat/9801035 — Wilson real modes / chirality |
