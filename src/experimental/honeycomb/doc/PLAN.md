# Reproduction plan — QEX implementation of the 16-cell honeycomb

Read [FORMULATION.md](FORMULATION.md) first; it is normative.  [SLIDES.md](SLIDES.md)
says what we are trying to reproduce.  Log progress in [STATUS.md](STATUS.md).

---

## 0. Ground rules for everyone working on this

* **Build**

  ```bash
  export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
  cd /Users/xjin/K/W/P003/qex/.claude/worktrees/reproduce-16-cell-slides-in-qex-d937db/build_mac
  make src/experimental/honeycomb/<file>.nim
  ./bin/<file>
  ```
  `make run <path>` builds and runs.  `SDKROOT` is mandatory — without it the C
  compiler cannot find system headers.

* **Style** — follow QEX: 2-space indent, `proc`/`template`, `import base, layout, field`
  style module imports, `threads:` blocks for field loops, `letParam:` +
  `installHelpParam(); echoParams()` for executables (see
  `src/examples/wflow_topo.nim` and `src/experimental/graph/pghmc.nim`).
  Reuse committed QEX code wherever it exists; do not re-invent
  `projectTAH`, `projectSU`, `newRNGField`, `saveGauge`, the RK integrators in
  `src/algorithms/rk.nim`, the solvers in `src/solvers/`, etc.

* **File ownership** — each task below owns a disjoint set of files.  Do not edit files
  owned by another task; if you need a change there, write it in STATUS.md instead.

* **Never run `git`.**  The worktree contains uncommitted work.

* **Every module gets a test** in `src/experimental/honeycomb/tests/` that runs
  standalone and prints `PASS`/`FAIL` lines, and exits non-zero on failure.

* When you finish, append a dated entry to `doc/STATUS.md`: what you built, the exact
  build/run commands, the numbers you got, and anything you could not resolve.

---

## 1. Architecture

The honeycomb is a hypercubic QEX `Layout` of **cells** (`N_s³ × N_t`) carrying a
2-site basis, plus 24 link fields per cell.  Nothing in QEX understands non-hypercubic
geometry, so the honeycomb layer is built on top of `Layout` + single-axis shifts.

```
  hcgeom.nim     pure combinatorics: 24 directions, 32 apex triangles, 16 hexagons,
                 point group.  No QEX Layout dependency.        [task G]
  hcfree.nim     free-field momentum space Dirac operator.      [task F]
  hclayout.nim   cell Layout + sublattices + the 16 delta-shifts.[task L]
  hcgauge.nim    24 gauge link fields, triangle loops, observables.[task L]
  hcaction.nim   triangle action, staples, force.               [task A]
  hchmc.nim      HMC for the triangle action.                   [task M]
  hcheatbath.nim SU(N) heatbath + overrelaxation (optional).    [task M]
  hcflow.nim     gradient flow (RK3) + t0.                      [task W]
  hctopo.nim     hexagon clover F_munu, E, Q.                   [task W]
  hcwilson.nim   interacting Wilson-Dirac operator.             [task D]
  hcarnoldi.nim  Arnoldi for lowest complex eigenvalues.        [task D]
```

Executables live next to the modules (`hcFreeSpectrum.nim`, `hcFreePressure.nim`,
`hcPureGauge.nim`, `hcMeasFlow.nim`, `hcSpectrum.nim`), tests in `tests/`.

### 1.1 The 16 diagonal shifts — do this the cheap way

The diagonal neighbours require shifts by `δ ∈ {0,1}⁴`.  QEX's `ShiftB`/`Shifter` only
supports single-axis displacements, and the general `comms/gather` path is slow.
**Build all 16 shifted copies with a binary tree: 15 single-axis shifts, not 64.**

```
  level 0:  f                                                  (δ = 0000)
  level 1:  shift(f, dir 0, +1)                                (δ = 0001)
  level 2:  shift each of the 2 by dir 1                       (4 fields)
  level 3:  ... dir 2                                          (8 fields)
  level 4:  ... dir 3                                          (16 fields)
```
Each level's shift is a standard `Shifter`, so this is SIMD-safe and MPI-safe.  Memory:
16 copies of a `ColorMatrix` field on 12⁴ cells ≈ 5.7 MB each — acceptable.

Alternative if memory becomes a problem: `src/comms/halo.nim`
(`makeHaloLayout`/`makeHaloMap` take an arbitrary list of offset vectors), which is the
only fully general SIMD-safe primitive in QEX.  Do not use `makeShiftSubQ` with a
general `disp`: its single `perm` field silently mis-permutes SIMD lanes when the
displacement crosses two or more SIMD-split directions.

### 1.2 Recommended first-pass simplifications

Correctness first.  Deliberately deferred optimisations, all recorded here so nobody
"fixes" them by accident:

* apply the full 4×4 `(γ·n − r)` instead of half-spinor projection (2× flops);
* recompute triangle products rather than caching shared sub-paths;
* use `newLayout(lat)` with default SIMD, but if a lane-permutation bug is suspected,
  fall back to `newLayout(lat, 1)` to isolate it.

---

## 2. Task breakdown

Dependency graph:

```
  G ──┬── F ────────────────────────────► Fig 3 (slide 14), Fig 4 (slides 12,13)
      └── L ── A ──┬── M ──┐
                   └── W ──┴── R ────────► Fig 1, Fig 2 (slide 10)
                            └── D ───────► Figs 5-8 (slides 16,18,19)
```

---

### Task G — geometry and symmetry  ⟶ *reproduces slides 4, 7, 8, and the counting on 9*

**Owns:** `hcgeom.nim`, `tests/tgeom.nim`

Deliver a compile-time/startup-time geometry table:

```nim
const
  nDirs* = 24            ## nearest-neighbour directions
  nAxis* = 8
  nDiag* = 16
  nTriPerSite* = 32      ## apex triangles
  nHexPerSite* = 16

type
  Dir* = array[4, float] ## one of the 24 unit neighbour vectors
  ApexTri* = object      ## triangle labelled by its apex
    d*, dp*: int         ## indices into the 16 diagonal directions
    mu*: int             ## the flipped component
  Hexagon* = object
    mu*: int             ## axis direction
    sigma*: array[4,int] ## +-1, sigma[mu] unused
    tri*: array[6, ...]  ## the 6 triangles, in cyclic order

proc dirs*(): array[24, Dir]
proc diagDirs*(): array[16, Dir]         ## d(delta)_mu = delta_mu - 1/2
proc apexTriangles*(): array[32, ApexTri]
proc hexagons*(): array[16, Hexagon]
proc omega*(h: Hexagon): array[4,array[4,float]]   ## the unit area 2-form, FORMULATION 2.3
proc pointGroupOrder*(vecs: openArray[Dir]): int
```

`tests/tgeom.nim` must check, and print, every ✅ item in FORMULATION §1–2:

1. 24 directions, all `|n| = 1`;
2. `Σ_i n_μ n_ν = 6 δ_μν` (and `2 δ_μν` for the 8 axis vectors alone);
3. exactly 32 unordered zero-sum triples; every one has 1 axis + 2 diagonal edges;
4. 96 triangles through a site; 32 have a unique apex there; 8 triangles per link;
5. 16 hexagons; each is 6 coplanar unit vectors at 60°; 16×6 = 96; 16×2 = 32;
6. `Σ_h Ω^{(h)}_{μν} Ω^{(h)}_{ρσ}`-type identity behind (4.2): verify numerically that
   `(3/8) Σ_h Ω^{(h)}_{μν} · (½ Ω^{(h)}_{αβ} F_{αβ}) = F_{μν}` for a random
   antisymmetric `F`;
7. `pointGroupOrder` = **1152** for the 24 directions and **384** for the 8 axis ones.

*(All seven have already been confirmed with an independent scratch calculation — the
numbers above are the answers, not guesses.  Item 6 is the only one not yet checked
numerically.)*

---

### Task F — free-field fermions  ⟶ *reproduces slide 14 (Fig 3) and slides 12–13 (Fig 4)*

**Owns:** `hcfree.nim`, `hcFreeSpectrum.nim`, `hcFreePressure.nim`, `tests/tfree.nim`
**Depends on:** G (direction list only).  **No QEX Layout needed** — this is 4×4 (cubic)
and 8×8 (16-cell, spin ⊗ sublattice) linear algebra in momentum space.

```nim
proc freeD*(p: array[4,float]; dirs: openArray[Dir]; pref, r, m: float): Mat4
  ## D(p) = M(p) + i gamma.K(p) ; see FORMULATION 5.3
proc freeEigs*(p: ...): array[2, Complex[float]]   ## M +- i|K|, exact - no eigensolver
```

**F1 — free spectrum (Fig 3).**  Enumerate all lattice momenta on 16⁴ (cubic: `p_μ =
2π k_μ/16`; 16-cell: `2 N⁴` momenta — see FORMULATION §5.3) and dump `Re λ, Im λ`.
Plot with gnuplot to match slide 14.
*Acceptance:* max `Re λ` = 8 (cubic) and 16/3 = 5.33333 (16-cell); max `|Im λ|` = 2 and
1.46789.

**F2 — free pressure (Fig 4).**  `O = (p(T)−p(0))/T⁴` for `N_t = 3..20`, `N_s → ∞`.
Recipe in FORMULATION §5.4.  Output `N_t, O_lat/O_cont` for both lattices.
*Acceptance:*
  * cubic ratio ≈ 3.83 at `N_t = 4`, 2.11 at 6, 1.47 at 8, → 1 as `N_t` grows;
  * 16-cell ratio ≈ 1.07 at `N_t = 4`, ≈ 1.00 from `N_t ≳ 6`;
  * fitting `O/O_cont − 1` to `c₂π²/N_t² + c₄π⁴/N_t⁴` (cubic) and
    `c₄π⁴/N_t⁴ + c₆π⁶/N_t⁶` (16-cell) reproduces
    `248/147, 635/147` and `127/980, 73/4158` to a few permille.
    **This fit is the sharpest single validation of the whole fermion formulation** —
    if it works, the operator and its normalisation are right.

`tests/tfree.nim`: the table in FORMULATION §5.3 plus the small-`p` limits
(`|K|/|p| → 1`, `M/(rp²/2) → 1`) for both lattices.

---

### Task L — lattice and gauge-field layer

**Owns:** `hclayout.nim`, `hcgauge.nim`, `tests/tgauge.nim`
**Depends on:** G

```nim
type
  HcLayout*[V:static[int]] = ref object
    lo*: Layout[V]              ## the CELL layout, N_s^3 x N_t
    ns*, nt*: int
  HcGauge*[F] = object
    uA*, uB*: array[4, F]       ## A(y)->A(y+e_mu), B(y)->B(y+e_mu)
    uD*: array[16, F]           ## B(y)->A(y+delta)

proc newHcLayout*(ns, nt: int): HcLayout
proc newHcGauge*(hl: HcLayout): HcGauge          ## all links = 1
proc random*(g: var HcGauge, r: var RNGField)
proc reunit*(g: var HcGauge)
proc gaugeTransform*(g: var HcGauge, v: ...)     ## needed for the invariance test
proc triangleSum*(g: HcGauge): float
  ## (1/(32 N_sites)) Sum_x Sum_{i=1}^{32} Re Tr P_i(x) / N   -- the "average triangle"
```

Also: the 16-way binary-tree shift helper (§1.1) as a reusable `HcShift16` object, since
both the action and the Dirac operator need it.

`tests/tgauge.nim`:
* unit gauge ⟹ `triangleSum = 1`;
* random gauge transformation leaves `triangleSum` invariant to `1e-12`;
* the 16-way shift reproduces a brute-force `lo.rankIndex`-based gather;
* every link appears in exactly 8 triangles (cross-check against `hcgeom`).

---

### Task A — action and force

**Owns:** `hcaction.nim`, `tests/taction.nim`
**Depends on:** L

```nim
proc hcAction*(beta: float, g: HcGauge): float
  ## S = (beta/2) Sum_x Sum_{i=1}^{32} (1 - Re Tr P_i / N)
proc hcForce*(beta: float, g: HcGauge, f: var HcGauge)
  ## traceless anti-Hermitian force, projectTAH convention as in QEX gaugeForce
```

The force needs the 8-staple sum per link (FORMULATION §3.2).  Reuse QEX
`projectTAH`/`contractProjectTAH`.

`tests/taction.nim`:
1. **finite-difference check**: for a random configuration and random algebra direction
   `H`, `d/ds S(e^{sH}U)|₀` from `hcForce` matches a central difference to `1e-8`;
2. **gauge invariance** of `hcAction`;
3. **classical continuum limit / β normalisation** (FORMULATION §3.1): put a smooth,
   weak Abelian field `A_μ(x)` (a single plane wave with small amplitude) on the
   honeycomb and on a cubic lattice of the same `a`, and check
   `S_16cell / S_cubic → 1` as the amplitude → 0 and `p a → 0`.
   Expect agreement to `O((pa)²)`;
4. strong-coupling: at small β, `⟨triangleSum⟩ ≈ β/(4N²)·(something)` — at minimum check
   that `⟨triangleSum⟩ → 0` as `β → 0` and `→ 1` as `β → ∞`.

---

### Task M — Monte Carlo

**Owns:** `hchmc.nim`, `hcPureGauge.nim`, `tests/thmc.nim`; optionally `hcheatbath.nim`
**Depends on:** A

**M1 — HMC.**  Momenta on all 24 link fields, leapfrog or 2MN (QEX has integrators in
`src/hmc/`; if they cannot be reused directly for a 24-field gauge object, write a
minimal leapfrog — it is ~60 lines).
*Acceptance:* reversibility `|ΔH_fwd+ΔH_bwd| < 1e-8`; `⟨exp(−ΔH)⟩ = 1` within errors
over ≥ 200 trajectories; acceptance > 70 %.

**M2 — heatbath + overrelaxation (strongly recommended for statistics).**
The paper used exactly this.  Each link's local action is `Re Tr(U Σ)` with `Σ` the sum
of 8 triangle staples — so the standard Cabibbo–Marinari SU(2)-subgroup heatbath and
microcanonical overrelaxation apply unchanged.  QEX has **no** SU(3) heatbath, so this
is new code (~200 lines).  It is 10–50× more efficient than HMC for quenched SU(3) and
is what makes the topological-susceptibility run feasible.
*Acceptance:* heatbath and HMC agree on `⟨triangleSum⟩` at the same β to within errors
on a small lattice.

**`hcPureGauge.nim`** — the generation executable: `letParam` for `ns, nt, beta, ntraj,
nsave, seed, outfile`; writes configurations (reuse `saveGauge` per link field, or a
simple raw dump — 24 fields).

---

### Task W — gradient flow and topology  ⟶ *reproduces Fig 1 and slide 10 (Fig 2)*

**Owns:** `hcflow.nim`, `hctopo.nim`, `hcMeasFlow.nim`, `tests/tflow.nim`, `tests/ttopo.nim`
**Depends on:** A

**W1 — flow.**  RK3 in the 2N form already used by `src/gauge/wflow.nim:21-66`; reuse
`src/algorithms/rk.nim` if convenient.  `wflow.nim` itself is hardwired to the
plaquette force (`gaugeForce(f,g)`), so this is a fork, not a call.

> **The flow-time normalisation must be calibrated numerically, not derived.**
> Required test: a weak Abelian plane wave must decay as `Â_μ(p,t) = e^{−t p²} Â_μ(p,0)`.
> Run the identical test on a cubic lattice with QEX's own `gaugeFlow` to validate the
> harness first, then fix the honeycomb constant.  See FORMULATION §4.4.

**W2 — topology.**  Hexagon clover (FORMULATION §4.1–4.3):
`F̂_Ω^{(h)} = (4/√3)·TAH[(1/6)C_h]`, then `F̂_μν = (3/8) Σ_h Ω^{(h)}_{μν} F̂_Ω^{(h)}`,
then `E` and `Q = ½ Σ_x q(x)`.
*(Resolved by task C: QEX's `topoQ` prefactor `−1/(4π²)` **is correct** — verified
exactly with an Atiyah–Singer constant-flux configuration; see FORMULATION §4.3.
Derive the honeycomb reduction independently anyway, and validate with the same
constant-flux test, which transfers directly and is sharper than integer-Q clustering.)*

**W3 — t₀.**  `t²⟨E⟩ = 0.3`; interpolate the flow-time series.  QEX has no t₀ finder.

`tests/tflow.nim`, `tests/ttopo.nim`:
* heat-kernel test above (this *is* the flow normalisation);
* on a smooth/heavily-flowed SU(3) configuration, `Q` must be within a few % of an
  integer, and `Q` must be stable under further flow;
* `F̂_μν` from the hexagon clover must reproduce the exact `F_μν` of a weak smooth
  Abelian field to `O(a²)`;
* the same for `E` — this validates the `4/√3`.

---

### Task R — physics runs  ⟶ *slide 10*

**Owns:** `hcAnalyze.nim`, `doc/RESULTS.md`, plot scripts under `doc/plots/`
**Depends on:** M, W

1. Generate quenched SU(3) ensembles.  The paper used 6⁴, 8⁴, 9⁴, 10⁴, 11⁴, 12⁴ with
   O(50 000) configurations; we will do what fits the machine — target ≥ 500–2000
   configurations per point, 3–5 β values, `N_s = N_t = 6..10`, at fixed physical
   volume `L/√t₀` where possible.
2. Measure `t₀` and `Q(t₀)` on each configuration.
3. `χ = ⟨Q²⟩/V` with `V = (a⁴/2)·N_sites = N_s³N_t a⁴`; plot `10⁴ t₀²χ` vs `a²/t₀`.
   Fit `O(a²)` and `O(a⁴)` forms.
   *Target:* 16-cell points flat out to `a²/t₀ ≈ 1`, extrapolating to `≈ 6.7×10⁻⁴`,
   with the `O(a²)` coefficient consistent with zero.
4. Histogram of `Q` (Fig 1).
5. Overlay the cubic reference points of Cè et al. (arXiv:1506.06052) — either from the
   paper's plot, or generate our own cubic data with QEX's existing `wflow_topo`
   machinery (cheaper to trust: **do both**, since QEX's cubic path is already tested).

**Reduced-scope fallback if statistics are the bottleneck:** run only 2–3 β on `6⁴`/`8⁴`
and demonstrate the *qualitative* result — the 16-cell `χ t₀²` is far flatter in `a²`
than the cubic one over the same range.  Say explicitly in RESULTS.md what statistics
were reached.

**Constraints learned from task C (cubic reference), binding for task R:**
* `t₀` does not exist below `L/√t₀ ≈ 6` (t²⟨E⟩ peaks below 0.3 and falls);
  aim for `L/√t₀ ≳ 9–10` on the points used for the χ fit.
* The `vlen:4` SIMD layout rejects `L = 10` and `L = 14` ("can't lay out inner
  geom"); usable sizes are `L ∈ {8, 12, 16, 20, 24}` (and small even L for tests).
* Thread oversubscription on this shared box is catastrophic (16 threads: 30×
  slowdown); run everything with `OMP_NUM_THREADS=4` and do not stack more than
  ~3 concurrent MC jobs.
* Fine `β` at fixed small `L` topologically freezes (cubic β=6.0/16⁴ froze);
  monitor the Q time series and exclude frozen points from fits.

---

### Task D — interacting Wilson–Dirac  ⟶ *slides 16, 18, 19 (stretch)*

**Owns:** `hcwilson.nim`, `hcarnoldi.nim`, `hcSpectrum.nim`, `tests/twilson.nim`
**Depends on:** L, W (for `Q_flow`), R (for configurations)

* **D1** — `hcwilson.nim`: `D ψ(x) = (m+4r)ψ + (1/6)Σ_i(γ·n_i − r)U_iψ(x+n_i)`.
  Needs 24 general `γ·n` matrices built from QEX's `gamma1..gamma4` — QEX's `spproj*`
  templates are axis-only and cannot be reused.
  *Test:* with `U = 1` on an `N⁴` cell lattice, the eigenvalues must match `hcfree.nim`'s
  momentum-space answer exactly.  **This is a strong end-to-end check that the
  interacting operator and the free analysis agree.**
* **D2** — stout smearing (ϱ = 0.05, 6 steps) on the honeycomb, plus tree-level clover.
  QEX's `src/gauge/stoutsmear.nim` is plaquette/axis-based; the smearing staples here
  are the triangle staples of task A, so this reuses `hcaction.nim`.
* **D3** — QEX has **no** non-Hermitian eigensolver (no Arnoldi/IRAM; LAPACK `zgeev` is
  bound at `src/eigens/linalgFuncs.nim:184` but never called).  Write an
  implicitly-restarted Arnoldi on top of that binding, or use `zgeev` on a small
  explicitly-formed matrix for tiny lattices as a first check.
* **D4** — measurements: distribution of `Re λ₀`; real modes and their chirality
  `⟨v|γ₅|v⟩`; `Q_Dirac = n₊ − n₋`; histogram of `Q_Dirac − Q_flow`.

---

## 3. What "reproduced" means, per slide

| slide | content | task | confidence |
|---|---|---|---|
| 4 | symmetry orders 384 / 1152 | G | **done analytically, ✅ verified** |
| 7, 8 | lattice + 24 neighbours | G | certain |
| 9 | 12 links, 32 triangles, β = 2N/g² | G, A | certain (analytic proof in FORMULATION §3.1) |
| 12 | free-pressure series | F | high — the series coefficients are a sharp test |
| 13 | free pressure vs N_t (Fig 4) | F | high |
| 14 | free Dirac spectrum (Fig 3) | F | **extremes already ✅ verified** |
| 10 | χ_top scaling (Fig 2) | R | medium — statistics-limited |
| — | Q distribution (Fig 1) | R | medium |
| 16 | λ₀ distribution (Figs 5, 6) | D | low — needs Arnoldi + ensembles |
| 18, 19 | Q_Dirac − Q_flow, chirality (Figs 7, 8) | D | low |
| 21 | 6× cost accounting | — | arithmetic, stated in FORMULATION §6 |

Priority order: **G → F** (two full slides, exact, cheap) **→ L → A → M/W → R → D**.
