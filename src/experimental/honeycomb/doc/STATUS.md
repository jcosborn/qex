# STATUS — running log

Append a dated entry when you finish a task.  Say what you built, the exact
build/run commands, the numbers you got, and anything unresolved.
Do **not** edit files owned by another task (see PLAN.md §0).

---

## 2026-08-21 — bootstrap (lead)

**Build recipe** (mandatory `SDKROOT`, otherwise clang cannot find system headers):

```bash
export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
cd /Users/xjin/K/W/P003/qex/.claude/worktrees/reproduce-16-cell-slides-in-qex-d937db/build_mac
make src/experimental/honeycomb/hcgeom.nim && ./bin/hcgeom
```

The build directory `build_mac/` was created inside this worktree with the standard
project options (clang-mp-22, `-Ofast -march=native`, QMP/QIO, `vlen:4`).

**Done**

* `doc/SLIDES.md` — full slide-by-slide record of the 22-slide deck plus what the
  authors did, with the paper's figures/parameters cross-referenced.
* `doc/FORMULATION.md` — normative formulation: lattice, link indexing, triangle and
  hexagon combinatorics, action normalisation proof, clover construction, the Dirac
  operator and its verified free-field numbers, conventions.
* `doc/PLAN.md` — task breakdown G/F/L/A/M/W/R/D with interfaces and acceptance tests.
* `hcgeom.nim` — geometry module (task **G**, module part).

**Verified numerically** (independent scratch calculation, then reproduced by
`hcgeom.nim`):

| quantity | value | status |
|---|---|---|
| nearest neighbours, all `\|n\|=1` | 24 | ✅ |
| `Σ_i n_μ n_ν` | `6 δ_μν` | ✅ |
| zero-sum triples = apex triangles / site | 32 | ✅ |
| triangles through a site | 96 | ✅ |
| hexagons per site | 16 | ✅ |
| **point group order, 16-cell** | **1152** | ✅ `./bin/hcgeom` |
| **point group order, cubic** | **384** | ✅ `./bin/hcgeom` |
| free Wilson–Dirac max `Re λ`, 16-cell (r=1) | 16/3 = 5.33333 | ✅ scratch |
| free Wilson–Dirac max `\|Im λ\|`, 16-cell | 1.46789 | ✅ scratch |
| free Wilson–Dirac max `Re λ` / `\|Im λ\|`, cubic | 8 / 2 | ✅ scratch |
| small-`p`: `\|K\|/\|p\|`, `M/(r p²/2)` both lattices | 1, 1 | ✅ scratch |

**Slide 4 is reproduced** (1152 vs 384).

**Open questions handed to the tasks**

1. *Flow-time normalisation* (task W).  The continuum flow rate depends on the link
   density, which differs from the cubic case.  Analytic derivation is
   factor-ambiguous; it **must** be calibrated by the Abelian heat-kernel test.
2. *QEX `topoQ` prefactor* (`src/gauge/gaugeUtils.nim:1274`, `−1/(4π²)`) appears to be a
   factor 2 off the textbook `1/(32π²)ε Tr FF`.  Derive ours independently; validate on
   a cooled configuration with near-integer Q.
3. *Slide 11's Wilson-term factor 2.*  `a(r/6)Σ∇*∇` and "`arp²/2` as usual" are
   inconsistent; the latter is right and is what FORMULATION §5.1 fixes (confirmed by
   `max Re λ = 16/3` matching the slide-14 plot).
4. *`Q = ½ Σ_x q(x)`* — the volume per site is `a⁴/2`.  Easy to forget.

---

## 2026-08-21 — task L (lattice + gauge layer) and the task G test

**Files added**

| file | contents |
|---|---|
| `tests/tgeom.nim` | task **G** test: all 9 acceptance items of PLAN.md §Task G |
| `hclayout.nim` | `HcLayout` (cell layout) + `HcShift16` (binary-tree 16-way shift) |
| `hcgauge.nim` | `HcGauge` (24 link fields), `triangleSum`, `gaugeTransform` |
| `tests/tgauge.nim` | task **L** test, incl. a brute-force reference for everything |

**Build / run**

```bash
export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
cd /Users/xjin/K/W/P003/qex/.claude/worktrees/reproduce-16-cell-slides-in-qex-d937db/build_mac
make src/experimental/honeycomb/tests/tgeom.nim  && ./bin/tgeom
make src/experimental/honeycomb/tests/tgauge.nim && OMP_NUM_THREADS=4 ./bin/tgauge
```

Both print `PASS:` lines and exit 0; `std/unittest` sets the exit code to 1 on any
failure (verified).  **`tgeom`: 42/42 PASS.  `tgauge`: 38/38 PASS.**

> ⚠️ **`OMP_NUM_THREADS` matters on this machine.**  It has 16 logical CPUs and a load
> average of 25–65 from other users.  QEX's `threadBarrier` is a *spin* wait, so the
> default 16 OpenMP threads livelock and `tgauge` appears to hang (I saw it stall for
> >60 s at random points).  With `OMP_NUM_THREADS` = 1, 2, 4 or 8 it finishes in ~0.1 s,
> 12/12 runs clean.  Results are **bit-identical** for all four thread counts.
> This is environmental, not a bug in the honeycomb code.

### hcgeom.nim: no changes needed

`hcgeom.nim` passed every check unmodified.  Two items that had *not* previously been
verified numerically are now verified:

* **FORMULATION (4.2), the `3/8`, is correct.**  `fromOmega(projectOmega(F)) = F` for
  random real antisymmetric `F` to `4.44e-16`, and the "gain" on a pure `F_01` is
  `1.000000000000000`.  (Analytically: for `(α,β)`, the 4 hexagons with axis `α` give
  `(4/3)F_αβ` and the 4 with axis `β` give another `(4/3)F_αβ`; `(3/8)(8/3) = 1`.)
* **`triPath` matches the hand-derived forms of FORMULATION §2.2 exactly**, as
  `LinkRef` sequences, for all 32 apex triangles on both sublattices.

### tgeom.nim numbers (all 9 PLAN items)

| check | result |
|---|---|
| 24 directions, `dot2(n,n) = 4`, closed under negation | ✅ |
| `Σ_i n_μ n_ν` over 24 / over 8 axis | `6 δ_μν` / `2 δ_μν` ✅ |
| unordered zero-sum triples | **32**, each 1 axis + 2 diagonal ✅ |
| triangles through a site (pairs at 60°) | **96** (192 oriented) ✅ |
| of those with the apex at the site | **32**, = `apexTris`, no duplicates ✅ |
| triangles per cell / links per cell (4⁴ cells) | 16384/256 = **64**, 6144 links ✅ |
| triangles containing a given link | min 8, max 8 → **exactly 8** ✅ |
| `triPath` closes, links chain, Σ displacement = 0 | ✅ (64 paths) |
| `triPath` == FORMULATION §2.2 forms | ✅ both sublattices |
| hexagons | **16**, unit vectors, consecutive `n·n' = 1/2`, opposite `−1`, coplanar ✅ |
| hexagon edge-pairs | **96** distinct, each in exactly one hexagon ✅ |
| apex triangles per hexagon | **2**, total `16×2 = 32` ✅ |
| `hexTriPaths` | 96 loops close, links chain, all 6 same orientation ✅ |
| clover reconstruction | `4.44e-16` ✅ |
| point group | **1152** (24 dirs) / **384** (8 axis) ✅ |

### tgauge.nim numbers

| check | result |
|---|---|
| `nCells`/`nSites`/`nLinks` on 4³×6 cells | 384 / 768 / 9216 ✅ |
| `HcShift16` forward vs brute force, SIMD (`V=4`, innerGeom `[1,1,2,2]`) | **0.0** (bit exact), all 16 δ |
| `HcShift16` backward, same | **0.0** |
| same on a `V=1` layout | **0.0** / **0.0** |
| same on `[2,4,4,6]` (L=2 ⇒ shifts wrap) | **0.0** / **0.0** |
| `run` after the source changed / `run(s, src)` rebind | **0.0** |
| unit gauge `triangleSum` | **1.0000000000000000** exactly, Im = 0 |
| field `triangleSum` vs brute-force `hcgeom.triPath` reference (random SU(3)) | `0.00181853981285684` vs `0.00181853981285684`, diff **2.2e-18** |
| same, `V=1` layout | diff **1.4e-18** |
| same, `[2,4,4,6]` | diff **1.5e-17** |
| **gauge invariance** under random SU(3) `vA`,`vB` | `2.8008866973819774e-4` → `2.8008866973818999e-4`, diff **8.7e-19** |
| ... and `Σ_links‖U'−U‖² = 5.5e4` (min per field 2.2e3) so the transform is non-trivial | ✅ |
| gauge-transformed **unit** config | `0.99999999999999134` (8.7e-15 from 1) |
| `triangleSum` of 8 random configs | mean `−3.8e-4`, max `|Re|` `3.2e-3`, vs `1σ = 1.5e-3` |
| every link in exactly 8 triangles (4⁴ cell torus, via `triPath`) | ✅ |
| `warm(0.3)` → `triangleSum` | `0.694132` |
| `checkSU` after `random`/`gaugeTransform` | max `9.6e-12` / `2.0e-11` (projectSU accuracy) |

The **brute-force reference** (`triangleSumRef` in `tests/tgauge.nim`) shares no code
with the field implementation: it walks `hcgeom.triPath` cell by cell, fetches every
link with `lo.coord`/`lo.rankIndex` single-site indexing and multiplies hand-written
3×3 complex matrices.  Agreement at the `1e-18` level is the strongest check in this
task, and it also validates the re-basing trick described below.

### API settled (task A builds on this)

```nim
# hclayout.nim
type HcLayout*[V: static[int]] = ref object
  lo*: Layout[V]; ns*, nt*: int; geom*: seq[int]
proc newHcLayoutX*(geom: openArray[int], V: static[int]): HcLayout[V]
template newHcLayout*(geom: openArray[int])                    # V = VLEN
template newHcLayout*(geom: openArray[int]; V: static[int])
template newHcLayout*(ns, nt: int)
template newHcLayout*(ns, nt: int; V: static[int])
template nCells*(hl: HcLayout): int        # global
template nSites*(hl: HcLayout): int        # 2*nCells, global
template nLinks*(hl: HcLayout): int        # 24*nCells, global

type HcShift16*[F, S] = object             # f[delta] = src(y + sign*delta)
  f*: array[16, F]; sh*: array[16, S]; sign*: int
proc newHcShift16*[F](src: F; sign: int = 1): auto   # outside threads:
proc run*(s: var HcShift16)                          # inside  threads:
proc setSrc*[F,S](s: var HcShift16[F,S], src: F)     # outside threads:
proc run*[F,S](s: var HcShift16[F,S], src: F)        # outside threads:
template `[]`*(s: HcShift16, delta: int)
func topBit*(delta: int): int

# hcgauge.nim
type HcGauge*[V: static[int], F] = object
  hl*: HcLayout[V]
  uA*, uB*: array[4, F]     # A(y)->A(y+e_mu), B(y)->B(y+e_mu)
  uD*: array[16, F]         # B(y)->A(y+delta)
proc newHcGauge*[V: static[int]](hl: HcLayout[V],
                                 nc: static[int] = getDefaultNc()): auto
proc newOneOf*[V,F](g: HcGauge[V,F]): HcGauge[V,F]
template lo*(g: HcGauge)
template eachLink*(g: HcGauge, u, body: untyped)   # allocation free
proc allLinks*[V,F](g: HcGauge[V,F]): seq[F]       # ALLOCATES: outside threads:
proc link*[V,F](g: HcGauge[V,F], k: LinkKind, idx: int): F
proc unit*(g: HcGauge)
proc random*(g: HcGauge, r: var RNGField)          # inside threads:
proc warm*(g: HcGauge, s: float, r: var RNGField)  # inside threads:
proc reunit*(g: HcGauge)                           # inside threads:
proc checkSU*(g: HcGauge): tuple[avg, max: float]  # inside threads:
proc `:=`*(a: HcGauge, b: HcGauge)                 # inside threads:
proc triangleTrace*[V,F](g: HcGauge[V,F]): tuple[re, im: float]  # opens threads:
proc triangleSum*(g: HcGauge): float                             # opens threads:
proc triangleSumIm*(g: HcGauge): float                           # opens threads:
proc gaugeTransform*[V,F](g: HcGauge[V,F], vA, vB: F)            # opens threads:
```

**Deviations from the interface sketched in PLAN.md §Task L / the task brief**

1. `HcGauge` and `HcShift16` each take **two** generic parameters, not one.  Nim cannot
   recover `V` from `Field[V,T]` in a field declaration, nor `Shifter[F,T]` from `F`.
   Callers use `auto`/`type(...)` so it is invisible in practice.
2. The mutators (`unit`, `random`, `warm`, `reunit`, `gaugeTransform`) take `g` by
   **value, not `var`** — the 24 fields are refs, so their contents are mutated anyway,
   and this matches QEX's own `random*[F:Field](g: openArray[F], ...)`.
3. Extra: `triangleTrace` (both parts at once), `triangleSumIm`, `newOneOf`, `:=`,
   `checkSU`, `eachLink`, `allLinks`, `link`, `setSrc`.

### Two things worth knowing

**(a) The A-apex loop is evaluated re-based, so it needs only one single-axis shift.**
`triangleTrace` computes the apex-B loops literally as
`uD[δ](y)·uA[μ](y+δ)·uD[δ'](y)†` — this is what `HcShift16` is for.  For apex A,
FORMULATION §2.2 gives `uD[db](y−db)†·uB[μ](y−db)·uD[db'](y−db')`, which would need
*three* diagonal shifts.  Since the trace is cyclic and `y ↦ z = y−db` is a bijection of
the torus, the same set of triangles is obtained from
`uB[μ](z)·uD[db'](z+e_μ)·uD[db](z)†` — one single-axis `+e_μ` shift.  Verified to
`1e-18` against the literal reference.  **Task A must not blindly copy this for the
force**: it is only valid for quantities summed over the whole lattice.

**(b) `triangleSumIm` does *not* vanish for a general configuration.**  The task brief
expected "small and real".  It is exactly 0 for unit and pure-gauge configurations, but
because the 32 apex triangles are enumerated with a **fixed orientation** there is no
reversed partner to cancel against (unlike the cubic plaquette sum, where `Σ Im Tr U` is
zero by pairing).  For random SU(3) it is the same statistical size as the real part
(`max|Im| = 3.5e-3` vs `1σ = 1.5e-3`).  Only `Re` enters the action, so this is
harmless — but do not use `Im ≈ 0` as a correctness check.

### Bug found and fixed while testing (in my own new code)

`hcgauge`'s helpers first built the 24-field list with `allLinks` (a `newSeq`) *inside*
the `threads:` block.  QEX runs `threads:` bodies on raw OpenMP threads, which have no
initialised Nim (`--mm:refc`) thread-local allocator, so this segfaulted/aborted at
random — roughly 1 run in 8.  Fixed by adding the allocation-free `eachLink` template
and rewriting `unit`/`random`/`warm`/`reunit`/`checkSU` on top of it.  I verified with
`awk` over the generated C that **no `newSeq`/`asgnRef`/`unsureAsgnRef` call remains
inside any `_Pragma("omp parallel")` region** of the whole `tgauge` binary.
*Rule for everyone downstream: never allocate inside a QEX `threads:` block.*

### Open issues / deferred

1. **`triangleSum` rebuilds all its shifters on every call** (4×15 tree shifts + 4
   single-axis shifters + 1 accumulator field, ≈65 `ColorMatrix` fields).  Fine for
   tests, hopeless for HMC.  Task A should hoist an `HcShift16` per `uA[μ]` (and
   whatever the staples need) into a persistent work object.  `HcShift16` is already
   designed for that: construct once, call `run` per evaluation.
2. **`triangleSum` accumulates into a `ColorMatrix` field and traces once** (the
   `plaq2` style, 2 matmuls + 1 matrix add per triangle).  A `redot`-based version would
   save the add; deliberately not done.
3. **Multi-rank is untested.**  The field implementation goes through ordinary QEX
   shifters and `threadRankSum`, so it should be MPI-correct, but the brute-force
   reference in `tests/tgauge.nim` is single-rank (`doAssert ri.rank == lo.myRank`) and
   the tests were only run on 1 rank.
4. **Nc ≠ 3 untested.**  `newHcGauge(hl, nc)` takes a static `Nc` and `random`/`warm`/
   `reunit` branch on `nrows == 1` (U(1)) vs SU(N), but only the default `Nc = 3` was
   exercised.
5. `HcLayout.ns` is just `geom[0]`; for a non-cubic spatial geometry use `geom`.

---

## 2026-08-21 — Task F, free-field fermions (slides 12, 13, 14)

**Files added**

* `hcfree.nim` — free momentum-space Wilson-type Dirac operator for both lattices:
  `freeM`/`freeK` (direct sum over the 8 / 24 neighbour vectors), `freeMClosed`/
  `freeKClosed` (the closed forms of FORMULATION §5.3), `freeEigs` (exact,
  `m+M ± i|K|`), the 4×4 matrix forms `freeD4`/`freeD4Direct`, the 8×8
  cell-momentum form `freeD8` + `gamma5x8`, DeGrand–Rossi `gammaMat`/`gamma5Mat`,
  small dense complex linear algebra (`mmul`/`adj`/`tr`/`det`), the momentum
  iterators `momenta`/`momentaAlt`/`cellMomenta`/`branchMomenta`, and the
  transfer-matrix machinery `polyB`/`cubicRoots`/`modeV`/`pressureIntegrand`.
  Dependencies: `std/[math, complex]` + `hcgeom` only.
* `hcFreeSpectrum.nim` — slide 14 / Fig. 3.
* `hcFreePressure.nim` — slides 12–13 / Fig. 4.
* `tests/tfree.nim` — 19 `unittest` cases, all passing, exit code 0.
* `doc/plots/freespec_{cubic,16cell}.dat`, `freespec.gnuplot`, `freespec.png`
* `doc/plots/pressure.dat`, `pressure_conv.dat`, `pressure_series.dat`,
  `pressure.gnuplot`, `pressure.png`

**Build / run**

```bash
export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
cd build_mac
make run src/experimental/honeycomb/tests/tfree.nim        # 0.2 s
make run src/experimental/honeycomb/hcFreeSpectrum.nim     # 1 s   -> freespec.png
make run src/experimental/honeycomb/hcFreePressure.nim     # 25 s  -> pressure.png
# optional arguments: hcFreeSpectrum [ns] [nt] [per|anti]
#                     hcFreePressure [ntMax] [maxGrid]      (defaults 32, 1024)
```

### Extreme eigenvalues (slide 14) — all ✅

| quantity | 16-cell | cubic |
|---|---|---|
| max `Re λ` | **5.333333333333** = 16/3 | **8.000000000000** |
| max `\|Im λ\|` | **1.467889825014** | **2.000000000000** |
| small-`p` `\|K\|/\|p\|` | 0.999999999 | 0.999999998 |
| small-`p` `M/(r p²/2)` | 0.999999994 | 0.999999994 |

Found by a 16⁴ grid scan over the full period box plus a shrinking local search.
`max Re λ` sits at `p=(π,π,π,π)` and equivalently at `(2π,0,0,0)`.

**New closed form.**  The 16-cell `max |Im λ|` is exactly

```
    3^(1/4) (1 + √3) / √6  =  1.46788982501387…
```

attained at `p = x·(1,1,1,1)/2` and at `p = (x,0,0,0)` with
`x = 2 arccos((√3−1)/2) = 2.392123788172`, i.e. `x/2 = 1.196062`.

### Free spectrum, 16⁴, r=1, m=0 (slide 14 = Fig. 3)

`doc/plots/freespec.png`.  Periodic in time by default (that is what makes the
cubic spectrum reach exactly `Re λ = 8`, as on the slide; `anti` is available).
131072 cubic modes (981 distinct to 1e-6) and 262144 16-cell modes (1196
distinct); the point group collapses the rest, the scatter plot is identical.
The blue 16-cell cloud is the compact near-ellipse `Re λ ∈ [0, 5.333]`,
`|Im λ| ≤ 1.4679`, **including the two small holes on the real axis near
`Re λ ≈ 4.3` and `4.7`** that the slide shows; red cubic fills `[0,8]×[−2,2]`
in the familiar three-lobed pattern.  Reproduced.

### Free pressure (slides 12–13 = Fig. 4)

**Method — the `ln p²` singularity is removed analytically, not numerically.**
At fixed spatial `q`, `M²+|K|²` is a polynomial `P(x)` in `x = cos θ`, with
`θ = p₃` (cubic, degree 1) or `θ = p₃/2` (16-cell, degree **3**; the `x⁴`
coefficient cancels identically).  With `x_i = cosh E_i`, `v_i = e^{−E_i}`,
`N' = N_t` / `2N_t`, the Matsubara product identity gives *exactly*

```
  O(N_t) = N_t³ · < Σ_i ln|1 + v_i(q)^{N'}| >_q ,   <·>_q = ∫d³q/(2π)³ over [0,2π)³
```

for **both** lattices.  The `T=0` subtraction and its logarithmic singularity
have cancelled; what is left is analytic on the torus except for the `|q|` cone
at the origin.  Because the odd part of `ln(1+e^{−x})` is exactly `−x/2`, the
uniform-grid error is a clean series `a/n_g⁴ + b/n_g⁶ + …`, killed by
Richardson extrapolation.  Two independent **geometric** grid families are used,
`A = 32,64,128,256,512,1024` and `B = 48,96,192,384,768` (ratio exactly 2 —
the Richardson recursion only eliminates the terms exactly for a constant grid
ratio); each is extrapolated on its own and their difference is the quoted
error.  The `n_g³` sums are 48-fold symmetry-reduced and the reduction is
verified against a brute-force sum (`1e-15`).  Result: `O/O_cont` is accurate
to `1e-13 … 1e-10` for `N_t ≤ 26`, growing to `≈1e-9` at `N_t = 32` (see
`pressure_conv.dat` and the "grid err" columns printed by the program).
Runtime 25 s.

For the 16-cell the BZ doubling is put in the **time** direction
(`p₃ ∈ [0,4π)`, `2N_t` Matsubara points, spatial `q ∈ [0,2π)³`) rather than in
`p₀`.  Both enumerations are complete coset representative sets of the momentum
lattice mod `2π D₄` and are checked to give identical spectral sums.

**`O/O_cont`** (`O_cont = 7π²/720 = 0.095954487233`):

| N_t | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 12 | 16 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| cubic | 4.9498 | **3.82813** | **2.80138** | **2.11549** | **1.70406** | **1.46437** | **1.32309** | **1.23683** | **1.14443** | 1.07266 | **1.04454** |
| 16-cell | 1.2279 | **1.07191** | 1.02584 | 1.01124 | 1.00573 | 1.00325 | 1.00200 | 1.00130 | 1.00062 | 1.00019 | 1.00008 |

Every acceptance target is hit: cubic 3.83 / 2.80 / 2.11 / 1.70 / 1.47 / 1.32 /
1.23 / 1.14 / 1.04, 16-cell ≈1.07 at N_t=4 and ≈1.00 from N_t≳6.  The paper's
statement "at N_t = 4 the correction is 7 % and 283 %" comes out as
**7.19 %** and **282.81 %**.

### Asymptotic coefficients — the sharp validation ✅

Fitting `O/O_cont − 1 = Σ_k c_{2k}(π²/N_t²)^k` (Neville extrapolation in
`1/N_t²` and least squares over `N_t = 12..32` with 4/5/6 terms; the quoted
error is the spread over those estimators):

| coefficient | this work | published | deviation |
|---|---|---|---|
| cubic `c₂` | **1.686996 ± 0.001521** | 248/147 = 1.6870748 | **0.005 %** |
| cubic `c₄` | **4.3499 ± 0.2073** (Neville alone: 4.31499) | 635/147 = 4.3197279 | 0.70 % (0.11 %) |
| 16-cell `c₄` | **0.1295918 ± 0.0000050** | 127/980 = 0.1295918 | **0.0000 %** |
| 16-cell `c₆` | **0.017787 ± 0.000052** | 73/4158 = 0.0175565 | 1.3 % |
| **16-cell `c₂`** | **(−2.1 ± 2.8)·10⁻⁸** | **0** | consistent with zero |

The 16-cell `1/N_t²` coefficient is **eight orders of magnitude smaller** than the
cubic one (`2e-8` vs `1.687`) — slide 12's headline claim (the leading 16-cell
correction is `O(a⁴)`) is confirmed at the `10⁻⁸` level.  It stays at that level
across all fit windows and term counts (`Nt = 6..32` / `8..32` / `12..32`,
4/5/6 terms).

Independent, assumption-free check (section C of the program output): with the
*published* coefficients inserted,
`R(N_t) = [O/O_cont − 1 − Σ published terms]·N_t^p/π^p` must tend to a finite
constant.  It does — cubic (p=6) settles towards ≈16.9 (so `c₆ ≈ 17`), and
16-cell (p=8) plateaus at **0.1896** around `N_t = 24` (so `c₈ ≈ 0.19`, matching
the fitted `c₈ ≈ 0.185`) before the `1e-9`-level grid noise, amplified by
`N_t⁸/π⁸ ≈ 10⁸`, takes over beyond `N_t ≈ 27`.  Had `c₂(16-cell)=0`, `c₄` or
`c₆` been wrong, `R` would have diverged as a power of `N_t`; it does not.

### Two errors found in FORMULATION.md §5.3

1. **BZ volume.**  "the site set `Z⁴ ∪ (Z+½)⁴` has reciprocal lattice `2π·D₄`, so
   the BZ has volume `(2π)⁴/2`" — the reciprocal lattice is right but the volume
   is wrong.  `D₄` has covolume 2, so `2πD₄` has covolume `2(2π)⁴`: the 16-cell
   BZ is **twice**, not half, the cubic one.  That is what makes `2N_s³N_t`
   momenta fit `2N_s³N_t` sites.  (The momentum counting elsewhere in the
   document is consistent with the correct value.)
2. **Location of `max |Im λ|`.**  The value 1.46789 is right and the
   `p ∝ (1,1,1,1)` statement is right, but the second quoted point
   `p = (2π/3,0,0,0)` is **not** a maximum: `|K|(2π/3,0,0,0) = 1.4433757`.
   The correct axis-direction maximum is at `p = (2.392124,0,0,0)`, i.e.
   `2 arccos((√3−1)/2)`, not `2π/3`.  Both are pinned in `tests/tfree.nim`.

Nothing else in FORMULATION §5 disagreed with the numerics.

### What `tests/tfree.nim` checks (19 cases, all pass)

Clifford algebra and `γ₅ = γ₀γ₁γ₂γ₃`; direct sum vs closed form for `M` and `K`
(1e-15) over 20000 random momenta and both lattices; `Σᵢ n_μ n_ν = 6δ / 2δ`;
`freeD4` vs the position-space assembly `Σᵢ(γ·nᵢ−r)e^{ip·nᵢ}`;
`det(D−λ)=0` for `λ = m+M±i|K|`; the small-`p` limits; the four extreme
eigenvalues plus the closed form; `γ₅ D₈ γ₅ = D₈†`; **8×8 cell-momentum spectrum
= union of the two true-momentum branches** (multiset equality via
`tr D⁸ⁿ, n=1..8`, plus `det(D₈−λ)=0`); momentum counts `= N_sites`; the two
16-cell momentum enumerations agree; exactly one zero mode; `P(x)` reproduces
`M²+|K|²`; the Matsubara identity per `q` (1e-14); an **end-to-end** check that a
direct `Σ_p ln det D(p)` momentum sum reproduces `O = N_t³⟨h⟩` including every
factor of 2 (1e-13); a fixed-grid pinning value of the pressure ratio; and
`ω/|q| → 1`.

### Open issues

* Nothing blocking.  `hcfree.nim` is self-contained (no QEX `Layout`), so task D
  can use `freeEigs`/`freeD8` directly as the `U = 1` reference for the
  interacting operator (PLAN task D1).
* The cubic `c₄` is the least well determined coefficient (0.2–0.7 %), because
  the cubic series has a large `c₆ ≈ 17`.  Pushing `N_t` beyond 32 would need
  `n_g > 1024` and is not worth it.
* `hcFreeSpectrum` defaults to **periodic** time; the paper's figure appears to
  be periodic too (its cubic band reaches exactly `Re λ = 8`).  Pass `anti` for
  the thermal spectrum.

---

## 2026-08-21 — Task C: cubic reference pipeline (flow → t₀ → Q → χ_top → continuum)

Full write-up with all numbers: **[RESULTS_CUBIC.md](RESULTS_CUBIC.md)**.

**Owns / added:** `hcanalysis.nim`, `refCubicGen.nim`, `refCubicMeas.nim`,
`tests/tanalysis.nim`, `doc/plots/cubic/*`, `doc/RESULTS_CUBIC.md`.
Nothing else was touched.

### ⚠️ The `topoQ` factor-2 question is settled: **QEX is right, FORMULATION §4.3 and PLAN task W2 are wrong**

`gaugeUtils.topoQ`'s prefactor `−1/(4π²)` on `(a−b+c)` is **correct**.
`hcanalysis.qexTopoQNormFix = 1.0`; the honeycomb clover charge must use the
same normalisation, i.e.

```
  q(x) = −(1/32π²) ε_{μνρσ} tr(F̂_μν F̂_ρσ)      F̂ traceless ANTI-hermitian
```

(FORMULATION §4.3 has the right magnitude `1/32π²`; only the sign convention
differs, and its "suspect by a factor 2" remark should be deleted.
**Task W: do not double the charge.**)

The suspicion came from pairing `Q = (1/32π²)∫εF^aF^a` with
`F^aF^a = −2 tr(FF)`. The first is a mis-remembered formula: the textbook
identity is `Q = (1/32π²)∫F^a F̃^a = (1/64π²)∫εF^aF^a` — the `εFF` form carries
`1/64π²`. (BPST check: `∫F^aF^a = 32π²` for a self-dual `Q=1` instanton, so
`1/32π²` on `εFF` would give 2.)

**Exact numerical proof, not statistical**: `refCubicMeas -abeliantest` builds a
constant-field-strength Abelian SU(3) configuration in the Cartan direction
`T = diag(1,−1,0)` with fluxes `n₁, n₂`, whose exact charge is
`Q = Σ_i q_i² n₁n₂ = 2n₁n₂` (Atiyah–Singer for a direct sum of U(1) bundles).

```
./bin/refCubicMeas -abeliantest:1 -lat:16,16,16,16 -n1:1 -n2:1
  exact 2      QEX topoQ 1.999598437023   ratio 0.99979922
```

The `1 − (φ₁²+φ₂²)/6` clover artefact is reproduced to all printed digits and
scales as `1/L⁴` (ratio 0.99679 / 0.99937 / 0.99980 / 0.99996 / 0.99999 at
L = 8/12/16/24/32); `densityE = 2sin²φ₁+2sin²φ₂` to 2e−16.
Table: `doc/plots/cubic/abeliantest.dat`.

**Warning for task W:** "Q is near an integer" is *not* a usable normalisation
test at coarse `a`. With the plain clover at `a²/t₀ ≳ 0.5`, `Q(t₀)` misses the
nearest integer by 0.2–0.45 routinely; at `a²/t₀ ≈ 1.5` it does not even
plateau under flow. Use the Abelian construction instead — it is exact, takes
0.3 s, and transfers directly to the honeycomb (a Cartan-valued
constant-field-strength configuration is equally easy to write down there).

### Two more things that will bite the honeycomb tasks

1. **`t₀` does not exist if `L/√t₀ ≲ 6`.** On 8⁴ at β=6.0 (`L/√t₀ = 4.5`)
   `t²⟨E⟩` rises only to 0.127 and then *falls* — the flow radius reaches `L/2`
   and the field is smoothed to a constant, so `t²⟨E⟩ = 0.3` is never reached.
   PLAN's implicit "`L/√t₀ ≳ 4`" is far too loose; **use `L/√t₀ ≳ 9–10`**.
2. **QEX's SIMD layout rejects `L = 10` and `L = 14`** (`vlen=4` splits two
   directions by 2 and then needs the outer geometry even):
   `Error: can't lay out inner geom`. Use `L ∈ {8,12,16,20,24}`.
3. Confirming the STATUS note above about threads, now quantified for HMC:
   one `refCubicGen` at 8⁴ takes 1.48 s (1 thread), 0.53 s (4), 0.54 s (8),
   **17.5 s (16)**. Never give a QEX process 16 threads on this box.

### Build / run

```bash
export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
cd build_mac
make src/experimental/honeycomb/tests/tanalysis.nim && ./bin/tanalysis   # 23/23 PASS
make src/experimental/honeycomb/refCubicGen.nim
make src/experimental/honeycomb/refCubicMeas.nim
./bin/refCubicMeas -abeliantest:1 -lat:16,16,16,16 -n1:1 -n2:1
P=../src/experimental/honeycomb/doc/plots/cubic
sh $P/run_gen.sh $TMPDIR/qexref;  sh $P/run_meas.sh $TMPDIR/qexref 6
python3 $P/harvest.py $TMPDIR/qexref --skip 3 --frozen 6.0 > $P/chitop.dat
sh $P/mkplot.sh
```

### HMC (acceptance criterion: `⟨exp(−ΔH)⟩ = 1`, acceptance > 70 %) — passed

Wilson plaquette action, 4MN5FV integrator, `τ=1`, `nsteps=4` (20 force
evaluations/trajectory), hot start.

| β, L | trajs | acceptance | `⟨exp(−ΔH)⟩` | `⟨plaq⟩` | literature |
|---|---|---|---|---|---|
| 5.6, 8⁴  | 6000 | 98.7 % | 0.99997(48)  | 0.524163(42) | 0.5242 |
| 5.7, 12⁴ | 4000 | 96.6 % | 1.00053(140) | 0.549085(21) | 0.5493 |
| 5.8, 12⁴ | 4000 | 96.1 % | 0.99826(152) | 0.567602(17) | 0.5670 |
| 5.9, 16⁴ | 1379 | 92.6 % | 1.00182(476) | 0.581855(15) | 0.5822 |
| 6.0, 16⁴ | 1537 | 92.6 % | 1.00175(499) | 0.593703(14) | 0.5937 |

Reversibility `|ΔH_fwd+ΔH_bwd| ≤ 5e−10` on every check.
`4MN5FV nsteps=4` beats `2MN nsteps=12` (100 % vs 55 % acceptance at 12⁴,
β=5.8, at lower cost) — worth copying for `hchmc.nim` (task M).

### Slide 10, cubic curve

`t₀/a²` measured: 0.6689(27), 0.9860(40), 1.5176(103), 2.3089(227), 3.3113(485)
at β = 5.6…6.0 — 2–6 % from the Necco–Sommer `r₀/a` + `t₀/r₀²=0.1107`
expectation, an independent validation of the flow and of the `t²⟨E⟩=0.3`
convention.

| β | L | n | a²/t₀ | 10⁴t₀²χ |
|---|---|---|---|---|
| 5.6 | 8  | 197 | 1.4951(59) | 1.902(201) |
| 5.7 | 12 | 133 | 1.0142(41) | 3.013(384) |
| 5.8 | 12 | 159 | 0.6589(45) | 4.079(496) |
| 5.9 | 16 | 55  | 0.4331(43) | 5.574(1112) |
| 6.0 | 16 | 59  | 0.3020(44) | 3.573(935) ← topology frozen, excluded |

O(a²) continuum extrapolation (β=5.6…5.9): **5.91 ± 0.68, χ²/dof = 0.38**;
O(a⁴): **7.65 ± 2.41, χ²/dof = 0.19**. Reference (Cè et al. 1506.06052):
**6.67 ± 0.07**. So the O(a²) fit lands 1.1 σ low and the O(a⁴) fit 1.4 σ high,
both consistent with 6.7 at the 12–31 % precision we reach.
The reason the O(a²) fit undershoots is that our `a²/t₀ ∈ [0.30, 1.50]` is
2–4× coarser than slide 10's red points (`0.15…0.35`), where the cut-off
effect is a factor 3.5 — exactly the regime the paper says the cubic lattice
cannot be extrapolated from. Plot: `doc/plots/cubic/chitop.png`.

**Honest statistics statement.** τ_int(Q) ≈ 17–22 trajectories at β ≤ 5.8
(60–170 effectively independent configurations, errors trustworthy); ≈ 24 at
β=5.9 (~14 independent, and ⟨Q⟩ = −0.48 is still 1.5 σ from 0 — marginal); at β=6.0 (a≈0.09 fm) the charge does a slow
random walk with `|Q| ≲ 2` over the whole stream, `⟨Q²⟩` was still rising when
the run stopped and the Madras–Sokal window reached a quarter of the sample —
**that point's error is not trustworthy** and it is excluded from the primary
fit.

### Open / handed on

* Task W should reuse `hcanalysis.findT0/findW0/jackknife/fitPoly` verbatim and
  `refCubicMeas`'s `-abeliantest` idea for the honeycomb clover normalisation.
* FORMULATION §4.3 and PLAN task W2 need their "factor 2" warnings corrected by
  their owners (I did not edit those files).
* Not done for lack of Monte Carlo time: β ≥ 6.0 with enough statistics to
  reach slide 10's `a²/t₀` window. That needs either far more trajectories or
  a topology-friendly algorithm (open boundaries / parallel tempering), or the
  SU(3) heatbath+overrelaxation of PLAN task M2 — which QEX still does not have.

---

## 2026-08-21 — Task A: triangle gauge action, staple derivative, force

**Files added** (nothing else touched):

| file | contents |
|---|---|
| `hcaction.nim` | staple recipe (auto-derived from `hcgeom.triPath`), `HcActionWork`, `hcAction`, `hcActionDeriv`, `hcForce`, `redot(HcGauge,HcGauge)`, `randomTAH(HcGauge,…)` |
| `tests/taction.nim` | 6 suites / 26 PASS: QEX-convention pin, FD force checks, invariance/covariance, continuum limit |

**Build / run**

```bash
export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
cd build_mac
make src/experimental/honeycomb/tests/taction.nim && OMP_NUM_THREADS=4 ./bin/taction
make src/experimental/honeycomb/hcaction.nim && OMP_NUM_THREADS=4 ./bin/hcaction  # smoke + timing
```

26/26 PASS, exit 0; identical results at `OMP_NUM_THREADS` = 1 and 4; three
repeat runs clean.  Verified (awk over the generated C, Task-L style) that no
`newSeq`/`asgnRef`/`unsureAsgnRef`/`rawNewObj` appears inside any
`omp parallel` region of `taction`/`hcaction` (16 + 8 regions scanned).

### Public API (task M consumes this)

```nim
type HcActionWork*[V:static[int],F,SH,SS] = ref object   # all shifters + scratch
proc newHcActionWork*[V,F](g: HcGauge[V,F]): auto        # outside threads:, once
proc hcAction*(w: HcActionWork, beta: float, g: HcGauge): float
proc hcAction*(beta: float, g: HcGauge): float           # convenience, allocates
proc hcActionDeriv*(w, beta, g, f: var HcGauge)          # staple sums, NOT TAH
proc hcActionDeriv*(beta, g, f: var HcGauge)             # convenience, allocates
proc hcForce*(w, beta, g, f: var HcGauge)                # the MD force
proc hcForce*(beta, g, f: var HcGauge)                   # convenience, allocates
proc redot*(a, b: HcGauge): float        # sum_l sum_x Re tr(a_l^dag b_l)
proc randomTAH*(g: HcGauge, r: var RNGField)             # inside threads:
```

All procs rebind the work object to the `g` they are given, so one work object
serves any number of same-shape gauges (that is what the FD test does).
`f`'s 24 fields must be distinct from `g`'s.  Everything after
`newHcActionWork` is allocation-free per call.

### Force convention (the exact equation; tasks M and W must keep it)

```
  S(U) = (β/2) Σ_x Σ_{i=1}^{32} ( 1 − (1/N) Re Tr P_i(x) )

  hcActionDeriv:  D_l = (β/2N) Σ_{k=1}^{8} V_k ,   V_k = W_k†,  P_k = U_l W_k
                  (P_k = triangle loop re-based to start with U_l un-daggered,
                   anchored at the cell where the link field lives)
  hcForce:        f_l = projectTAH( U_l D_l† ) = (β/2N) projectTAH( Σ_k U_l W_k )

  d/ds S( exp(s P_l) U_l ) |_{s=0}  =  Σ_l redot( P_l , f_l )        …(★)
```

with `redot(a,b) = Re tr(a†b)` summed over sites (QEX field `redot`), P_l any
traceless anti-Hermitian momenta.  (★) is **the same convention QEX's
`gaugeForce` satisfies** — test 1 verifies that for QEX itself on a cubic
lattice (rel 7.3e−13), then tests 3–4 pin `hcForce` to it.  Hence the
`puregaugehmc.nim` leapfrog is drop-in:
`p −= dt*f;  U := exp(dt*p)*U` conserves `H = S + ½ Σ_l redot(p_l,p_l)`.
Note the **plus** sign in `f_l = +(β/2N) projectTAH(U_l Σ W_k)`: with QEX's
`redot(P,f) = −Re tr(P f)` for anti-Hermitian P, a minus sign here (which a
naive reading of `dS/ds = −(β/2N) Re tr(P U W)` suggests) would violate (★);
the FD tests would have caught it instantly.

### How the staples are generated (and why they can't drift)

At module load `hcaction.nim` walks all 64 triangles per cell via
`hcgeom.triPath`, re-bases every loop at each of its 3 links (reversing
orientation where the link enters daggered), converts to QEX staples
`V = W†`, and `doAssert`s the resulting 192 terms equal the closed-form six
families the executor implements (3 per apex-B triangle, 3 per apex-A
triangle; module docs list them).  Per `apexTris` entry the executor needs the
forward 16-tree `uA[μ](x+δ)` (4 `HcShift16`s), one `+e_μ` shift of `uD[dbp]`,
one backward chain shift by `−δ` of the 2-link product `uD[δ]†uD[δ']`
(≤3 single-axis shifts), and one `−e_μ` shift of `uB[μ]†uD[db]` — 172
single-axis shifts and ~200 field matmuls per deriv call, all preallocated in
`HcActionWork`.  `hcAction` reuses the same recipe through per-pair `redot`s
(each triangle once), an independent code path from Task L's `triangleTrace`.

### Numbers (β = 5.7)

Finite-difference force vs `Σ redot(P,f)` (warm config, `ndiff` ordMax 7):

| momentum | [4,4,4,4] rel | [2,4,4,6] rel |
|---|---|---|
| all 24 fields | **8.8e−12** | **1.6e−11** |
| uA only | 3.3e−11 | — |
| uB only | 6.3e−11 | — |
| uD only | 8.0e−12 | — |
| uD[5] only | 6.9e−12 | 1.7e−10 |
| all, V=1 layout | 1.6e−11 | — |

* `hcAction` vs `(β/2)·32·nSites·(1 − triangleSum)` (Task L's brute-force-
  verified path): rel **3.9e−16** (warm), **0.0** (random).
* Unit gauge: `hcAction = 3.9e−12` absolute (scale `β/2·nTri = 7.0e4`, so
  ~5e−17 relative — pure `nTri − Σ` cancellation noise); force **exactly 0.0**.
* Gauge invariance of `hcAction` under random SU(3) `vA,vB`: rel diff **0.0**
  (bit-identical doubles).  Force covariance `f → V_start f V_start†`:
  rel **1.6e−27**.
* Random config: S = 70058.25 > 0; warm(0.35): 27819.89; ordering sane.

**Classical continuum limit / β-normalisation** (test 6): weak Abelian plane
wave `A_μ = ε_μ cos(p·x)` embedded via `T = diag(1,−1,0)`, links from exact
straight-line integrals on an `Ns⁴`-cell honeycomb and an `Ns⁴` cubic lattice,
same β; cubic side cross-checked against its analytic small-ε value (≤1e−4).

| Ns | p dir | k | ε₀ | S_16cell | S_cubic | ratio | \|ratio−1\|/p² |
|---|---|---|---|---|---|---|---|
| 12 | 0 | 1 | 1e−3 | 7.77156015e−3 | 7.65365943e−3 | 1.01540449 | 0.05619 |
| 12 | 0 | 2 | 1e−3 | 3.03857354e−2 | 2.85638430e−2 | 1.06378317 | 0.05816 |
| 12 | 0 | 1 | 3e−2 | 6.99439002e+0 | 6.88821060e+0 | 1.01541466 | 0.05623 |
| 12 | 3 | 1 | 1e−3 | 7.03459033e−3 | 6.92783749e−3 | 1.01540926 | 0.05621 |
| 8  | 0 | 1 | 1e−3 | 3.42135233e−3 | 3.30514711e−3 | 1.03515887 | 0.05700 |

Ratio → 1 with `O((pa)²)` deviations exactly as FORMULATION §3.1 predicts:
momentum scaling k=1→2 gives dev ratio **4.14** (≈4), volume scaling 12→8
gives **2.28** (≈2.25), amplitude dependence 1e−5 at ε₀=3e−2 (O(ε²)), and the
time-direction case equals the spatial one to 5 digits (isotropic artifact).
Measured artifact coefficient: `S_hc/S_cubic = 1 + 0.0562(2)·p² + O(p⁴)`
(relative to the cubic Wilson action at the same β).  **β means the same thing
on both lattices** — slide 9's normalisation confirmed numerically.

### Cost (for task M planning)

8⁴ cells, `OMP_NUM_THREADS=4`, this (loaded) box: `hcForce` **49 ms/call**,
`hcAction` **6.7 ms/call** (`./bin/hcaction` prints these).  A τ=1, 4MN5FV,
nsteps=4 trajectory ⇒ ~20 force calls ≈ 1 s at 8⁴.

### Findings / no FORMULATION changes needed

* FORMULATION §2.2/§3.2 are exactly consistent with `hcgeom.triPath`: the
  startup scan reproduces "axis link: 8 B-apexes at `y−δ`, `δ_μ=0`; diagonal
  link: 4 apexes at each end" verbatim.  No Task L file needed changes.
* The apex-A re-basing trick of `triangleTrace` (STATUS Task L note (a)) is
  *not* used for the force; each staple is attributed at the correct cell and
  the per-kind FD tests confirm the attribution link-kind by link-kind.

### Open issues / deferred

1. **Multi-rank untested** (same status as Task L): every shift is a plain
   single-axis `Shifter`, incl. the ≤3-step backward chains, so it should be
   MPI-correct, but only 1 rank was run.
2. **Nc ≠ 3 untested** (recipe and force are Nc-generic).
3. Performance headroom if HMC needs it: the 32 `+e_μ` shifts of `uD[dbp]`
   could become 4 backward trees of `uB` products; `hcAction` could fuse its
   64 `redot` reductions into one accumulator field.  Correctness-first kept.
4. `hcActionDeriv` always overwrites `f` (no `accumulate` flag yet — task W
   can add one if the flow wants to sum several terms).

---

## 2026-08-21 — Task D1: interacting Wilson–Dirac operator (`hcwilson.nim`)

**Files added** (only these two; plus this entry)

| file | contents |
|---|---|
| `hcwilson.nim` | `HcFermion`, `HcWilson`, `D`/`Ddag`, `setBC`, `applyGamma5`, helpers |
| `tests/twilson.nim` | 22 checks in 10 `unittest` cases, all passing, exit 0 |

**Build / run**

```bash
export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
cd build_mac
make src/experimental/honeycomb/tests/twilson.nim && OMP_NUM_THREADS=4 ./bin/twilson
make src/experimental/honeycomb/hcwilson.nim && ./bin/hcwilson   # smoke only
```

12/12 repeat runs clean at `OMP_NUM_THREADS` = 1, 2, 4, 8; all physics numbers
**bit-identical** across thread counts.

### What was built

`D psi(x) = (m + 4 r) psi(x) + (1/6) Σ_{i=1}^{24} (γ·n_i − r) U_i(x) psi(x+n_i)`
(FORMULATION §5.1, `a = 1`), full 4×4 `(γ·n − r)` per hop, **no half-spinor
projection** (deferred by design, PLAN §1.2).

```nim
# hcwilson.nim
type HcFermion*[F] = object
  a*, b*: F                        # one DiracFermion field per sublattice
type SpinMat* = typeof(gamma0)     # 4x4 complex constant spin matrix
proc newDiracField*(hl: HcLayout, nc: static[int] = getDefaultNc()): auto
   # any V incl. V=1 (built like ColorMatrix(l,n) from lo.newDComplexV;
   # lo.DiracFermion() itself is VLEN-only)
proc newHcFermion*(hl: HcLayout): auto            # ALLOCATES: outside threads:
proc newOneOf*(x: HcFermion): HcFermion           # ALLOCATES
proc `:=`*(r: HcFermion, x: HcFermion|SomeNumber) # inside threads: ok
proc gaussian*(x: HcFermion, r: var RNGField)     # inside threads:
proc norm2*(x: HcFermion): float                  # global, both sublattices
proc dot*(x, y: HcFermion): auto                  # conj on FIRST arg (QEX conv.)
proc norm2diff*(x, y: HcFermion): float
proc applyGamma5*(r: HcFermion, x: HcFermion)     # inside threads: ok
proc gammaDotDir*(dir: int): SpinMat              # γ·n for hcgeom dir 0..23
proc hopMat*(dir: int, rw: float): SpinMat        # (γ·n − r)/6

type HcWilson*[...] = ref object                  # ref: threads: can capture it
proc newHcWilson*(g: HcGauge): auto               # ALLOCATES: outside threads:
proc gaugeRefresh*(w: HcWilson)                   # call after links change!
proc D*(w: HcWilson, r: var HcFermion, x: HcFermion, m: float, rw = 1.0)
proc Ddag*(w: HcWilson, r: var HcFermion, x: HcFermion, m: float, rw = 1.0)
proc setBC*(g: HcGauge)                           # antiperiodic time; inside threads:
```

Machinery: 16 QEX `Transporter`s for the axis hops (fused `U*shift`, both
signs); two `HcShift16` trees per application (`psiA(y+δ)` forward,
`psiB(y−δ)` backward, 30 single-axis shifts); and 16 **pre-shifted link
copies** `uDsh[δ](y) = uD[δ](y−δ)` (needed because the A-row diagonal link
lives at the far cell, `uD[δbar](y−δbar)†`), rebuilt by `gaugeRefresh` with 32
chained single-axis backward shifts, **once per gauge configuration**.

The hop tables are **derived from `hcgeom.step` at construction** and every
`LinkRef` field is `doAssert`ed there — a mismatch between the operator and
`hcgeom`'s link convention aborts `newHcWilson`.

`Ddag` is a **direct implementation**, not `γ5 D γ5`: since the γs are
Hermitian and the reverse of hop `n` is hop `−n` through the same link,
`D†` is the same hopping sum with the spin factor of the *opposite* direction
(`gm[opposite(dir)]`).  γ5-hermiticity is therefore a genuine test.

### Conventions the D2/D3/D4 owners need

1. **B-sublattice plane-wave phase.**  `hcfree.freeD8(k)`'s 8×8 block (index
   `2*spin + sub`, sub 0 = A, 1 = B) corresponds to
   `psiA(y) = a e^{i k·y}`, `psiB(y) = b e^{i k·y}` — **both phased with the
   integer CELL coordinate `y`; NO extra half-site phase
   `e^{i k·(½,½,½,½)}` on B.**  (Verified entrywise to 4e-15, tests 2/3/3b.)
2. **QEX `dot` conjugates its FIRST argument** (`dot(x, ix) = +i|x|²`,
   established empirically in test 1).  So adjointness reads
   `dot(x, D y) == dot(Ddag x, y)`.
3. **Antiperiodic time** (`setBC`, sign flip on `uA[3]`, `uB[3]`, and the 8
   `uD[δ]` with bit 3 set, at `y3 = Nt−1`): the operator on the cell-phase
   plane wave at `k` equals `freeD8` at `k3 → k3 + π/Nt` **exactly** (test 7,
   4.4e-15).  Call `setBC` before `newHcWilson`, or `gaugeRefresh` after.
4. **`D`/`Ddag`/`gaugeRefresh` open their own `threads:` block — do NOT call
   them from inside one** (unlike `wilsonD`'s conventions).  `r` must not
   alias `x` (asserted).  `rw` is a per-call parameter; the 24 hop matrices
   are cached and rebuilt only when it changes.
5. `applyGamma5` gives `γ5 ⊗ 1_sublattice` for chirality measurements.

### Test results (`twilson`: 22/22 PASS)

| check | result |
|---|---|
| QEX dot convention: `dot(x, ix) = +i\|x\|²` | ✅ conj-first |
| spinOld `gamma1..4` == `hcfree.gammaMat` entrywise | **0.0** |
| `hopMat(dir, r)` == `(γ·n − r)/6` from hcfree, 24 dirs | 2.8e-17 |
| **free field vs `freeD8`, 4⁴, 10 momenta × (m, r) ∈ {(0,1),(.1,1),(.3,.8),(−.2,1),(0,.7)}** | worst entry **4.44e-15** |
| `Ddag` block == `freeD8†` entrywise | 4.44e-15 |
| free field, asymmetric `[4,4,2,6]`, 5 momenta | 3.55e-15 |
| free field on a `V = 1` (non-SIMD) layout | 1.60e-14 |
| point source (both sublattices, src at wrap corner `[0,3,0,3]`): support | **exactly 25 sites** (1 + 8 axis + 16 diag) |
| ... every neighbour value == `(γ·n − r)/6` column, colors ≠ 0 vanish | 2.78e-17 |
| **gauge covariance** `D[U^g](Vψ) = V(D[U]ψ)`, random SU(3) + random V, `[4,4,4,6]` | rel **1.83e-14** |
| adjointness `<x, D y> = <Ddag x, y>`, random config | 8.45e-16 |
| **γ5-hermiticity**: `γ5 D γ5 x == Ddag x` as fields | **0.0 (bit-identical)** |
| `<x, γ5 D γ5 y> = <D x, y>` | 5.0e-16 |
| antiperiodic blocks == `freeD8(k3 + π/Nt)`, 4 momenta | 4.44e-15 |
| ... and setBC really changes the periodic-k block | diff 0.67 |
| mass linearity `D(m1) − D(m0) = (m1−m0)·1`, random config | rel 1.06e-15 |

(γ5 D γ5 = Ddag being *bitwise* zero is expected, not suspicious: the γ5
sandwich flips exactly the entries that `opposite(dir)` flips, and the
summation order of the two code paths is identical.)

### Timing (8⁴ cells = 4096, random SU(3), best of 5 batches; box load 36–58)

| `OMP_NUM_THREADS` | 1 | 2 | 4 | 8 |
|---|---|---|---|---|
| `D` per application | 4.8 ms | 3.7 ms | **3.8–4.4 ms** | 2.5 ms |

`gaugeRefresh` (once per configuration): 0.5–1.6 ms.  Budget number for the
D3 eigensolver: **~1 µs per cell per D application at 4 threads**, i.e.
~20 ms per D on 12⁴; a 300-eigenvalue Arnoldi at a few·10³ D's ⇒ tens of
seconds per 12⁴ configuration.  (Timings on this shared box fluctuate ~2×
with load; the earlier 78 ms figure seen once at load 59 was interference.)

### Core-QEX gap worked around (in-module, no core files touched)

`Transporter`'s forward apply needs a fused `mul(DiracFermion, ColorMatrix,
DiracFermion)`; `spinOld.mul` only has the all-`Spin` form, so `wilsonD`-style
DiracFermion transport was never instantiated in QEX.  `hcwilson.nim` adds

```nim
template mul*(r: var Spin, x: Color, y: Spin2) = mul(r[], x, y[])
```

(the color matrix then acts as a `Sca2` on the spin `VectorArray`).  Picked up
by `transporterApply`'s `mixin mul` at instantiation; nothing else in QEX is
affected unless it imports `hcwilson`.

### Open issues / deferred

1. **Callers must not be inside `threads:`** (see item 4 above).  If D3's
   Arnoldi wants an inside-threads operator, split `applyDirac` into
   setSrc/threads-body parts — trivial refactor, not needed yet.
2. **Multi-rank untested** (same status as Tasks L/A).  All comms are plain
   single-axis `Shifter`/`Transporter`/`ShiftB` patterns, so it should be
   MPI-correct; only 1 rank was run.
3. **Nc = 3 only tested.**  `newDiracField(hl, nc)` takes a static `nc`, but
   `newHcWilson` currently pairs it with the gauge field's default-Nc type.
4. Performance headroom, deliberately not taken (correctness first):
   half-spinor projection (§5.2; needs 24 custom projectors), fusing the spin
   matrix into the transporter gather, merging the A/B axis loops (46
   barrier-carrying shifts per application today), and comm/compute overlap
   via `startSB`/`boundarySB` directly.
5. `HcFermion` is two independent QEX fields; nothing enforces that `a` and
   `b` share a layout (constructors always do).

---

## 2026-08-21 — Task D3, non-Hermitian eigensolver (Krylov–Schur Arnoldi)

QEX had no non-Hermitian eigensolver (the `zgeev` binding in
`src/eigens/linalgFuncs.nim` is never called, and `import qex` does not even
link LAPACK unless `-d:lapackLib`).  This task adds one as a standalone,
operator-agnostic module, tested on the existing cubic Wilson operator.

**Files added**

* `hcarnoldi.nim` — restarted Arnoldi, **Krylov–Schur** (thick restart in the
  Schur basis of the projected matrix).  Pure Nim + LAPACK, **no QEX
  dependency**: the operator enters through closures, the vector space
  through six mixin ops (`vcopy/vzero/vscale/vaxpy/vdot/vnorm2`).  Contains
  its own `zgeev_`/`zgees_`/`ztrsen_` prototypes; links
  `-framework Accelerate` **by default on macOS** (override:
  `-d:hcLapackLib="..."`, e.g. `:hcLapackLib=-llapack` through the build
  system; elsewhere the default is `-llapack -lblas`).  `when isMainModule`
  is a QEX-free self-test against dense LAPACK (400×400, all four `which`).
* `tests/tarnoldi.nim` — the QEX test (37 checks, all PASS, exit 0) plus the
  reusable ~20-line QEX-`Field` mixin adapter that any later caller
  (`hcSpectrum`) can copy.

**Public API**

```nim
type ArnoldiOp*[V] = object
  apply*: proc (r: var V; x: V)   ## r = A x
  newVec*: proc (): V             ## fresh zeroed vector
  start*: proc (v: var V)         ## deterministic (seeded) start vector

proc arnoldi*[V](op: ArnoldiOp[V]; nev, ncv: int; tol: float;
                 maxRestarts: int; which: string; verb = 0):
    tuple[vals: seq[Complex64]; vecs: seq[V]; resids: seq[float]; nApply: int]
  ## which: "LM" | "SM" | "LR" | "SR".  resids are DIRECT
  ## |A v - lambda v| / max(|lambda|, floor), one extra apply per pair.
  ## verb>=3 prints basis-orthonormality and factorization-error diagnostics.

# also exported: ZMat (column-major complex dense), newZMat, zeig (zgeev),
# zschur (zgees), ztrsenReorder, refineEigvec, eigOrder, frobNorm
```

Algorithm: CGS2 orthogonalization extended adaptively (DGKS) to ≤4 passes;
`zgees`+`ztrsen` restart keeping `kKeep = nev + (ncv-nev)/2` best Schur
vectors (converged pairs soft-lock: their spike entries are ~0); convergence
by the Arnoldi estimate `beta|y_m| ≤ tol·max(|λ|, eps^(2/3)|G|_F)`; final
Ritz vectors from `zgeev` **plus inverse-iteration refinement** on the small
matrix (`refineEigvec`) — necessary because soft-locking re-seeds converged
Ritz values, and LAPACK eigenvectors of the resulting clustered/near-defective
projected matrix have residuals ~eps/gap (measured 6e-7; refinement restores
eps).  Measured basis orthonormality 2e-14, factorization error 1e-13 after
~100 restarts.  Basis lives in one `seq[V]` (ncv+1 vectors); restart rotation
uses a kKeep-vector scratch pool (total memory ncv+kKeep+nev+3 vectors).

**Build / run** (no LAPACK flags needed on this machine — Accelerate is the
in-module default):

```bash
export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
cd /Users/xjin/K/W/P003/qex/.claude/worktrees/reproduce-16-cell-slides-in-qex-d937db/build_mac
make src/experimental/honeycomb/hcarnoldi.nim        && ./bin/hcarnoldi   # self-test, no QEX
make src/experimental/honeycomb/tests/tarnoldi.nim   && OMP_NUM_THREADS=4 ./bin/tarnoldi
```

`tarnoldi` runs in ~6 s.  All field ops in the tests and in the driver are
**serial** (no `threads:` blocks) — correct QEX usage that sidesteps the
loaded-machine spin-barrier issue; only QEX's own `cgSolve` (test 3) threads
internally.  Consequence: tests 0–2 are **bit-identical** across
`OMP_NUM_THREADS` ∈ {2,4} and across reruns; test 3 varies at the 1e-14 level
with thread count (cgSolve reduction order), bit-identical at fixed count.

### Test numbers (all 37 checks PASS)

* **Test 0, conventions verified**: QEX field `dot(x,y) = Σ conj(x)·y`
  (conjugate-linear in the FIRST argument), `norm2 = Σ|x|²`, complex axpy via
  `y += newComplex(a.re,a.im)*x` — all against manual site loops (≤1e-13).
* **Test 1, diagonal operator** (Complex field, 4·4·8·8 = 1024 sites, known
  values, angles restricted to a sector — see caveat below): LM 8/8 in 3
  cycles / 56 applies, worst |ritz−exact| 6.8e-14, worst direct resid
  7.3e-15; SM 8/8 in 90 cycles / 752 applies, worst 3.0e-15 / 7.6e-13.
  (Spec asked 1e-10 / 1e-8.)
* **Test 2, free Wilson–Dirac** (`newWilson`, unit gauge, `setBC`
  antiperiodic time, m = 0.3).  Exact spectrum extracted from the operator
  itself: plane waves (half-integer time momenta) → apply D → 4×4 projection
  → `zgeev`; the projection **completeness defect** (certifies momentum set,
  BC handling, and that D is translation invariant) is 3.5e-15 (4⁴) /
  6.2e-14 (6⁴).  Informational: the extracted 4×4 eigenvalues match the
  closed form `m + Σ(1−cos p) ± i|sin p|` to **1.2e-14** — so the `D` entry
  point is the *standard* Wilson normalization (the "normalized to 2*D_w"
  comment refers to the internal kernels, which the wrapper rescales).
  | run | nev/ncv | cycles | applies (per pair) | worst \|ritz−exact\| | conj-pairing | shells |
  |---|---|---|---|---|---|---|
  | 4⁴ LM | 16/48 | 10 | 224 (14) | 9.8e-14 | 2.0e-14 | 1.6e-14 |
  | 4⁴ SM | 16/48 | 10 | 224 (14) | 1.5e-14 | 1.4e-14 | 4.4e-15 |
  | 6⁴ LM | 16/48 | 15 | 304 (19) | 9.1e-14 | 1.0e-13 | 1.8e-14 |
  | 6⁴ SM | 16/48 | 21 | 400 (25) | 2.6e-14 | 8.4e-15 | 1.4e-14 |
  "shells" = worst distance from any distinct exact eigenvalue among the
  first 16 (in `which` order) to a converged Ritz value — i.e. nothing was
  missed, including the 12-fold degenerate copies (the eigenvalues here come
  in 12+12 conjugate multiplets; the solver reproduces multiple copies).
* **Test 3, shift-invert** (4⁴): apply = `D⁻¹x = (D†D)⁻¹D†x` by CG on the
  normal equations (`cgSolve`, r2req 1e-24, ~17.8 iterations/solve at unit
  gauge; inversion verified directly, |D y−x|/|x| = 6.8e-14).  LM on D⁻¹:
  16/16 in **2 cycles**, 96 outer applies, 1708 CG iterations ⇒ ≈3512 Wilson
  D applies, **220 per converged pair**; worst |1/λ − exact| = 5.7e-15.
  Cost model: shift-invert ≈ 16× the applies of direct SM here, but its
  cycle count is essentially independent of the spectral gap — it is the
  right tool once smallest modes are interior/clustered (rough gauge fields,
  m → 0).

### Things downstream users must know

1. **"SM" without shift-invert only works when the origin is not enclosed by
   the spectrum.**  Restarted Arnoldi builds polynomial filters; no
   polynomial is small on a region surrounding the target.  Wilson at m>0
   (spectrum in Re>0) is fine — direct SM was even cheap here.  For interior
   targets use test 3's shift-invert mode.  Both the module self-test and
   `tarnoldi` document this (the diagonal test's values are restricted to an
   angular sector for exactly this reason; full-circle angles stagnate at
   nconv=0 — verified).
2. **Completeness caveat** (all restarted Krylov methods): with ncv=3·nev on
   a densely clustered synthetic ring, 1 of 8 wanted SM values was crowded
   out for one seed (returned values were still correct eigenvalues with
   tiny residuals — the *set* was short by one).  ncv=4·nev fixed it and was
   ~3× cheaper (`./bin/hcarnoldi 24 0` reproduces; default self-test uses
   ncv=32).  The Wilson tests' first-16-shells checks all pass at ncv=3·nev.
3. `lo.DiracFermion()`/`newWilson(g)` are **VLEN-only**; 6⁴ needs V=1
   (`newLayout(lat,1)`: VLEN=4 cannot lay out odd outer dims).  The generic
   pattern that works for any V (element types must stay Simd-wrapped even
   for V=1 or `ShiftB`'s `pack` breaks on bare float64):
   `type DF = Spin[VectorArray[4, Color[VectorArray[3, type(lo.newDComplexV)]]]]; lo.newField(DF)`,
   plus the `newWilson(@g, proto)` overload.  See `makeWilson` in the test.
4. Nim/QEX codegen gotcha: passing a `seq[cint]` across a generic
   instantiated in a qex-importing module collides with `seq[int32]`
   instances (incompatible C struct names, clang errors).  That is why
   `ztrsenReorder` takes `openArray[int]` and builds the LOGICAL array
   internally.
5. Driver field ops are serial by design; the operator closure owns
   threading (wrap `s.D` in `threads:` for production sizes).  Per-restart
   basis rotation costs kKeep·ncv field axpys; `refineEigvec` costs
   O(ncv³) per reported pair (negligible ≤ 48; ~seconds at ncv≈600).
6. For the eventual "lowest 300 of the honeycomb D" (task D4): point
   `ArnoldiOp.apply` at `hcwilson`'s D via the same adapter, use
   shift-invert, and expect memory = (ncv + kKeep + nev + 3) vectors — at
   ncv=600, kKeep=450, nev=300 that is ~1350 × 40 MB ≈ 54 GB, so plan on
   fewer nev per run (e.g. 4×75 with different shifts) or a leaner rotation
   (Householder-applied Q, 1 temp vector) — noted as future work.

### Open issues

* None blocking.  The completeness caveat (2) is inherent to restarted
  Krylov methods; mitigate with ncv ≥ 4·nev when the target window is
  clustered.
* `hcarnoldi` is deliberately QEX-free, so the ~20-line mixin adapter lives
  in `tests/tarnoldi.nim`; when task D grows `hcSpectrum.nim`, move the
  adapter into a shared module instead of copying it a third time.

---

## 2026-08-21 — Task M: Monte Carlo for the triangle action (HMC, heatbath+OR, I/O, generator)

**Files added** (plus this entry; nothing owned by other tasks touched):

| file | contents |
|---|---|
| `hchmc.nim` | HMC: `HcHmc`, `newHcHmc`, `trajectory`, `revCheck`, `integrate`, `mdSchedule` (leapfrog / 2MN / 4MN5FV, merged stage lists) |
| `hcheatbath.nim` | Cabibbo–Marinari SU(N) heatbath (Kennedy–Pendleton + small-α rejection) and SU(2)-subgroup microcanonical overrelaxation; `HcHeatbath`, `hbSweep`, `orSweep`, `update`, `axisStaples`, `dStaple`, `sampleA0`, `su2HeatbathQuat`, `su2OverrelaxQuat` |
| `hcio.nim` | configuration save/load: `saveHcGauge`, `loadHcGauge`, `hcFileGeom`, `HcMeta` |
| `hcPureGauge.nim` | generation executable: `-algo:hb|hmc|scan`, letParam/installHelpParam/echoParams |
| `tests/thmc.nim` | 9 suites / 28 PASS, exit 0 (exit 1 on any failure — verified with a failing build) |
| `doc/plots/hb_scan.dat` | `<triangleSum>` vs β ∈ [1,12], 23 points, 4⁴ cells (for task R) |

**Build / run**

```bash
export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
cd build_mac
make src/experimental/honeycomb/tests/thmc.nim && OMP_NUM_THREADS=4 ./bin/thmc   # ~11 s, 28/28 PASS
make src/experimental/honeycomb/hcPureGauge.nim
OMP_NUM_THREADS=4 ./bin/hcPureGauge -geom:8,8,8,8 -beta:8 -algo:hb -norsweeps:3 \
  -nwarm:300 -ntraj:2000 -savefreq:10 -outdir:<dir>            # production pattern
# smoke + timing mains: ./bin/hchmc, ./bin/hcheatbath, ./bin/hcio
```

28/28 PASS at `OMP_NUM_THREADS` = 4 **and** 1, bit-identical numbers.  A short
heatbath stream is bit-identical at 1 vs 4 threads (per-site RNG streams).
Verified (awk over the generated C, house rule): **no allocation call inside
any `omp parallel` region** of `hchmc` (9), `hcheatbath` (6), `hcio` (0),
`hcPureGauge` (22), `thmc` (18 regions).

### HMC conventions (M1)

* Momenta: one TAH field per link field (`HcGauge` of TAH matrices), drawn by
  `gaugeUtils.randomTAH` ⇒ P(p) ∝ exp(−|p|²/2).
* `T(p) = ½ Σ_l redot(p_l, p_l)`; `H = T + hcAction`.  With Task A's
  convention `dS/ds|₀ = Σ_l redot(P_l, f_l)` this conserves H — pinned by the
  eps-scaling and `<exp(−dH)> = 1` tests, exactly the `puregaugehmc.nim`
  normalisation (`0.5*p2`).
* Updates: `U ← exp(dt·p) U` (`mdt`), `p ← p − dt·f` (`mdv`, `f = hcForce`).
* Integrators (stage lists identical to QEX's mdevolve; adjacent same-kind
  stages merged across step boundaries): `leapfrog` (nsteps force calls per
  τ), `2MN` λ=0.1931833275037836 (2·nsteps), `4MN5FV` (4·nsteps+1).
* Metropolis: one **lattice-global** uniform from a caller-supplied
  `MRG32k3a` seeded with the broadcast `seed(seed,ix)` — identical decision
  on all ranks (puregaugehmc pattern).  Accept ⇒ `reunit` (projectSU);
  reject ⇒ restore saved copy.

### Heatbath / overrelaxation (M2) — the paper's algorithm

Local action pinned from `hcActionDeriv`, not prose:
`P(U_l) ∝ exp( c·ReTr(U_l Σ_l†) )`, `c = β/(2N)`, `Σ_l =` the raw 8-staple
sum = `hcActionDeriv` output at `beta = 2N` (its prefactor is β/2N).  Pinned
two independent ways (tests 2, 3): the heatbath staple fields equal
`hcActionDeriv(2N)` to **1.7e−32** (relative, squared) on all 24 fields, and
single-link replacements satisfy `dS = −(β/2N)·ΔReTr(UΣ†)` to **3.9e−11**.

**Update schedule** (structural facts verified from `hcgeom.triPath` in test
1, not assumed): every one of the 64 triangles/cell has exactly 1 axis + 2
diagonal edges, and the 2 diagonal edges always carry different δ.  Hence
(a) axis staples contain only uD links ⇒ all 8 axis fields updated in one
pass with staples fixed; (b) each `uD[δ]` field never shares a triangle with
itself ⇒ updated in one pass.  One sweep = axis staples → update
uA[0..3],uB[0..3] → refresh the 4 uA shift-trees → for δ in 0..15: staple of
uD[δ] → update it.  Staple algebra reuses `HcActionWork` (task A) shifters;
per-site SU(2) work happens in plain scalar 3×3 arrays via the `f{i}` lane
proxies with per-site RNG streams (`r{i}`), so results are thread-count
independent.

* SU(2) sampler `sampleA0`: Kennedy–Pendleton for α > 1, plain rejection for
  α ≤ 1 (both exact; split is efficiency only).  Test 6:
  `<a0>` vs I₂(α)/I₁(α) at α = 0.5 / 2 / 8 / 20 with 4M samples:
  diffs 3.2e−5, 9.7e−5, 8.1e−5, 3.9e−6 (all ≤ 1.1σ); the full quaternion
  update on a rotated staple reproduces I₂/I₁ to 1.0e−4.  (I₂/I₁ computed
  two independent ways — series and cosθ-Simpson — agreeing < 1e−8.)
* Overrelaxation: microcanonical reflection `a = (V†)²` per SU(2) subgroup
  (involution, preserves the local exponent identically).  Test 7: one OR
  sweep changes hcAction by **8.0e−16** relative while moving the links by
  per-link |U′−U|² = **2.79**; SU(3) drift without reunit: max 1.4e−12
  (`hcPureGauge` reunits after each hb stream; long streams should reunit
  periodically).

### Test results (tests/thmc.nim, 4⁴/[4,4,4,6] cells, β = 8.0)

| test | result |
|---|---|
| 1 schedule structure | 64 triangles: 1 axis + 2 diag, δ's differ ✅ |
| 2 staples ≡ hcActionDeriv(2N) | worst rel² **1.7e−32** ✅ |
| 3 single-link dS pin | worst rel **3.9e−11** (uA/uB/uD) ✅ |
| 4 reversibility (leapfrog/2MN/4MN5FV) | \|dH_f+dH_b\| = **1.5e−11**, per-link \|U_back−U₀\| ≤ **1.5e−15** ✅ |
| 5 leapfrog dH ∝ eps² | dH(16/32/64 steps) = 0.971/0.193/0.045; ratios **5.04, 4.28** → 4 from above (dH = c₂eps²+c₄eps⁴) ✅ |
| 6 KP moments | table above ✅ |
| 7 OR invariance | rel ΔS **8e−16**, links move 2.79/link ✅ |
| 8 `<exp(−dH)>`, acceptance, HB↔HMC | acceptance **88.7 %** (300 traj, 2MN τ=1 nsteps=10); `<exp(−dH)> = 0.9897 ± 0.0137`; `<triangleSum>`: HMC **0.589676(802)** (τ_int 8.3) vs HB **0.589297(148)** (τ_int 0.98, 800 upd) → **0.46σ** ✅ |
| 9 I/O round trip | bit exact (Σ\|U−U′\|² = 0.0), metadata + geometry round trip ✅ |

Longer independent cross-check (hcPureGauge, 4⁴, β=8): HMC 500 traj
`<ts> = 0.588937(369)`, `<exp(−dH)> = 0.998(12)`, acc 88.6 % vs heatbath
2000 upd `<ts> = 0.589209(94)` — 0.7σ.  Volume consistency at β=8 (HB, 500
upd): 4⁴ 0.589297(148), 6⁴ **0.589088(89)**, 8⁴ **0.588949(50)**.

### β scan (doc/plots/hb_scan.dat, 4⁴ cells, hot starts, 150 warm + 250 meas)

Smooth strong-coupling rise to `ts ≈ 0.32` at β = 6.5, a **steep crossover
between β = 6.5 and 7.5** (`ts` jumps 0.325 → 0.549; τ_int spikes to ~10 at
β = 7.0), then a slow approach to 1.  The requested window
`<triangleSum> ∈ [0.5, 0.65]` is **β ∈ [7.5, 9.0]**:
7.5→0.5491(4), 8.0→0.5893(3), 8.5→0.6213(3), 9.0→0.6483(2).
⚠️ For task R: the jump near β ≈ 7.0 may be a bulk transition — check
hot/cold hysteresis before simulating β ≤ 7.5; the scan used hot starts only.

### Configuration format (hcio.nim)

One SciDAC/LIME file per configuration via QEX's `Writer` (`saveGauge`
machinery): file metadata `<hcGauge><version>1</version><geom>…</geom>
<beta>…</beta><traj>…</traj><info>…</info></hcGauge>`; **one** binary record
with `datacount = 24` — the link fields in canonical order
`uA[0..3], uB[0..3], uD[0..15]`, one 3×3 colour matrix per **cell** per
field, native double ⇒ bit-exact round trip (test 9).  The stored lattice
size is the cell geometry (`getFileLattice`/`hcFileGeom` return it).  RNG
state is not stored (runs log seed + update number).  8⁴ file: 14.2 MB.

### Timings (this loaded box, OMP_NUM_THREADS=4, β = 8)

| lattice | update | cost/update | acceptance |
|---|---|---|---|
| 6⁴ cells (simdlen 1) | HB: 1 heatbath + 3 OR sweeps | **~10 ms** | — |
| 8⁴ cells | HB: 1 heatbath + 3 OR sweeps | **~22–28 ms** | — |
| 6⁴ cells (simdlen 1) | HMC 4MN5FV τ=1 nsteps=4 (17 force) | ~55 ms | 100 % (30 traj) |
| 6⁴ cells (simdlen 1) | HMC 2MN τ=1 nsteps=10 (20 force) | ~56 ms | 70 % |
| 8⁴ cells | HMC 4MN5FV τ=1 nsteps=4 | ~162 ms | 100 %, `<e^{−dH}>`=1.0005(74) |
| 8⁴ cells | HMC 2MN τ=1 nsteps=12 | ~194 ms | 57 % |

(Single sweeps at 8⁴: heatbath 9.2 ms, OR 6.1 ms — `./bin/hcheatbath`
prints these.)  **Heatbath+OR is ~7× cheaper per update than the best HMC
and has τ_int(ts) ≈ 1 vs ≈ 4–8 ⇒ ≳30× cheaper per independent
configuration.  Use `-algo:hb` for production (it is also what the paper
used); keep HMC as the cross-check and for future dynamical work.**
10k updates at 8⁴ ≈ 4–5 min; task R's O(10k+) updates per point are easy.

### Recommended production settings (task R)

| β | `<triangleSum>` (4⁴) | lattice | algorithm | cost/update |
|---|---|---|---|---|
| 7.5 | 0.549 | 6⁴–8⁴ | hb, 1 HB + 3 OR | 10–28 ms |
| 8.0 | 0.589 | 6⁴–10⁴ | hb, 1 HB + 3 OR | 10–50 ms |
| 8.5 | 0.621 | 8⁴–12⁴ | hb, 1 HB + 3 OR | 22–90 ms |
| 9.0 | 0.648 | 8⁴–12⁴ | hb, 1 HB + 3 OR | 22–90 ms |

(β↔a mapping must come from task W's t₀; these are the scan-based starting
points.  Save every 5–10 updates; `-savefreq`, `-loadcfg` support chained
streams.)

### Deviations from the brief / notes

1. Param name `measfreq` (brief wrote `meashfreq`); `ntraj` counts measured
   updates for both algorithms.
2. `hcPureGauge -simdlen:{0,1,2}`: QEX's default `vlen 4` layout **rejects
   6⁴, 9⁴, 10⁴, 11⁴** cell geometries ("can't lay out inner geom", cf. the
   task C note).  `-simdlen:1` runs them (6⁴ verified); 8⁴/12⁴ use the
   default.  Results across simdlen are statistically equivalent, not
   bit-identical (site→RNG-stream pairing changes).
3. `-loadcfg:` added (restart from a saved configuration) — task R will need
   it; `revCheck` keeps the RNG stream aligned across ranks by construction
   (all draws are per-site or global).
4. RNG state is not saved in configurations (documented in hcio.nim).
5. Nim gotcha hit twice: `strformat`'s `&""` cannot interpolate identifiers
   created inside the *same* template body (fails "undeclared identifier" at
   definition); use a helper proc or plain `echo` in templates.
6. Multi-rank and Nc ≠ 3 untested (same status as tasks L/A).  The heatbath
   site kernels are Nc-generic (all N(N−1)/2 SU(2) subgroups) but only
   Nc = 3 was exercised.
7. `trajectory` draws momenta with `gaugeUtils.randomTAH` (non-var RNGField)
   directly rather than `hcaction.randomTAH` (var param cannot be captured
   by the `threads:` closure).

---

## 2026-08-21 — Task W: gradient flow, hexagon clover, E, Q  ⟶ *feeds slide 10*

**Files added** (plus this entry; nothing owned by other tasks touched):

| file | contents |
|---|---|
| `hcflow.nim` | triangle-action Wilson/Lüscher flow, Lüscher RK3 (fork of `gauge/wflow.nim:4-19`), `hcGaugeFlow` template, calibrated `hcFlowCflow` |
| `hctopo.nim` | hexagon clover `hcFmunu`, `hcEQ` (avgE, Q), `HcTopoWork`, `hcCloverSign`, `pairIndex`; full q-prefactor derivation in module docs |
| `hcMeasFlow.nim` | measurement executable: load (task M `hcio`) or warm start, flow, columns `t E t²E t·d/dt(t²E) Q`, stops, local t₀/w₀ interpolation |
| `tests/tflow.nim` | 9 PASS: cubic heat-kernel harness, honeycomb cflow calibration + rational pin, calibrated verification, eps-independence, E-monotonicity, SU(3) preservation |
| `tests/ttopo.nim` | 22 PASS: brute-force reference, weak-field F̂/E, **Atiyah–Singer constant-flux Q (exact)**, gauge invariance, sector recovery, integer clustering |

**Build / run**

```bash
export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
cd build_mac
make src/experimental/honeycomb/tests/tflow.nim && OMP_NUM_THREADS=4 ./bin/tflow  # 63 s, 9/9 (230 s at 1 thread)
make src/experimental/honeycomb/tests/ttopo.nim && OMP_NUM_THREADS=4 ./bin/ttopo  # 170 s, 22/22
make src/experimental/honeycomb/hcMeasFlow.nim
OMP_NUM_THREADS=4 ./bin/hcMeasFlow -geom:8,8,8,8 -warms:0.35 -tmax:2          # warm demo
OMP_NUM_THREADS=4 ./bin/hcMeasFlow -gaugefile:<cfg.lime> -eps:0.02            # hcio config
# smoke + timing mains: ./bin/hcflow, ./bin/hctopo
```

Both tests exit 0 (unittest exit 1 on failure), at OMP_NUM_THREADS = 4 and
(tflow) 1.  House rule verified: a brace-balanced scan of all 204 generated
C files found **~280–292 `omp parallel` regions (across build states) and
zero allocation calls** (`newSeq|newObj|rawNewObj|asgnRef|unsureAsgnRef|
newSeqPayload|allocImpl`) inside any of them.

### W1 — flow normalisation: **cflow = 6, pinned exact** (headline number)

Definition (exact structural parallel of QEX's cubic `gaugeFlow`):

```
  Z = -eps * cflow * nc * (redot-gradient of S_tri at beta = 1)
    = -eps * hcForce(beta = cflow*nc)          [task A convention (★)]
```

* **Cubic harness first** (QEX `gaugeFlow`, weak Abelian plane wave from exact
  line integrals, S(t) of the mode = e^{-2λt}): measured λ/p̂² with
  p̂² = 4sin²(p/2) — the exact linearised lattice prediction — is
  **1.00000199 / 1.00000084 / 1.00001364** at (Ns,k) = (12,1)/(8,1)/(12,2).
  λ/p² = 0.977363 / 0.949642 / 0.911903, 3-point extrapolation
  **1.000000**, O(p²) coefficient **−0.0833** = the exact −1/12.  The
  harness and QEX's flow-time convention are thereby pinned.
* **Honeycomb with provisional cflow = 1**: λ/p² = 0.1660338 / 0.1652466 /
  0.1641535 at the same momenta (plateau drift ≤ 2.2e−5).  Removing the
  O(p²) artifact: c_HC = **0.1666673 ± 3.8e−6** (3-pt quad vs 2-pt linear
  spread), i.e. `cflow = 1/c_HC = 5.9999771 ± 1.4e−4` →
  **cflow = 6 exactly** (|dev| = 2.3e−5, the only simple rational within
  1e−3).  FORMULATION §4.4's rough estimate ("1/6" for the rate constant)
  was right.
* **Verification at cflow = 6**: λ/p² = 0.9962001 (12,1), 0.9914820 (8,1);
  artifacts −0.0038/−0.0085; extrapolated **0.999974** — |dev| = 2.6e−5,
  acceptance `< 1e-2` met with ~400× margin.  RK3 step-size independence:
  rate(eps=0.05) vs rate(0.025) differ by **2.1e−6**.  Linearity in cflow
  checked at cflow=4.5 (λ/p² extrapolated 0.749899 = 4.5/6 to 1.3e−4).
* Mode amplitude 3e−3: at 1e−3 the extracted rate sat close to the
  cancellation-noise floor of the action sum and moved by ~5e−5 (relative)
  between OMP thread counts (1-thread cflow dev from 6 was 1.7e−3, outside
  the 1e−3 pin gate); at 3e−3 both thread counts pin cleanly
  (dev 2.3e−5 at 4 threads, 2.0e−4 at 1), nonlinear bias O(ε²) ~ 1e−5.
* Analytic footnote (in `hcflow.nim` docs): the Rayleigh quotient of the
  naively sampled continuum mode gives (2/9)·cflow·p² — an *upper bound
  only* (the true slow eigenvector has an O(p) axis/diagonal optical
  admixture).  Do not "re-derive" 9/2; the heat kernel says 6.
* **Bonus**: the honeycomb flow-rate O(p²) artifact coefficient is −0.0024,
  ~35× smaller than the cubic −1/12.
* `flowScaleSetup`: with `cflow = hcFlowCflow` the `wflowT` handed to
  `measure` is continuum flow time in units a², 1:1 comparable to cubic
  flow time (t₀, √(8t) smoothing radius, etc.).
* Flow preserves SU(3) to **2.3e−15** over 60 RK steps (no reunit); E from
  the clover decreases strictly monotonically (1.4678 → 1.6e−5 over t=0→3).

### W2 — hexagon clover, E, Q

Implementation: the 2×16×6 = 192 triangle-loop recipes are generated at
module load from `hcgeom.hexTriPaths` (never hand-derived), each evaluated as
an anchored field product with ≤2-hop single-axis shift chains (association
order chosen per loop to minimise shifts), then
`F̂_ab = Σ_h [hcCloverSign·Ω^{(h)}_ab/(4√3)]·TAH[C_h]` (weights ±1/12).
All shifters/scratch live in a persistent `HcTopoWork`; calls are
allocation-free.

**q-prefactor derivation** (independent, documented in `hctopo.nim`; agrees
with QEX's cubic `topoQ` as the lead's update predicted):

```
  q(x) = -(1/32π²) ε_{μνρσ} Tr[F̂_μν F̂_ρσ]          (F̂ anti-Hermitian = i a² F_H)
       = -(1/4π²) [ t(F̂01F̂23) - t(F̂02F̂13) + t(F̂03F̂12) ]     (ε-contraction = 8 terms/partition)
  Q    = ½ Σ_x q(x)                                  (a⁴/2 per site)
```

**Sign finding (FORMULATION (4.1) correction needed):** hcgeom's hexagon
ring `d⁻, d⁺, ê_μ, …` runs **clockwise** with respect to `omega`'s 2-form
(angles 120°,60°,0°,…), so each triangle loop encloses flux −(√3/4)a²F_Ω and
(4.1) must read `F̂_Ω = −(4/√3)·TAH[(1/6)C_h]` for F̂ = +i a²F·T with
U = exp(+i∫A·dl).  Implemented as `hcCloverSign = -1`; pinned by the two
exact tests below (a wrong sign flips Q and the site-wise F̂ ratio to −1).
The magnitudes 4/√3, 3/8 and the ½ in Q are all confirmed exactly.
*(FORMULATION owner: (4.1) needs the minus given hcgeom's ring order — or
equivalently hexTriPaths' "positive orientation" comment is w.r.t. −Ω.)*

Validation (tests/ttopo.nim, all PASS):

1. **Brute force**: independent single-site implementation (hand 3×3 complex
   algebra, literal (4/√3)/(3/8)/Ω sums) vs `hcFmunu`/`hcEQ` on random
   configs: max site-wise |F̂ diff| **3.3e−16**; avgE rel **3.3e−16**;
   Q **1.9e−16**; identical on V=1 and on [2,4,4,6] (wrapping shifts).
2. **Weak-field F̂, site by site** (plane wave, exact line integrals, 12⁴):
   max|F̂ − iF_exact T|/max|F| = **0.0189 (k=1) → 0.0738 (k=2)**, ratio
   **3.905** (O(p²) ✓, incl. the [1,1] = −iφ and off-diagonal structure).
   E/E_exact = **0.96258 / 0.85792**, (1−E) ratio **3.797** (O(p²) ✓).
3. **Atiyah–Singer constant-flux test (the decisive one, exact)**: Cartan
   T = diag(1,−1,0) configuration with fluxes (n1,n2); construction = task
   C's cubic recipe generalised to the 24 link types by exact line integrals
   + the boundary transition function at the *endpoint* (−f1·L0·x1(end);
   reduces to refCubicMeas's on the cubic lattice; needed because diagonal
   links cross the seam with nonzero transverse displacement).
   `Q/(2n1n2)` = **0.99679136 (L=4), 0.99979922 (L=8), 0.99996034 (L=12)**;
   artifact scaling **15.98 ≈ (8/4)⁴** and **5.06 ≈ (12/8)⁴** — pure 1/L⁴;
   sign-odd under n2 → −n2; (n1,n2) = (2,1) gives Q/4 = 0.999498;
   E/E_cont = 0.999799 at L=8.  A missing ½ would read 2× exactly.
   **Bonus:** dev·L⁴ = 0.8224 = (2π)²/48 — the hexagon clover's Q artifact
   is exactly **16× smaller than the cubic 1×1 clover at equal L** (cubic
   dev·L⁴ = 2(2π)²/6 = 13.16; our L=8 number coincides digit-for-digit with
   task C's cubic 16⁴ value 0.99979922).
4. Gauge invariance: E rel **2.5e−14**, Q abs **1.3e−14** under random
   SU(3) vA,vB.
5. **Sector recovery (sharp integer-Q test)**: exact flux configs kicked by
   U → exp(0.35·X)U with random TAH X (raw clover Q then reads 0.55–1.43
   for the Q=2 sectors), flowed to t=2: Q returns to the exact sector
   integer, max |Q − Q_exact| = **0.0029** over sectors {0, ±2, 4},
   plateau max |Q(2)−Q(1)| = **0.0013**.
6. **Integer-Q clustering on generic rough configs (secondary, loose)**:
   12 × warm(0.85) 8⁴ configs flowed to t=10: histogram of round(Q):
   **{−1:4, 0:6, 1:1, 2:1}**, mean |Q−int| = **0.094** (uniform: 0.25;
   nearest-half-odd: 0.406), 6/12 with |Q|>0.5, max dist 0.304 (one config
   still mid-anneal).  Exactly task C's cubic finding: cutoff-scale lumps
   keep dislocating, so snapshot clustering is loose — warm(0.7) configs
   instead all anneal to Q = 0.0000 exactly, and warm(0.85) at t=6 is still
   turbulent (mean dist 0.27).  The sharp statements are tests 3 and 5.

### W3 — hcMeasFlow

Pattern of `examples/wflow_topo.nim`; task M's `hcio` landed in time and is
integrated (`-gaugefile:` loads, geometry from the file; `-simdlen:1` for
6⁴-type geometries; warm start otherwise).  Prints `HCFLOW t E t²E
t·d/dt(t²E) Q`, stops on tmax/t2Emax/tdt2Emax, then `T0 <t0> <Q(t0)>` and
`W0SQ <w0²>` from **local linear interpolators** (deliberately not importing
task C's concurrently-edited `hcanalysis.nim`; task R should consolidate on
`hcanalysis.findT0/findW0`).

Demo on a real β=8 heatbath configuration (task M `hcPureGauge`, 8⁴):
E(0) = 3.03, flow stops at t·d/dt(t²E) = 0.35 at t = 1.88, w0²/a² = 1.73,
t²E = 0.24 still rising ⇒ t₀ not reachable on 8⁴ at β=8 (t₀ ≈ 2.2 would give
L/√t₀ ≈ 5.4 < 9 — the finite-volume bound task C established; the tool says
so explicitly).  Q plateaus at −2.596 for t ∈ [1.4, 2.4], then a dislocation
near t ≈ 4.5 moves it towards −2 — the coarse-clover behaviour task C
documented on the cubic side, reproduced here.
**Task R note: 8⁴ at β = 8 is too small/coarse for a t₀-based χ_top point;
push to 10⁴–12⁴ (simdlen 1 for 10⁴) and/or larger β.**

### Timings (idle box, OMP_NUM_THREADS=4; task A's 49 ms force figure was
taken at load average 25–65)

| operation | 8⁴ cells |
|---|---|
| `hcForce` | 3.6 ms |
| `hcAction` | 1.3 ms |
| `hcEQ` (E and Q, full clover) | **7.8 ms** |
| RK3 flow step (3 forces) + per-step (E,Q) | **31.4 ms** |
| full flow+measure to t = 1.9 (94 steps, per-step meas) | **2.95 s** |

12⁴ scales ×5.06 (volume).  A per-configuration measurement (eps 0.02,
measevery 2–5, t ≤ 3) costs 3–10 s at 8⁴, 15–50 s at 12⁴ — task R's
ensembles are cheap to measure.

### Public API

```nim
# hcflow.nim
const hcFlowCflow* = 6.0                 # calibrated; t is continuum flow time (a²)
template hcGaugeFlow*(g: HcGauge; steps: int; eps, cflow: float; measure)
template hcGaugeFlow*(g: HcGauge; eps, cflow: float; measure)
template hcGaugeFlow*(g: HcGauge; eps: float; measure)   # cflow = hcFlowCflow
  # injects wflowT; `break` in measure stops; allocates its own HcActionWork once

# hctopo.nim
const hcCloverSign* = -1.0               # see the (4.1) sign finding above
func pairIndex*(a, b: int): int          # (1,0)=0 (2,0)=1 (2,1)=2 (3,0)=3 (3,1)=4 (3,2)=5
type HcTopoWork*[V,F,SS] = ref object    # f*: array[2, array[6, F]]  = F̂_ab per sublattice
proc newHcTopoWork*[V,F](g: HcGauge[V,F]): auto          # outside threads:, once
proc hcFmunu*(w: HcTopoWork, g: HcGauge)                 # fills w.f, allocation free
proc hcEQ*(w: HcTopoWork, g: HcGauge): tuple[e, q: float] # (avgE, Q = ½Σq)
```

### Open issues / deferred

1. **Multi-rank untested** (same status as L/A/M; every shift is a plain
   single-axis `Shifter`).
2. **Nc ≠ 3 untested**; `hctopo` itself is Nc-generic, the tests' Abelian
   embeddings assume Nc = 3.
3. `hcEQ` recomputes all 192 loop products per call; sharing radial-link
   subproducts across a hexagon could save ~2×.  Correctness-first kept
   (7.8 ms is cheap enough).
4. The t₀/w₀ interpolators in `hcMeasFlow` are local and linear — task R:
   consolidate on `hcanalysis.findT0/findW0` (cubic-validated, cubic
   interpolation).
5. `tests/ttopo.nim` test 6 is seeded and deterministic, but its assertions
   are statistical by nature (mean dist < 0.2 with measured 0.094); if a
   future RNG change trips it, loosen with the histogram in view.
6. FORMULATION §4.1 sign and §4.4 constant: reported above for the owner;
   only this file and my own code state them.

---

## 2026-08-21 — Task D2: stout smearing + tree-level clover improvement (honeycomb AND cubic)

**Files added** (plus this entry; nothing owned by other tasks touched):

| file | contents |
|---|---|
| `hcstout.nim` | Morningstar–Peardon stout on the honeycomb (`HcStout`, triangle staples via task A's `HcActionWork`) + thin uniform wrapper `CubicStout` around QEX `gauge/stoutsmear` (REUSED, not reimplemented); `hcStoutKappa`, `cubicStoutKappa` |
| `hcclover.nim` | clover-improved Wilson–Dirac operators: `HcCloverWilson` (wraps task D1 `HcWilson` + task W `HcTopoWork`), `CubicCloverWilson` (wraps QEX `physics/wilsonD` + QEX `fmunu(g,1)`); `hcCloverGam`, `cubicFmunuSign` |
| `tests/tstout.nim` | 16 PASS: SU(3) exactness, gauge covariance, rho=0 identity, monotone smoothing, **heat-kernel kappa measurement both lattices**, timing |
| `tests/tclover.nim` | 30 PASS: cSW=0 bit-reduction, gauge covariance, gamma5-hermiticity + Hermiticity/gamma5-evenness of the added term, on-site locality, **constant-flux normalisation pin both lattices** (+ chirality splitting), free-field identity, timing |

**Build / run**

```bash
export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
cd build_mac
make src/experimental/honeycomb/tests/tstout.nim  && OMP_NUM_THREADS=4 ./bin/tstout   # ~55 s, 16/16
make src/experimental/honeycomb/tests/tclover.nim && OMP_NUM_THREADS=4 ./bin/tclover  # ~1.5 s, 30/30
# smoke + timing mains: ./bin/hcstout, ./bin/hcclover
```

Both exit 0 at `OMP_NUM_THREADS` = 4 and 1 (kappa shifts only at the 4e-6
level between thread counts — the tflow-documented action-sum noise floor);
4-thread reruns bit-identical.  House rule verified (brace-balanced scan of
all 242 generated C files): **350 omp parallel regions, zero allocation calls
inside any**.

### D2.1 — stout smearing

**Honeycomb = MP with the 8-triangle staples.**  `Sigma_l = hcActionDeriv` at
`beta = 2N` (raw staple sum, prefactor beta/2N verified in the code and by
task M's 1.7e-32 pin), so the MP step `U' = exp(projectTAH(rho Sigma U^dag))U`
is implemented as `U' = exp(-rho*projectTAH(U Sigma^dag))U` — the
gradient-descent direction of the triangle action.  **Cubic = QEX
`StoutSmear` reused**: its `exp(-alpha*nc*projectTAH(U ds^dag))U` with
`ds = gaugeActionDeriv(plaq=1) = (1/nc)*staple-sum` **is exactly the isotropic
MP step with rho = alpha** (shown in hcstout.nim docs); `CubicStout` only
gives it the same `smear`/`smearN` interface.

**Heat-kernel constant kappa (t_eff = n·rho·kappa per n steps), the headline
numbers** — measured with the exact-line-integral Abelian plane wave
(taction/tflow recipe), per-step action decay, 3 momenta, O(p²) artifact
removed by quadratic extrapolation:

| lattice | analytic | measured (rho=0.05) | rho=0.025 | per-step pin |
|---|---|---|---|---|
| cubic | **1** (stout step ≡ Euler `gaugeFlow` step at eps=rho; Luscher 2010) | **1.000013** | — | factor = `-ln(1-rho*phat^2)` to **4.8e-6** |
| honeycomb | **1/3** (= 2/cflow, exact given task W's cflow=6: stout step = Euler flow step of eps=rho·2/cflow) | **0.3333308** (rational scan: 1/3, dev 2.5e-6) | 0.3333320 | — |

Raw rates (rho=0.05): cubic rate/p² = 0.983966/0.963825/0.935473 at
p² = 0.274/0.617/1.097, O(p²) coeff −0.0584 (= flow's exact −1/12 + Euler
rho/2 = −0.0583 ✓); honeycomb 0.3328247/0.3321896/0.3312962, O(p²) coeff
**−0.0018** — the honeycomb stout kernel is ~32× more rotationally clean,
echoing task W's flow-artifact finding.

> **rho equivalence for task D4.**  t_eff per step: cubic **rho**, honeycomb
> **rho/3**.  The paper's "6 steps at rho = 0.05 on both lattices", read in
> the plain MP staple-sum convention (C = rho·Sigma with Sigma the unweighted
> 8-triangle staple sum — the only convention stated), gives smearing radius
> sqrt(8 t_eff) = **1.55 a (cubic)** vs **0.89 a (honeycomb)**: NOT equal
> smearing.  If D4 wants equal physical smearing instead, use rho_hc = 3
> rho_cubic (or 3× the steps).  Since a is ~2× coarser on the honeycomb at
> equal physics in the paper's programme, the physical radii are closer than
> the lattice-unit ones; D4 should simply state which convention it runs.

Other tstout numbers: smeared links SU(3) to **6.4e-16** (hc) / **4.9e-12**
(cubic, dominated by the input `random` config's own projectSU accuracy);
covariance `smear(U^g) = (smear U)^g` to **9.7e-15 / 7.6e-15**; rho=0
identity **exact (diff 0.0)**; warm-8⁴ smoothing sequences (rho=0.05):
triangleSum 0.602948 → 0.800535 → 0.912025 → 0.960301 → 0.979971 → 0.988693 →
0.993064, plaq 0.508376 → 0.703829 → 0.839190 → 0.916270 → 0.955984 →
0.975961 → 0.986206 (both strictly monotone).

**Analytic remark** (also in the module docs): stout smearing commutes
exactly with `setBC`'s Z2 time twist (every staple has the same crossing
parity as its link), so smear-then-setBC and setBC-then-smear agree; use the
canonical order smear → `setBC` → build/refresh the operator anyway.

### D2.2 — clover term: the coefficient and its pin

Implemented operator (both lattices, SAME definition; `pairIndex`/Fhat
conventions of hctopo):

```
  D_c psi = D psi − (cSW·r_w/2) Σ_{a>b} γ_a γ_b Fhat_ab psi
          =  D psi − (cSW·r_w/4) Σ_{μν} σ_μν F_μν psi          (continuum form)
  σ_μν = (i/2)[γ_μ,γ_ν],   Fhat_ab = +i a²F_ab (TAH clover field)
```

**Derivation of the coefficient** (hcclover.nim docs): both operators are
`m + γ·D − (r a/2)D² + O(a²)` — the naive term's O(a) piece is odd in n and
cancels; the Wilson term is `−(r a/2)·(1/6)Σ_24 (n·D)² = −(r a/2)D²` with NO
F-term of its own (n_μ n_ν symmetric ⇒ only {D_μ,D_ν}) — the honeycomb 1/6
and cubic 1/2 normalisations both land on the same `a r p²/2` (pinned by
tasks F/D1).  Adding −(r a/4)σ·F (cSW=1) makes the O(a) term
`−(r a/2)(γ·D)²`, a pure on-shell mass shift: standard tree-level improvement,
**identical coefficient on both lattices** (unlike the flow's cflow=6, no
lattice-dependent factor appears — the 1/6 is already inside D²).

Field strength: honeycomb `hctopo.hcFmunu` (Fhat = a²F pinned by task W);
cubic QEX `fmunu(g,1)`.  **Measured (new pin): QEX `fmunu(g,1)[a][b]` (a>b)
= +i F_ab T on the constant-flux background — the SAME +iF_ab convention as
hctopo pair storage; `cubicFmunuSign = +1`** (this was previously unpinned:
topoQ/densityE are quadratic in F and blind to it).

**The acceptance pin** (tclover test 5): on the exact Atiyah–Singer
constant-flux background (T = diag(1,−1,0), F_01 = f1 = 2πn1/L², F_23 = f2),
the on-site 12×12 matrix (D_c − D)(x), extracted from 12 point sources, must
equal `(i/2)[f1 s1 γ2γ1 + f2 s2 γ4γ3] ⊗ T` with the **exactly known clover
artifact factors s_hc = 4 sin(f/4)/f and s_cubic = sin(f)/f**:

| lattice | L | s_meas (= measured/continuum-analytic) | s_exact | \|meas − exact\| | structure resid |
|---|---|---|---|---|---|
| honeycomb | 4 | 0.998394393036 | 0.998394393036 | 5.6e-16 | 1.1e-15 |
| honeycomb | 8 | 0.999899604216 | 0.999899604216 | 3.3e-16 | 1.1e-16 |
| honeycomb | 12 | 0.999980168255 | 0.999980168255 | 2.0e-14 | 6.3e-16 |
| cubic | 8 | 0.998394393036 | 0.998394393036 | 6.2e-15 | 1.7e-16 |
| cubic | 12 | 0.999682720392 | 0.999682720392 | 2.7e-15 | 4.7e-17 |
| cubic | 16 | 0.999899604216 | 0.999899604216 | 1.8e-14 | 1.9e-16 |

measured/lattice-analytic = 1 to **2e-14** (worst case, all L, both planes,
both lattices); measured/continuum = 1 − O(1/L⁴) with ratios
**16.0 (L=4→8) / 5.06 (8→12)** hc and **5.06 (8→12) / 3.16 (12→16)** cubic —
the artifact scaling demonstrated.  (Note the brief anticipated 1 + O(1/L²);
on this constant-flux background f ∝ 1/L² so the artifact is O(f²) = O(1/L⁴).)
Bonus symmetry: the hexagon-clover argument is f/4 vs the cubic f, so
**s_hc(L) = s_cubic(2L) digit-for-digit** — the honeycomb clover term is 16×
closer to continuum at equal L, matching task W's factor-16 Q-artifact
finding.  The hc matrix is constant across sublattices and across the
transition-function seam (diffs ≤ 9e-17).

**Chirality-splitting cross-check** (relative-sign/factor errors downstream):
for f1 = f2 ("self-dual") the added term vanishes on the **γ5 = +1** spin
sector (|P₊CP₊| ~ 1e-17 vs |P₋CP₋| = 0.196), for f1 = −f2 on γ5 = −1,
**identically on both lattices** (σ_01σ_23 = −γ5 in the DeGrand–Rossi basis).
A wrong overall sign flips s to −1 (excluded at 1e-14), a wrong relative
plane sign or an hc/cubic mismatch breaks the null sector (excluded at 1e-17).

Other tclover numbers: cSW=0 → bare and unit-gauge (F̂ exactly 0) → bare both
give **norm2diff = 0.0**; gauge covariance of D_c **4.1e-14 / 1.4e-14**;
adjointness ⟨x,D_c y⟩=⟨D_c†x,y⟩ **5.0e-16 / 4.4e-16**; γ5-hermiticity
⟨x,γ5D_cγ5 y⟩=⟨D_c x,y⟩ **4.9e-16 / 1.5e-15**; support of D_c δ unchanged
(25 sites hc / 9 cubic), added term strictly on-site (off-site diff exactly
0.0 on random gauge).

> **Brief correction, verified numerically**: the added term C = D_c − D is
> **HERMITIAN** (⟨x,Cy⟩=⟨Cx,y⟩ to 1.9e-15; the anti-Hermitian combination is
> O(1), measured 2.00), and γ5-even (1.1e-15).  It must be: γ5 C γ5 = C (σ
> commutes with γ5), so γ5-hermiticity of D_c ⇔ C† = C.  The brief's
> "anti-Hermitian" expectation is a slip (the *naive γ·D part* of D is the
> anti-Hermitian piece).

### Public API (task D4 consumes this)

```nim
# hcstout.nim
const hcStoutKappa* = 1.0/3.0        # t_eff = n*rho*kappa
const cubicStoutKappa* = 1.0
type HcStout*[V,F,W] = ref object    # rho*, w* (HcActionWork), d* (staple sums)
proc newHcStout*(g: HcGauge, rho: float): auto            # outside threads:
proc smear*(s: HcStout, gin: HcGauge, gout: var HcGauge)  # opens threads:; gout may alias gin
proc smearN*(s: HcStout, g: HcGauge, gout: var HcGauge, n: int)
type CubicStout*[G] = ref object     # wraps gauge/stoutsmear.StoutSmear
proc newCubicStout*(lo: Layout, rho: float): auto
proc smear*(s: CubicStout, gin, gout: G)                  # in-place ok
proc smearN*(s: CubicStout, g, gout: G, n: int)

# hcclover.nim  (exports hcwilson, hctopo, wilsonD)
const cubicFmunuSign* = 1.0
let hcCloverGam*: array[6, SpinMat]  # gamma_a gamma_b, pairIndex order
type HcCloverWilson*[W,TW] = ref object   # w* = bare HcWilson, tw* = HcTopoWork, cSW*
proc newHcCloverWilson*(g: HcGauge, cSW: float): auto     # computes Fhat; outside threads:
proc gaugeRefresh*(c: HcCloverWilson)     # uDsh AND Fhat, after links change in place
proc D*(c: HcCloverWilson, r: var HcFermion, x: HcFermion, m: float, rw = 1.0)
proc Ddag*(c: HcCloverWilson, r: var HcFermion, x: HcFermion, m: float, rw = 1.0)
type CubicCloverWilson*[W,MF] = ref object  # s* = physics/wilsonD Wilson (r_w = 1), f*[6], cSW*
proc newCubicCloverWilson*(g: seq[G], cSW: float): auto        # VLEN layouts
proc newCubicCloverWilson*(g: seq[G], v: T, cSW: float): auto  # generic-V (proto fermion)
proc gaugeRefresh*(c: CubicCloverWilson)  # QEX fmunu: ALLOCATES, once per config
proc D*(c: CubicCloverWilson, r: var auto, x: auto, m: SomeNumber)
proc Ddag*(c: CubicCloverWilson, r: var auto, x: auto, m: SomeNumber)
```

Conventions/pattern for D4: load config → `smearN(st, g, g, 6)` →
`threads: g.setBC` → `newHcCloverWilson(g, 1.0)` (or `gaugeRefresh` an
existing one — it recomputes both the pre-shifted links and F̂); D/Ddag open
their own `threads:` (do not call from inside one; serial drivers like
tarnoldi's work as-is); `r` must not alias `x`; cSW=0 short-circuits to the
bare operator; the honeycomb clover term scales with the per-call `rw`
(tree-level term ∝ r), the cubic wrapper has r_w = 1 fixed by QEX wilsonD.
Chirality for Q_Dirac: `applyGamma5` (hc) / `gamma5 * x[e]` (cubic) — same
DeGrand–Rossi γ5 on both lattices, and test 5 fixed the σ·F ↔ γ5 correlation
identically on both.

### Timings (8⁴, OMP_NUM_THREADS=4, this box near-idle)

| operation | honeycomb (8⁴ cells) | cubic (8⁴ sites) |
|---|---|---|
| one stout step | **10.96 ms** (21.6 ms under load) | **3.36 ms** |
| D_c apply | **1.56 ms** (bare D 1.32 ms ⇒ +18%) | **0.30 ms** |
| gaugeRefresh (once per cfg) | 8.05 ms | 4.66 ms (QEX fmunu, allocates) |

### Open issues / deferred

1. **Multi-rank untested** (same status as L/A/M/W/D1); everything is built
   from the same single-axis Shifter/ShiftB machinery, so it should be
   MPI-correct.
2. **Nc = 3 only tested** (the operators and hcCloverGam are Nc-generic; the
   flux-background tests embed via T = diag(1,−1,0)).
3. `CubicCloverWilson.gaugeRefresh` allocates (QEX `fmunu` work fields +
   `GC_fullCollect`) — fine once per configuration; a persistent cubic clover
   workspace would need a QEX fmunu refactor (not worth it).
4. Each D_c application launches two `threads:` regions (bare D + clover
   add); fusing them needs an hcwilson refactor — the measured 18% overhead
   does not justify it.
5. No stout force/chain rule on the honeycomb (quenched measurements only,
   per brief); the cubic wrapper's inner `StoutSmear` still has QEX's
   `smearDeriv` if ever needed.

---

## 2026-08-21 — Task R: physics runs — topological susceptibility scaling (slide 10 / Figs. 1–2)

Full write-up with all numbers and plots: **[RESULTS.md](RESULTS.md)**.
**Owns / added:** `doc/RESULTS.md`, `doc/plots/honeycomb/*` (run scripts,
harvester, fitter, plot scripts, per-point data, chitop.png, qhist.png).
Nothing owned by other tasks was touched.  No compiled analyzer was needed
(`hcChiAnalyze.nim` not created): hcMeasFlow's per-invocation overhead is
only ~30 ms, so shell loops + a pure-python harvester
(`doc/plots/honeycomb/harvest.py`, conventions mirroring `hcanalysis`)
sufficed.

### Headline

Six quenched ensembles (heatbath + 3 OR, `hcPureGauge`; flow + hexagon
clover, `hcMeasFlow`) at β = 6.90–7.20 on 8⁴/12⁴ cells, a²/t0 = 2.64 → 0.60,
n_eff(Q²) = 92–525 per point:

* O(a⁴) fit: **10⁴t0²χ|_{a=0} = 6.38(31)** (all 6 points, χ²/dof 0.88);
  **6.60(42)** on the slide-10 window x ≤ 2.1 — vs Cè et al. **6.67(7)**
  and slide 10's ≈ 6.78. ✅
* Freeing the O(a²) term: **c2 = −0.51 ± 1.39, consistent with zero**
  (slide 10's headline), vs the cubic reference's −2.70(51) (task C, same
  pipeline). ✅
* At a²/t0 ≈ 1.5 the honeycomb retains ≈ 83 % of its continuum χ, the
  cubic ≈ 32 % — factor 2.8(4) smaller cut-off effect at equal a²/t0. ✅
* Q(t0) histograms symmetric, ⟨Q⟩ = 0 within errors, integer peaks visible
  at the finest point, washed out at a²/t0 = 2.6 (expected for the plain
  clover at coarse a). ~paper Fig. 1.

### Things later tasks / reruns should know

1. **The scale moves fast**: t0/a² = 0.23 at β = 6.7 but ≳ 2.5 by β = 7.3
   (d ln t0/dβ ≈ 5–6 through the crossover).  The whole slide-10 window
   lives at **β ∈ [6.9, 7.2]**; the pre-run guess window β ≥ 7.5 has
   t0 ≳ 4 (no t²E = 0.3 crossing even on 12⁴ up to t = 4.2).  hb_scan's
   β ≈ 7.0 feature is a **smooth crossover**: hot/cold hysteresis clean at
   β = 6.8…7.3 on 8⁴ (second-half ⟨ts⟩ diffs ≤ 6·10⁻⁴), τ_int(ts) peaks
   ≈ 13 updates at β = 6.8.
2. **Honeycomb flow RK3 stability**: eps = 0.05 accurate (1e−4 in t²E vs
   eps 0.02), **eps = 0.10 blows up**.  Boundary between 0.05 and 0.10.
3. τ_int(Q) with 1 HB + 3 OR: ≈ 2 updates at a²/t0 = 2.6 → ≈ 43 updates at
   0.60 (12⁴).  Volume rule confirmed again: β = 7.20 on 8⁴ has
   L/√t0 = 6.0 and 7/20 configs never cross t²E = 0.3 (t0 unmeasurable);
   on 12⁴ the same β is healthy (L/√t0 = 9.3).
4. The two harness kills (15:03, 16:11) cost wall time but no data:
   `resume_point.sh` re-flows surviving configs and continues chains from
   the last saved configuration; measurement series stay contiguous.
   `run_recover*_slot*.sh` document the exact restarts.
5. Reproduction: `sh doc/plots/honeycomb/run_calib_{A..E}.sh` (calibration,
   3 concurrent), `sh doc/plots/honeycomb/run_prod_slot{A,B,C}.sh`
   (production), `sh doc/plots/honeycomb/harvest_all.sh` (analysis+plots).
   All jobs `OMP_NUM_THREADS=4`, ≤ 3 concurrent; configs live under
   `$TMPDIR/hcR` and are deleted after measurement (12⁴ files are 72 MB).

### Open / not done

* a²/t0 < 0.6 needs L ≥ 16 cells and ~10× the decorrelation budget — out
  of scope here; the O(a⁴)-vs-O(a²) discrimination below x = 0.6 rests on
  the intercept comparison with Cè et al., not on our own points.
* Multi-rank still untested anywhere in the honeycomb stack (inherited).
* The b7.07 point sits at L/√t0 = 8.1 (edge of the volume rule); flagged
  in RESULTS.md §6.

## 2026-08-21 — Task D4 close-out (lead harvest after user-requested stop)

The D4 agent was stopped at 17:57 on user request ("stop if the handoff explains how
to resume").  Harvest of the logs as they stood (hc 66 cfgs ρ=0.05 + 8 cfgs ρ=0.15,
cubic 60 cfgs): see `doc/RESULTS_FERMIONS.md` for all numbers and the verdicts
(slide 16 ✅, slide 19 ✅, slide 18 volume-limited ◐).

**Known bug left open**: `doc/plots/fermions/harvest.py` does not join the honeycomb
`FLOWQ` lines into the per-config Q_flow (the CFG line's field 11 is 0 for hc) — the
agent was mid-fix when stopped.  The external join used for the results:
`Q_flow(cfg i) = FLOWQ field 3` (first plateau sample at the per-config t0), applied
when CFG field 11 is 0.  `data/fig7_hc*.dat` were regenerated with this fix;
`fig7_qdiff.png` replotted.  Anyone rerunning `harvest.py` must re-apply the fix or
port it into the script.

The unkillable generation jobs (`hc722b`, target 76; `cub586`) may append past the
harvest cut (n = 66/60/8 used); ignore the tail or re-harvest with more configs.
Resume instructions: HANDOFF.md §3.
