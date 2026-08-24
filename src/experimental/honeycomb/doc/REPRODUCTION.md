# Reproduction of "Lattice QCD on the 16-cell honeycomb" — final report

Target: S. D. Katz, D. Nogradi, *Lattice QCD on the 16-cell honeycomb*, Lattice 2026
(slides, 22 pp.) and arXiv:2512.10604 (paper).  Everything below was implemented in QEX
under `src/experimental/honeycomb/` in this worktree; every number was produced by code
in this tree.  Details: [SLIDES.md](SLIDES.md) (what the talk shows),
[FORMULATION.md](FORMULATION.md) (normative formulation), [PLAN.md](PLAN.md) (work
breakdown), [STATUS.md](STATUS.md) (per-task logs with all numbers),
[RESULTS_CUBIC.md](RESULTS_CUBIC.md), [RESULTS.md](RESULTS.md),
RESULTS_FERMIONS.md (fermion runs).

**Statistics disclaimer.** All Monte-Carlo results here are *preliminary, laptop-scale*
(hundreds of configurations per point, single shared machine, one afternoon), intended
to establish correctness of the implementation and visibility of the signals.  The
paper used O(50 000) configurations per point.  A production rerun on a cluster only
needs bigger numbers in `hcPureGauge`/`hcSpectrum` inputs — no code changes.

---

## 1. Slide-by-slide verdict

| slide | content | our result | verdict |
|---|---|---|---|
| 4 | point-group orders: cubic 384, 16-cell 1152 | `hcgeom`: **384 / 1152** by explicit enumeration | ✅ exact |
| 7–8 | D₄* lattice, 24 unit neighbours (8 axis + 16 diagonal) | `hcgeom` + `tgeom` (42 checks) | ✅ exact |
| 9 | 12 links/site, 32 triangles/site, `β = 2N/g²` unchanged | 32 apex triangles verified; classical-limit test: `S₁₆/S_cubic = 1 + 0.056 p²` → 1, O(p²) scaling confirmed | ✅ |
| 10 | `10⁴t₀²χ` vs `a²/t₀`: 16-cell flat to ~1–2, O(a²) coeff ≈ 0, common continuum ≈ 6.7 | 6 points, `a²/t₀ = 0.60–2.64`: O(a⁴) fit **c₀ = 6.60(42)**, curvature −0.577(150) vs slide's −0.59; freed **c₂ = −0.51(1.39) ≡ 0** vs cubic **−2.70(51)** (5.3σ); factor 2.8(4) less suppression at `a²/t₀ = 1.5` | ✅ within 7–12 % errors |
| (Fig 1) | Q(t₀) histogram, integer peaks | integer peaks at the fine point, washed out at `a²/t₀ = 2.6`, symmetric, no freezing | ✅ qualitative |
| 11 | `D₀ = (1/6)Σ₂₄ γ·n ∇`, Wilson term, 1 physical mode | implemented (`hcfree`, `hcwilson`); note: slide's `a(r/6)Σ∇*∇` is a factor-2 slip vs its own `arp²/2` statement — we use the latter (see §3) | ✅ |
| 12 | free pressure series: cubic O(1/N_t²), 16-cell leading **O(1/N_t⁴)** | fitted coefficients: cubic c₂ = 1.686996 (248/147 to 0.005 %), c₄ = 4.35(21) (635/147); 16-cell c₄ = 0.1295918 (127/980 to 7 digits), c₆ = 0.01779(5) (73/4158); 16-cell **c₂ = (−2.1±2.8)·10⁻⁸ ≡ 0** | ✅ exact |
| 13 | p/p_cont vs N_t: 3.83 vs 1.07 at N_t=4 | **3.82813 vs 1.07191** (283 % vs 7.19 %); full curve reproduced | ✅ exact |
| 14 | free Wilson–Dirac spectrum: compact near-ellipse, Re ≤ 16/3, \|Im\| ≤ 1.47 | max Re λ = **16/3** exactly, max \|Im λ\| = **1.467890 = 3^{1/4}(1+√3)/√6** (closed form found); full 16⁴ scatter plotted | ✅ exact |
| 16 | quenched λ₀: 16-cell arc lower & ~4× narrower than cubic at a≈0.12 fm | matched-radius smearing: ⟨λ₀⟩ **0.0347(11) vs 0.0999(14)** (2.9× lower), width **3.6× narrower** (paper ≈3× / ≈4×); 66+8 vs 60 configs | ✅ preliminary |
| 18 | `Q_Dirac − Q_flow` much narrower on 16-cell | not resolvable at (0.9 fm)⁴: both tight (RMS 0.185 vs 0.183); honeycomb flow-Q plateau 5× flatter; needs cluster-scale volume | ◐ inconclusive (volume-limited) |
| 19 | real-mode chirality piles at ±1 on 16-cell | matched radius: mean\|χ\| **0.987**, 100 % \|χ\|>0.9, vs cubic 0.871 / 62 % | ✅ preliminary (small subset) |
| 21 | 6× cost/site-volume; ~10× net HMC gain | arithmetic verified (FORMULATION §6); measured cost ratios consistent (hcForce 49 ms vs cubic staples; 2 sites × 3 neighbours) | ✅ |

## 2. The validation chain (why the numbers can be trusted)

Every layer was gated on an independent check before anything was built on it:

1. **Geometry** (`tgeom`, 42 checks): 24/32/96/16 combinatorics, clover reconstruction
   identity, point groups — all exact.
2. **Gauge layer** (`tgauge`): 16-way diagonal shift bit-exact vs coordinate gather;
   triangle sum vs an independent hand-rolled 3×3-matrix reference to 2e-18; gauge
   invariance to 9e-19.
3. **Action/force** (`taction`): finite-difference force check 9e-12 (per link kind);
   β-normalisation vs QEX's own Wilson action on the same continuum field.
4. **Monte Carlo** (`thmc`): reversibility 1.5e-11; ⟨e^{−ΔH}⟩ = 1; KP heatbath vs exact
   SU(2) moments; **HMC and heatbath agree at 0.46σ** on ⟨triangleSum⟩.
5. **Flow** (`tflow`): harness validated on QEX's cubic flow (rate/p² = 1.000000,
   artifact −1/12 exact); honeycomb constant **cflow = 6 pinned** (5.9999771(14e-4)).
6. **Topology** (`ttopo`): weak-field F̂ site-wise vs exact; **Atiyah–Singer
   constant-flux test**: Q/2n₁n₂ → 1 with pure 1/L⁴ artefacts — pins the prefactor
   *and* the ½ site-volume factor (an error would read exactly 2×).
7. **Dirac operator** (`twilson`): entrywise vs the momentum-space blocks at 4e-15;
   gauge covariance 2e-14; γ₅-hermiticity bit-exact; 25-site locality.
8. **Stout/clover** (`tstout`, `tclover`): stout heat-kernel constants κ_cubic = 1,
   **κ_hc = 1/3 exact**; clover matrix elements on constant flux to 2e-14 including the
   *predicted* artefact factors; self-dual flux nulls exactly one chirality sector.
9. **Eigensolver** (`tarnoldi`): Krylov–Schur vs exact Wilson spectra to 1e-13;
   shift-invert ≈ 220 applies/pair.
10. **Analysis** (`tanalysis`): t₀ finder, jackknife, fits on synthetic data.
11. **Cubic end-to-end reference** (task C): the entire flow→t₀→Q→χ→extrapolation chain
    run on the ordinary lattice with existing QEX code reproduced the literature
    (t₀/a² within 2–6 % of Necco–Sommer; `10⁴t₀²χ → 5.91(68)` vs 6.67(7)).

## 3. Things we learned that are not (clearly) in the paper

* **Slide 11's Wilson term normalisation is internally inconsistent by a factor 2**
  (`a(r/6)Σ∇*∇` vs "`arp²/2` as usual").  The `arp²/2` reading is correct — fixed
  independently by `max Re λ = 16/3` matching slide 14.
* **Closed forms**: max \|Im λ\| = `3^{1/4}(1+√3)/√6 = 1.467890`, attained on the
  (1,1,1,1) diagonal and on-axis at `p = 2 arccos((√3−1)/2)` (not at 2π/3).
* **The hexagon clover is dramatically better than the cubic clover**: its topological
  artefact at size L equals the cubic one at 2L *digit for digit*
  (`s_hc(L) = s_cubic(2L)`; both 1 − O(1/L⁴), ratio 16) — the "2× coarser lattice"
  claim visible inside a single measurement operator.
* **Flow and stout normalisations**: the honeycomb flow constant is exactly **1/6**
  (cflow = 6 in our convention), and one MP-stout step with plain triangle staples is a
  flow step of **ρ/3** (vs ρ on the cubic lattice).  Hence "ρ = 0.05, 6 steps" is a
  ~3× smaller smearing radius on the honeycomb if taken in the plain staple
  convention — a genuine convention ambiguity when comparing to the paper.
* **The β window**: with the triangle action, the slide-10 scaling window sits at
  β = 6.9–7.2 (⟨triangleSum⟩ ≈ 0.45–0.52); the β ≈ 7.0 crossover is hysteresis-free.
* **Sampler structure**: triangles have 1 axis + 2 different-index diagonal edges, so
  all 8 axis link fields can be heatbath-updated simultaneously and each of the 16
  diagonal fields sweeps whole — 17 staple builds per full sweep.

## 4. What a cluster run adds (no code changes needed)

* points below `a²/t₀ ≈ 0.5` (needs L ≥ 16 cells) → the O(a²)-vs-O(a⁴) discrimination
  without the external continuum anchor;
* per-point χ errors at the few-% level (needs 10⁴–10⁵ decorrelated Q measurements;
  watch τ_int(Q) which grew from 2 → 43 updates over our β range);
* the paper-scale fermion statistics (1215 configs, 12⁴, 300 eigenvalues/config);
* multi-rank: the code paths are MPI-clean by construction (QEX shifts/gathers), but
  only single-rank was exercised here — validate `tgauge`/`taction` on ≥ 2 ranks first.

## 5. Build & run quickstart

```bash
export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
cd build_mac
make src/experimental/honeycomb/tests/tgeom.nim && ./bin/tgeom          # geometry
make src/experimental/honeycomb/hcFreePressure.nim && ./bin/hcFreePressure  # slides 12-13
make src/experimental/honeycomb/hcPureGauge.nim                        # generation
make src/experimental/honeycomb/hcMeasFlow.nim                         # flow+t0+Q
# full test suite: tgeom tfree tanalysis tgauge taction thmc tflow ttopo twilson tstout tclover tarnoldi tspectrum
```

Use `OMP_NUM_THREADS=4` on shared machines (QEX spin barriers livelock when
oversubscribed).  Lattice sizes with the default `vlen:4` SIMD layout: L ∈ {4, 8, 12,
16, 20, 24}; other sizes need `-simdlen:1`.
