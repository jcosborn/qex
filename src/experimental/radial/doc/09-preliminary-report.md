# Preliminary report — QED3 in radial quantization on S²×ℝ, reproduced in QEX

> **Convention notice (2026-08-26):** Tier-1/massless results remain directly applicable.
> The Tier-2 files under `output/radial/t2` were generated with the retired additive-mass
> convention and are preserved as legacy results. The active standard-overlap campaign writes
> to `output/radial/t2-standard-overlap` and must be generated afresh; do not combine the two.
> **As of 2026-09-09 that directory has never been created: every Tier-2 number in §3 is
> legacy data.** The m=0 ensembles are convention independent (D(0) is the same operator; the
> Hasenbusch mass 0.5 is a preconditioner); only the condensate rows depend on it.

> **Review corrections (2026-09-09), applied throughout:** see §6 at the end. The most
> consequential: the σ_FS and axial "exact" statements of §3.1 were Ginsparg–Wilson tautologies of
> the connected contraction; the paper's (IV.12) is *not* printed inverted; the flat area identity
> was never in the paper's (IV.2); configuration counts below are now the numbers actually measured.

2026-08-21, corrected 2026-09-09. Code: `src/experimental/radial/` (worktree
`qed3-slides-reproduction-plan-0e70b6`, branch `experimental/radial-quant`, three commits on
devel plus the review commits). Targets: the Lattice 2026 talk (doc/01) and arXiv:2510.03085
(doc/02). Everything below is reproducible: 13 test suites (`tests/t*.nim`, all green at the
review), deterministic free-limit apps (`rfree`, `rspec`), and a checkpointed HMC campaign
(`campaign/t2.sh`). **Restart instructions: doc/08-handoff.md.** Audit trail: doc/06-status.md.

Statistics disclaimer, per the user's direction: Tier 2 ran at *preliminary* statistics
(1 to 19 measured configurations per ensemble, one L=2 attempt cut before thermalization) to
establish signal-over-noise; production statistics are a cluster campaign (same script,
raise the targets).

---

## 1. Tier 1 — the free-limit paper (arXiv:2510.03085): **reproduced, deterministic** (except the n_max rows, below)

| quantity | published | ours | deviation |
|---|---|---|---|
| fermion Δ₀ (L=1, L_t=168, T=16) | 0.953918 | **0.953918** | −1.4e−7 |
| fermion Δ₀^cont (V.7 fit, paper's grid) | 0.999998(34) | **0.999999(44)** | exact = 1 |
| gauge Δ₀ (L=1, L_t=120, T=16) | 1.33242 | **1.332430** | +7.5e−6 |
| gauge Δ₀^cont | 1.41409(18) | **1.414208(70)** | exact = √2, dev −5.5e−6 |
| Table I (doubler condition), 20 cells | — | **all 20 reproduced** | exact |
| Fig 10 propagator T-fold: Wilson / overlap | violated / protected | 0.87 / **8.8e−14** | 13 orders |
| Figs 4–12 | — | remade, TSV+PNG in `output/radial/free/` | visual agreement incl. normalization |
| slide-8 free legend, L=2 / L=4 | 1.010 / 0.965 | 1.010466 / 0.964330 (flat-κ) | see §3 |
| slide-8 free legend, L=1 | 1.154 | 1.158156 (flat-κ, a_t=0.1333) | the L=1 panel used the *free-limit paper's* a_t, not the campaign's 0.2 |

Fig-7/Fig-11 normalizations agree with the analytic correlators including the overall constant
(max rel. dev. in t∈[2,6]: 9.9→1.6 % and 55.6→1.2 % for L=1→8, pure O(a²)).

**Not reproduced (reported, not tuned):** the Eq. (V.9) n_max integers and the gauge residuals
at L=4, 8. Ours (least-squares, any weighting): fermion 3/7/14/28, gauge 3/7/16/30; published
6/10/19/32 and 3/8/18/35; gauge residual/DOF at the published n_max 0.0049/0.0063 vs the
published 0.0031/0.0037. The published fermion *residual column* equals our relative
residual/DOF evaluated **at their n_max** (0.0279 vs 0.028 at L=1; 0.0040 vs 0.0039 at L=4), so
the correlators agree; the paper states only "a least-squares fit to determine the cutoff n_max
that minimizes the residual" (DOF = L_t − 2), and the norm, weighting and the treatment of t=0
are unstated, so the selection rule cannot be adjudicated. Gauge L=1 matches exactly (3 @ 0.0031).

Cross-validation: every headline number was obtained by ≥2 independent implementations
(Nim per-Matsubara-mode dense pipeline; Nim real-space CG pipeline; a pure-Python oracle with
independently constructed geometry and spin connection — the oracle scripts lived in /tmp and are
gone; their results are transcribed in doc/06).

## 2. The coupling convention, as settled

The paper's couplings use the **exact spherical kite area**
\(A_e=\sum_\pm 4\arctan[\tan(\ell/4)\tan(\ell^*_\pm/2)]\) in both
\(\kappa_e=2A_e/(\bar a_s\ell)\) and \(\beta_\ell=2A_e/(g^2\ell^2a_t)\). That is in fact how the
paper *defines* \(A_{y_1y_2}\) (Sec. III: the sum of two spherical triangles); the flat form
\(\tfrac12\ell(\ell^*_1{+}\ell^*_2)\) is the form in which the paper's per-prism decomposition
(IV.4) and gauge derivation (IV.33)/(IV.35) are written and which it calls "equally possible" at
O(a²). An earlier version of this report attributed the flat identity to Eq. (IV.2) itself; that
was our transcription error. What the reproduction settled is the paper's own O(a²) ambiguity:
both published Δ₀ values pin the exact form to six digits (flat gives 0.921250 and 1.356697 —
off by 3.4 % and 1.8 %; the fermion flat value is from the Python oracle and is not on file), and
only the exact diamonds tile the sphere (ΣA_e = 4π to 1e−12; flat misses by 3.6 % at L=1).
Exception: the slide-8 Wilson-spectrum legends match the *flat* convention — evidently older
diagnostics.

Also derived en route: the generalized-eigenvalue weight of (IV.12) is
\(\overline{\delta V}/\delta V\), which is what the paper prints (an earlier claim that "(IV.12) as
printed is inverted" was a misreading of the generalized problem; the code was always right); and
the free gauge tower for slide 14's right panel is \(\Delta_{\ell}=\sqrt{\ell(\ell+1)}\) from
(C.37)/(V.14), so the continuum free ratio is √3 = 1.732, while the slide draws its free line
near 1.22 (= √6/2 = Δ₂^free/Δ_A is a plausible explanation, not an established one).

## 3. Tier 2 — interacting system (preliminary statistics)

Ensembles measured (L=1, N_f=2, a_t=0.2, M=1, exact-area convention, rational orders 31/11,
**legacy additive-mass binaries**): `L1g15m00` g²R=1.5 m=0 (17 configurations saved, of which
10 measured for scalars/gluonic/Wilson spectra and 11 for currents; the measurement job stopped
mid-configuration), `L1g10m00` g²R=1.0 m=0 (21 saved, 19 measured; currents unmeasured —
expensive), `L1g15m01` μ=0.1 (15 cfgs, condensate), `L1g15m04` μ=0.4 (1 cfg), plus pure-gauge
`pureL1`, `pureL2` (256 cfgs each, exact heatbath). HMC health across all (post-warmup, from the
hmc logs): acceptance 97–98 %, mean |ΔH| 0.02–0.03 with maxima up to 0.09, ⟨e^{−ΔH}⟩ = 1.000(3),
kernel windows inside.

### 3.1 Exact (statistics-independent) results — the deck's structural claims
| claim (slide) | result |
|---|---|
| ℓ=1,2 multiplets protected by I_h (13) | the *exact free* correlator matrices are ∝ 𝟙 to 1e−9 and the 60-rotation group average of one configuration is ∝ 𝟙 to 1e−16 (a Schur identity); ℓ=3 splits exactly 3+4 (T₂+G). The ensemble-mean residual from a pure block is what rmeas now reports (`l{1,2}_block_residual`); it was not measured on the dynamical data |
| σ_PS and σ_FS identical spectra (16) | the **connected** correlators are identical at every dt ≠ 0 at m=0 as a Ginsparg–Wilson identity of the contraction, and the dt=0 contact cancels to 3e−14. This is an algebraic consequence of `(1−D†)S† = −S`, **not** a check of the slide: σ_FS is the flavor singlet, and its fermion-disconnected (hairpin) piece was never computed in the legacy data. rmeas now measures it (`Delta_FS_full`, `fs_hairpin_fraction`); no dynamical number exists yet |
| conserved current (7, 11) | Ward: charge plateau flat to 1.6e−7 on dynamical configs; insertion = i×propagator jump ✓ (an exact gauge-variance statement at any mass, not conservation of the current–current correlator) |
| m=0 condensate | exactly 0 by Ginsparg–Wilson (measured 1e−14, smoke ensemble) — slide 10's content is entirely in m>0 |

### 3.2 Slide 8 — Wilson spectrum on dynamical configurations (T2.1)
min|D_W − 1|, ours (exact-κ, 19 and 10 cfgs) vs published legend (flat-κ, 3 cfgs):

| | g²R=1.0 | g²R=1.5 |
|---|---|---|
| ours | 0.8394(63) | 0.6606(171) |
| published | 0.814 | 0.682 |

The **additive mass shift** (monotone decrease with g²R) is clearly resolved; the ~3 % offsets
have the size and sign pattern of the flat-vs-exact κ convention difference. **Signal: yes.**

### 3.3 Slide 10 — condensate (T2.3)
⟨σ_PS⟩ (our normalization: per site, per 2-component flavor pair, contact-subtracted; the deck's
normalization is undefined on the slide and differs by an overall constant). **These rows are
legacy additive-mass data: the mass column is the additive μ of `D_ov + μ`, and the
standard-convention mass that labels the same operator is m = μ/(1+μ/2).**

| additive μ | standard m | ⟨σ_PS⟩ | configs |
|---|---|---|---|
| 0 | 0 | 0 exactly (GW) | — |
| 0.1 | 0.0952 | 0.004604(44) | 15 |
| 0.4 | 0.333 | 0.01347 (provisional) | 1 |

Nonzero, small, decreasing to zero with m — **consistent with no SSB** (the deck's conclusion);
2τ ≈ 1 (well-decorrelated). The m=0.2, 0.3 points, real statistics for m=0.4, and the whole scan
in the standard convention are the `cond_scan` entries of `t2.sh` (a fresh campaign; nothing
resumes). **Signal: yes; slope precision: cluster.**

### 3.4 Slides 14–15 — gluonic sector (GEVP over loop shapes at flow time s=0.6)
| quantity | pureL1 (free, 256 cfgs) | pureL2 (free, 256) | g²R=1.0 (19) | g²R=1.5 (10) | reference |
|---|---|---|---|---|---|
| Δ_{F,ℓ=1} | 1.204(78) | 1.375(84) | 0.99(20) | 2.1(1.0) | continuum free √2 = 1.414 (the exact L=1 lattice value is O(a²) below it; tmeas prints it as the `jtopCorrExact` reference); CFT: 2 |
| **R_{F,ℓ2/ℓ1}** | **1.633(145)** | 1.423(168) | 1.62(44) | 0.98(50) | continuum free √3 = 1.732; CFT 3/2 |
| R_{F²/F} | 1.418(135) | 1.644(168) | 1.76(45) | 0.88(54) | CFT large-N_f: 2 |

The free-ensemble ℓ2/ℓ1 ratio is 0.7σ from √3 and 0.9σ from 3/2 and does not discriminate them;
it is 2.8σ from 1.225, the value drawn on the slide. The correct comparison for a finite-L
pure-gauge ensemble is the exact lattice ratio at the same L from `jtopCorrExact`, which is not
yet tabulated. Dynamical points exist but are trend-level at ≤19 configs; the deck's Δ_F/Δ_A
panel needs Δ_A (below). **Signal: free yes; interacting: cluster.**

### 3.5 Slide 16 — scalars
Connected-only Δ_PS = 2.284(21) at g²R=1.0 and 2.272(37) at g²R=1.5 (local-effmass estimator at
t=1; free value 2); the connected Δ_FS is the same number by the GW identity of §3.1, so it is
not an independent measurement. Combined with the deck's Δ_A ≈ 2.2–2.4 the connected Δ_PS sits
in the published Δ_PS/Δ_A ≈ 0.88–0.98 band. The full σ_FS spectrum, including its hairpin, has
not been measured on any ensemble (rmeas now produces it). Our own Δ_A did not resolve (below),
so the ratio itself is a cluster deliverable. **Signal: yes (the absolute connected Δ_PS).**

### 3.6 Slides 11–13 — current spectroscopy: **not resolvable at this size** (honest failure)
The stochastic axial/vector correlators at ≤17 configs × few noise sources give NaN plateau fits
(Δ_A, Δ_V, the ℓ-resolved ratios). This is a statistics wall, not a code defect (the same
estimators pass their exact tests, and the deck itself flags the disconnected noise on slide 12).
Requirements for the cluster pass are estimated in doc/06 (≥100–300 configs, more noise hits,
and/or point-source dominance).

### 3.7 Slide 9 — gradient flow scale scan
Pure-gauge flow curves (256 cfgs, L=1 and L=2, flow times 0–1.6) are in
`output/radial/t2/pure*/analysis/flow.tsv`; dynamical E_s(t) measured on both m=0 ensembles.
Known open item (doc/06 WP-G): for pure gauge the nine-curve g²-splitting of the slide must come
from the fermions, and raw E_s collapses across L while E_s·√L does not — plot both when the
dynamical statistics exist.

## 4. Costs (measured) and the cluster plan
41.5 s/trajectory at L=1, n_t=60, N_f=2 (serial, this Mac; 3 concurrent runs → ~2.5–4×).
Measurement: currents 115 s, condensate 21 s (8 noises), scalars 28 s, gluonic 0.03 s,
wspec 3.9 s per config. L=2 ≈ 5–8× L=1 per trajectory. The full deck grid
(g²a ∈ {0.5,1,1.5} × L ∈ {1,2,4} × N_f ∈ {2,4,6}, ≥10³ trajectories each) is a
straightforward ensemble-parallel cluster campaign with the existing `t2.sh` (raise the
targets; re-enable nf4/nf6; add L=4). No code changes required; the code is single-rank by
design and parallelizes over ensembles.

## 5. Deliverables index
`output/radial/free/` (Tier-1 TSVs + 9 figure PNGs) · `output/radial/t2/<ens>/{meas,analysis}/`
(legacy additive-mass data; old TSV format, rejected by the current rmeas) · `doc/01–08` · tests
`src/experimental/radial/tests/` · apps `rgeom rfree rspec rgauge rhmc rmeas` · campaign
`campaign/t2.sh` (+ `campaign/free/*.gp`).

## 6. Review corrections (2026-09-09)

A full review of the code, outputs and documents (own derivations plus paper verification
against the published PDF) found the following; all are applied above and in the code.

1. **σ_FS hairpin.** The singlet scalar's fermion-disconnected piece was never computed; the
   "identical correlators, stronger than the slide" claim was a GW identity of the connected
   contraction. Fixed: `scalarSample`/`scalarConn`/`scalarOnePoint` and the `Delta_*_full` rows.
2. **Axial hairpin.** The block-τ₃ current has the one-point function −2i Im tr[KS] per pair; the
   doc/07 argument from tr τ₃ = 0 assumed a flavor-proportional propagator, which diag(S, S†) is
   not. Fixed: `Delta_A_full_l1` alongside the connected rows.
3. **Vector hairpin weight.** rmeas weighted the disconnected piece by N_f; per-pair normalization
   requires N_f/2. Fixed.
4. **ℓ=3 splitting.** The per-m spread is frame dependent; the invariant T₂/G block analysis with
   the icosahedral projectors is now the reported quantity (`Delta_A_conn_l3_{T2,G}`,
   `l3_block_split`), with `l{1,2}_block_residual` as the protection check.
5. **Paper attributions.** (IV.12) is printed correctly; the flat area identity is not in (IV.2);
   (B.4) has η=1 in both cases; the paper's Fig. 10 text/caption disagree on L. Corrected in
   doc/02, doc/06 (dated note), doc/08 and here.
6. **Polyakov mode.** Freezing the uniform temporal mode with fermions is a physical choice
   (the temporal twist), not gauge fixing; documented in doc/02 §5 and hmc/trajectory.nim.
7. **Bookkeeping.** Configuration counts, the additive-μ labels of §3.3, the never-created
   standard-convention output directory, the unsourced "1.3299", test counts and branch names.
8. **Code structure.** One generic CG (`ops/solve.nim`) for spinor, gauge and site-scalar
   fields with reusable scratch; one dense-helper module (`core/dense.nim`); one FNV hash; one
   `levels` ladder; one effective-mass routine; shared test fixtures (`tests/thelpers.nim`); rmeas
   uses the tested estimators of `meas/observables` instead of private copies; kernel-window
   violations stop rmeas; the summary is written atomically; Metropolis hooks take the sampler by
   value (devel's `hmc/metropolis.nim` API).
