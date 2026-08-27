# Preliminary report — QED3 in radial quantization on S²×ℝ, reproduced in QEX

> **Convention notice (2026-08-26):** Tier-1/massless results remain directly applicable.
> The Tier-2 files under `output/radial/t2` were generated with the retired additive-mass
> convention and are preserved as legacy results. The active standard-overlap campaign writes
> to `output/radial/t2-standard-overlap` and must be generated afresh; do not combine the two.

2026-08-21. Code: `src/experimental/radial/` (worktree
`qed3-slides-reproduction-plan-0e70b6`, branch `claude/qed3-slides-reproduction-plan-0e70b6`,
uncommitted). Targets: the Lattice 2026 talk (doc/01) and arXiv:2510.03085 (doc/02).
Everything below is reproducible: 12 test suites (~196 tests, all green), deterministic
free-limit apps (`rfree`, `rspec`), and a checkpointed HMC campaign (`campaign/t2.sh`).
**Current restart instructions: doc/08-handoff.md.** Full audit trail: doc/06-status.md.

Statistics disclaimer, per the user's direction: Tier 2 ran at *preliminary* statistics
(10–21 measured configurations per ensemble, one L=2 attempt cut before thermalization) to
establish signal-over-noise; production statistics are a cluster campaign (same script,
raise the targets).

---

## 1. Tier 1 — the free-limit paper (arXiv:2510.03085): **fully reproduced, deterministic**

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

**Known discrepancy (reported, not tuned):** the Eq. (V.9) n_max integers. Ours (least-squares,
any weighting): fermion 3/7/14/28, gauge 3/7/16/30; published 6/10/19/32 and 3/8/18/35. But the
published *residual column* equals our relative residual/DOF evaluated **at their n_max**
(0.0279 vs 0.028 at L=1; 0.0040 vs 0.0039 at L=4) — the correlators agree perfectly; only the
paper's unstated n_max selection rule differs. Gauge L=1 matches exactly (3 @ 0.0031).

Cross-validation: every headline number was obtained by ≥2 independent implementations
(Nim per-Matsubara-mode dense pipeline; Nim real-space CG pipeline; a pure-Python oracle with
independently constructed geometry and spin connection).

## 2. The physics discovery of the reproduction: the coupling convention

The paper's couplings use the **exact spherical kite area**
\(A_e=\sum_\pm 4\arctan[\tan(\ell/4)\tan(\ell^*_\pm/2)]\) in both
\(\kappa_e=2A_e/(\bar a_s\ell)\) and \(\beta_\ell=2A_e/(g^2\ell^2a_t)\) — **not** the flat form
\(\tfrac12\ell(\ell^*_1{+}\ell^*_2)\) that Eq. (IV.2) presents as an equivalent identity (it is
equivalent only at O(a²)). Evidence: both published Δ₀ values pin the exact form to six digits
(flat gives 0.921250 and 1.356697 — off by 3.4 % and 1.8 %); and only the exact diamonds tile
the sphere (ΣA_e = 4π to 1e−12; flat misses by 3.6 % at L=1). Exception: the slide-8
Wilson-spectrum legends match the *flat* convention — evidently older diagnostics.
Also derived and verified en route: the generalized-eigenvalue volume weight is
\(\overline{\delta V}/\delta V\) (the paper's (IV.12) as printed is inverted relative to its own
(IV.11)); and the free gauge tower reference for slide 14's right panel is
\(\Delta_{\ell}= \sqrt{\ell(\ell+1)}\), so the free ratio is **√3 = 1.732**, not the drawn
√(3/2) = 1.225 (which is Δ₂^free/Δ_A). §4.4 below confirms √3 numerically.

## 3. Tier 2 — interacting system (preliminary statistics)

Ensembles measured (L=1, N_f=2, a_t=0.2, M=1, exact-area convention, rational orders 31/11):
`L1g15m00` g²R=1.5 m=0 (17 cfgs, full observables), `L1g10m00` g²R=1.0 m=0 (21 cfgs; scalars,
gluonic, Wilson spectra; currents unmeasured — expensive), `L1g15m01` m=0.1 (15 cfgs, condensate),
`L1g15m04` m=0.4 (1 cfg), plus pure-gauge `pureL1`, `pureL2` (256 cfgs each, exact heatbath).
HMC health across all: acceptance ≈ 97 %, |ΔH| ≲ 0.05, ⟨e^{−ΔH}⟩ ≈ 1, kernel windows inside.

### 3.1 Exact (statistics-independent) results — the deck's structural claims
| claim (slide) | result |
|---|---|
| ℓ=1,2 multiplets protected by I_h (13) | correlator matrices ∝ 𝟙 at 1e−16 per configuration; ℓ=3 splits exactly 3+4 |
| σ_PS and σ_FS identical spectra (16) | **identical correlators at every dt**: max dev 0.0 exactly on dynamical ensembles; GW contact at dt=0 verified to 3e−14 (stronger statement than the slide's) |
| conserved current (7, 11) | Ward: charge plateau flat to 1.6e−7 on dynamical configs; insertion = i×propagator jump ✓ |
| m=0 condensate | exactly 0 by Ginsparg–Wilson (measured 1e−14) — slide 10's content is entirely in m>0 |

### 3.2 Slide 8 — Wilson spectrum on dynamical configurations (T2.1)
min|D_W − 1|, ours (exact-κ, 17–21 cfgs) vs published legend (flat-κ, 3 cfgs):

| | g²R=1.0 | g²R=1.5 |
|---|---|---|
| ours | 0.8394(63) | 0.6606(171) |
| published | 0.814 | 0.682 |

The **additive mass shift** (monotone decrease with g²R) is clearly resolved; the ~3 % offsets
have the size and sign pattern of the flat-vs-exact κ convention difference. **Signal: yes.**

### 3.3 Slide 10 — condensate (T2.3)
⟨σ_PS⟩ (our normalization: per site, per 2-component flavor pair, contact-subtracted; the deck's
normalization is undefined on the slide and differs by an overall constant):

| mR | ⟨σ_PS⟩ | configs |
|---|---|---|
| 0 | 0 exactly (GW) | — |
| 0.1 | 0.004604(44) | 15 |
| 0.4 | 0.01347 (provisional) | 1 |

Nonzero, small, decreasing to zero with m — **consistent with no SSB** (the deck's conclusion);
2τ ≈ 1 (well-decorrelated). The m=0.2, 0.3 points and real statistics for m=0.4 are queued in
`t2.sh` (resume per doc/08). **Signal: yes; slope precision: cluster.**

### 3.4 Slides 14–15 — gluonic sector (GEVP over loop shapes at flow time s=0.6)
| quantity | pureL1 (free, 256 cfgs) | pureL2 (free, 256) | g²R=1.0 (19) | g²R=1.5 (11) | reference |
|---|---|---|---|---|---|
| Δ_{F,ℓ=1} | 1.204(78) | 1.375(84) | 0.99(20) | 2.1(1.0) | free exact (L=1 lattice): 1.3299; CFT: 2 |
| **R_{F,ℓ2/ℓ1}** | **1.633(145)** | 1.423(168) | 1.62(44) | 0.98(50) | **free √3 = 1.732** ✓ (2.8σ from the slide's √(3/2)); CFT 3/2 |
| R_{F²/F} | 1.418(135) | 1.644(168) | 1.76(45) | 0.88(54) | CFT large-N_f: 2 |

The free-ensemble ℓ2/ℓ1 ratio **numerically confirms √3** and refutes the slide's free line —
closing doc/06 open question 1. Dynamical points exist but are trend-level at ≤19 configs;
the deck's Δ_F/Δ_A panel needs Δ_A (below). **Signal: free yes; interacting: cluster.**

### 3.5 Slide 16 — scalars
Δ_PS = Δ_FS = 2.284(21) at g²R=1.0 and 2.272(37) at g²R=1.5 (local-effmass estimator at t=1;
free value 2). Combined with the deck's Δ_A ≈ 2.2–2.4 this sits exactly in the published
Δ_PS/Δ_A ≈ 0.88–0.98 band. Our own Δ_A did not resolve (below), so the ratio itself is a
cluster deliverable. **Signal: yes (the absolute Δ_PS).**

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
· `doc/01–08` · tests `src/experimental/radial/tests/` · apps `rgeom rfree rspec rgauge rhmc
rmeas` · campaign `campaign/t2.sh` (+ `campaign/free/*.gp`).
