# Quenched Wilson–Dirac spectra — preliminary results (slides 16/18/19, paper Figs. 5–8)

**Status: PRELIMINARY.**  The measurement agent (task D4) was stopped on user request
before its planned statistics completed; the lead harvested the logs as they stood
(procedure: HANDOFF.md §3; join fix for the honeycomb FLOWQ column applied — the raw
`fig7_hc*.dat` were regenerated, see `data/hcall.log` note below).  A cluster-scale
rerun is planned.  Figures: `fig{5,6,7,8}_*.png`; raw data `data/*.dat`,
`data/{hcall,cub586}.log`; drivers `hcSpectrum.nim` / `cubicSpectrum.nim`;
`harvest.py` (NB: its in-file hc Q_flow join was the bug the agent was fixing when
stopped — the numbers below use the corrected external join).

## Setup

| | cubic | 16-cell honeycomb |
|---|---|---|
| β | 5.86 | 7.22 |
| volume | 8⁴ sites | 8⁴ cells (2·8⁴ sites) |
| t₀/a² (ensemble) | 2.534 (a²/t₀ = 0.395) | per-config ⟨t₀⟩ ≈ 1.8–2.3 |
| configs harvested | **60** | **66** (ρ=0.05) + **8** (ρ=0.15 subset) |
| operator | Wilson + tree-level clover c_SW = 1, m = 0, r = 1, antiperiodic time |
| smearing | 6 stout, ρ = 0.05 | 6 stout, ρ = 0.05; subset ρ = 0.15 |
| eigenvalues | 16 lowest per config (shift-invert Krylov–Schur), residuals ≲ 6e-6 |
| cost | 12.6 s/cfg | 60.4 s/cfg (4 threads) |

**Smearing-radius caveat (important).**  In the plain Morningstar–Peardon staple
convention used here, one stout step is a flow step of ρ (cubic) but only ρ/3
(honeycomb) — task D2's measured κ.  So "ρ = 0.05 on both" (the paper's stated
parameters, taken literally) smears the honeycomb 3× *less*; the ρ = 0.15 subset is
the matched-radius run.  The results below strongly suggest the paper's ρ = 0.05
means matched radius (our ρ = 0.15).

## Results

### Re λ₀ — additive mass fluctuation (slide 16 / Figs. 5–6)

| ensemble | ⟨Re λ₀⟩ | width σ |
|---|---|---|
| cubic | 0.0999(14) | 0.0109 |
| 16-cell, ρ=0.05 | 0.1625(8) | 0.0068 |
| **16-cell, ρ=0.15 (matched radius)** | **0.0347(11)** | **0.0031** |

At matched smearing the honeycomb λ₀ sits **2.9× lower** and its distribution is
**3.6× narrower** than cubic — quantitatively matching the slide (paper: ≈0.031 vs
≈0.09, ~4× narrower), with only 8 configs in the subset.  Even at the 3×-weaker
ρ=0.05 the honeycomb width is 1.6× narrower, though the mean is then *larger* than
cubic — direct evidence that the smearing convention, not the lattice, drives the mean,
and that the paper's convention is matched-radius.

### Chirality of real low modes (slide 19 / Fig. 8)

| ensemble | n real modes | mean \|χ\| | frac \|χ\| > 0.9 |
|---|---|---|---|
| cubic | 29 | 0.871 | 0.62 |
| 16-cell, ρ=0.05 | 53 | 0.916 | 0.77 |
| **16-cell, ρ=0.15** | 9 | **0.987** | **1.00** |

The matched-radius honeycomb modes pile sharply at ±1 exactly as slide 19 shows;
the cubic distribution is broad.  (Statistics on the subset are small; the trend
ρ=0.05 → 0.15 is monotone.)

### Q_Dirac − Q_flow (slide 18 / Fig. 7)

| ensemble | RMS(Q_D−Q_f) | frac in [−0.5, 0.5] | flow-plateau spread (max) |
|---|---|---|---|
| cubic | 0.183 | 0.95 | 0.258 |
| 16-cell, ρ=0.05 | 0.185 | 0.97 | **0.051** |
| 16-cell, ρ=0.15 | 0.153 | 1.00 | 0.022 |

**The dramatic slide-18 contrast is NOT resolved at our parameters** — both lattices
already agree well — because our volumes (≈(0.9 fm)⁴, sd(Q) ≈ 0.6–0.9) are far
smaller than the paper's 12⁴ at a≈0.12 fm (sd(Q) ≈ 2–3): with so little topology,
integer-matching is easy on both lattices; also our cubic point is slightly finer
(a²/t₀ = 0.40 vs the paper's ≈0.52).  What *is* visible: the honeycomb flow-Q plateau
is 5× flatter (last column), consistent with task W's finding that the hexagon-clover
topological artefacts are 16× smaller.  A cluster run at 12⁴-scale volumes is
expected to expose the histogram contrast.

### Verdict vs the slides

* **slide 16 (λ₀ mean & width): reproduced** at matched smearing (2.9× / 3.6× vs the
  paper's ≈3× / ≈4×), preliminary statistics.
* **slide 19 (chirality pileup at ±1): reproduced**, small subset.
* **slide 18 (Q_Dirac−Q_flow narrowing): not resolvable at (0.9 fm)⁴ volumes** —
  both already tight; indirect indicator (flow-plateau flatness, 5×) favours the
  honeycomb.  Needs the cluster-scale volume.

## Reproduction notes

* `data/hcall.log` = `hc722.log` (configs 1–44) + `hc722b.log` (45–66, plus the
  ρ=0.15 subset); the leftover generation jobs may have appended more configs after
  the harvest cut — the numbers above used exactly n = 66/60/8.
* Q_flow per config: honeycomb from the `FLOWQ` line (Q at the per-config t₀; three
  plateau samples recorded), cubic from the driver's in-line join at the ensemble t₀.
* Rerun everything: `run_gen_cubic.sh`, `run_cubic.sh`, `run_hc.sh`,
  `run_hc_rho15.sh`, then `harvest.py` (fix its hc FLOWQ join first, or reuse the
  external join snippet recorded in STATUS.md task-D4-closeout).
