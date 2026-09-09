# Results: reproducing "QCD on the 16-cell honeycomb"

Target: S. D. Katz, D. Nogradi, arXiv:2512.10604 (6 pages, 9 figures) and the
22-slide Lattice 2026 talk.  Every number below comes from code in this tree
or from the archived logs of the first campaign (`archive/plots/`, re-analysed
by `archive/analysis/`).  All Monte-Carlo results are laptop-scale
(hundreds of configurations per point against the paper's O(50 000)); they
establish correctness of the implementation and the visibility of the signals,
not the paper's precision.

## 1. Verdict by slide

| slide | paper content | this tree | verdict |
|---|---|---|---|
| 4 | point groups 384 / 1152 | exhaustive enumeration: 384 / 1152 | exact |
| 7–9 | D₄* lattice, 24 neighbours, 12 links, 32 triangles, β = 2N/g² unchanged | geometry tests; classical limit S₁₆/S_cubic = 1 + 0.056p² → 1 | exact |
| 10 | 10⁴t₀²χ vs a²/t₀ flat, O(a⁴) fit, O(a²) coefficient consistent with zero, continuum ≈ 6.77 | six points at a²/t₀ = 0.60–2.64; the five inside the paper's range sit on the paper's O(a⁴) curve (χ² = 2.7 for 5 points, no free parameter); O(a⁴) fit 6.43(32) (all) / 6.68(44) (x ≤ 2.1) vs paper 6.771(40) and Cè et al. 6.67(7) | consistent; the O(a²)-vs-O(a⁴) shape is not discriminated by these data (see §3) |
| 11 | D₀, Wilson term, one physical mode | implemented; the paper's eq. (4) has a factor-2 and sign slip against its own a r p²/2 statement | exact |
| 12 | free pressure series | all four coefficients derived analytically from the dispersion relation and reproduced numerically to 7 digits | exact |
| 13 | p/p_cont: 283 % vs 7 % at N_t = 4 | 3.828128 vs 1.071914, from an exact transfer-matrix closed form | exact |
| 14 | free spectrum: Re λ ≤ 16/3, Im λ ≤ 1.47 | closed forms 16/3 and 3^{1/4}(1+√3)/√6 = 1.467890 | exact |
| 16 | λ₀ mean lower and narrower on the honeycomb (factor ≈ 3 in the paper) | matched smearing, 60 configs each: ⟨Re λ₀⟩ 0.0364(7) vs 0.0999(14), 2.75× lower; width 0.0055 vs 0.0109, 2.0× narrower | reproduced at low statistics |
| 18 | Q_Dirac − Q_flow narrower | RMS 0.126 (100 % within ±½) vs 0.183 (95 %) | direction reproduced; small volumes make both tight |
| 19 | real-mode chirality piles at ±1 | mean\|χ\| 0.964, 93 % above 0.9, vs 0.871, 62 % | reproduced at low statistics |
| 21 | 6× cost per site volume | arithmetic and measured cost ratios | exact |

## 2. Validation chain (tests/)

1. `tgeom`: 24/32/96/16 counts, triangle paths, hexagon rings, the (3/8)
   reconstruction identity, point groups.
2. `tgauge`: 16-way shift tree bit-exact against a coordinate gather;
   triangle sum against a brute-force 3×3 reference; gauge invariance.
3. `taction`: finite-difference force per link class; β normalisation against
   QEX's own Wilson action on the same continuum field.
4. `thmc`: staple sums equal actionDeriv(2N); single-link action; reversibility
   for three integrators; mdevolve leapfrog equals an explicit TVT composition;
   ⟨e^{−ΔH}⟩ = 1; HMC and heatbath agree on ⟨triangleSum⟩; file round trip.
5. `tflow`: harness validated on QEX's cubic flow (rate = p̂², artefact −1/12);
   honeycomb constant 1/c_HC = 5.9996(3), pinned to the analytic 6.
6. `ttopo`: site-wise F̂ against the exact weak-field F; Atiyah–Singer
   Q/2n₁n₂ → 1 with pure 1/L⁴ artefacts, sensitive to the ½ site volume.
7. `twilson`: operator entrywise against the momentum-space blocks; gauge
   covariance; γ₅-hermiticity; antiperiodic time.
8. `tstout`, `tclover`: stout heat-kernel constants κ_cubic = 1, κ_hc = 1/3;
   clover matrix elements on the constant-flux background with the predicted
   artefact factors; chirality sector nulling on self-dual flux.
9. `tarnoldi`: Krylov–Schur against exact Wilson spectra and a dense
   non-normal 400×400 zgeev reference; shift-invert.
10. `tanalysis`, `tfree`, `tspectrum`: t₀ finder, jackknife, fits; free
    operator identities; the qDiracSign pin and the flux-background real modes.

Independent checks outside the tests: the geometry counts and point groups by
a brute-force script; the flow constant from a linearised 24×24 Bloch matrix
(slow eigenvalue λ/p² = 1/6 − 0.0023p², identical on all 24 link classes);
the pressure series by symbolic expansion (cubic 248/147 and 635/147 recovered
as a check, honeycomb a₂ = 0, ⟨a₄⟩ = −1/210, ⟨a₆⟩ = −1/11340); the finite-N_t
pressure by the transfer-matrix roots.

## 3. Topological susceptibility (slide 10)

Honeycomb ensembles: heatbath + 3 overrelaxation sweeps per update, Q = Q_flow
at the ensemble t₀ from the hexagon clover, χ = ⟨Q²⟩/V with V = N_cells.
Errors: delete-block jackknife with block ≥ 2τ_int of Q and of Q²; burn-in
10τ_int(Q).  Re-analysed from the archived series (`archive/analysis/`).

| β | cells | n | t₀/a² | a²/t₀ | ⟨Q²⟩ | 10⁴t₀²χ | τ_int(Q) [updates] | n_eff(Q) |
|---|---|---|---|---|---|---|---|---|
| 6.90 | 8⁴ | 522 | 0.3787(6) | 2.641(4) | 8.98(62) | 3.15(22) | 2.2 | 465 |
| 6.95 | 8⁴ | 390 | 0.5092(14) | 1.964(6) | 6.87(51) | 4.35(32) | 2.4 | 390 |
| 7.00 | 8⁴ | 382 | 0.6726(23) | 1.487(5) | 4.81(37) | 5.32(40) | 4.5 | 214 |
| 7.07 | 8⁴ | 341 | 0.968(6) | 1.033(7) | 2.92(26) | 6.68(56) | 9.6 | 107 |
| 7.15 | 12⁴ | 251 | 1.359(8) | 0.736(4) | 6.39(74) | 5.69(65) | 15.6 | 81 |
| 7.20 | 12⁴ | 284 | 1.653(11) | 0.605(4) | 4.89(65) | 6.44(82) | 48.8 | 35 |

The 8⁴ errors plateau with the block size; the two 12⁴ points do not (their
errors at blocks of 16–32 are 20–45 % larger than quoted), and n_eff(Q) = 35
at β = 7.20.  Every configuration's own t²E crossed 0.3; ⟨Q⟩ is consistent
with zero at every point; hot and cold starts agree (no bulk transition in
β = 6.8–7.3).

Comparison with the paper's Fig. 2, digitised from the PDF vector paths
(16-cell points at a²/t₀ = 0.50–2.0 with errors 0.03–0.04, continuum
6.771(40), curve y = 6.771 − 0.588x²):

| β | a²/t₀ | 10⁴t₀²χ | paper curve | pull |
|---|---|---|---|---|
| 6.90 | 2.641 | 3.15(22) | 2.67 (outside the paper's range) | +2.2 |
| 6.95 | 1.964 | 4.35(32) | 4.50 | −0.5 |
| 7.00 | 1.487 | 5.32(40) | 5.47 | −0.4 |
| 7.07 | 1.033 | 6.68(56) | 6.14 | +1.0 |
| 7.15 | 0.736 | 5.69(65) | 6.45 | −1.2 |
| 7.20 | 0.605 | 6.44(82) | 6.56 | −0.1 |

χ² = 2.7 for the five points inside the paper's range with no free parameter.
Since a wrong flow, E or Q normalisation would displace both axes, this
confirms the normalisation chain (β convention, cflow = 6, hexagon clover,
½ site volume, cell volume) end to end.

Fits, y = 10⁴t₀²χ, x = a²/t₀ (weighted least squares):

| form | points | c₀ | c₂ | c₄ | χ²/dof |
|---|---|---|---|---|---|
| c₀ + c₂x | all 6 | 7.78(48) | −1.74(22) | | 0.84 |
| c₀ + c₄x² | all 6 | 6.43(32) | | −0.478(61) | 0.77 |
| c₀ + c₂x + c₄x² | all 6 | 6.99(1.19) | −0.71(1.45) | −0.29(40) | 0.94 |
| c₀ + c₄x² | x ≤ 2.1 | 6.68(44) | | −0.598(156) | 0.79 |
| c₀ + c₂x | x ≤ 2.1 | 7.60(68) | −1.60(43) | | 1.07 |
| cubic c₀ + c₂x (β = 5.6, 5.7, 5.8) | 3 | 5.33(86) | −2.29(62) | | 0.30 |
| Cè et al. 1506.06052 | | 6.67(7) | | | |
| paper Fig. 2 refit | | 6.772(19) | | −0.588(9) | 0.31 |

What these data do and do not show: the pure O(a²) and pure O(a⁴) forms
describe the honeycomb points equally well, and the three-parameter fit's
parameters are 91–99 % correlated, so "c₂ consistent with zero" here is a
statement of low power, not a measurement.  The O(a⁴) form is singled out
only through the external continuum anchor (the O(a²) intercept is 2.1–2.3σ
above both references while the O(a⁴) intercepts agree with them) and through
the agreement of the points with the paper's curve.  The paper's 1σ bound on
c₂, an order of magnitude below the cubic slope, needs its precision.
Same-model comparisons of the two lattices: O(a²) slopes −1.74(22) vs
−2.29(62), 0.8σ apart; the qualitative contrast is that the honeycomb still
carries 5.3(4) at a²/t₀ ≈ 1.5 where the cubic lattice has 1.9(2).

Cubic reference (Wilson action, HMC, QEX `gaugeFlow`/`fmunu`/`topoQ`, clover
Q at the per-configuration t₀): β = 5.6, 5.7, 5.8 on 8⁴, 12⁴, 12⁴ are usable;
β = 5.9 (16⁴, 35 decorrelated configurations after burn-in, ⟨Q⟩ = −1.4(5))
and β = 6.0 (n_eff(Q) = 17) are not.  t₀/a² agrees with the Necco–Sommer
r₀/a parametrisation (valid for 5.7 ≤ β ≤ 6.92) combined with
t₀/r₀² = 0.1108(17) (Cè et al. eq. 5.2) within 2–6 %.

## 4. Quenched Wilson–Dirac spectra (slides 16, 18, 19)

Setup: tree-level clover c_SW = 1, m = 0, r = 1, antiperiodic time, 6 stout
steps, 16 eigenvalues nearest σ = −0.8 by shift-invert Krylov–Schur, real
modes |Im λ| < 10⁻⁵, Q_flow at fixed t = 1.917.  Honeycomb β = 7.22 on 8⁴
cells (ensemble t₀/a² = 1.73), cubic β = 5.86 on 8⁴ (t₀/a² = 2.53): the
spacings differ by 21 % and both volumes are ≈ (0.9 fm)⁴, far below the
paper's 12⁴ at a ≈ 0.12 fm.  With ρ = 0.05 in the plain staple convention the
honeycomb is smeared three times less than the cubic lattice; ρ = 0.15 is
the matched-radius run and reproduces the paper's λ₀, so the paper's
"same parameters" is read as equal smearing radius.

| ensemble | n | ⟨Re λ₀⟩ | sd(λ₀) | real modes | mean\|χ\| | frac \|χ\| > 0.9 | RMS(Q_D − Q_f) | frac \|ΔQ\| ≤ ½ |
|---|---|---|---|---|---|---|---|---|
| cubic, ρ = 0.05 | 60 | 0.0999(14) | 0.0109 | 29 | 0.871 | 0.62 | 0.183 | 0.95 |
| honeycomb, ρ = 0.05 | 66 | 0.1625(8) | 0.0068 | 53 | 0.916 | 0.77 | 0.185 | 0.97 |
| honeycomb, ρ = 0.15 | 60 | 0.0364(7) | 0.0055 | 43 | 0.964 | 0.93 | 0.126 | 1.00 |

Completeness of the real-mode count: the converged window reaches
σ + max|λ − σ| = 0.22 (cubic, minimum over configurations), 0.25 (honeycomb
ρ = 0.05) and 0.16 (ρ = 0.15); cubic real modes extend to the window edge
(7 of 29 above 0.15), honeycomb ρ = 0.15 modes do not.  With a common cut at
0.15 the cubic Q_Dirac changes on 5 configurations (RMS 0.226); the honeycomb
numbers are unchanged.  `spectrum` now reports the reach per configuration and
counts real modes only below `-recut`.  The cubic ensemble has
⟨Q_flow⟩ = −0.22(7), three standard deviations from zero: its topology is not
well sampled.

## 5. Free fermions (slides 12–14)

`freepressure` evaluates O = (p(T) − p(0))/T⁴ exactly (transfer-matrix
roots, Richardson-extrapolated 3D integral, grid error ≤ 10⁻⁹) and fits the
1/N_t series: cubic 1.686996 (248/147 to 5·10⁻⁵) and 4.35(21) (635/147);
honeycomb 0.1295918 (127/980 to 7 digits), 0.01779(5) (73/4158), c₂ =
(−2.1 ± 2.8)·10⁻⁸.  The same coefficients follow analytically from the
dispersion relation (FORMULATION §6).  `freespectrum` reproduces the 16⁴
scatter of slide 14 with the closed-form extremes.

## 6. Findings about the formulation

* cflow = 6 and stoutKappa = 1/3 are exact (FORMULATION §5), not calibration
  constants; the honeycomb flow's O(p²) artefact is 35 times smaller than the
  cubic one.
* The paper gives no β values, flow equation, E definition or clover
  normalisation; all conventions here are derived and then confirmed by the
  Fig. 2 comparison.
* The paper's eq. (4) Wilson term is a factor 2 (and a sign, depending on the
  meaning of ∇*) away from its stated a r p²/2; the free spectrum fixes the
  latter.
* The hexagon clover's topological artefact at size L equals the cubic
  clover's at 2L (projected triangle area ¼ of a plaquette).
* F₄ symmetry removes only the O(4)-breaking dimension-6 operator; the
  vanishing O(a²) term in t₀²χ is an empirical statement (FORMULATION §7).

## 7. Reproduction

```bash
# generation
hcpuregauge -geom:8,8,8,8 -beta:7.0 -algo:hb -nwarm:200 -ntraj:2000 -savefreq:5 -outdir:/tmp/hc70
refcubicgen -lat:8,8,8,8 -beta:5.7 -ntraj:2000 -savefreq:10 -outdir:/tmp/cub57
# flow, t0, Q, chi_top (ENSEMBLE line), either lattice
measflow -lattice:hc -cfgs:'/tmp/hc70/hc.*.lime' -bin:8
measflow -lattice:cubic -cfgs:'/tmp/cub57/cfg.*.lime' -bin:4
measflow -chifit:chitop.dat            # continuum fits of a table
# spectra, chirality, Q_Dirac vs Q_flow
spectrum -lattice:hc -cfgs:'/tmp/hc70/hc.*.lime' -rho:0.15 -sigma:-0.8 -nev:16 -recut:0.15
spectrum -lattice:cubic -cfgs:'/tmp/cub57/cfg.*.lime' -rho:0.05 -sigma:-0.8 -nev:16 -recut:0.15
# free fermions
freepressure -outdir:/tmp/free ; freespectrum -outdir:/tmp/free
```

What a production run adds: points below a²/t₀ ≈ 0.5 (L ≥ 16 cells) and
few-per-cent errors per point (10⁴–10⁵ decorrelated charges; τ_int(Q) grows
from 2 to 50 updates over β = 6.9–7.2), which is what the O(a²) bound
requires; fermion spectra at matched spacing on 12⁴ with the paper's 1215
configurations and 300 modes; multi-rank validation (all runs here were
single-rank).

Archive of the first campaign: `archive/plots/{honeycomb,cubic,fermions}`
(per-configuration series, logs, run scripts, plots), `archive/analysis/`
(the re-analysis scripts and tables behind §3–4), `archive/docs/` (the
session documents).
