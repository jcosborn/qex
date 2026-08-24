# Implementation plan — work packages

Read [`02-formulation.md`](02-formulation.md) for the physics, [`03-targets.md`](03-targets.md) for
what counts as done, and [`04-interfaces.md`](04-interfaces.md) for the signatures you must
implement. Record progress in [`06-status.md`](06-status.md).

## Dependency graph

```
              WP-A geom/lattice ────┬──► WP-E wilson ──┬──► WP-F overlap ──┬──► WP-H hmc ──┐
                                    │                  │                   │               │
              WP-D spinor/solve ────┴──────────────────┘                   ├──► WP-I meas ─┼──► WP-L
                                                                            │               │   interacting
              WP-B zolotarev ──────────────────────────┘                    │               │   campaign
                                                                            │               │
              WP-G gaugeact/flow ◄── WP-A ─────────────────────────────────┘               │
                                                                                            │
              WP-C analytic ───┬──► WP-K free-limit campaign (Tier 1) ──────────────────────┘
              WP-J fit/stats ──┘
```

Wave 1 (no dependencies, start together): **A, B, C, D, J**
Wave 2: **E, G**
Wave 3: **F**
Wave 4: **H, I, K**
Wave 5: **L**

---

## WP-A — Geometry: icosahedral S² and S²×ℝ
Files: `core/types.nim`, `core/geom.nim`, `core/lattice.nim`, `tests/tgeom.nim`, app `rgeom.nim`.

1. Vec3/Mat2/Spinor value types and Pauli algebra.
2. Icosahedron → level-L refinement → radial projection → dedup. Build `nbr`, `nbe`, `edges`,
   `faces` with consistent outward orientation.
3. Spherical geometry: triangle areas (excess), circumcenters, signed \(\ell^*\), dual polygon
   areas \(A_y\), diamond areas \(A_e\), sub-areas \(\tilde A_i\), \(\bar a_s\).
4. Tangent frames \(e^a\) at both ends of every edge.
5. **Spin connection** by the coordinate-free route (SO(3) transport → SO(2) angle → Spin(2) lift
   → \(\mathbb F_2\) solve for the \(\mathbb Z_2\) signs → verify every face holonomy), plus the
   polar-chart oracle `omegaChart`.
6. `Lat`: \(\kappa\), \(\kappa'\), volume weights, \(\beta_\triangle\), \(\beta_\ell\), doubler check.
7. `rgeom` app: dump a table of \(L, N_V, N_E, N_F, \bar a_s, \min/\max A_y, \bar a_s/a_t\), and run
   every identity check.

Acceptance: **T1.1a–T1.1h**. Also the flat-equilateral unit test (§2.2 of the formulation) and
gauge (tilt) independence of all scalars.

## WP-B — Zolotarev
Files: `ops/zolotarev.nim`, `tests/tzolo.nim`.

AGM `ellipticK`, descending-Landen `jacobiSn`, `newRat` with poles/zeros/residues in **σ² units**,
partial fractions computed in log space with sign tracking, endpoint-equioscillation
normalization, measured max error over 20001 log-spaced samples, FNV-1a coefficient hash.

Acceptance: **T1.3a, T1.3b**; `ellipticK`/`jacobiSn` against known values; orders 11 and 31 both
built and their error compared.

## WP-C — Analytic continuum formulas
Files: `core/analytic.nim`, `tests/tanalytic.nim`.

(V.3) closed form and truncated sums, (V.14), the flat lattice spectrum (IV.8), the two S²
checks (C.54)/(C.55) and (C.56)/(C.57), image sums for finite antiperiodic T, and the
\(\Delta_{\rm eff}\) definition (V.4)–(V.5).

Acceptance: **T1.6a, T1.6b**; \(G^{(1,1)}(t)\to1/(4\pi t^2)\) as \(t\to0\); truncated sum → closed form.

## WP-D — Spinor fields, CG, multishift CG
Files: `core/spinor.nim`, `ops/solve.nim`, `tests/tspinor.nim`, `tests/tsolve.nim`.

Plain-array vector ops; CG; Jegerlehner multishift CG with per-shift **true-residual
recomputation and single-shift refinement fallback**; zero initial guess always; the 1.001
roundoff guard.

Acceptance: **T1.3h**; CG against a dense LAPACK solve on a random SPD matrix; allocation
regression.

## WP-E — Wilson-Dirac operator
Files: `ops/wilson.nim`, `tests/twilson.nim`.

Exactly (IV.1) with the free-limit \(\kappa,\kappa'\); the volume-normalized \(\hat D_W\) wrappers;
explicit adjoint; tangent `applyDwDeriv`; `dwPullback`; `denseDw`.

Acceptance: **T1.2a–T1.2e**, test-ladder steps 1–5.

## WP-G — Non-compact U(1) gauge action, zero modes, flow
Files: `ops/gaugeact.nim`, `ops/flow.nim`, `tests/tgauge.nim`, `tests/tflow.nim`, app `rgauge.nim`.

Gaussian action over spatial and temporal plaquettes with \(\beta_\triangle,\beta_\ell\); analytic
force validated against a dense incidence oracle and finite differences; divergence/gradient and
the CG zero-mode projector; exact Gaussian heatbath; the double-CG pseudo-inverse propagator
(V.16)–(V.17); gradient flow on `algorithms/rk`.

Acceptance: **T1.5g**; \(S=0\) on pure gauge; \(\sum_f\Theta_f=0\) per time slice; flow vs exact
matrix exponential on a small lattice.

## WP-F — Overlap operator
Files: `ops/overlap.nim`, `tests/toverlap.nim`.

\(X=\hat D_W-M\), \(H=X^\dagger X\), \(D_{\rm ov}=1+X\,R(H)\) with the frozen Zolotarev window;
persistent workspace pool; `kernelWindow` monitor that **stops the run** when the window is
violated; the single shared `ovGradient` pullback; `denseOv` polar-factor oracle.

Acceptance: **T1.3c–T1.3h**, test-ladder steps 1–6.

## WP-K — Free-limit campaign (Tier 1 headline)
Files: apps `rfree.nim`, `rspec.nim`; `meas/fit.nim` from WP-J.

Reproduce, with numbers written to TSV and figures scripted:
Figs. 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 and Table I; \(\Delta_0\) values and \(n_{\max}(L)\).

Acceptance: **T1.2f–T1.2h, T1.3g, T1.4a–T1.4i, T1.5a–T1.5f**.

## WP-H — HMC
Files: `hmc/pseudofermion.nim`, `hmc/trajectory.nim`, `tests/thmc.nim`, app `rhmc.nim`.

Pseudofermions with a Hasenbusch ladder; `mdevolve` nested Omelyan; `hmc/metropolis`
accept/reject; **projection of the committed field, the refreshed momentum, and every MD force**;
trajectory-addressed RNG; versioned checkpoint with full manifest validation.

Acceptance: test-ladder steps 7–12.

## WP-I — Measurements
Files: `meas/harmonics.nim`, `meas/observables.nim`, `meas/gevp.nim`, `tests/tmeas.nim`.

Propagators and effective masses; \(\sigma_{PS},\sigma_{FS}\) and the condensate; conserved link
currents (vector and axial, connected and disconnected); \(\ell,m\) projection on the icosahedral
lattice; the 7 Wilson-loop shapes, \(J_{\rm top}\) and \(F^2\); GEVP through
`eigens/linalgFuncs.zeigsgv`.

Acceptance: cold-field zeros; stochastic condensate vs the dense-spectrum oracle; \(\ell=1,2\)
degeneracy exact under \(I_h\) (**T2.6 protection claim**); \(\sigma_{PS}\) and \(\sigma_{FS}\)
spectra identical (**T2.9 structural claim**).

## WP-J — Statistics and fitting
Files: `meas/fit.nim`, `meas/dataio.nim`, `tests/tfit.nim`, app `ranalyze.nim`.

Effective mass, plateau fit (V.6), \(O(a^2)\) fit (V.7), \(n_{\max}\) least-squares fit (V.9),
jackknife + integrated autocorrelation via `utils/resample`, TSV read/write, gnuplot scripts.

Acceptance: recovers exact answers on synthetic data; AR(1) autocorrelation within tolerance;
\(n_{\max}\) fit reproduces **T1.4f/T1.5e**.

## WP-L — Interacting campaign (Tier 2)
Files: `campaign/` scripts, app driver wiring.

Grid \(g^2a\in\{0.5,1,1.5\}\times L\in\{1,2\}\times N_f\in\{2,4,6\}\) at \(a_t=0.2\)
(L=4 for the pure-gauge-only observables). Produce T2.1–T2.9 to whatever statistics fit, and
**state the statistics honestly in every plot footer**.

---

## Rules for every work package

1. **Work only in the worktree** `/Users/xjin/K/W/P003/qex/.claude/worktrees/qed3-slides-reproduction-plan-0e70b6`.
   Never write to `/Users/xjin/K/W/P003/qex/src/experimental/qed3` or anywhere else in the main
   checkout — it holds the user's uncommitted work.
2. **Never run a state-changing git command.** No `checkout`, `reset`, `clean`, `stash`, `commit`.
3. Build and test from `<worktree>/build_mac` with
   `make run experimental/radial/tests/t<name>`.
4. Follow `CLAUDE.md`: short names, equation-first comments, no defensive wrappers, no accessors
   for plain data, no NaN/Inf special-casing.
5. Every new proc that implements an equation carries the equation and its paper tag in a comment.
6. If a target number does not come out, **report the discrepancy** — do not tune tolerances to
   make a test pass.
