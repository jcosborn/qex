# Reproduction targets

Two tiers.

* **Tier 1 — deterministic.** No Monte Carlo. Every number is a linear-algebra result and must
  match the published value. These are the *hard* acceptance criteria.
* **Tier 2 — Monte Carlo.** Requires HMC ensembles. At the scale available here the numbers are
  statistics-limited; the acceptance criterion is that the machinery is correct and the central
  values are consistent within errors.

---

## Tier 1 — free limit (arXiv:2510.03085), fully deterministic

### T1.1 Geometry
| # | target | criterion |
|---|---|---|
| T1.1a | \(N_V,N_E,N_F=10L^2{+}2,\,30L^2,\,20L^2\) for L=1,2,4,8 | exact integers |
| T1.1b | \(\sum_\triangle A_\triangle=\sum_yA_y=4\pi\) | \(<10^{-12}\) |
| T1.1c | \(\sum_i 4\arctan[\tan(\ell_i/4)\tan(\ell^*_i/2)]=A_\triangle\) per triangle (the **exact spherical** kite identity) | \(<10^{-12}\) |
| T1.1c′ | its flat form \(\sum_i\tfrac12\ell_i\ell^*_i=A_\triangle\) | only \(O(\bar a_s^2)\): 3.56e-2 at L=1, falling like \(\bar a_s^2\) — assert the scaling, not the value |
| T1.1d | closure relation (IV.6) \(\sum_i\ell_i\ell^*_ie^ae^b=A_\triangle\delta^{ab}+O(\bar a_s^2)\) | off-diagonal / diagonal \(=O(\bar a_s^2)\), decreasing with L |
| T1.1e | triangle holonomy \(\prod\Omega=\exp(i\sigma_3A_\triangle/2)\) for **every** face | \(<10^{-12}\) |
| T1.1f | global holonomy over all faces = 1 | \(<10^{-11}\) |
| T1.1g | icosahedral symmetry: \(\bar a_s\), \(A_y\) multiplicities match \(I_h\) orbits | exact orbit counts |
| T1.1h | \(\bar a_s(L)\) table | report; \(\bar a_s(1)=2\arcsin\sqrt{(5-\sqrt5)/10}=1.10715\) |

### T1.2 Wilson operator
| # | target | criterion |
|---|---|---|
| T1.2a | \(C\) antihermitian, \(B\) hermitian | \(<10^{-13}\) |
| T1.2b | explicit adjoint: \(\langle x,D_Wy\rangle=\langle D_W^\dagger x,y\rangle\) | \(<10^{-13}\) rel |
| T1.2c | \(U(1)\) gauge covariance incl. the \(\varphi\) seam and the temporal wrap | \(<10^{-20}\) rel |
| T1.2d | matrix-free apply == dense basis assembly | \(<10^{-12}\) |
| T1.2e | antiperiodic temporal mode: half-integer Matsubara | \(<10^{-12}\) |
| T1.2f | **Fig. 4** — \(D_W\) spectrum, T=4, L=1,2,4 at \(L_t=24\); \(L_t=16,24,48\) at L=2 | remake figure |
| T1.2g | **Fig. 5** — flat-limit spectrum (IV.8) with matching \(\bar a_s,a_t\) | remake figure |
| T1.2h | free spatial spectrum → \(\pm i(\ell+1)\) after volume rescaling | error decreasing in L |

### T1.3 Overlap operator
| # | target | criterion |
|---|---|---|
| T1.3a | Zolotarev: pole interlacing \(0<p_i<z_i<p_{i+1}\), positive residues | exact |
| T1.3b | equioscillation: \(\ge n{+}1\) alternating extrema at \(\pm\)maxRelError | \(<10^{-6}\) rel |
| T1.3c | \(D_{\rm ov}\) vs exact dense polar factor of \(X\) | \(<5\times\)maxRelError |
| T1.3d | Ginsparg–Wilson \(D+D^\dagger-D^\dagger D=0\) (2-component circle identity) | \(<10^{-13}\) — the floor for a rational+iterative operator is \(\approx2\,\)maxRelErr + solver residual, so the original \(10^{-17}\) was unreachable (WP-F; achieved 1–7e-14) |
| T1.3e | \(\Gamma\mathcal D+\mathcal D\Gamma=\mathcal D\Gamma\mathcal D\) for \(\Gamma=\gamma_4,\gamma_5\) (IV.22) | \(<10^{-13}\), same floor; reduces exactly to T1.3d + its adjoint (reduction stated in toverlap.nim) |
| T1.3f | parity relation (IV.18) | \(<10^{-13}\) |
| T1.3g | **Fig. 6** — \(D_{\rm ov}\) vs \(D_W\) spectrum and generalized eigenvalues, T=4, L=4, \(L_t=24\), M=1 | remake figure |
| T1.3h | multishift CG == independent single-shift CG per pole | \(<10^{-18}\) rel |

### T1.4 Free fermion propagator — **the headline number**
| # | target | published |
|---|---|---|
| T1.4a | **Fig. 7** \(G^{(1,1)}(t)\), L=1,2,4, T=12, \(L_t=168\) vs \(\frac{1}{16\pi\sinh^2(t/2)}\) | agreement incl. normalization |
| T1.4b | **Fig. 8** \(\Delta_{\rm eff}(t)\), \(L_t=168\), L=1,2,4,8 | remake figure |
| T1.4c | \(\Delta_0\) at L=1, \(L_t=168\), T=16, fit \(4\le t<8\) | **0.953918** |
| T1.4d | **Fig. 9** \(O(a^2)\) scaling in \(\bar a_s^2\) and \(a_t^2\) | remake figure |
| T1.4e | \(\Delta_0^{\rm cont}\) from (V.7) | **0.999998(34)**, exact 1 |
| T1.4f | \(n_{\max}(L)\) at T=12, \(L_t=168\) | **6, 10, 19, 32** (L=1,2,4,8) |
| T1.4g | residual/DOF for T1.4f, DOF=166 | **0.028, 0.012, 0.0039, 0.038** |
| T1.4h | **Fig. 10** T-symmetry: \(D_W\) violates, \(D_{\rm ov}\) preserves | remake figure |
| T1.4i | Table I doubler table, T=16, \(\bar a_s/a_t\ge4/3\) | exact membership |

### T1.5 Free gauge current correlator
| # | target | published |
|---|---|---|
| T1.5a | **Fig. 11** \(G_g(t)\), \(1/g^2=20\), T=12, \(L_t=120\), L=1,2,4,8 vs (V.14) | agreement incl. normalization |
| T1.5b | \(\Delta_0\) at L=1, \(L_t=120\) | **1.33242** |
| T1.5c | **Fig. 12** \(O(a^2)\) scaling | remake figure |
| T1.5d | \(\Delta_0^{\rm cont}\) | **1.41409(18)**, exact \(\sqrt2\) |
| T1.5e | \(n_{\max}(L)\) at T=12, \(L_t=120\) | **3, 8, 18, 35** |
| T1.5f | residual/DOF, DOF=118 | **0.0031, 0.0023, 0.0031, 0.0037** |
| T1.5g | zero-mode projection: \(\tilde M^{-1}\) via double CG (V.16-17) works, result gauge invariant | \(<10^{-10}\) |

### T1.6 Two-dimensional analytic checks
| # | target |
|---|---|
| T1.6a | **Fig. 13** — \(S^2\) fermion propagator (C.55) → \(\sigma_1/(4\pi\sin(\theta/2))\) as \(n_{\max}\to\infty\) |
| T1.6b | **Fig. 14** — \(S^2\) \(J^tJ^t\) delta-function formation (C.57), \(n_{\max}=20,40\) |

---

## Tier 2 — interacting system (the slides)

Parameter grid: \(g^2a\in\{0.5,1.0,1.5\}\), \(L\in\{1,2,4\}\), \(g^2R=(g^2a)L\), \(N_f\in\{2,4,6\}\),
\(M=1\), \(a_t=0.2\).

### T2.1 Slide 8 — Wilson spectrum on configurations *(cheap, high value)*
Eigenvalues of \(D_W\) near the domain-wall height on 3 configurations per point, plus
\(\min|D_W-1|\). Published legend values:

| L | free \(U{=}0\) | \(g^2R\) values and \(\min\|D_W-1\|\) |
|---|---|---|
| 1 | 1.154 | 0.5 → 0.993 · 1.0 → 0.814 · 1.5 → 0.682 |
| 2 | 1.010 | 1.0 → 0.858 · 2.0 → 0.650 · 3.0 → 0.610 |
| 4 | 0.965 | 2.0 → 0.849 · 4.0 → 0.742 · 6.0 → 0.604 |

The **free** column is deterministic — reproduce it exactly. The interacting columns are
3-configuration averages and will differ; the trend (monotone decrease with \(g^2R\), i.e. the
additive mass shift) is the target.

### T2.2 Slide 9 — gradient-flow scale scan
\(E_s(t)\sqrt L\) vs \(r/t\) on log-log, nine curves, reference slope 3/2.
Targets: (a) the free-Maxwell \(t^{-3/2}\) slope at small \(r/t\); (b) collapse of all nine curves in
the "weak coupling; universal" window; (c) \(g^2a\)-dependent splitting at large \(r/t\).

### T2.3 Slide 10 — PS condensate
\(\langle\sigma_{PS}\rangle\) vs \(mR\), \(N_f=2\), at \((L,g^2R)=(1,1.5),(2,3.0),(4,6.0)\).
Published points:

| ensemble | \(mR\) → \(\langle\sigma_{PS}\rangle\) |
|---|---|
| L1 g1.5 | 0.1→0.073, 0.2→0.143, 0.3→0.210, 0.4→0.271 |
| L2 g3.0 | 0.1→0.040, 0.2→0.079, 0.3→0.116, 0.4→0.153 |
| L4 g6.0 | 0.1→0.0195, 0.5→0.097, 1.0→0.190, 1.5→0.277 |

Acceptance: linear through the origin (no SSB), slope within errors.

### T2.4 Slide 11 — pseudovector conserved current
Effective mass of the axial link-current correlator; \(\Delta_{\ell=2}/\Delta_{\ell=1}\) vs \((a/R)^2\).
CFT = 3/2. Published band: ≈1.40–1.44 at L=2, ≈1.28–1.36 at L=1.

### T2.5 Slide 12 — vector conserved current
\(\Delta_V/\Delta_A\) vs \(g^2R\) at L=1, \(N_f=2,4,6\). CFT = 1. Requires the disconnected diagram.

### T2.6 Slide 13 — \(\ell=3\) spherical-symmetry breaking
Per-\(m\) effective masses for \(\ell=3\); \(\Delta_{\ell=3}/\Delta_{\ell=1}\) vs \((a/R)^2\), CFT = 2.
Key quantitative claim: **~3 % spread among the \(m\) components at L=2**, and \(\ell=1,2\) protected
by \(I_h\) (their spread must be ~0). This one is testable even at modest statistics.

### T2.7 Slide 14 — gluonic sector
\(J^\mu_{\rm top}=\epsilon^{\mu\nu\rho}F_{\nu\rho}\) via a GEVP over 7 Wilson-loop shapes.
\(\Delta_F/\Delta_A\) vs \(g^2R\) (CFT 1, free \(1/\sqrt2\)); \(\Delta_{F,\ell=2}/\Delta_{F,\ell=1}\)
(CFT 3/2, slide's free line \(\sqrt{3/2}\) — **see the open question in
[`02-formulation.md`](02-formulation.md) §9; compute the free value ourselves**).

### T2.8 Slide 15 — \(F^2\)
\(\Delta_{F^2}/\Delta_F\) vs \(g^2R\); CFT large-\(N_f\) = 2. Published band 2.0–2.3.

### T2.9 Slide 16 — scalars
\(\Delta_{PS}/\Delta_A\) and \(\Delta_{FS}/\Delta_A\) vs \(g^2R\); free = 1; published 0.88–0.98.
Key structural claim: **\(\sigma_{PS}\) and \(\sigma_{FS}\) give identical spectra** — that is a sharp,
cheap test that does not need good statistics.

---

## Priority order for a laptop-scale campaign

1. **T1.1 → T1.6** in full. These are exact and are the real proof the formulation is right.
2. **T2.1 free column** (deterministic) and **T2.6 \(\ell=1,2\) protection** (cheap).
3. **T2.2** gradient flow — needs only gauge configurations, and for the Gaussian action those can
   be generated by exact heatbath, no fermions. Nine curves are cheap.
4. **T2.3** condensate at L=1, \(N_f=2\) — needs dynamical overlap HMC but only a few hundred trajectories.
5. **T2.9** scalar-spectrum identity, then **T2.4** — L=1, \(N_f=2\).
6. **T2.7/T2.8** gluonic GEVP — expensive smearing + GEVP; attempt at L=1.
7. **T2.5** vector current with disconnected — noisiest; attempt last.
