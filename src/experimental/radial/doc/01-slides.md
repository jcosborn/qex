# Slide deck: "Studying QED3 in radial quantization: Interacting system on coarse lattice"

Source: <https://indico.global/event/16565/contributions/161926/attachments/74534/144803/qed3_v10.pdf>
Nobuyuki Matsumoto (Boston U, RBRC), Lattice 2026, 18 pages (17 numbered + "Thanks!").
Collaboration: P. Boyle (BNL), R. C. Brower (BU), G. Fleming (FNAL), A. Katz (BU), N. Matsumoto, R. Misra (BU).

Companion paper (the *free limit*, and the source of every formula we implement):
**arXiv:2510.03085v2** [hep-lat], "Studying QED3 with radial quantization on the lattice: Free limit",
Boyle, Brower, Fleming, Katz, Matsumoto, Misra; to appear in Phys. Rev. D; FERMILAB-PUB-25-0712-T.
See [`02-formulation.md`](02-formulation.md) for the complete formulation extracted from it.

This document records *what the slides show*. The reproduction targets extracted from
it are collected in [`03-targets.md`](03-targets.md).

---

## 1. Narrative of the deck

### Slide 1 — Title
Free limit already published (PRD, arXiv:2510.03085). This talk = the **interacting** system
on coarse lattices. Resources: USQCD Type B+C (Fermilab), BU Shared Computing Cluster.

### Slide 2 — "Lattice field theory beyond QCD/SM"
Motivation: lattice QCD succeeded; other QFT branches would benefit. Example: **CFTs in higher
dimensions** — motivated by UV completion of the SM, composite Higgs, and structureless simplicity.

### Slide 3 — "CFT in higher D from lattice"
- A torus IR regulator breaks scale invariance.
- Good IR regulator: the **cylinder** / radial picture \(S^{D-1}\times\mathbb R\).
  Dilatation = translation along \(\mathbb R\). Conformal map introduces curvature.
- **Quantum Finite Element (QFE)** project: QFT on curved lattices
  (Brower–Fleming–Neuberger, Phys. Lett. 2013; Brower, Matsumoto et al., PoS 2024).
- Ambitions: curvature-triggered phase transitions in QCD; strongly-interacting QFT in
  astrophysical/cosmological backgrounds.

### Slide 4 — "QED3 — as a prototype"
Continuum action
\[
S_{\rm cont}=\int d^3x\sqrt{g}\Big[\tfrac{1}{4g^2}F_{\mu\nu}F^{\mu\nu}
+\sum_{f=1}^{N_f}\bar\psi_f\,\sigma^\mu(\partial_\mu+A_\mu+\omega_\mu)\,\psi_f\Big].
\]
- Conformal window for \(N_f\gg1\) (even).  Refs: Deser–Jackiw–Templeton (1982);
  Redlich PRL/PRD (1984).
- Sketch: \(\beta(\mu)\) vs \(1/\mu\) — super-renormalizable at high energy, CFT in the IR.
- Conformality breaking may be tied to SSB of \(SU(N_f)\) flavor symmetry, \(\langle\bar\psi_f\psi_f\rangle\ne0\).
  Parity cannot break in a vector-like theory (Vafa–Witten 1984).
- Also the effective theory of a superconducting system (Dorey–Mavromatos, PLB 1990).

### Slide 5 — "Theory landscape"  (scale hierarchy)
| energy | meaning | dimensionless coupling | regime |
|---|---|---|---|
| \(a^{-1}\) | cutoff | — | lattice junk |
| \(g^2\) | theory scale | \(\hat g_{UV}^2 \equiv g^2 a\) | perturbative for \(g^2\ll a^{-1}\) |
| \(R^{-1}\) | curvature | \(\hat g_{IR}^2 \equiv g^2 R\) | well described by CFT for large \(N_f\) (Appelquist–Nash–Wijewardhana '88, Nash '89) |

**CFT limit: \(g^2R\to\infty\)** (strongly interacting). E.g. Chester–Pufu '16.

### Slide 6 — "Lattice action"
- **Non-compact \(U(1)\)**, Gaussian gauge action.
- **Overlap fermions** (Zolotarev): \(N_f=2,4,6\).
  - HMC, Hasenbusch-preconditioned.
  - **15 poles (n = 31)** for the accept/reject action.
  - **6 poles (n = 11)** for the force.
  - Exact flavor + parity symmetry on the lattice (Karthik–Narayanan PRD 2016a,b).
- Four-component construction \( \mathcal L_{N_f}=(\xi^\dagger,\eta^\dagger)\gamma_4\mathcal D_{\rm ov}\binom{\xi}{\eta}=\sum_f\bar\psi_f D_{\rm ov}\psi_f \).
  Parity: \(\xi\leftrightarrow\eta\). Flavor basis = boundary modes \(\psi_f\).
- Two chirality operators \(\gamma_4,\gamma_5\).
- Two scalar operators:
  - \(\sigma_{PS}\equiv\eta^\dagger\xi+\xi^\dagger\eta\)  (parity symmetric, flavor non-singlet)
  - \(\sigma_{FS}\equiv\eta^\dagger\xi-\xi^\dagger(1-D_{\rm ov}^\dagger)\eta\)  (flavor symmetric, parity non-singlet)

### Slide 7 — "\(S^2\) lattice and anisotropy"
- Refined icosahedron projected onto \(S^2\); refinement levels **L = 1, 2, 4**; \(N_{\rm site}=10L^2+2\).
- Uses the **free-limit couplings** derived in the previous paper — a "handwaving but plausible"
  working hypothesis that they need no non-perturbative retuning. Supporting arguments:
  - Spin connection fixes the curvature per plaquette:
    \(\exp\oint dx^\mu\omega_\mu=\tfrac12\int dx^{\mu\nu}R_{\mu\nu}=\tfrac{i\sigma_3}{2}A_\triangle\)
    ⇒ **\(S^2\) geometry is protected**.
  - Exactly conserved vector currents on links, \( \sum_{y:\,nn}J_{xy}=0,\ J_{xy}=\kappa\bar\Psi K_{xy}\Psi\)
    ⇒ **fermion hopping \(\kappa\) is protected** (does not renormalize).
  - Super-renormalizable ⇒ **gauge coupling \(\beta\) scales trivially**.
- **But** the space-time anisotropy *can* renormalize ⇒ only **ratios** of energy levels
  (= operator dimensions \(\Delta\)) are physical. This is why every plot in the deck is a ratio.

### Slide 8 — "Lattice parameters"
Phase-diagram cartoon in the plane \((1/(g^2R),\ a/R\equiv1/L)\):
- Lines of constant \(g^2a=\hat g^2_{UV}=0.5,\,1.0,\,1.5\) run through the origin.
- Origin = **CFT limit**; horizontal = continuum limit; leftwards = stronger coupling.
- Grey wedge at small \(1/(g^2R)\) = **Aoki phase**; further left "no finite measure".
- Systematics ansatz: \(\langle O\rangle_{\rm CFT}=\langle O\rangle_{\rm lat}+c_{a^2}(a/R)^2+c_{g^2}/(g^2R)+\cdots\)

Bottom: **Wilson spectrum** (complex plane \(\mathrm{Re}\,D_W\) vs \(\mathrm{Im}\,D_W\), by Arnoldi,
around the domain-wall height, fixed \(a_t\)), for L = 1, 2, 4, each with the free case \(U=0\)
and three couplings, plus \(\min|D_W-1|\) in the legend (see [`03-targets.md`](03-targets.md)).
An "additive mass shift" is annotated for L=4.

### Slide 9 — "Scale scan with the gradient flow"
Gradient flow \(\dot A_\mu(t)=-\delta S/\delta A_\mu(x)\); plot of \(E_s(t)\sqrt L\) vs \(r/t\)
on log–log, nine curves (L=1,2,4 × three \(g^2R\)). Reference dashed line of slope 3/2
(free Maxwell in 3D, \(E\propto t^{-3/2}\)).
Three regions: small \(r/t\) (large flow time) = lattice-UV-clean; mid = **weak coupling, universal**
(all nine curves collapse); large \(r/t\) (small flow time) = **lattice UV contaminated,
\(g^2a\)-dependent**.
Cf. A. Hasenfratz, C. Rebbi, O. Witzel 2019; A. Hasenfratz, C. Peterson 2024; Robert, Harlander, Mason 2026.

### Slide 10 — "Conformal window study: Condensate"
Soft mass term with \(\sigma_{PS}\) (breaks flavor symmetry); measure \(\langle\sigma_{PS}\rangle\) vs \(mR\).
**No sign of SSB even at \(N_f=2\)** — the curves extrapolate linearly to zero.
Cross-checked with the overlap spectrum; consistent with Karthik–Narayanan PRD 2016a,b.

### Slide 11 — "Cross-check 1: conserved current (pseudovector)"
CFT formula in flat space
\(\langle J_\mu(x)J_\nu(y)\rangle = I_{\mu\nu}(x-y)/|x-y|^{2(D-1)}\),
\(I_{\mu\nu}(x)=\delta_{\mu\nu}-2x_\mu x_\nu/x^2\), so \(\Delta_J=D-1=2\).
Effective-mass plot for the axial (pseudovector) two-point function, then
\(\Delta_{\ell=2}/\Delta_{\ell=1}\) (first descendant / primary) vs \((a/R)^2\), CFT value **3/2**.

### Slide 12 — "Cross-check 2: conserved current (vector)"
Vector correlator = connected − disconnected; the disconnected piece is noisy.
\(\Delta_V/\Delta_A\) vs \(g^2R\) at L=1 for \(N_f=2,4,6\); CFT (= free) value **1.0**.

### Slide 13 — "More check 1: spherical symmetry breaking, \(\ell=3\)"
Icosahedral symmetry protects \(\ell=1,2\); \(\ell=3\) splits. Six effective-mass panels
(L=1 × g=0.5,1.0,1.5; L=2 × g=1.0,2.0,3.0) resolved by \(m=-3\ldots3\).
- Visible but small (~3 %) systematic spread at L = 2 ⇒ strength of the sphere-breaking terms.
- Anticipated; scaling toward the continuum is what matters; the average scales correctly.
\(\Delta_{\ell=3}/\Delta_{\ell=1}\) vs \((a/R)^2\), CFT value **2**.

### Slide 14 — "More check 2: fermionic vs gluonic sector"
Topological current \(J^\mu_{\rm top}\equiv\epsilon^{\mu\nu\rho}F_{\nu\rho}\) (topologically conserved).
Generalized eigenvalues from a **GEVP over 7 Wilson-loop shapes**.
- \(\Delta_F/\Delta_A\) vs \(g^2R\): CFT = 1, free = \(1/\sqrt2\). Data approach the CFT value at
  strong coupling, **faster for larger \(N_f\)**.
- \(\Delta_{F,\ell=2}/\Delta_{F,\ell=1}\) vs \(g^2R\): CFT = 1.5, free line drawn at \(\sqrt{3/2}\approx1.2247\).
  (See the open question in [`03-targets.md`](03-targets.md): the naive free-Maxwell value on
  \(S^2\times\mathbb R\) is \(\sqrt{6}/\sqrt2=\sqrt3\approx1.732\). We recompute this ourselves.)

### Slide 15 — "Prediction 1: \(F^2\) — the relevant term in the action"
Effective masses of the \(F^2\) (0++) channel at L=1, 2 vs flow time \(t\); then
\(\Delta_{F^2}/\Delta_F\) vs \(g^2R\), with the CFT large-\(N_f\) value **2.0**.

### Slide 16 — "Prediction 2: scalar correlators"
Nontrivial lattice prediction, comparable to the conformal bootstrap.
\(\sigma_{PS}\) and \(\sigma_{FS}\) give **identical spectra** (shown by the overlapping effective masses).
\(\Delta_{FS}/\Delta_A\) and \(\Delta_{PS}/\Delta_A\) vs \(g^2R\), free value **1.0**;
data fall to ~0.88–0.98.

### Slide 17 — Summary / discussion
- Analysis of QED3 on \(S^2\times\mathbb R\) on coarse lattices.
- Conformal-window prediction: \(N_f\ge2\).
- Consistent with analytic CFT predictions for conserved currents.
- Nontrivial predictions deliverable to the CFT community.
- Discussion: more statistics, better GEVP; global \(((a/R)^2, g^2R)\) fit toward the CFT limit;
  finite \(a_t/R\) corrections. Order of limits matters — the Aoki phase can have SSB;
  compactness of \(U(1)\) may be relevant (a finite measure exists with \(g^2a\to\infty\));
  monopole proliferation (Polyakov 1977, cf. Karthik–Narayanan PRD 2019); nontrivial scaling with
  \(g^2a=\)const; zero modes/topological excitations even though \(\pi_1(S^2)=0\)
  (Banks–Casher '80, Marinari–Parisi–Rebbi '81, Leutwyler–Smilga '92).
- Finer lattices: large-scale code with domain-wall fermions (with P. Boyle, in progress).
- If the non-renormalization hypothesis holds: QFT on curved spacetime; chiral fermions on curved
  manifolds (Aoki–Fukaya 2023, Kaplan–Sen 2024, Golterman–Shamir 2024, Clancy–Kaplan 2025,
  Yamamoto–Kan–Fukaya 2025).

---

## 2. What the author actually did (method summary)

1. Built the \(S^2\times\mathbb R\) lattice: refined icosahedron (L = 1, 2, 4) projected onto the
   sphere, \(\times\, L_t\) copies with temporal spacing \(a_t\); geometric weights from spherical
   (geodesic) formulas; lattice spin connection \(\Omega_{y_1y_2}\) from the parallel-transport
   integral, with a pole cut giving the antiperiodic \(\varphi\) boundary condition.
2. Used the **free-limit couplings** \(\kappa,\kappa',\beta_\triangle,\beta_\ell\) of the companion
   paper unchanged in the interacting theory (working hypothesis, argued from the three protection
   mechanisms of slide 7).
3. Gauge sector: **non-compact** \(U(1)\), Gaussian action over spatial (triangular) and temporal
   plaquettes.
4. Fermion sector: 2-component Wilson kernel \(D_W\) ⇒ overlap \(D_{\rm ov}=1+X(X^\dagger X)^{-1/2}\),
   \(X=D_W-M\), \(M=1\); Zolotarev rational approximation, order 31 (15 poles) for the action,
   order 11 (6 poles) for the force; multi-shift CG.
5. Multi-flavor sign-definite action pairing \(D_{\rm lat}\) and \(D^\dagger_{\rm lat}\)
   (four-component \(\mathcal D_{\rm lat}\)) so that \(N_f = 2,4,6\) is a positive measure.
6. **HMC** with Hasenbusch preconditioning; nested integrator.
7. Diagnostics: Wilson (Arnoldi) spectrum around the domain-wall height on gauge configurations;
   gradient-flow scale scan of the smeared energy.
8. Measurements:
   - \(\langle\sigma_{PS}\rangle\) vs soft mass \(mR\) (conformal-window / SSB test),
   - conserved link currents \(J_{xy}\): pseudovector (axial) and vector (conn − disc), projected on
     spherical harmonics \(\ell=1,2,3\) (and \(m\)),
   - topological current \(J_{\rm top}=\epsilon F\) via a **GEVP of 7 Wilson-loop shapes** with
     gradient-flow smearing,
   - \(F^2\) (0++),
   - scalars \(\sigma_{PS},\sigma_{FS}\).
9. Extracted \(\Delta\) from plateaux of effective masses, formed **ratios** (anisotropy-independent),
   plotted vs \((a/R)^2\) (continuum approach) and vs \(g^2R\) (CFT approach).

## 3. Parameter grid used in the deck

Bare coupling \(g^2a=\hat g_{UV}^2\in\{0.5,\,1.0,\,1.5\}\) held fixed across refinements, with
\(a/R\equiv1/L\) so that
\[
g^2R \;=\; (g^2a)\cdot L .
\]

| \(g^2a\) | L=1 | L=2 | L=4 |
|---|---|---|---|
| 0.5 | \(g^2R=0.5\) | 1.0 | 2.0 |
| 1.0 | 1.0 | 2.0 | 4.0 |
| 1.5 | 1.5 | 3.0 | 6.0 |

This reproduces every legend in the deck: slide 8 (L1: 0.5/1/1.5; L2: 1/2/3; L4: 2/4/6),
slide 9 (same nine), slide 10 (L1 g1.5, L2 g3.0, L4 g6.0 — all at \(g^2a=1.5\)),
slides 12/14/15/16 (L1: 0.5/1/1.5; L2: 1/2/3).

Flavors \(N_f\in\{2,4,6\}\). Domain-wall height \(M=1\). \(a_t\) held fixed across L
("fixed \(a_t\)", slide 8); \(a_t=0.2\) is consistent with the doubler condition
\(\bar a_s/a_t\ge4/3\) at all three refinements (see [`02-formulation.md`](02-formulation.md) §6).
