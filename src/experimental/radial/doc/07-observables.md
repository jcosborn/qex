# Observable definitions

The talk names its observables but does not define the lattice contractions, and the free-limit
paper only covers the fermion propagator and \(J^tJ^t\). **This document is my derivation of the
rest.** Where a choice was necessary I say so and give the test that validates it, so a wrong
choice fails loudly rather than producing a plausible number.

Companion to [`02-formulation.md`](02-formulation.md) (which is purely what the papers state).

---

## 1. Conserved link currents

### 1.1 The kernel
The gauge field enters the fermion action only through \(U_l=e^{i\theta_l}\). Gauge invariance of
\(S_F=\bar\Psi\mathcal D(\theta)\Psi\) under \(\theta_l\to\theta_l+(d\Lambda)_l\) says
\(\sum_l(\partial S/\partial\theta_l)(d\Lambda)_l=0\) for every \(\Lambda\), i.e. \(d^\dagger J=0\).
So the **exactly conserved current is simply**
\[
J_l=\frac{\partial S_F}{\partial\theta_l}=\bar\Psi K_l\Psi,\qquad
K_l\equiv\frac{\partial\mathcal D}{\partial\theta_l},
\]
which reproduces slide 7's \(\sum_{y:\,nn}J_{xy}=0\) and \(J_{xy}=\kappa\bar\Psi K_{xy}\Psi\)
(the \(\kappa\) is inside \(K\)).

**Consequence for the code:** `ovGradient` — the pullback already needed for the HMC force — *is*
the current kernel. There must not be a second, separately derived current routine; if there
were, the Ward test would only validate one of them. (This is the one design rule the prior
attempt got right and wrote down.)

### 1.2 Flavor structure: which current is which
The gauge field couples identically to all flavors, so a current for a flavor generator \(T\) uses
the *same* link kernel with \(T\) inserted: \(J^T_l=\bar\Psi\,T\,K_l\,\Psi\). For two bilinears
\(O_i=\bar\Psi T_i K_i\Psi\),
\[
\langle O_1O_2\rangle=-\,{\rm tr}(T_1T_2)\,{\rm tr}\big[K_1SK_2S\big]
\;+\;{\rm tr}(T_1)\,{\rm tr}(T_2)\,{\rm tr}[K_1S]\,{\rm tr}[K_2S],
\qquad S=\mathcal D^{-1}.
\]

| slide | current | \(T\) | contractions |
|---|---|---|---|
| 12 "vector" | flavor **singlet** — this is the current the photon couples to | \(T=1\) | connected **and** disconnected ("conn − disc" on the slide) |
| 11, 13 "axial"/pseudovector | flavor **non-singlet**, the \(\gamma_{4,5}=1\otimes\tau_3\) relative current between the two 2-component blocks | \({\rm tr}\,T=0\) | **connected only** |

Both are conserved and both have \(\Delta=D-1=2\) in the CFT — hence \(\Delta_V/\Delta_A=1\)
(slide 12). The axial one is clean because it has no disconnected piece; that is exactly why the
deck uses it for the precision ratios on slides 11 and 13 and flags "noise from the disc diagram"
on slide 12.

### 1.3 The operator that is actually correlated
Use the **temporal** component, i.e. the current on the temporal link at site \(y\), time \(t\)
— this is the charge density, the natural radial-quantization operator:
\[
O^{J}_{\ell m}(t)\;=\;\sum_{y} A_y\,Y_{\ell m}(\hat y)\,J^{\,t}_{y,t},
\qquad
C^{J}_{\ell m}(\Delta t)=\big\langle O^{J}_{\ell m}(t+\Delta t)\,O^{J\,*}_{\ell m}(t)\big\rangle .
\]
\(A_y\) is the dual area (the lattice measure). Use **real** spherical harmonics so everything
stays real.

---

## 2. Spherical projection and the icosahedral selection rule

Under \(SO(3)\to I_h\), the spin-\(\ell\) representation decomposes as

| \(\ell\) | dim | \(I\) content | degenerate? |
|---|---|---|---|
| 1 | 3 | \(T_1\) | **yes — irreducible** |
| 2 | 5 | \(H\) | **yes — irreducible** |
| 3 | 7 | \(T_2\oplus G\) | **no — splits 3 + 4** |

This is precisely slide 13's "icosahedral symmetry protects \(\ell=1,2\)". It gives two free tests
of the geometry and the measurement code:

* **exact test**: the \((2\ell+1)\times(2\ell+1)\) correlator matrix for \(\ell=1\) and \(\ell=2\)
  must be proportional to the identity to machine precision on *any* configuration (it is a
  symmetry statement, not a statistical one). If it is not, the lattice or the harmonics are wrong.
* **physics measurement**: for \(\ell=3\) the matrix splits into a 3-fold and a 4-fold block; the
  spread between them is the "~3 % at L=2" number the deck quotes as the strength of the
  sphere-breaking terms.

Report \(\ell=3\) as the two block eigenvalues **and** their multiplicity-weighted average — the
deck's claim is that the *average* scales correctly toward the continuum.

---

## 3. Scalar operators

With \(\Psi=\binom{\xi}{\eta}\) and \(\mathcal L=(\xi^\dagger,\eta^\dagger)\gamma_4\mathcal D_{\rm ov}\binom{\xi}{\eta}\),
\(\gamma_4=\begin{psmallmatrix}0&1\\1&0\end{psmallmatrix}\), \(\mathcal D_{\rm ov}={\rm diag}(D,D^\dagger)\):
\[
\mathcal L=\xi^\dagger D^\dagger\eta+\eta^\dagger D\,\xi
\;\Longrightarrow\;
\langle\xi\eta^\dagger\rangle=D^{-1},\quad\langle\eta\xi^\dagger\rangle=(D^\dagger)^{-1}.
\]
Hence, from slide 6,
\[
\langle\sigma_{PS}\rangle=\langle\eta^\dagger\xi+\xi^\dagger\eta\rangle={\rm tr}\,D^{-1}+{\rm tr}\,D^{-\dagger}=2\,{\rm Re}\,{\rm tr}\,D_{\rm ov}^{-1},
\]
\[
\langle\sigma_{FS}\rangle=\langle\eta^\dagger\xi-\xi^\dagger(1-D^\dagger_{\rm ov})\eta\rangle
={\rm tr}\,D^{-1}-{\rm tr}\big[(1-D^\dagger)D^{-\dagger}\big]
={\rm tr}\,D^{-1}-{\rm tr}\,D^{-\dagger}+{\rm tr}\,\mathbb 1 .
\]
The \((1-D^\dagger)\) factor is the Ginsparg–Wilson contact subtraction; the residual \({\rm tr}\,\mathbb 1\)
is a field-independent constant and drops out of every connected correlator.

At finite standard-overlap mass, with
\(D(m)=(1-m/2)D_{\rm ov}+m\), \(S=D(m)^{-1}\), the FS factor must be evaluated as
\[
(1-D_{\rm ov}^\dagger)S^\dagger
=\frac{(1+m/2)S^\dagger-\mathbb 1}{1-m/2}.
\]
This identity, including both the propagator coefficient and the contact coefficient, is used
by the dense and point-source contractions.

**Soft mass (slide 10).** We adopt the standard lattice-QCD overlap convention above. The mass
derivative is \(\partial_mD=1-D_{\rm ov}/2\), so measure the corresponding GW-improved
condensate
\[
\frac1N{\rm Re\,tr}\left[(1-D_{\rm ov}/2)D(m)^{-1}\right]
=\frac1N\sum_k{\rm Re}\frac{1-\lambda_k/2}
{(1-m/2)\lambda_k+m}
\]
with Gaussian volume noise and cross-check it against the exact dense-spectrum value on a
small lattice. The slides' raw \(m\,\sigma_{PS}\) wording and this improved insertion agree in
the continuum; at finite cutoff the latter is the convention deliberately chosen here.

**Acceptance test for this whole section (slide 16):** \(\sigma_{PS}\) and \(\sigma_{FS}\) must give
**identical spectra**. That is a sharp structural check that needs no statistics — if the two
effective-mass curves do not lie on top of each other, one of the two contractions above is wrong.

---

## 4. Gluonic sector

### 4.1 Operators
Non-compact \(U(1)\): a "Wilson loop" is the real flux \(\Theta_C=\sum_{l\in C}\eta_l\theta_l\).
From (IV.28) \(\Theta_C\simeq A_C\cdot\tfrac12\epsilon^{\mu\nu}F_{\mu\nu}\), so \(\Theta\) is linear
in \(F\) (parity-odd, the \(J_{\rm top}\) channel) and \(\Theta^2\) is the \(F^2\) channel.

Topological current \(J^\mu_{\rm top}=\epsilon^{\mu\nu\rho}F_{\nu\rho}\), temporal component
(V.11)/(V.12): \(J^t_{\rm lat}=\Theta_\triangle/A_\triangle\). Projected:
\[
O^{\rm top}_{\ell m}(t)=\sum_\triangle A_\triangle\,Y_{\ell m}(\hat c_\triangle)\,\frac{\Theta_\triangle(t)}{A_\triangle}
=\sum_\triangle Y_{\ell m}(\hat c_\triangle)\,\Theta_\triangle(t),
\]
\(\hat c_\triangle\) = the triangle's dual point. \(F^2\) (the \(0^{++}\) of slide 15):
\[
O^{F^2}_{\ell m}(t)=\sum_\triangle A_\triangle\,Y_{\ell m}(\hat c_\triangle)\Big(\frac{\Theta_\triangle(t)}{A_\triangle}\Big)^2
\;+\;\sum_{\ell\text{-links}} \big(\text{temporal-plaquette analogue}\big),
\]
**vacuum-subtracted** for \(\ell=0\) (the \(\ell=0\) mode carries \(\langle O\rangle^2\); forgetting
this is a classic way to get a flat, meaningless correlator).

### 4.2 The GEVP basis — our choice of "7 Wilson loop shapes"
The deck says "generalized eigenvalues (from 7 Wilson loop shapes)" and plots \(m_{\rm eff}\) versus
**gradient-flow time** \(t\in[0.2,1.6]\), so the correlators are measured on flowed configurations
and the 7 shapes form the variational basis at each flow time. The shapes themselves are not
stated. **Ours** (record any change here):

| # | shape | links | channel |
|---|---|---|---|
| 1 | elementary spatial triangle | 3 | \(F_{\theta\varphi}\) |
| 2 | rhombus: two triangles sharing an edge | 4 | \(F_{\theta\varphi}\) |
| 3 | vertex star: the link ring around a site (5 or 6) | 5–6 | \(F_{\theta\varphi}\) |
| 4 | \(L{\ge}2\) only: the "quadruple" triangle (4 elementary triangles) | 6 | \(F_{\theta\varphi}\) |
| 5 | temporal plaquette on a spatial link | 4 | \(F_{\cdot t}\) |
| 6 | temporal rectangle, extent 2 in \(t\) | 6 | \(F_{\cdot t}\) |
| 7 | temporal plaquette on a *next*-neighbour path (two-link spatial side) | 6 | \(F_{\cdot t}\) |

At L=1 shape 4 degenerates (the sphere is a bare icosahedron) — drop it and run a 6×6 GEVP there,
and say so in the output.

### 4.3 GEVP
\(C(t)v=\lambda(t,t_0)C(t_0)v\), solved with QEX's committed `eigens/linalgFuncs.zeigsgv`
(Hermitian generalized, `zhegv` underneath, with its automatic diagonal-regularisation retry).
\(\Delta_n(t)=-\frac{1}{a_t}\ln\frac{\lambda_n(t+a_t,t_0)}{\lambda_n(t,t_0)}\).

---

## 5. Reference values and what they test

| ratio | CFT | free | note |
|---|---|---|---|
| \(\Delta_{\ell=2}/\Delta_{\ell=1}\), current | 3/2 | 3/2 | descendant / primary; same in both limits, so a pure discretization test |
| \(\Delta_{\ell=3}/\Delta_{\ell=1}\), current | 2 | 2 | ditto |
| \(\Delta_V/\Delta_A\) | 1 | 1 | tests the disconnected diagram |
| \(\Delta_F/\Delta_A\) | 1 | \(1/\sqrt2\) | \(\Delta_F^{\rm free}=\sqrt{1\cdot2}\), \(\Delta_A=2\) |
| \(\Delta_{F,\ell=2}/\Delta_{F,\ell=1}\) | 3/2 | \(\sqrt3\) — **not** the \(\sqrt{3/2}\) drawn on slide 14; see [`06-status.md`](06-status.md) open question 1 | |
| \(\Delta_{F^2}/\Delta_F\) | 2 (large \(N_f\)) | — | |
| \(\Delta_{PS}/\Delta_A\), \(\Delta_{FS}/\Delta_A\) | — | 1 | free scalar bilinear has \(\Delta=2\) |

The free gauge tower on \(S^2\times\mathbb R\) is \(\Delta_\ell=\sqrt{\ell(\ell+1)}\), \(\ell\ge1\),
with degeneracy \(2\ell+1\) — read off from (C.37) with \(\ell\equiv n+|m|\) and the \((2n+1)\)
weight in (V.14). We compute it ourselves in the free limit rather than trusting the slide.
