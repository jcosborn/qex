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
| 11, 13 "axial"/pseudovector | the \(\gamma_{4,5}=1\otimes\tau_3\) relative current between the two 2-component blocks of one \(\Psi\) | \(T=\tau_3\) | connected **and** a disconnected piece of its own (below) |

Both are conserved and both have \(\Delta=D-1=2\) in the CFT — hence \(\Delta_V/\Delta_A=1\)
(slide 12).

**The factorized formula above needs a flavor-proportional propagator, and the block
propagator is not one.** In the \((\xi,\eta)\) basis \(\mathcal S={\rm diag}(S,S^\dagger)\) with
\(S\ne S^\dagger\) (no \(\gamma_5\)-hermiticity in 3D), so \({\rm tr}\,\tau_3=0\) does not kill the
one-point function of the block current:
\[
\langle J^{\tau_3}_l\rangle=-{\rm tr}[K_lS]+{\rm tr}[K_l^\dagger S^\dagger]=-2i\,{\rm Im\,tr}[K_lS]
\equiv -i\,\tau_l ,
\qquad
\langle J^{1}_l\rangle=-2\,{\rm Re\,tr}[K_lS]\equiv -T_l .
\]
\(\tau_l\) is parity odd (it averages to zero over the ensemble) but fluctuates configuration by
configuration: it is the lattice form of the mixed Chern–Simons response of the parity-invariant
pair, and its correlator is a genuine piece of the \(\tau_3\) current two-point function. With
\(n_p=N_f/2\) four-component pairs (all with the same propagator) the full correlators, normalized
per pair, are
\[
C_V(t_2,t_1)=-2\,{\rm Re\,tr}[K_2SK_1S]+n_p\big(\langle T_2T_1\rangle-\langle T\rangle^2\big),
\qquad
C_A(t_2,t_1)=-2\,{\rm Re\,tr}[K_2SK_1S]-n_p\big(\langle\tau_2\tau_1\rangle-\langle\tau\rangle^2\big),
\]
where \(\langle\cdot\rangle\) is the gauge average and the fermion-connected trace is common to
both. The weight is \(n_p=N_f/2\), **not** \(N_f\): the connected piece is a sum over pairs, the
disconnected one a double sum. A connected-only axial channel is therefore a *choice* (the deck
may well have made it; it is not derivable from \({\rm tr}\,\tau_3=0\)); rmeas reports both, as
`Delta_A_conn_l*` and `Delta_A_full_l1`, and the vector as `Delta_V_full`. A flavor-adjoint
current between *different* four-component copies has no hairpin, but needs \(N_f\ge4\).

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

**How the blocks are read off (meas/harmonics `icosaProjectors`).** The real \(Y_{3m}\) of our
chart are not adapted to \(I_h\) (the \(z\) axis is a 2-fold axis of the icosahedron), so the
diagonal entries \(C_{mm}\) mix the \(T_2\) and \(G\) eigenvalues and their spread over \(m\) is
frame dependent (rmeas keeps it as `l3_spread_over_m_frame`, the deck's presentation). The
invariant statement uses the projectors \(P_a=\frac{d_a}{60}\sum_{g\in I}\chi_a(g)D^{(3)}(g)\)
built from the 60 rotations: for an \(I\)-invariant matrix \(C=c_{T_2}P_{T_2}+c_GP_G\) and
\(c_a={\rm tr}[P_aC]/d_a\). rmeas measures the full \((2\ell+1)^2\) cross-\(m\) correlator
matrices (`currmat`), reports `Delta_A_conn_l3_T2`, `Delta_A_conn_l3_G`, their weighted mean
`Delta_A_conn_l3_blocks` and `l3_block_split`, and for \(\ell=1,2\) the residual
`l{1,2}_block_residual` of the ensemble mean from a pure block (the "protection", a statement
about the *ensemble average* — a single stochastic-source configuration is not \(I_h\)-symmetric).

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

### 3.1 Connected correlators, and why their equality is not a check

With \(\Pi_t\) the slice projector, \(z_1(dt)={\rm tr}[\Pi_{t_2}S\Pi_{t_1}S]\) and
\(F=(1-D_{\rm ov}^\dagger)S^\dagger\), the fermion-connected pieces are
\[
C^{\rm conn}_{PS}=-2\,{\rm Re}\,z_1,\qquad
C^{\rm conn}_{FS}=-{\rm Re}\,z_1-{\rm Re\,tr}[\Pi_{t_2}F\Pi_{t_1}F].
\]
At \(m=0\), \(F=-S\) by (IV.17), so the two are **identical at every \(dt\ne0\) as an algebraic
identity**, and at \(dt=0\) the contact term cancels because \(S+S^\dagger=1\) gives
\({\rm Re\,tr}[\Pi_tS]=n_V\). This equality is a Ginsparg–Wilson tautology of the connected
contraction (rmeas records it as `psfs_conn_maxdev`, `psfs_contact_m0`); it is **not** a test of
the slide's statement, which is about the full correlators.

### 3.2 The hairpin of the singlet

\(\sigma_{FS}\) is the flavor **singlet** (slide 6: "flavor symmetric"), so its correlator has the
fermion-disconnected piece of §1.2. The one-point functions on a configuration follow from the
propagators of §3 with \(T_t\equiv{\rm tr}[\Pi_tS]\):
\[
\langle\sigma_{PS}(t)\rangle=-2\,{\rm Re}\,T_t,\qquad
\langle\sigma_{FS}(t)\rangle=-T_t+\frac{(1+m/2)\,\overline{T_t}-2n_V}{1-m/2}
\;\xrightarrow{m=0}\;-2n_V-2i\,{\rm Im}\,T_t .
\]
At \(m=0\) the \(\sigma_{PS}\) one-point function is the constant \(-2n_V\) (no hairpin, the GW
statement behind the vanishing condensate), while \({\rm Im}\,T_t=\tfrac12{\rm tr}[\Pi_t(1-V)(1+V)^{-1}]\)
is unconstrained and parity odd: exactly the structure of \(\tau_l\) above. Per pair,
\[
C_{FS}=C^{\rm conn}_{FS}+n_p\Big({\rm Re}\big\langle O_{FS}(t_2)O_{FS}(t_1)\big\rangle-{\rm Re}\,\langle O_{FS}\rangle^2\Big),
\qquad
C_{PS}=C^{\rm conn}_{PS}+n_p\big(\langle O_{PS}(t_2)O_{PS}(t_1)\rangle-\langle O_{PS}\rangle^2\big),
\]
with the plain (unconjugated) product, \({\rm Re}[O_2O_1]=a_2a_1-b_2b_1\) for \(O=a+ib\). At
\(m\ne0\) \({\rm Re}\,T_t\) fluctuates too, so \(\sigma_{PS}\) acquires a hairpin as well (the
non-singlet flavors carry masses \(\pm m\) and their propagators differ).

**Estimators (meas/observables `scalarSample`, `scalarConn`, `scalarOnePoint`).** A noise pair
\((\eta,\xi)\) with \(y=S\xi\), \(z=S\eta\) gives \(a_t=\sum_{x\in t}\eta_x^\dagger y_x\),
\(b_t=\sum_{x\in t}\xi_x^\dagger z_x\), \(d_t=\sum_{x\in t}\eta_x^\dagger z_x\) with
\(E[a_{t_2}b_{t_1}]=z_1\) and \(E[d_t]=T_t\); the cross-sample product of the per-sample
one-point functions is the unbiased hairpin. rmeas writes the ingredients (`scalarvol`,
`scalardisc`) and assembles `Delta_PS_full`, `Delta_FS_full` at analysis time; the point-source
connected correlators remain the precise connected-only reference (`Delta_PS_conn`,
`Delta_FS_conn`). The size of the omitted piece is reported as `fs_hairpin_fraction`.

**Acceptance test for this section:** the stochastic estimators reproduce the dense traces and
`scalarConn` reproduces `scalarCorrDense` exactly (tmeas). Whether the *full* \(\sigma_{PS}\) and
\(\sigma_{FS}\) spectra coincide is a physics question about the interacting theory, answered by
`psfs_full_maxdev` with statistics, not by construction.

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

**As implemented (doc/06 WP-I; `LoopShape` in meas/observables):** shapes 5–7 are realized as the
temporal plaquettes (rectangles) *combined round a face* with the face orientation signs, in the
time-reflection-even second-difference form — the raw per-edge temporal plaquette carries the
arbitrary canonical orientation of its edge and is not \(I_h\)-covariant. Those combinations equal
exact second time differences of the spatial fluxes, so they enlarge the variational space but are
not an independent \(F_{\cdot t}\) channel. At L=1 every spatial shape collapses onto ONE operator
after \(\ell\)-projection, so the basis there is {triangle, temporal plaquette, extent-2 rectangle}
with rank truncation in the GEVP; the full seven need \(L\ge2\).

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
| \(\Delta_V/\Delta_A\) | 1 | 1 | tests the disconnected diagrams (both channels have one, §1.2) |
| \(\Delta_F/\Delta_A\) | 1 | \(1/\sqrt2\) | \(\Delta_F^{\rm free}=\sqrt{1\cdot2}\), \(\Delta_A=2\) |
| \(\Delta_{F,\ell=2}/\Delta_{F,\ell=1}\) | 3/2 | \(\sqrt3\) from (C.37)/(V.14); slide 14 draws its free line at \(\approx1.22\), which would be \(\sqrt6/2=\Delta_2^{\rm free}/\Delta_A\) — a reading, not an established fact | |
| \(\Delta_{F^2}/\Delta_F\) | 2 (large \(N_f\)) | — | |
| \(\Delta_{PS}/\Delta_A\), \(\Delta_{FS}/\Delta_A\) | — | 1 | free scalar bilinear has \(\Delta=2\) |

The free gauge tower on \(S^2\times\mathbb R\) is \(\Delta_\ell=\sqrt{\ell(\ell+1)}\), \(\ell\ge1\),
with degeneracy \(2\ell+1\) — read off from (C.37) with \(\ell\equiv n+|m|\) and the \((2n+1)\)
weight in (V.14), both verified against the published PDF. The lattice free values at finite L
differ from these by O(a²) (the exact L=1 \(\ell=1\) value is the `jtopCorrExact` reference in
tmeas, 6 % below \(\sqrt2\)), so a pure-gauge Monte-Carlo ratio must be compared with the exact
lattice ratio at the same L, not with \(\sqrt3\); pureL1's 1.63(15) is 0.7σ from \(\sqrt3\) and
0.9σ from 3/2 and decides nothing.
