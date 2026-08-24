# QED3 on \(S^2\times\mathbb R\): complete formulation

Everything here is extracted from **arXiv:2510.03085v2** (Boyle, Brower, Fleming, Katz,
Matsumoto, Misra, "Studying QED3 with radial quantization on the lattice: Free limit")
plus the interacting-system slides (see [`01-slides.md`](01-slides.md)).
Equation tags `(IV.1)`, `(C.29)`, … refer to that paper.

This is the **normative reference** for every implementation work package.
If code and this document disagree, one of them is wrong — fix it, don't work around it.

Conventions: sphere radius \(R=1\) throughout, so all dimensionful quantities are in units of \(R\).
\(a_s\) means a spatial (geodesic) edge length, \(\bar a_s\) their average, \(a_t\) the temporal
spacing, \(L_t\) the number of time slices, \(T=a_t L_t\).

---

## 1. Continuum theory

\[
S_{\rm cont}=\int dV\Big[\frac{1}{4g^2}F_{\mu\nu}F^{\mu\nu}
+\sum_{f=1}^{N_f}\bar\psi_f\,\sigma^a e_a^\mu(\nabla^S_\mu+iA_\mu)\,\psi_f\Big]
\equiv S_g+S_{N_f}.
\tag{II.1}
\]
\(N_f\) two-component Dirac fermions, \(N_f\) even.
Coordinates \((x^\mu)=(\theta,\varphi,t)\), metric \(ds^2=dt^2+d\theta^2+\sin^2\theta\,d\varphi^2\) (II.3).
Vierbein (II.4):
\[
e^1_\mu=\partial_\mu\theta,\qquad e^2_\mu=\sin\theta\,\partial_\mu\varphi,\qquad e^3_\mu=\partial_\mu t .
\]
\(\nabla^S_\mu=\partial_\mu+\omega^S_\mu\); \(\sigma^a\) are the Pauli matrices.

### 1.1 Four-component packaging (needed for \(N_f\) and for the scalars)
For \(f=1,\dots,N_f/2\) (II.5):
\(\Psi_f=\binom{\psi_f}{\psi_{f+N_f/2}}\), \(\bar\Psi_f=(\bar\psi_f,\,-\bar\psi_{f+N_f/2})\),
\(\gamma^a={\rm diag}(\sigma^a,-\sigma^a)\) (II.7),
\(\gamma_4=\begin{pmatrix}&1_2\\1_2&\end{pmatrix}\),
\(\gamma_5=-\gamma_1\gamma_2\gamma_3\gamma_4=\begin{pmatrix}&-i1_2\\ i1_2&\end{pmatrix}\) (II.8),
\(\gamma_4=1_2\otimes\tau_1,\ \gamma_5=1_2\otimes\tau_2,\ \gamma_{4,5}\equiv i\gamma_4\gamma_5=1_2\otimes\tau_3\) (II.9).

### 1.2 Discrete symmetries
\(P:(\theta,\varphi,t)\to(\pi-\theta,\varphi+\pi,t)\) (II.10);
\(T:(\theta,\varphi,t)\to(\theta,\varphi,-t)\) (II.11).
Two-component: \(P:\psi\to\sigma_1\psi(x_P),\ \bar\psi\to-\bar\psi(x_P)\sigma_1\) (II.16);
\(T:\psi\to\sigma_3\psi(x_T),\ \bar\psi\to-\bar\psi(x_T)\sigma_3\) (II.17).
Continuum Dirac operator identity (IV.13): \(-\sigma_1 D(x_P;A_P)\sigma_1=D(x;A)\).

---

## 2. The lattice: refined icosahedron × \(\mathbb R\)

### 2.1 Construction
Take the regular icosahedron; for refinement level \(L\), split each of the 20 faces into
\(L^2\) small triangles by dividing the edges into \(L\) equal segments and drawing lines parallel
to the edges; then **project every vertex radially onto the unit sphere**. Result:
\[
N_V=10L^2+2,\qquad N_E=30L^2,\qquad N_F=20L^2 .
\]
(Euler: \(V-E+F=2\) ✓.) 12 sites have 5 neighbours, the remaining \(N_V-12\) have 6.
The isometry group of the lattice is the full icosahedral group \(I_h\), which contains \(P\).

Temporal direction: \(L_t\) identical copies with spacing \(a_t\); \(T=a_t L_t\).

### 2.2 Geometric weights (all intrinsic / spherical)
For a link \(y_1y_2\), \(\gamma_{y_1y_2}\) is the geodesic; \(\ell_{y_1y_2}\) its length;
\(e^\alpha_{y_1y_2}=dy^\alpha/ds\) the unit tangent. \(\bar a_s\equiv\langle\ell\rangle_{\rm links}\).

- **Dual point** of a triangle = its circumcenter projected onto \(S^2\).
  For vertices \(a,b,c\) (unit vectors), \(\hat c=\pm\,\mathrm{normalize}\big((b-a)\times(c-a)\big)\);
  pick the sign that puts it on the same hemisphere as \(a+b+c\).
- \(A_\triangle\) = spherical area of the triangle = spherical excess.
- \(A_y\) = area of the **dual polygon** around site \(y\) (vertices = circumcenters of the
  incident triangles).
- \(A_{y_1y_2}\) = area of the **diamond** \(y_1\,c_1\,y_2\,c_2\), the two spherical triangles
  formed by the link and the two adjacent dual points (Fig. 2).
- \(\ell^*_{i,i+1}\) = **signed** geodesic distance from the triangle's dual point to the edge
  \(y_iy_{i+1}\). With \(\hat n=\mathrm{normalize}(y_i\times y_{i+1})\), \(\ell^*=\arcsin(\hat c\cdot\hat n)\)
  with the sign convention that \(\ell^*>0\) when the dual point is inside the triangle.
  (Obtuse triangles give \(\ell^*<0\); the sum \(\sum_i\ell^*_{i,i+1}\ell_{i,i+1}/2=A_\triangle\) still holds.)
- \(\tilde A_i\) = area of the sub-region of \(\triangle\) belonging to vertex \(i\)
  (bounded by the two half-edges at \(y_i\) and the two dual lines), so
  \(\sum_{i=1}^3\tilde A_i=A_\triangle\) and \(\sum_{\triangle\ni y}\tilde A_{(\triangle,y)}=A_y\).

**Identities that MUST be tested** (they are the correctness criteria for the geometry module):
| identity | why |
|---|---|
| \(\sum_\triangle A_\triangle = 4\pi\) | closure of the sphere |
| \(\sum_y A_y = 4\pi\) | dual tiling closure |
| \(\sum_{\rm links} A_{y_1y_2} = 2\cdot 4\pi/\ldots\) → use \(A_{y_1y_2}=\tfrac12\ell_{y_1y_2}(\ell^{*}_{1}+\ell^{*}_{2})\) | consistency of diamond and dual-length definitions |
| \(\sum_i 4\arctan\!\big[\tan(\ell_i/4)\tan(\ell^*_i/2)\big]=A_\triangle\) | dual decomposition, **exact on the sphere** |
| \(\sum_{i}\tfrac12\ell_{i,i+1}\ell^*_{i,i+1}=A_\triangle+O(\bar a_s^2)\) | the *flat* form of the same identity — only \(O(a^2)\), **do not assert it exactly** |
| \(\sum_{\triangle\ni y}\tilde A_{(\triangle,y)}=A_y\) | dual-cell decomposition |
| \(\sum_{i}\ell_{i,i+1}\ell^*_{i,i+1}e^a_{i,i+1}e^b_{i,i+1}=A_\triangle\delta^{ab}(1+O(\bar a_s^2))\) | **(IV.6)**, the simplicial (Christ–Friedberg–Lee) closure relation — the deepest test |
| \(\sum_{\rm links}\ell=N_E\bar a_s\) | definition of \(\bar a_s\) |

**Flat equilateral cross-check** (exact, use as a unit test with a flat triangular patch):
edge \(a\), \(\ell^*=a/(2\sqrt3)\), \(A_\triangle=\sqrt3a^2/4\), \(A_y=\sqrt3a^2/2\),
\(A_{y_1y_2}=a^2/(2\sqrt3)\), hence \(\kappa=1/\sqrt3\), \(\kappa'=(\sqrt3/2)(a/a_t)\) — matching (IV.8).

### 2.3 Lattice spin connection \(\Omega_{y_1y_2}\)
\[
\Omega_{y_1y_2}=\cos\frac{\omega_{y_1y_2}}{2}+i\sigma_3\sin\frac{\omega_{y_1y_2}}{2}
= \exp\!\big(i\sigma_3\,\omega_{y_1y_2}/2\big),
\tag{III.1}
\]
\[
\omega_{y_1y_2}=\int_{\gamma_{y_1y_2}}\!ds\;e^\alpha_{y_1y_2}\omega^{12}_\alpha
=-\int_{\gamma_{y_1y_2}}\!ds\;\frac{d\varphi(s)}{ds}\cos\theta(s).
\tag{III.2}
\]
Under parity \(\omega_{y_1^Py_2^P}=-\omega_{y_1y_2}\) (III.3), i.e. \(\Omega_{y_1^Py_2^P}=\Omega^\dagger_{y_1y_2}\) (IV.16).

**Closed-form construction (preferred; no quadrature).** From the tetrad hypothesis (A.13)
\(\gamma_ae^a_{y_1y_2}(y_1)\Omega_{y_1y_2}=\Omega_{y_1y_2}\gamma_ae^a_{y_1y_2}(y_2)\) one derives
\[
\boxed{\;\omega_{y_1y_2}=\alpha_2-\alpha_1\;}
\]
where \(\alpha_k\) is the angle of the **outgoing geodesic tangent** at \(y_k\) measured in the local
orthonormal frame \((\hat e_\theta,\hat e_\varphi)\):
\(\hat t=\cos\alpha\,\hat e_\theta+\sin\alpha\,\hat e_\varphi\), i.e.
\(\alpha=\mathrm{atan2}(\hat t\cdot\hat e_\varphi,\ \hat t\cdot\hat e_\theta)\), with
\(\hat t(y_1)\) pointing from \(y_1\) toward \(y_2\) and \(\hat t(y_2)\) pointing away from \(y_1\).
Derivation: \(e^{-i\sigma_3\omega/2}(\cos\alpha_1\sigma_1+\sin\alpha_1\sigma_2)e^{i\sigma_3\omega/2}
=\cos(\alpha_1+\omega)\sigma_1+\sin(\alpha_1+\omega)\sigma_2\).

The \(2\pi\) branch of \(\alpha_2-\alpha_1\) is *not* free — \(\Omega\) flips sign under
\(\omega\to\omega+2\pi\). Resolve it by **continuous tracking**: sample \(\alpha(s)\) at
\(N_q\simeq32\) points along the geodesic, unwrap, and take \(\omega=\alpha(\ell)-\alpha(0)\).
This is exactly \(-\int\cos\theta\,d\varphi\) and is the reference implementation.

**Pole cut / antiperiodicity in \(\varphi\).** For a closed curve \(\gamma\) (B.4):
\[
\oint_\gamma dy^\alpha\omega^S_\alpha=\frac{i\sigma_3}{2}\Big[\,{\rm Area}({\rm int}\,\gamma)-2\pi\eta(\gamma)\Big],
\]
with \(\eta=\pm1\) if \(\gamma\) encircles the north (south) pole in the positive (negative)
\(\varphi\) direction, else 0. Introduce a **cut** \(\gamma_c\): a meridian from N to S pole placed so
it touches no site; every link whose geodesic crosses \(\gamma_c\) gets an extra factor \(-1\) on
\(\Omega\). This implements the antiperiodic boundary condition for \(\psi\) in \(\varphi\)
(\(m\in\mathbb Z+\tfrac12\)) and exactly cancels the \(\eta\) term.

**Practical gauge choice.** Tilt the polar axis so that (i) no site sits near a pole and (ii) no
link geodesic passes near a pole. This is a choice of local Lorentz gauge and is physically
irrelevant; observables must be invariant under it — **test this**.

**Validation test (the strongest one available).** For every triangle, with edges traversed
counterclockwise as seen from outside,
\[
\Omega_{12}\Omega_{23}\Omega_{31}=\exp\!\big(i\sigma_3 A_\triangle/2\big)\times(\pm1)_{\rm cut},
\]
i.e. \(\omega_{12}+\omega_{23}+\omega_{31}\equiv \pm A_\triangle \pmod{2\pi}\)
— this is the "\(S^2\) geometry is protected" statement on slide 7.
Also check that \(\sum_\triangle A_\triangle=4\pi\) makes the global holonomy trivial.

### 2.4 Tangent-frame components
\(e^a_{y_1y_2}\equiv e^\mu_{y_1y_2}e^a_\mu\), i.e. at \(y_1\):
\(e^{1}=\cos\alpha_1,\ e^{2}=\sin\alpha_1\) (and \(e^3=0\) for spatial links).
Only \(e^a_{y_1y_2}(y_1)\) enters (IV.1).

---

## 3. Wilson fermion (two-component)

\[
\begin{aligned}
S_W=\;&\sum_t\sum_{y_1,y_2}\kappa_{y_1y_2}\,\bar\psi_{y_1,t}
\Big[-\tfrac12\big(1-e^a_{y_1y_2}(y_1)\sigma_a\big)U_{y_1,t;y_2,t}\Omega_{y_1y_2}\psi_{y_2,t}
+\tfrac12\psi_{y_1,t}\Big]\\
+&\sum_t\sum_{y}\kappa'_y\,\bar\psi_{y,t}
\Big[-\tfrac12(1-\sigma_3)U_{y,t;y,t+1}\psi_{y,t+1}
-\tfrac12(1+\sigma_3)U_{y,t;y,t-1}\psi_{y,t-1}+\psi_{y,t}\Big]
\end{aligned}
\tag{IV.1}
\]
with **free-limit couplings** (IV.2)
\[
\boxed{\ \kappa_{y_1y_2}=\frac{2A_{y_1y_2}}{\bar a_s\,\ell_{y_1y_2}}
=\frac{\ell^{*}_{1}+\ell^{*}_{2}}{\bar a_s},\qquad
\kappa'_y=\frac{A_y}{\bar a_s\,a_t}\ }
\]
\(\kappa_{y_1y_2}\) is nonzero only for nearest neighbours and is symmetric.
The second form follows from \(A_{y_1y_2}=\tfrac12\ell(\ell^*_1+\ell^*_2)\).

**These same couplings are used unchanged in the interacting theory** (slide 7 working hypothesis).

Gauge links: \(U_m=\exp(i\theta_m)\) with \(\theta_m\in\mathbb R\) (non-compact).
\(U_{y,t;y,t-1}=U^*_{y,t-1;y,t}\).

**Boundary conditions**: antiperiodic in \(t\); antiperiodic in \(\varphi\) via the cut (§2.3).

**Structure.** Rearranged per prism (IV.3): \(S_W=\sum_{(\triangle,t)}\bar\psi(C^{(\triangle,t)}+B^{(\triangle,t)})\psi\)
with the naive (antihermitian) part
\[
\bar\psi C^{(\triangle,t)}\psi=\sum_{i}\frac{\ell^*_{i,i+1}}{2\bar a_s}
\Big[\bar\psi_{i,t}\Omega_{i,i+\frac12}e^a_{i,i+1}\sigma_a\Omega_{i+\frac12,i+1}\psi_{i+1,t}+({\rm h.c.\ direction})\Big]
+\sum_i\frac{\tilde A_i}{4\bar a_s a_t}\Big[\bar\psi_{i,t}\sigma_3(\psi_{i,t+1}-\psi_{i,t-1})-(\bar\psi_{i,t+1}-\bar\psi_{i,t-1})\sigma_3\psi_{i,t}\Big]
\tag{IV.4}
\]
and \(B^{(\triangle,t)}\) the (hermitian) Wilson term. Consequence: **\(C\) antihermitian, \(B\) hermitian**
⇒ imaginary and real parts of the spectrum respectively. Use this as a numerical test.

Continuum limit (IV.7): \(C^{(\triangle,t)}=\frac{1}{\bar a_s a_t}A_\triangle a_t\,\sigma^ae^\mu_a\overleftrightarrow\nabla^S_\mu(1+O(a^2))\)
using the simplicial closure relation (IV.6).

### 3.1 Flat-limit spectrum (validation target)
On the equilateral triangular lattice with radial extension (IV.8):
\[
\lambda_{\rm flat}=\kappa\big(3-\cos k_1-\cos k_2-\cos(k_1{+}k_2)\big)+\kappa'\big(1-\cos k_t\big)
\pm i\sqrt{\kappa^2\big(e_{12}\sin k_1+e_{23}\sin k_2+(e_{12}{+}e_{23})\sin(k_1{+}k_2)\big)^2+\kappa'^2\sin^2k_t}
\]
with \(\kappa=1/\sqrt3\), \(\kappa'=(\sqrt3/2)(\bar a_s/a_t)\), \(e_{12}=(1,0)\), \(e_{23}=(-1/2,\sqrt3/2)\).
Reproducing Fig. 5 of the paper is a cheap standalone check.

### 3.2 Volume factor / generalized eigenvalues
The lattice operator carries a volume factor (IV.11):
\((\bar a_sa_t)D_{\rm lat,xx'}\Leftrightarrow\sqrt{g}\,d^3x\,D(x)\delta^3(x-x')d^3x'\).
Define \(\delta V\equiv{\rm diag}(A_ya_t)\) and \(\overline{\delta V}\) its average. Then the
**generalized eigenvalues** that approach the continuum spectrum are
\[
D_{\rm lat}\psi=\tilde\lambda\,(\overline{\delta V})^{-1}\delta V\,\psi
\quad\Longleftrightarrow\quad
\tilde\lambda\ \text{= eigenvalues of}\ \ {\rm diag}\!\big(\overline{A_ya_t}/(A_ya_t)\big)\,D_{\rm lat}.
\]

> **Corrected 2026-08-21 (WP-E).** This section originally carried the weight the other way up,
> \({\rm diag}(A_y a_t/\overline{A_ya_t})\), transcribing (IV.12) as printed. That is inconsistent
> with (IV.11): per prism, \(C^{(\triangle,t)}\sim\frac{1}{\bar a_sa_t}A_\triangle a_t\,\sigma\!\cdot\!\nabla\),
> so summing the prisms around a site gives \((D_{\rm lat}\psi)_y\simeq(A_y/\bar a_s)(D\psi)(y)\),
> i.e. \(D_{\rm cont}=\bar a_s\,{\rm diag}(1/A_y)\,D_{\rm lat}\) — the compensating weight is
> \(\overline{\delta V}/\delta V\). Verified analytically on the flat equilateral lattice
> (\([C\psi]=(A_y/\bar a)\sigma\!\cdot\!\partial\psi\) exactly) and numerically: with this weight the
> free spectrum converges to \(\pm i(\ell+1)\) (max dev 0.084 → 0.024 → 0.006 for L = 1, 2, 4)
> while the inverted weight moves away (0.134 → 0.508 → 0.745). Note the correct
> \(\tilde\lambda\) satisfies \({\rm eig}(\hat D_W)=\bar a_s\,{\rm eig}(D_{\rm cont})\) up to the
> Wilson term.

---

## 4. Overlap fermion

\[
X\equiv D_W-M,\qquad
D_{\rm ov}\equiv 1+X\,\frac{1}{\sqrt{X^\dagger X}} .
\tag{IV.9}
\]
Note: in 3D the two-component \(D_W\) has **no \(\gamma_5\)-hermiticity**, so this is the *unitary
polar factor* of \(X\), not a hermitian sign function. \(M=1\) throughout the papers.

**Allowed range of \(M\)** (IV.10):
\[
0<M<\alpha M_0,\qquad M_0=\min\Big(\frac{4}{\sqrt3},\ \frac{\sqrt3\,\bar a_s}{a_t}\Big),\qquad \alpha=0.9 .
\]
Under \(\bar a_s/a_t\ge4/3\) this reduces to the conventional \(0<M<2\).

**Ginsparg–Wilson** (IV.22): \(\Gamma\mathcal D_{\rm ov}+\mathcal D_{\rm ov}\Gamma=\mathcal D_{\rm ov}\Gamma\mathcal D_{\rm ov}\)
for \(\Gamma=\gamma_4,\gamma_5\), with \(\mathcal D_{\rm lat}={\rm diag}(D_{\rm lat},D^\dagger_{\rm lat})\) (IV.21).
Also \(D_{\rm ov}^{-1}=1-(D_{\rm ov}^\dagger)^{-1}\) (IV.17), and the parity relation (IV.18)
\(-\sigma_1(D_{\rm ov}|_{x\to x_P})^{-1}\sigma_1=D_{\rm ov}^{-1}-1\).
**All three are exact numerical identities to test.**

### 4.1 Zolotarev approximation
\((X^\dagger X)^{-1/2}\) is approximated by the optimal rational function of \(1/\sqrt{z}\) on
\([\varepsilon,1]\) after rescaling, with \(\varepsilon=\lambda_{\min}/\lambda_{\max}\) of \(X^\dagger X\):
\[
\frac{1}{\sqrt z}\;\approx\;A\prod_{i=1}^{m}\frac{z+c_{2i}}{z+c_{2i-1}}
\;=\;\sum_{i=1}^{m}\frac{b_i}{z+c_{2i-1}}\ (+\,b_0),
\qquad
c_i=\frac{{\rm sn}^2\!\big(iK'/n;\,k'\big)}{1-{\rm sn}^2\!\big(iK'/n;\,k'\big)} .
\]
Reference: Chiu, Hsieh, Huang, Huang, PRD 66 (2002) 114502 [hep-lat/0206007];
van den Eshof, Frommer, Lippert, Schilling, van der Vorst, CPC 146 (2002) 203.
Requires Jacobi \({\rm sn}\) and the complete elliptic integral \(K(k)\) — **neither exists in QEX**,
so we implement them (AGM for \(K\), descending Landen / AGM for \(\rm sn\)).

Slide 6: **order n = 31 (15 poles)** for the accept/reject action,
**order n = 11 (6 poles)** for the force. Make the order a parameter and reproduce both.

Applying \(D_{\rm ov}\) then needs one **multi-shift CG** solve of \((X^\dagger X+c_{2i-1})\chi_i=\psi\).

### 4.2 Multi-flavor action
Sign-definite form (IV.19)/(IV.20):
\[
S_{N_f,{\rm lat}}=\sum_{f=1}^{N_f/2}\Big[\bar\psi_fD_{\rm lat}\psi_f-\bar\psi_{f+N_f/2}D^\dagger_{\rm lat}\psi_{f+N_f/2}\Big]
=\sum_{f=1}^{N_f/2}\bar\Psi_f\mathcal D_{\rm lat}\Psi_f,\qquad
\mathcal D_{\rm lat}=\begin{pmatrix}D_{\rm lat}&\\&D^\dagger_{\rm lat}\end{pmatrix}.
\]
\(\det\mathcal D_{\rm lat}=|\det D_{\rm lat}|^2\ge0\) ⇒ positive measure for even \(N_f\).
For HMC, one pseudofermion pair per \(N_f=2\); \(N_f=4,6\) = 2, 3 pairs (or Hasenbusch-split).

### 4.3 Scalar operators (slide 6)
With \(\Psi=\binom{\xi}{\eta}\):
\[
\sigma_{PS}\equiv\eta^\dagger\xi+\xi^\dagger\eta \quad\text{(parity symmetric, flavor non-singlet)},
\]
\[
\sigma_{FS}\equiv\eta^\dagger\xi-\xi^\dagger\big(1-D^\dagger_{\rm ov}\big)\eta
\quad\text{(flavor symmetric, parity non-singlet)} .
\]
The soft mass term used for the condensate scan is \(m\,\sigma_{PS}\).

---

## 5. Gauge action (non-compact \(U(1)\), Gaussian)

Degrees of freedom: real link angles
\(\theta_m=\int_{\gamma_m}dx^\mu A_\mu\) (IV.25),
one per **spatial link per time slice** (\(N_E L_t\)) and one per **site per temporal link**
(\(N_V L_t\)). \(U_m=\exp(i\theta_m)\).

\[
S_{g,{\rm lat}}=\sum_{\triangle,t}\frac{\beta_\triangle}{2}\Big(\sum_{m\in p(\triangle,t)}\eta_m\theta_m\Big)^2
+\sum_{\ell,t}\frac{\beta_\ell}{2}\Big(\sum_{m\in p(\ell,t)}\eta_m\theta_m\Big)^2
\tag{IV.24}
\]
with \(\eta_m=\pm1\) the relative orientation, and **free-limit couplings** (IV.26)
\[
\boxed{\ \beta_\triangle=\frac{1}{g^2}\frac{a_t}{A_\triangle},\qquad
\beta_\ell=\frac{1}{g^2}\frac{2A_\ell}{\ell^2_\ell\,a_t}\ }
\]
(\(A_\ell\equiv A_{y_1y_2}\), the diamond area of the link.)

> **Convention finding (WP-G, 2026-08-21).** \(A_\ell\) must be the **exact spherical diamond
> area** \(\sum_{\pm}4\arctan[\tan(\ell/4)\tan(\ell^*_\pm/2)]\), *not* the flat form
> \(\tfrac12\ell(\ell^*_1+\ell^*_2)\): only the exact form tiles the sphere
> (\(\sum_eA_e=4\pi\) to 1e-12; the flat form misses by 3.6/1.0/0.25 % at L=1/2/4), and only it
> reproduces the published \(\Delta_0=1.33242\) (exact: 1.332430; flat: 1.356697; chordal:
> 1.659221). Whether \(\kappa\) in (IV.2) uses the same exact area or the flat identity
> \(\kappa=(\ell^*_1+\ell^*_2)/\bar a_s\) is discriminated by T1.4c (published 0.953918) — see
> [`06-status.md`](06-status.md).

- \(p(\triangle,t)\): the spatial triangular plaquette, 3 links on slice \(t\), oriented consistently.
- \(p(\ell,t)\): the temporal plaquette for spatial link \(\ell=y_1y_2\):
  \(\theta_{\ell,t}+\theta^{(t)}_{y_2,t}-\theta_{\ell,t+1}-\theta^{(t)}_{y_1,t}\).

The action is **exactly Gaussian**: \(S_{g,\rm lat}=\tfrac12\theta^TM\theta\). Consequences:
- pure-gauge configurations can be generated by an **exact heatbath** (no HMC needed) once the
  zero modes are handled;
- the free gauge two-point function is a linear-algebra problem, no Monte Carlo (§7.2);
- the gradient flow \(\dot\theta=-M\theta\) is **linear**.

**Gauge zero modes.** \(M\) is singular; its kernel is the gauge orbit
\(\theta_m\to\theta_m+(\Lambda_{y_2}-\Lambda_{y_1})\). For gauge-invariant observables use the
pseudo-inverse \(\tilde M^{-1}\), realized with CG plus projection (V.16)/(V.17):
```
b' = P b = Mtilde^{-1} (M b)     # first CG: project the source
chi =      Mtilde^{-1} b'        # second CG: the actual solve
```
starting CG from \(\chi_0=0\), so the Krylov space never contains the kernel.

---

## 6. Doubler condition and parameter constraints

Flat-space estimate of the first doubler: \(\min(4\kappa,2\kappa')=\min(4/\sqrt3,\ \sqrt3\bar a_s/a_t)\).
The spatial doubler position is fixed by the construction; the temporal one moves with the
anisotropy. In the flat limit they coincide at \(\bar a_s/a_t=4/3\), and as \(\bar a_s/a_t\to\infty\)
the temporal doubler freezes out. **Require \(\bar a_s/a_t\ge4/3\).**

Approximate \(\bar a_s\) (to be computed exactly by the geometry module):
| L | \(N_V\) | \(N_E\) | \(N_F\) | \(\bar a_s\) (geodesic, measured) | \(\bar a_s/a_t\) at \(a_t=0.2\) |
|---|---|---|---|---|---|
| 1 | 12 | 30 | 20 | 1.107149 | 5.536 |
| 2 | 42 | 120 | 80 | 0.590946 | 2.955 |
| 4 | 162 | 480 | 320 | 0.299474 | 1.497 |
| 8 | 642 | 1920 | 1280 | 0.150227 | 0.751 ✗ (needs smaller \(a_t\)) |

\(a_t=0.2\) is the deck's "fixed \(a_t\)" (consistent with the \({\rm Im}\,\lambda\) ranges of slide 8
and with \(\bar a_s/a_t\ge4/3\) up to L = 4). The exact \(\bar a_s\) values must come from the code.

---

## 7. Analytic continuum correlators (the free-limit targets)

### 7.1 Fermion, temporal, coincident spatial point
\[
G(t)=\sigma_3\,{\rm sign}(t)\,\frac{1}{4\pi}\sum_{n\ge0}(n+1)e^{-(n+1)|t|}.
\tag{V.3}
\]
Closed form for the (1,1) component:
\[
G^{(1,1)}(t)=\frac{{\rm sign}(t)}{4\pi}\frac{e^{-|t|}}{(1-e^{-|t|})^2}
=\frac{{\rm sign}(t)}{16\pi\sinh^2(|t|/2)} \;\xrightarrow{t\to0}\;\frac{{\rm sign}(t)}{4\pi t^2},
\]
the correct 3D flat-space normalization. Lowest dimension \(\Delta_0=1\) (integer tower ⇒ CFT with
a primary of \(\Delta=1\)).
Truncated version used for the \(n_{\max}\) fits:
\(G^{(1,1)}(t;n_{\max})=\frac{1}{4\pi}\sum_{n=0}^{n_{\max}}(n+1)e^{-(n+1)|t|}\).

Lattice counterpart (V.2): \(G_{\rm lat}(x,x')=\frac{1}{\bar a_sa_t}[D_{\rm lat}^{-1}]_{xx'}\).

### 7.2 Gauge current, temporal
\[
J_\rho=\tfrac12\epsilon^{\rho\mu\nu}F_{\mu\nu},\qquad J^t=\frac{1}{\sqrt g}F_{\theta\varphi},
\tag{V.10-11}
\]
lattice operator (V.12): \(J^t_{\rm lat}=\frac{1}{A_\triangle}\sum_{m\in p(\triangle,t)}\eta_m\theta_m\).
\[
G_g(t)=\frac{1}{g^2}\langle J^t(0,0,t)J^t(0)\rangle
=\frac{1}{8\pi}\sum_{n\ge1}\sqrt{n(n+1)}\,(2n+1)\,e^{-\sqrt{n(n+1)}\,|t|}.
\tag{V.14}
\]
Lowest dimension \(\Delta_0=\sqrt2\). Non-integer exponents ⇒ pure gauge (\(N_f=0\)) is not conformal.

Eigenvalues behind these: fermion \(\lambda_{|k|,|m|,n}=\sqrt{k^2+(n+|m|+\tfrac12)^2}\) (C.22);
gauge \(\lambda_{|k|,|m|,n}=k^2+(n+|m|)(n+|m|+1)\) (C.37), zero mode \((m,n)=(0,0)\) excluded.
So on \(S^2\times\mathbb R\) the free gauge tower is \(\Delta_\ell=\sqrt{\ell(\ell+1)}\), \(\ell\ge1\).

### 7.3 Effective dimension and fits
\[
f(t)=\cosh^{-1}\!\frac{G^{(1,1)}(t)}{G^{(1,1)}(T/2)},\qquad
\Delta_{\rm eff}(t)=-\frac{1}{a_t}\big[f(t+a_t)-f(t)\big].
\tag{V.4-5}
\]
Plateau fit \(\Delta_{\rm eff}(t)\simeq\Delta_0+c\,e^{-\Delta't}\) (V.6);
continuum extrapolation \(\Delta_0(\bar a_s,a_t)=\Delta_0^{\rm cont}+c_s\bar a_s^2+c_ta_t^2\) (V.7)
(linear terms forbidden by P and T).

### 7.4 Two-dimensional analytic checks (cheap, exact)
Free fermion propagator on \(S^2\) (C.54): \(\langle\psi(\theta,0)\bar\psi(0,0)\rangle=\dfrac{\sigma_1}{4\pi\sin(\theta/2)}\),
which must equal (C.55) \(\sigma_1\sum_{n\ge0}\frac{(-1)^n}{\pi\sqrt2}\xi_{1/2,n}(-z)\) as \(n_{\max}\to\infty\).
\(J^tJ^t\) on \(S^2\) (C.56-57): \(\frac{1}{g^2}\langle J^t(\theta,0)J^t(0,0)\rangle=\frac{1}{2\pi}[\delta(z-1)-\tfrac12]\).

---

## 8. Published free-limit numbers we must reproduce

| quantity | setup | published value |
|---|---|---|
| \(\Delta_0^{\rm cont}\), fermion | fit over L=2,4,8 × \(L_t\)=120,144,168, T=16, fit range \(4\le t<8\) | **0.999998(34)** (exact 1) |
| \(\Delta_0\), fermion, L=1 | \(L_t=168\), T=16 | **0.953918** |
| \(n_{\max}(L)\), fermion | T=12, \(L_t=168\) | **6, 10, 19, 32** for L=1,2,4,8 |
| residual/DOF, fermion | DOF = 168−2 | **0.028, 0.012, 0.0039, 0.038** |
| \(\Delta_0^{\rm cont}\), gauge | L=2,4,8 × \(L_t\)=48,64,96,120, T=16 | **1.41409(18)** (exact \(\sqrt2=1.41421\)) |
| \(\Delta_0\), gauge, L=1 | \(L_t=120\) | **1.33242** |
| \(n_{\max}(L)\), gauge | T=12, \(L_t=120\), \(1/g^2=20\) | **3, 8, 18, 35** for L=1,2,4,8 |
| residual/DOF, gauge | DOF = 120−2 | **0.0031, 0.0023, 0.0031, 0.0037** |
| doubler table (Table I) | T=16, \(\bar a_s/a_t\ge4/3\) | L=1: all \(L_t\); L=2: \(L_t\ge64\); L=4: \(L_t\ge96\); L=8: \(L_t\ge144\) |

Figures to remake: Fig. 4 (\(D_W\) spectrum, L and \(L_t\) scans), Fig. 5 (flat spectrum),
Fig. 6 (\(D_{\rm ov}\) vs \(D_W\), generalized eigenvalues), Fig. 7 (overlap propagator),
Fig. 8 (\(\Delta_{\rm eff}\)), Fig. 9 (\(O(a^2)\) scaling), Fig. 10 (T-symmetry: \(D_W\) vs \(D_{\rm ov}\)),
Fig. 11 (\(J^tJ^t\)), Fig. 12 (gauge \(O(a^2)\) scaling), Figs. 13–14 (\(S^2\) checks).

---

## 9. CFT / free reference values for the interacting-run ratios

| ratio | CFT | free | source |
|---|---|---|---|
| \(\Delta_J\) | 2 = \(D-1\) | 2 | slide 11 |
| \(\Delta_{\ell=2}/\Delta_{\ell=1}\) (current) | **3/2** | 3/2 | descendant/primary |
| \(\Delta_{\ell=3}/\Delta_{\ell=1}\) (current) | **2** | 2 | slide 13 |
| \(\Delta_V/\Delta_A\) | **1** | 1 | slide 12 |
| \(\Delta_F/\Delta_A\) | **1** | \(1/\sqrt2\) | slide 14 (\(\Delta_F^{\rm free}=\sqrt{1\cdot2}\), \(\Delta_A=2\)) |
| \(\Delta_{F,\ell=2}/\Delta_{F,\ell=1}\) | **3/2** | slide draws \(\sqrt{3/2}\approx1.225\) | **see open question below** |
| \(\Delta_{F^2}/\Delta_F\) | **2** (large \(N_f\)) | — | slide 15 |
| \(\Delta_{PS}/\Delta_A\), \(\Delta_{FS}/\Delta_A\) | — | **1** | slide 16 |

> **Open question (resolve numerically, do not guess).** The free Maxwell tower on
> \(S^2\times\mathbb R\) is \(\Delta_\ell=\sqrt{\ell(\ell+1)}\) (§7.2), giving
> \(\Delta_{\ell=2}/\Delta_{\ell=1}=\sqrt6/\sqrt2=\sqrt3\approx1.732\), not the \(\sqrt{3/2}\approx1.225\)
> drawn on slide 14. Our own free-limit computation (WP-B, WP-F) settles this; report both.

## 10. Open items to settle in code, not by assumption

1. \(a_t\) used in the interacting runs (best evidence: 0.2; see §6).
2. \(L_t\) / \(T\) for the interacting runs (evidence: slide 13 x-axis `dt` to 30 with fit window
   8–16 ⇒ \(L_t\approx60\), \(T\approx12\); slide 11 shows \(t\) to 8 ⇒ \(T=16\)).
3. Exact definition of the 7 Wilson-loop shapes in the slide-14 GEVP.
4. The precise Zolotarev pole counts behind "n=31 → 15 poles" and "n=11 → 6 poles".
5. Sign/branch conventions in \(\Omega\) — fixed by the triangle-area test (§2.3).
6. The free reference in §9 row 6.
