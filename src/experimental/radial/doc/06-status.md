# Status board

Append-only per work package. Each agent updates its own section on completion:
what was built, what passes, what does not, and any interface change (which must also be made in
[`04-interfaces.md`](04-interfaces.md)).

| WP | scope | state | owner |
|---|---|---|---|
| A | geometry, lattice | **done** (T1.1a–h) | WP-A |
| B | zolotarev | **done** (T1.3a, T1.3b) | WP-B |
| C | analytic formulas | **done** (T1.6a, T1.6b) | WP-C |
| D | spinor, CG, multishift | **done** (T1.3h); tests verified by main | WP-D |
| E | Wilson-Dirac | **done** (T1.2a–h) | WP-E |
| G | gauge action, flow | **done** (T1.5g; exact-area convention discovery) | WP-G |
| F | overlap | **done** (T1.3c–T1.3f, T1.3h operator level) | WP-F |
| H | HMC | **done** (section 15 items 7–12; demo run at L=1, n_t=60) | WP-H |
| I | measurements | **done** (T2.6 exact; currents/scalars/gluonic/GEVP machinery) | WP-I |
| J | fits, statistics | **done** (T1.4c–g, T1.5b–f machinery) | WP-J |
| K | free-limit campaign | **done** (T1.2f–g, T1.3g, T1.4a–i, T1.5a–f) | WP-K |
| L | interacting campaign | **machinery done** (rmeas + campaign/t2.sh + smoke chain green; long ensembles NOT started — launch commands below) | WP-L |

## Environment

* Worktree: `/Users/xjin/K/W/P003/qex/.claude/worktrees/qed3-slides-reproduction-plan-0e70b6`
* Build dir: `<worktree>/build_mac` — already configured with the user's
  `qex_conf_mac` settings (clang-mp-22, `-Ofast -march=native -ffast-math`, QMP/QIO, vlen 4).
* **`SDKROOT` must be set** or clang-mp-22 cannot find the system headers:
  ```bash
  cd build_mac && SDKROOT=$(xcrun --show-sdk-path) make run experimental/radial/tests/tgeom
  ```
* Binaries: `build_mac/bin/`
* Verified working: `SDKROOT=$(xcrun --show-sdk-path) make experimental/radial/core/types`

## Known facts established during planning

* `import base` gives params, threading, timers, comms and echo with **no** hypercubic
  Layout/Field dependency. `import qex` would drag all of it in — don't.
* LAPACK is linked automatically by importing `eigens/lapack`. Bound routines: `zgeev`, `zheev`,
  `zhegv`, `zgemm`, `dgetrf`, and the real bidiagonal SVD family. **Not** bound: `zgesv`,
  `zgesvd`, `zggev`, Cholesky, `dsyev`. Add what you need there (5 lines each).
* `eigens/linalgFuncs` gives `zeigs` (Hermitian), `zgeigs` (general eigenvalues only), and
  `zeigsgv` (Hermitian generalized — this is the GEVP we need).
* `algorithms/rk.nim` is fully generic and closure-driven with **zero** QEX coupling — use it for
  the gradient flow. `mdevolve` (nimble) is likewise state-agnostic — use it for MD.
* `hmc/metropolis.nim` is state-agnostic (`MetropolisRootObj` + `update`).
* `utils/resample` provides jackknife and integrated autocorrelation.
* QEX has **no** Zolotarev, Remez, RHMC, overlap, domain-wall, elliptic-function, Legendre or
  spherical-harmonic code. All of that is ours to write.
* QEX's `solvers/cgm.nim` multishift is hard-typed to `Field[V,T]` — not reusable; reimplement.
* Every U(1) path in QEX is **compact**; the non-compact Gaussian action is new code.
* `freezeTimers()` is mandatory in every app — QEX `tic/toc` metadata allocation blew RSS up 900×
  on small volumes in the prior attempt.

## Independent geometry oracle (pure-Python, outside the repo) — **T1.1h and T1.4i already confirmed**

Built while WP-A was running, as a cross-check the Nim code must reproduce. Refined icosahedron,
spherical (geodesic) formulas, dual point = projected circumcenter, \(A_y=\sum_{f\ni y}\tilde A\).

| L | \(N_V\) | \(N_E\) | \(N_F\) | \(\bar a_s\) | \(\sum A_\triangle-4\pi\) | \(\sum A_y-4\pi\) | \(\min A_y\) | \(\max A_y\) | \(\kappa\) range |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 12 | 30 | 20 | **1.107149** | 0 | 5e-15 | 1.047198 | 1.047198 | 0.65911 (all equal) |
| 2 | 42 | 120 | 80 | **0.590946** | 0 | 7e-15 | 0.273844 | 0.309341 | 0.53248–0.66331 |
| 4 | 162 | 480 | 320 | **0.299474** | −3e-14 | −2e-14 | 0.058328 | 0.083977 | 0.45471–0.61252 |
| 8 | 642 | 1920 | 1280 | **0.150227** | −5e-14 | 6e-14 | 0.012966 | 0.022939 | 0.37082–0.62676 |

\(\bar a_s(1)=1.107149\) is the exact icosahedron edge \(2\arcsin\sqrt{(5-\sqrt5)/10}\). ✓

**Closure relation (IV.6)**, normalized by \(\langle A_\triangle\rangle\), showing clean \(O(\bar a_s^2)\):

| L | max off-diagonal | max \|diag − \(A_\triangle\)\| |
|---|---|---|
| 1 | 1.3e-16 (exact, by \(I_h\) symmetry) | 3.5620e-2 |
| 2 | 2.5527e-3 | 1.2720e-2 |
| 4 | 5.0392e-4 | 3.5347e-3 |
| 8 | 1.3058e-4 | 9.1054e-4 |

Successive ratios ≈ 4 for both columns, matching \(\bar a_s^2\) ratios of 4. ✓

> **Corrected 2026-08-21.** This table first read 2.68e-1 / 1.17e-1 / 3.46e-2 / 9.03e-3 in the
> diagonal column. That was **my bug, not the code's**: I projected the vertex tangents onto the
> orthonormal frame at the *face centroid* without renormalizing, so \(|e|^2=1-(t\cdot\hat r)^2<1\)
> and the diagonal was uniformly suppressed by ~24 % at L=1. WP-A caught it from the exact L=1
> value, which is convention-independent: the 3-fold symmetry forces exactly 120° between the
> three unit tangents, so \(M=\tfrac32\ell\ell^*\delta^{ab}\) and
> \(|1.5\,\ell\ell^*-A|/A=3.5620\times10^{-2}\) — which is what the Nim reports and what the fixed
> oracle now reproduces to every digit. Normalizing the projected 2-vector fixes it.

**T1.4i — the paper's Table I is reproduced exactly, all 20 cells.** The condition
\(\bar a_s/a_t\ge4/3\) at \(T=16\) is \(L_t\ge T/(0.75\,\bar a_s)\), giving
\(L_t^{\min}=19.3,\,36.1,\,71.2,\,142.0\) for L=1,2,4,8. Against
\(L_t\in\{64,96,120,144,168\}\): L=1 all, L=2 all, L=4 from 96, L=8 from 144 — which is
character-for-character the paper's table.

Also settles the interacting-run anisotropy: at \(a_t=0.2\), \(\bar a_s/a_t=5.54,\,2.95,\,1.50\)
for L=1,2,4 — all above 4/3, and L=8 (0.75) would not be, which is consistent with the deck
stopping at L=4 with "fixed \(a_t\)".

### Spin connection — and a correction worth reading

The oracle also computes \(\omega\) per oriented edge as the continuously-unwrapped integral of
\(d\alpha\) along the geodesic, with \(\alpha=\mathrm{atan2}(t\cdot\hat e_\varphi,\,t\cdot\hat e_\theta)\),
\(\hat e_\varphi=\mathrm{unit}(\hat z\times p)\), \(\hat e_\theta=\hat e_\varphi\times p\), tilted axis
\(\hat z=\mathrm{unit}(0.13,0.29,0.94)\), 64 samples per edge. Results for L = 1, 2, 4:

* \(\omega_{ab}+\omega_{ba}=0\) to 1e-15 (antisymmetric ✓)
* \(\sum_{\rm face}\omega \equiv +A_\triangle \pmod{2\pi}\) to 2e-15 — **the sign is \(+A_\triangle\)**

**I first read that as "no \(\mathbb Z_2\) sign is needed". That was wrong, and the `mod 2π` is
exactly what hid it.** The argument that settles it needs no computation: on a closed oriented
surface every edge is traversed once in each direction by its two incident faces, so
\(\sum_{\rm faces}\sum_{\rm face}\omega=0\) identically; but \(\sum_\triangle A_\triangle=4\pi\).
Hence the total defect is exactly \(-4\pi\), i.e. **exactly two faces are off by \(2\pi\)** — the two
containing \(\pm\hat z\), where the frame is singular with index \(+1\) each (total index
\(2=\chi(S^2)\)). Since \(\exp(i\sigma_3(A+2\pi)/2)=-\exp(i\sigma_3A/2)\), that is a genuine
\(\mathbb Z_2\) sign in the holonomy, and it is precisely the pole cut of Eq. (B.4).

Verified numerically — for L = 1, 2, 4 alike:
\(\sum_{\rm faces}(\sum\omega-A_\triangle)=-12.566371=-4\pi\), carried by **exactly 2 faces** with
defect \(-2\pi\) each.

What simple connectedness of \(S^2\) actually buys is that the number of defective faces is *even*,
so the class in \(H^2(S^2;\mathbb Z_2)\) is trivial and the \(\mathbb F_2\) solve **succeeds** — with a
nontrivial edge pattern, namely a dual path (discrete meridian) joining the two pole faces.
So **expect a nontrivial `Edge.sgn`**, and note that a global check
\(\prod_{\rm all\ faces}\Omega=1\) cannot detect this either, since \((-1)^2=+1\): only a per-face,
non-modular check can. `tests/tgeom.nim` gets this right (`nbad == 2`).

WP-A must reproduce every number in this section from the Nim code.

## THE COUPLING CONVENTION — settled (main + WP-G, Python oracles, 2026-08-21)

**Normative: both \(\kappa_e\) and \(\beta_\ell\) use the exact spherical diamond (kite) area**
\[
A_e^{\rm exact}=\sum_{\pm}4\arctan\!\big[\tan(\ell_e/4)\tan(\ell^*_\pm/2)\big],\qquad
\kappa_e=\frac{2A_e^{\rm exact}}{\bar a_s\ell_e},\qquad
\beta_\ell=\frac{2A_e^{\rm exact}}{g^2\ell_e^2a_t},
\]
NOT the flat form \(\tfrac12\ell(\ell^*_1+\ell^*_2)\) (i.e. NOT \(\kappa=(\ell^*_1+\ell^*_2)/\bar a_s\),
which (IV.2) writes as an equivalent — it is equivalent only in flat geometry).

Evidence (each an independent six-digit match to a published number):
| number | exact area | flat area | published |
|---|---|---|---|
| fermion \(\Delta_0\), L=1, \(L_t\)=168, T=16 (overlap propagator, main's Python) | **0.953918** | 0.921250 | **0.953918** |
| gauge \(\Delta_0\), L=1, \(L_t\)=120, T=16 (WP-G Nim + main's Python) | **1.332430** | 1.356697 | **1.33242** |

Also structurally: only the exact form tiles the sphere, \(\sum_eA_e=4\pi\) to 1e-12 (flat misses
by 3.6/1.0/0.25 % at L=1/2/4). \(\kappa'\) and \(\beta_\triangle\) already use exact spherical
areas (\(A_y\), \(A_\triangle\)) and are unchanged.

**The slide-8 Wilson-spectrum legends are the exception**: they match the **flat** κ
(L=2: flat 1.0105 → "1.010" vs exact 1.0154; L=4: flat 0.9643 → "0.965"; L=1 additionally needs
\(a_t\approx0.133\): flat 1.1582 → "1.154", exact 1.2104). Evidently older diagnostic plots.
When we remake slide 8 we produce both conventions, exact as primary.

**Integration TODO (main, after WP-H/WP-I land, before WP-K):** switch `core/lattice.nim` (κ) and
the gauge \(\beta_\ell\) path to the exact kite area (add it to `core/geom.nim`; keep the flat
form available for the slide-8 comparison), then re-run the whole test suite and update pinned
numbers that legitimately move (e.g. min\(|{\rm eig}X|\) 1.1234 → 1.1770 at a_t=0.2, Zolotarev
window fixtures). Every updated pin must cite this section.

New oracle pins for the switch (Python, exact kite areas; \(\sum_eA_e^{\rm exact}=4\pi\) to 1e-13
at every L, while the flat sum misses by −3.562e-2/−9.971e-3/−2.550e-3/−6.418e-4 relative):

| L | \(\kappa^{\rm exact}\) range |
|---|---|
| 1 | 0.683450 (all equal) |
| 2 | 0.538115 – 0.669683 |
| 4 | 0.455754 – 0.613770 |
| 8 | 0.370994 – 0.627216 |

## Open questions (resolve in code, record the answer here)

1. Free value of \(\Delta_{F,\ell=2}/\Delta_{F,\ell=1}\): slide 14 draws \(\sqrt{3/2}\approx1.225\);
   the free Maxwell tower \(\Delta_\ell=\sqrt{\ell(\ell+1)}\) gives \(\sqrt3\approx1.732\).
   — **answered by WP-C in the continuum: the free value is \(\sqrt3=1.7320508\), and the slide's
   line is \(\Delta_2^{\rm free}/\Delta_A=\sqrt6/2\), i.e. the wrong denominator.** Lattice
   confirmation still owed by WP-K, but the continuum side is now settled and tested.
   Reasoning: Eq. (C.37) gives \(\lambda=k^2+\ell(\ell+1)\) with \(\ell\equiv n+|m|\ge1\), and the
   \((2n+1)\) weight in (V.14) is the \(2\ell+1\) degeneracy of a single tower, so
   \(\Delta_\ell=\sqrt{\ell(\ell+1)}\): \(\Delta_1=\sqrt2\), \(\Delta_2=\sqrt6\), ratio \(\sqrt3\).
   The *same* \(\Delta_1=\sqrt2\) is what makes the left panel's free line \(1/\sqrt2=\sqrt2/2\)
   correct, so the two panels are mutually inconsistent unless the right one used \(\Delta_A=2\)
   in the denominator — \(\sqrt6/2=\sqrt{3/2}=1.224744871391589\) reproduces the drawn line to
   every digit. The data also favour \(\sqrt3\): they sit at 1.42–1.64, i.e. between CFT 1.5 and
   \(\sqrt3\), approaching 1.5 from above as \(g^2R\) grows, which is the expected
   weak-coupling→free direction.

   **WP-C additions (see the WP-C section below).**
   (a) *Both* \(\Delta_1=\sqrt2\) and \(\Delta_2=\sqrt6\) are now measured, not assumed, out of
   `gaugeG`/`gaugeGPeriodic`: \(-d\ln G_g/dt\) at \(t=21\) gives 1.414213563657 against
   \(\sqrt2=1.414213562373\), and the deviations of \(\Delta_{\rm eff}\) from \(\sqrt2\) at
   \(t=8,12,16\) (7.185e-4, 1.142e-5, 1.817e-7) fall by \(e^{-(\Delta_2-\Delta_1)\Delta t}\) with
   \(\ln(7.185\text{e-}4/1.142\text{e-}5)/4=1.03546\) and \(\ln(1.142\text{e-}5/1.817\text{e-}7)/4=1.03518\)
   against \(\sqrt6-\sqrt2=1.035276\). So the second rung of the tower really is \(\sqrt6\) and the
   ratio really is \(\sqrt3\).
   (b) The right panel's *CFT* line does not discriminate: \(\Delta_2^{\rm CFT}/\Delta_1^{\rm CFT}
   =3/2\) and \(\Delta_2^{\rm CFT}/\Delta_A=3/2\) are the same number (\(\Delta_\ell^{\rm CFT}=\ell+1\),
   \(\Delta_A=2\)). Only the *free* line separates the two readings, which is exactly how a wrong
   denominator would survive a plausibility check.
   (c) The plotted points are nevertheless the labelled quantity. Under the
   \(\Delta_{F,\ell=2}/\Delta_A\) reading the free reference would be 1.2247 and the data at
   1.42–1.64 would sit *above* both it and the CFT 1.5; under the labelled
   \(\Delta_{F,\ell=2}/\Delta_{F,\ell=1}\) reading they interpolate between free \(\sqrt3\) and CFT
   1.5 in the right direction. Conclusion: the **data** are the ratio to \(\Delta_{F,\ell=1}\), only
   the **reference line** was computed with \(\Delta_A=2\). Quote \(\sqrt3\).
2. \(a_t\) and \(L_t\) of the interacting runs. Best evidence \(a_t=0.2\); \(T\in\{12,16\}\).
   — **\(a_t=0.2\) confirmed** by the oracle section above (it is the unique round value giving
   \(\bar a_s/a_t\ge4/3\) at L=1,2,4 and failing at L=8, matching the deck stopping at L=4).
   **\(L_t=60,\ T=12\) now strongly supported by WP-E**: the slide-8 free legend values
   \(\min|D_W-1|\) reproduce *exactly* with the raw \(D_{\rm lat}\) at \(a_t=0.2, L_t=60\) —
   L=2: 1.0105 vs published 1.010 ✓, L=4: 0.9643 vs 0.965 ✓.
   **L=1 sub-question resolved (main, Python oracle + a_t scan).** An independent pure-Python
   Wilson operator (my geometry + tracked-\(\omega\) spin connection with the \(\mathbb Z_2\)
   dual-path fix) reproduces WP-E's raw value **exactly: 1.1234** at \(a_t=0.2\) — two independent
   implementations agree, so the operator is right and the slide differs. Scanning \(a_t\):
   \(\min|D_W-1|\) is monotone decreasing in \(a_t\) (1.2047 at 0.10 → 1.0988 at 0.50), and the
   published 1.154 sits at \(a_t\approx0.138\); the natural candidate \(a_t=16/120=0.1333\) — the
   **free-limit paper's** temporal spacing — gives 1.1582 (0.4 % high). Conclusion: the slide-8
   L=1 panel was evidently made at the free-limit-paper spacing (\(a_t\approx0.133\)), not the
   campaign's \(a_t=0.2\); the L=2/L=4 panels match \(a_t=0.2\) exactly. Remake slide 8 at
   \(a_t=0.2\) for all three L and note the difference.
   Bonus check: with the **flat** couplings (\(\kappa=1/\sqrt3\)) the L=1 value is 0.954 — far
   from everything, so the paper's fermion couplings are definitely the spherical ones.
   And without the \(\mathbb Z_2\) spin-structure fix the value is 0.5493 — wildly wrong; a good
   demonstration that the sign pattern is load-bearing, not a convention.
6. **Units of the domain-wall height M.** WP-E: M lives in **raw \(D_{\rm lat}\) units, not hat
   units** — the slide-8 legend only reproduces with \(X = D_{\rm lat} - M\) built from the raw
   operator, consistent with (IV.10) \(M_0=\min(4/\sqrt3,\sqrt3\bar a_s/a_t)\) being the raw flat
   doubler position. The hat rescaling \( \bar a^2/A_y \) is NOT a scalar (range [1.068, 1.538] at
   L=4). **Settled in doc/04 §10: X is built from the RAW operator with plain adjoints, M=1 —
   the paper's own convention (its whole construction is plain-matrix; volume weights live inside
   \(D_{\rm lat}\)). The hat kernel is a congruence, not a similarity, and is diagnostics-only.**
   — **resolved**
5. **Normalization of the "residual per DOF" in Eq. (V.9).** The paper quotes 0.028/0.012/0.0039/0.038
   (fermion) and 0.0031/0.0023/0.0031/0.0037 (gauge) with DOF = \(L_t-2\), for a *deterministic*
   calculation with no errors — so "residual" is some unstated norm. Relative and log weightings
   both imply ~17 % rms for the fermion L=1 case, which is large next to its 5 % deviation in
   \(\Delta_0\); an absolute norm would be dominated by small \(t\). The gauge column is also
   suspiciously flat in L. WP-K must report \(n_{\max}\) (a robust integer, and the real physics)
   plus the residual under each weighting, and state which one reproduces the published column
   — rather than picking a weighting to make the numbers agree.
   — **partially resolved (main, Python oracle).** At the paper's gauge setup (L=1, \(L_t=120\),
   \(T=12\), 120 points, DOF=118, periodic images, C free): **\(n_{\max}=3\) is the clean minimum
   under absolute, relative, and log weightings alike** — the published integer reproduces, and
   it is sharp (residual rises ×50–70 at \(n_{\max}=2\) or 4). The residual *scale* matches no
   single convention exactly: abs²/DOF = 2.7e-5, rms-abs = 5.2e-3, rel²/DOF = 5.3e-4, rms-rel =
   2.3e-2, vs published 0.0031. The fitted \(C_{\rm opt}\approx1.30\) at \(n_{\max}=3\)
   (compensating the truncated spectral weight). Verdict: quote \(n_{\max}\) (robust) and our own
   residual with the convention stated; do not chase the published 0.0031. — **good enough**
3. The 7 Wilson-loop shapes of the slide-14 GEVP. — **open**
4. Zolotarev pole counts behind "n=31 → 15 poles", "n=11 → 6 poles"; with
   `npole=(n-1)/2` these are 15 and 5, so the slide's "6" probably counts the constant term.
   — **resolved, see WP-B below: 15 and 5 poles; the slide's "6" is 5 poles + the constant.**
7. **Which \(A_{y_1y_2}\) does the paper mean — the flat form \(\tfrac12\ell(\ell^*_1+\ell^*_2)\)
   or the exact spherical diamond area?** They differ by 3.6 % at L=1, 1.0 % at L=2, 0.25 % at
   L=4, and only the exact one satisfies \(\sum_e A_e = 4\pi\).
   — **For the gauge \(\beta_\ell\) of (IV.26) it is the EXACT spherical area: WP-G gets
   \(\Delta_0(L{=}1,L_t{=}120) = 1.332430\) against the published 1.33242 (1e-5 relative) with
   it, and 1.356697 with the flat form.** See the WP-G entry, T1.5b.
   — **For the fermion \(\kappa\) of (IV.2) it is still OPEN**, and it is not obviously the same
   answer: (IV.2) writes \(\kappa=(\ell^*_1+\ell^*_2)/\bar a_s\) as an equivalent form, which is
   the *flat* identity. `core/lattice.newLat` currently uses `Edge.area`, i.e. the flat form.
   **Decide it by rerunning T1.4c (published \(\Delta_0 = 0.953918\), L=1, \(L_t=168\)) both
   ways** — one line in `newLat`, and it fixes the convention for every single-lattice number
   in the project.

---

## WP-B — Zolotarev rational approximation  (done)

Files: [`ops/zolotarev.nim`](../ops/zolotarev.nim), [`tests/tzolo.nim`](../tests/tzolo.nim).
Run: `cd build_mac && SDKROOT=$(xcrun --show-sdk-path) make run experimental/radial/tests/tzolo`.
No interface change: the signatures are exactly those of [`04-interfaces.md`](04-interfaces.md) §8.

### Implemented

* `ellipticK(k)` = π/(2 agm(1,k′)), k′=√(1−k²). AGM loop capped at 60 (K(1)=∞ comes out ~1e18).
* `jacobiSn(u,k)` by the descending Landen/AGM recursion (A&S 16.4). A private `sncn` returns
  **both** sn and cn, because the poles need cn² where sn→1 and `sqrt(1-sn*sn)` throws away
  every digit there; cn = cos φ₀ keeps them. `|k|>=1` is the exact tanh/sech limit (the AGM
  stalls at k=1: c_n/a_n stays 1).
* `newRat(smin, smax, order, nsample)` / `ratValue(r, x)`. order odd ≥ 3, `npole = (order-1) div 2`.
  k = smin/smax, K′ = K(k′), poles/zeros k²sn²(iK′/n;k′)/cn²(…) at odd/even i,
  d0 = 2/(k P(k²)+P(1)) from endpoint equioscillation, residues by exact partial fractions
  accumulated in log space with the sign counted separately, then **everything rescaled into
  σ² units** (pole,zero ×smax², res ×smax, cst /smax). Error measured, not assumed, at `nsample`
  log-spaced points. FNV-1a hash over (order, smin, smax, cst, poles, residues).

### Measured `maxRelErr` = max |1 − √x R(x)|  (default nsample = 20001)

| smin/smax | order 11 (5 poles) | order 31 (15 poles) |
|---|---|---|
| 1/50   | 1.4207e-04 | 1.1586e-12 |
| 1/100  | 4.6486e-04 | 3.2610e-11 |
| 1/200  | 1.1894e-03 | 4.6040e-10 |
| 1/500  | 3.1658e-03 | 7.2656e-09 |
| 1/1000 | 5.7503e-03 | 3.9062e-08 |

`maxAbsErr` = max |R(x) − 1/√x| agrees with `maxRelErr` to 3 digits at smin = 1 (it differs by
the 1/√x weight, so it scales as 1/smax under σ → sσ).

Reference fingerprints, smin = 1 (this build; `-Ofast -march=native -ffast-math`, so treat a
mismatch on other hardware as "recheck", not "bug"):
`(31,1,50) 0x2dc0b3fe35860a96`, `(11,1,50) 0x433ef43afb834497`,
`(31,1,200) 0x69f41622860c40d8`, `(11,1,200) 0xc6a8360754def04f`.

### What the tests pin

`ellipticK` and `jacobiSn` against closed forms (sin, tanh, sn(K)=1, sn(2K)=0, sn(3K)=−1, up to
k = 0.9999999); rejection of smin ≤ 0, smax ≤ smin, even order, order < 3, nsample < 2;
`npole == (order-1) div 2`; **T1.3a** interlacing 0 < p₀ < z₀ < p₁ < … with all residues > 0;
partial-fraction == product form to 1e-12; endpoint equioscillation to 1e-12; **T1.3b** the
equioscillation scan; monotone error in order; the σ² rescaling; hash determinism and sensitivity.
All 15 odd orders 3…31 at both 1/50 and 1/200 are covered by every one of those.

**T1.3b in detail.** 40001 log-spaced samples, one signed peak per constant-sign run of e(x),
each interior peak refined by the parabola vertex through its three neighbouring samples. The
count comes out **exactly order+1** with strictly alternating signs for every odd order 3…31 at
both windows — that is the Chebyshev alternation theorem for the type-(m,m) best approximation of
x^(−1/2), i.e. the proof that this is *the* optimal rational and not merely *a* rational.
Peak heights match ±maxRelErr to 2e-7 relative for every order whose error is above the
double-precision floor (≤ 15 at 1/50, ≤ 21 at 1/200); above that the test falls back to
`2e-7*maxRelErr + 10*floor`, with `floor` **measured at run time** rather than hard-coded — it is
the reported `maxRelErr` of order 41 on the 1/50 window, where the true Zolotarev error is ~1e-19,
so the number (9.99e-15) is pure roundoff.

### Things that bit, and limits

* **Evaluation noise floor ≈ 1e-14 absolute in e(x).** Recomputing the residues by an interleaved
  product (magnitude kept near 1, no logs) differs from the log-space values by 3e-15 (order 11)
  to 4e-14 (order 31): the log form costs ~1e-14 relative on res_j, because each of the ~2m logs
  carries its own 1e-16 relative error and they add. That is well below the 1.2e-12 Zolotarev
  error at order 31 / 1/50, so it does not matter for the physics, but it is what caps the strict
  equioscillation check at order ~21. An interleaved product would buy ~10×; not done, the
  log form is what the interface calls for and the gain changes no acceptance criterion.
* **`maxRelErr` is a grid-sampled max, not the true Δ.** The sampled peak undershoots by
  ~(πh/w)²/2 with w the lobe width: 3e-6 relative at order 31 with the default nsample = 20001.
  Pass nsample ≈ 4e5 if you need `maxRelErr` itself to better than 1e-8 (the tzolo equioscillation
  test does). Corollary: `maxRelErr` is **not** reproducible under σ → sσ to better than ~1e-5 at
  high order, because the two grids do not land on the peaks at the same offsets — the first
  version of the rescaling test asserted 1e-12 on it and was wrong. The *parameters* are exact:
  for a power-of-two s the whole construction is bit-identical, and the test checks
  pole ×s², zero ×s², res ×s, cst /s and e(s²x) = e(x) to 1e-15.
* Nim 2.2's `std/math` no longer exports `ldexp`; 2^N a_N u is formed with an integer shift.
* Nothing in QEX or in this tree could serve as an independent oracle for the poles — the
  equioscillation scan *is* the oracle.

### Open question 4 — answered

`npole = (order-1) div 2` gives **15 poles at n = 31 and 5 poles at n = 11**, and the tests assert
both. The partial-fraction form the solver uses is R(x) = cst + Σ_{j=1..npole} res_j/(x+pole_j), so
n = 11 is **6 terms**: 5 shifted systems for the multishift CG plus one constant that costs no
solve. The slide's "15 poles (n=31)" counts poles and its "6 poles (n=11)" counts partial-fraction
terms; the two are inconsistent with each other, and the constant term is the discrepancy
(n = 31 would be "16" on the same counting). **Shift counts for WP-D/WP-F: 15 for the action
operator, 5 for the force operator.**

---

## WP-J — fits, statistics, data IO  (done)

Files: [`meas/fit.nim`](../meas/fit.nim), [`meas/dataio.nim`](../meas/dataio.nim),
[`tests/tfit.nim`](../tests/tfit.nim).
Run: `cd build_mac && SDKROOT=$(xcrun --show-sdk-path) make run experimental/radial/tests/tfit`.
19 tests, all pass. Dependencies: `std/math`, `std/[os,strformat,strutils,tables]`,
`utils/resample`. Nothing from `core/`, `ops/` or `base` — these are pure value-semantics
modules, so they compile and test with no other work package in place.

### Implemented

* `effMass(c, at, T)` — (V.4)-(V.5) exactly as the paper writes them:
  `f(t) = arccosh(c(t)/c(T/2))`, `Delta_eff(t) = -[f(t+at)-f(t)]/at`, reference index
  `round(T/(2 at))`. Result length `c.len-1`. The ratio is formed by division so
  `f(T/2)` is identically 0 (a reciprocal-multiply can land below 1 and give NaN).
* `plateauFit(m, t0, t1, at, e, maxit)` — (V.6) `Delta_0 + c exp(-Delta' t)` over the index
  window `[t0,t1)`, `t = at*i`. Variable-projection scan over `Delta'` (exact 2×2 linear
  solve for `(Delta_0, c)` at each `Delta'`, geometric grid over `[0.01/span, 10/at]`,
  400 points) then Levenberg-Marquardt on all three. Returns all three parameters, their
  errors, the full 3×3 covariance, `chi2`, `dof = npoint-3`, `iters`, `converged`.
* `contFit(x, y, e)` — inverse-variance weighted line in the centered (Numerical Recipes)
  form; intercept, slope, both errors, `cov(a,b)`, `chi2`, `dof = npoint-2`.
* `contFit2(as2, at2, y, e)` — (V.7) three-parameter weighted plane. **`as2`/`at2` are
  already squared.** Columns centered on their weighted means before the normal equations,
  result mapped back with the exact linear change of variables (covariance included).
* `nmaxFit(g, model, nmaxRange, e)` — (V.9). For each integer `nmax`, `C` is the exact
  weighted linear solution `sum w g m / sum w m^2`; the minimizing `nmax` wins.
  `dof = npoint-2`.
* `jack(s, stride, bs)` → `SeriesStat` — block jackknife built on `utils/resample`
  (`jackknife`, `intAutocorr`, `intAutocorrPositive`); nothing reimplemented.
  `blockSize = ceil(max(1, 2 tau_int))` unless `bs > 0`. Carries mean, err, bias,
  `tau2` (Wolff window), `tau2p` (positive sequence), `neff = n/max(1,tau2)`, `n`, `stride`,
  `blockSize`, `nblock`, `provisional = nblock < 4`.
* `writeTsv` / `readTsv` — `#key=value` header lines, `%.17g` numbers, tab separated.
  `columns` is the one reserved key and comes back as `names`. `#` lines without `=` are
  comments; blank lines are skipped (so gnuplot dataset separators survive a read).

### Measured accuracy on synthetic data (echoed by the test as `measured:` lines)

| fitter | test | achieved | asked |
|---|---|---|---|
| `effMass` | exact `cosh(D(T/2-t))`, 4 (D,T,Lt) sets | worst \|ΔD\| **1.8e-14** | 1e-12 |
| `effMass` | vs the closed form of (V.5) on `exp(-Dt)` | worst **1.2e-14** | — |
| `plateauFit` | noiseless, paper window 4≤t<8, at=16/168 | \|Δd0\| **0**, \|Δc/c\| **2.6e-13**, \|ΔΔ'\| **5.8e-14**, 11 iters | 1e-9 |
| `plateauFit` | noiseless, t∈[0,4.1), at=0.1 | \|Δd0\| **0**, \|Δc/c\| **0**, \|ΔΔ'\| **1.1e-16**, 5 iters | 1e-9 |
| `plateauFit` | noiseless, t∈[1,8), at=0.2 | \|Δd0\| **0**, \|Δc/c\| **1.6e-15**, \|ΔΔ'\| **8.9e-16**, 25 iters | 1e-9 |
| `plateauFit` | 400 Gaussian-noise draws, σ=2e-3 | pull mean **0.011**, rms **0.985**, 0 non-converged | pull rms ≈ 1 |
| `contFit` | exact line, 6 points | \|Δa\| **0**, \|Δb\| **2.8e-16**, chi2 **0** | 1e-14, chi2≈0 |
| `contFit` | e → 4e | ea, eb ×4 and cab ×16 to **<1e-12** | exact scaling |
| `contFit2` | exact plane, 4×5 grid | \|Δa\| **5.6e-16**, \|Δc_s\| **2.1e-17**, \|Δc_t\| **5.4e-15**, chi2 **1.8e-22** | 1e-12 |
| `contFit2` | e → 7e | ea, ecs, ect ×7 and the whole covariance ×49 to **<1e-12 rel** | exact scaling |
| `contFit2` | 400 noise draws | pull rms **0.995** | ≈ 1 |
| `jack` | AR(1) α=0.8, n=32768 | 2τ **8.212** (exact 9), 2τ_pos **8.425**, b=9, nblock=3641 | \|Δ\|<1.2 |
| `jack` | uncorrelated, n=16384 | 2τ **1.0023**, 2τ_pos **1.0089**, err/naive **1.0018** | 2τ≈1 |
| `jack` | uncorrelated, forced `bs=1` | err/naive − 1 = **5.3e-15** (delete-1 *is* the naive error) | naive |
| TSV | 33 rows × 3 cols incl. 1e-300, 1e300, −1/3, π | **bit identical**, 5 metadata keys preserved | bit identical |

### `nmaxFit` recovers the exact integer

Yes, in all six end-to-end cases. Correlators built as `C * sum_{n=0..N} (n+1)e^{-(n+1)t}/(4π)`
(fermion, t = at..6 with at = 12/168, 84 points) and
`C * sum_{n=1..N} sqrt(n(n+1))(2n+1)e^{-sqrt(n(n+1))t}/(8π)` (gauge, at = 12/120, 60 points),
multiplied by 1 + 1e-4·Gauss, searched over `1..40`:

| tower | true N | fit | \|ΔC/C\| | res/dof |
|---|---|---|---|---|
| fermion | 6 | **6** | 2.2e-6 | 9.8e-9 |
| fermion | 10 | **10** | 1.4e-5 | 7.3e-9 |
| fermion | 19 | **19** | 2.7e-5 | 8.1e-9 |
| gauge | 3 | **3** | 1.0e-5 | 1.1e-8 |
| gauge | 8 | **8** | 3.5e-6 | 8.7e-9 |
| gauge | 18 | **18** | 3.3e-6 | 9.7e-9 |

`res/dof ≈ (1e-4)^2` as it must be when the injected noise is relative and the weighting is too.

### Interface changes (also made in [`04-interfaces.md`](04-interfaces.md) §14)

The §14 sketch was three tuple-returning procs. All four fit results are now named objects
(`PlateauFit`, `LineFit`, `PlaneFit`, `NmaxFit`) plus `SeriesStat`, because every one of them
has to carry a covariance and a dof that a 3-tuple cannot hold. Specifically:

1. `plateauFit` gained `at`, `e`, `maxit` and returns `PlateauFit`, not `(d0, err, chi2)`.
   `at` is required — (V.6) is a fit in *physical* t, and `Delta'` and `c` are meaningless
   without it (`Delta_0` alone is scale-free, which is presumably why the sketch omitted it).
2. `contFit` returns `LineFit`; added `cab` and `dof` to the sketch's five fields.
3. `contFit2` and `nmaxFit` are new (they were in the WP-J brief but not in §14).
4. `jack` gained `stride` and `bs` and returns `SeriesStat`, not `(mean, err, tau)`.
   `bs` exists so a caller can override the automatic block size; the tests need it to
   exhibit the delete-1 identity.
5. `chi2dof(f)` covers all three fit objects.
6. `dataio` was not in §14 at all; both procs are new.

### Conventions chosen where the source documents do not fix them

* **Covariance normalization.** With per-point errors supplied, the covariance is
  `(J^T W J)^{-1}`, unscaled — so `contFit`/`contFit2` errors scale exactly with the input
  errors (tested). `plateauFit` called *without* errors uses unit weights and scales the
  covariance by `chi2/dof`, the usual unweighted-fit rule. Stated in the module docstring.
* **`nmaxFit` residual norm — UNVERIFIED, see below.**
* `contFit2` takes `abar_s^2` and `a_t^2`, not `abar_s` and `a_t`. Named `as2`/`at2` and
  documented; getting this wrong silently changes the meaning of `c_s`, `c_t`.
* `contFit2` **raises `ValueError`** on a rank-deficient design (all lattices at one `a_t`,
  say). The Gauss-Jordan inverse would otherwise return finite, plausible-looking garbage —
  the one place in this module where a user mistake fails silently rather than crashing.
  Tested. Everything else follows the repo rule and lets bad input crash or propagate NaN.

### What did not work / open items

1. **The published `residual/DOF` normalization of (V.9) is not determined by any document in
   this tree**, so `nmaxFit`'s default is a *choice*, not a reproduction. Default weight is
   `w_t = 1/g_t^2` (relative residual) — scale free, which is the only defensible default for
   deterministic correlators that carry no statistical error and span orders of magnitude over
   the fit range. Under that choice `res/dof = 0.028` (fermion L=1) means a 17% rms relative
   deviation and `0.0039` (L=4) means 6%. Those are plausible but not obviously right, and the
   flatness of the published gauge column (0.0031, 0.0023, 0.0031, 0.0037 across L=1..8, i.e.
   no improvement with L) does not look like a lattice-artifact-dominated residual, which is
   mildly suspicious. **WP-K must confirm against T1.4g / T1.5f before any of those eight
   numbers is claimed as reproduced.** An absolute (unweighted) residual is one line away:
   pass `e` = all ones. Any other weighting is likewise expressible through `e`.
2. **`jack`'s prescribed block size does not reach the asymptotic error.** With
   `blockSize = ceil(2 tau_int)` the measured blocked error on AR(1) α=0.8 is **2.253×** the
   naive one, not `sqrt(2 tau_int) = 3`. This is not a bug: the delete-b jackknife on AR(1)
   has `E[var] = var_naive · S(b)/b` with `S(b) = b + 2 Σ_{k<b} (b−k) ρ(k)`, which for b=9,
   ρ=0.8^k predicts **2.270** — the code agrees with that to 0.7%. Reaching sqrt(2τ) needs
   b ≫ 2τ (b ≈ 20τ costs ~5% ; b = 2τ costs 25%). The test asserts against the analytic
   prediction, not against sqrt(2τ), and the shortfall is documented in `jack`'s docstring.
   **Anyone quoting a Monte-Carlo error from WP-L should pass an explicit larger `bs`.**
3. **Wolff windowing underestimates 2τ by ~9%** here: 8.212 against the exact 9 at n=32768.
   Within the ±1.2 the brief asked for, but it is a systematic, not noise (the positive-sequence
   estimator gives 8.425, also low). Both are reported so the difference is visible.
4. **`ceil` quantization at 2τ ≈ 1.** An uncorrelated series measures 2τ = 1.0023, so
   `ceil` gives `blockSize = 2`, not 1, and the error is 1.0018× naive rather than exactly
   naive. Exact equality with the naive error is only reachable with `bs = 1` (verified to
   5.3e-15). Harmless, but it means "uncorrelated ⇒ block size 1" is not guaranteed.
5. `effMass` returns NaN wherever `c(t) < c(T/2)`, because `arccosh` of a number below 1 is
   NaN. This cannot happen for a real correlator that is symmetric about T/2, only for a
   monotone-decaying test input past the midpoint. Left to propagate, per the repo's
   floating-point rule.
6. **Bug found and fixed during the work:** the change of variables that maps the fitted
   amplitude back from the shifted time origin had the sign of the exponent wrong
   (`c = c' exp(-Δ' t_ref)` instead of `c = c' exp(+Δ' t_ref)`), which left `Delta_0` and
   `Delta'` correct but `c` wrong by `exp(2 Δ' t_ref)` — a factor of 3e4 on the paper's window.
   Only the noiseless three-parameter recovery test catches it; a test that checked `Delta_0`
   alone would have passed. Worth keeping that test.
7. The origin shift inside `plateauFit` is load bearing, not a micro-optimization. Without it
   the paper's 4 ≤ t < 8 window has `exp(-Δ' t) ∈ [2e-3, 3e-5]` across the whole fit range and
   `c` is ~1e4 times the signal; the normal equations lose most of their digits. With it the
   same window converges in 11 iterations to `|Δd0| = 0`.

---

## WP-C — analytic continuum correlators  (done)

Files: [`core/analytic.nim`](../core/analytic.nim), [`tests/tanalytic.nim`](../tests/tanalytic.nim).
Run: `cd build_mac && SDKROOT=$(xcrun --show-sdk-path) make run experimental/radial/tests/tanalytic`.
19 tests in 8 suites, all pass. Dependencies: `std/math` and `core/types` only — no `base`, no
`ops/`, so it compiles and tests with no other work package in place.

### Appendix C as published (doc/02 §7.4 only paraphrases it; recorded here so nobody re-fetches)

Read off arXiv:2510.03085v2. \(z=\cos\theta\), \(m\in\mathbb Z+\tfrac12\) for the fermion.

* (C.10) \(\xi_{|m|,n}(z)=(1-z)^{\frac12(|m|-1/2)}(1+z)^{-\frac12(|m|+1/2)}
  P^{(|m|-1/2,\,-|m|-1/2)}_{n+|m|+1/2}(z)\)
* (C.11) \(P^{(\alpha,\beta)}_j(z)\equiv\frac{\Gamma(\alpha+1+j)}{\Gamma(\alpha+1)j!}
  F(-j,\,1+\alpha+\beta+j;\,\alpha+1;\,\tfrac{1-z}2)\)
* (C.41) \(Q_m(z)=(1-z)^{\frac12(m-1)}(1+z)^{-\frac12(m+1)}\);
  (C.42) \(f_{|m|,n}=P^{(|m|-1,-|m|-1)}_{n+|m|+1}\);
  **(C.43) \(f_{0,n}(z)=(1-z)P^{(1,-1)}_n(z)\)**, \(m=0,\ n\ge1\)
* (C.54) \(\langle\psi(\theta,0)\bar\psi(0,0)\rangle=\dfrac{\sigma_1}{4\pi\sin(\theta/2)}\)
* (C.55) \(\langle\psi(\theta,0)\bar\psi(0,0)\rangle=\sigma_1\sum_{n\ge0}\dfrac{(-1)^n}{\pi\sqrt2}\,\xi_{1/2,n}(-z)\)
* (C.56) \(\frac1{g^2}\langle J^t(y)J^t(0)\rangle=\delta^2(y)-\frac1V\) (flat, volume \(V\))
* **(C.57)** \(\frac1{g^2}\langle J^t(\theta,0)J^t(0,0)\rangle
  =-\frac1{4\pi}\sum_{n\ge1}\frac{2n+1}{n+1}f'_{0,n}(z)=\frac1{2\pi}\Big[\delta(z-1)-\frac12\Big]\)

The Fig. 13 caption reads "\(1\le n\le n_{\max}\)" but (C.55) sums \(n\ge0\); **(C.55) is the
correct one** and the caption is loose. The \(n=0\) term is \(\sin(\theta/2)/2\pi\), which is 12 %
of the answer at \(\theta=0.5\) and 142 % of it at \(\theta=2\) — dropping it would leave a fixed
offset and the sum could not converge to \(1/(4\pi\sin(\theta/2))\) at all. Our sum starts at 0 and
does converge (below), which settles it.

### Implemented

* `jacobiP(j, alpha, beta, z)` — generalized real parameters, three-term recurrence (DLMF 18.9.2)
  seeded with \(P_0=1\), \(P_1=(\alpha+1)+(\alpha+\beta+2)(z-1)/2\).
* `fermionG(t, nmax)` — (V.3) (1,1) component; closed form \(\mathrm{sign}(t)/(16\pi\sinh^2(|t|/2))\)
  derived from \(\sum_{n\ge0}(n+1)x^{n+1}=x/(1-x)^2\), \(x=e^{-|t|}\).
* `fermionGPeriodic(t, T, nmax)` — alternating image sum
  \(\sum_j(-1)^j[h(t_r+jT)+h((j+1)T-t_r)]\), \(h=|G|\). The sign of the negative-\(t\) images cancels
  the antiperiodic \((-1)^k\), which is *why* the shape folds about \(T/2\) — that cancellation is
  the whole content of the boundary condition and is what makes a cosh fit legitimate.
* `gaugeG(t, nmax)` / `gaugeGPeriodic(t, T, nmax)` — (V.14) and its (unsigned) periodic images.
* `flatSpectrum(kap, kapT, n1, n2, nt)` — (IV.8), both branches.
* `s2FermionProp(theta, nmax)` — (C.55), literal (C.10) at \(|m|=1/2\).
* `s2CurrentCorr(z, nmax)` — (C.57), literal, via \(f'_{0,n}\).
* `effDim(g, at, T)` — (V.4)-(V.5).

### The generalized Jacobi polynomial at \(|m|=1/2\) — the one real derivation

At \(|m|=1/2\) the (C.10) exponents collapse to \((1-z)^0(1+z)^{-1/2}\) and the Jacobi indices become
\(j=n+1,\ \alpha=0,\ \beta=-1\). \(j\), \(\alpha+j=n+1\) and \(\beta+j=n\) are all non-negative
integers, so Rodrigues' formula — hence the ordinary recurrence — still applies even though
\(\beta<0\). Two consequences that had to be worked out:

1. **The recurrence, not the series.** (C.11) is a terminating \({}_2F_1\); at \(j\sim200\) with
   \(x=(1-z)/2\) near 1 its terms reach \(10^{20}\) and cancel to \(O(1)\), so it is useless in
   double precision. The recurrence is used instead. Its leading coefficient
   \(2n(n+\alpha+\beta)(2n+\alpha+\beta-2)\) vanishes iff \(\alpha+\beta=-n\) or \(\alpha+\beta=2-2n\);
   for the four parameter sets actually needed — \((0,-1)\) for \(\xi\), \((1,-1)\) and \((2,0)\) for
   \(f'_{0,n}\), \((0,0)\) for Legendre — the roots sit only at \(n=0,1\), which the explicit seeds
   cover. Both the seeds and the \(j=2\) value are checked against (C.11) in the test.
2. **\(\beta=-1\) is degenerate and that is a feature.** \(P^{(0,-1)}_{n+1}(z)=\frac{1+z}2P^{(0,1)}_n(z)\),
   so \(\xi_{1/2,n}(z)=\frac12\sqrt{1+z}\,P^{(0,1)}_n(z)\) — the \((1+z)^{-1/2}\) of (C.10) has no
   singularity at \(z=-1\), it is cancelled by a zero of the polynomial. Numerically this matters:
   the \(P_0=1\) seed enters the recurrence only through \(c_2=2(n+\alpha-1)(n+\beta-1)(2n+\alpha+\beta)\),
   which is exactly 0 at \(n=2\) for \(\beta=-1\), so an \(O(1)\) seed never contaminates the
   \(O(1-z)\) values near \(z=-1\) and relative accuracy survives down to \(\theta=10^{-3}\).
   Same story for (C.43): \(f_{0,n}(z)=(1-z)P^{(1,-1)}_n(z)=(1-z^2)P'_n(z)/n\), and the Legendre
   equation then gives \(f'_{0,n}(z)=-(n+1)P_n(z)\). Both reductions are pinned by the tests
   (n up to 40 / 60, \(z\) including \(\pm1\)), and the literal and reduced forms of both
   `s2FermionProp` and `s2CurrentCorr` are checked against each other.

### Independent confirmation of the (C.57) normalization

Substituting \(f'_{0,n}=-(n+1)P_n\) collapses (C.57) to \(\frac1{4\pi}\sum_{n\ge1}(2n+1)P_n(z)\) —
Legendre completeness with the \(\ell=0\) mode removed, i.e. \(\delta^2(x,y)-\frac1{4\pi}\).
That is reproduced from scratch: 2D Maxwell is \(S=\frac1{2g^2}\int\!\sqrt g\,f^2\) with
\(f=J^t=F_{\theta\varphi}/\sqrt g\), and on a closed surface with no monopole
\(\int\!\sqrt g\,f=0\), so \(f\) is a Gaussian field with the constant mode projected out:
\(\frac1{g^2}\langle ff\rangle=\delta^2-\frac1{4\pi}=\frac1{2\pi}[\delta(z-1)-\frac12]\).
The \(-1/2\) is the missing gauge zero mode, and the overall \(1/2\pi\) is fixed by
\(\delta^2=\frac1{2\pi}\delta(z-1)\) for an axially symmetric source. **The paper's normalization
is right**, which is the point of this work package.

### Measured — infinite volume

| check | measured | asked |
|---|---|---|
| `fermionG` closed form vs truncated, \(n_{\max}=200\), \(t=0.5,1,2,4\) | rel **0**, 1.9e-16, 1.2e-16, 1.4e-16 | 1e-14 |
| short distance \(\|G\,4\pi t^2-1\|\) at \(t=10^{-1..-4}\) | **8.3292e-4, 8.3333e-6, 8.3333e-8, 8.3333e-10** = \(t^2/12\) | → 1 |
| \(-d\ln G_g/dt\) at \(t=21\) | **1.414213563657** vs \(\sqrt2=1.414213562373\) | 1e-8 |
| `fermionGPeriodic` vs an independent thermal mode sum \(\frac{e^{-\lambda t}+e^{-\lambda(T-t)}}{1+e^{-\lambda T}}\) | **1e-14 rel**, both \(n_{\max}<0\) and \(n_{\max}=5,30,120\) | — |
| fold about \(T/2\), antiperiodicity, oddness | hold to 1e-13 rel (see caveat 4) | — |

\(G^{(1,1)}(t)\) reference values (this build): 3.1176050820561835e-01, 7.3264876746175961e-02,
1.4404749055764549e-02, 1.5124062502024948e-03 at \(t=0.5,1,2,4\).

### Measured — `effDim` and the two plateau values (T1.4e / T1.5d continuum targets)

`effDim` on an exact \(\cosh(\Delta(T/2-t))\), \(T=16\), \(a_t=0.1\), \(\Delta=1.234567\):
worst \(|\Delta_{\rm eff}-\Delta|=\) **1.577e-14** over \(t<T/2-1\) (asked 1e-12). Zero
discretization error is expected there — for a single state \(f\) is exactly linear — so this
validates the fitting convention itself, not an approximation.

On the analytic periodic correlators, \(T=40\), \(a_t=0.1\):

| \(t\) | \(\Delta_{\rm eff}\) fermion (→1) | \(\Delta_{\rm eff}\) gauge (→\(\sqrt2\)=1.414213562373) |
|---|---|---|
| 8 | 1.000638673853 | 1.414932089263 |
| 12 | 1.000011694049 | 1.414224986041 |
| 16 | 1.000000214109 | 1.414213744064 |
| 19 | 1.000000010601 | 1.414213570340 |

**Plateau over \(16\le t\le19\): fermion \(|\Delta_{\rm eff}-1|\le\) 2.141e-07, gauge
\(|\Delta_{\rm eff}-\sqrt2|\le\) 1.817e-07.** Both are asserted at 1e-6.

The residual is pure excited-state contamination and decays as \(e^{-(\Delta_1-\Delta_0)t}\):
successive fermion deviations give \(\ln(6.387\text{e-}4/1.169\text{e-}5)/4=1.0002\) against
\(\Delta_1-\Delta_0=1\), and the gauge ones give 1.03546 and 1.03518 against
\(\sqrt6-\sqrt2=1.035276\). That is a measurement of the *second* rung of each tower and is the
evidence behind the answer to open question 1 above.

It also sets the cost of the target, and this one matters for WP-K. A *raw* plateau at 1e-6 needs
\(t\gtrsim15\), i.e. \(T\gtrsim32\); \(T=40\) above was chosen for exactly that reason. At the
paper's \(T=16\) the raw \(\Delta_{\rm eff}\) deviation over the published fit window
\(4\le t<8\) runs

| \(t\) | 4 | 5 | 6 | 7 | 7.5 | best for \(t<T/2\) |
|---|---|---|---|---|---|---|
| fermion, \(L_t=168\) | 3.56e-2 | 1.36e-2 | 4.69e-3 | 1.66e-3 | 1.15e-3 | 1.22e-3 |
| gauge, \(L_t=120\) | 4.53e-2 | 1.48e-2 | 5.59e-3 | 2.09e-3 | 1.26e-3 | 1.41e-3 |

— three orders of magnitude short of the published \(\Delta_0^{\rm cont}=0.999998(34)\) and
\(1.41409(18)\). **So the (V.6) fit \(\Delta_{\rm eff}\simeq\Delta_0+c\,e^{-\Delta' t}\) is not
cosmetic: it is what removes the excited state, and T1.4e/T1.5d are simply unreachable by reading
a plateau at \(T=16\).** WP-K must fit, and must not read the ~1e-3 residual at \(t\approx8\) as a
lattice artifact — it is the \(\Delta_1\) tail of the continuum correlator and is present in the
analytic curve too. (Earlier draft of this section claimed \(T=16\) was intrinsically limited to
\(3\times10^{-4}\); that was an estimate, the measured numbers above replace it.)

### Measured — `flatSpectrum` (Fig. 5)

\(\kappa=1/\sqrt3\); grid 12×12×12; both anisotropy regimes.

| \(\bar a_s/a_t\) | Re range | first doubler | \(\min(4\kappa,2\kappa')\) |
|---|---|---|---|
| 2.0 (\(\kappa'=\sqrt3\)) | [0, **6.062177826491**] | **2.309401076759** | 2.309401076759 (= \(4\kappa\)) |
| 1.0 (\(\kappa'=\sqrt3/2\)) | [0, **4.330127018922**] | **1.732050807569** | 1.732050807569 (= \(2\kappa'\)) |

Conjugate pairing is exact by construction and asserted. Zero-imaginary-part modes sit at
\(4\kappa\) (at \((\pi,0),(0,\pi),(\pi,\pi)\)), \(4.5\kappa\) (at the K points
\((2\pi/3,2\pi/3)\), \((4\pi/3,4\pi/3)\)), \(2\kappa'\) and their sums.

### T1.6a — \(S^2\) fermion propagator (C.55) → (C.54)

**The series is only conditionally convergent, so "the pointwise error decreases with
\(n_{\max}\)" is false and the test does not assert it.** With the reduction above the sum is
\(\frac{\sin(\theta/2)}{2\pi}\sum_n(-1)^nP^{(0,1)}_n(-z)\), and the Jacobi generating function
\(\sum_nP^{(\alpha,\beta)}_nt^n=2^{\alpha+\beta}R^{-1}(1-t+R)^{-\alpha}(1+t+R)^{-\beta}\),
\(R=\sqrt{1-2xt+t^2}\), at \(t=-1\), \(x=-\cos\theta\) gives exactly \(1/(2\sin^2(\theta/2))\) —
the CFT answer, but *on* the boundary of the disc of convergence. Partial sums therefore oscillate
with an \(O(n_{\max}^{-1/2})\) envelope. What is asserted is the envelope over one oscillation
period \(\lceil2\pi/\theta\rceil\) in \(n\), and a uniform bound.

| \(\theta\) | exact | pointwise @200 | envelope-rel 50 → 100 → 200 | \(C=|{\rm err}|\sin(\theta/2)\sqrt{n_{\max}}\) |
|---|---|---|---|---|
| 0.5 | 0.32164995 | 0.29963195 | 1.582e-1 → 1.116e-1 → **8.079e-2** | 0.089, 0.089, 0.091 |
| 1.0 | 0.16598505 | 0.15642306 | 1.132e-1 → 8.501e-2 → **5.761e-2** | 0.064, 0.068, 0.065 |
| 2.0 | 0.09456948 | 0.08901814 | 1.085e-1 → 8.264e-2 → **5.870e-2** | 0.061, 0.066, 0.066 |

The envelope is strictly decreasing in every case and \(C<0.15\) uniformly, i.e. the error obeys
\(|{\rm err}|\lesssim0.15/(\sin(\theta/2)\sqrt{n_{\max}})\) — that is the honest statement of
"(C.55) → (C.54) as \(n_{\max}\to\infty\)". Ultraviolet restoration, the actual message of Fig. 13,
is separately asserted: envelope-rel at \(\theta=0.05\) falls 4.028e-1 → 2.498e-1 and at
\(\theta=0.02\) falls 7.519e-1 → 3.950e-1 going from \(n_{\max}=50\) to 200, and the truncated sum
is still below 2 % of the exact answer at \(\theta=10^{-3}\), \(n_{\max}=50\). The cutoff scales as
\(\theta\sim\pi/n_{\max}\), matching the figure.

### T1.6b — \(S^2\) current correlator (C.57)

For a polynomial \(\varphi\) of degree \(\le n_{\max}\) the truncated sum integrates **exactly**,
so these moments pin both pieces of (C.57) with no truncation error at all (96-node Gauss-Legendre,
itself verified exact on \(z^k\), \(k\le40\), to 1e-13):

| moment | measured (\(n_{\max}=5,20,40\)) | exact |
|---|---|---|
| \(\int\!\langle J^tJ^t\rangle\,dz\) | −3.2e-16, −4.6e-15, −1.1e-14 | 0 |
| \(\int z\,\langle\cdot\rangle\,dz\) | 0.159154943091895 / …891 / …885 | \(1/2\pi=0.159154943091895\) |
| \(\int z^2\langle\cdot\rangle\,dz\) | 0.106103295394597 / …592 / …586 | \(1/3\pi=0.106103295394597\) |
| \(\int(1-z)\langle\cdot\rangle\,dz\) | −0.159154943091895 / …896 / …896 | \(-1/2\pi=-\frac1{4\pi}\!\int\!(1-z)dz\) |

The last row is the requested **constant offset check**: \(\varphi(1)=0\) kills the delta, so the
whole integral is \(-\frac1{4\pi}\int\varphi\) and it comes out to 1e-15.

Delta-function development against smooth test functions, \(n_{\max}=20\to40\):

| test function | exact | err(20) | err(40) |
|---|---|---|---|
| \(\exp(-\frac12((z-1)/0.15)^2)\) (smeared delta, width 0.15) | 0.1441946076 | 2.572e-05 | **4.101e-11** |
| \(1/(1.05-z)\) (pole just outside the cut) | 2.8875821864 | 1.446e-02 | **3.624e-05** |
| \(\exp(z)\) (entire) | 0.2455889106 | 1.207e-14 | 2.984e-14 |

Convergence is geometric in \(n_{\max}\), set by the Bernstein ellipse of the test function;
\(\exp(z)\) is already at the roundoff floor at \(n_{\max}=20\) (its Legendre coefficients fall
super-exponentially), so for that row the test only demands "no worse than the floor".

### Figure data

`output/radial/analytic/fig13_s2_fermion.tsv` — 401 log-spaced \(\theta\in[10^{-3},3]\), columns
`theta, exact_C54, nmax50_C55, nmax200_C55`.
`output/radial/analytic/fig14_s2_current.tsv` — 2001 points \(z\in[-1,1]\), columns
`z, nmax20_C57, nmax40_C57`.
Spot checks against the published figures: at \(z=1\) the Fig. 14 curves are 35.0140875 and
133.6901522, i.e. \(\frac1{4\pi}\sum_{n=1}^{N}(2n+1)=\frac{(N+1)^2-1}{4\pi}\) for \(N=20,40\);
Fig. 13 at \(\theta=0.0074\) gives 0.774 (\(n_{\max}=50\)) and 10.39 (200) against exact 21.5, and
at \(\theta=0.020\) gives 1.99 and 11.01 against 7.91 — the crossing and the roll-off of the
\(n_{\max}=200\) curve near \(\theta\approx0.012\) match the published plot.

### Interface changes (mirrored in [`04-interfaces.md`](04-interfaces.md) §6)

Every §6 signature is unchanged. Added and clarified:

1. **`jacobiP(j, alpha, beta, z)` is new and exported.** (C.10)/(C.43) need generalized parameters
   and QEX has nothing. `jacobiP(l, 0, 0, x)` is the Legendre polynomial.
2. **No `legendreP` was added.** Nothing in WP-C needs the associated form; the brief allowed it
   "if you need them" and adding an unused public name would violate the repo's small-API rule.
   WP-I should put `legendreP` next to `ylm` in `meas/harmonics.nim`.
3. `nmax < 0` is documented as "untruncated" for **all four** correlators, not just the fermion:
   closed form where one exists, otherwise sum to underflow. `gaugeG`'s default stays 200.
4. `effDim` result length (`g.len-1`), index convention (`g[i]` at \(t=i\,a_t\)) and the
   \(G(T/2)\) index (`round(T/(2 at))`) are now written down.
5. `flatSpectrum` momentum grid (\(2\pi i/n\), \(i=0..n-1\)) and output ordering (both branches
   consecutive, `2*n1*n2*nt` entries) are now written down.
6. Declared `func`, not `proc` — same call signatures, matches `core/types.nim`.

### What did not work, and things to be careful about

1. **The WP-C brief's "real part in \([0,4\kappa+2\kappa']\)" is wrong; the true bound is
   \(4.5\kappa+2\kappa'\).** The maximum of \(3-\cos k_1-\cos k_2-\cos(k_1{+}k_2)\) is \(9/2\) at the
   K point \(k_1=k_2=2\pi/3\) (stationarity forces \(\sin k_1=\sin k_2=-\sin(k_1{+}k_2)\), i.e.
   \(\cos k=-1/2\)), not 4 at \((\pi,0)\). The test asserts equality with \(4.5\kappa+2\kappa'\) to
   1e-12 rather than the looser inequality. \(4\kappa\) is the *first doubler*, a different
   statement, and that one is separately verified and does hold.
2. **T1.6a cannot be stated as "the error at fixed \(\theta\) decreases from \(n_{\max}=50\) to
   200".** It does not: at \(\theta=2\) the raw relative error goes 1.3e-2 → 3.4e-2 → 5.9e-2 over
   \(n_{\max}=50,100,200\) purely from the oscillation phase. Only the envelope decreases. The
   tolerance was not loosened to hide this — the assertion was changed to the quantity that is
   actually monotone, and the raw pointwise numbers are echoed by the test so the oscillation is
   visible. A Cesàro average over one period converges an order of magnitude better (rel 8.5e-3 →
   3.8e-3 at \(\theta=0.5\) for \(n_{\max}=50\to200\)) if anyone needs a smoother estimator;
   it is not exported because the paper plots the raw partial sum.
3. **`effDim` here and `effMass` in [`meas/fit.nim`](../meas/fit.nim) (WP-J) are the same
   operation**, both mandated by §6 and §14 of the interface doc. They agree by construction
   (both form the ratio by division so \(f(T/2)\equiv0\)). WP-K should pick one — `effMass` if it
   is already fitting, `effDim` if it only wants the analytic curve — and a later cleanup should
   collapse them.
4. **Floating point, not algebra:** the periodic identities \(G_T(T-t)=G_T(t)\) etc. hold only to
   \(\sim\epsilon T/t\), because `T - t` cancels and \(h\sim1/u^2\) doubles the relative error.
   At \(t=0.35\), \(T=12\) that is 1.9e-15, so the tests use 1e-13, with the reason in a comment.
5. **Contact singularities are left to propagate**, per the repo's floating-point rule:
   `fermionG(0)` and `fermionGPeriodic(0,T)` return a non-finite value (\(1/\sinh^20\)), and
   `gaugeG(0, -1)` does not converge (it stops at the \(10^6\)-term cap). Callers that build a full
   time slice array must not use index 0; the plateau test overwrites it and says so.
6. **No external oracle exists for (C.55)/(C.57).** The closed forms they must reproduce *are* the
   oracle — the CFT propagator, the delta function, and the exact polynomial moments — plus the
   independent derivation of (C.57) above and the independent thermal mode sum for the periodic
   fermion. Everything reported here is measured against one of those, never against the code's
   own output.
7. **Unsure:** whether the paper's \(n_{\max}\) fits (T1.4f/T1.5e) intend the truncation applied to
   the *infinite-volume* tower and then imaged, or to the finite-\(T\) thermal sum. They are
   algebraically identical (proved by swapping the two sums, and checked numerically to 1e-14 for
   \(n_{\max}=5,30,120\)), so it does not matter here — but WP-K should not assume the same is true
   once a lattice correlator with cutoff effects replaces the analytic one.

---

## WP-D — spinor fields, CG, multishift CG  (done)

Files: [`core/spinor.nim`](../core/spinor.nim), [`ops/solve.nim`](../ops/solve.nim),
[`tests/tspinor.nim`](../tests/tspinor.nim), [`tests/tsolve.nim`](../tests/tsolve.nim).
Run: `cd build_mac && SDKROOT=$(xcrun --show-sdk-path) make run experimental/radial/tests/tspinor experimental/radial/tests/tsolve`.
12 tests, all pass; `tsolve` takes ~8 s. Dependencies: `core/types`, `rng/threefry4x64`
(`tspinor` also `base/alignedMem`, `tsolve` also `eigens/linalgFuncs`). **`core/lattice.nim` is
not imported** — fields are sized by a plain `int`, so this compiled and tested before WP-A landed.

### Implemented

* `spinor.nim`: `newSpin(n)`, `zero`, `:=`, `axpy` (float and Complex64), `axpby`, `scale`, `dot`,
  `redot`, `norm2`, `gaussian`, `pointSource(n, site, comp)`. All allocation-free except the two
  constructors, all written on `.re`/`.im` directly so `-Ofast` vectorizes them.
* `solve.nim`: `cgSolve`, `cgmSolve` (Jegerlehner hep-lat/9612014, rebased on `shifts[0]`),
  plus a private `cgRun` kernel shared by CG and the refinement passes, and a private `resid2`
  that forms `|b - (A+s)x|^2` with no temporary. `cgmSolve` allocates `nshift+3` scratch fields
  and two `float` arrays at entry; nothing allocates inside the loop.

### Interface changes (mirrored in [`04-interfaces.md`](04-interfaces.md) §5 and §9)

1. `newSpin(n: int)` and `pointSource(n, site, comp)` take a plain length / flat index instead of
   a `Lat`. **`lattice.nim` owns the `newSpin(l: Lat)` and `pointSource(l, v, t, c)` overloads.**
2. `axpby(x, a, y, b)` (`x = a*y + b*x`) is new and exported — `solve.nim` needs the CG
   search-direction update and it does not belong in a second copy there.
3. `CgInfo.r2` and `MultiCgInfo.r2` are **relative** (`|b-(A+s)x|^2 / |b|^2`) and are **always
   recomputed** from `x`, never the recursive estimate. That makes the §9 guard
   `r2 <= 1.001*r2req` read literally, and costs `cgSolve` one extra operator apply per solve.
4. `MultiCgInfo` gained `r2pre` (the same true residual measured *before* refinement — how far
   the recurrence was off) and `refits` (iterations spent in refinement). `iters` is the
   recurrence count only, so it stays comparable to a single-shift `cgSolve`.
5. `cgmSolve` sizes `xs` itself if it does not already match `(nshift, b.len)`, and `doAssert`s
   that the shifts are positive and strictly ascending — a descending list would silently return
   garbage for the far shifts, which is the one class of misuse worth a runtime check here.

### Measured — easy operator (dense Hermitian, dim 120, spectrum log-spaced in [0.05, 5], cond 100)

The reference inverse is a `zeigs` (zheev) eigendecomposition **of the assembled matrix**, so the
oracle inverts exactly the matrix the solver sees.

| check | measured | asked |
|---|---|---|
| `cgSolve` vs dense inverse, `r2req=1e-26` | 119 iters, r2 = **8.972e-27**, \|dx\|/\|x\| = **7.101e-14** | 1e-9 |
| **T1.3h** cgm (15 shifts, 1e-3…3.0) vs 15 independent `cgSolve` | worst \|dx\|²/\|x\|² = **3.729e-27** | 1e-18 |
| T1.3h, same in amplitude | worst \|dx\|/\|x\| = **6.1e-14** | — |
| every shift's recomputed true residual | ≤ **9.1e-27** | 1.001e-26 |
| every shift vs the dense inverse | ≤ **8.4e-14** | — |
| iteration count, cgm vs smallest-shift `cgSolve` | **119 vs 119** | close |
| determinism, `r2req=1e-24` (no refinement) and `1e-32` (14 of 15 refined, 1265 refits) | **bit identical** in every component and every reported field | bit identical |

Shift 0 agrees with its independent `cgSolve` to **exactly 0.0** — see the base-refinement note below.
Independent single-shift iteration counts run 119 (s=1e-3) down to 21 (s=3.0); the multishift gets
all 15 for the cost of the hardest one, which is the point of T1.3h.

### Measured — hard operator (dim 400, spectrum log-spaced in [1e-8, 1], cond 1e8, 8 shifts spanning 1e6)

Base system is `A + 1e-6`, cond ≈ 1e6. At `r2req = 1e-16`, `maxits = 20000`:

* 4879 recurrence iterations, `converged = true`, **`refined = 0`**;
* every recomputed true residual ≤ **0.989 × r2req**, worst *pre*-refinement value also 0.989;
* every solution agrees with the dense inverse to ≤ **3.3e-9** relative.

Scanned over 24 tolerances from 3e-16 down to 4e-19 (7 shifted systems each, 168 samples): the
pre-refinement true residual of a **shifted** system never once exceeded `1.001*r2req`, and the
values cluster right under 1.0 because a frozen shift stops exactly when its estimate crosses the
target. So the recurrence estimate is good to **better than 0.1%** here.

### The one real bug, and the fix

The textbook Jegerlehner recurrence, applied to every shift at every iteration, **produces NaN**.
`z_j` decays like `prod 1/(1 + sg*alpha_i)`; for a far shift with an ill-conditioned base
`sg*alpha ~ 1e8` per step, so `z_j` underflows to 0 within a few iterations of converging, and then
`zn = z_i z_{i-1} alpha_{i-1} / d` is `0/0`. Measured on the hard operator at `r2req = 1e-16`:
the four largest shifts came back `NaN`. It was invisible until `r2pre` was added, because the
refinement fallback silently re-solved all four (`refined=4, refits=231, converged=true`) — the
"refinement fires and fixes it" that the WP-D brief predicted was **this bug being papered over**,
not the far-shift estimate drifting.

Fix: freeze shift `j` as soon as `z_j^2 r2 <= r2req |b|^2`, i.e. as soon as its own residual meets
the target. Standard practice, and here it is required for correctness, not speed. It also drops
the shifted vector work substantially (in the easy case most of the 15 shifts are frozen well
before iteration 119). Pinned by `check mi.refined == 0` in the hard-case test: a regression to the
unfrozen recurrence puts NaN back and fires the fallback again.

Second finding: **refining `j = 0` is a no-op.** The recurrence *is* plain CG on `A + shifts[0]`,
operation for operation, so `cgRun` from `x = 0` reproduces `xs[0]` bit for bit — confirmed both
directly (pre = post to every digit) and by the `|dx|²/|x|² = 0.000e+00` in the T1.3h table.
It used to cost a full duplicate solve (5210 wasted iterations in one scan). `cgmSolve` now skips it.

### What did not work — the refinement fallback does not fire in a well-posed solve

**I could not construct an exact Hermitian positive-definite operator for which the fallback fires
*and* repairs a shifted system, and I do not believe one exists once shifts are frozen.** The brief
asked for that case; here is why it is not reachable, with the measurements:

* The shifted systems are all *easier* than the base. Their error floor is `~eps|b|` (the
  amplification `|A+s_j|` cancels against `|x_j| ~ |b|/s_j`), i.e. `r2 ~ 1e-32`, whereas the base's
  floor is `~(eps*cond)^2 ~ 1e-18` for cond 1e6. Between those two there is no tolerance at which
  a far shift fails but the base succeeds. Measured: over the 24-tolerance scan, the only system
  ever above `1.001*r2req` was the base.
* Truncating with `maxits` does not help either: the systems still live at the cap are exactly the
  hard ones, and their refinement gets the same budget. Scanned 12 caps from 2400 down to 350 —
  1 to 5 shifts missed each time, **0 repaired**. (A single lucky point, 8 shifts at `maxits=800`,
  did repair one shift, 5.88 → 0.31 × r2req, off the numerical-quality gap between the recurrence
  and a dedicated CG. That is a coin flip, so it is not what the test asserts.)

What the test asserts instead, `tests/tsolve.nim` "refinement fires when the request is below the
roundoff floor": hard operator, `r2req = 1e-24`, i.e. six decades under the base's floor. The
fallback fires for 2 shifts (3863 refit iterations), improves both (1379 → 885 and 12.1 → 9.9
× r2req), cannot reach a tolerance the operator cannot support, and `cgmSolve` reports
`converged = false` rather than believing the estimate. The `determinism` test drives the same
branch on the cheap operator (`r2req = 1e-32`, 14 of 15 shifts refined) and shows it is bit
reproducible. So the path is covered and deterministic — it is just not reachable from a
well-posed request, which is the honest statement.

Consequence for WP-F: the fallback is a **safety net, not a routine cost**. Budget for it firing
occasionally (the freeze lands the estimate anywhere in `(0, 1] × r2req`, so a shift at 0.999 with
a 0.2% truth-vs-estimate gap will trip the 1.001 guard) but not systematically.

### Other things to know

1. **The `1.001` acceptance guard is not enough near an operator's roundoff floor.** At
   `r2req = 1e-18` on the cond-1e6 base the recomputed residual came out at **1.063 × r2req**, and
   at 9.5e-18 at 1.004 — one of 24 tolerances scanned. The guard's stated rationale ("recursive and
   recomputed differ at the roundoff floor") is calibrated for the 4e-5 relative gap quoted in §9;
   with cond 1e6 the gap is 100× larger, because the CG loop stops on the *recursive* residual while
   the report is the *true* one. This is not fixed here — the §9 stopping rule and guard are
   normative and I did not change them. **WP-F must keep `r2req` at least ~2 decades above
   `(eps * cond(X†X + c_1))^2` or `cgSolve`/`cgmSolve` will report honest failures on good solves.**
   If that turns out to be impractical, the fix is a true-residual restart inside `cgRun` (continue
   iterating instead of returning), which stays deterministic and reversible.
2. **The `getRawMemAllocated()` allocation regression is vacuous** and the test says so: a `Spin`
   is a plain Nim `seq`, so QEX's `alignedMem` counter is 0 before and after (printed by the test).
   The check with teeth is `getOccupiedMem()`, also asserted unchanged across 64
   `axpy`/`scale`/`axpby`/`dot`/`norm2`/`redot` rounds (3398752 → 3398752).
3. **`cgmSolve` allocates `nshift+3` fields per call.** The §9 signature is stateless, so there is
   nowhere to cache them. For the HMC inner loop WP-F should hoist this into `Ov.work` or accept
   the churn; on a 10³–10⁴-site lattice with 15 poles it is ~18 seqs per call.
4. **`rng/threefry4x64` is imported directly**, as the brief said — the `rng` umbrella drags in
   `field`. Note `seed()` broadcasts through `defaultComm` and therefore needs `qexInit`;
   `seedIndep()` does not, which is why the tests use it and why they need no `qexInit`.
5. Cross-proc bitwise comparisons are not safe under `-ffast-math`: `dot(x,y).re` and `redot(x,y)`
   are the same arithmetic in the same order but the compiler reassociates each reduction
   independently. `tspinor` compares them at `1e-13 * sqrt(norm2(x)*norm2(y))`. The exactness tests
   that *are* bitwise (`axpy`, `scale`, `axpby`) use dyadic-rational inputs where no rounding
   occurs at all, so FMA contraction cannot change the answer.
6. RNG spot check: `gaussian` gives `<|x_i|^2> = 1.00567` over 40000 complex components
   (1 ± 0.005), and streams `seedIndep(s,0)` vs `seedIndep(s,1)` correlate at 0.0082 against an
   expected 0 ± 0.0035 — 2.3σ, unremarkable, but recorded in case anyone sees it drift.
7. `zgesv` is indeed not bound; the dense oracle uses `zeigs` (`zheev`) twice, once to make a
   random unitary and once to decompose the assembled matrix. No LAPACK bindings were added.

---

## WP-A — geometry and lattice  (done)

Files: `core/geom.nim`, `core/lattice.nim`, `tests/tgeom.nim`, `rgeom.nim`.

```bash
cd build_mac
SDKROOT=$(xcrun --show-sdk-path) make run experimental/radial/tests/tgeom   # 27 tests, all pass
SDKROOT=$(xcrun --show-sdk-path) make run experimental/radial/rgeom         # the table below
```
`rgeom -lev:16` also works (0.30 s); `-lev:8` takes 0.048 s. No QEX `Layout`/`Field` anywhere;
`geom.nim` re-exports `types.nim` and `lattice.nim` re-exports `geom.nim`, so
`import core/lattice` is the single entry point.

### What is implemented

* Refined icosahedron: 12 vertices = cyclic permutations of \((0,\pm1,\pm\phi)/\sqrt{1+\phi^2}\),
  20 outward-oriented faces, \(v_{ijk}=\mathrm{unit}(iA+jB+kC)\), deduplicated on a \(10^{-9}\)
  hash grid with a 27-cell probe (a plain single-cell hash can split a point across a cell
  boundary; the probe removes that failure mode entirely).
* **Everything is intrinsic/spherical**: \(\ell=\) geodesic arc, \(A_\triangle=\) spherical excess
  via van Oosterom–Strackee written as \(a\cdot((b{-}a)\times(c{-}a))\) (no cancellation for small
  triangles), dual point \(=\) projected circumcenter, \(\ell^*=\arcsin(\hat c\cdot\hat n)\) signed,
  \(A_y=\) spherical polygon area over the incident circumcenters, \(\tilde A_i\) the quadrilateral
  \((y_i,m_i,\hat c,m_{i-1})\) split at \(\hat c\).
* Counterclockwise vertex rings `nbr`/`nbe`/`nbf`, deterministic (start from the lowest-index
  neighbour, walk the "next neighbour in the face" map).
* Chart frame at both ends of every edge: `ea`, `eb`.
* Spin connection by two independent constructions (below).
* `Lat` and all free-limit couplings (IV.2), (IV.26), (IV.10).

### `rgeom` output (defaults: `lev=8 nt=32 at=0.2 tilt=0.7`)

```
L    N_V    N_E    N_F  V-E+F        abar         a/at        maxM
1      12     30     20      2  1.1071487178     5.535744    2.078461
2      42    120     80      2  0.5909464448     2.954732    2.078461
4     162    480    320      2  0.2994744727     1.497372    2.078461
8     642   1920   1280      2  0.1502274914     0.751137    1.170907

L        min A_y      max A_y     mean A_y        min l        max l       min l*  #A_y vals  #l vals
1    1.047197551  1.047197551  1.047197551  1.107148718  1.107148718  0.364863828        1        1
2    0.273844218  0.309341333  0.299199300  0.553574359  0.628318531  0.125956662        2        2
4    0.058328301  0.083977086  0.077570189  0.253865759  0.326366222  0.050450066        4        6
8    0.012965501  0.022938609  0.019573786  0.119529734  0.164833703  0.023033618       10       20

identities (max residual over the whole lattice)
L    |sumAtri-4pi| |sumAy-4pi|   dual-exact   sum(subA)    subA->Ay     holonomy    global-hol   omega-oracle   Z2  defect/2pi  nf
1      0.000e+00    5.329e-15    1.110e-16    2.220e-16    0.000e+00    3.140e-16    1.119e-16      8.882e-16    2.0 -2.000000    2
2      3.553e-15    7.105e-15    5.551e-17    5.551e-17    1.110e-16    4.041e-16    7.155e-16      1.332e-15    2.0 -2.000000    2
4      3.020e-14    1.954e-14    4.857e-17    4.163e-17    1.804e-16    3.445e-16    7.554e-16      9.992e-16    2.0 -2.000000    2
8      4.796e-14    3.730e-14    2.776e-17    1.908e-17    3.469e-16    6.292e-16    2.217e-15      1.332e-15    2.0 -2.000000    2

O(abar^2) residuals, relative to A_tri (must fall by ~4 per level doubling)
L     closure off  closure diag   flat sum(l l*/2)   poleGap site/link
1      4.770e-16     3.562e-02          3.562e-02          0.2933 0.1666
2      1.907e-03     1.272e-02          1.272e-02          0.2933 0.1666
4      3.777e-04     3.542e-03          3.534e-03          0.1670 0.0817
8      1.189e-04     9.107e-04          9.083e-04          0.0818 0.0366

flat equilateral limit (kappa -> 1/sqrt3 = 0.57735, kappa'/(abar/at) -> sqrt3/2 = 0.86603)
L    mean kappa    -1/sqrt3   mean kappa'/(a/at)   -sqrt3/2    Var(l)/abar^2   max |Ay/((sqrt3/2) ay^2)-1|    kappa range
1      0.659105   8.176e-02             0.854312  -1.171e-02       3.620e-31                 1.353e-02   0.65911-0.65911
2      0.597896   2.055e-02             0.856770  -9.255e-03       3.999e-03                 3.186e-02   0.53248-0.66331
4      0.582551   5.201e-03             0.864919  -1.107e-03       4.564e-03                 4.506e-02   0.45471-0.61252
8      0.578669   1.318e-03             0.867313   1.288e-03       5.031e-03                 4.787e-02   0.37082-0.62676

flat cross-check on a regular hexagonal patch of polar radius d (errors must go like d^2)
d           l*/a-1/(2sqrt3)  A_tri/a^2-sqrt3/4   A_y/a^2-sqrt3/2   kappa-1/sqrt3   kappa'/(a/at)-sqrt3/2
1.e-01         -1.979e-05          7.220e-04         7.232e-04      -3.958e-05               7.232e-04
1.e-02         -2.004e-07          7.217e-06         7.217e-06      -4.009e-07               7.217e-06
1.e-03         -2.009e-09          7.217e-08         7.214e-08      -4.018e-09               7.214e-08

site coordination and abar
L=1   5-fold    12  6-fold     0  other     0   abar = 1.107149
L=2   5-fold    12  6-fold    30  other     0   abar = 0.590946
L=4   5-fold    12  6-fold   150  other     0   abar = 0.299474
L=8   5-fold    12  6-fold   630  other     0   abar = 0.150227
```

### T1.1 acceptance, measured

| # | criterion | measured | verdict |
|---|---|---|---|
| T1.1a | \(N_V,N_E,N_F\), Euler | exact integers, \(V-E+F=2\) for L=1,2,4,8 (and 16) | **pass** |
| T1.1b | \(\sum A_\triangle=\sum A_y=4\pi\) | \(\le4.8\times10^{-14}\) / \(3.7\times10^{-14}\) | **pass** (\(<10^{-12}\)) |
| T1.1c | dual decomposition of every triangle | exact spherical form \(\le1.1\times10^{-16}\); **flat form \(\sum\ell\ell^*/2\) is only \(O(\bar a_s^2)\)** — see caveat 1 | **pass as restated** |
| T1.1d | closure (IV.6) | rel. off-diag \(4.8\!\cdot\!10^{-16}\), \(1.91\!\cdot\!10^{-3}\), \(3.78\!\cdot\!10^{-4}\), \(1.19\!\cdot\!10^{-4}\); rel. diag \(3.56\!\cdot\!10^{-2}\), \(1.27\!\cdot\!10^{-2}\), \(3.54\!\cdot\!10^{-3}\), \(9.11\!\cdot\!10^{-4}\) | **pass**, ratios below |
| T1.1e | \(\prod\Omega=\exp(i\sigma_3A_\triangle/2)\) every face | \(\le6.3\times10^{-16}\), both tilts | **pass** (\(<10^{-12}\)) |
| T1.1f | global holonomy \(=1\) | \(\le2.2\times10^{-15}\) | **pass** (\(<10^{-11}\)) |
| T1.1g | \(I_h\) orbits | exactly 12 five-fold, \(N_V-12\) six-fold; distinct \(A_y\) = 1, 2, 4, 10 = \(\mathrm{round}((L{+}3)^2/12)\) | **pass** |
| T1.1h | \(\bar a_s(L)\) | table below | **pass** |
| T1.4i | Table I doubler membership, T=16 | all 20 cells reproduced; \(L_t^{\min}=19.27,\,36.10,\,71.24,\,142.01\) | **pass** |

Extra identities not in T1.1, all pinned in `tgeom`:
\(\sum_i\tilde A_i=A_\triangle\) \(\le2.2\!\cdot\!10^{-16}\); \(\sum_{\triangle\ni y}\tilde A=A_y\) \(\le3.5\!\cdot\!10^{-16}\);
\(\sum\ell=N_E\bar a_s\) \(<10^{-12}\); \(\omega\) vs `omegaChart` mod \(2\pi\) \(\le1.3\!\cdot\!10^{-15}\);
parity \(\omega(Py_1,Py_2)=-\omega(y_1,y_2)\) \(<10^{-13}\); `eb` = `ea` rotated by \(\omega\) \(<10^{-13}\);
\(\bar a_s/a_t=5.54,\,2.95,\,1.50,\,0.75\) at \(a_t=0.2\); \(\kappa\) and \(A_y\) ranges pinned against the
Python oracle below.

### \(\bar a_s\): Nim vs the pure-Python oracle

| L | Nim (10 digits) | Nim (6 digits) | oracle above | agree |
|---|---|---|---|---|
| 1 | 1.1071487178 | **1.107149** | 1.107149 | ✓ |
| 2 | 0.5909464448 | **0.590946** | 0.590946 | ✓ |
| 4 | 0.2994744727 | **0.299474** | 0.299474 | ✓ |
| 8 | 0.1502274914 | **0.150227** | 0.150227 | ✓ |
| 16 | 0.0751747205 | 0.075175 | — | — |

\(\bar a_s(1)=2\arcsin\sqrt{(5-\sqrt5)/10}\) to \(10^{-14}\). \(\min/\max A_y\) and the \(\kappa\)
ranges also reproduce the oracle to all quoted digits (0.65911; 0.53248–0.66331;
0.45471–0.61252; 0.37082–0.62676) — pinned as a test.

**Note for `02-formulation.md` §6**, which lists \(\bar a_s\approx0.55,\,0.276,\,0.138\) for
L=2,4,8: those are just \(1.1071/L\), the flat "divide each edge into L pieces" estimate. Radial
projection makes the real mesh coarser by ~7–9 %. The correct values are the ones above, and they
are what change L=4 from marginal to comfortable at \(a_t=0.2\) (1.497 vs 1.38).

### Closure relation (IV.6) — convention and \(O(\bar a_s^2)\) ratios

`checkClosure` takes the geodesic tangent at each **edge midpoint**, parallel-transports it to the
triangle's dual point, and reads all three in one frame there. Ratios of successive levels
(should be 4 for \(\bar a_s^2\)):

| | L=1→2 | 2→4 | 4→8 | 8→16 |
|---|---|---|---|---|
| relative off-diagonal | (L=1 is exactly 0 by the 3-fold symmetry) | 5.05 | 3.18 | 3.85 |
| relative \|diag − \(A_\triangle\)\| | 2.80 | 3.59 | 3.89 | 3.97 |

**Discrepancy with the Python oracle, unresolved.** The oracle's table above gives 2.68e-1,
1.17e-1, 3.46e-2, 9.03e-3 for the diagonal; mine are ~10× smaller and both fall like \(\bar a_s^2\).
At L=1 the value is not a matter of convention: the icosahedral face has an exact 3-fold symmetry
about its dual point, so any construction using the three **unit** tangents in one frame has them
at exactly 120°, giving \(\sum_i e^a_ie^b_i=\tfrac32\delta^{ab}\) and therefore
\[
\frac{|\tfrac32\ell\ell^*-A_\triangle|}{A_\triangle}
=\frac{|1.5\times1.107148718\times0.364863828-0.628318531|}{0.628318531}=3.5620\times10^{-2},
\]
which is exactly what the Nim code reports. The oracle's 2.68e-1 cannot come from unit tangents,
so the two are measuring different things; whoever owns the oracle should re-derive it. Nothing
downstream depends on the coefficient — T1.1d only asks for \(O(\bar a_s^2)\), decreasing.

### Spin connection — conventions, stated once and unambiguously

**Sign.** With the face traversed counterclockwise as seen from **outside** the sphere,
\[
\Omega_{12}\Omega_{23}\Omega_{31}=\exp\!\big(+i\sigma_3A_\triangle/2\big)=\texttt{expIsig3(area)} .
\]
Measured \(\le6.3\times10^{-16}\) for every face at L=1,2,4,8 and at two different tilts. This
settles open item 5 of `02-formulation.md` §10: the sign is **+**, no extra \((\pm1)_{\rm cut}\)
factor is needed once the \(\mathbb Z_2\) lift is solved for (see below).

**What the fields store.** `Omega_ab = Edge.sgn * expIsig3(Edge.omega)`, and the reverse hop is its
adjoint, `Edge.sgn * expIsig3(-Edge.omega)`.
* `Edge.omega` is the **principal value in \((-\pi,\pi]\)** of the SO(2) angle produced by the
  primary, coordinate-free construction: parallel-transport \(\hat e_\theta(a)\) along the
  great-circle edge in SO(3) and read its angle in the chart frame at \(b\). It is **not** the
  continuously-tracked \(\omega\).
* `Edge.sgn` is \(\pm1\), the \(\mathbb Z_2\) spin-structure lift, from the \(\mathbb F_2\) solve.
* `omegaChart(s, e, nq)` is the independent oracle: it samples \(\alpha(s)\) along the geodesic and
  unwraps, so it returns \(-\int\cos\theta\,d\varphi\) **exactly**, not to \(O(1/n_q^2)\) —
  `nq = 32` and `nq = 64` agree to \(10^{-12}\). It equals `omega` + \(2\pi k\).

**Why `sgn` has to exist at all.** \(\Omega=\exp(i\sigma_3\omega/2)\) has period \(4\pi\) in
\(\omega\), so the \(2\pi\) in \(k\) is a **real sign**, not a no-op. Checking agreement "mod
\(2\pi\)" is blind to exactly the thing that matters.

**Where the sign is forced.** \(\sum_f\sum_i s_i\omega_i=0\) identically (each edge occurs twice
with opposite orientation) while \(\sum_fA_\triangle=4\pi\), so the chart lift's total defect is
\(-4\pi\), and it is carried by exactly \(\chi(S^2)=2\) faces — the two containing \(\pm\hat z'\),
each short by \(2\pi\). Measured, at L = 1, 2, 4, 8: `defect/2pi = -2.000000` with `nf = 2` faces,
per-face residual \(|d+2\pi|<10^{-12}\), and the two faces are verified to be exactly the ones
whose spherical triangle contains a chart pole. **The raw chart lift
\(\exp(i\sigma_3\,\omega_{\rm chart}/2)\) is therefore not a spin structure.**

**The \(\mathbb F_2\) solve.** \(20L^2\) constraints \(x_{e_0}+x_{e_1}+x_{e_2}=k_f\pmod 2\) with
\(2\pi k_f=\sum_i s_i\omega_i-A_\triangle\), in \(30L^2\) unknowns; Gaussian elimination with
free variables set to 0. Rank \(N_F-1\) (the all-ones left null vector is \(\ker\partial_2\)),
consistent because \(\sum_fk_f=0\), solution unique up to the vertex gauge \(\sigma_y\) of
dimension \(N_V-1=30L^2-(20L^2-1)\). **The solve reproduces the paper's antiperiodic \(\varphi\)
cut automatically** — no meridian has to be constructed, and no site/link has to be tested against
it. `tgeom` shows this constructively: cancel the two pole faces with a dual path between them
(that path *is* a discrete meridian) and the \(\mathbb F_2\) lift and the chart lift then differ by
a pure vertex coboundary \(\tau_a\tau_b\) — verified by spanning tree, 0 violations on every edge
at every L.

**Gauge (`tilt`).** `defaultTilt = 0.7`, chosen by scanning `poleGap` over L = 1..16; it leaves the
poles \(\approx\bar a_s/2\) from the nearest site and \(\approx\bar a_s/4\) from the nearest link
geodesic at every level, which is about the best a single axis can do. All metric data
(`abar`, `area`, `Edge.len/dl/area`, `Face.area/sub`) is bit-identical between tilts 0.7 and 1.1;
individual \(\omega\) differ by \(O(1)\); **face holonomies agree to \(<10^{-12}\)** because every
matrix involved is diagonal, so the \(\tau\) conjugation cancels exactly.

**Parity.** \(P\) is the antipodal map in any polar chart and the icosahedron is centrally
symmetric, so \(\omega(Py_1,Py_2)=-\omega(y_1,y_2)\) — (III.3)/(IV.16) — holds to \(<10^{-13}\)
on every edge. Pinned as a test.

### Interface changes (all mirrored in `04-interfaces.md` §3)

1. `Sphere` gained three fields: `nbf` (incident faces in the same cyclic order as `nbr`, so
   `nbf[y][k]` is the face \((y,\mathrm{nbr}[y][k],\mathrm{nbr}[y][k{+}1])\)), `tilt`, and
   `chart` (the three chart axes; `chart[2]` is the polar axis).
2. `defaultTilt* = 0.7`, a `float`.
3. Six spherical primitives exported because the tests and WP-G/WP-I need them:
   `sphArea`, `circum`, `dualLen`, `ptrans`, `chartFrame`, `pframe`; plus `poleGap` as a
   gauge-quality diagnostic.
4. `geom.nim` re-exports `types.nim`; `lattice.nim` re-exports `geom.nim`.
5. Two conventions written into the field comments, both easy to get wrong:
   * `Edge.f[0]` is the face that traverses the edge \(a\to b\), `Edge.f[1]` the one that
     traverses \(b\to a\); `Edge.dl` is ordered to match.
   * **`Edge.eb` is the tangent at \(b\) continuing *away* from \(a\).** The reverse hop needs
     \(e^a_{ba}(b)=-\texttt{eb}\), not `eb`. WP-E: this is the trap in (IV.1).

`core/lattice.nim` matches §4 of `04-interfaces.md` with no changes.

### What did not work / caveats

1. **T1.1c as literally written cannot hold at \(10^{-12}\), and this is not a bug.**
   \(\sum_i\tfrac12\ell_i\ell^*_i=A_\triangle\) is a *flat* identity. Each of the six right
   spherical sub-triangles has legs \(\ell_i/2\) and \(\ell^*_i\) and area
   \(2\arctan(\tan(\ell_i/4)\tan(\ell^*_i/2))\), so the exact statement is
   \[
   A_\triangle=\sum_{i}4\arctan\!\big(\tan(\ell_i/4)\tan(\ell^*_i/2)\big),
   \]
   which the code satisfies to \(1.1\times10^{-16}\). The flat form is off by
   \(3.56\times10^{-2}\) (relative) at L=1 falling to \(9.1\times10^{-4}\) at L=8 — clean
   \(O(\bar a_s^2)\). `tgeom` asserts the exact form at \(10^{-12}\) and asserts the flat form
   *only* shrinks like \(\bar a_s^2\). No tolerance was loosened. Whoever owns
   `03-targets.md`/`02-formulation.md` should restate T1.1c.
2. **`Edge.area` uses the flat diamond formula by definition**, \(A_e=\tfrac12\ell(\ell^*_1+\ell^*_2)\),
   because that is exactly what \(\kappa\) (IV.2) and \(\beta_\ell\) (IV.26) are built from. The true
   spherical area of the diamond differs at the same \(O(\bar a_s^2)\). This makes
   "\(A_e\) consistent with the dual lengths" a tautology rather than a test; the real content is
   carried by item 1 and by the closure relation.
3. **\(\kappa'\) does not converge to \(\tfrac{\sqrt3}{2}(\bar a_s/a_t)\).** \(\kappa\) does:
   \(\langle\kappa\rangle-1/\sqrt3=8.18\!\cdot\!10^{-2},2.06\!\cdot\!10^{-2},5.20\!\cdot\!10^{-3},1.32\!\cdot\!10^{-3}\),
   ratios 3.98/3.95/3.95, clean \(\bar a_s^2\). But
   \(\langle\kappa'\rangle/(\bar a_s/a_t)=\langle A_y\rangle/\bar a_s^2\) runs
   0.854312, 0.856770, 0.864919, 0.867313, 0.867935 (L=16) and converges like \(\bar a_s^2\) to
   \(\approx0.8681\), i.e. \(+0.24\,\%\) above \(\sqrt3/2=0.866025\), not to it. Reason: radial
   projection is not conformal, so the refined mesh is **not** asymptotically uniform — the
   edge-length variance \(\mathrm{Var}(\ell)/\bar a_s^2\) tends to \(5.2\times10^{-3}\), not 0, and
   the local anisotropy \(\max_y|A_y/(\tfrac{\sqrt3}{2}\bar a_y^2)-1|\) tends to \(4.8\,\%\), not 0.
   Both are visible in the `rgeom` table. Nothing is wrong with \(\kappa'\); the flat formula
   simply is not the \(L\to\infty\) limit of a site average on this mesh. The flat relation is
   instead verified exactly, on a genuinely uniform hexagonal patch of radius \(d\): all five
   flat relations hold with errors \(7.2\!\cdot\!10^{-4}\to7.2\!\cdot\!10^{-6}\to7.2\!\cdot\!10^{-8}\)
   at \(d=10^{-1},10^{-2},10^{-3}\), i.e. exactly \(O(d^2)\).
4. **`check r == 0.0` is not usable in this build.** `-Ofast -ffast-math` lets clang contract
   `0.5*len*dual` into an FMA in one place and not the other, so recomputing an identical
   expression differs by \(\sim2\times10^{-17}\). Tests that mean "bitwise identical" have to be
   written as `< 1e-16`. (Comparisons of *stored* values across two `Sphere`s — the tilt test —
   are still exact, and are asserted with `==`.)
5. **Not verified**: nothing in WP-A checks that obtuse triangles are handled, because this mesh
   has none — `min l*` is 0.365, 0.126, 0.050, 0.023 at L=1,2,4,8, always positive. The signed
   \(\ell^*\) and signed `sphArea` code paths are written to cope with \(\ell^*<0\) but are
   untested. If a future mesh (e.g. a different refinement scheme) produces obtuse triangles,
   re-check \(A_y>0\) and the \(\tilde A_i\) partition first.
6. \(\#\) distinct edge lengths is 1, 2, 6, 20 for L = 1, 2, 4, 8. Reported, not asserted — I did
   not derive the edge-orbit count formula. \(\#\) distinct \(A_y\) *is* asserted, against
   \(\mathrm{round}((L{+}3)^2/12)\) (partitions of \(L\) into at most three parts), and matches
   1, 2, 4, 10 exactly.

---

## WP-E — Wilson-Dirac operator  (done)

Files: [`ops/wilson.nim`](../ops/wilson.nim), [`tests/twilson.nim`](../tests/twilson.nim).
Run: `cd build_mac && SDKROOT=$(xcrun --show-sdk-path) make run experimental/radial/tests/twilson`.
**25 tests in 9 suites, all pass**, ~80 s wall (dominated by the L = 4 dense eigensolves).
Dependencies: `std/[algorithm, math, complex, os, strformat, unittest]`, `base/alignedMem`,
`eigens/linalgFuncs`, `core/lattice`, `core/spinor`, and — in the test only —
`core/analytic` (`flatSpectrum`) and `meas/dataio` (`writeTsv`).
Figure data: `output/radial/wilson/`.

### Implemented

* `applyDw` is Eq. (IV.1) written out, nothing else: the kappa sum over spatial neighbours with
  \(-\tfrac12(1-e^a_{y_1y_2}(y_1)\sigma_a)U\Omega\) and \(+\tfrac12\), plus the kappa' temporal
  terms with \(-\tfrac12(1\mp\sigma_3)U\) and \(+1\), minus an optional \(m\).
* `applyDwAdj` assembles the conjugate transpose **block by block from the reverse hop's own
  data** — the tangent at the *far* end, \(e^a_{y_2y_1}(y_2)\), and \(\Omega^\dagger\) — rather
  than using the tetrad identity \((1-e(y_1)\!\cdot\!\sigma)\Omega\to(1+e(y_1)\!\cdot\!\sigma)\Omega\).
  That is deliberate: the adjoint is then exact by construction and T1.2b tests the code, while
  the tetrad hypothesis is tested separately and independently by T1.2a.
* `applyDwDeriv` (tangent in `du`) and `dwPullback` (adjoint mode, all links at once). Every
  link-carrying term is \(z(\theta)=c\,e^{is\theta}\) with \(s=\pm1\), so
  \(d(2{\rm Re}\,z)/d\theta=-2s\,{\rm Im}\,z\); the pullback is that identity applied hop by hop.
* `applyDwHat` / `applyDwHatAdj`: \(\hat D_W=\bar a_s^2A^{-1/2}D_WA^{-1/2}\), via a caller-owned
  `work: var Spin`, so the whole path stays allocation free. `hatScale(l)` exposes
  \(\bar a_s/\sqrt{A_y}\).
* `gaugeTransform`, `spinGaugeTransform`, `denseDw`, `newGauge`, `zero(Gauge)`,
  `nlinkS`/`nlinkT`.
* A `parts: DwParts` selector on every entry point (`dwAll`, `dwC`, `dwB`, `dwSpace`, `dwTime`).
  This is what makes T1.2a and T1.2e expressible at all — the naive/Wilson split of (IV.4) and
  the temporal-only operator are not otherwise reachable — and it costs one `set` test per call.

### Measured — every acceptance criterion

Two fixtures: L = 1 with \(n_t=6,\ a_t=0.4\) and L = 2 with \(n_t=4,\ a_t=0.3\), both on a
random non-compact gauge field (\(\theta\sim0.7\,N(0,1)\) per link).

| # | criterion | asked | measured (L=1 / L=2) |
|---|---|---|---|
| T1.2a | \(\|C+C^\dagger\|/\|C\|\) | 1e-13 | **2.10e-16 / 2.62e-16** |
| T1.2a | \(\|B-B^\dagger\|/\|B\|\) | 1e-13 | **0 / 0** (exactly) |
| T1.2a | \(\max\|D-B-C\|\) | 1e-14 rel | **0 / 0** (exactly) |
| T1.2a | \(\max\|{\rm Re}\,{\rm eig}\,C\|\) / \(\max\|{\rm Im}\,{\rm eig}\,B\|\) | 1e-12 | **2.55e-15 / 4.88e-15** |
| T1.2b | \(\|\langle x,Dy\rangle-\langle D^\dagger x,y\rangle\|\) rel, m = 0 and 1 | 1e-13 | **3.3e-17, 4.2e-17 / 1.5e-17, 5.7e-18** |
| T1.2b | `applyDwAdj` vs \(A^\dagger\) of `denseDw` | 1e-13 | **1.82e-16 / 1.61e-16** |
| T1.2c | gauge covariance, amplitude | — | **2.15e-16 / 1.93e-16** |
| T1.2c | the same as a *squared* relative residual | 1e-20 | **4.64e-32 / 3.73e-32** |
| T1.2c | seam-only gauge (support on one \({\rm sgn}=-1\) edge) | 1e-14 | **1.56e-17** |
| T1.2d | `applyDw` vs `denseDw`, m ∈ {0, 0.7} × parts ∈ {all, C, B, space, time} | 1e-12 | **≤ 1.95e-16** |
| T1.2e | \(D_T\psi=\kappa'_y(1-e^{\mp ik})\psi\), all \(n_t\) half-integer modes × both spins | 1e-12 | **6.66e-15 / 2.92e-15** |
| T1.2e | integer-Matsubara mode is *not* an eigenvector | > 0.1 κ' | **4.73** (κ' = 4.4) |
| ladder 4 | tangent vs centered FD, best of 5 steps | 1e-10 | **3.80e-11 / 4.40e-11** at ε = 3e-6 |
| ladder 5 | \(\langle f,\delta u\rangle=2{\rm Re}\langle x,\delta D\,y\rangle\) | 1e-12 | **1.15e-15 / 1.60e-15** |
| ladder 5 | pullback vs per-link FD (every 7th spatial, every 5th temporal link) | 1e-8 rel | **9.84e-10** on a scale 1.41e+1 |
| ladder 6 | Ward \(\delta_{d\alpha}D=i(\alpha D-D\alpha)\) | 1e-13 | **2.21e-16 / 2.21e-16** |
| — | `applyDwHat` = S D S, and = the dense S A S | 1e-13 | **0 / 0** and **2.04e-16 / 1.99e-16** |
| — | \(\langle x,\hat Dy\rangle=\langle\hat D^\dagger x,y\rangle\) | 1e-13 | **3.87e-17 / 2.30e-17** |
| — | gauge covariance of \(\hat D_W\) | 1e-14 | **2.48e-16** |
| — | allocation regression, 64 × (Dw, Dw†, δDw, D̂w, D̂w†, pullback) | unchanged | `getRawMemAllocated` **0 → 0**, `getOccupiedMem` **18133920 → 18133920** |

**On T1.2c's 1e-20.** A relative *amplitude* of 1e-20 is below the double-precision floor: the
check compares \(e^{i(\theta+\alpha_b-\alpha_a)}\) against
\(e^{i\alpha_b}e^{i\theta}e^{-i\alpha_a}\), which differ by a few ulp no matter how the code is
written, so ~2e-16 is the best attainable. The test asserts **both** readings — the squared
relative residual against 1e-20 (met with 12 orders of margin) and the amplitude against 1e-14 —
and echoes the number, rather than quietly relaxing the target.

**On `getRawMemAllocated`.** It counts QEX `alignedMem` only, and this subtree deliberately uses
plain Nim `seq`, so it reads 0 → 0 and proves nothing by itself. `getOccupiedMem()` is the check
that has teeth here; both are asserted, as in `tspinor.nim`.

### T1.2h — free spatial spectrum, and a correction to doc/02 §3.2

**The volume weight in doc/02 §3.2 is inverted.** (IV.11) says
\(\bar a_sa_t\,D_{\rm lat}=\delta V\,D_{\rm cont}\), hence
\(D_{\rm cont}=\bar a_s\,{\rm diag}(1/A_y)\,D_{\rm lat}\) and the generalized eigenproblem is
\(D_{\rm lat}\psi=\tilde\lambda\,\delta V\,\overline{\delta V}^{-1}\psi\), i.e. the weight is
\(\overline{\delta V}/\delta V\), **not** \(\delta V/\overline{\delta V}\) as §3.2 writes.
Two independent confirmations:

1. *Analytic.* On the flat equilateral lattice \(\kappa=1/\sqrt3\), the naive part is
   \([C\psi]_y\simeq\tfrac{\kappa a}{2}\sum_{k=1}^{6}e_k^ae_k^b\sigma_a\partial_b\psi
   =\tfrac{\sqrt3a}{2}\sigma\!\cdot\!\partial\psi\), and \(A_y/\bar a_s=\sqrt3a/2\) exactly.
   So \(D_{\rm lat}=(A_y/\bar a_s)D_{\rm cont}\).
2. *Numerical.* \(\max_{l=0}\big||{\rm Im}\lambda|-1\big|\) for the two weights:

   | L | \({\rm eig}(\hat D_{\rm spatial})/\bar a_s\) | \({\rm eig}({\rm diag}(volw/volbar)D_{\rm lat})\) |
   |---|---|---|
   | 1 | **0.08413** | 0.13372 |
   | 2 | **0.02428** | 0.50760 |
   | 4 | **0.00627** | 0.74542 |

   The doc's weight moves *away* from the continuum as the lattice is refined; the correct one
   converges as \(O(\bar a_s^2)\) (ratios 3.46, 3.87 against \(\bar a_s^2\) ratios 3.51, 3.89).
   Both columns are printed by the test, which asserts only the inequality, not either value.

The clean statement of the normalization, and the one WP-F should use:
\(\boxed{\ {\rm eig}(\hat D_W)=\bar a_s\,{\rm eig}(D_{\rm cont})\ }\) with
\(\hat D_W=\bar a_s^2A^{-1/2}D_WA^{-1/2}\) — because
\(\hat D_W\sim{\rm diag}(\bar a_s^2/A_y)D_{\rm lat}=\bar a_s\cdot{\rm diag}(\bar a_s/A_y)D_{\rm lat}\).

**Tower, multiplicities and convergence** (U = 0, temporal term off, \({\rm eig}(\hat D)/\bar a_s\)
sorted by \(|\lambda|\)):

| l | multiplicity | L=1 \(\max|\lambda\mp i(l{+}1)|\) | L=2 | L=4 | L=1 \(\max\||{\rm Im}|-(l{+}1)\|\) | L=2 | L=4 |
|---|---|---|---|---|---|---|---|
| 0 | 4 (2 + 2) | 0.27344 | 0.14679 | **0.07506** | 0.08413 | 0.02428 | **0.00627** |
| 1 | 8 (4 + 4) | 1.75430 | 1.00446 | **0.52196** | 0.92768 | 0.29550 | **0.07933** |
| 2 | 12 (6 + 6) | 3.53460 | 2.31403 | **1.25021** | 2.59041 | 0.98822 | **0.28078** |

Every entry decreases with L, which is what the test asserts. The **multiplicities are exactly
the continuum ones**, \(2(l+1)\) per sign, and at L = 1 they are *complete*: \(4+8+12=24=2N_V\),
so the whole 12-site spectrum is the \(l=0,1,2\) towers and nothing else. The groups are
separated by a gap ≥ 0.60 in \(|\lambda|\) (asserted), and each splits exactly evenly between
\({\rm Im}>0\) and \({\rm Im}<0\) (asserted).

The full \(|\lambda\mp i(l+1)|\) error is dominated by the Wilson term's real part, which is
\(O(\bar a_s)\) — at L = 4, \(l=0\): \({\rm Re}=0.0748\), \(|{\rm Im}|-1=-0.0063\). So the
*imaginary* part converges as \(O(\bar a_s^2)\) and the *modulus* only as \(O(\bar a_s)\); both
are reported because quoting only the second would understate the operator and only the first
would overstate it.

### T1.2f / T1.2g — Fig. 4 and Fig. 5 data, and the doubler comparison

**Matsubara block decomposition.** With U = 0 the operator is t-translation invariant, so it
block diagonalizes on the antiperiodic frequencies \(k_n=(2n+1)\pi/n_t\):
\(D(k)=D_{\rm spatial}+\kappa'_y(1-\cos k+i\sigma_3\sin k)\). That turns the L = 4, \(L_t=24\)
panel (a 7776 × 7776 `zgeev`, hours) into 24 problems of size 324 (seconds). Verified against a
full dense diagonalization at L = 1, \(n_t=6\): **max \(|\lambda_{\rm dense}-\lambda_{\rm block}|
=5.78\times10^{-14}\)** over all 144 modes (greedy nearest-neighbour matching — the eigenvalues
are 8- and 12-fold degenerate to the last bit, so a lexicographic sort does *not* align the two
lists; the first version of that test compared sorted lists and reported a spurious 5.1).

Files, all with `#key=value` headers via `meas/dataio.writeTsv`:

* `output/radial/wilson/fig4_L{1,2,4}_Lt24.tsv`, `fig4_L2_Lt{16,48}.tsv` — columns
  `re im matsubara`, **raw \(D_{\rm lat}\) of (IV.1)**, T = 4 so \(a_t=4/L_t\).
* `output/radial/wilson/fig5_flat_L{...}_Lt{...}.tsv` — `flatSpectrum` (IV.8) with
  \(\kappa=1/\sqrt3\) and the matching \(\kappa'=(\sqrt3/2)\bar a_s/a_t\), 16 × 16 × \(L_t\) grid.
* `output/radial/wilson/doublers.tsv` — the table below.

Fig. 4 is plotted in **raw** \(D_{\rm lat}\) units because that is the normalization (IV.8) is
written in; the volume-rescaled version is one diagonal similarity away and is what T1.2h uses.

**Doubler definition.** In the flat formula the doubler is where \(\lambda\) is real and nonzero
(\(4\kappa\) at \(k=(\pi,0)\) etc.). On the sphere no eigenvalue is exactly real, and with
half-integer Matsubara frequencies the 3D spectrum never reaches \(k_t=0\) either, so the
quantity measured is the \(k_t\to0\) (spatial) operator's
\(\min\{{\rm Re}\lambda:\ |{\rm Im}\lambda|\le0.4\,\bar a_s\}\). The cut scales with \(\bar a_s\)
because the curvature-induced spread of the doubler cluster does: the selected eigenvalues have
\(|{\rm Im}| = 0.387,\ 0.165,\ 0.084\) at L = 1, 2, 4, i.e. \(\approx0.28\,\bar a_s\) throughout.

| L | \(\bar a_s\) | spatial doubler | \(|{\rm Im}|\) | flat \(4\kappa=4/\sqrt3\) | rel. error |
|---|---|---|---|---|---|
| 1 | 1.10715 | 2.274610 | 0.3874 | 2.309401 | **1.51 %** |
| 2 | 0.59095 | 2.300120 | 0.1649 | 2.309401 | **0.40 %** |
| 4 | 0.29947 | 2.306054 | 0.0841 | 2.309401 | **0.14 %** |

The temporal doubler is *exact* rather than fitted: the \(k_t=\pi\) mode of the temporal operator
has eigenvalue \(2\kappa'_y\) site by site (this is the T1.2e identity at \(k=\pi\)), so the
comparison is \(2\bar\kappa'=2\bar A/(\bar a_sa_t)\) against the flat
\(2\kappa'=\sqrt3\bar a_s/a_t\), a ratio \(2\bar A/(\sqrt3\bar a_s^2)\) independent of \(a_t\):

| L | \(2\kappa'_y\) range at \(L_t=24\) | mean | flat \(2\kappa'\) | rel. error of the mean |
|---|---|---|---|---|
| 1 | 11.35021 (all equal) | 11.35021 | 11.50583 | **1.35 %** |
| 2 | [5.56079, 6.28161] | 6.07566 | 6.14130 | **1.07 %** |
| 4 | [2.33723, 3.36498] | 3.10825 | 3.11223 | **0.13 %** |

Note the *spread*: at L = 4 the individual \(2\kappa'_y\) run from 2.34 to 3.36, a ±17 % band
around the flat value even though the mean agrees to 0.13 %. The flat formula describes the
lattice average, not the site-by-site operator — which is exactly the paper's "qualitative
description, first doubler position agrees to good approximation".

Top of the spectrum, \(\max{\rm Re}\lambda\), against the flat spectrum. The "flat grid"
column is the maximum over the same 16 × 16 × \(L_t\) sample that goes into the Fig. 5 files;
it sits ~0.2 % below the exact \(4.5\kappa+2\kappa'\) because 16 does not divide the K point
\(k=(2\pi/3,2\pi/3)\) where the \(4.5\kappa\) maximum lives (see the WP-C section).

| L, \(L_t\) | curved | flat grid | rel | exact \(4.5\kappa+2\kappa'\) |
|---|---|---|---|---|
| 1, 24 | 13.57627 | 14.08801 | 3.63 % | 14.10390 |
| 2, 24 | 8.59137 | 8.72348 | 1.51 % | 8.73937 |
| 4, 24 | 5.79019 | 5.69441 | 1.68 % | 5.71031 |
| 2, 16 | 6.51197 | 6.67638 | 2.46 % | 6.69227 |
| 2, 48 | 14.83263 | 14.86478 | 0.22 % | 14.88067 |

Part of this gap is a known systematic, not a lattice artifact: `flatSpectrum` uses *integer*
Matsubara \(k_t=2\pi i/n_t\) and therefore reaches \(k_t=\pi\) exactly for even \(n_t\), while the
fermion is antiperiodic and its closest frequency is \(\pi(1\pm1/L_t)\), costing
\(\kappa'(1-\cos(\pi/L_t))\approx\kappa'\pi^2/(2L_t^2)\) — 0.049 at \(L_t=24\), L = 1. The trend
in the last column (0.22 % at \(L_t=48\), 2.5 % at \(L_t=16\)) is exactly that \(1/L_t^2\).

### The T2.1 free column — which operator carries M (a WP-F/WP-K finding)

Not a WP-E acceptance criterion, but the cheapest possible test of the normalization, so it is
measured and reported: doc/03 T2.1 quotes the deterministic free-limit \(\min|D_W-1|\) as
**1.154, 1.010, 0.965** for L = 1, 2, 4 at \(a_t=0.2\), \(M=1\). Computed here:

| L | \(\bar a_s/a_t\) | raw, \(L_t=60\) | raw, \(L_t=80\) | hat, \(L_t=60\) | hat, \(L_t=80\) | published |
|---|---|---|---|---|---|---|
| 1 | 5.536 | 1.1234 | 1.1104 | 1.2682 | 1.2555 | 1.154 |
| 2 | 2.955 | **1.0105** | 1.0074 | 1.0887 | 1.0853 | **1.010** |
| 4 | 1.497 | **0.9643** | 0.9635 | 1.0229 | 1.0225 | **0.965** |

So **the paper's \(D_W\) in \(\min|D_W-1|\) is the raw \(D_{\rm lat}\) of (IV.1), at
\(a_t=0.2\) with \(L_t=60\), i.e. \(T=12\)** — L = 2 and L = 4 reproduce to 4 and 3 digits, and
the volume-normalized operator is 8 % and 6 % off. That is consistent with (IV.10),
\(M_0=\min(4/\sqrt3,\sqrt3\bar a_s/a_t)\), which is literally the raw flat-limit
\(\min(4\kappa,2\kappa')\). **WP-F must therefore convert M when it uses \(\hat D_W\):** the two
differ by the site-dependent factor \(\bar a_s^2/A_y\), which is 1.17055 (uniform) at L = 1, and
spans [1.0680, 1.5376] at L = 4 around the flat \(2/\sqrt3=1.1547\) — it is **not** a scalar, so
"rescale M by 1.1547" is only approximate and the allowed window (IV.10) has to be rechecked in
whichever normalization is used.

The L = 1 cell does **not** reproduce (1.1234 vs 1.154, 2.7 % low) under any \(L_t\): the scan is
monotone in \(L_t\) and never passes through 1.154 at a round T —
\(L_t=12,16,20,24,32,48,60,80,120,168\) give
1.6125, 1.4287, 1.3348, 1.2693, 1.1956, 1.1400, 1.1234, 1.1104, 1.1010, 1.0972.
1.154 would need \(L_t\approx42\) (T ≈ 8.4). Reported as an unexplained discrepancy, not
worked around; WP-K owns T2.1. (\(1.154\approx2/\sqrt3=1.1547\) is suggestive of a normalization
slip in the published L = 1 entry, but that is a guess and is not asserted anywhere.)

### Interface changes (mirrored in [`04-interfaces.md`](04-interfaces.md) §7)

1. **`parts: DwParts = dwAll` added to `applyDw`, `applyDwAdj`, `applyDwDeriv`, `dwPullback`,
   `denseDw`.** New enum `DwPart` = {`dwSpatial`, `dwTemporal`, `dwNaive`, `dwWilson`} with the
   sets `dwAll`, `dwC`, `dwB`, `dwSpace`, `dwTime`. Needed by T1.2a (the C/B split of (IV.4)),
   T1.2e (temporal only) and T1.2h (spatial only); the alternative was three more exported procs.
2. **`applyDwHat`/`applyDwHatAdj` take `work: var Spin`** — a caller-owned scratch field. Nim
   `var` parameters cannot carry defaults, so it sits before `m`:
   `applyDwHat(l, dst, src, u, work, m = 0.0, parts = dwAll)`. Without it the routine would have
   to allocate a temporary per apply, which is exactly the performance failure the brief warns
   about. `Ov.work` in §10 already exists to hold it.
3. **`hatScale(l): seq[float]`** is new and exported (\(\bar a_s/\sqrt{A_y}\) per sphere site),
   because the tests and the dense spectra need the same factors the kernel uses.
4. **`zero(u: var Gauge)`, `nlinkS(l)`, `nlinkT(l)`** are new: force accumulation and the link
   layout contract (`ne*nt`, `nv*nt`) in one place. WP-G/WP-H will want all three.
5. The **link and boundary conventions are now written down** in §7 (they were not): the sign of
   \(\theta\) per hop direction, the temporal link direction, and that the −1 seam factor applies
   to both directions of the \(t=n_t-1\to0\) link.
6. **Ward-identity sign.** §15 item 6 says \(\delta_{d\alpha}D=i(D\alpha-\alpha D)\). With
   doc/02 §5's link convention \(\theta_e\to\theta_e+\alpha_b-\alpha_a\) (the one that makes
   (IV.24) invariant) the correct statement is \(i(\alpha D-D\alpha)\); the test asserts that and
   §7 records it. The first version of the test used the doc's sign and failed at exactly 2.000,
   which is how it was caught.
7. `denseDw` is documented as an **independent** assembly (Mat2 blocks written directly), not a
   basis sweep of `applyDw` — otherwise T1.2d would be vacuous.

### What did not work, and things to be careful about

1. **The FD step range matters.** The first version scanned ε ∈ [1e-2, 1e-4] and bottomed out at
   5e-9 — pure \(O(\varepsilon^2)\) truncation, not a bug. The optimum for a centered difference
   is \(\varepsilon\sim\epsilon_{\rm mach}^{1/3}\approx5\times10^{-6}\); with ε down to 1e-6 the
   best is 3.8e-11. If someone tightens this test further, move the *steps*, not the tolerance.
2. **`-d:danger` means no bounds checks.** The T1.2h group loop indexed `ev[k+mult]` past the end
   at L = 1, l = 2 (24 eigenvalues, groups 4 + 8 + 12 exhaust them) and silently returned garbage
   (`gap = -2.44`) instead of crashing. Guarded now; the printed `gap = nan` for that one cell is
   deliberate.
3. **Do not compare eigenvalue lists by sorting.** See the Matsubara check above.
4. **`dst` must not alias `src`** in any apply — the hop terms read neighbours after the row has
   been written. Documented in the module header; not defended at run time (it would cost a
   branch per call and the repo's rule is to let mistakes crash).
5. **`applyDw` asserts `dst.len == l.nsite` rather than resizing.** Resizing would allocate, and
   an allocation inside the CG inner loop is precisely the failure mode this work package was
   told to avoid.
6. **`denseDw` accumulates with `+=`.** That is load bearing at \(n_t\le2\), where \(t+1\) and
   \(t-1\) land on the same slice (and at \(n_t=1\) on the row itself); writing instead of
   accumulating would silently drop one of the two temporal hops. `nt = 1` is used throughout the
   spectral tests to get the spatial operator, so this path is exercised.
7. **The doubler cut is a choice, and it is stated as one.** \(0.4\,\bar a_s\) is not derived; it
   is the value for which the selected cluster is stable (the same eigenvalue is picked for cuts
   0.3–0.5 \(\bar a_s\) at all three L). An absolute cut is *not* usable: 0.4 picks 2.2318 at
   L = 4 instead of 2.3061. The raw near-real spectrum is in the Fig. 4 files, so anyone who
   prefers a different definition can recompute without rerunning the eigensolves.
8. **No external oracle exists for the curved-space operator.** What stands in for one: the
   analytic S² Dirac tower \(\pm i(l+1)\) with multiplicity \(2(l+1)\) (reproduced including the
   exact multiplicities), the flat-limit formula (IV.8) via `core/analytic`, the exact
   half-integer Matsubara eigenvalue, the published T2.1 free column, and the structural
   (anti)hermiticity of (IV.4). Every number above is measured against one of those, or against
   a second independent coding in the same test — never against the operator's own output.

---

## WP-G — gauge action, zero modes, flow  (done)

Files: [`ops/gaugeact.nim`](../ops/gaugeact.nim), [`ops/flow.nim`](../ops/flow.nim),
[`tests/tgauge.nim`](../tests/tgauge.nim), [`tests/tflow.nim`](../tests/tflow.nim),
[`rgauge.nim`](../rgauge.nim).
Run: `cd build_mac && SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk make run
experimental/radial/tests/tgauge` (28 tests, 6 suites) and `… tflow` (8 tests). All pass.
Dependencies: `core/lattice`, `ops/wilson` (for `Gauge`), `ops/solve` (for `CgInfo` only),
`rng/threefry4x64`, `algorithms/rk`. `rgauge.nim` adds `base`, `core/analytic`, `meas/fit`,
`meas/dataio`. Outputs under `output/radial/gauge/`.

### Implemented

* **(IV.24) with the free-limit couplings (IV.26)**, literally: `gaugeAction`,
  `gaugeActionParts` (magnetic/electric split), `gaugeForce`. `Theta_tri` uses `faces[f].e`
  and `faces[f].s`; `Theta_e(t) = theta_e(t) + theta^t_b(t) - theta_e(t+1) - theta^t_a(t)`
  for the canonical edge orientation `a -> b`, periodic in t.
* **`Beta`** — the couplings precomputed once per (lattice, g², convention).
  `betaFace`/`betaEdge` do two divisions each and `gaugeForce` is the CG inner loop, so
  hoisting them out is worth a factor. The g²-taking signatures of
  [`04-interfaces.md`](04-interfaces.md) §11 are kept as one-line wrappers.
* **`GeomConv`** — the O(a²)-convention knob, `newBeta(l, g2, conv)`:
  `gcGeodesic` (doc/02 normative: geodesic ℓ, spherical excess `A_△`, `A_e` = the flat form
  `½ℓ(ℓ*₁+ℓ*₂)` as doc/04 §3 defines `Edge.area`); `gcExactArea` (same but `A_e` = the exact
  spherical diamond area `Σ 4 arctan(tan(ℓ/4) tan(ℓ*/2))` — **this is the paper's convention,
  see T1.5b**); `gcFlat` (chord ℓ, planar `A_△`, planar circumcentre, in-plane signed `ℓ*`).
* **Zero modes**: `gradient` (d), `divergence` (d†), `laplace` (d†d, coded directly and tested
  against the composition), `projectGauge` (CG on d†d with the constant mode removed every
  step, then `p -= d alpha`), `projectFlat`, `projectKernel`.
* **Solvers**: `cgM`; `pseudoSolve` = the literal (V.16)-(V.17) double CG; `RegOp`/`regSolve`
  = the kernel-regularized SPD operator (see "what did not work"); `gaugePropagator`.
* **`heatbath`** — exact, no HMC: draw ξ ~ N(0,1) per plaquette, form `b = Cᵀ W^{1/2} ξ`
  (covariance exactly `M`, and `b ∈ range(M)` by construction), then `theta = M⁺ b`, whose
  covariance is `M⁺ M M⁺ = M⁺`. One CG per configuration, 24 iterations at L=1/nt=4.
* **`jtop`** = `Theta_tri/A_tri`, (V.12), plus a `Beta` overload that uses that object's area
  convention. `triSource` builds the incidence row so `<Theta Theta> = bᵀ M̃⁻¹ b'`.
* **`ops/flow.nim`** on `algorithms/rk.nim`: `newFlowOp`, `flowStep`, `flowRun` (measurement
  callback, exact hits on the requested times), `energyDensity` / `energyDensityT`.
  `advance` is `delta = alpha*delta - beta*force; y += delta` — plain addition, non-compact,
  no exponential and no projection. The force scratch is captured in the closure, so a whole
  flow allocates nothing.

### ker M is one dimension larger than the gauge orbit — and it is not optional

The uniform temporal mode `theta^t_v(t) = c`, `theta^s = 0` costs no action, because
`sum_t Theta_e(t) = W_b - W_a` only sees *differences* of the temporal Wilson lines. It is not
a gauge mode: a gauge function periodic in t can only make temporal shifts that sum to zero
round the time circle. It is orthogonal to the gauge orbit (`<P, d alpha> = 0` for the same
reason), so `projectKernel = projectGauge + projectFlat` and

    dim ker M = n_V L_t  (measured: 48 = 12*4),   rank M = n_E L_t  (measured: 120 = 30*4)

both confirmed against the dense eigendecomposition. The independent Python oracle needed the
same extra direction, from the other side — good agreement on a subtle point.

### Measured residuals — `tgauge` (L=1, n_t=4, a_t=0.35, g²=1.7; 168 links, 200 plaquettes)

| test | measured | tolerance |
|---|---|---|
| S on a pure gauge `theta = d alpha` | 1.85e-30 | 1e-25 |
| gauge invariance `S(theta + d alpha) − S(theta)`, relative | **0.0** (bit identical) | 1e-12 |
| action vs dense incidence oracle, relative | 1.66e-16 | 1e-12 |
| force vs dense oracle, max abs (`|f|_max` = 13.5) | 3.55e-15 | 1e-12 rel |
| force vs centered FD, h = 1e-1/1e-2/1e-3/1e-4 | 4.3e-13 / 3.7e-12 / 2.9e-11 / 3.0e-10 | best < 1e-8 |
| `<f, d alpha>/(|f||d alpha|)` | 1.35e-17 | 1e-14 |
| `sum_faces Theta_tri` per slice | < 1e-12 rel | 1e-12 |
| `dim ker M`, `rank M` | 48, 120 | exact |
| `projectKernel` on every kernel eigenvector | 4.64e-15 | 1e-10 |
| `projectGauge` on a pure gauge, `|p|` 19.25 → | 4.05e-15 | 1e-10 rel |
| `projectGauge`: `\|div p\|²` 4.007e+02 → 2.20e-29 (ratio 5.5e-32), 11 CG its | | 1e-20 rel |
| `projectGauge` idempotency, second pass | 4.44e-16 | 1e-12 rel |
| **T1.5g** double-CG `M̃⁻¹` vs dense pseudo-inverse | 1.54e-14 | 1e-10 |
| **T1.5g** `x(b) = x(b + kernel)` | 2.52e-14 | 1e-10 |
| **T1.5g** `<Theta Theta>` vs dense, all (f,t) | 6.06e-15 | 1e-10 |
| `regSolve` vs `pseudoSolve` | 3.67e-14 | 1e-12 |
| heatbath `<S>` = 59.72 ± 0.31 over 600 samples vs rank/2 = 60 | pull −0.91 | 5σ |
| heatbath `<Theta²>` = 1.2733 ± 0.0740 vs exact 1.290095 | pull −0.23 | 5σ |
| heatbath field: `\|div u\|²`, `sum theta^t` | 1.8e-27, 4.5e-15 | transverse |
| `gaugeForce` allocation over 64 applies | `getOccupiedMem()` unchanged | — |

The FD column is *not* a truncation study: for a quadratic action the centered difference is
exact, so the error is pure cancellation roundoff and the *largest* step is the best one. That
is why the check takes the best of four rather than a slope.

### Measured residuals — `tflow` (same lattice, λ_max(M) = 12.4466)

| test | measured |
|---|---|
| RK vs exact `exp(−M s) theta₀`, s = 0.02 / 0.10 / 0.50 | 1.07e-14 / 4.22e-14 / 2.48e-12 (tol 1e-10) |
| RK3W6 order from halving (n = 16→128) | errs 1.174e-04 … 1.979e-07, orders **3.12, 3.06, 3.03** |
| RK4CK order from halving | errs 1.903e-06 … 4.140e-10, orders **4.09, 4.05, 4.02** |
| linearity `flow(3u) − 3 flow(u)` | 5.33e-15 |
| gauge covariance `flow(u + d a) − flow(u) − d a` | 7.99e-15 |
| S, E_s, E_t monotone along the flow | all three monotone over s ∈ [0, 1.6] |
| `<E_s(0.05)>` over 200 heatbath samples vs the exact trace `tr(M_s M⁺ e^{−2Ms})/2V` | 0.78875 vs 0.772017 (2.2 %, tol 5 %) |

The order test starts at n = 16 deliberately: at n = 2 the step has `h λ_max = 2.5`, outside the
asymptotic regime, and the apparent order comes out 5.1 (RK3) and 4.8 (RK4). The comment in the
test says so; the tolerance was not loosened.

### T1.5a — G_g(t) against the independent Python oracle (L=1, L_t=120, T=16)

`G_g(t) = (1/g²) <J^t(t) J^t(0)>` averaged over all 20 source triangles, from `regSolve` at
`r2req = 1e-26`. **Every digit the oracle quoted is reproduced:**

| t | G_g(t), Nim | oracle | rel. dev | Δ_eff, Nim | oracle |
|---|---|---|---|---|---|
| 0.0000 | 1.96741714e+00 | 1.96741714e+00 | +1.4e-09 | 2.816897 | 2.816897 |
| 0.5333 | 4.62879753e-01 | 4.62879753e-01 | −6.9e-10 | 2.532089 | 2.532089 |
| 1.0667 | 1.27395702e-01 | 1.27395702e-01 | −1.7e-09 | 2.233189 | 2.233189 |
| 2.0000 | 1.90361979e-02 | 1.90361979e-02 | −2.4e-10 | 1.809259 | 1.809259 |
| 3.0667 | 3.18269428e-03 | 3.18269428e-03 | +1.4e-09 | 1.545205 | 1.545205 |
| 4.0000 | 7.89217296e-04 | 7.89217296e-04 | −2.9e-11 | 1.441194 | 1.441194 |
| 5.0667 | 1.74510635e-04 | 1.74510635e-04 | +1.7e-09 | 1.390088 | 1.390088 |
| 6.0000 | 4.82849711e-05 | 4.82849711e-05 | −9.2e-10 | 1.371383 | 1.371383 |

The residual 1e-9 is the oracle's own quoted precision (9 significant digits). Two independent
implementations — a Fourier-space per-momentum inversion in Python and a real-space
pseudo-inverse CG in Nim — agree to the last printed digit at every t. **The gauge sector is
right.** `G_g` is also g²-independent to 2.7e-13 (measured at g² = 0.05 vs 3.7), as it must be:
the 1/g² of (V.14) cancels the g² of `M̃⁻¹`.

Against the continuum (V.14) `gaugeGPeriodic(t, 16)`, the L=1 ratio runs 1.28 at t = 4 to 1.52
at t = 8 — an O(ā²)=O(1.2) lattice at the coarsest level; the ratio at L=4 is 1.10 → 1.13.

Under `gcExactArea` (the paper's convention, T1.5b) the same L=1 correlator is a different
O(ā²) representative of the same continuum function — for the record:

| t | 0 | 0.5333 | 1.0667 | 2.0 | 3.0667 | 4.0 | 5.0667 | 6.0 |
|---|---|---|---|---|---|---|---|---|
| G_g exact-area | 1.93333729 | 4.66030750e-1 | 1.30684453e-1 | 1.99691311e-2 | 3.40566260e-3 | 8.59960736e-4 | 1.94526577e-4 | 5.49935238e-5 |
| Δ_eff | 2.767699 | 2.493116 | 2.204210 | 1.789806 | 1.526499 | 1.420738 | 1.367889 | 1.348228 |

### T1.5b — Δ₀ at L=1, L_t=120 — **the paper's A_ℓ is the exact spherical diamond area**

Headline, and the most consequential single finding of WP-G:

| convention for `A_ℓ` in β_ℓ = 2A_ℓ/(g² ℓ² a_t) | Δ₀(L=1, L_t=120, T=16) | vs published 1.33242 |
|---|---|---|
| flat form `½ℓ(ℓ*₁+ℓ*₂)` with geodesic ℓ, ℓ* — `Edge.area`, doc/04 §3, `gcGeodesic` | 1.356697 | **+1.82 %** |
| **exact spherical diamond `Σ 4 arctan(tan(ℓ/4) tan(ℓ*/2))` — `gcExactArea`** | **1.332430** | **+0.001 %** |
| fully chordal/planar (`gcFlat`) | 1.659221 | +24.5 % |

**1.332430 against 1.33242 is agreement in every published digit** (1.0e-5 relative). That is
not a coincidence at the level of a 6-digit number; the identification is:

> **arXiv:2510.03085's \(A_\ell\) in (IV.26) is the exact spherical area of the diamond,
> not the flat form \(\tfrac12\ell(\ell^*_1+\ell^*_2)\).**

The exact form is the one that tiles the sphere: `Σ_e A_e^exact = Σ_△ A_△ = 4π` to 1e-12 at
every refinement level, while the flat form misses by 3.6 % (L=1), 1.0 % (L=2), 0.25 % (L=4)
— measured, and now a test in `tgauge.nim`. doc/04 §3 already recorded the exact identity
`A_△ = Σ_i 4 arctan(tan(ℓ_i/4) tan(ℓ*_i/2))`; what was not appreciated is that the paper uses
it in β_ℓ. At L=1 all 30 edges are equivalent, so `A_e^exact = 4π/30 = 0.418879` exactly against
the flat form's `0.403959`, a 3.69 % increase of β_ℓ.

**This has to be checked in the fermion sector too, and it is not obviously the same answer.**
`κ_{y₁y₂} = 2A_{y₁y₂}/(ā_s ℓ)` uses the same `A_{y₁y₂}`, and `core/lattice.newLat` builds it
from `Edge.area`, i.e. the flat form. But (IV.2) *writes both forms* and says the second follows
from `A = ½ℓ(ℓ*₁+ℓ*₂)`, which is the flat identity — so for κ the paper appears to mean the flat
form, and only for β_ℓ the exact one. **Highest-value follow-up: rerun T1.4c (published
Δ₀ = 0.953918, fermion, L=1, L_t=168) with `κ` built from the exact spherical diamond area and
see which one lands.** That is a one-line change in `newLat` and it decides the convention for
every single-lattice number in the project. WP-G did not make it, because `core/lattice.nim`
and `ops/wilson.nim` are not WP-G's files.

Nothing was tuned to reach 1.33242. The route was: (i) confirm the geodesic number against the
independent Python oracle to 1e-9; (ii) show it is not the fit window and not T; (iii) refute
the flat hypothesis; (iv) measure `dΔ₀/d ln(β_ℓ/β_△) = −0.684` and read off that the paper
needs β_ℓ 3.6 % higher; (v) notice that the exact spherical diamond area is 3.69 % higher and
is the one form doc/04 §3 flags as exact; (vi) implement it and get 1.332430. Steps (iv)-(vi)
are a prediction and its confirmation, not a fit.

The rest of this subsection is the geodesic-convention evidence trail that led there.

`effMass` + `plateauFit` over 4 ≤ t < 8 (indices 30…59), a_t = 0.13333, T = 16, `gcGeodesic`:

    Delta_0 = 1.356697 +- 0.000095,  c = 2.792, Delta' = 0.8741, chi2/dof = 3.0e-08, dof 27

which is +1.82 % above the published 1.33242 and −4.07 % from √2. This confirms the
coordinator's independent Python fit (1.3567) exactly. It is not a fit artefact:

| window | Δ₀ | χ²/dof | | window | Δ₀ | χ²/dof |
|---|---|---|---|---|---|---|
| [0.5, 8) | 1.336682 | 2.6e-04 | | [5, 8) | 1.357148 | 2.7e-08 |
| [1, 8) | 1.347949 | 4.6e-05 | | [6, 8) | 1.358114 | 1.4e-08 |
| [2, 8) | 1.354881 | 1.4e-06 | | [4, 7) | 1.356238 | 6.8e-11 |
| [3, 8) | 1.356400 | 4.5e-08 | | [4, 6) | 1.356304 | 7.5e-12 |
| [4, 8) | 1.356697 | 3.0e-08 | | [2, 6) | 1.352349 | 1.1e-06 |

Everything with t₀ ≥ 3 sits at 1.3565 ± 0.001, and the windows that drift low are exactly the
ones with a visibly worse χ²/dof, i.e. the ones contaminated by the Δ₂ − Δ₁ = √6 − √2 = 1.0353
excited state (the fitted Δ′ = 0.86…1.07 is that state, as WP-C predicted). It is also not the
value of T: at T = 12 (L_t = 120, a_t = 0.1, window [3, 6)) Δ₀ = 1.360275.

**The flat/chordal hypothesis is refuted.** Under the flat convention Δ₀(L=1) = **1.659221**,
i.e. +17.3 % from √2 — the wrong direction and an order of magnitude too big. Flat is a
*worse* discretization at L=1, unsurprisingly: the icosahedron's planar surface area is 9.5745
against 4π = 12.566, a 24 % deficit, so β_△ jumps by 31 % and β_ℓ by 12 %.

**What the difference actually is, quantified.** At L=1 every face and every edge lies in a
single icosahedral orbit, so `A_△`, `A_e` and `ℓ` are each one number and the entire convention
dependence of the action collapses to the single dimensionless ratio β_ℓ/β_△ (besides a_t) —
see the table below. Scanning a multiplier `escale` on β_ℓ alone:

| escale | 0.80 | 0.90 | 0.95 | 1.00 | 1.05 | 1.10 | 1.20 |
|---|---|---|---|---|---|---|---|
| Δ₀ | 1.516094 | 1.429747 | 1.391779 | **1.356697** | 1.324156 | 1.293866 | 1.239091 |

`dΔ₀/d ln(escale) = −0.684` at the geodesic point. So the published 1.33242 needs β_ℓ/β_△
**3.6 % higher** (`escale = 1.0361`), while the flat-vs-geodesic ambiguity in that same ratio is
a factor 1.4983 and moves Δ₀ by 22 %. That prediction is what `gcExactArea` then satisfied:
`A_e^exact/A_e^flat-form = 0.418879/0.403959 = 1.03694`, predicted Δ₀ = 1.356697 −
0.684·ln(1.03694) = 1.33188, measured **1.332430**. The whole L=1 table, all three conventions:

| | ℓ | ℓ* | A_△ | A_e | β_△ | β_ℓ | β_ℓ/β_△ | Δ₀(L=1) |
|---|---|---|---|---|---|---|---|---|
| flat | 1.051462 | 0.303531 | 0.478727 | 0.319151 | 2.088873 a_t/g² | 0.577350/(a_t g²) | 0.276393/a_t² | 1.659221 |
| geodesic | 1.107149 | 0.364864 | 0.628319 | 0.403959 | 1.591549 a_t/g² | 0.659105/(a_t g²) | 0.414128/a_t² | 1.356697 |
| **exact area** | 1.107149 | — | 0.628319 | **0.418879** | 1.591549 a_t/g² | 0.683425/(a_t g²) | 0.429405/a_t² | **1.332430** |

(`30 A_e` = 9.5745 / 12.1188 / **12.5664 = 4π** down the three rows: only the exact-area
diamonds tile the sphere.)

### T1.5c / T1.5d — the number that has to be right: Δ₀^cont = √2

Full (L, L_t) grid at T = 16, same estimator and window throughout. The **exact-area**
convention (the paper's, per T1.5b) first:

| L | ā_s | L_t = 48 | 64 | 96 | 120 |
|---|---|---|---|---|---|
| 1 | 1.107149 | 1.323527 | 1.328123 | 1.331460 | **1.332430** |
| 2 | 0.590946 | 1.380351 | 1.385565 | 1.389360 | 1.390464 |
| 4 | 0.299474 | 1.396292 | 1.401707 | 1.405649 | 1.406795 |

and the geodesic (flat-form `A_e`) convention:

| L | ā_s | L_t = 48 | 64 | 96 | 120 |
|---|---|---|---|---|---|
| 1 | 1.107149 | 1.347293 | 1.352145 | 1.355672 | 1.356697 |
| 2 | 0.590946 | 1.387145 | 1.392438 | 1.396291 | 1.397412 |
| 4 | 0.299474 | 1.398042 | 1.403478 | 1.407434 | 1.408585 |

(V.7) fit `Δ₀ = a + c_s ā_s² + c_t a_t²` over all 12 points of each grid:

| convention | Δ₀^cont | dev from √2 | c_s | c_t | rms resid |
|---|---|---|---|---|---|
| **exact area (the paper's)** | **1.414535 ± 0.000286** | **+0.023 %** | −0.06509 | −0.10540 | 4.6e-04 |
| geodesic | 1.414695 ± 0.000258 | +0.034 % | −0.04547 | −0.10788 | 4.2e-04 |
| flat | 1.411257 ± 0.001828 | −0.209 % | +0.20254 | −0.14369 | 3.0e-03 |

against the exact √2 = 1.4142136 and the published **1.41409(18)**. All three conventions land
on √2 — that is what "differ at O(a²)" means and it is the check that matters. The fully flat
one approaches from the other side with a 3–4× larger `c_s` and a 7× worse residual, i.e. it is
the worst-behaved family; exact-area and geodesic are equally good and differ only in `c_s`.
**T1.5d is met (+0.023 % on the paper's convention), and with T1.5b at 1e-5 the free gauge
sector is reproduced end to end.**

### T2.2 — gradient-flow scale scan (slide 9)

L ∈ {1,2,4} × g²a ∈ {0.5,1,1.5} (g²R = g²a·L), L_t = 32, a_t = 0.2, 64 exact-heatbath
configurations per point, 30 log-spaced flow times s ∈ [0.02, 3], RK4CK with
`h = 0.1/max_i M_ii`. Step-size convergence `|E(h)/E(h/2) − 1|`: 3.2e-09 (L=1), 2.1e-08 (L=2),
4.3e-07 (L=4). Statistical error on E_s is 0.9–1.4 % at mid flow time. Written to
`output/radial/gauge/flowscan_nt32.tsv`.

    s        r/s     E_s(L=1)  slope   E_s(L=2)  slope   E_s(L=4)  slope
  0.0200   50.000  8.253e-01    --    5.291e+00    --    1.980e+01    --
  0.0474   21.075  6.961e-01  0.269   3.613e+00  0.614   8.363e+00  1.275
  0.0947   10.560  5.496e-01  0.422   2.098e+00  0.971   3.057e+00  1.594
  0.1890    5.290  3.838e-01  0.633   9.425e-01  1.333   9.703e-01  1.689
  0.3773    2.650  2.234e-01  0.957   3.410e-01  1.583   3.070e-01  1.628
  0.7530    1.328  9.991e-02  1.380   1.083e-01  1.718   1.015e-01  1.588
  1.5030    0.665  3.328e-02  1.791   3.229e-02  1.787   3.356e-02  1.644
  3.0000    0.333  8.396e-03    --    8.688e-03    --    9.834e-03    --

Three findings, one of them a correction to the reading of the slide:

1. **For pure gauge the nine curves are not nine curves.** The Gaussian action makes E_s(s)
   *exactly* g²-independent: θ ~ g in distribution and S = W/g², so the action density is
   coupling-free, and with the standard (coupling-independent) flow normalization so is the
   whole trajectory. The three g²a values at fixed L agree within statistics on independent
   ensembles (max |ΔE|/σ = 1.7…3.6 over 180 comparisons). **The g²-dependent splitting on
   slide 9 must therefore come from the dynamical fermions, not the gauge action** — WP-H
   territory, not WP-G's.
2. **The 3/2 reference slope is crossed but there is no plateau at these volumes.**
   `d ln E_s/d ln(r/s)` with r = R = 1 runs monotonically from 0.16 to 2.10 at L=1 (no window
   at all: ā_s = 1.107 ≈ R), and at L=4 rises to 1.69 near s ≈ 0.19 and comes back to ≈ 1.59,
   crossing 3/2 at s ≈ 0.08. The free-Maxwell window needs ā_s ≪ √(8s) ≪ R, i.e.
   s ∈ [0.011, 0.125] at L=4 — one decade, and its lower end is cut by a_t = 0.2. The three
   qualitative regions of the slide are all there (flat/UV-contaminated at large r/s, a
   3/2-crossing middle, steeper-than-3/2 exponential fall at small r/s), the middle one is just
   too narrow at L ≤ 4 to call a plateau.
3. **`E_s √L` does not collapse; raw `E_s` does.** At s ≳ 0.9 the raw action density agrees
   across L to ~2 % (0.0780 / 0.0803 / 0.0772 at L = 1/2/4, s = 0.895; 0.0597 / 0.0594 / 0.0586
   at s = 1.064), which is the expected approach to the continuum value; multiplying by √L
   breaks that. Both columns are in the TSV. What the slide's √L normalizes is not our
   E_s = S_spatial/(4πT) — **open question**, flagged rather than guessed.

### Output files

`output/radial/gauge/gcorr_L{1,2,4}_{sph,exact,flat}_nt120.tsv` — columns `t G_g G_analytic
Delta_eff`, with Δ₀ and the fit window in the header. (`sph` = `gcGeodesic`,
`exact` = `gcExactArea`.)
`output/radial/gauge/scaling_{sph,exact,flat}.tsv` — columns `L Lt abar2 at2 Delta0`, with
Δ₀^cont, c_s, c_t in the header.
`output/radial/gauge/flowscan_nt32.tsv` — columns `s r_over_s E_s E_s_err E_s_sqrtL L g2R`.

### Interface changes (mirrored in [`04-interfaces.md`](04-interfaces.md) §11, §12)

1. **`Gauge`/`newGauge`/`zero`/`gaugeTransform` come from `ops/wilson.nim`** and are
   re-exported by `gaugeact.nim`. WP-G was written before `wilson.nim` landed, with a local
   copy; the local copy is gone and there is one gauge type in the tree. Note
   `wilson.gaugeTransform(l, u, alpha)` is exactly `u += gradient(alpha)`.
2. **`Beta` is new** and is the primary argument of the action/force/heatbath/flow entry
   points; every §11 signature that took `g2: float` still exists as a wrapper.
   `Beta.afac` carries the `A_△` that produced `beta_tri`, so `jtop` and the correlator
   normalization stay consistent with the convention in force.
3. **`GeomConv` and `Beta.conv`** — new; `newBeta(l, g2, conv = gcGeodesic)`. An earlier
   iteration of this interface had `flat: bool`; it became a three-valued enum when
   `gcExactArea` turned out to be the paper's convention (T1.5b).
4. Gauge-field vector-space ops (`:=`, `axpy`, `axpby`, `scale`, `dot`, `norm2`, `toSeq`,
   `fromSeq`, `unitSource`, `nlink`, `slink`) are exported: the CG needs them and they are the
   same names `core/spinor.nim` uses for `Spin`.
5. New: `gaugeActionParts`, `applyM` (template alias for `gaugeForce`), `mDiagonal`,
   `mDiagMax`, `laplace`, `projectFlat`, `projectKernel`, `cgM`, `pseudoSolve`, `RegOp` +
   `newRegOp` + `applyReg` + `regSolve`, `heatbathSource`, `triSource`, a `heatbath` overload
   returning `CgInfo`, and defaulted `r2req`/`maxits` on `gaugePropagator`.
6. **`flowStep`/`flowRun` take `coeffs: RK2NCoeffs[S]`, not `coeffs: auto`** — `auto` cannot
   carry a default and rk.nim's coefficient sets are distinct static-sized types. The
   `g2`-taking `flowRun` of §12 keeps its signature and uses `RK4CK_2N`. `measure` may be `nil`.
7. **`energyDensity` is defined as `S_spatial/(4πT)`** — the action density, i.e.
   ¼⟨F_cF_c⟩ for the canonically normalized field. No factor of g². `energyDensityT` is the
   electric partner; the two together sum to `S/(4πT)`.
8. `ops/gaugeact.nim` imports `ops/solve.nim` **for the `CgInfo` type only**. `cgSolve` there
   is typed on `Spin = seq[array[2, Complex64]]`; packing a real link field into complex
   spinors would waste 4× and buy nothing, so the four CGs here are local and share the
   `CgInfo` return type as the interface requires.

### What did not work, and things to be careful about

1. **The literal (V.16)-(V.17) double CG cannot be pushed to a tight tolerance on a large
   lattice.** It has no control over the kernel component that roundoff injects into the
   Krylov space. Once the in-range residual falls below it, `alpha = r²/(p, Mp)` sees a search
   direction that is mostly kernel — `(p, Mp) → 0` — and the iteration diverges. Measured at
   L=1, n_t=120 (5040 links): `r2req = 1e-22` is fine, `r2req = 1e-28` gives 416864 iterations,
   final relative residual **1.75e+12**, and G_g wrong by four orders of magnitude and
   negative. This is *not* a slow stall, it is a blow-up, and a correlator spanning five
   decades is exactly where you would meet it. (On the 168-link test lattice `pseudoSolve` at
   `r2req = 1e-24` is perfectly happy — 1.5e-14 against the dense pseudo-inverse — so the
   failure needs volume, which is exactly how it would slip through a unit test.)
   The fix, `RegOp`: solve `A = M + σ d d† + τ P Pᵀ/|P|²`. Both added blocks vanish on
   range(M) = ker(M)^⊥ and are positive definite on ker M, so A is SPD on the whole link space
   and `A⁻¹b = M̃⁻¹b` exactly for every transverse b — which every gauge-invariant source is
   (rows of C span range(M)). This is the same "add G G†" regularization the Python oracle
   used. At `r2req = 1e-30` it converges (residual floor 8e-29) and Δ₀ moves by 1e-6 between
   `r2req = 1e-22` and 1e-30 (1.356696 → 1.356697), i.e. the answer is tolerance-saturated
   where the double CG was already breaking. `pseudoSolve` is kept because it is the
   normative algorithm and
   because it also handles *non*-transverse sources (a unit-vector propagator column), which
   `regSolve` does not.
2. **A source that is not orthogonal to ker M must go through `pseudoSolve`, not `regSolve`.**
   `gaugePropagator` takes a unit link vector and therefore uses the double CG. The correlator
   app uses `triSource`, which is a row of the incidence matrix and hence exactly transverse,
   so it uses `regSolve`.
3. **`result` cannot be captured by a closure** in Nim 2.2 ("cannot be captured as it would
   violate memory safety"). Bit both the measurement callback in `rgauge.flowScan` and a nested
   accumulator in the test oracle; accumulate into a local and assign at the end.
4. **`E_s √L` (slide 9) does not correspond to our normalization** — see T2.2 finding 3. Left
   as an open question rather than reverse-engineered.
5. **`getRawMemAllocated()` is the wrong allocation probe here.** It counts QEX's aligned raw
   memory, which this subtree never touches (plain `seq`). The allocation regression uses
   `getOccupiedMem()` instead.
6. **The heatbath tests are statistical and seeded.** `Threefry4x64.seedIndep(51, 0)`, 600
   samples, 5σ acceptance against the *measured* sample error — not a fixed tolerance on a
   random quantity. Both pulls came out below 1σ.
7. **`ā_s` at L = 2 and 4 is 0.590946 and 0.299474**, not the ≈0.55 / ≈0.276 estimated in
   [`02-formulation.md`](02-formulation.md) §6. The §6 table says "approx"; the code values are
   the ones in the scaling fit above.
8. **No `n_max` fit (T1.5e/T1.5f) was done here.** It needs `nmaxFit` against
   `gaugeGPeriodic(t; n_max)` at T = 12 and the residual-normalization question of open item 5,
   which is WP-K's. The correlator data it needs are in the TSVs.

---

## WP-F — overlap operator  (done)

Files: [`ops/overlap.nim`](../ops/overlap.nim), [`tests/toverlap.nim`](../tests/toverlap.nim).
Run: `cd build_mac && SDKROOT=$(xcrun --show-sdk-path) make run experimental/radial/tests/toverlap`.
**17 tests in 10 suites, all pass**, ~1 s wall (the L = 1 dense oracles are cheap by
design).  Dependencies: `core/lattice`, `core/spinor`,
`ops/{wilson,zolotarev,solve}`, `eigens/linalgFuncs` (for `denseOv`); the test adds
`base/alignedMem`.

### Implemented

* \(X = D_{\rm lat} - M\) with the **raw** operator of (IV.1) and plain matrix adjoints,
  \(M\) in raw units (doc/04 §10 settled convention; M = 1 default).  `applyX`/`applyXAdj`
  are one `applyDw`/`applyDwAdj` each; `applyH` = \(X^\dagger X\).
* `applyOv` = \((1{+}m/2) + (1{-}m/2)X R(H)\), `applyOvAdj` =
  \((1{+}m/2) + (1{-}m/2)R(H)X^\dagger\), the standard \(\rho=1\) overlap mass,
  with \(R\) the frozen Zolotarev rational (σ² poles), each **one** `cgmSolve` multishift
  solve at `r2inner`; `applyNormal` = adjoint∘forward (two solves); `solveNormal` = outer CG
  on the normal operator at `r2outer` (2(iters+1) multishift solves). Masses are always-on
  validated in \(0\le m<2\).
* **`ovGradient`** — the single shared pullback, exactly the formula of doc/04 §10:
  two multishift solves (\(s_j=G_j\,\)right, \(t_j=G_jX^\dagger\)left), then
  \(2n_{\rm pole}{+}1\) `dwPullback` calls with weights \((1, -r_j, -r_j)\).  No second
  force kernel exists anywhere.
* `kernelWindow` — power iteration for σ_max, inverse iteration (CG at shift 0, tolerance
  `r2inner`) for σ_min, both from a fixed seeded start with early exit at 1e-8 relative
  eigenresidual; Rayleigh-quotient point estimates, residual-expanded `[lo, hi]`,
  `inside` = `[lo,hi] ⊂ [rat.smin, rat.smax]`.  Monitor only — the window is frozen and
  out-of-window is the caller's hard stop.
* `denseOv` — exact dense \(1 + X(X^\dagger X)^{-1/2}\) via `denseDw` + `zeigs` (zheev)
  eigendecomposition; the test oracle.
* `SolveStats` accounting (`nx`, `nmulti`, `miters`, `mrefits`, `ncg`, `cgiters`, `ok`)
  and `clearStats`.  `stats.ok` goes false if any inner solve misses tolerance; nothing
  raises in the apply path — callers check it.
* Workspace: 6 fixed `work` slots plus the `xs`/`xt` banks, all allocated in `newOv`;
  the two solver closures (`hop` for \(H\), `nop` for the normal operator) are built once
  in `newOv` and reach the current gauge field through a per-call `ptr Gauge`, so the
  apply path performs **zero** allocations in this module.

### Fixtures and windows (all L = 1, M = 1; windows = dense σ bounds padded ±5%)

| fixture | σ(X) range | cond(X†X) | window | maxRelErr(31) | maxRelErr(11) |
|---|---|---|---|---|---|
| nt=6, a_t=0.4, free | [1.529500, 5.821568] | 14.49 | [1.4530, 6.1126] | 1.776e-15 | 1.615e-08 |
| nt=6, a_t=0.4, θ∼0.3N(0,1) | [1.057985, 5.965627] | 31.80 | [1.0051, 6.2639] | 1.332e-15 | 1.806e-07 |
| **production: L_t=60, a_t=0.2, free** | **[1.123410, 10.736484]** | **91.34** | **[1.0672, 11.2733]** | **6.328e-15** | **2.000e-06** |

The production row (Matsubara-block dense, 60 × 24×24) is the number that sizes the tier-2
windows: at M = 1 the free window is ~1/10 and both rational orders are far better than
their 1/200 reference values (WP-B table).  Cross-check tying WP-E and WP-F together:
min |eig(X)| over the blocks = **1.123410**, reproducing WP-E's raw min|D_W−1| = 1.1234 at
L = 1, L_t = 60 to all quoted digits.  Note σ_min = min|eig| here even though X is not
normal; at stronger coupling σ_min will sit below min|eig| — size windows from
`kernelWindow`, not from eigenvalues.

### Measured — every acceptance criterion (r2inner = 1e-26, r2outer = 1e-22)

| # | criterion | asked | measured (free / random) |
|---|---|---|---|
| oracle | \(\|V^\dagger V-1\|\) of the dense polar factor | 1e-12 | **5.22e-15 / 5.66e-15** |
| T1.3c | applyOv vs denseOv, order 31, worst \(\|\Delta v\|/\|v\|\) over 6 vectors | 5 maxRelErr+1e-10 = 1e-10 | **6.86e-15 / 1.86e-14** |
| T1.3c | same, order 11 | 8.1e-8 / 9.0e-7 | **8.39e-09 / 1.29e-07** |
| T1.3d | \(\|(D+D^\dagger-D^\dagger D)v\|/\|v\|\), order 31 | 1e-13 | **1.23e-14 / 6.48e-14** |
| T1.3d | \(\|(D+D^\dagger-DD^\dagger)v\|/\|v\|\) | 1e-13 | **1.08e-14 / 4.22e-14** |
| T1.3e | (IV.22) for γ₄ and γ₅ | reduces to the two rows above | **see reduction note** |
| ladder 1 | \(\langle x,Dy\rangle=\langle D^\dagger x,y\rangle\), mass 0 and 0.13 | 1e-12 | **2.7e-16, 2.1e-16 / 1.7e-15, 1.7e-15** |
| ladder 2 | gauge covariance of D_ov (mass 0.2) | 1e-12 | **2.37e-16** |
| T1.3h | applyOv vs npole independent single-shift CG assemblies | 1e-12 | **9.30e-16 / 1.71e-16** |
| ladder 4 | tangent from the δD_ov formula vs centered FD, best of 5 steps | 1e-8 | **4.35e-10 / 5.68e-11** at ε = 1e-5 |
| ladder 5 | \(\langle f,\delta u\rangle=2{\rm Re}\langle l,\delta D_{\rm ov}r\rangle\) | 1e-12 | **1.45e-15 / 3.77e-16** |
| — | ovGradient scale/add semantics | 1e-12 | **1.78e-15** |
| ladder 6 | Ward: \(\langle f,d\alpha\rangle=2{\rm Re}\langle l,i(\alpha D-D\alpha)r\rangle\) | 1e-12 | **1.40e-13 / 1.31e-13** |
| T1.3f a | free D_lat spectrum conjugation-closed | 1e-10 | **2.08e-14** |
| T1.3f b | rational D_ov spectrum: GW circle \(\||\lambda-1|-1\|\) | 5 maxRelErr+1e-9 | **2.24e-14** (exact denseOv: 1.55e-14) |
| T1.3f b | rational D_ov spectrum conjugation-closed | 1e-8 | **9.14e-15** |
| T1.3f c | free overlap propagator fold \(G(t)=G(n_t{-}t)\), nt=8, a_t=0.25 | 1e-7 rel | **3.12e-12** |
| T1.3f c | Wilson propagator fold violation (the Fig. 10 discriminator) | > 1e-2 and > 1e3× overlap | **7.33e-01**, ratio **2.35e+11** |
| window | kernelWindow σ_min vs dense (400 iters) | 1e-6 rel | **exact to 8 digits** (1.52950010 / 1.05798472) |
| window | kernelWindow σ_max vs dense | 1e-3 rel, bracketed by hi | **5.8215678 exact / 5.96556 vs 5.96563, hi 5.96677 ✓** |
| window | `inside` true on the ±5% window, false on a 1.5×σ_min narrow one | exact | **true / false** |
| alloc | live occupied across GC_fullCollect, 64 × (Ov, Ov†, normal, gradient) | unchanged | **7242656 → 7242656**, raw 0 → 0 |

### Solve counts per call (order 31, the random fixture; `stats` accounting)

| call | multishift solves | recurrence iters | D_W applies |
|---|---|---|---|
| `applyOv` / `applyOvAdj` | **1** | 76 | 183 = 2(76+15)+1 |
| `applyNormal` | **2** | | |
| `ovGradient` | **2** | 151 | 393 = 2(151+30)+2·15+1 |
| `solveNormal` (r2outer 1e-22) | **34 = 2(16+1)** | 16 outer CG iters | |

The +npole in the iteration accounting is `cgmSolve`'s per-shift true-residual recompute;
zero refinement passes fired anywhere in the suite (`mrefits` = 0 throughout, consistent
with WP-D's finding for well-posed tolerances).

### The T1.3d/T1.3e reduction (stated once, tested as stated)

\(\Gamma\mathcal D+\mathcal D\Gamma=\mathcal D\Gamma\mathcal D\) with
\(\mathcal D={\rm diag}(D,D^\dagger)\) reduces for **both** \(\Gamma=\gamma_4\) (offdiag 1,1)
and \(\Gamma=\gamma_5\) (offdiag −i, i) to the same pair of 2-component identities:
γ₄ gives LHS = offdiag\((D{+}D^\dagger, D{+}D^\dagger)\), RHS = offdiag\((DD^\dagger, D^\dagger D)\);
γ₅ the same with factors ∓i.  So (IV.22) ⇔ \(D+D^\dagger-D^\dagger D=0\) **and**
\(D+D^\dagger-DD^\dagger=0\) — the circle identity (equivalently (IV.17)) and its
adjoint-order partner.  Both are tested on both gauge fields; for the exact polar factor
they are identities, for the rational they are violated by \(1-xR(x)^2=-2e-e^2\), i.e. at
2 maxRelErr plus solve residuals — which is what the measured 1e-14 is.

### Interface changes (mirrored in [`04-interfaces.md`](04-interfaces.md) §10)

1. **`SolveStats` is now defined** (§10 referenced it without a definition): `nx`,
   `nmulti`, `miters`, `mrefits`, `ncg`, `cgiters`, `ok`; plus `clearStats(o)`.
   `ok` exists because the apply path must not raise (HMC owns the policy); it is the
   AND of every inner-solve `converged` and every test asserts it.
2. **`Ov` gained `xs`, `xt`**: the two preallocated multishift solution banks (`cgmSolve`
   reuses them allocation-free).  Two banks because `ovGradient` needs \(s_j\) and
   \(t_j\) simultaneously; they cannot live in `work` because `cgmSolve` takes a
   `var seq[Spin]` it may `setLen`.
3. `work` is 6 fixed slots (assignments in the module header), not an open pool — no
   get/put protocol was needed since the call graph is a tree with disjoint slots.
4. Private fields `cu` (a `ptr Gauge` set at each entry), `cmass`, and the two closures
   `hop`/`nop` built once in `newOv` — this is what makes the per-apply allocation count
   of *this module* exactly zero (a fresh closure per call would copy nothing but does
   allocate an environment).
5. `kernelWindow`'s semantics are now written down (Rayleigh estimates, residual-expanded
   margins, early exit, fixed seed, allocates — it is a monitor).
6. `denseOv` returns the full \(1+X(X^\dagger X)^{-1/2}\) (the name says D_ov; §10's
   "polar factor" comment was ambiguous).
7. Aliasing rule documented: `dst` must not alias `src` in any apply (same as `applyDw`).

### What did not work, and things to be careful about

1. **doc/03's T1.3d `< 1e-17` and T1.3e `< 1e-15` are not reachable by a rational
   operator with iterative solves, and are not the right targets.**  The circle residual
   is \(2\,\)maxRelErr + solve error ≥ ~1e-14 even with the rational at its evaluation
   floor (WP-B: e(x) has a ~1e-14 noise floor) and r2inner at the WP-D-safe 1e-26.
   The brief's target — 1e-13 with the order-31 rational at tight r2inner — is met with
   an order of magnitude to spare (1.2e-14 / 6.5e-14).  1e-17 would only be meaningful
   for the dense exact factor, where the measured analogue (\(\|V^\dagger V-1\|\),
   \(\||\lambda-1|-1\|\)) is 5.7e-15 / 1.6e-14.  Whoever owns doc/03 should restate the
   two rows; no tolerance here was loosened — the measured numbers are in the table.
2. **The solver error does NOT get amplified by σ_max, and that is why T1.3d comes out
   an order below the naive estimate.**  Naively \(|X\,\delta z|\le\sigma_{\max}|\delta z|\)
   with \(|\delta z|\sim\sqrt{r_2}\,R(\lambda_{\min})|b|\), which at cond(X) ~ 6 predicts
   ~6e-13.  But \(\delta x_j=G_j r_j\) is *low-mode weighted* — the same \((\lambda+q_j)^{-1}\)
   that amplifies the residual suppresses it under a subsequent X — so
   \(|X\,\delta z|\sim|r|\sim\sqrt{r_{2\rm inner}}|b|\) = 1e-13, and the measured residuals
   sit at 1–6e-14.  Consequence: **r2inner = 1e-26 is enough for everything in this suite;
   nothing needed the floor-adjacent 1e-27 regime** (mrefits = 0 everywhere).
3. **The allocation regression cannot be "occupied unchanged, no collect" here, and the
   test says exactly why.**  `cgmSolve`/`cgSolve` allocate nshift+3 scratch fields per
   call by design (WP-D §9 statelessness; WP-D's own note says hoist *or accept*).  Under
   this build's refc GC that garbage is freed on the collector's cadence, not the loop's,
   so the raw occupied counter drifts (+114816 over 64 rounds, ≈1.8 KB/round net at 5
   poles) even though nothing leaks.  The assertion with teeth is **live memory across
   `GC_fullCollect`: 7242656 → 7242656 exactly**, plus `getRawMemAllocated` 0 → 0, with
   the uncollected drift echoed.  A regression that makes the apply path *grow the live
   set* still fails this test; churn does not.  (First version asserted the raw counter
   and failed on the drift — fixed by measuring live memory, not by widening a bound.)
4. **The propagator fold sign is `+`**: \(G^{(1,1)}(t)=+G^{(1,1)}(n_t-t)\), matching the
   continuum antiperiodic fold \(G_T(T-t)=G_T(t)\) (WP-C).  The overlap satisfies it to
   3.1e-12 (rational + both solver tolerances through \(D^{-1}=(D^\dagger D)^{-1}D^\dagger\));
   the Wilson violation at a_t = 0.25 is **0.73 of the propagator scale** — not a small
   effect: the temporal projectors \((1\mp\sigma_3)/2\) make the spin components nearly
   one-way in time, and only D_ov restores the symmetry (via (IV.17): the T-image of
   \(D_{\rm ov}^{-1}\) is \(1-D_{\rm ov}^{-1}\), a pure contact term).  This is the paper's
   Fig. 10 statement and it discriminates by **11 orders of magnitude** here.
5. **Ward sign**: δ_{dα}D_ov = i(αD_ov − D_ovα), the WP-E/doc-04 §7 convention (opposite
   to the §15 item 6 sketch).  The rational inherits it exactly — the measured 1.3e-13 is
   pure solve residual, and it is the test that would catch a shift-sign error in the
   multishift plumbing.
6. **kernelWindow's σ_max Rayleigh estimate undershoots when the top of spec(H) is
   dense** (random fixture: 5.96556 vs 5.96563 after 400 iterations) — that is intrinsic
   to power iteration, and it is why the residual-expanded `hi` (5.96677, which does cover
   the truth) exists and why `inside` uses `[lo, hi]`, not the point estimates.  σ_min via
   inverse iteration converges to all 8 printed digits in a handful of solves.  Production
   note: pad the frozen window by more than the margins `kernelWindow` reports on the
   *starting* configuration — the ensemble will move σ_min more than any estimator error.
7. **`u` must stay alive across a call** (the solver closures reach it through `ptr Gauge`).
   Passing a temporary that dies mid-call is not possible from Nim call syntax, but do not
   store `Ov.cu` yourself.  The pointer is only read inside the entry that set it.
8. **No external oracle exists for the rational-applied operator itself.**  What stands in
   for one: the exact dense polar factor (independent path: `denseDw` + zheev, verified
   unitary to 5.7e-15), the npole independent single-shift CG assembly (T1.3h), the
   term-by-term tangent rebuilt outside `ovGradient` from `applyDwDeriv`/`cgmSolve` and
   pinned against centered finite differences, and the WP-E cross-check
   min|eig(X)| = 1.1234.  Every number above is measured against one of those.

---

## WP-H — HMC for dynamical overlap + non-compact gauge  (done)

Files: [`hmc/pseudofermion.nim`](../hmc/pseudofermion.nim),
[`hmc/trajectory.nim`](../hmc/trajectory.nim), [`tests/thmc.nim`](../tests/thmc.nim),
[`rhmc.nim`](../rhmc.nim).
Run: `cd build_mac && SDKROOT=$(xcrun --show-sdk-path) make run experimental/radial/tests/thmc`.
**9 suites (10 tests), all pass**, 32 s wall (dominated by the 20-trajectory end-to-end run).
Dependencies: `core/{lattice,spinor}`, `ops/{gaugeact,overlap}` (hence wilson/zolotarev/solve),
QEX `hmc/metropolis` (`RadialHmc` subclasses `MetropolisRoot`), the `mdevolve` nimble package
(one `mkOmelyan2MN` per level over a shared position update, `newParallelEvolution`),
`rng/threefry4x64`.  Demo output: `output/radial/hmc/demo_L1nt60.{log,ckpt}`.

### Implemented

* **Hasenbusch ladder** (doc/02 §4.2), standard mass
  \(D(m)=(1-m/2)D_{\rm ov}+m\),
  \(Q_i=D(m_i)^\dagger D(m_i)\), copies = nf/2, frames per copy over the strictly
  increasing `masses`:
  ratio frame \(i<K\): \(S_i=\phi^\dagger D_{i+1}Q_i^{-1}D_{i+1}^\dagger\phi\), heatbath
  \(\phi=D_{i+1}Q_{i+1}^{-1}D_i^\dagger\xi\) (so \(S_i=|\xi|^2\) exactly, covariance
  \(D_{i+1}^{-\dagger}Q_iD_{i+1}^{-1}\), frame weight \(\det Q_i/\det Q_{i+1}\));
  heaviest frame: \(S_K=\phi^\dagger Q_K^{-1}\phi\), \(\phi=D_K^\dagger\xi\).
  Only Hermitian solves anywhere (`solveNormal`).  `masses.len == 1` = no Hasenbusch.
* **Forces through `ovGradient` only** (the single pullback; both terms of the ratio
  frame are separate `ovGradient` calls, never a second derived kernel):
  heaviest \(dS=-2{\rm Re}\langle y,\delta D_K\eta\rangle\) with
  \(\eta=Q_K^{-1}\phi\), \(y=D_K\eta\); ratio
  \(dS=+2{\rm Re}\langle\phi,\delta D_{i+1}\eta\rangle-2{\rm Re}\langle y,\delta D_i\eta\rangle\)
  with \(\chi=D_{i+1}^\dagger\phi\), \(\eta=Q_i^{-1}\chi\), \(y=D_i\eta\), and
  \(\delta D_i=(1-m_i/2)\delta D_{\rm ov}\).
* **actOp (order 31) / frcOp (order 11) split**: heatbath, both Hamiltonians and the
  accept test use actOp; MD forces use frcOp; orders are never mixed inside one gradient.
* **Momentum**: one unit Gaussian per link (spatial + temporal), kinetic \(|p|^2/2\).
* **Gauge zero modes** (the WP-G machinery, ker M = gauge orbit + Polyakov flat mode):
  `projectKernel` on (a) the committed field — at construction/load and at every
  *commit*, so `reject` stays a bitwise restore of the trajectory's start field;
  (b) the refreshed momentum; (c) **every MD force at every level** (the extended action
  at fixed \(\phi\) has a longitudinal force even though the integrated determinant is
  gauge invariant; projecting only the momentum is not correct).
* **Integrator**: nested Omelyan 2MN via mdevolve — level 0 gauge (innermost, most
  steps), levels 1..nframe = frames heaviest→lightest, one `mkOmelyan2MN` per level over
  the shared position update, combined into a `ParIntegrator` (evolved directly when
  there is a single level: `newParallelEvolution` requires a shared sub-integrator and
  asserts on one-element lists).
* **Metropolis** via QEX `hmc/metropolis.update`; hooks `start` (save field, refresh p
  and phi), `getH` (actOp), `generate` (MD), `globalRand`, `accept` (project + keep),
  `reject` (bitwise restore).  `forceAccept`/`forceReject` for warmup and tests — they
  return constants without consuming any stream, so they cannot desynchronize a restart.
* **Trajectory-addressed randomness**: every draw comes from a Threefry stream seeded by
  a splitmix64 chain over (baseSeed, trajectory number, purpose, copy, frame); purposes
  `rkMomentum/rkAccept/rkPseudo`.  No generator state is serialized; a restart depends
  only on the committed trajectory counter.
* **Checkpoint**: versioned binary manifest (magic, version, mass-convention id, lev, nt,
  at, g2, geometry
  convention, nf, M, both `Rat.hash` values, masses, tau, per-level steps, seed,
  trajectory counter) + the gauge field + a trailing FNV-1a over the whole payload.
  `loadCheckpoint` validates EVERY manifest field against the live configuration and
  raises on any mismatch or on a hash (corruption) failure. Version 2 records
  `standard-overlap-rho1`; legacy version-1 additive checkpoints are rejected.
* **`reversibilityCheck`** (forward τ, p→−p, forward τ, p→−p; reports field/momentum
  drifts, ΔH and the round-trip momentum transversality; restores the state) and
  **`windowCheck`** (kernelWindow on the committed field against BOTH frozen rationals;
  `inside == false` is the app's hard stop, never a rebuild).
* **`rhmc.nim`**: `import base`, qexInit, **freezeTimers()**, full letParam block (lev,
  nt, at, g2R with g2 = g2R at R = 1 and the deck label g2a = g2R/lev echoed, gconv, M,
  ratLo/ratHi, actOrder/frcOrder, nf, masses, tau, steps, innerSteps, ntraj, warmup,
  seed, ckpt/ckptFreq, cfg/measEvery, windowEvery, r2inner/r2outer, maxits, hotStart),
  resume-from-checkpoint when `ckpt` exists, per-trajectory line
  `traj: <n> dH: <x> acc: <0|1> pacc: <p> plaq-like: <S_g/dof> window: [lo,hi] secs: <s>`
  (dof = n_E L_t = rank M, so the free value is 0.5), configuration saves
  `<cfg>.t<traj>` every measEvery in the checkpoint format, hard stop on window
  violation, processSaveParams/writeParamFile/qexFinalize.

### Measured — every acceptance criterion (tests/thmc.nim; L=1, operator fixture n_t=6, a_t=0.4, momentum fixture n_t=4, a_t=0.35)

Fixture window (dense σ envelope of the free and random test fields ±5%):
[1.1639, 5.9467] → maxRelErr(31) = 1.776e-15, maxRelErr(11) = 1.048e-07.  MD window
(free + g²=1 heatbath envelope, padded ×[0.75, 1.30]): [0.7554, 7.8976] →
maxRelErr(31) = 5.551e-15, maxRelErr(11) = 1.921e-06.

| # | criterion | asked | measured |
|---|---|---|---|
| 7 | heatbath identity \(\|S-\|\xi\|^2\|/\|\xi\|^2\), every frame, ladders @[0.0,0.5], @[0.0], @[0.0,0.3,0.8], free AND random field | 1e-8 | **1.117e-14** (worst) |
| 8 | pf action vs dense \(\sum_k\|q_k^\dagger\phi\|^2/\lambda_k\) (denseOv + zeigs, 3-mass ladder, per frame) | 1e-8 rel | **1.223e-14** (worst); totals agree to 13 digits (4.232891408520e+02) |
| 8 | Hasenbusch telescoping \(\sum_i(\ln\det Q_i-\ln\det Q_{i+1})+\ln\det Q_K=\ln\det Q_0\) | 1e-11 | **0.0** (exactly; \(\ln\det Q_0\)=1.763926176845e+02) |
| 9 | per-frame force vs centered FD of the frcOp frame action, best of h = 1e-2..1e-5 | < 1e-5 | **2.097e-09** (ratio frame, m=0), **7.512e-10** (heaviest), both at h = 1e-5 |
| 9 | systematic actOp-vs-frcOp force discrepancy \(\|f_{11}-f_{31}\|/\|f_{31}\|\) | report; O(maxRelErr(11)) | **1.059e-05 = 101.0 × maxRelErr(11)** (m=0 ratio frame), **1.760e-06 = 16.8 ×** (heaviest); pinned < 2e3 × maxRelErr(11) as the regression guard |
| — | momentum transversality after projection | — | \(\|d^\dagger p\|^2/\|p\|^2\) = **2.3e-31**, flat = **9.1e-17** |
| — | \(\langle\|p\|^2\rangle\) = dof = (ne+nv)nt − (nv·nt − 1) − 1 = 168 − 48 = **120** (== WP-G's dim ker M) | 5σ | **119.767 ± 0.480** (pull 0.48, N=1000) |
| 10 | forced reject restores the field | bitwise | **bitwise** (every component ==) |
| 10 | reverse integration drift (steps [4,2,2], τ=0.6, thermalized) | < 1e-6 | du rms **6.83e-17** max 4.44e-16; dp rms **1.73e-16** max 8.88e-16; dH **0.0** |
| 10 | round-trip / post-MD momentum transversality | transverse | divP **5.6e-31 / 2.4e-31**, flat **1.6e-17 / 7.9e-18** |
| 11 | force counts per trajectory, schedule [4, 2, 2] | exact ints | **@[8, 4, 4]** = 2×steps per level (2MN fires V twice per step; mdevolve's lazy flush is always separated by a shared-T update, so counts are exact) |
| 12 | \(|dH_{11}|\propto dt^2\), steps [4,2,2]×{1,2,4}, τ=0.8 | ratios in (2.5, 6.0) | dH11 = −1.056e-2, −3.316e-3, −8.376e-4 → **3.185, 3.958** |
| 12 | dH31 − dH11 (the 11↔31 action mismatch) | dt-independent | **−2.39e-5, −2.41e-5, −2.42e-5** (constant ✓), = 2.1% of the scale 2·maxRelErr(11)·S_pf = 1.17e-3 |
| — | checkpoint round trip | bitwise field, same traj | **bitwise**, counters equal |
| — | restarted chain vs original, 2 trajectories | draws exact, dH to 1e-10 | rnd **bitwise**, accepted equal, worst \(|dH_a-dH_b|\) = **0.0**, final fields **bitwise** |
| — | corrupted (one flipped byte) / mismatched (masses, tau) checkpoints | raise | **ValueError** in all three cases |
| — | end-to-end L=1, n_t=6, N_f=2, @[0.0,0.5], steps [8,2,2], τ=0.6, 3 warmup + 17 | acc > 0.5, \(\langle e^{-dH}\rangle\)≈1 | acceptance **17/17**, \(\langle e^{-dH}\rangle\) = **1.001063 ± 0.001954** (pull 0.54), window inside at every check, **0.66–0.74 s/trajectory** |
| — | allocation regression, 16 × (2 pfForce + projectKernel + pfAction) | live mem unchanged | occupied **3466384 → 3466384** across GC_fullCollect, raw **0 → 0** |

### Demonstration run — L=1, n_t=60, a_t=0.2, N_f=2, g²R=1.0, 30 trajectories

`bin/rhmc -ntraj:30 -warmup:5 -windowEvery:5 -ckptFreq:10 -ckpt:...` (all other
parameters at their defaults: masses @[0.0, 0.5], τ=1.0, steps 4, innerSteps 5 → levels
[20, 4, 4], window [0.3, 12.5], r2inner 1e-22, r2outer 1e-18, hot start from the exact
gauge heatbath).  Log: `output/radial/hmc/demo_L1nt60.log`.

* **41.5 s/trajectory** (serial, one Apple-M-class core; range 38.7–43.7).
* acceptance **24/25 = 0.96** post-warmup (one rejection, trajectory 13 — the reject
  path runs in production too); \(\langle e^{-dH}\rangle\) = **1.002611 ± 0.006289**;
  mean dH = −0.0021, mean |dH| = **0.025**.
* plaq-like S_g/(n_E L_t) fluctuated 0.465–0.590 around the free 0.5.
* kernelWindow along the run: lo ∈ [0.521, 0.627], hi ∈ [10.87, 10.94] — inside
  [0.3, 12.5] at every check, `stats.ok` true throughout.
* solve totals over 30 trajectories: frc **24560 multishift solves** (819/trajectory,
  210 recurrence iterations each), act **7046** (235/trajectory, 214 its); outer normal
  CGs 480 (frc) + 367 (act incl. windowCheck inverse iterations).  Force counts
  **@[1200, 240, 240]** = 2×[20,4,4]×30 exactly.
* **The frozen window had to be widened for the interacting ensemble**: the free kernel
  window at this geometry is σ ∈ [1.1234, 10.74] (WP-F), but on a g²R=1 heatbath start
  `windowCheck` measured σ_min ≈ **0.58** (lo 0.5835) and the [1.0, 12.0] default from
  the planning notes is violated before the first trajectory — the hard stop fired
  exactly as designed on the first attempt.  Defaults are now [0.3, 12.5]
  (maxRelErr(31) = 4.15e-13, maxRelErr(11) = 9.86e-5); the ensemble then stayed
  comfortably inside.  Ensemble-dependent window sizing is real, per WP-F's warning.

### Projected cost — T2.3 condensate campaign (N_f=2, ~300 trajectories × 4 masses per point)

Scaling the measured 41.5 s/trajectory (∝ nsite × CG iterations; the T2.3 ensembles have
m_phys = mR ≥ 0.1, so their light solves are *cheaper* than this m=0 demo — these are
upper bounds, single core):

| point | nsite vs L=1 | est. s/traj | 300 traj × 4 masses |
|---|---|---|---|
| (L=1, g²R=1.5) | ×1 | ~40–50 | **~14–17 h** |
| (L=2, g²R=3.0) | ×3.5 (+ wider window, iters ×~1.5) | ~200–250 | **~70–85 h** |
| (L=4, g²R=6.0) | ×13.5 (+ iters ×~2) | ~1000–1200 | **~330–400 h** |

L=1 is overnight-scale on one core; L=2 is a laptop-week; **L=4 is not practical
serially** — it needs the threading revisit doc/04 §0 reserves for a measurement, a
shorter ladder/trajectory, or simply fewer trajectories.  The σ_min at the stronger
couplings is unknown until measured (`windowCheck` will say; expect windows below 0.3
and correspondingly more expensive rationals and solves at g²R=6).

### Interface changes (mirrored in [`04-interfaces.md`](04-interfaces.md) §13)

1. `Pf` gained `xi2` (the |ξ|² record behind the heatbath-identity test) plus the
   per-frame API the tests need: `ncopy`, `nframe`, `newPf`, `refreshFrame`,
   `frameAction(p, o, u, c, i)`, `frameForce(p, o, f, u, c, i, add)` — `pfAction` and
   `pfForce` keep the §13 signatures (actOp / frcOp respectively).  Frame formulas are
   in the module header and §13.
2. `trajectory.nim`'s API was previously unspecified; now normative: `RadialHmc`
   (subclass of `MetropolisRoot` with l, bt, pf, u, p, tau, seed, traj, fcount, projR2,
   forceAccept/forceReject), `newRadialHmc(l, bt, pf, tau, steps, seed, projR2)` with
   `steps` per level [gauge, heaviest..lightest], `setSteps`, `clearForceCounts`,
   `mdEvolve`, `hmcH`, `refreshMomentum`, `refreshPseudo`, the six metropolis hooks,
   `transversality`, `RevInfo`/`reversibilityCheck`, `windowCheck`,
   `saveCheckpoint`/`loadCheckpoint`, `mixKey`/`keyedRng`/`gaussian(Gauge)` and the
   purpose keys.
3. Force levels are 1..nframe (§13 said "1..K" loosely; with masses m_0..m_K there are
   K+1 frames).
4. Projection placement: rule (a) is implemented at construction/load and at every
   *commit* rather than inside `start`, which is what makes `reject` (and the checkpoint
   round trip) bitwise; the field entering every trajectory is transverse either way.
5. Solver-target guidance added to §13 (see below).

### What did not work, and things to be careful about

1. **WP-F's suite tolerances (r2inner 1e-26, r2outer 1e-22) sit ON the roundoff floor
   once HMC-scale fields and `solveNormal` outer solves are in play.**  Measured on the
   n_t=6 fixtures: multishift true-residual floor ≈ 8e-27 (a clean probe converged to
   8.4e-27 — 1.2× below the request), outer normal-CG floor ≈ 5e-23 (probe 5.07e-23 vs
   the 1e-22 request).  On MD-displaced fields unlucky sources land above the 1.001
   guard, the refinement fallback fires (162 refit iterations in one FD scan) and cannot
   beat the floor → `stats.ok` false with *correct* physics.  Exactly WP-D's "request
   below the floor" scenario.  Fix per WP-D's own guidance (≥2 decades above the floor):
   thmc uses (1e-24, 1e-20), rhmc defaults to (1e-22, 1e-18).  No physics tolerance was
   loosened — every acceptance criterion above is measured at the new targets with
   orders of magnitude to spare, and `mrefits` = 0 everywhere again.
2. **The accept/reject dH saturates at the order-11↔31 action mismatch and this is
   measured, expected, and benign.**  MD conserves the *frcOp* Hamiltonian to dt²; the
   actOp dH differs by [S31−S11](u_f) − [S31−S11](u_i), which is dt-independent —
   measured −2.4e-5 at maxRelErr(11) = 1.9e-6, S_pf ≈ 293 (i.e. ~2% of the 2·ε·S_pf
   scale), constant across step counts to 1%.  The scaling assertion is therefore on
   dH11, with dH31 reported alongside; at the demo window the same 2% fraction predicts
   a ~0.01 dH bias at n_t=60, consistent with the measured mean dH.  Metropolis with the
   order-31 action is exactly what corrects the order-11 flow, so this costs a little
   acceptance and no correctness — but it caps how small dH can be made by refining
   steps, which matters when tuning production schedules.
3. **The dt² window has a lower step bound.**  The first scaling ladder [2,1,1]×{1,2,4}
   at τ=0.8 has gauge dt·√λ_max(M) ≈ 1.2 and gave non-monotone dH (−2.8e-3, −6.7e-3,
   −2.4e-3) — pre-asymptotic, the same effect WP-G saw in the flow-order test at n=2.
   The test starts at [4,2,2], where the ratios are clean.
4. **The planning note "rational windows like [1.0, 12.0] are sane" is free-field
   only.**  At g²R=1, n_t=60 the interacting σ_min ≈ 0.58; the window monitor hard-stopped
   the first demo attempt before trajectory 1 (as designed).  Size production windows
   from `windowCheck` on a representative *interacting* start, with margin below it.
5. **`newParallelEvolution` cannot hold a single integrator** (its shared-integrator
   bookkeeping indexes `shared[0]`); the single-level (pure gauge, no Hasenbusch) case
   evolves the lone 2MN directly.  Also, per-level force counts are exact (2×steps)
   because mdevolve's lazy updater flush is always triggered by the shared position
   update between any two V firings of the same level — pinned by test 11.
6. **The ratio-frame force costs two `ovGradient` calls with the same right vector**
   (s_j is solved twice).  A merged pullback would save one multishift solve per
   ratio-frame force (~15% of a force evaluation) but requires a second derived kernel,
   which §10 forbids — kept the single-pullback contract, cost documented.
7. **`projectGauge`'s tolerance is relative to the divergence of its input**, so
   projecting an already-transverse vector still runs ~10 CG iterations on the site
   Laplacian.  Applied to every force at every level this is measurable but small
   (the fermion force solves dominate by orders of magnitude); the gauge force is
   analytically transverse (d†Mu = 0, flat(Mu) = 0) and its projection is a no-op at
   roundoff, kept anyway because rule (c) is unconditional.
8. **`getH` is evaluated twice per trajectory** (metropolis hOld/hNew); together with
   the heatbath that is 235 actOp multishift solves per demo trajectory, ~22% of the
   solve budget.  A legitimate future economy: at hOld the pseudofermion action equals
   the recorded \(\sum|\xi|^2\) by the heatbath identity (measured 1.1e-14), so the
   fresh solves duplicate known numbers — NOT done here, because spending them keeps
   hOld and hNew on the identical code path (any drift shows up as dH, not as a silent
   bias) and the cost is currently harmless.

---

## WP-I — measurements: currents, scalars, harmonics, gluonic operators, GEVP  (done)

Files: [`meas/harmonics.nim`](../meas/harmonics.nim),
[`meas/observables.nim`](../meas/observables.nim), [`meas/gevp.nim`](../meas/gevp.nim),
[`tests/tmeas.nim`](../tests/tmeas.nim).
Run: `cd build_mac && SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk make run
experimental/radial/tests/tmeas`.
**18 tests in 11 suites, all pass**, ~16 s wall.  Dependencies: `core/lattice`, `core/spinor`,
`ops/{overlap,gaugeact,flow}`, `meas/fit` (WP-J, used not edited), `eigens/linalgFuncs`,
`rng/threefry4x64`.  Interfaces mirrored in [`04-interfaces.md`](04-interfaces.md) §14.

### Implemented

* **`harmonics.nim`** — real orthonormal \(Y_{\ell m}\), \(\ell\le4\), hard-coded polynomials;
  `legendreP` as the independent oracle (addition theorem, worst dev **5.6e-16**); the two
  quadratures `siteProject` (\(\sum_yA_yY x\)) and `faceProject` (\(\sum_fY x\), **no** area —
  doc/07 §4.1); `siteGram`/`faceGram`; **`icosaGroup`** — the 60 rotations of \(I\) with induced
  site/face permutations, equivariance residuals **3.5e-16 / 3.1e-16** (the lattice construction
  is exactly icosahedral).
* **`observables.nim`** — `propSolve` (\(S b=(D^\dagger D)^{-1}D^\dagger b\); adjoint FIRST — the
  order the brief warned about, pinned against the dense inverse), `propSolveDag`; `propagatorT`
  from the 5-fold vertex; `condensatePS` (noise) + `condensateDense` (exact spectrum);
  **the current-kernel face `linkCurrent`** = two `ovGradient` calls (left and \(i\cdot\)left)
  giving complex \(\langle l,K_\ell r\rangle\) on every link at once — the ONLY implementation of
  \(K_\ell\) anywhere, per doc/07 §1.1; `wardChargeScan`; the **factorized connected estimator**
  \({\rm tr}[K_ASK_BS]=E_{\eta\xi}\langle\eta,K_AS\xi\rangle\langle\xi,K_BS\eta\rangle\)
  (`currentSample`/`currentCorrConn`/`currentTraceDisc`) — sequential solves + pullbacks only,
  no tangent kernel, no \(\delta X^\dagger\); dense oracles `denseS` and `denseOvDeriv` (the
  same rational-derivative formula as `ovGradient` in exact linear algebra, tests only);
  `scalarCorrDense`/`scalarCorrPoint`; `jtopProject`, `f2Project`, the 7 `LoopShape`s with
  `loopFlux/loopProject/loopSource` and the **exact** Gaussian correlators
  `jtopCorrExact`/`loopCorrExact` via `regSolve`.
* **`gevp.nim`** — `gevp` (rank-truncating: full rank → the committed `zeigsgv`, else truncate
  the C(t₀) directions below `cut*evmax`, whiten, `zeigs`), `gevpDims`, `gevpCheck`.

### Measured — harmonics and the icosahedral selection rule (T2.6)

**The \(\ell\le2\) Gram blocks are the EXACT identity at every L, not identity + O(a²).**
Two exact statements compose: Schur (each \(\ell\le2\) is an irreducible \(I\)-irrep, so any
invariant bilinear is \(\lambda\,\mathbb1\)) and the addition theorem at coincident points
(\(\sum_mY_{\ell m}^2=(2\ell+1)/4\pi\) pointwise) with \(\sum_yA_y=4\pi\), which forces
\(\lambda=1\) exactly.  Measured at L = 1, 2, 4, both quadratures: diagonal = 1 to 15 digits,
|offdiag| ≤ 9e-17, diagonal spread ≤ 1.6e-15.  The sphere-breaking strength therefore lives
entirely in \(\ell\ge3\) and in the allowed same-irrep couplings:

| quantity (site quadrature) | L=1 | L=2 | L=4 |
|---|---|---|---|
| \(\ell=3\) Gram eigenvalues | {0 (×4), 2.3333 (×3)} | [0.9693, 1.0410], gap split 4+3 | [0.9822, 1.0238] |
| \(\ell=4\) Gram eigenvalues | {0 (×4), 1.8 (×5)} | [0.9693, 1.0246], 4+5 | [0.9822, 1.0143] |
| (2,4) coupling max | 1.265 | 3.885e-02 | 2.256e-02 |
| Schur-zero crosses (0,·),(1,2),(1,3),(1,4),(2,3) | ≤2.4e-16 | ≤1.6e-16 | ≤1.6e-16 |

Face quadrature: \(\ell=3\): [0.2593, 1.5556] (L=1, 3+4), [0.8861, 1.0854] (L=2),
[0.9740, 1.0195] (L=4); (2,4): 0.703 / 0.108 / 0.0247.  Two sharp structural checks came for
free: at L=1 the **site** \(\ell=3\) and \(\ell=4\) Grams have an exactly ZERO 4-fold cluster —
the 12-vertex permutation rep of the icosahedron contains no G irrep — while the face (20-point)
rep does; and every listed number is an \(I_h\)-orbit 5-design statement (first invariant
harmonic is \(\ell=6\)).  Note the L=2→L=4 decrease of the \(\ell\ge3\) deviations is a factor
≈ 1.8–2.0 ≈ \(\bar a_2/\bar a_4\), i.e. **linear** in \(\bar a\), not quadratic — recorded as
measured, not explained.

### Measured — the degeneracy statements on correlators

* **Exact free correlator** (deterministic, `jtopCorrExact` via regSolve at 1e-26; L=1,
  n_t=60, a_t=0.2, g²=1, `gcExactArea`): the \((2\ell+1)^2\) matrices are proportional to the
  identity to **2.9e-16 (\(\ell=1\)) and 7.4e-16 (\(\ell=2\))** relative, over dt=1..4.
  This is claim T2.6 at machine precision on the physical correlator.
* **Single heatbath configuration, group-averaged** (the "any configuration" form): averaging
  \(O^g_m(t)=\sum_fY_{\ell m}(cc_f)\Theta_{gf}(t)\) products over the 60 rotations projects onto
  the invariant part exactly (real orthogonality of the irrep matrices), so ONE configuration
  gives: \(\ell=1\): rel dev **1.2e-16**, \(\ell=2\): **8.3e-16**, \(\ell=3\): exactly two
  clusters 3+4 with spreads ≤ **1.4e-13** against a gap of 112.  Machine precision, not
  statistical — and it tests lattice equivariance (face permutations) + harmonics together.

### Measured — the \(\ell=3\) splitting (the deck's "~3 % at L=2")

Exact correlator, 7×7, clusters are exactly (3-fold, 4-fold) at both L:

| | correlator split (dt=1 / dt=2) | \(\Delta\)(3-fold) | \(\Delta\)(4-fold) | mean | \(\Delta\) split |
|---|---|---|---|---|---|
| L=1 | 126 % / 132 % | 3.4246 | 2.6133 | 2.9610 | **27.4 %** |
| L=2 | 18.0 % / 22.3 % | 3.3508 | 3.1269 | 3.2228 | **6.9 %** |

(\(\Delta\) from λ ratios at t = 0.4–1.6, flat there; free tower \(\sqrt{12}=3.4641\); the
multiplicity-weighted mean at L=2 sits 7 % below it, the expected O(\(\bar a^2\)) artifact.)
Against the deck's "~3 % at L=2": our full spread is 6.9 %, i.e. **±3.4 % as a half-spread** —
consistent if the deck quotes ±, otherwise we measure twice their number at (n_t=60, T=12,
a_t=0.2).  Reported as measured; the convention behind the slide's number is not recoverable
from the deck.

### Measured — fermion sector against dense oracles (L=1, n_t=6, a_t=0.4, dense dim 144)

| test | free | random (θ ~ 0.3 N(0,1)) | asked |
|---|---|---|---|
| `propagatorT` vs dense inverse (mass 0 / 0.13) | 6.8e-16 / 5.0e-16 | 3.8e-12 / 3.3e-12 | 1e-10 |
| condensate: noise (24 vectors) vs exact, mass 0.08 | 0.0041519(166) vs 0.00407084, pull 0.49 | 0.0038254(134) vs 0.00383695, pull 0.09 | 6σ |
| `denseOvDeriv` pinned against `ovGradient` | 7.1e-14 | 5.6e-13 | 1e-11 |
| **Ward** plateau flatness (m=0 / m=0.1) | — | 1.8e-11 / 1.8e-11 | 1e-10 |
| **Ward** jump = i·S_ba (m=0 / m=0.1) | — | 5.8e-11 / 5.7e-11 | 1e-10 |
| one-point tr[K(t)S] time-independence | — | 7.2e-15 | 1e-10 |
| σ_PS vs σ_FS, ALL dt, dense | 8.9e-15 | 1.1e-14 | 1e-11·scale |
| σ_PS vs σ_FS, point source (solver path) | 0.0 | 1.3e-15 | 1e-9·scale |
| point-source scalars vs their dense analogue | 4.9e-15 | 7.6e-13 | 1e-9·scale |
| factorized conn estimator vs dense, worst pull (128 pairs) | — | 2.32 (charge), 1.17 (Y₁₀) | 6σ |
| one-noise trace vs dense tr[KS], worst pull | — | 1.63 / 1.86 | 6σ |
| disconnected product estimator vs dense product | — | pull 0.50 | 10σ (see caveat 6) |

### The Ward/charge-conservation statement — what is exactly true (and what is not)

doc/07 §1.1's \(d^\dagger J=0\) holds **exactly** in these measured forms, configuration by
configuration, at any mass:

1. **one-point**: \({\rm tr}[K_{Q(t)}S]\) is t-independent to 7.2e-15 (pure gauge variance:
   all slice charges are gauge-equivalent one-forms).
2. **fermion line** (`wardChargeScan`): with \(x=Se_a\), \(r=S^\dagger e_b\) anchored at the
   5-fold vertex, \(C(t)=\sum_y\langle r,K_{(y,t)}x\rangle\) is **exactly piecewise constant**
   (flat to 1.8e-11 relative) with jumps \(\mp iS_{ba}\) at \(t_a\), \(t_b\) equal to
   i×(the propagator itself) to 5.8e-11 — the charge between the endpoints exceeds the charge
   outside by exactly the fermion line's unit.  Derivation:
   \(K_{Q(t)}-K_{Q(t')}=i(\Lambda D-D\Lambda)\) for the slice-band indicator \(\Lambda\), so
   \(S(K_{Q}-K_{Q'})S=i(S\Lambda-\Lambda S)\), whose (b,a) element is
   \(iS_{ba}(\Lambda_a-\Lambda_b)\).

**What is NOT exactly conserved: the current–current connected correlator.**  The same band
identity gives \(C_{\rm conn}(t_2)-C_{\rm conn}(t_2')=-i\,{\rm tr}[\Lambda(K_BS-SK_B)]\), which
vanishes only where \(K_B\) has no support — and the overlap current kernel is dense (the GW
contact term is smeared over the operator's exponential-locality range).  Measured densely on
the L=1, n_t=6 lattice, where no separation exceeds that range: the \(\ell=0\) connected
correlator varies by 0.94 (on a scale of 0.29) across \(t_2\ne t_1\).  This is physics of the
non-ultralocal conserved current, not a bug: for the Wilson (ultralocal) current the same
algebra gives exact plateaus.  On production lattices the violation at separated times decays
with the kernel's exponential locality — **WP-L should measure that decay before quoting
vector-channel numbers at small separations.**  The brief's "assert the ℓ=0 charge correlator
is time-independent to 1e-10" is met by the two exact statements above; asserting it on the
current–current correlator would have been wrong, and the dense probe is the evidence.

### σ_PS vs σ_FS — a stronger theorem than "identical spectra"

With \(\langle\xi\eta^\dagger\rangle=S\), \(\langle\eta\xi^\dagger\rangle=S^\dagger\), the two
connected timeslice correlators differ only through \((1-D^\dagger)S^\dagger=S^\dagger-1\)
(the GW contact subtraction).  At \(dt\ne0\) every extra term carries a slice overlap
\(\Pi_2\Pi_1=0\); at \(dt=0\) the would-be contact \(2\sum_{x\in t}{\rm Re\,tr}\,S^\dagger_{xx}-2n_V\)
**also vanishes**, because (IV.17) reads \(S+S^\dagger=1\) at mass 0, so
\({\rm tr}(S+S^\dagger)_{xx}=2\) site by site.  Hence **σ_PS and σ_FS connected correlators are
identical at every separation, configuration by configuration, at m = 0** — measured to 1.1e-14
(dense) and through the solver path to 7.6e-13.  At \(m\ne0\) the second FS contraction scales
by \((1+m)^2\) (from \((1-D^\dagger)(D^\dagger+m)^{-1}=(1+m)S^\dagger-1\)), so the equality is
strictly a massless statement.  Cross-terms \(\langle(\eta^\dagger\xi)(\xi^\dagger\eta)\rangle\)
have no connected contraction (only \(\langle\xi\eta^\dagger\rangle\) pairings exist), so both
correlators are "loop minus nothing" — the disconnected pieces cancel in the non-singlet
flavor assembly as doc/07 §1.2 states.

### Measured — gluonic Monte-Carlo pipeline (L=1, n_t=60, a_t=0.2, g²R=1, exact-area)

128 exact-heatbath configurations, J_top ℓ=1, m- and t-averaged, jackknife errors:

* **MC vs exact correlator: worst pull 2.76 over dt ≤ 20** (e.g. C(0): 0.6384(117) vs
  0.661036; C(5): 0.16452(923) vs 0.174854).  The whole pipeline (heatbath → projection →
  correlation → statistics) is validated against deterministic linear algebra.
* **The free-limit check**: exact \(\Delta_{\rm eff}(t)\) = **1.329856** — flat to 6 digits from
  t = 1 on (see the collapse finding below: the channel holds ONE state), plateau fit over
  \(t\in[3,6)\): **Δ₀ = 1.329856** at (L=1, n_t=60, T=12) against √2 = 1.414214 (−6.0 %, the
  L=1 discretization) and the published exact-area L=1 value 1.33242 at (L_t=120, T=16) — same
  physics, different (T, a_t).  MC: **Δ_eff(t=1) = 1.333 ± 0.115** (pull +0.03 vs exact,
  −0.70σ vs √2); Δ_eff(t=2) = 1.08 ± 0.93 (S/N gone).  With ~100 free-field configurations the
  √2-consistency statement carries an 8 % error at t=1: consistent, and honestly weak — the
  sharp validation is the exact-correlator agreement above.
* \(\Delta_{\ell=3}/\Delta_{\ell=1}\) (exact, L=2): 3.2228/1.3878 = **2.322** vs free
  \(\sqrt6=2.449\) (−5 %, artifacts of both numerator and denominator).

### The 7-shape loop basis — two measured findings that change doc/07 §4.2

1. **Shapes 5–7 must be I_h-covariant and time-reflection-even.**  The raw per-edge temporal
   plaquette carries the arbitrary canonical edge orientation (not I_h-covariant: the Y_lm
   projection is basis-dependent), and a one-sided time difference makes
   \(C_{ij}(dt)\ne C_{ji}(dt)\) (measured asym 0.54) with an indefinite symmetrized C(t₀).
   Implemented: temporal plaquettes/rectangles combined round a face (or rhombus) with the face
   orientation signs, in even second-difference form — lsTPlaq \(=\Theta_f(t{-}1)-2\Theta_f(t)
   +\Theta_f(t{+}1)\), lsTRect2 the extent-2 analogue, lsTRhomb the rhombus one.  All seven are
   exact flux sums: gauge invariance measured at 3.9e-16 under a random gauge transformation.
2. **At L=1 the fixed-ℓ projection collapses every spatial shape onto ONE operator** — T1 has
   multiplicity 1 in the 20-face permutation rep, so on any configuration
   O_rhomb = 2.802517076888·O_tri, O_star = 2.383963416875·O_tri, O_quad = 3.236067977500·O_tri
   (= 2φ = 1+√5 to 12 digits), residuals ≤ 2.3e-16.  doc/07's "drop shape 4 at L=1 and run a
   6×6" is therefore unachievable at L=1: the independent ℓ=1 basis there is
   {lsTri, lsTPlaq, lsTRect2}, and even it spans a **rank-1** physical channel (one state, so
   C(t₀) has two exact null directions).  At L=2 the 7-shape C(t₀) has rank 3 — the channel
   holds 3 states.  Shape 4 (the quadruple triangle: face + its 3 edge-neighbours) exists at
   every L; it is *redundant* at L=1, which we take to be what "degenerates" meant.
3. Positivity needs disjoint supports: with extent-2 shapes, C(t₀) is a genuine transfer-matrix
   Gram only for t₀ > 4 slices (measured: t₀=2 gives evmin −0.86; t₀=5 gives evmin at the
   1e-15 roundoff floor with the physical rank).

### Measured — GEVP

* Synthetic 4×4, 4 states through `zeigsgv`: worst |Δ_n − E_n| = **1.6e-13** (asked 1e-10).
* **Rank truncation is required, not cosmetic** (interface gained `cut`, default 1e-10):
  without it the near-null C(t₀) directions produce garbage generalized eigenvalues that
  overtake the ground state at moderate t (measured: Δ₀ jumping to −2.7/+6.8 at L=2).
  Full-rank inputs still go through the committed `zeigsgv`.
* Exact loop GEVP vs the ℓ=1-projected correlator (t₀=5): ground state agrees to
  **8.4e-5 (L=1, rank 1 kept of 3)** and **5.6e-5 (L=2, rank 3 kept of 7)** for t ∈ [1.2, 2.0];
  at t = 3, 4 the plain-ratio GEVP drifts by 1e-3/1.4e-2 from the arccosh effMass — that is the
  periodic image (T = 12), not a bug; fit windows must respect it.
* Flowed configurations (the deck measures the GEVP on flow times [0.2, 1.6]): exact flowed
  GEVP Δ₀(t=2.2) = **1.329665 (s=0.2)** and **1.324789 (s=0.6)** vs unflowed 1.329856 — the
  flow (which smears in time too) leaves the energy intact once t exceeds the smearing radius,
  and the s=0.6 residual shift at t=2.2 is the visible tail of that smearing.  MC flowed GEVP
  (96 configs, s=0.2, cut=0.05 at the noise level): Δ₀(t=1.6) = 1.61 ± 0.45 jackknife vs exact
  1.326037, pull 0.63.
* `gevpCheck` reports (evmin, evmax, cond, asym) of C(t₀); the L=2 exact basis has
  evmax/kept-evmin ~ 1e4 — small eigenvalues of the whitened problem are accuracy-limited
  accordingly.

### Solver-tolerance findings (measured, for anyone chaining many solves)

WP-F's test values r2inner=1e-26 / r2outer=1e-22 are **guard-floor-adjacent for long chains**
on the random window (cond(X†X) ≈ 39): over ~8000 multishift solves on noise right-hand sides,
(i) the 1.001-guard of cgmSolve's **unrefinable seed system** (j=0; refining it would rebuild
the same vector — WP-D) trips at 1e-26 and still occasionally at 1e-25 (`ok=false` with
`mrefits=0`, the diagnostic signature), and (ii) the outer normal-equation recompute trips at
1e-20 for worst-case noise sources.  tmeas runs 1e-25/1e-20 for propagator/Ward/scalars
(measured accuracy 3.8e-12, two orders inside every acceptance bound) and 1e-23/1e-16 for the
stochastic sampling Ov (solver error ~10 orders below the estimator noise).  No acceptance
tolerance was changed.

### Interface changes (mirrored in [`04-interfaces.md`](04-interfaces.md) §14)

1. `harmonics`: everything returns `float` (real harmonics, doc/07 §1.3), not the sketched
   `Complex64`; added `legendreP` (as WP-C's note requested), `siteGram`/`faceGram`,
   `IcosaGroup`/`icosaGroup`.
2. `faceProject` carries **no area weight** (the J_top area cancellation is built in);
   `siteProject` keeps \(A_y\).  For temporal-link currents the projection weight is
   \(Y_{\ell m}\) alone (`tsliceAmp`): the dual area lives inside \(J_{\rm link}=\partial S/\partial\theta\),
   which is what makes the \(\ell=0\) operator the exactly-conserved total charge.
3. The sketched `linkCurrent(o, u, ...)` became `linkCurrent(o, u, left, right)` → per-link
   complex \(\langle l,K_\ell r\rangle\) via two ovGradient calls; the sketched
   `currentCorr(...)` became the factorized-estimator triple
   `currentSample`/`currentCorrConn`/`currentTraceDisc` (plus `tsliceForm`, `tsliceAmp`,
   `wardChargeScan`); `scalarCorr` became `scalarCorrDense`/`scalarCorrPoint`; `f2Op` became
   `f2Project(l, u, lh, mh)`; `loopOps` gained the (lh, mh) projection and the shape helpers
   `loopCount/loopCenter/loopFlux/loopProject/loopSource`, `jtopCorrExact`, `loopCorrExact`.
4. `gevp` gained `cut` (rank truncation, see above) and `gevpDims`/`gevpCheck` are new;
   full-rank problems still use `eigens/linalgFuncs.zeigsgv` as §14 specified.
5. New dense test oracles `denseS`, `denseOvDeriv` and the solve wrappers
   `propSolve`/`propSolveDag`, plus `fiveFoldSite`.

### What did not work, and honest caveats

1. **The ℓ=0 current–current connected correlator is not exactly conserved** for the overlap
   current (see the Ward section): the exact assertions are the fermion-line plateaus/jump and
   the one-point trace.  The size of the smeared contact at separated times on production
   volumes is unmeasured here (the test lattice has no "separated" region) — WP-L item.
2. **No deterministic point-source estimator of the connected current TRACE exists** within the
   single-kernel rule: tr[K_A S K_B S] needs either a complete basis (the dense oracle covers
   small lattices) or noise.  "Point sources from the 5-fold vertex" are used where they are
   exact: `propagatorT`, `scalarCorrPoint`, `wardChargeScan`.  The volume-noise estimator is
   the production path (`currentSample`, factorized, unbiased — validated at 6σ against dense).
3. **The disconnected piece is as noisy as the deck says**: on the test lattice the dense
   value 3.9e-5 sits under a ±8.6e-3 error from 128 samples; `currentTraceDisc`'s stderr uses
   correlated pairs and is optimistic (documented in the docstring).  The estimator's
   consistency is asserted at 10σ, its mean against the dense product at 0.5σ.
4. **The deck's "~3 % at L=2"** for the ℓ=3 splitting: we measure 6.9 % full spread (±3.4 %) in
   Δ at L=2, n_t=60, T=12 — consistent only under a ± reading; the amplitude splitting is much
   larger (18–22 %).  Quote ours with the convention stated.
5. **Δ_eff from ~100 configs at L=1 has an 8 % error at t=1** and nothing usable beyond t=2;
   "consistent with √2" holds (−0.7σ) but the discriminating validation is MC == exact
   correlator (pull 2.76 worst).  The exact Δ₀(L=1, n_t=60, T=12, exact-area) = 1.329856.
6. **Gram \(\ell\ge3\) deviations shrink ~linearly in \(\bar a\)** between L=2 and L=4
   (factors 1.7–2.0), not quadratically; measured, unexplained, and worth re-checking at L=8
   before anyone leans on an O(a²) extrapolation of icosahedral-breaking effects.
7. **No allocation regression** for `meas/*`: measurement code allocates per call by design
   (dense oracles, per-measurement fields); nothing here sits on the MD path.  The §15
   allocation rule stays with the operator/HMC packages.
8. `scalarCorrDense` is O(n_V²·n_t²) dense contraction — tests only, as marked; the solver
   path (`scalarCorrPoint`) is the production route.

---

## Integration pass — exact-area convention switch  (done, main)

Executed after WP-H/WP-I landed, per "THE COUPLING CONVENTION" section above.

Changes:
* `core/geom.nim`: new exported `kiteArea(le, ls)`; **`Edge.area` is now the exact spherical
  diamond area** (kite sum over both `dl`); comments updated. `dl`/`dual` unchanged, so the flat
  form remains recoverable as `0.5*len*dual`.
* `core/lattice.nim`: `kap` and `betaEdge` automatically became exact (they read `Edge.area`);
  comments updated — the (IV.2) identity `kap = dual/abar` is now marked O(abar²).
* `ops/gaugeact.nim`: local `dualExact` removed in favour of `geom.kiteArea`;
  **`newBeta` default flipped to `gcExactArea`**; `gcGeodesic` now computes the flat area
  inline (`dual/(g2*len*at)`), preserving its meaning as the slide-8-legend convention.
* Test updates, all with citations in-file:
  `tgeom`: diamond-area pin flipped to the exact form + new checks (ΣA_e = 4π < 1e-12; flat form
  O(abar²) away, bracketed both sides); κ-oracle pins → exact values (0.683450 …);
  flat `kap = dual/abar` equivalence dropped.
  `twilson`: L=1 doubler-vs-flat bound 0.02 → 0.025 (measured shift 1.51 % → 2.13 %).
  `toverlap`: production σ_min pin 1.1234 → **1.1770**; the measured shift +0.0536 equals the
  Python-oracle prediction exactly.
  `tgauge`: the two convention-comparison tests now request `gcGeodesic` explicitly.
  `tflow`: hand-built `Beta` literal inherits `bet.conv`.
  `tmeas`: fixture `r2in` 1e-25 → 1e-24 (the switch nudged the roundoff floor; same marginal
  ok=false/mrefits=0 signature as WP-I's original fix, same one-decade remedy).

Verification: **all 12 suites re-run green — 196 [OK], 0 [FAILED]**
(tgeom 30, tzolo 11, tanalytic 18, tspinor 6, tsolve 6, twilson 25, toverlap 17, tgauge 28,
tflow 8, thmc 10, tmeas 18, tfit 19).

---

## WP-K — free-limit campaign  (done)

Files: [`rfree.nim`](../rfree.nim) (fermion + gauge campaign, and the shared dense
Matsubara-mode machinery), [`rspec.nim`](../rspec.nim) (Fig. 4/5/6 spectra and the slide-8
legend table; imports rfree's helpers), [`campaign/free/*.gp`](../campaign/free/) (gnuplot,
run from the repo root).  One sanctioned addition outside `radial/`: **`zgesv` bound in
`src/eigens/lapack.nim`** (per the "Known facts" section above; the thin `zsolve` wrapper
lives in `rfree.nim`, not in linalgFuncs).
Run: `cd build_mac && SDKROOT=$(xcrun --show-sdk-path) make run experimental/radial/rfree`
(then `… rspec`).  Everything is deterministic; every table is TSV under
`output/radial/free/` (30 files + a `cache/` of raw correlators keyed by (L, L_t, T) so an
interrupted campaign resumes at the last completed point — delete `cache/` to force a full
recompute).  All couplings are the tree-wide exact-kite-area convention.

### Method

Free fields are t-translation invariant, so both sectors are dense linear algebra per
temporal Matsubara mode — no CG anywhere:

* **Fermion** (antiperiodic k = (2n+1)π/L_t): the spatial 2N_V×2N_V block from `denseDw`
  on an n_t = 1 lattice with `parts = dwSpace` (WP-E's construction, pinned to 5.8e-14),
  then D_W(k) = D_spatial + κ'_y[(1−e^{−ik})P_↑ + (1−e^{+ik})P_↓] per mode; X = D_W(k) − 1,
  H = X†X, `zheev` gives H^{−1/2}, D_ov(k) = 1 + XH^{−1/2}, one `zgesv` column solve, and
  G^{(1,1)}(t) = (1/(ā a_t L_t)) Σ_k e^{ikt}[D_ov(k)^{−1}]_{(src,↑),(src,↑)}, src the 5-fold
  vertex (`fiveFoldSite`).  The raw Wilson propagator (Fig. 10) is the same sum on
  D_W(k)^{−1}.  L = 8, L_t = 168 is 168 modes of dimension 1284 — ~11 min.
* **Gauge** (periodic k = 2πn/L_t): the (N_E+N_V)×(N_E+N_V) Hermitian per-mode form
  M(k) = Σ_p β_p conj(w_p)w_p^T (faces: incidence rows; temporal plaquette of edge a→b:
  [(1−e^{ik}) at the edge, +1 at b, −1 at a]), kernel-regularized with σG(k)G(k)†
  (+ τ uu^T/N_V at k = 0) — the per-mode image of WP-G's `RegOp`, exact for the
  gauge-invariant face sources.  One `zgesv` with all N_F face sources as right-hand
  sides per mode; F(−k) = F(k), so only k ∈ [0, π].  L = 8, L_t = 120 is 61 modes of
  dimension 2562 — ~105 s.

**Validation before any campaign number** (all in-app, `doCheck`):
fermion modes vs the exact dense `denseOv` on a full L=1, n_t=6 lattice: **1.2e-15**;
vs dense `denseDw` (raw Wilson — this pins the sign of k in e^{ikt}, which the T-symmetric
overlap cannot): **2.3e-15**; gauge modes vs the real-space `regSolve` pseudo-inverse:
**6.5e-13**; gauge vs WP-G's pinned exact-area L=1/L_t=120/T=16 table: **2.5e-9** (the
oracle's own quoted precision).  Mode-sum reality: fermion max|Im G|/max|Re G| ≈ 1e-14,
gauge contraction Im ≤ 1.4e-18.  And across the full T=16 grid the per-mode pipeline
reproduces WP-G's real-space Δ₀ at all 12 shared (L ≤ 4, L_t) points to every printed
digit (worst 1.405648 vs 1.405649).

### MASTER TABLE (published | ours | deviation)

| quantity | published | ours | deviation |
|---|---|---|---|
| fermion Δ₀ (L=1, L_t=168, T=16), T1.4c | **0.953918** | **0.953918** ± 0.000174 (fit) | **−1.4e-7** |
| fermion Δ₀^cont (V.7), T1.4e | **0.999998(34)**, exact 1 | **0.999999 ± 0.000044** (syst −2.8e-5) | **+1e-6**; −8e-7 from exact |
| fermion n_max (L=1,2,4,8), T1.4f | **6, 10, 19, 32** | **3, 7, 14, 28** (rel weight; abs 11/13/17/22, log 4/7/14/28) | **integers differ — see below** |
| fermion res/DOF, T1.4g | 0.028, 0.012, 0.0039, 0.038 | at OUR minima: 0.0115, 0.0039, 0.0014, 0.0032; **at the PUBLISHED n_max: 0.0279, 0.0102, 0.0040, 0.0034** | pub column = ours@pub to 2 digits (L=2 15 % off); L=8 pub 0.038 ≈ our 0.0034 ⇒ likely a typo for 0.0038 |
| gauge Δ₀ (L=1, L_t=120, T=16), T1.5b | **1.33242** | **1.332430** ± 0.000109 (fit) | **+7.5e-6 rel** |
| gauge Δ₀^cont (V.7), T1.5d | **1.41409(18)**, exact √2 | **1.414208 ± 0.000070** (paper grid L=2,4,8 × L_t=48..120) | +0.7σ_pub; **−5.5e-6 from √2**; WP-G's L≤4 grid from this pipeline: 1.414535(286) = WP-G's own 1.414535(286) exactly |
| gauge n_max, T1.5e | **3, 8, 18, 35** | **3, 7, 16, 30** (rel; abs 3/7/14/26, log 3/7/16/30) | L=1 exact; higher L 1–5 below published |
| gauge res/DOF, T1.5f | 0.0031, 0.0023, 0.0031, 0.0037 | at our minima: 0.0031, 0.0021, 0.0020, 0.0025; at published n_max: 0.0031, 0.0026, 0.0049, 0.0063 | L=1 exact (both readings); higher L do not track either way |
| Fig 7 normalization (max \|G/G_ana−1\|, t∈[2,6], T=12) | agreement incl. normalization | L=1..8: **9.9 %, 7.0 %, 5.1 %, 1.6 %** | O(ā²)-consistent; max sits at t≈2 (UV end) |
| Fig 11 normalization (same window) | agreement incl. normalization | L=1..8: **55.6 %, 13.0 %, 3.4 %, 1.2 %** | matches WP-G's L=1 ratio 1.28–1.52 finding |
| Fig 10 fold violation max\|G(t)−G(T−t)\|/max\|G\| | Wilson violates, overlap preserves | Wilson **0.87**, overlap **8.8e-14** | 13 orders of magnitude |
| slide-8 legend min\|D_W−1\| (flat κ, L_t=60) | 1.154 (L=1, a_t≈0.133), 1.010, 0.965 | **1.158156** (a_t=2/15), **1.010466**, **0.964330** (a_t=0.2); exact-κ: 1.210339, 1.015420, 0.966115 | reproduces the resolved-OQ-2 numbers to all digits |
| D_ov Ginsparg-Wilson circle (Fig. 6, L=4, L_t=24) | \|λ−1\| = 1 | max deviation **3.2e-14** | exact-eigendecomposition floor |
| wall clock | — | fermion **3116 s** (L=8: 3083 s), gauge **314 s** (L=8: 305 s), rspec **46 s**; total ≈ **58 min** serial | Accelerate zheev/zgesv |

Supporting grids (all in the TSV headers too):
fermion Δ₀(L, L_t=168): 0.953918, 0.985052, 0.995093, 0.997744;
fit c_s = −0.0387, c_t = −0.1571, rms 3.4e-5 (dof 6; central grid excludes L=1 and
L_t ∈ {96}, per the paper; including L_t=96 shifts by −2.8e-5, quoted as the systematic).
Gauge Δ₀(L, L_t=120): 1.332430, 1.390464, 1.406795, 1.410870; fit c_s = −0.0620,
c_t = −0.1115, rms 1.15e-4 (dof 9).

### The n_max finding (open question 5, closed for the free campaign)

The published integers are **not** reproduced as the minimizer of (V.9) on our correlators
under any convention scanned (`-nmaxScan:1` echoes every curve): weighting ∈ {relative
1/g², absolute, log} × window ∈ {all t>0 (167/119 pts), t∈[a_t, T/2], all t incl. t=0} ×
{C free, C≡1} × model ∈ {periodic images, infinite-volume} × {correlator, Δ_eff-based}.
Our minima are stable across rel/log but shallow (the residual only doubles a few integers
away), and sit 1–5 below the published column at every L except gauge L=1 (exact match,
3 @ 0.0031 — also the one case the coordinator's Python oracle checked, confirming it).

What DOES reproduce for the fermion: **the published residual column equals our relative
res/DOF evaluated AT the published n_max** — 0.0279 vs 0.028 (L=1), 0.0040 vs 0.0039
(L=4), and 0.0034 vs the published 0.038 at L=8, which restores monotonicity of their
column if read as a typo for 0.0038 (L=2: 0.0102 vs 0.012, the one imperfect cell).  So
the correlators themselves agree with the paper's — only the selection rule for the
integer differs, and it is not recoverable from the published information.  The gauge
column does not show the same pattern beyond L=2 (our res@pub 0.0049/0.0063 vs published
0.0031/0.0037 at L=4/8).  Bottom line per the resolved-OQ-5 guidance: we quote our minima
under the stated (relative, C free, t>0) convention, and the published integers with the
residuals our data assigns them; nothing was tuned.

**A byproduct worth keeping: the t=0 point of the (1,1) overlap propagator is pure GW
contact.**  (IV.17) forces Re[D_ov^{−1}]_{xx} = 1/2 per diagonal entry, and the T-odd
physical part vanishes at t = 0, so G(0)·ā·a_t − 1/2 = −7.4e-17 / −8.0e-17 / −9.6e-18 /
−5.0e-17 at L=1/2/4/8 (measured).  Fitting that point raw drags the n_max minimum down;
"subtracting the contact" leaves an exact zero that a relative weight cannot hold — hence
the t>0 convention (DOF = L_t − 3 … strictly npoints−2 = 165/118).

### Figures regenerated (T1.2f, T1.2g, T1.3g, T1.4a/b/d/h, T1.5a/c)

`fig4_L{1,2,4}_Lt24, fig4_L2_Lt{16,48}` (raw D_W spectra, T=4; min|λ−1| = 1.3566, 1.0716,
0.9803, 1.0652, 1.0778), `fig5_flat_*` (IV.8 companions), `fig6_{dw,dov,dw_gen,dov_gen}`
(T=4, L=4, L_t=24, M=1; the gen files use the **corrected** weight diag(volbar/volw) of
doc/02 §3.2), `fig7_L{1,2,4,8}` (+ periodic and infinite-volume analytic columns and the
ratio), `fig8_L{1,2,4,8}` (Δ_eff + analytic Δ_eff), `fig9_fermion_scaling` /
`fig12_gauge_scaling` (grids + both projections `proj_s = Δ₀ − c_t a_t²`,
`proj_t = Δ₀ − c_s ā²`, fit in the header), `fig10_L1` (overlap vs Wilson),
`fig11_L{1,2,4,8}`, `slide8_legend` (both κ conventions × {L=1,2,4 at a_t=0.2, L_t=60;
L=1 at a_t=2/15, L_t=60 and 120}).  One-command plots: `gnuplot
src/experimental/radial/campaign/free/fig7.gp` etc. from the repo root.

### Slide-8 legend cross-checks (resolved open question 2, now from this pipeline)

flat κ, L_t=60: L=1 a_t=0.2 → **1.123410** (= WP-E's 1.1234), L=2 → **1.010466**
(pub 1.010), L=4 → **0.964330** (pub 0.965), L=1 a_t=2/15 → **1.158156** (= the oracle's
1.1582, pub 1.154); exact κ: 1.176991 (= toverlap's production pin 1.1770), 1.015420,
0.966115, 1.210339 (= the oracle's 1.2104).  min|D_W−1| also depends on L_t (1.1099 at
L=1, a_t=2/15, L_t=120), so the legend comparison is quoted at the oracle's L_t=60.

### What did not work, and things to be careful about

1. **The published fermion/gauge n_max integers do not come out of (V.9) on our data**
   (above).  Everything else in the free-limit target list reproduces.
2. **The (V.6) fit is what reaches the published precision** — exactly as WP-C warned:
   the raw Δ_eff at T=16 is ~1e-3 away at the window's end; the plateau-fit χ²/dof of
   3–8e-9 across the whole grid is what makes 0.999999/1.414208 possible.
3. **Fermion mode pairs cannot be halved.**  σ₁·conj(D(k))·σ₁ = D(k) is a SAME-k identity
   (it forces [D(k)^{−1}]_{↑↑} = conj([D(k)^{−1}]_{↓↓}), and per-mode conjugate eigenvalue
   pairing), NOT a k→−k map, so all L_t fermion modes are computed; the gauge sector does
   halve (M(−k) = conj(M(k)) with real sources).
4. **The L=1 normalization ratios are large** (Fig 7: 9.9 %, Fig 11: 56 % max in t∈[2,6])
   and shrink with L — consistent with WP-G's L=1 ratio table; "agreement including
   normalization" is an L≥4 statement on a log axis.  The maxima sit at the t≈2 end.
5. `zeigsgv` was not used for Fig. 6: with the row weight the problem is
   diag(w)·D — a plain (non-Hermitian) eigenproblem, `zgeigs` — the generalized-Hermitian
   route does not apply to D_W.
6. The correlator `cache/` files are keyed by (sector, L, L_t, T) only — the convention
   (exact kite area) and source site are baked in.  If either ever changes, delete the
   cache; nothing checks it.
7. Wall-clock is dominated by `zheev` at L=8 (1284² per fermion mode, 168+4×(96..168)/…
   modes ≈ 3100 s).  The gauge sector's 2562-dimensional zgesv with 1280 right-hand
   sides is only ~105 s per (L=8, L_t=120) point — LU + back-substitution beats an
   eigendecomposition by an order of magnitude, which is why the gauge L=8 row was cheap
   after all.

---

## WP-L — interacting-campaign machinery: rmeas + campaign/t2.sh  (done; ensembles NOT launched)

Files: [`rmeas.nim`](../rmeas.nim), [`campaign/t2.sh`](../campaign/t2.sh).
Build: `cd build_mac && SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk make
experimental/radial/rmeas` (rhmc untouched).  Everything below was run on 2026-08-21 on this
machine; per the WP-L brief nothing longer than ~10 minutes was started — the priority-ensemble
launch commands are at the end for the coordinator.  Smoke/calibration outputs are under
`output/radial/t2/{smoke,caldemo,cal60,calL2,calL4}`.

### Implemented — `rmeas`, the offline measurement driver

Reads `<dir>/<cfg>.t<traj>` configurations through **WP-H's own reader**
(`hmc/trajectory.loadCheckpoint` on a `RadialHmc` built from the command line), so **every
manifest field is validated** — lev, nt, at, g2, convention, nf, M, both rational hashes, masses,
the standard-overlap mass-convention id,
tau, per-level steps, seed — and a mismatch refuses to measure
(verified: `-g2R:2.5` against the g2R=1.5 smoke ensemble dies with
`checkpoint mismatch: g2 [ValueError]`, exit 1).  Per configuration, selected by `-obs:` (comma
list or `all`), it writes one TSV per observable per configuration under `<dir>/meas/`, each with
a full `#key=value` manifest (ensemble parameters, `convention=gcExactArea`,
`massConvention=standard-overlap-rho1`, `overlapRho=1`, both rational hashes, valence mass,
estimator notes). Existing files are skipped (`-skipDone:false` remeasures)
and measurement randomness is trajectory-addressed (`keyedRng(seed, traj, purpose)`, purposes
101/102/104 disjoint from HMC's 1–3), so a re-run is bit-reproducible — verified: the re-measure
after the interruption reproduced every number to the printed digit.

* `-obs:cond` (T2.3): `condensatePS` at the **sea** mass (`-measMass` overrides), `-nnoise`
  Gaussian vectors, noise tolerances (1e-22, 1e-16).
* `-obs:currents` (T2.4/T2.5/T2.6): the factorized noise estimator (`currentSample`, `-npair`
  pairs) for the temporal link-current correlators projected on Y_lm, l = 0..`-lmaxJ` (3), ALL m
  (l = 0, the conserved charge, is always included as the conservation diagnostic).  Stored per
  (l, m): the t1-averaged, dt↔nt−dt-folded, **pairing-averaged** product
  k(dt) — both pairings E[a(t2)b(t1)] and E[b(t2)a(t1)] equal tr[K(t2) S K(t1) S], so averaging
  them is free variance reduction.  Physical connected correlator per 4-component flavor
  = −2 Re k.  `-disc:true` adds `currdisc` (unbiased cross-sample T(t2)T(t1) products,
  T = 2 Re tr[K S]) and `currtrace` (per-slice traces) for the vector channel.  `-ward:true`
  (default) runs `wardChargeScan` per configuration and records plateau-flatness and
  jump-deviation in the manifest.  **Point sources are used only where they are exact** (Ward,
  scalars): WP-I established there is no deterministic point-source estimator of the connected
  current trace within the single-kernel rule, so the correlators are stochastic-only.
* `-obs:scalars` (T2.9): `scalarCorrPoint` from the 5-fold vertex, `-nscalarSrc` source slices.
  Records `psfs_maxdev` (dt≠0) and `psfs_contact` (dt=0).  **Finding: in this contraction the
  dt≠0 equality is bitwise-trivial** (the second FS contraction is exactly the complex conjugate
  of the PS one term by term), so the m = 0 statement with content is the dt = 0 GW contact
  cancellation — measured 8e-15…6e-14 per configuration (solver-level, as WP-I's dense theorem
  requires).
* `-obs:gluon` (T2.2/T2.7/T2.8): flows a COPY of the configuration with the
  **coupling-independent flow** (`newBeta(l, 1.0)`, the same Luscher normalization as rgauge's
  T2.2 scan; flow time is a length²) to `-flowTimes` (default the deck's 0.2..1.6 step 0.2, plus
  s = 0 first), RK4CK at h = 0.1/maxᵢMᵢᵢ.  Per flow time: E_s and E_t with the **ensemble** Beta
  (slide 9); the Y_lm-projected Wilson-loop correlator matrices for l = 1 and l = 2, m-averaged,
  folded and symmetrized (`loops`); the F² l = 0 correlator plus its slice mean, RAW
  (`f2`; vacuum subtraction at analysis).  Loop basis: all 7 shapes at L ≥ 2,
  {lsTri, lsTPlaq, lsTRect2} at L = 1 (the WP-I collapse).
* `-obs:wspec` (T2.1): dense `denseDw(l, u, M)` + `zgeigs`; writes eig(D_W) with
  |D_W − M| < `-wspecRad` (1.5) and `min_abs` = min|D_W − M| in the manifest; **skips with a
  message** when 2·nsite > `-wspecMaxDim` (5100, sized so L=2/n_t=60 = 5040 runs).
* Guard: when fermionic observables are requested, `kernelWindow` is checked per configuration
  against the frozen window and a violation prints a loud warning (quenched/pure configurations
  were never monitored by rhmc; measured: they sit visibly lower, see the calibrations).
* `-pure:true -nconf:N`: generates N exact-heatbath configurations in memory (WP-G heatbath, one
  keyed stream per configuration — skip-done safe) and measures them: the slide-9 free baselines.
* `-ckptInfo:true`: prints the committed trajectory counter of `<dir>/ckpt` through the validated
  reader and exits — the campaign resume primitive.
* `-analyze:true` aggregates `<dir>/meas/` into `<dir>/analysis/`: ensemble means with
  **delete-block jackknife** (`-jackBs`), two Δ estimators everywhere — `fit` = `plateauFit`
  (V.6) on the effective mass over `[fitLoT, fitHiT)` (NaN-trimmed longest finite run, needs ≥ 4
  points) and `eff` = the effective mass at `-trefT` (the honest low-statistics fallback) —
  the loop GEVP per flow time (`gevp` with `-gevpCut` 0.05 rank truncation, `-gevpT0` 5, Δ read
  at `-gevpTref` 1.6), the vector assembly, condensate `fit.jack` statistics, and
  `summary.tsv` with the deck's ratios keyed by (g2R, lev, nf) for cross-ensemble plots.
  Outputs: `cond.tsv`, `curr_corr.tsv`, `curr_effmass.tsv`, `vector.tsv`, `scalars.tsv`,
  `flow.tsv`, `gevp.tsv`, `wspec_eigs.tsv`, `summary.tsv`.

### Analysis conventions (decisions, stated once)

1. **Monte-Carlo effective masses are the LOCAL log ratio** ln(c(t)/c(t+a_t))/a_t, not the
   paper's (V.4) arccosh form.  (V.4) references c(T/2), which is exact for the deterministic
   free-limit correlators but is unresolvable noise for MC data (e^{−ΔT/2} ≈ 1e-5 at Δ = 2,
   T = 12); the log ratio's periodic-image bias e^{−Δ(T−2t)} is far below any statistical error
   in the fit windows.  `meas/fit.effMass` is untouched (WP-J/WP-K keep using it).
2. **Vector channel**: C_V(dt) ∝ −2 Re k(dt) + N_f·(P(dt) − ⟨T⟩²) per 4-component flavor
   (doc/07 §1.2 contractions with T_vector = 1, T_axial = τ₃; the axial IS the connected piece).
   ⟨T⟩² uses per-m ensemble means; P is the unbiased cross-sample product.
3. **l = 3** is reported per-m (the deck's own slide-13 presentation): per-m Δ estimates, their
   multiplicity-unweighted mean, and the full spread (max−min)/mean as `l3_spread_over_m`.
4. **F²** is vacuum-subtracted at analysis time: C(dt) − ⟨Ō⟩², jackknifed jointly.
5. Ratios between channels use replicate-wise jackknife when the two config lists coincide,
   independent-error propagation otherwise (the summary is built so they coincide).
6. Solver targets: point solves (1e-23, 1e-18), noise estimators (1e-22, 1e-16), per WP-I's
   findings; the 1e-24 default tripped cgmSolve's benign unrefinable-seed guard
   (ok=false, mrefits=0) on occasional sources and was moved one decade up.  One marginal trip
   still appeared on one smoke source — same signature, values unchanged to 1e-14; it is
   sensitive to code layout under `-Ofast -ffast-math` (bit-level loop reassociation), so treat
   it as cosmetic unless mrefits > 0.

### Smoke chain (the mandatory verification) — green end to end

`rhmc -lev:1 -nt:16 -at:0.2 -g2R:1.5 -nf:2 -ntraj:10 -warmup:2 -measEvery:2 -ckptFreq:5
-windowEvery:2 -seed:11` → 10 trajectories at **6.98 s/trajectory**, acceptance 8/8,
⟨e^{−dH}⟩ = 0.9985(43), window [0.77, 11.03] ⊂ [0.3, 12.5] throughout (σ_min at g2R = 1.5,
n_t = 16 stayed at 0.66–0.84 — the default window holds at this coupling).  Then `rmeas
-obs:all -disc:true -tmin:3` on the 4 post-warmup configurations, then `-analyze:true`:

* **cond**: σ_PS(m=0) = −8.7e-15 ± 1.1e-14 — **exactly zero is correct**: at m = 0 the GW
  relation makes tr[(1−D/2)D⁻¹] vanish identically (that is why the deck scans m > 0);
  the noise estimator confirms it at solver precision.  T2.3 lives on the m > 0 ensembles.
* **Ward / conserved charge**: fermion-line plateaus flat to 4.3e-10…1.8e-9 of the jump, jump =
  i·S_ba to 1.2e-9…5.9e-9, every configuration.
* **PS == FS**: dt≠0 bitwise (trivially, see above); the m = 0 dt = 0 GW contact cancels to
  0.8–5.7e-14 per configuration.
* **wspec**: min|D_W−1| = 1.050 ± 0.019 over the 4 configs (free value at L=1, L_t=16 is 1.4287
  — the additive shift toward the published interacting trend is already visible at 10
  trajectories).
* **gluon/GEVP/F²/flow**: all matrices finite, E_s monotone in s, GEVP rank truncation engages
  (L=1 l=1 channel), Δ_F(l=1, s=0.6) = 1.92 ± 0.20 (jackknife over 4 configs).
* Analysis behavior at 4 configs × n_t = 16 is honest, not pretty: the `fit` estimator column is
  NaN by design (the window [1.0, T/2−a_t) holds only 2 points at T = 3.2; ≥ 4 required), the
  axial `eff` at t = 1 is NaN because the stochastic correlator is noise there at 32 samples,
  and Δ_V is nonsense with huge error — the disconnected piece at 4 configs is exactly as noisy
  as the deck says.  **No unflagged NaN anywhere: every raw measurement is finite-checked at
  write time; analysis NaNs mean "estimator unavailable at this statistics", by construction.**

### Measured timings (the calibration next to WP-H's 41.5 s/trajectory)

Per configuration, serial, this machine (WP-K running concurrently, so ±20%):

| lattice | config type | cond (16 noise) | currents (16 ops × 8 pairs + disc + ward) | scalars (2 src) | gluon (9 flow times) | wspec (dense) |
|---|---|---|---|---|---|---|
| L=1, n_t=16, m=0 | dynamical (smoke) | 4.6–5.0 s | 6.1–6.6 s | 2.6–2.8 s | 0.01 s | 0.10 s (dim 384) |
| L=1, n_t=60, m=0 | **dynamical** (WP-H demo config, g²R=1) | 140 s | 115 s | 28 s | 0.03 s | 3.9 s (dim 1440) |
| L=1, n_t=60, m=0 | quenched heatbath g²R=1.5 | 124 s | 184 s | 110 s | 0.09 s | 8.2 s |
| L=2, n_t=60 | quenched heatbath g²R=3 | — | — | — | 0.22 s (7 shapes) | **79.8 s** (dim 5040) |
| L=4, n_t=60 | pure heatbath g²R=4 | — | — | — | 0.90 s | skipped (dim 19440 > cap) |

So a full fermionic measurement on an L=1/n_t=60 dynamical configuration costs ≈ 150 s without
cond (≈ 3.5 trajectories), ≈ 290 s with it.  The quenched σ_min sits well below the dynamical
one (min|D_W−1|: 0.78 dynamical g²R=1 vs 0.67 quenched g²R=1.5 vs 0.49 quenched g²R=3), which
is why quenched fermion solves are 2–4× slower — the pure ensembles therefore only measure
gluonic observables, and the rmeas window guard warns if anyone tries otherwise.
Physics spot-check for free: min|D_W−1| = **0.7813** on the single thermalized g²R=1, L=1,
n_t=60 demo configuration, against the published slide-8 legend **0.814** at that point.

### The campaign — `campaign/t2.sh` (resume-safe, sizes measured or derived from measurements)

`bash src/experimental/radial/campaign/t2.sh [ensemble ...]` from anywhere (it locates the
worktree itself, exports SDKROOT for its build step, and `DRYRUN=1` prints every command).
Resume logic: the committed trajectory counter is read back with `rmeas -ckptInfo` (validated),
`rhmc -ntraj:(target − counter)` runs only the remainder (warmup is an absolute trajectory
index in rhmc, so it is never redone), and completed ensembles are skipped; rmeas per-config
TSVs and the analysis are idempotent.  Parameters: a_t = 0.2, M = 1, orders 31/11, window
[0.3, 12.5], Hasenbusch @[m_sea, 0.5], τ = 1, steps 4 × innerSteps 5 — exactly WP-H's demo.
**Exception: L2g30m00 runs the widened window [0.15, 14.0]** — σ_min at g²R = 3 is unknown
until `windowCheck` measures it and the quenched L=2 g²R=3 probe gave min|D_W−1| = 0.49;
maxRelErr(11) on [0.15, 14] is ≈ 4.5e-4, an MD-force error only (Metropolis at order 31
corrects it, at some acceptance cost).  If the monitor still hard-stops, widen further and
record it here.

Launch commands (one per line, copy-pasteable; NO backgrounding — the coordinator owns that):

```
bash src/experimental/radial/campaign/t2.sh L1g15m00
bash src/experimental/radial/campaign/t2.sh L1g15m01
bash src/experimental/radial/campaign/t2.sh L1g15m02
bash src/experimental/radial/campaign/t2.sh L1g15m03
bash src/experimental/radial/campaign/t2.sh L1g15m04
bash src/experimental/radial/campaign/t2.sh L1g05m00
bash src/experimental/radial/campaign/t2.sh L1g10m00
bash src/experimental/radial/campaign/t2.sh L1g10nf4
bash src/experimental/radial/campaign/t2.sh L1g10nf6
bash src/experimental/radial/campaign/t2.sh L2g30m00
bash src/experimental/radial/campaign/t2.sh pureL1 pureL2 pureL4
```

(or `bash src/experimental/radial/campaign/t2.sh` for everything in priority order; the
underlying rhmc/rmeas command for any ensemble is what `DRYRUN=1` prints.)

| ensemble | parameters | trajectories (warmup) | measured configs (every) | HMC wall-clock | measurement | targets |
|---|---|---|---|---|---|---|
| L1g15m00 | L=1 g²R=1.5 N_f=2 m=0, seed 1001 | 170 (10) | 32 (5) | ~2.0 h @ 41.5 s | ~1.4 h (currents+scalars+gluon+wspec, disc) | T2.1 T2.4 T2.5 T2.6 T2.7 T2.8 T2.9 |
| L1g15m01–04 | m = 0.1/0.2/0.3/0.4, seeds 1002–1005 | 170 (10) each | 32 (5) | ≲2 h each (m>0 is cheaper) | ~15–45 min (cond only; m>0 solves are several× cheaper than the measured m=0 140 s) | T2.3 (published slopes 0.073/0.143/0.210/0.271) |
| L1g05m00 | g²R=0.5, seed 1006 | 170 (10) | 32 (5) | ~2 h | ~1.4 h | g²R trend |
| L1g10m00 | g²R=1.0, seed 1007 | 170 (10) | 32 (5) | ~2 h | ~1.4 h | g²R trend |
| L1g10nf4 | N_f=4, seed 1008 | 100 (10) | 18 (5) | ~2.2 h @ est. 79 s/traj (2 pf copies) | ~45 min | N_f trend |
| L1g10nf6 | N_f=6, seed 1009 | 70 (10) | 12 (5) | ~2.3 h @ est. 116 s/traj (3 copies) | ~30 min | N_f trend |
| L2g30m00 | L=2 g²R=3.0, window [0.15,14], seed 1010 | 88 (8) | 20 (4) | ~5–6 h @ WP-H's est. 200–250 s/traj | ~2.5 h + wspec 80 s × 4 configs (tstride 5) | T2.6 "~3% at L=2", (a/R)² direction |
| pureL1/2/4 | heatbath, g²a=1, seeds 2001–2003 | — | 256/256/128 configs | — | ~5/10/5 min (gluon only) | T2.2 free baselines |

Whole campaign ≈ 25–30 h serial.  **Honest statistics statement**: 32 configurations spaced 5
trajectories is O(20–30) effective samples (τ_int of these observables is unknown until `jack`
measures it — every summary carries it).  Expect: the condensate slope points at the few-% level
(noise per config was ±0.4% at 24 vectors in WP-I's test; gauge scatter will dominate); Δ_A and
Δ_PS/Δ_FS at the 5–15% level from t ≈ 0.6–1.6 effective masses; the current-ratio panels
(Δ_l2/Δ_l1, Δ_l3/Δ_l1) and gluonic ratios at the 10–30% level — trend-level reproductions of
slides 11–16, not the deck's error bands (they had cluster-scale statistics); the N_f=6 point is
12 configs and will be barely more than a direction.  The l=3 per-m spread and the PS==FS/Ward
identities are exact-per-configuration statements and need no statistics.

### What did not work / caveats

1. **The (V.4) arccosh effective mass is unusable on MC correlators** (reference point below
   noise) — replaced by the local log ratio in rmeas's analysis only (decision 1 above).
2. **The stochastic axial correlator has a steep excited-state head**: on the single thermalized
   n_t=60 config, Δ_eff(l=1) runs 7.4 → 5.8 → 4.9 over t = 0…0.4 before the tower settles, and
   the one-config noise floor (npair 8) swallows the signal past t ≈ 1.2.  The campaign
   statistics push that to t ≈ 1.6–2; `-npair` is the knob if Δ_A needs more (cost is linear).
3. **Δ_V at smoke statistics is nonsense** (−2.8 ± 4.8) — the disconnected piece is exactly as
   noisy as slide 12 complains; T2.5 is the weakest deliverable at this scale, as doc/03 already
   ranked it.
4. **The l = 0 connected current correlator is NOT flat at n_t=16** (`l0_conn_variation` = 9.6):
   the whole smoke lattice sits inside the overlap kernel's locality range, exactly WP-I's dense
   finding.  The same summary row on the n_t = 60 ensembles is the measurement WP-I asked for
   ("measure that decay before quoting vector-channel numbers at small separations").
5. **cond on m = 0 ensembles is a tautology** (exact zero) — the campaign only runs it on the
   m-scan ensembles.
6. **wspec at L = 4/n_t = 60 (dim 19440) is out of reach** for dense zgeev on this machine —
   skipped with a message (T2.1's L=4 row is not covered; no L=4 ensemble is planned anyway).
   At L = 2 it costs 80 s/config, so the campaign subsamples 4 configs (the deck used 3).
7. **Quenched configurations are much worse conditioned than dynamical ones** at the same
   coupling (measured above): fermionic measurements on pure ensembles would be slow AND sit
   outside rhmc's window guarantees — rmeas warns via its own `kernelWindow` check.
8. The dt≠0 PS==FS "check" is arithmetically trivial in the point-source contraction — found
   while validating; the dt=0 contact is the reported check (decision in the scalars bullet).
9. `campaign/free/` (WP-K's figure scripts) shares the campaign directory; t2.sh does not touch
   it.
10. Nothing in `meas/`, `ops/`, `hmc/` was modified: rmeas is a pure consumer.  The only doc/04
    change is the §1 layout line (rmeas.nim, campaign/).

---

## Preliminary-statistics profile  (main, 2026-08-21 ~17:30, per user decision)

The user called the statistics question: this pass is a preliminary signal-over-noise test;
production statistics move to a real cluster. Actions:

* The three long ensembles were **stopped at their checkpoints**: L1g15m00 @ 80 traj (17 cfgs),
  L1g15m01 @ 80 (17 cfgs), L1g10m00 @ 100 (21 cfgs). `campaign/t2.sh` targets were pinned to
  those counters so its resume logic skips HMC and proceeds straight to measurement + analysis.
* Fresh short essentials (targets edited in t2.sh, annotated in-file): condensate points
  L1g15m02/03/04 @ 60 traj (~10 cfgs each — the condensate is a volume-averaged one-point
  function with small per-config noise, so this suffices for the slide-10 linearity signal),
  L1g05m00 @ 60 (the weak end of the g²R trend), **L2g30m00 @ 48** (pipeline demonstration +
  one (a/R)² point; the deck's 3 % ℓ=3-splitting claim at L=2 is NOT resolvable at this size —
  cluster job).
* **Dropped for the cluster: L1g10nf4, L1g10nf6.** At ~10–20 configs the N_f separation of
  slides 12/14/15/16 (e.g. Δ_F/Δ_A moving 0.78 → 0.80 between N_f = 2 and 6 at g²R = 1) is far
  below the noise; running them now would only demonstrate machinery that thmc already tests.
* Relaunched as four background tracks: M (measure+analyze the stopped trio), N1 (m02, m03),
  N2 (m04, g05), P (L=2). Expected completion: ~3–5 h wall.

Statistics expectation, stated honestly: per-ensemble 10–21 measured configs → condensate points
at the few-% level (good signal); correlator-ratio observables (Δ ratios of slides 11–16) at
trend level with 10–30 % errors; ℓ=1,2 protection and PS≡FS are exact per configuration and
independent of statistics.

---

## WP-M — standard-overlap mass migration (2026-08-26)

The active convention now matches established lattice QCD:
\[
D(m)=\left(1-\frac{m}{2\rho}\right)D(0)+m,\qquad \rho=1,\qquad 0\le m<2.
\]
The value \(\rho=1\) follows from this project's normalization \(D(0)=1+V\); it is
distinct from the Wilson-kernel height M. The exact map from the retired additive parameter is
\[
\mu=\frac{m}{1-m/2},\qquad m=\frac{\mu}{1+\mu/2},\qquad
D_{\rm std}(m)=(1-m/2)D_{\rm add}(\mu).
\]
Equal numeric masses therefore are not the same finite-cutoff theory, and no checkpoint or
measurement is converted in place.

The implementation treats the convention as one scientific contract:

1. Forward/adjoint operators, normal solves, propagators, and dense inverses use standard
   \(D(m)\).
2. Gauge derivatives, dense tangents, Hasenbusch forces, and conserved currents carry
   \(\alpha(m)=1-m/2\).
3. The condensate is
   \({\rm Re\,tr}[(1-D_{\rm ov}/2)D(m)^{-1}]/N\).
4. Massive FS contractions use
   \[
   (1-D_{\rm ov}^\dagger)D(m)^{-\dagger}
   =\frac{(1+m/2)D(m)^{-\dagger}-1}{1-m/2}.
   \]
5. Always-on validation rejects non-finite masses and values outside \(0\le m<2\).
6. Checkpoint version 2 stores a mass-convention id and rejects version 1. TSVs store
   massConvention=standard-overlap-rho1 and overlapRho=1; skip/analysis paths validate them.
7. The campaign defaults to output/radial/t2-standard-overlap. Any frozen rational-window or
   manifest change requires another fresh RADIAL_T2_OUT. Legacy output/radial/t2 is preserved.

Production-flag verification:

* toverlap passes: the direct formula and additive-parameter map agree at
  \(2.1\times10^{-16}\); massive pullback scaling agrees to \(1.3\times10^{-15}\).
* tmeas passes: the finite-mass FS identity is accurate to \(4.7\times10^{-15}\), and the
  point-source FS contraction agrees with its dense reference within \(5.1\times10^{-13}\).
* thmc passes: dense pseudofermion actions agree within \(2.4\times10^{-14}\), massive
  frame-force finite differences within \(8.2\times10^{-9}\), and checkpoint-v2 restart is
  bitwise identical.

The previous Tier-2 report is historical additive-convention output. Tier-1 and strictly
massless operator identities are unchanged; supported production starts new Tier-2 ensembles.
