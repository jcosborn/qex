# The 16-cell honeycomb: formulation and conventions

Normative reference for `src/experimental/honeycomb`.  Every statement here is
either derived below or verified by a test named next to it.  If code and this
document disagree, one of them is wrong.

## 1. Lattice

D₄ = {n ∈ Z⁴ : Σn_μ even} and D₄* = Z⁴ ∪ (Z+½)⁴ are similar (ratio 1/√2).
The code uses D₄* with nearest-neighbour distance a = 1: a has the cubic
meaning and volumes are comparable with a cubic N_s³×N_t lattice.

Cell encoding: a hypercubic QEX layout of cells y ∈ Z⁴ with two sites per cell,
A(y) = y and B(y) = y + ½(1,1,1,1).  Sites 2N_cells, links 12 per site,
volume per site a⁴/2.  Periodic in the three (1,0,0,0)-type directions with
period N_s and (anti)periodic in (0,0,0,1) with period N_t; both map A→A and
B→B, so any N_s, N_t is allowed.

The 24 unit neighbour vectors:

    axis      ±ê_μ                       8   A↔A, B↔B
    diagonal  d(δ) = δ − ½(1,1,1,1)     16   A↔B,  δ ∈ {0,1}⁴,  δ = δ₀ + 2δ₁ + 4δ₂ + 8δ₃

    Σ_{i=1}^{24} n_iμ n_iν = 6 δ_μν,     Σ_i (p·n_i)⁴ = 3 (p²)²        (tgeom, tfree)

Link fields per cell (`HcGauge`): `uA[μ]` A(y)→A(y+ê_μ), `uB[μ]` B(y)→B(y+ê_μ),
`uD[δ]` B(y)→A(y+δ) with displacement d(δ); A(y)→B(y−δ) is `uD[δ]†` at cell
y−δ.  `links` holds the same 24 fields flat in that order.

Point groups: 1152 (Weyl group of F₄) and 384 (B₄), by exhaustive enumeration
of orthogonal maps of the neighbour set onto itself (tgeom).

## 2. Triangles and hexagons

n_i + n_j + n_k = 0 has 32 unordered solutions, each with exactly one axis edge
and two diagonal edges; every triangle is equilateral of side a, area √3/4.
Each triangle has a unique apex, the corner between its two diagonal edges,
which gives the bijection

    triangle  ⟷  (site x, diagonal δ, direction μ),   δ' = δ xor 2^μ,   δ_μ = 0

so `for x in sites: for (δ, μ) in apexTris` visits every triangle once:
32 per site, 64 per cell, 96 through each site, 8 containing each link.

Apex on B at cell z:  `uD[δ](z) · uA[μ](z+δ) · uD[δ'](z)†`.
Apex on A, re-based at z = apex − δ̄ (δ̄ = δ xor 15, δ̄' = δ' xor 15):
`uD[δ̄](z)† · uB[μ](z) · uD[δ̄'](z+ê_μ)`.

Hexagons: six coplanar neighbour vectors {±u, ±v, ±(u+v)} with u·v = −½.  Every
hexagon contains one axis pair ±ê_μ and two diagonal pairs ±d⁻, ±d⁺ with
d⁺ − d⁻ = ê_μ; 16 per site (4 axes × 4 sign patterns modulo overall sign).
Ring order in `hexTriPaths`: d⁻, d⁺, ê_μ, −d⁻, −d⁺, −ê_μ.  Unit 2-form
Ω^(h) = ê_μ ∧ f with f = (1/√3) Σ_{ν≠μ} σ_ν ê_ν, i.e. Ω_μν = σ_ν/√3.

## 3. Gauge action

    S = (β/2) Σ_x Σ_{i=1}^{32} (1 − Re Tr P_i(x)/N),     β = 2N/g₀²

Classical limit: for the apex triangle (d, μ) the flux is
Φ = ½ S_αβ F_αβ = d_μ Σ_{ν≠μ} d_ν F_μν, and with Σ_{16 d} d_ν d_ρ = 4 δ_νρ

    Σ_{i=1}^{32} Tr Φ_i² = ½ Σ_{μ≠ν} Tr F_μν²,     S → ¼ Σ_x Σ_μν Tr F² = ½ ∫ Σ_μν Tr F_μν F_μν

using Σ_x = 2∫d⁴x.  This is the cubic Wilson action at the same β
(taction test 6: S₁₆/S_cubic = 1 + 0.056 p² on a plane wave).

Force (`hcaction`): with P = U_l V_k† for the 8 triangles containing link l,

    D_l = (β/2N) Σ_k V_k,   f_l = projectTAH(U_l D_l†),   d/ds S(e^{sP}U)|₀ = Σ_l redot(P_l, f_l)

with redot(a,b) = Re tr(a†b); HMC uses p −= dt f, U ← exp(dt p) U and
T = ½ Σ_l redot(p_l, p_l) (taction tests 1, 3; thmc).

## 4. Field strength, E, Q

Hexagon clover at site x: C_h = Σ_{k=1}^{6} P_k, the six triangle loops of
hexagon h based at x in ring order.  Each loop encloses flux (√3/4) a² F_Ω with
F_Ω = ½ Ω_αβ F_αβ, so

    F̂_Ω^(h) = s (4/√3) TAH[C_h/6],       F̂_μν = (3/8) Σ_h Ω^(h)_μν F̂_Ω^(h)

since Σ_h Ω_μν F_Ω = (8/3) F_μν (four hexagons with axis μ give (4/3)F_μν, four
with axis ν give the same).  The whole thing collapses to weights
c_h,μν = s Ω_μν/(4√3) = ±1/12 on TAH[C_h].  The ring order is clockwise
with respect to Ω, so s = `cloverSign` = −1 and F̂_μν = +i a² F_μν T for
U = exp(+i∫A·dl).  E and Q are even in F̂; the sign is pinned by the
site-wise weak-field test (ttopo test 2) and matters only for the clover term
of the Dirac operator.

    E(x) = −½ Σ_μν Tr F̂_μν F̂_μν,               ⟨E⟩ = (1/N_sites) Σ_x E(x)        (intensive)
    q(x) = −(1/32π²) ε_μνρσ Tr F̂_μν F̂_ρσ = −(1/4π²)[tr F̂₀₁F̂₂₃ − tr F̂₀₂F̂₁₃ + tr F̂₀₃F̂₁₂]
    Q    = ½ Σ_x q(x)                                                          (a⁴/2 per site)

The reductions are those of QEX `densityE`/`topoQ`, so honeycomb and cubic t₀
and Q share one normalisation.  Atiyah–Singer check (ttopo test 3): the
constant-flux Cartan background T = diag(1,−1,0) with fluxes n₁, n₂ gives
Q = 2n₁n₂ with a pure 1/L⁴ clover artefact 4 sin(f/4)/f, f = 2πn/L²; the
cubic clover artefact is sin(f)/f, so the honeycomb artefact at size L equals
the cubic one at 2L (the triangle's projected area is ¼ of a plaquette's).

## 5. Gradient flow and stout smearing

Lüscher's flow V̇ = −g₀² ∂S_W V with S_W over oriented plaquettes gives
Ȧ_μ = D_νG_νμ with t in units of a².  On the honeycomb, linearising the eight
triangle fluxes of a link around its midpoint gives for axis and diagonal
links alike, in U(1)-angle language,

    −Σ_{8 tri ∋ l} ε_tri φ_tri = (1/3) n_l · ∂_ν F_ν·        (cubic: Σ_{6 plaq} → ∂_ν F_ν μ)

and Re tr(1 − P) ≈ ¼ (Φ^a)² for tr T^aT^b = −½ δ^ab, so the naive transcription
V̇ = −∂_l[Σ_tri Re tr(1−P)] V runs at 1/6 of Lüscher's rate.  Equivalently:
the least-squares projection Σ_{12 links/site} n_μ n_ν = 3δ_μν and the half
site volume give 3Ȧ = ½ ∂F.  Hence

    V̇_l = −6 ∂_l[Σ_tri Re tr(1 − P)] V_l,        `cflow = 6`  (exact; tflow test 2 measures 5.9996(3))

The sampled plane wave is the exact leading-order eigenmode; the O(p²)
artefact of the honeycomb flow rate is −0.0023 p² against −p²/12 on the cubic
lattice.  The code uses Z = −ε force(β = cflow N) in the RK3 scheme of
`gauge/wflow.nim` (Lüscher App. C):
W₁ = e^{Z₀/4}W₀, W₂ = e^{8/9 Z₁ − 17/36 Z₀}W₁, V(t+ε) = e^{3/4 Z₂ − 8/9 Z₁ + 17/36 Z₀}W₂.

t₀: t²⟨E⟩ = 0.3; w₀²: t d/dt(t²⟨E⟩) = 0.3.  Scale in the paper: √(8t₀) = 0.47 fm.

Stout (Morningstar–Peardon) with the raw 8-triangle staple sum Σ_l:

    U' = exp(−ρ projectTAH(U_l Σ_l†)) U_l,      θ' = θ + (ρ/3) n·∂F

so one step is a flow step of size ρ·stoutKappa, stoutKappa = 2/cflow = 1/3
(cubic: 1; tstout test 5).  Equal smearing radius √(8t) needs ρ_hc = 3ρ_cubic.
The paper states "6 steps, ρ = 0.05, the same parameters on both lattices"
without a normalisation; the runs below use both readings.

## 6. Fermions

Naive term D₀ = (1/6) Σ_{i=1}^{24} (γ·n_i) ∇_i with ∇_i f(x) = (U_i f(x+n_i) − f(x))/a;
in momentum space D₀(p) = (i/6a) Σ_i (γ·n_i) sin(a p·n_i) → iγ·p by Σn n = 6δ.

Wilson term: the paper's eq. (4), D = D₀ + a(r/6) Σ_{i=1}^{24} ∇*_i∇_i, is
stated to give a r p²/2 at small p; as written it gives a r p² with the sign of
−∇*∇ (the 24 terms count each ±n pair twice).  The a r p²/2 normalisation is
the one used, fixed independently by max Re λ = 16/3 of the free spectrum:

    W(p) = (r/6a) Σ_{i=1}^{24} (1 − cos a p·n_i) → a r p²/2,
    D ψ(x) = (m + 4r/a) ψ(x) + (1/6a) Σ_i (γ·n_i − r) U_i(x) ψ(x+n_i)

(γ·n_i)² = 1 for all 24 directions, so (γ·n_i − 1)/2 is minus a rank-2
projector as on the cubic lattice.  D† is the same hopping sum with the spin
factor of the opposite direction; γ₅Dγ₅ = D† is tested (twilson).  The Wilson
term is −(r a/2) D² in the continuum on both lattices, so the tree-level
clover coefficient is the standard one:

    D_c = D − (c_SW r/4) σ_μν F_μν = D − (c_SW r/2) Σ_{a>b} γ_a γ_b F̂_ab,   c_SW = 1

with σ_μν = (i/2)[γ_μ,γ_ν], [D_μ,D_ν] = iF_μν and F̂ the anti-Hermitian
clover field of §4 (tclover test 5 pins the coefficient on the constant-flux
background; tspectrum test 4a is the physical pin: with c_SW = 1 the
would-be zero modes of the flux background sit at |Re λ| < 0.05 while the
Wilson term alone puts them at (r/2)(|f₁|+|f₂|) ≈ 0.39).

Free spectrum (a = 1, r = 1): D(p) = M + iγ·K, eigenvalues M ± i|K|, each twice,

    M(p)   = 4r − (r/6)[2 Σ_μ cos p_μ + 16 Π_μ cos(p_μ/2)]
    K_μ(p) = (1/6)[2 sin p_μ + 8 sin(p_μ/2) Π_{ν≠μ} cos(p_ν/2)]
    small p:  K = p(1 − p²/12),  M = r p²/2 − r(p²)²/16,  D†D = p² + (r²/4 − 1/6)(p²)² + O(p⁶)

| | 16-cell | cubic |
|---|---|---|
| max Re λ | 16/3 at p = (π,π,π,π) and (0,0,0,2π) | 8 |
| max Im λ | 3^{1/4}(1+√3)/√6 = 1.467890 on p ∝ (1,1,1,1) and on-axis at cos(p₀/2) = (√3−1)/2 (related by F₄) | 2 |

Brillouin zone: reciprocal lattice 2π D₄ (covolume 2(2π)⁴), fundamental domain
[0,2π)³ × [0,4π), 2N_s³N_t momenta; p₃ = (2n+1)π/N_t antiperiodic.

Plane-wave convention: both sublattices are phased with the integer cell
coordinate (no half-site phase); the 8×8 (spin ⊗ sublattice) free operator
`freeD8` carries the diagonal hops with e^{±ik·δ} (tfree, twilson).

Free pressure O = (p(T) − p(0))/T⁴ per fermionic degree of freedom, continuum
7π²/720.  For r = 1 the time dependence of M² + K² at fixed spatial momentum is
a polynomial in cos(p₃/2) (linear in cos p₃ on the cubic lattice) whose roots
w_k = e^{E_k} give the exact thermal sum Σ_k 2 ln(1 + w_k^{−2N_t}) per ln det;
the 1/N_t series follows from the physical dispersion relation
E(p⃗) = |p⃗|(1 + a₂|p⃗|² + a₄|p⃗|⁴ + …) and Γ(n)η(n) integrals:

    O_cubic   = (7π²/720)[1 + (248/147)(π²/N_t²) + (635/147)(π⁴/N_t⁴) + …]      ⟨a₂⟩ = −4/15
    O_16-cell = (7π²/720)[1 + (127/980)(π⁴/N_t⁴) + (73/4158)(π⁶/N_t⁶) + …]      a₂ = 0, ⟨a₄⟩ = −1/210, ⟨a₆⟩ = −1/11340

The O(a²) term is absent on the honeycomb because D†D is O(4) invariant at
that order (paper eq. (6)), so E² = p⃗² + O(a⁴).  Values at N_t = 4:
3.828128 and 1.071914 (cubic, 16-cell); at N_t = 6: 2.115496 and 1.011240
(freepressure; independent derivation in RESULTS.md).

Q_Dirac = qDiracSign (n₊ − n₋) over real low modes (|Im λ| < ε) with
Re λ below a window `recut` that must not exceed the reach of the converged
set, σ + max|λ_i − σ|; `qDiracSign = −1` on the constant-flux background
(Q = +2n₁n₂, zero modes have chirality −sign(n₁n₂)) in the DeGrand–Rossi
convention (tspectrum).

## 7. Symmetry and cut-off effects

For free fermions the absence of O(a²) is a theorem (§6).  In the gauge sector
F₄ symmetry forbids only the O(4)-breaking dimension-6 operator
Σ_μν tr(D_μF_μν D_μF_μν); the O(4)-invariant on-shell operator survives, so a
vanishing O(a²) term in t₀²χ is an empirical fit statement, as the paper says,
not a consequence of the symmetry.

## 8. Cost

Per unit volume: 2 sites × 3 more neighbours = 6 times the cubic work.  With
HMC cost ∝ V^{3/2}, a factor 2 in lattice spacing at fixed physics buys 64.

## 9. Code conventions

μ = 0..3, μ = 3 time (QEX `setBC`, `densityE`).  Sublattice 0 = A, 1 = B.
γ matrices: QEX `gamma1..gamma4` (DeGrand–Rossi), γ₅ = diag(1,1,−1,−1).
TAH(M) = (M − M†)/2 − tr(M − M†)/2N = QEX `projectTAH`.  QEX `dot(x, y)`
conjugates x.  Extensive gluonic observables carry a⁴/2 per site; χ_top uses
V = N_cells a⁴.
