# QCD on the 16-cell honeycomb

Lattice gauge theory and Wilson-type fermions on the four-dimensional
16-cell honeycomb {3,3,4,3} (vertex set D₄, equivalently the 4D bcc lattice
Z⁴ ∪ (Z+½)⁴) inside QEX, after

> S. D. Katz and D. Nogradi, *QCD on the 16-cell honeycomb*, arXiv:2512.10604
> [hep-lat]; talk *Lattice QCD on the 16-cell honeycomb*, Lattice 2026.

Every site has 24 nearest neighbours and 12 links, the shortest closed loop is
an equilateral triangle, and the point group has 1152 elements (Weyl group of
F₄) against 384 for the hypercubic lattice.  The claim under test is that the
larger symmetry gives smaller cut-off effects.

| document | contents |
|---|---|
| [doc/FORMULATION.md](doc/FORMULATION.md) | normative: lattice, link indexing, action, flow, clover, Dirac operator, free-fermion results, all conventions |
| [doc/RESULTS.md](doc/RESULTS.md) | what was reproduced, the validation chain, the numbers, and what remains |

Raw data, logs, plots, scripts and the session documents of the first campaign
are archived outside the tree (`archive/` in the worktree, see RESULTS.md).

## Modules

| file | contents |
|---|---|
| `honeycomb.nim` | barrel |
| `hcgeom.nim` | geometry: 24 directions, 32 apex triangles, 16 hexagons, point group |
| `hcgauge.nim` | `HcGauge` (24 link fields per cell), 16-way shift tree, triangle loops, gauge transformations, configuration files |
| `hcaction.nim` | triangle action, 8-staple sums, force (`ActionWork`) |
| `hcflow.nim` | gradient flow (RK3, `cflow = 6`) and the stout step (`stoutKappa = 1/3`) |
| `hctopo.nim` | hexagon-clover F̂_μν, energy density, topological charge (`TopoWork`) |
| `hcwilson.nim` | Wilson-Dirac operator on all 24 directions, optional clover term, antiperiodic time |
| `hchmc.nim`, `hcheatbath.nim` | HMC (mdevolve integrators, `hmc/metropolis`) and Cabibbo–Marinari heatbath + overrelaxation |
| `hcarnoldi.nim` | Krylov–Schur Arnoldi for non-Hermitian operators (LAPACK zgeev/zgees/ztrsen) |
| `hcspec.nim` | shift-invert operator, chirality, real modes, Q_Dirac bookkeeping |
| `hcanalysis.nim` | t₀/w₀ finders, autocorrelation, jackknife (`utils/resample`), polynomial fits, flow bookkeeping |
| `hcfree.nim` | free momentum-space operator on both lattices (spectrum, pressure) |
| `cubic.nim` | cubic-lattice counterparts: repeated `StoutSmear`, clover Wilson operator |

Executables: `hcpuregauge` (generation), `refcubicgen` (cubic HMC generation),
`measflow` (flow, t₀, Q, χ_top; `-lattice:hc|cubic`), `spectrum` (low modes,
chirality, Q_Dirac vs Q_flow; `-lattice:hc|cubic`), `freespectrum` and
`freepressure` (free-fermion slides).  Tests in `tests/`, one per module.

## Build and test

```bash
export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk   # macOS clang headers
cd build_mac
make src/experimental/honeycomb/tests/tgeom.nim && OMP_NUM_THREADS=4 ./bin/tgeom
```

All suites: `for t in tgeom tgauge taction thmc tflow ttopo twilson tstout
tclover tarnoldi tanalysis tfree tspectrum; do make
src/experimental/honeycomb/tests/$t.nim && OMP_NUM_THREADS=4 ./bin/$t; done`.
QEX thread barriers spin; keep `OMP_NUM_THREADS` at or below the free cores.
The default SIMD layout (`VLEN = 4`) accepts cell sizes 4, 8, 12, 16, 20, 24;
other sizes need `-simdlen:1`.

## Conventions in one place

`a = 1` is the cubic-sublattice spacing (the nearest-neighbour distance);
`μ = 3` is time.  Sites per cell 2, volume per site `a⁴/2`, so extensive
gluonic observables carry a factor ½ per site and χ_top uses the cell count as
the volume.  `β = 2N/g²` with the action `(β/2) Σ_x Σ_{32}(1 − Re Tr P/N)`.
Flow time is continuum-normalised (`cflow = 6`, analytic).  The stout step is
`ρ/3` of flow time, so equal smearing radii need `ρ_hc = 3 ρ_cubic`.  The
Wilson term is `a r p²/2`; the clover coefficient is `c_SW = 1` at tree level
on both lattices.  QEX `dot` conjugates its first argument.
