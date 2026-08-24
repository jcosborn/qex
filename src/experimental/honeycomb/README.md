# QCD on the 16-cell honeycomb

An implementation of lattice gauge theory and Wilson-type fermions on the four
dimensional **16-cell honeycomb** ({3,3,4,3}, vertex set = the D₄ lattice) inside QEX,
reproducing the Lattice 2026 talk

> S. D. Katz and D. Nogradi, *Lattice QCD on the 16-cell honeycomb*,
> Lattice 2026; companion paper arXiv:2512.10604 [hep-lat].

Instead of the hypercubic lattice, sites are `Z⁴ ∪ (Z+½)⁴` (4D body-centred cubic).
Each site has 24 nearest neighbours and 12 links, the shortest closed loop is an
equilateral triangle, and the point symmetry group has 1152 elements instead of 384.
The claim under test is that the larger symmetry gives much smaller cut-off effects.

## Documents

| file | contents |
|---|---|
| [doc/REPRODUCTION.md](doc/REPRODUCTION.md) | **final report**: slide-by-slide verdict, validation chain, findings |
| [doc/SLIDES.md](doc/SLIDES.md) | slide-by-slide record of the talk and what the authors did |
| [doc/FORMULATION.md](doc/FORMULATION.md) | **normative**: lattice, link indexing, action, clover, Dirac operator, all conventions and verified numbers |
| [doc/PLAN.md](doc/PLAN.md) | work breakdown, module interfaces, acceptance tests |
| [doc/STATUS.md](doc/STATUS.md) | per-task logs: APIs, commands, all measured numbers |
| [doc/RESULTS.md](doc/RESULTS.md) | honeycomb χ_top production runs (slide 10 / Figs 1–2) |
| [doc/RESULTS_CUBIC.md](doc/RESULTS_CUBIC.md) | cubic reference pipeline results |
| doc/RESULTS_FERMIONS.md | quenched Wilson–Dirac spectra (slides 16/18/19 / Figs 5–8) |

## Modules

| file | contents |
|---|---|
| `hcgeom.nim` | geometry: 24 directions, 32 apex triangles, 16 hexagons, point group |
| `hcfree.nim` | free-field momentum-space Dirac operator (spectrum, pressure) |
| `hclayout.nim` | cell `Layout` + 2-site basis + the 16-way diagonal shift (`HcShift16`) |
| `hcgauge.nim` | 24 gauge link fields (`HcGauge`), triangle loops, observables |
| `hcaction.nim` | triangle action, 8-staple sums, force (`hcAction`, `hcForce`) |
| `hchmc.nim`, `hcheatbath.nim` | HMC (leapfrog/2MN/4MN5FV) and Cabibbo–Marinari KP heatbath + OR |
| `hcio.nim` | configuration save/load (single LIME file, 24 fields) |
| `hcflow.nim` | Lüscher flow for the triangle action, RK3, calibrated `cflow = 6` |
| `hctopo.nim` | hexagon-clover F̂_μν, energy density, topological charge |
| `hcwilson.nim` | interacting Wilson-type Dirac operator on all 24 directions |
| `hcstout.nim`, `hcclover.nim` | stout smearing + tree-level clover term (honeycomb *and* cubic) |
| `hcarnoldi.nim` | Krylov–Schur Arnoldi (non-Hermitian, LAPACK zgeev/zgees) |
| `hcanalysis.nim` | t₀ finder, jackknife, autocorrelation, continuum fits |

Executables: `hcFreeSpectrum`, `hcFreePressure` (slides 12–14), `hcPureGauge`
(generation), `hcMeasFlow` (flow + t₀ + Q), `hcSpectrum`/`cubicSpectrum` (fermion
spectra), `refCubicGen`/`refCubicMeas` (cubic reference).  Tests live in `tests/`
(one per module; all standalone, exit non-zero on failure).

## Building

```bash
export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
cd build_mac
make run src/experimental/honeycomb/hcgeom.nim
```

`build_mac/` is configured with the project's standard options
(`cctype:clang cc:clang-mp-22 cflagsspeed:"-Ofast -march=native -ffast-math"
qmpdir:... qiodir:... simd: vlen:4`).  `SDKROOT` is required or clang will not find the
system headers.
