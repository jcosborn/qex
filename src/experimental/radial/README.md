# radial — QED3 in radial quantization on \(S^2\times\mathbb R\)

Reproduction of

* **arXiv:2510.03085v2**, Boyle, Brower, Fleming, Katz, Matsumoto, Misra,
  *Studying QED3 with radial quantization on the lattice: Free limit* (PRD, to appear), and
* the Lattice 2026 talk *Studying QED3 in radial quantization: Interacting system on coarse
  lattices* (N. Matsumoto).

The lattice is a refined icosahedron projected onto \(S^2\) (\(N_V=10L^2{+}2\)) times a regular
temporal direction. The gauge field is a **non-compact** \(U(1)\) one-cochain with a Gaussian
action; the fermion is a two-component **overlap** operator built from a simplicial Wilson kernel.
Finite masses use the standard lattice-QCD convention
\(D(m)=(1-m/2)D_{\rm ov}+m\) at \(\rho=1\). Legacy additive-mass checkpoints and
measurements are rejected rather than silently reused.

Because the lattice is simplicial, this code deliberately does **not** use QEX's hypercubic
`Layout`/`Field` machinery. It imports `base` for parameters, timers and comms, and reuses
`algorithms/rk` (gradient flow), `mdevolve` + `hmc/metropolis` (HMC), `eigens/lapack` and
`eigens/linalgFuncs` (dense spectra, GEVP), `rng/threefry4x64`, and `utils/resample`.

## Continuing this work

**Start with [`doc/08-handoff.md`](doc/08-handoff.md)** — the self-contained handoff report
(state, running jobs, exact resume commands, remaining steps, gotchas).

## Documents

| file | contents |
|---|---|
| [`doc/01-slides.md`](doc/01-slides.md) | slide-by-slide transcription of the talk and what the author did |
| [`doc/02-formulation.md`](doc/02-formulation.md) | **normative** formulation: geometry, spin connection, Wilson, overlap, gauge action, analytic correlators |
| [`doc/03-targets.md`](doc/03-targets.md) | every reproduction target with its published number and acceptance tolerance |
| [`doc/04-interfaces.md`](doc/04-interfaces.md) | **normative** module layout and proc signatures |
| [`doc/05-plan.md`](doc/05-plan.md) | work packages and dependency graph |
| [`doc/06-status.md`](doc/06-status.md) | living status board — every measured number, the audit trail |
| [`doc/07-observables.md`](doc/07-observables.md) | derived observable definitions (currents, scalars, GEVP basis) |
| [`doc/08-handoff.md`](doc/08-handoff.md) | **handoff report — read first when resuming** |
| [`doc/09-preliminary-report.md`](doc/09-preliminary-report.md) | **the preliminary results report** (Tier 1 complete; Tier 2 signal-level) |

## Build

```bash
cd build_mac && SDKROOT=$(xcrun --show-sdk-path) make run experimental/radial/tests/tgeom
```

`SDKROOT` is required — without it `clang-mp-22` cannot find the system headers.
Binaries land in `build_mac/bin/`.
