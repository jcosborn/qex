# HANDOFF REPORT — QED3 radial-quantization reproduction

Written 2026-08-21 ~17:50 by the session that built this project, for a successor model.
**Read this file first, then `06-status.md` for details. Everything you need is on disk.**

---

## 2026-08-26 convention migration — read before any Tier-2 command

The active code now uses the standard lattice-QCD overlap mass
\[
D(m)=(1-m/2)D_{\rm ov}+m\qquad(\rho=1),
\]
not the legacy additive \(D_{\rm ov}+\mu\). Binary checkpoints are version 2 and carry an
explicit mass-convention id; all version-1 checkpoints are intentionally rejected. Measurement
TSVs carry `massConvention=standard-overlap-rho1` and `overlapRho=1`, and restart/analysis
rejects files without matching metadata.

Everything below that refers to `output/radial/t2`, its running jobs, results, or resume commands
is a **historical 2026-08-21 record**. Do not resume those checkpoints or append measurements
with the current binaries. They remain useful only as legacy additive-convention artifacts.
The active campaign defaults to:

```
output/radial/t2-standard-overlap/
```

Start it with `bash src/experimental/radial/campaign/t2.sh <ensemble>`. To change a frozen
rational window, mass ladder, or other ensemble manifest field, set `RADIAL_T2_OUT` to another
fresh directory; never alter parameters and resume an existing checkpoint. The exact parameter
map is \(\mu=m/(1-m/2)\), \(m=\mu/(1+\mu/2)\), but no automatic data migration is performed.

---

## 0. What this project is

Reproduce, inside QEX, the Lattice 2026 talk *"Studying QED3 in radial quantization:
Interacting system on coarse lattices"* (N. Matsumoto) and its companion paper
**arXiv:2510.03085** (free limit). Lattice = refined icosahedron projected on S² (N_V = 10L²+2)
× regular time direction; non-compact U(1) Gaussian gauge action; two-component **overlap**
fermion (Zolotarev order 31 action / 11 force); HMC with Hasenbusch; radial-quantization
spectroscopy (operator dimensions Δ from temporal correlators).

* Worktree (ALL work happens here):
  `/Users/xjin/K/W/P003/qex/.claude/worktrees/qed3-slides-reproduction-plan-0e70b6`
* All new code: `src/experimental/radial/` (that directory of the worktree).
* Docs, in reading order: `doc/01-slides.md` (what the talk shows), `doc/02-formulation.md`
  (every equation, normative), `doc/03-targets.md` (every published number + tolerance),
  `doc/04-interfaces.md` (module APIs, normative), `doc/05-plan.md` (work packages),
  `doc/06-status.md` (**the audit trail — every measured number lives here**),
  `doc/07-observables.md` (derived observable definitions), this file.

### Absolute safety rules (violating these destroys the user's uncommitted work)
1. NEVER write to `/Users/xjin/K/W/P003/qex/src/...` or anything in the main checkout.
   The main checkout contains ~5000 lines of the user's UNTRACKED prior work
   (`src/experimental/qed3/`, several `docs/qed3_*.md`). Work only in the worktree above.
2. NEVER run state-changing git commands (no checkout/reset/clean/stash/commit/add) unless the
   user explicitly asks. Current branch: `claude/qed3-slides-reproduction-plan-0e70b6`.
   Everything under `src/experimental/radial/` and `output/` is intentionally uncommitted.
3. Subagents you spawn must be given rules 1–2 verbatim. A past project lost work this way.

### Build (memorize this)
```bash
cd /Users/xjin/K/W/P003/qex/.claude/worktrees/qed3-slides-reproduction-plan-0e70b6/build_mac
SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk make run experimental/radial/tests/tgeom
```
* `SDKROOT` is MANDATORY (clang-mp-22 cannot find system headers without it; `$(xcrun ...)` may
  fail under the sandbox — use the literal path).
* `make <suffix-of-path>` finds targets case-insensitively; binaries land in `build_mac/bin/`.
* `make` can exit 0 even when the build failed — check the binary exists / rerun and grep Error.
* Every app calls `freezeTimers()` right after `qexInit()` — do not remove (QEX tic/toc metadata
  allocation once blew RSS up 900× on these small volumes).

---

## 1. State: what is DONE and verified (do not redo)

**Tier 1 — the free-limit paper is fully reproduced.** 12 test suites, ~196 tests, all green
(`tests/t{geom,zolo,analytic,spinor,solve,wilson,overlap,gauge,flow,hmc,meas,fit}.nim`).
Master table (published | ours):

| quantity | published | ours |
|---|---|---|
| fermion Δ₀ (L=1, L_t=168, T=16) | 0.953918 | 0.953918 (−1.4e−7) |
| fermion Δ₀^cont | 0.999998(34) | 0.999999(44) |
| gauge Δ₀ (L=1, L_t=120, T=16) | 1.33242 | 1.332430 |
| gauge Δ₀^cont | 1.41409(18) | 1.414208(70) |
| Fig 10 T-fold: Wilson / overlap | visible / protected | 0.87 / 8.8e−14 |

Figures 4–12 of the paper: TSVs + rendered PNGs in `output/radial/free/`, gnuplot scripts in
`src/experimental/radial/campaign/free/`. Apps: `rfree.nim` (correlators + fits),
`rspec.nim` (spectra). Logs: `build_mac/rfree_full2.log`, `rspec_full2.log`.

**One known, honest discrepancy (do NOT try to fix by tuning):** the Eq. (V.9) n_max integers.
Ours: fermion 3/7/14/28, gauge 3/7/16/30. Published: 6/10/19/32 and 3/8/18/35. BUT the published
residual column equals our relative residual/DOF evaluated AT their n_max — the correlators agree
perfectly; only the paper's (unstated) n_max selection rule differs. Recorded in 06-status WP-K.

**The single most important physics fact discovered here** (doc/06 section "THE COUPLING
CONVENTION"): the paper's couplings use the **exact spherical kite area**
A_e = Σ± 4·arctan(tan(ℓ/4)·tan(ℓ*±/2)) for both κ_e and β_ℓ — NOT the flat ½ℓ(ℓ*₁+ℓ*₂) that
Eq. (IV.2) writes as an equivalent. Both published Δ₀ values pin this to six digits. `Edge.area`
in `core/geom.nim` IS the exact form. Exception: the slide-8 Wilson-spectrum legends match the
FLAT convention (`gcGeodesic` in `ops/gaugeact.nim`), and the L=1 legend additionally used
a_t ≈ 0.1333 instead of the campaign's 0.2. All reproduced; see `output/radial/free/slide8_legend.tsv`.

**Tier-2 exact results already in hand** (statistics-independent, from tests + smoke runs):
ℓ=1,2 correlator multiplets exactly degenerate under I_h (machine precision — slide 13's
"protection"); ℓ=3 splits exactly 3+4 (splitting 27.4% at L=1, 6.9% at L=2 free-field);
σ_PS ≡ σ_FS exactly at every dt (stronger than the slide's "identical spectra"); Ward charge
conservation ~1e−9 on dynamical configs; m=0 condensate exactly 0 by Ginsparg–Wilson.

---

## 2. Historical 2026-08-21 state (legacy additive campaign; do not resume)

At handoff time (~17:50, 2026-08-21) four background shell jobs were running the preliminary
Tier-2 campaign. **They may or may not have survived the session handoff** (they are OS
processes started by the previous session). CHECK FIRST:

```bash
cd /Users/xjin/K/W/P003/qex/.claude/worktrees/qed3-slides-reproduction-plan-0e70b6
ls -lT output/radial/t2/*/ckpt 2>/dev/null       # ckpt mtimes advancing => still running
tail -2 output/radial/t2/*/hmc.log 2>/dev/null    # if the script logs there; else check dirs
```
These commands are retained only to document the old run. They must not be used with the
current standard-overlap binaries:
```bash
cd /Users/xjin/K/W/P003/qex/.claude/worktrees/qed3-slides-reproduction-plan-0e70b6
export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
bash src/experimental/radial/campaign/t2.sh <ensemble>
```

The preliminary campaign (user decision, 17:30: low statistics, signal-over-noise; production
runs move to a real cluster later — targets already edited into `t2.sh` with annotations):

| ensemble | target traj | state at handoff | measures |
|---|---|---|---|
| pureL1, pureL2 (pure gauge, exact heatbath) | 256 cfgs | **done** | flow scan (slide 9) |
| L1g15m00 (L=1, N_f=2, g²R=1.5, m=0) | 80 (pinned at stop) | HMC done, **measuring** | currents, scalars, gluonic, wspec (slides 11–16, T2.1) |
| L1g15m01 (m=0.1) | 80 (pinned) | HMC done, meas queued | condensate (slide 10) |
| L1g10m00 (g²R=1.0, m=0) | 100 (pinned) | HMC done, meas queued | full set |
| L1g15m02 / 03 / 04 (m=0.2/0.3/0.4) | 60 each | HMC running (~traj 6+) | condensate |
| L1g05m00 (g²R=0.5, m=0) | 60 | queued behind m04 | full set |
| L2g30m00 (L=2, g²R=3.0, m=0) | 48 | HMC just started | full set (subsampled wspec) |
| L1g10nf4, L1g10nf6 | — | **deliberately skipped** | N_f trend is below laptop noise; cluster job |

Expected completion ≈ 21:00–23:00 local. Each ensemble dir gets
`output/radial/t2/<ens>/meas/*.tsv` (per-config) and `analysis/summary.tsv` (jackknifed).

The four job loops were: M = {L1g15m00, L1g15m01, L1g10m00} (measure+analyze only);
N1 = {L1g15m02, L1g15m03}; N2 = {L1g15m04, L1g05m00}; P = {L2g30m00}. Rerunning
`t2.sh <ens>` for each name in that order reproduces them exactly; running them one at a time
is ~2× faster per trajectory than 4 concurrent (CPU contention).

---

## 3. What remains — the finish line

### Step 1 — babysit / finish the campaign
When all ensembles above have `analysis/summary.tsv`, the data collection is complete.
If an ensemble aborts with `window: ... inside false`, the kernel window monitor fired: widen
`-ratLo`/`-ratHi` for that ensemble in `t2.sh` (documented in-file; L2 already uses [0.15, 14])
and rerun. Do NOT rebuild rationals mid-ensemble any other way.

### Step 2 — aggregate into the slide-by-slide preliminary report
Write (or extend `ranalyze.nim` if it exists — check) a short aggregation that reads every
`output/radial/t2/*/analysis/summary.tsv` and produces `output/radial/t2/report/`:
1. **Slide 10**: ⟨σ_PS⟩ vs mR from L1g15m01..04 (published points 0.073/0.143/0.210/0.271 at
   m=0.1..0.4 — ours are lower-statistics; the SIGNAL is linearity through the origin = no SSB).
2. **Slide 11/13 analogues**: axial-current Δ_{ℓ=2}/Δ_{ℓ=1} (CFT 3/2) and Δ_{ℓ=3}/Δ_{ℓ=1}
   (CFT 2) from the m00 ensembles; per-m ℓ=3 splitting.
3. **Slide 12**: Δ_V/Δ_A (CFT 1) — disconnected piece is noisy; report with honest errors.
4. **Slide 14/15**: gluonic Δ_F/Δ_A (free 1/√2 → CFT 1) and Δ_{F²}/Δ_F (CFT 2) from the
   GEVP summaries vs g²R ∈ {0.5, 1.0, 1.5} (+3.0 at L=2). NOTE (recorded in 06-status): the free
   reference for Δ_{F,ℓ=2}/Δ_{F,ℓ=1} is **√3**, not the slide's √(3/2) (slide used Δ_A in the
   denominator); at L=1 all spatial loop shapes collapse to ONE operator after ℓ-projection, so
   the L=1 GEVP is 3 temporal shapes with rank truncation.
5. **Slide 16**: Δ_PS/Δ_A, Δ_FS/Δ_A (free 1; published 0.88–0.98) + the PS≡FS exactness.
6. **Slide 9**: E_s(t)√L vs r/t from pureL1/pureL2 + the dynamical m00 ensembles. Known open
   item: raw E_s collapses across L for pure gauge, E_s·√L does NOT — the slide's √L
   normalization is unexplained (06-status WP-G); plot both and say so.
7. **Slide 8**: `rmeas -obs:wspec` outputs vs the published legend (free column already fully
   explained — see slide8_legend.tsv; interacting min|D_W−1| decreases with g²R = the
   "additive mass shift").
Use `meas/fit.nim` (jackknife `jack`, `effMass`, `plateauFit`) — the analysis step of rmeas
already produces most of this; the aggregation is mostly collation + gnuplot.
State statistics honestly in every plot footer (10–21 configs; trend-level).

### Step 3 — final write-up
Append a "Tier-2 preliminary results" master section to `doc/06-status.md` mirroring the WP-K
master table: slide | published | ours (± jackknife) | verdict (signal seen / trend / needs
cluster). Update `README.md` with a results summary. The user then decides about committing.

### Cluster scaling (user's stated plan; prepare, don't run)
`t2.sh` targets back to 170+ (annotations show the original values), re-enable nf4/nf6,
add L=4 ensembles and the g²a-constant lines of doc/01 §3. Code is serial; parallelism comes
from running ensembles concurrently. MPI is wired in QEX but radial code is single-rank
(checkpoint refuses nRanks≠1) — that is fine for ensemble-parallel running.

---

## 4. Gotchas that will bite you (each cost hours; all are in 06-status)

1. **σ² units**: Zolotarev poles live in σ² (spec(X†X) ⊂ [σmin², σmax²]). Windows are FROZEN
   per ensemble; the monitor hard-stops on violation. Never rebuild mid-ensemble.
2. **Solver floors**: keep `r2` targets ≥ 2 decades above (eps·cond)² or good solves report
   failure (`stats.ok=false, mrefits=0` is the marginal-floor signature; the fix is one decade
   of headroom, never a physics tolerance).
3. **The overlap kernel X = D_W − M uses the RAW operator, plain adjoints, M=1** (paper
   convention). The "hat" volume-normalized kernel is diagnostics-only.
4. **Generalized eigenvalues** use weight diag(volbar/volw) — doc/02 §3.2 as corrected (the
   paper's (IV.12) as printed is inverted vs its own (IV.11)).
5. `ovGradient` is THE only force/current kernel. Never write a second one — the Ward test
   only protects one.
6. Zero-mode projection hits the committed field, refreshed momentum, AND every MD force.
   ker M = gauge orbit + the uniform temporal (Polyakov) mode: dim = n_V·L_t (not n_V·L_t − 1).
7. Nim: never assign `result` inside `threads:`; `tFor` collides in one scope; the radial code
   is deliberately serial (OpenMP was 60× slower at these volumes in a prior attempt).
8. Only `import base` (+ explicit small modules), never `import qex`, in radial code.
9. LAPACK bindings live in `src/eigens/lapack.nim`; `zgesv` was added there (3 lines) — the
   only change outside `radial/`. `zeigs`/`zgeigs`/`zeigsgv` wrappers: `eigens/linalgFuncs`.
10. The effective-dimension convention is the paper's arccosh form (V.4–V.5) for deterministic
    data; rmeas's MC analysis uses a local log-ratio instead (arccosh is unusable on noise).
11. `rmeas` re-measurement is bit-reproducible (trajectory-addressed RNG); skip-done means you
    can always just rerun the whole `t2.sh <ens>` line.

## 5. Independent oracles (for re-verification)
Pure-Python cross-checks used to validate the Nim live in `/tmp/claude-502/qed3/*.py` (may be
gone — /tmp; the important results are all transcribed into 06-status). The strongest in-repo
checks: `tests/` (196 tests), the WP-K four-way validation (dense `denseOv`, dense `denseDw`,
real-space `regSolve`, WP-G pinned table), and `rgeom`/`rfree`/`rspec` which re-verify on
every run. The papers: slides + free-limit PDF summaries are fully transcribed in
`doc/01-slides.md` / `doc/02-formulation.md`; source PDFs were at
`https://indico.global/event/16565/contributions/161926/attachments/74534/144803/qed3_v10.pdf`
and arXiv:2510.03085.
