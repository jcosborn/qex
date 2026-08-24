# HANDOFF — 16-cell honeycomb reproduction (read this first)

Written 2026-08-21 ~17:50 CDT by the session lead, for a successor (human or model)
continuing this work.  Everything is on disk in this worktree; nothing lives only in
the previous session's memory.  **One background job set is still running** (§3) — that
is your first order of business.

---

## 0. Mission and current state in one paragraph

We reproduced the Lattice 2026 talk "Lattice QCD on the 16-cell honeycomb"
(Katz & Nogradi, arXiv:2512.10604) inside QEX, in `src/experimental/honeycomb/`.
Status: slides 4, 7–9, 12, 13, 14, 21 reproduced exactly; slide 10 (topological
susceptibility scaling, the headline) reproduced at preliminary statistics with the
O(a²) coefficient consistent with zero; slides 16/18/19 (quenched Wilson–Dirac spectra)
were **in flight** when this handoff was written — two measurement jobs were running
and an autonomous agent was driving them; it may or may not have finished by the time
you read this.  Read [REPRODUCTION.md](REPRODUCTION.md) for the verdict table, then
come back here.

## 1. Reading order

1. This file.
2. [REPRODUCTION.md](REPRODUCTION.md) — final report, slide-by-slide verdicts, findings.
3. [FORMULATION.md](FORMULATION.md) — normative math/conventions.  If code and this
   document disagree, one of them is wrong; investigate, don't paper over.
4. [STATUS.md](STATUS.md) — per-task logs (11 sections): exact APIs, commands, every
   measured number.  Search it before re-deriving anything.
5. [PLAN.md](PLAN.md) — original work breakdown; [SLIDES.md](SLIDES.md) — what the talk
   shows; [RESULTS.md](RESULTS.md) / [RESULTS_CUBIC.md](RESULTS_CUBIC.md) /
   `RESULTS_FERMIONS.md` (if present) — physics run writeups.

## 2. Environment and build (memorize these five things)

```bash
export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk   # MANDATORY
cd <worktree>/build_mac                                              # already configured
make src/experimental/honeycomb/<file>.nim   # builds to build_mac/bin/<file>
OMP_NUM_THREADS=4 ./bin/<file> [-key:val ...]
```

1. **`SDKROOT` is mandatory** — without it clang cannot find `<string.h>` and every
   build fails with hundreds of header errors.
2. **`OMP_NUM_THREADS=4` always.**  The box has 16 logical CPUs shared with other
   users; QEX thread barriers are spin waits and livelock at 16 threads (30× slowdown
   or apparent hangs).  Results are bit-identical for 1/2/4/8 threads.
3. **Lattice sizes**: the default `vlen:4` SIMD layout accepts L ∈ {4, 8, 12, 16, 20,
   24} per dimension; L = 6, 9, 10, 11, 14 need the executable run with `-simdlen:1`
   (hcPureGauge supports it) or a `newLayout(lat,1)` code path.
4. **Never allocate QEX fields / touch the GC inside a `threads:` block** — it
   segfaults intermittently (1 in ~8 runs).  Allocate everything up front; see the
   persistent work-object pattern (`HcActionWork`, `HcTopoWork`, `HcWilson`).
5. **Background shell jobs cannot be killed** in this sandbox (`kill`/`pkill` denied)
   and `ps` is blocked.  Detect activity by file mtimes.  BSD `find -newermt` does NOT
   understand `'-90 minutes'` relative syntax (it silently matches nothing) — compute
   an absolute timestamp: `find . -newermt "$(date -v-10M '+%Y-%m-%dT%H:%M:%S')"`.

Full test suite (all should pass; ~5 min total):
```bash
for t in tgeom tfree tanalysis tgauge taction thmc tflow ttopo twilson tstout tclover tarnoldi tspectrum; do
  make src/experimental/honeycomb/tests/$t.nim && OMP_NUM_THREADS=4 ./bin/$t || echo "FAIL $t"
done
```
(`tspectrum` exists only if the in-flight task D4 finished writing it — check.)

## 3. ⚠️ THE IN-FLIGHT WORK — task D4 (fermion spectra, slides 16/18/19)

> **UPDATE (18:05 CDT):** D4 was stopped on user request and the lead harvested the
> logs (66 + 8 honeycomb / 60 cubic configs) — results are in
> `RESULTS_FERMIONS.md`, REPRODUCTION.md §1 is fully filled in, and the harvest.py
> FLOWQ-join bug + its fix are recorded in STATUS.md ("Task D4 close-out").
> Nothing below is *required* anymore; it remains valid if you want to extend the
> statistics (the unkillable jobs may have appended more configs — just re-harvest).

An autonomous agent ("D4") was running two measurement jobs when this was written.
Its goal: reproduce paper Figs. 5–8 at preliminary statistics.  The user's directive:
**low statistics are fine; harvest as soon as the signal is visible; a real cluster run
comes later.**

### What was running (state at 17:45 CDT)

| job | log (append-mode) | progress | rate | target |
|---|---|---|---|---|
| honeycomb spectrum, β=7.22, 8⁴ cells | `doc/plots/fermions/data/hc722.out` (configs ≤45) **continued by** `hc722b.out` (relaunched, counts "n/76") | 50/76 | ~49 s/cfg | 76 configs |
| cubic spectrum, β=5.86, 8⁴ | `doc/plots/fermions/data/cub586.out` | 26 done | ~11.6 s/cfg | ~60 |
| cubic ensemble generation | `cubgen586.out` | **done** (112 configs) | — | — |

Setup: 6 stout steps ρ=0.05 both lattices, tree-level clover c_SW=1, bare m=0, r=1,
antiperiodic time BC, shift-invert Krylov–Schur for the low spectrum (pilot logs
`pilot_s050/s065/s080_*.log` record the shift choice), per-config gradient flow for
Q_flow.  Configs live under `$TMPDIR` and are disposable; **the logs contain all
measurements** and live in the repo.

### What to do, by scenario

* **D4 finished** (check: `doc/RESULTS_FERMIONS.md` exists, STATUS.md has a
  "Task D4" section, `doc/plots/fermions/*.png` exist): verify the four figures and
  numbers, fill the three ⏳ rows of REPRODUCTION.md §1 with the measured values, done.
* **D4 died mid-way** (logs stale, no RESULTS_FERMIONS.md): everything needed is in
  `doc/plots/fermions/`: `harvest.py` (parses the logs — read it for the exact line
  format), `figs.gnuplot`, `run_hc.sh`, `run_cubic.sh`, `pilot2.sh`.  Procedure:
  1. Let/leave the unkillable jobs be; pick a cutoff config number per log
     (≥50 honeycomb, ≥40 cubic is enough) and note it.
  2. `python3 harvest.py` (read its usage first) → per-config tables of
     eigenvalues, chiralities, λ₀, Q_Dirac, Q_flow.
  3. Make the four figures (gnuplot): (a) scatter Im λ vs Re λ, both lattices;
     (b) normalised histogram of Re λ₀ — *the headline*: expect the honeycomb
     distribution ~3–4× narrower and centred lower than cubic; (c) histogram of
     Q_Dirac − Q_flow with a [−0.5, 0.5] band — expect honeycomb visibly narrower;
     (d) histogram of real-mode chirality — expect honeycomb piled at ±1.
  4. Write `doc/RESULTS_FERMIONS.md` (setup table, n_configs used, the four
     comparisons with numbers, "preliminary statistics" caveat), append a STATUS.md
     section, update REPRODUCTION.md §1.
  5. If a signal is NOT visible at this sample size, write that down as the finding —
     do not silently run longer.
* **Logs unusable / want a clean restart**: regenerate with
  `run_hc.sh`/`run_cubic.sh` (they wrap `bin/hcSpectrum` / `bin/cubicSpectrum`;
  binaries exist in `build_mac/bin/`).  Matched-point parameters as above.  Costs:
  honeycomb ~50–80 s/cfg, cubic ~12 s/cfg at 4 threads.

### Physics cross-check values for the fermion runs
Honeycomb β=7.22 corresponds to t₀/a² ≈ 1.7–1.9 (task R's calibration:
β=7.15 → 1.359, β=7.20 → 1.653, growing with d ln t₀/dβ ≈ 5–6); cubic β=5.86
interpolates to t₀/a² ≈ 1.9 (task C: β=5.8 → 1.518, β=5.9 → 2.309).  Both ≈ the
paper's a ≈ 0.12 fm point (t₀/a² = 1.92 for √(8t₀) = 0.47 fm).

## 4. What is DONE and trustworthy (do not redo)

All in STATUS.md with numbers; headline gates:

| layer | test | key number |
|---|---|---|
| geometry `hcgeom` | `tgeom` 42 checks | point groups 1152/384; 32 triangles; 16 hexagons |
| gauge `hclayout`/`hcgauge` | `tgauge` 38 | shift16 bit-exact; gauge invariance 9e-19 |
| action/force `hcaction` | `taction` 26 | FD force 9e-12; β-normalisation vs Wilson: 1+0.056p² |
| MC `hchmc`/`hcheatbath`/`hcio` | `thmc` 28 | HMC↔heatbath agree 0.46σ; ⟨e^{−ΔH}⟩=1 |
| flow `hcflow` | `tflow` 9 | **cflow = 6 exact** (heat-kernel calibrated) |
| topology `hctopo` | `ttopo` 22 | Atiyah–Singer Q/2n₁n₂ → 1, 1/L⁴; **hcCloverSign=−1** |
| Dirac `hcwilson` | `twilson` 22 | entrywise vs momentum space 4e-15 |
| stout/clover `hcstout`/`hcclover` | `tstout`+`tclover` 46 | κ_hc = **1/3** (ρ_hc=3ρ_cubic for equal smearing); clover pinned 2e-14 |
| eigensolver `hcarnoldi` | `tarnoldi` 37 | vs exact Wilson spectra 1e-13 |
| analysis `hcanalysis` | `tanalysis` 23 | t₀ finder, jackknife, fits |
| cubic reference | RESULTS_CUBIC.md | 10⁴t₀²χ → 5.91(68), lit. 6.67(7); **QEX topoQ is CORRECT** (no factor 2) |
| slide 10 production | RESULTS.md | c₀(O(a⁴)) = 6.60(42); c₂ = −0.51(1.39) ≡ 0 vs cubic −2.70(51) |
| free fermions `hcfree` | `tfree` 19 | pressure series coefficients exact; slides 12–14 |

Conventions you must not re-litigate (all verified, see FORMULATION.md):
Wilson term is `arp²/2` (slide 11 has a factor-2 slip); Q carries the ½ site-volume
factor; force convention `d/ds S(e^{sP}U)|₀ = Σ redot(P,f)` with a **plus** sign;
QEX `dot` conjugates its FIRST argument; B-sublattice plane waves use the integer cell
coordinate (no half-site phase); the honeycomb β window for physics is 6.9–7.2.

## 5. Next steps, in priority order

1. **Close out D4** (§3) and fill REPRODUCTION.md §1's three ⏳ rows.
2. Re-run the full test suite (§2) once, to confirm the tree is green end-to-end.
3. **Nothing is committed.**  The whole of `src/experimental/honeycomb/` (and
   `build_mac/`, which should NOT be committed) is untracked.  Ask the user before
   committing; suggested shape: one commit adding `src/experimental/honeycomb/`
   (code + tests + doc + small .dat/.png plot files; the `doc/plots/*/data/` raw
   per-config files are a judgement call — they are small here).
4. Small cleanups (optional): consolidate `hcMeasFlow`'s local t₀ interpolator with
   `hcanalysis.findT0`; add an `accumulate` flag to `hcActionDeriv`; multi-rank (MPI)
   validation of `tgauge`/`taction` — everything is single-rank-tested only.
5. Cluster-scale production (the real run the user plans): REPRODUCTION.md §4 lists
   exactly what more statistics buys.  The executables take all parameters on the
   command line; no code changes needed.  Remember τ_int(Q) grows steeply with β
   (2 → 43 updates over β = 6.90 → 7.20).

## 6. File map

```
src/experimental/honeycomb/
  hcgeom.nim hclayout.nim hcgauge.nim hcaction.nim        # geometry → force
  hchmc.nim hcheatbath.nim hcio.nim hcPureGauge.nim       # Monte Carlo + I/O + generator
  hcflow.nim hctopo.nim hcMeasFlow.nim                    # flow, topology, measurement exe
  hcfree.nim hcFreeSpectrum.nim hcFreePressure.nim        # free fermions (slides 12–14)
  hcwilson.nim hcstout.nim hcclover.nim hcarnoldi.nim     # interacting fermions
  hcSpectrum.nim cubicSpectrum.nim                        # fermion measurement drivers (D4)
  hcanalysis.nim refCubicGen.nim refCubicMeas.nim         # analysis + cubic reference
  tests/  t{geom,gauge,action,hmc,flow,topo,wilson,stout,clover,arnoldi,free,analysis,spectrum}.nim
  doc/    HANDOFF.md REPRODUCTION.md FORMULATION.md SLIDES.md PLAN.md STATUS.md
          RESULTS.md RESULTS_CUBIC.md [RESULTS_FERMIONS.md]
  doc/plots/  {freespec,pressure}.*  cubic/  honeycomb/  fermions/
build_mac/            # configured build dir (bin/, nimcache/) — do not commit
$TMPDIR (/tmp/claude-502)/hcR/   # task R scratch: configs (disposable; results harvested)
```

## 7. Known open items / rough edges

* D4 closure (§3) — the only unfinished deliverable.
* `hc722.out` vs `hc722b.out`: the first honeycomb job stalled around config 45 and
  was superseded by `hc722b` (which resumed the chain and reports "n/76") — when
  harvesting, avoid double-counting configs present in both logs (harvest.py should
  handle it; verify).
* Multi-rank and Nc≠3 untested everywhere (single-rank SU(3) only).
* The unkillable leftover jobs will keep appending to their logs until they hit their
  own config targets; ignore output past your chosen harvest cutoff.
* An unrelated Go build (another user/session) shares `/tmp/claude-502` — ignore
  gocache/go-build noise when scanning for activity.
