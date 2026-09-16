# Graph validation and benchmarks

See [Graph design](../src/experimental/graph/DESIGN.md) for runtime and numerical
contracts.

## Optional test generation

`make tests experimental` builds optional fixtures and generates
`testscript-experimental.sh`; `make tests experimental/graph` selects this group.
The generated script uses `RUNJOB` for MPI-capable tests and `RUN1` for tests
intended for one rank. From an MPI-enabled configured build directory, generate
and run the optional script explicitly:

```sh
make tests experimental/graph
RUNJOB='mpiexec -n 2' RUN1='mpiexec -n 1' ./testscript-experimental.sh
```

The default GitHub pipeline and its `runtests.sh` MPI pass run the normal test
script; they do not run optional graph tests automatically. QEX graph test
drivers using `tests/helpers` disable unittest's automatic command-line name
filtering. Their QEX options use `-OPTION:VALUE`, and all unittest cases execute;
no trailing `'*'` selector is needed. Check nonzero executed case counts as well
as exit status. Generated commands do not set `-expectRanks`, so explicit MPI
gates should supply the intended rank count and an appropriate geometry.

## Graph storage benchmark

From a configured build directory, build with `make graph/benchGraphStorage` and
run `./bin/benchGraphStorage`. Each invocation is one fresh process for one case
and mode. Use the same executable, lattice, rank count/geometry, threads and
options for both modes. Record the compiler flags, precision, colour count and
SIMD configuration with the logs. Deterministic HMC inputs and gauge fingerprint
weights depend on the local layout, so changing the decomposition changes the
comparison data.

| Option | Default | Applies to |
| --- | --- | --- |
| `-case:chain`, `branch`, `gradient`, `hmc` | `chain` | All runs |
| `-mode:direct` or `planned` | `planned` | All cases |
| `-lat:4,4,4,4` | Global lattice derived from local `8^4` for chain, `4^4` otherwise | All cases; an explicit value is the global lattice |
| `-threads:1` | `1` | All cases; overrides `OMP_NUM_THREADS` |
| `-reps:N` | `20`, or `3` for HMC | Calls in each phase after the one-call warm phase; positive |
| `-delta:0.0001` | `0.0001` | Deterministic increments in changed phases |
| `-steps:24` | `24` | Chain only; positive |
| `-inplace:true` or `false` | `false` | Chain only; permits the scale node to reuse its gauge input |
| `-alpha:1.01` | `1.01` | Chain, branch and gradient; HMC starts with beta `6` |

Nonapplicable `steps`, `inplace` and `alpha` options are accepted and ignored;
matched cases print `ignored_options`. Branch, gradient and HMC use the stock
operators' declared alias/reuse contracts. The benchmark has no unittest name
filters.

Both modes evaluate the same requested results. The `roots` header is the
ordered root contract. Branch, gradient and HMC wrap every slot in one
`multiValues` evaluation root, including the carrier's evaluation in both modes.

| Case | Requested results |
| --- | --- |
| `chain` | `chain.final.gauge`, after `steps` scalar-times-gauge nodes |
| `branch` | `application`, `independent`, `application.norm2` |
| `gradient` | `value`, `loss`, `dLoss/dAlpha`, `dLoss/dBeta`, `dLoss/dX`, `dLoss/dY` |
| `hmc` | `initial.H`, `initial.S`, `initial.T`, `force[0].norm2`, `force[1].norm2`, `final.gauge`, `final.momentum`, `final.H`, `final.S`, `final.T`, `dH`, `loss`, `dLoss/ddt`, `dLoss/dlambda` |

Branch evaluates `z=alpha*x*x+y` or `z=beta*x+y*y` through a conditional lambda;
the independent root is `beta*y`. Gradient uses `z=alpha*x+beta*y` and
`loss=norm2(z)`. Both start with identity-scaled fields `x=1`, `y=0.75`, beta
`0.625` and the requested alpha.

HMC measures a Wilson proposal graph with 2MN, one integration step, initial
`dt=0.05`, force norm roots and training expressions. Its construction timing
includes the trajectory constructor's copies of numerical inputs in both modes.
The complete `runHmc` sampler has additional RMS/min/max diagnostics, snapshots,
acceptance, callbacks, reverse checks and measurements; those operations are
outside this benchmark's workload.

Every case runs `warm`, `unchanged`, `all-inputs-changed` and `partly-changed` in
that order. Chain inserts `updated` before `all-inputs-changed` to change alpha
alone. Its all-input phase changes alpha and the gauge; its partial phase changes
the gauge. Branch and gradient change x, y, alpha and beta together; branch also
alternates its selector. Their partial phase changes x. HMC changes gauge,
momentum, beta, dt and lambda together, then momentum alone. Direct evaluation
can retain unaffected subgraphs in partial phases, so compare the reported work
and time without requiring equal forward counts across modes.

`tic`/`toc` bracket graph construction and each evaluation. Input changes,
fingerprint reductions, analytic/identity checks and explicit phase collections
are outside evaluation timing. A raw allocator threshold can still trigger a
collection during evaluation; the log records `raw_gc_threshold`.
`-d:graphPlanMemory` adds collections and printing during plan construction or
rebuild, and `-d:nimAllocStats` adds allocator instrumentation. Either flag marks
the log `timing_run_kind=instrumented`; keep these logs separate from ordinary
timing comparisons and use the same flags in a direct/planned pair.

Memory baselines are captured after input setup and a full collection, before
graph construction. The `*_above_inputs` fields subtract that whole baseline,
which also contains layout and runtime state. `input_raw_bytes` describes only
the numerical input gauge payloads.

| Counter | Meaning |
| --- | --- |
| `construction_raw_bytes`, `eval_raw_allocated` | Cumulative raw allocation deltas during construction and timed evaluations, respectively |
| `phase_raw_allocated` | Raw allocation delta over the phase loop, including work outside evaluation timing |
| `allocated_above_inputs`, `cumulative_raw_allocated` | Cumulative raw allocation since the baseline and since process start; frees do not reduce these counters |
| `current_above_inputs`, `current_after_gc_above_inputs` | Current raw occupancy minus the baseline, before and after the explicit collection |
| `process_raw_peak_above_baseline` | `max(0, process lifetime raw high water - baseline raw occupancy)`; earlier phases and pre-baseline allocations can determine it |
| `managed_occupied*`, `managed_heap` | Nim managed occupancy and reserved heap, separate from raw gauge allocations |
| `process_rss_peak_bytes` | OS process lifetime RSS high water in bytes, or `-1` on an unsupported OS |
| Arena `bytes`, `peak_live_bytes`, `workspace_bytes` | Planner buffer/workspace accounting; excludes other graph metadata and process memory |

The raw/RSS maxima are never reset per phase. In diagnostic plan reports,
`nodePoolBytes` counts one pool-index integer per plan node, not all plan metadata.
QEX prints these counters from rank zero; they describe that process. Gauge
fingerprints use reductions across all ranks.
This benchmark keeps its graph, plan and results alive until process exit and
does not measure teardown or establish reclamation.

### Paired HMC validation

Chain checks an analytic norm; branch and gradient check analytic fingerprints.
HMC checks finite root fingerprints, preserved inputs and the internal identities
`H=S+T`, `dH=Hfinal-Hinitial` and `loss=-min(1,exp(-dH))*dt^2`. These identities
use the same evaluated bundle, so a successful HMC process alone does not
establish agreement between modes. Its `fingerprint_comparison=external-required`
marks the emitted data for the following separate comparison.

Run each command to completion and require exit status zero. These commands use
the same executable in separate processes; retain both complete logs:

```sh
make ARGS="--assertions:on" graph/benchGraphStorage
mkdir -p graph-storage-logs
./bin/benchGraphStorage -case:hmc -mode:direct -lat:4,4,4,4 -threads:1 -reps:3 -delta:0.0001 >graph-storage-logs/hmc-direct.log 2>&1
./bin/benchGraphStorage -case:hmc -mode:planned -lat:4,4,4,4 -threads:1 -reps:3 -delta:0.0001 >graph-storage-logs/hmc-planned.log 2>&1
```

The paired-log check must:

1. Match all run settings except `mode`, including instrumentation, lattice,
   ranks, threads, reps, delta, ordered root names, root count and fingerprint
   layout. For these commands require `ranks=1`, `root_slots=14` and four phases
   with call counts `1,3,3,3`.
2. Require exactly one fingerprint record for every `(phase, root index)` in
   each log: 56 records for this fixed HMC graph. Reject missing or duplicate
   records, changed root names, and wrong component counts. Scalars have one
   component; gauges have three: norm squared, weighted real sum and weighted
   imaginary sum.
3. Compare every component of both `fingerprint_checksum` and `fingerprint_last`.
   Reject nonfinite values. For the default double precision build require
   `abs(planned-direct) <= 1e-10 * max(1, abs(direct))`, and report the largest
   normalized difference. Checksums weight each call by its one-based index;
   `fingerprint_last` retains the final call in that phase.
4. Require the internal assertions to pass and every reported HMC
   `max_relative_error` to be finite and below `1e-10`. Compare performance and
   memory only after the fingerprint comparison passes. Do not require equal
   allocation or forward counts as an equivalence condition.

This check establishes agreement of the emitted summaries for these requested
roots and phases. The three gauge summaries are not an elementwise proof; the
graph/HMC test suites supply independent derivative and trajectory coverage.
Repeat fresh-process pairs for other cases, lattices or instrumentation builds
and retain their settings with each result.
