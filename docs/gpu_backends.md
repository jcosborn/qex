# GPU backends

[`backend/accel`](../src/backend/accel.nim) provides `gpuFor`, memory
operations and reductions. Select a backend with `-d:Backend=...` in
`nimargs`; the default is CPU. The staggered operator, gauge operations,
HYP smearing, RNGs and HMC use this interface.

| Backend | Compiler | Kernel execution | IPC between GPUs on one host |
|---|---|---|---|
| CPU | Host compiler | Loop on the calling host thread | None (QMP messages) |
| OpenMP | oneAPI icx | `target teams distribute parallel for` | Level Zero |
| SYCL | oneAPI icpx | `parallel_for` | Level Zero |
| CUDA | clang or nvcc | Device lambdas in `qexFor` | CUDA IPC |
| HIP | amdclang | Device lambdas in `qexFor` | HIP IPC |

## Kernel calls and synchronization

Place `gpuFor` launches inside procs. Kernel bodies can use templates,
which expand in place, and procs that compile for the device. Captured
pointers refer to device storage. Keep their allocations alive until all
kernels and exchanges using them have completed.

CUDA and HIP require device annotations. The `gpuInline` and `gpuCall`
pragmas in [`base/basicOps`](../src/base/basicOps.nim) supply `QEX_HD`,
which expands to `__host__ __device__` in device compilation. A host/device
proc must also be valid device code: host globals and CPU SIMD types do
not belong in its body.

| Proc mark | OpenMP / SYCL kernels | CUDA / HIP kernels | Inlining |
|---|---|---|---|
| Template | Yes | Yes | Expansion in Nim |
| `{.inline.}` | Yes | No | Compiler decides |
| `{.alwaysInline.}` | Yes | No | Required |
| `{.gpuInline.}` | Yes | Yes | Required |
| `{.gpuCall.}` | Yes | Yes | Compiler decides |

OpenMP and SYCL need callable definitions in the generated C/C++ file;
Nim places inline procs in each file that uses them. Use `gpuInline` for
small arithmetic helpers and `gpuCall` for larger helpers such as RNG
draws. `alwaysInline` is suitable for host SIMD wrappers. The older
kernel-block API uses `inlineProcs:` from `base/metaUtils` to copy called
bodies into the kernel; `gpuFor` does not perform that transformation.

`gpuFor` completes its kernel before returning. `gpuForAsync` submits work;
`gpuWaitAsync` waits for queued kernels. OpenMP synchronous launches are
not ordered after outstanding `nowait` kernels: wait before mixing the
two forms or consuming their results on the host.

With the CPU backend, both forms execute on the calling host thread and
complete on return. `gpuWaitAsync` is a no-op, device allocations are host
allocations, and every halo message uses QMP. Increasing the host thread
count does not parallelize a `gpuFor` loop called by one thread.

## Halo exchange

[`GpuHaloEx`](../src/comms/halogpu.nim) owns gather maps, buffers and
communication handles. `pack` queues stores into send slots. `start`
waits for the producing kernels before starting messages; `wait` completes
the exchange before its receive data is consumed. With no messages,
`wait(sync=true)` still waits for queued kernels.

For peers on the same host, IPC maps the receiving GPU's memory into the
sender. Kernels store boundary values there directly; QMP messages signal
completion. Other messages carry device buffers through QMP and require
GPU-aware MPI. `haloIpc=false` (`-ipc:0` in the programs) uses that message
path for peers on the same host as well.

Ranks sharing an exchange must agree on the sequence of sends and receives.
Finish pending kernels and messages before reusing or freeing their buffers.

## Staggered operator and storage

[`StagGpu`](../src/physics/stagGpu.nim) implements the four-dimensional
staggered operator with three colors. It owns links, neighbor maps, halo
exchanges and solver scratch. Construct it from phased gauge fields with
compatible layouts. `upload` and `download` require vector fields with
matching layout, SIMD width and scalar precision.

A site `k` and real component `c` have vector offset
`(k div V)*(6*V) + c*V + k mod V`. Even sites precede odd sites. The
`solveEE` field interface uploads and downloads only the even part.
The pointer interface to `solveM` requires storage for all `6*n` reals of
the solution, even when the source is restricted to even sites.

`solveM` source and destination device buffers are disjoint from each other
and from the operator's work vectors. An instance's scratch and exchanges serve one
solve at a time. Copies of a `StagGpu` share its allocations; release the
owner with `free` once pending work has completed. Construct with
`batch=true` for calls containing several systems.

### Link representation

`newStagGpu(g, T, reals=12, fwd=-1, ...)` selects storage from the requested
format and the gauge. Its `nl` field gives the actual number of stored
reals per link.

| Reals | Link condition | Stored data and reconstruction |
|---|---|---|
| 18 | General 3x3 complex matrix | All entries |
| 12 | `U = s W`, `W` in SU(3), `s = ±1` | First two rows; row 2 is `s conj(row 0 × row 1)` |
| 14 | Unitary 3x3 matrix | First two rows and `det U`; row 2 is `det U conj(row 0 × row 1)` |

A request for 12 reals falls back to 14 when the sign representation is
insufficient, then to 18 when determinant reconstruction is insufficient.
A request for 14 falls back to 18. The construction checks the third row
at relative norm tolerance 1e-12 across all ranks.

Both forward and backward caches store links in forward orientation.
A backward hop takes the adjoint of the reconstructed forward matrix;
independently compressing its adjoint would define a different operator.
With `fwd=1`, only forward links are stored, and backward hops fetch the
neighbor's link. `fwd=0` stores both caches. Automatic selection (`-1`)
uses forward storage when both caches would exceed 128 MiB.

Compressed FP32 links have two optional reconstruction modes:

- `reconFma=true` (`-reconFma:1`) uses explicit FP32 fused multiply-add
  operations for the cross product and complex determinant factor.
- `recon64=true` (`-recon64:1`) reads FP32 entries, reconstructs in FP64,
  and rounds the result to FP32 before the matrix-vector product. It takes
  precedence over `reconFma`.

Both default to false. Reconstruction changes arithmetic, not storage
precision. Scalar and batched hops use the same matrix for local and halo
links. These options do not affect FP64 operators or 18-real links.

`setLinks` updates caches from a `GpuGauge` and `stagSigns`; it preserves
the selected format and geometry. The gauge must satisfy that format's
link condition, including SU(3) for 12-real storage before the signs.
Mixed solves require both precision objects to represent the same current
gauge; update both after changing the links.

### Equation and stopping contract

With `M = m + D/2` and `A = 4m² - D_eo D_oe`, FP64 and mixed GPU solves
start from zero and require:

- `solveEE`: `||b_e - A x_e||² <= r2req ||b_e||²`.
- `solveM`: `||b - M x||² <= r2req ||b||²`. Unless `full=true`, the odd
  source is taken as zero and is not read. The output includes both parities.

`solveM` uses a nonzero mass and the even system

```text
q   = 4m b_e - 2 D_eo b_o
A x_e = q
x_o = b_o/m - D_oe x_e/(2m)
```

In exact arithmetic, the full residual is `(q-Ax_e)/(4m)` on even sites
and zero on odd sites. The full solve uses its absolute target
`stop = r2req ||original source||²`. For each correction source `b`,
the even CG target is
`16m² (stop-||b_o||²)` when `||b_o||² <= stop/2`, otherwise
`0.99 * 16m² * stop`. The latter leaves a margin for reconstruction.

Convergence is decided by the residual recomputed with the precise
operator. A trial correction replaces the accepted solution only if it
reduces that residual or meets the target. A failed mixed correction
retries in FP64 on the device; a non-improving FP64 correction ends the
solve. All attempts share `maxits`. The returned statistics report the
accepted solution's actual residual, including when the target cannot be
reached. A zero source returns zero with no CG iterations and relative
residual zero.

Full solves check both parities after reconstruction and refine the full
equation. This matters at small mass: the even residual has arithmetic
error of order `epsilon * ||A|| * ||x_e||`, amplified by `1/(4m)` in the
reconstructed full solution. A small correction `M y = b-Mx` has a smaller
rounding error than further refinement of the large `x_e` alone.
Reconstruction may also make progress without a CG iteration.

`SolverParams` accumulates calls and iteration counts. Its `r2` statistic
uses relative **squared** residuals; `reliable` counts inner reliable
updates separately from outer refinement checks. The precise operator is
the one defined by the stored links, including any link reconstruction.
The FP32-only `solveEE` overload recomputes residuals in FP32 and stops
when a restart no longer improves them.

### CG modes and communication

With an FP32 inner operator, `gpuCg` selects:

| `gpuCg` | Inner iteration | Partial solution | Batch scheduling |
|---|---|---|---|
| 0 (default) | Chronopoulos–Gear, one global reduction per step | FP32 | Refill completed slots |
| 1 | Conventional CG with reliable residual updates | `gpuAcc64` selects FP32/FP64 | Complete groups of up to `nBatch` |
| 2 | Chronopoulos–Gear | `gpuAcc64` selects FP32/FP64 | Refill completed slots |

The residual, direction and stencil vectors are FP32 in mixed solves;
dot products use FP64 products and reductions. Mode 2 with `gpuAcc64=1`
uses FP64 alpha for correction accumulation while keeping the recurrence
in FP32. A nonpositive curvature ends an inner correction so the outer
solve can check its result and recover.

The one-reduction recurrence uses `w = A r` and

```text
gamma = r.r; delta = w.r
beta  = gamma/gamma_previous
alpha = gamma/(delta - beta*gamma/alpha_previous)
p = r + beta*p; v = w + beta*v
x += alpha*p; r -= alpha*v
```

The first step has `beta=0`, `alpha=gamma/delta`. Each rank updates halo
copies of `r` and `v` by the same formulas. The first hop of the next
`A r` therefore reuses those halo copies; exchange of `w` overlaps the
global reduction. `hopSplit` (`-split`) controls overlap of the second
hop's local work with the first hop's halo exchange: `1` always splits,
`0` does not, and `-1` splits when there are neighbors on another node.

`nBatch` is a compile-time integer, default 3. Batched hops share each
loaded link across systems. Each system has its own residual history and
iteration budget. A failed mixed system uses the separate scalar FP64
workspace for recovery. `StagGpu.fixed` selects fixed-order reductions;
atomic reductions permit variation in the order of additions. The
programs expose this as `-fixed`.

### Reliable updates and controls

Reliable CG periodically accumulates the partial solution into the precise
total and replaces its recursive residual with `b-Ax`. Retaining the
direction projects it orthogonal to this new residual before continuing.
The optional modified beta numerator is `Re(r_new†(r_new-r_old))`, with
`||r_new||²` used if that numerator is negative; both divide by `||r_old||²`.

A recursive residual predicting convergence, the arithmetic limit below,
or a reduction by `Delta` relative to the last precise norm requests a check. A reduction
relative to the largest recursive norm since that check, or reaching
`Period`, also requests an update. Growth counts toward `MaxInc` and
`MaxTotal` only when the recurrence predicted progress: CG residual norms
can increase at ordinary periodic checks.

For an inner source `q`, GPU modes 0 and 2 use the squared stopping target
`max(requested, max(r2in, epsilon(FP32)) * ||q||²)`. Thus `r2in=0` selects
the FP32 epsilon bound, about 1.19e-7. Mode 1 uses
`max(requested, r2in * ||q||²)` and can resolve smaller targets through
reliable updates.

The reliable inner correction also estimates its arithmetic accuracy as

```text
eta = Floor * epsilon(outer precision) * (||q|| + Anorm * ||y||)
Anorm = ||A q|| / ||q||
```

`Anorm` is estimated from the first FP32 operator application, not a
rigorous operator bound. The correction can stop at this estimate once
it improves on its source. GPU norms share existing paired reductions.
The outer solve still enforces the original requested tolerance. A zero
`Floor` disables the estimate. Disabling both `Delta` and `Period` leaves
checks at candidate convergence, the arithmetic estimate, and the GPU
iteration limit.

CPU controls live in `SolverParams.cg`; GPU controls live in the FP32
`StagGpu.ctrl`. Their command-line arguments and defaults are:

| CPU argument | GPU argument | CPU / GPU default | Meaning |
|---|---|---|---|
| `-cpuCg` | `-gpuCg` | -1 / 0 | CPU: automatic, ordinary CG (0), or reliable CG (1); GPU modes above |
| `-cpuDelta` | `-gpuDelta` | 0.01 / 0.1 | Norm reduction factor; 0 disables this trigger |
| `-cpuPeriod` | `-gpuPeriod` | 1000 / 1000 | Maximum iterations between precise updates; 0 disables this trigger |
| `-cpuFloor` | `-gpuFloor` | 1 / 1 | Arithmetic error estimate factor |
| `-cpuAcc64` | `-gpuAcc64` | 0 / 1 | FP64 partial accumulation; reliable CG also widens fused direction arithmetic |
| `-cpuBeta` | `-gpuBeta` | 0 / 1 | Modified beta numerator when enabled |
| `-cpuKeep` | `-gpuKeep` | 1 / 1 | Retain and reorthogonalize the direction |
| `-cpuMaxInc` | `-gpuMaxInc` | 1 / 1 | Allowed consecutive unexpected residual increases |
| `-cpuMaxTotal` | `-gpuMaxTotal` | 10 / 10 | Allowed total unexpected residual increases |
| `-cpuR2in` | `-r2in` | 0 / 1e-6 | Minimum relative squared inner target |

The reliable controls affect GPU mode 1; `gpuAcc64` also affects mode 2.
CPU automatic selection uses reliable CG for `solveReconL`, and ordinary
CG corrections for `solveReconR` and standalone `solveEE`/`solveOO`.
The reconstruction choice follows the current source's parity norms.
CPU mixed precision is enabled with `sloppySolve=1`. GPU mixed precision
uses a separate FP32 operator; the programs expose `-mixed:1`.

## Build configuration

Use [configure](../INSTALL.md#configuration-examples) in a separate build
directory. Each build needs its own `nimcache`. Set `QEX_SOURCE`,
`QMP_PREFIX` and `QIO_PREFIX` to the source directory and library prefixes
for your installation. The examples use x86 hosts and SIMD width 8; layout geometry
must support that width. `-march=native` requires building on the same
CPU architecture as the compute nodes; otherwise name that architecture
explicitly. Static QMP/QIO libraries built without PIE require `-no-pie`
in the link flags.

### OpenMP, Intel PVC

Use MPI wrappers for icx/icpx and the oneAPI runtime:

```sh
"$QEX_SOURCE/configure" \
  qmpdir:"$QMP_PREFIX" qiodir:"$QIO_PREFIX" \
  cctype:clang cc:mpicc cpp:mpicxx \
  cflagsspeed:"-O3 -march=native" cppflagsspeed:"-O3 -march=native" \
  ldflags:"-g -O3 -Xs '-device pvc' -ftarget-register-alloc-mode=pvc:large -ldl" \
  ldppflags:"-g -O3 -Xs '-device pvc' -ftarget-register-alloc-mode=pvc:large -ldl" \
  env:"OMPFLAG=-fiopenmp -fopenmp-targets=spir64_gen" \
  simd:auto vlen:8 nimargs:'@["-d:Backend=OpenMP"]'
```

Device compilation occurs at link time and needs the optimization flag
there as well. `-ftarget-register-alloc-mode=pvc:auto` permits the compiler
to select 128 or 256 registers per kernel.

### SYCL, Intel PVC

```sh
"$QEX_SOURCE/configure" \
  qmpdir:"$QMP_PREFIX" qiodir:"$QIO_PREFIX" \
  cctype:clang ccdef:cpp cc:mpicc cpp:mpicxx \
  cflagsspeed:"-O3 -march=native" cppflagsspeed:"-O3 -march=native" \
  cppflagsalways:"-g -fsycl -fsycl-targets=spir64_gen" \
  ldppflags:"-g -O3 -fsycl -fsycl-targets=spir64_gen -Xsycl-target-backend '-device pvc' -ftarget-register-alloc-mode=pvc:large -ldl" \
  env:"OMPFLAG=-fiopenmp" simd:auto vlen:8 nimargs:'@["-d:Backend=SYCL"]'
```

The C++ MPI wrapper must invoke icpx.

### CUDA with clang

Set `CUDA_HOME` to a toolkit supported by clang and `GPU_ARCH` to the
device architecture, such as `sm_90` for H100 or `sm_100` for B200. These
wrapper settings use OpenMPI; MPICH uses `MPICH_CC` and `MPICH_CXX`.

```sh
"$QEX_SOURCE/configure" \
  qmpdir:"$QMP_PREFIX" qiodir:"$QIO_PREFIX" \
  cctype:clang ccdef:cpp cc:mpicc cpp:mpicxx \
  env:"OMPI_CC=clang" env:"OMPI_CXX=clang++" \
  cflagsspeed:"-O3 -march=native" cppflagsspeed:"-O3 -march=native" \
  cppflagsalways:"-g -x cuda --cuda-gpu-arch=$GPU_ARCH --cuda-path=$CUDA_HOME -Xarch_device -mllvm=-disable-machine-sink" \
  ldppflags:"-g -L$CUDA_HOME/lib64 -Wl,-rpath,$CUDA_HOME/lib64 -lcudart -ldl" \
  simd:"SSE,AVX,AVX512" vlen:8 nimargs:'@["-d:Backend=CUDA"]'
```

`backend/cuda/nimbase.h` supplies `QEX_HD` and `qexFor`. The device flag
`-disable-machine-sink` works with the `[[clang::always_inline]]` call in
`qexFor`: together they keep batched-hop arithmetic near its loads and
limit register spills. To link where the NVIDIA driver is unavailable,
add `-L$CUDA_HOME/lib64/stubs`; the driver is required at runtime.

### CUDA with nvcc

[`build/nvcc.sh`](../build/nvcc.sh) translates Nim's clang command line to
nvcc for files containing kernels and uses `QEX_HOSTCXX` for other files.
Set `CUDA_HOME`, `GPU_ARCH` and the MPI installation prefix `MPI_HOME`.
Use a host compiler supported by the CUDA toolkit.

```sh
"$QEX_SOURCE/configure" \
  qmpdir:"$QMP_PREFIX" qiodir:"$QIO_PREFIX" \
  cctype:clang ccdef:cpp cc:mpicc cpp:"$QEX_SOURCE/build/nvcc.sh" \
  env:"OMPI_CC=gcc" env:"OMPI_CXX=g++" env:"QEX_HOSTCXX=g++" \
  env:"CUDA_HOME=$CUDA_HOME" env:"CPATH=$MPI_HOME/include" \
  cflagsspeed:"-O3 -march=native" cppflagsspeed:"-O3 -march=native" \
  cppflagsalways:"-g -x cuda --cuda-gpu-arch=$GPU_ARCH" \
  ldpp:mpicxx \
  ldppflags:"-g -L$CUDA_HOME/lib64 -Wl,-rpath,$CUDA_HOME/lib64 -lcudart -ldl" \
  simd:"SSE,AVX,AVX512" vlen:8 nimargs:'@["-d:Backend=CUDA"]'
```

`CPATH` supplies headers to the host compiler and nvcc's preprocessor;
the MPI wrapper performs the link.

### HIP, AMD MI300X / MI300A

Use MPI wrappers for amdclang/amdclang++ and a ROCm environment:

```sh
"$QEX_SOURCE/configure" \
  qmpdir:"$QMP_PREFIX" qiodir:"$QIO_PREFIX" \
  cctype:clang ccdef:cpp cc:mpicc cpp:mpicxx \
  env:"OMPI_CC=amdclang" env:"OMPI_CXX=amdclang++" \
  cflagsspeed:"-O3 -march=native" cppflagsspeed:"-O3 -march=native" \
  cppflagsalways:"-g -x hip --offload-arch=gfx942" \
  ldppflags:"-g -ldl --hip-link --offload-arch=gfx942" \
  env:"STATIC_UNROLL=1" simd:"SSE,AVX,AVX512" vlen:8 nimargs:'@["-d:Backend=HIP"]'
```

`gfx942` covers both MI300X and MI300A. `backend/hip/nimbase.h` supplies
`QEX_HD`, `qexFor` and the HIP declarations. Generate C++ for all modules
with `ccdef:cpp` so host/device declarations have compatible linkage.

## Device placement and programs

CUDA, HIP and SYCL choose global rank `myRank mod` the number of visible devices;
OpenMP uses the default device. Set visibility per rank before `qexInit`
when using several nodes or a specific GPU/tile mapping. Bind each rank's
host threads near its GPU. For Intel tiles, visibility can be restricted
with `ZE_AFFINITY_MASK`.

For MPICH on PVC, `MPIR_CVAR_CH4_IPC_GPU_P2P_THRESHOLD=0` permits GPU IPC
for small messages and `MPIR_CVAR_CH4_OFI_ENABLE_HMEM=1` enables the device
buffer path between nodes. Cray MPICH uses `MPICH_GPU_SUPPORT_ENABLED=1`.
Set `OMP_TARGET_OFFLOAD=MANDATORY` for OpenMP device execution. Large local
volumes may need a larger host OpenMP stack, such as `OMP_STACKSIZE=256M`.

| Program | Purpose |
|---|---|
| `backend/examples/bestream` | Device copy, scale, add and triad bandwidth |
| `backend/examples/bestagcg` | Even-system CG throughput; `-nb` selects source count |
| `backend/examples/bestagm` | Complete `M` solves; `-masses`, `-nb`, `-full`, `-mixed`, `-r2req`, `-r2in` |
| `backend/examples/bestagres` | Recursive and recomputed residuals, iteration counts and source/statistics checks |
| `backend/examples/berng`, `bernt` | RNG output and kernel throughput |
| `examples/staghmcgpu_sh` | Staggered HMC with GPU smearing and force evaluation |
| `prod/lsd/eightFlavorSMGgpu` | GPU eight-flavor production application |

Build `bestagm` with `-d:cpuOnly` for the production CPU variant. Compiling with
`-d:stagWorkCount` enables per-rank stencil-site, launch, CG-start and
reduction counters. Reconstruction with FP64 intermediates still counts
as FP32 stencil work when the fields and matrix-vector products are FP32.

## Representative PVC performance

Complete `solveM` calls for three full sources on random SU(3) links,
masses 0.05/0.06/0.07, relative squared tolerance 1e-16, 14-real links,
fixed reductions and SIMD width 8. Both reconstruction options are off.
Mode 0 uses `r2in=1e-6`; mode 1 uses `r2in=0` and the reliable defaults above.
Times are medians in milliseconds per three-source call.

Machine configuration: Intel PVC on Sunspot, six GPUs / twelve tiles per
node, Xeon CPU Max 9470C hosts, one rank per tile and eight bound CPU cores
per rank. Compilers are oneAPI 2026.1.0 icx/icpx with MPICH 5.0.0,
`-O3 -march=native` and `-ftarget-register-alloc-mode=pvc:large`.

| Global lattice | Nodes / ranks | Rank geometry | SYCL mode 0 | SYCL mode 1 | OpenMP mode 0 | OpenMP mode 1 |
|---|---|---|---:|---:|---:|---:|
| 16^4 | 1 / 8 | 2x2x2x1 | 82.68 | 101.58 | 86.68 | 104.99 |
| 32^4 | 1 / 8 | 2x2x2x1 | 281.05 | 269.92 | 288.74 | 280.64 |
| 32^4 | 1 / 8 | 1x2x2x2 | 249.68 | 269.15 | 257.70 | 280.90 |
| 48^4 | 1 / 12 | 3x2x2x1 | 632.79 | 586.19 | 638.00 | 600.85 |
| 48^4 | 2 / 24 | 2x2x2x3 | 497.92 | 501.64 | 505.08 | 506.13 |

Local volume and rank geometry affect the cost of reliable updates and
halo exchanges. The `gpuCg` modes make that tradeoff explicit.
