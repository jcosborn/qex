# GPU backends

The GPU code (backend/accel gpuFor kernels, physics/stagGpu, gauge/hypGpu,
gauge/gaugeGpu, comms/halogpu, hmc/hmcActionGpu, rng/rngGpu) builds with four
GPU backends, selected by `-d:Backend=...` in `nimargs`, and with the CPU
backend, the default:

| Backend | Compilers | Kernels | Halo stores into peers on the node |
|---|---|---|---|
| CPU | the host compiler | loops on the calling host thread | none, QMP messages |
| OpenMP | icx (oneAPI), OpenMP offload | `target teams distribute parallel for` | Level Zero IPC (comms/zeipc) |
| SYCL | icpx (oneAPI) | `parallel_for` lambdas | Level Zero IPC (comms/zeipc) |
| CUDA | clang `-x cuda`; nvcc through build/nvcc.sh | device lambdas of the qexFor<<<>>> template | CUDA IPC (comms/cudaipc) |
| HIP | amdclang `-x hip` | device lambdas of the qexFor<<<>>> template | HIP IPC (comms/hipipc) |

The tested builds (September 2026) pass both CI suites,
tests/extra/tstaghmc_sh and tests/extra/teightFlavorSMG, on 1, 2 and 4
GPUs (PVC tiles):

| GPU | System | Host CPU | Compiler | MPI |
|---|---|---|---|---|
| Intel PVC | Sunspot | Xeon CPU Max 9470C | oneAPI 2026.1.0: icx (OpenMP), icpx (SYCL) | MPICH 5.0.0 |
| NVIDIA H100 SXM | JLSE hopper00 | Xeon Platinum 8468 | clang 22.1.8 with CUDA 12.9.1; nvcc 13.3 with g++ 14.3 | OpenMPI 4.1.1 |
| NVIDIA B200 | JLSE blackwell00 | Xeon 6960P | clang 22.1.8 with CUDA 12.9.1 | OpenMPI 4.1.1 |
| AMD MI300X | JLSE amdgpu00 | EPYC 9654 | ROCm 10.0.0 amdclang | OpenMPI 4.1.1 |
| AMD MI300A | Tuolumne | MI300A (Zen 4 cores) | ROCm 10.0 amdclang | Cray MPICH 9.1.0 |

On September 30 devel, with the solver stopping rules below, passed on
all of these builds: both CI suites, tests/base/taccel and tgaugegpu,
backend/examples/bestagres (also with `-split:0` and `-split:1`, links of
12 reals and an 8^4 lattice, and with `-ipc:0` where MPI reads device
memory, on Sunspot and Tuolumne), berng and staghmcgpu `-check:1`, on 1, 2
and 4 GPUs and on 12 PVC tiles.
The stopping rules first ran with the CPU backend, on a FreeBSD 15.1 host
(Xeon E5-2687W v2, clang 19.1.7, MPICH 5.0.1, `vlen:16`).

## Procs that kernels call

A gpuFor body may use templates freely: they expand in place.  The procs it
calls must also compile for the device.

- OpenMP and SYCL compile for the device the functions a kernel calls that
  are defined in the same C or C++ file.  Nim writes each inline proc into
  every file that uses it, so kernels call inline procs without marks.
- CUDA and HIP compile for the device only the functions marked
  `__device__`; clang, amdclang and nvcc reject a kernel's call to any other
  function.  The procs that kernels call carry gpuInline or gpuCall
  (base/basicOps), which put `QEX_HD`, `__host__ __device__` from
  backend/cuda/nimbase.h and backend/hip/nimbase.h, in front of their
  declarations.  The other procs stay host only: nvcc compiles every host
  device function for the device as well and rejects host globals and x86
  vector types there.  build/nvcc.sh gives the files without kernels to the
  host compiler, where `QEX_HD` is empty.

The declarations Nim writes for a proc `f` returning `T`:

| Mark | CUDA and HIP | Other backends | Inlined | Kernels may call it |
|---|---|---|---|---|
| template | expanded in place | same | always | with every backend |
| `{.inline.}` (Nim) | `static N_INLINE(T, f)`, that is `static inline T f` | same | as the compiler decides | with OpenMP and SYCL |
| `{.alwaysInline.}` (basicOps) | `inline __attribute__((always_inline)) T f` | same | always | with OpenMP and SYCL |
| `{.gpuInline.}` (basicOps) | `QEX_HD inline __attribute__((always_inline)) T f` | as alwaysInline | always | with every backend |
| `{.gpuCall.}` (basicOps) | `QEX_HD static N_INLINE(T, f)` | as `{.inline.}` | as the compiler decides | with every backend |

Which to use:

- gpuInline: small procs that kernels call, where inlining always pays: the
  float `+=`, `-=` and `*=` of basicOps and the 3x3 matrix procs of
  gauge/hypGpu.
- gpuCall: larger procs that kernels call, left to the compiler: the draws
  of the RNGs in rng/.  With gpuInline, clang's Philox4x64 and Threefry4x64
  kernels ran 13-33% slower on H100 (berng 32^4, kernel times from nsys);
  with gpuCall, clang and amdclang compile the same kernels as for Nim's
  inline procs, and nvcc inlines the draws either way.
- alwaysInline: host code that must be inlined, such as the SIMD wrappers.
- `{.inline.}`, `{.noinline.}` and procs without a mark: host code.  OpenMP
  and SYCL kernels can call the inline ones, CUDA and HIP kernels none.
- The macro `inlineProcs:` (base/metaUtils) copies into its block the
  bodies of the procs called there, in Nim.  The kernels of the older API
  (the cuda, hip, sycl and OpenMP kernel blocks of backend/, and
  backend/vectorized) use it; gpuFor bodies don't.

## CPU backend

Without `-d:Backend`, backend/cpu runs each gpuFor and gpuForAsync kernel
as a loop over its indices on the calling host thread, complete when it
returns; gpuWaitAsync does nothing, gpuMalloc and gpuMallocHost allocate
host memory, and GpuHaloEx sends every message through QMP.  The GPU
programs build and run there as the reference of the device code, on one
core per rank.  tests/base/taccel and tests/base/tgaugegpu run in the CI
suite this way.

## Solver stopping rules

physics/stagGpu solves with M = m + D/2 and A = 4m^2 - D_eo D_oe, as
stag.solve and solveEE:

- solveEE until |b_e - A x_e|^2 <= r2req |b_e|^2;
- solveM until |b - M x|^2 <= r2req |b|^2, with b_o = 0 unless `full`
  (then b_o is not read), through A x_e = q = 4m b_e - 2 D_eo b_o and
  x_o = b_o/m - D_oe x_e/(2m).  (b - M x)_e = (q - A x_e)/(4m) and
  (b - M x)_o = 0 in exact arithmetic, so the CG stops at
  |q - A x_e|^2 <= 16m^2 (r2req |b|^2 - |b_o|^2) where
  |b_o|^2 <= r2req |b|^2/2, as solveReconR and solveEE of the CPU, and
  else at 0.99 16m^2 r2req |b|^2, as solveReconL, whose margin covers the
  rounding of x_o.  HMC trajectories then follow those of hmcAction to
  rounding.

The CG iterations stop on the recursive residual.  solveEE then computes
the true residual b - A x, and while it misses the request and drops, with
iterations left, restarts the CG from it with the same x.  solveM computes
b - M x instead, and while it misses the request and drops, with
iterations left, adds to x the solution y of M y = b - M x, as stag.solve
does.  solveEE and solveM record the true relative squared residual in
SolverParams, 0 for b = 0, and accumulate the calls and iterations.
The true |q - A x_e| can't drop below about eps |A| |x_e|, which for
small m and r2req lies above the CG target (m = 0.001 and r2req = 1e-24 on
a 4^3x8 lattice): rounding x_e by d changes b - M x by A d/(4m) through x_o,
while the correction y is small and so is its rounding.  The restarts in
double of the mixed precision solver also stop once the residual no longer
drops.  A solve out of iterations stops there and records the residual it
reached.

The September 30 production change, ac1a9f3b, caches the single precision
operator for a public solve, keeps the inner solve approximate, and refines
against the residual of the requested equation.  stagGpu already keeps
both precisions in its StagGpu objects: HMC refreshes both with setLinks
after smearing.  Its mixed solveM now takes one single precision CG per
correction and checks b - M x in double, including the reconstructed odd
sites.  In several systems, those corrections, reconstructions and full
residuals share the links and halo exchanges of up to nBatch systems.
The slot takes another system after its full residual meets the request,
exhausts its iterations, or stops dropping.

For both mixed solveEE and solveM, an inner CG with source q stops at
max(requested squared residual, max(r2in, epsilon(float32)) |q|^2).
r2in defaults to 1e-6.  The GPU programs accept `-r2in`; zero selects the
FP32 epsilon floor, about 1.19e-7.  The final tolerance still applies to the true
residual in double.  solveM does this refinement on M directly, which
avoids spending inner iterations refining A at its rounding floor before
checking the full equation.

The CG programs bestagcg, bestagm and bestagres also accept `-recon64:1`.
For compressed FP32 links, this evaluates the cross product and the sign
or determinant factor in FP64, rounds the reconstructed row to FP32,
and applies the matrix to the vector in FP32. The default is zero.
`newStagGpu(..., recon64 = true)` selects the same mode through the API.
Storage stays at 12 or 14 FP32 reals per link. FP64 operators and 18-real
links use their usual kernels. The host selects the kernel specialization;
scalar and batched hops use the same reconstruction for local and halo
links. Batched hops share each reconstructed matrix across their RHS.
The optional work counters classify stencils by field precision, so this
mode still contributes to the FP32 site count.

For example, compare `-recon64:0` and `-recon64:1` with the same arguments:

```
bin/bestagcg -mixed:1 -reals:12 -recon64:1 -r2in:1e-6
bin/bestagm -mixed:1 -reals:14 -recon64:1 -nb:3 -r2in:1e-6
```

The CPU production solver and GPU backend also have an optional conventional
mixed CG with reliable residual updates. The following controls can be set
independently for each backend. CPU controls are initialized in SolverParams;
GPU controls are stored in the FP32 StagGpu object's `ctrl` record.

| CPU argument | GPU argument | Meaning |
|---|---|---|
| `-cpuCg` | `-gpuCg` | CPU: -1 automatic, 0 CG corrections, 1 reliable CG; GPU: 0 one-reduction CG, 1 reliable CG, 2 widened correction accumulation |
| `-cpuDelta` | `-gpuDelta` | Reliable-update factor in the residual norm; 0 disables the norm trigger |
| `-cpuPeriod` | `-gpuPeriod` | Maximum iterations between precise residual checks, default 1000; 0 uses only the norm trigger |
| `-cpuFloor` | `-gpuFloor` | Factor for the inner arithmetic error estimate; default 1, 0 disables |
| `-cpuAcc64` | `-gpuAcc64` | Partial-solution precision: 1 uses FP64, 0 uses FP32; reliable CG also widens fused direction arithmetic when 1 |
| `-cpuBeta` | `-gpuBeta` | 1 uses the modified residual-difference numerator for beta; 0 uses the residual norm squared |
| `-cpuKeep` | `-gpuKeep` | Retain and reorthogonalize the search direction at reliable updates |
| `-cpuMaxInc` | `-gpuMaxInc` | Allowed consecutive non-improving reliable residuals |
| `-cpuMaxTotal` | `-gpuMaxTotal` | Allowed total non-improving reliable residuals |
| `-cpuR2in` | `-r2in` | Minimum relative squared target for an inner correction |

The reliable solver applies A to the search direction and periodically
recomputes b-Ax with the precise operator. A candidate convergence always
triggers a precise check. The partial solution is accumulated into the
precise total before replacing the residual; direction retention projects
the old direction orthogonal to the new residual before continuing CG.
The modified beta uses Re(r_new^+ (r_new-r_old))/|r_old|^2, with the usual
norm numerator if it is negative. Dot-product reductions use FP64.

The CPU default is `cpuCg=-1`, `cpuAcc64=0`, `cpuBeta=0`, `cpuDelta=0.01`,
`cpuPeriod=1000` and `cpuKeep=1`. These choices apply when mixed precision
is enabled (`sloppySolve=1` in production CPU programs, `mixed=1` in
bestagm). Automatic selection follows the reconstruction
choice: solveReconL uses reliable CG, while solveReconR and standalone
solveEE/solveOO use CG corrections. The choice is made
from the current source's parity norms, using the same criterion as the
reconstruction. `cpuCg=0` selects CG corrections; precise residual checks and recovery
apply to both algorithms.

The GPU numerical defaults are `gpuCg=0`, `r2in=1e-6`, `reconFma=0`
and `recon64=0`. For the optional conventional mode, the controls default
to `gpuDelta=0.1`, `gpuPeriod=1000`, `gpuAcc64=1`, `gpuBeta=1` and
`gpuKeep=1`. The defaults are `cpuFloor=1` and `gpuFloor=1`; `cpuR2in=0` leaves the
CPU inner target unrestricted. Both allow one consecutive and ten total
unexpected increases of the precise residual before ending a reliable
correction. An increase counts only when the recursive residual predicted
progress; periodic checks can observe ordinary nonmonotone CG residuals. The complete equation's residual
and total iteration budget still decide whether the solve succeeds.

GPU mode 2 uses the one-reduction recurrence, coupled halo updates and
system queue.
With `gpuAcc64=1`, it accumulates each correction in FP64, using the FP64
alpha coefficient for that accumulation. The residual, direction and
cached A-direction recurrence use their existing FP32 arithmetic. It
uses outer corrections with one global reduction per inner iteration.
`gpuDelta`, `gpuPeriod`, `gpuBeta`, `gpuKeep` and the residual-increase
limits select behavior in GPU mode 1.

`delta` acts on the norm: 0.1 represents a tenfold reduction during
monotonic convergence. The `r2in` controls act on squared norms. The
ordinary CG corrections clamp an inner target to FP32 epsilon. Conventional
mixed CG can check a smaller target with reliable updates, so an inner
floor of zero permits the requested target. In both GPU solveEE and solveM,
`r2in` selects the accuracy of each even-system correction. `gpuFloor`
scales the estimated outer-precision arithmetic error
`epsilon(FP64) * (norm(source) + Anorm * norm(correction))`, with
`Anorm = norm(A source) / norm(source)` estimated from the first sloppy
operator application. Setting `gpuFloor=0` disables this estimate.
Its norms share existing paired reductions. This floor ends an inner
correction; it does not relax the requested final equation tolerance.

GPU solvers evaluate trial corrections with the precise operator and accept
only a smaller residual (or one meeting the target). Failed mixed corrections
retry from the accepted solution in FP64 on the device. All attempts share
the original iteration budget; exhaustion returns the best accepted solution
and its actual residual. The batched paths recover each failing system using
the separate scalar workspace. Full solves check both parities after odd-site
reconstruction. Compressed links in the precise operator still define that
operator; comparisons with uncompressed CPU links require the documented
roundoff allowance.
The period bound ensures that a stalled recursive residual still reaches
the precise residual's progress checks. Disable both `Delta` and `Period`
to check the precise residual only at candidate convergence (and, on the
GPU, the iteration limit).

Conventional GPU CG shares links and halo exchanges within groups of up
to nBatch systems. Each slot has its own reliable-update history and
iteration budget. Groups finish before taking the next group; the
one-reduction solver replaces individual completed systems in its queue.

`-reconFma:1` uses explicit FP32 fused operations for the cross product. The phase multiplication
of a 14-real link also uses FP32 fused operations. Both link orientations
and halo paths use the same reconstructed matrix. `recon64` takes
precedence when both reconstruction options are enabled. The HMC GPU
programs expose both reconstruction options as well.

## Configuration examples

As in [INSTALL.md](../INSTALL.md#configuration-examples), `<configure>` is
the configure script, including path, in the QEX source directory; it
writes the options into `qexconfig.nims` of the build directory.  Each build
directory needs its own `nimcache` (the default, `nimcache` in the build
directory): builds sharing one corrupt each other.  The Nim phase of the
large programs (bestagcg, bestagres, staghmcgpu_sh, eightFlavorSMGgpu)
takes 13-26 minutes and 30-40 GB of memory, with any backend.

`-march=native` fits a build on a node of the kind that runs the program.
On a login node with another CPU, name the CPU of the compute nodes
(`sapphirerapids` for the Xeon 8468 and Xeon CPU Max, `graniterapids` for
the Xeon 6, `znver4` for the EPYC 9654 and MI300A, `znver2` for the EPYC
7742) and set `simd:` to match.  All tested builds use `vlen:8`.  Where
QMP and QIO are static libraries built without `-fPIE`, add `-no-pie` to
the link flags.

### OpenMP offload, Intel PVC (Aurora, Sunspot)

```
<configure> \
  qmpdir:"$HOME/lqcd/install/qmp" \
  qiodir:"$HOME/lqcd/install/qio" \
  cctype:"clang" \
  cc:"mpicc" \
  cflagsspeed:"-O3 -march=native" \
  ldflags:"-g -O3 -Xs '-device pvc' -ftarget-register-alloc-mode=pvc:large -ldl" \
  cpp:"mpicxx" \
  cppflagsspeed:"-O3 -march=native" \
  ldppflags:"-g -O3 -Xs '-device pvc' -ftarget-register-alloc-mode=pvc:large -ldl" \
  env:"OMPFLAG=-fiopenmp -fopenmp-targets=spir64_gen" \
  nimargs:'@["-d:Backend=OpenMP"]'
```

with the oneapi and mpich modules, whose `mpicc` runs icx.  Nim generates C;
the device code is compiled ahead of time at the link, which needs `-O3`:
icx treats `-g` alone as `-O0` for the device code.

### SYCL, Intel PVC

```
<configure> \
  qmpdir:"$HOME/lqcd/install/qmp" \
  qiodir:"$HOME/lqcd/install/qio" \
  cctype:"clang" \
  ccdef:"cpp" \
  cc:"mpicc" \
  cflagsspeed:"-O3 -march=native" \
  cpp:"mpicxx" \
  cppflagsalways:"-g -fsycl -fsycl-targets=spir64_gen" \
  cppflagsspeed:"-O3 -march=native" \
  ldppflags:"-g -O3 -fsycl -fsycl-targets=spir64_gen -Xsycl-target-backend '-device pvc' -ftarget-register-alloc-mode=pvc:large -ldl" \
  env:"OMPFLAG=-fiopenmp" \
  nimargs:'@["-d:Backend=SYCL"]'
```

Nim generates C++, compiled by icpx.  For both Intel backends
`-ftarget-register-alloc-mode=pvc:auto` lets the compiler choose 128 or 256
registers per kernel: mixed precision staghmcgpu_sh runs 5-7% faster,
double precision within 2%.

### CUDA with clang, NVIDIA H100, B200

```
<configure> \
  qmpdir:"$HOME/lqcd/install/qmp" \
  qiodir:"$HOME/lqcd/install/qio" \
  cctype:"clang" \
  ccdef:"cpp" \
  cc:"mpicc" \
  env:"OMPI_CC=clang" \
  cflagsspeed:"-O3 -march=native" \
  cpp:"mpicxx" \
  env:"OMPI_CXX=clang++" \
  cppflagsalways:"-g -x cuda --cuda-gpu-arch=sm_90 --cuda-path=$CUDA_HOME -Xarch_device -mllvm=-disable-machine-sink" \
  cppflagsspeed:"-O3 -march=native" \
  ldppflags:"-g -L$CUDA_HOME/lib64 -Wl,-rpath,$CUDA_HOME/lib64 -lcudart -ldl" \
  simd:"SSE,AVX,AVX512" vlen:8 \
  nimargs:'@["-d:Backend=CUDA"]'
```

with `CUDA_HOME` a CUDA toolkit that clang supports (clang 22 with CUDA
12.9) and OpenMPI wrappers (`OMPI_CC`, `OMPI_CXX`; MPICH takes `MPICH_CC`
and `MPICH_CXX`).  `--cuda-gpu-arch` names the GPU: `sm_90` for H100,
`sm_100` for B200.  backend/cuda/nimbase.h, found first through `-iquote`,
defines `QEX_HD` and the qexFor kernel template.  Only programs that import
backend/accel find it, so CPU programs, such as tests/base/trngfieldv, don't
build with these flags, which compile every file as CUDA.

- The device flag `-disable-machine-sink` and the `[[clang::always_inline]]`
  call in qexFor belong together: LLVM's machine sinking moves the
  arithmetic of a hop past the loads of the later directions, and LLVM leaves
  large lambdas uninlined, so the batched hop of 3 systems spilled (255
  registers, up to 4 KB of stack; either change alone does not help).  With
  both, 164 registers and no stack, as with nvcc, and the default nBatch 3
  fits.
- To link on a node without the NVIDIA driver, add
  `-L$CUDA_HOME/lib64/stubs` (the stub libcuda; the driver's at run time).

### CUDA with nvcc, NVIDIA H100

```
<configure> \
  qmpdir:"$HOME/lqcd/install/qmp" \
  qiodir:"$HOME/lqcd/install/qio" \
  cctype:"clang" \
  ccdef:"cpp" \
  cc:"mpicc" \
  env:"OMPI_CC=gcc-14" \
  cflagsspeed:"-O3 -march=native" \
  cpp:"$PWD/qex/build/nvcc.sh" \
  env:"OMPI_CXX=g++-14" \
  env:"QEX_HOSTCXX=g++-14" \
  env:"CUDA_HOME=$CUDA_HOME" \
  env:"CPATH=$MPI_HOME/include" \
  cppflagsalways:"-g -x cuda --cuda-gpu-arch=sm_90" \
  cppflagsspeed:"-O3 -march=native" \
  ldpp:"mpicxx" \
  ldppflags:"-g -L$CUDA_HOME/lib64 -Wl,-rpath,$CUDA_HOME/lib64 -lcudart -ldl" \
  simd:"SSE,AVX,AVX512" vlen:8 \
  nimargs:'@["-d:Backend=CUDA"]'
```

run in the build directory, with `CUDA_HOME` and `MPI_HOME` set (`$PWD/qex`
is the link that configure makes to the source).  build/nvcc.sh turns the
compile commands Nim writes for clang into nvcc commands (`-x cu
-arch=sm_XY --extended-lambda`, the gcc options through `-Xcompiler`) for
the files with kernels and gives the other files to the host compiler
`QEX_HOSTCXX`.  `CPATH` holds the MPI headers for the host compiler and
nvcc's preprocessing; `ldpp` links with the MPI wrapper.  On one H100 the
nvcc build runs within 1-6% of the clang build.

### HIP, AMD MI300X, MI300A

```
<configure> \
  qmpdir:"$HOME/lqcd/install/qmp" \
  qiodir:"$HOME/lqcd/install/qio" \
  ccdef:"cpp" \
  cc:"mpicc" \
  cflagsspeed:"-O3 -march=native" \
  cpp:"mpicxx" \
  cppflagsalways:"-g -x hip --offload-arch=gfx942" \
  cppflagsspeed:"-O3 -march=native" \
  ldppflags:"-g -ldl --hip-link --offload-arch=gfx942" \
  simd:"SSE,AVX,AVX512" vlen:8 \
  env:"STATIC_UNROLL=1" \
  nimargs:'@["-d:Backend=HIP"]'
```

with MPI wrappers running amdclang: on Tuolumne the modules
rocmcc/10.0-magic, rocm/10.0 and cray-mpich/9.1.0; with OpenMPI,
`env:"OMPI_CC=amdclang" env:"OMPI_CXX=amdclang++"` and ROCm's bin in
`PATH`.  backend/hip/nimbase.h includes hip_runtime.h (for
`__launch_bounds__`) and defines `QEX_HD` and the qexFor kernel template;
as with CUDA, CPU programs don't build with these flags.  All modules must
be C++ (`ccdef:"cpp"`): in C, `__host__ __device__` does not compile.
gfx942 is both MI300X and MI300A.

## Runs

- PVC: one rank per tile (gpu_tile_compact.sh sets each rank's
  `ZE_AFFINITY_MASK`), 8 cores per rank,
  `MPIR_CVAR_CH4_IPC_GPU_P2P_THRESHOLD=0`, `MPIR_CVAR_CH4_OFI_ENABLE_HMEM=1`
  off the node, `OMP_STACKSIZE=256M` for large local volumes.
- CUDA, HIP and SYCL ranks take GPU `myRank mod` the number of GPUs they see
  (gpuInit); OpenMP ranks take the default device.
- JLSE: one rank per GPU on the cores of the GPU's NUMA domain, 16 threads
  per rank (12 on blackwell00, where two GPUs share each 24-core domain),
  `OMP_STACKSIZE=256M`.  Its OpenMPI 4.1.1 is not GPU aware; the halos
  between GPUs of a node go through CUDA or HIP IPC.
- Tuolumne: `flux run -N1 -n4 -c24` with `MPICH_GPU_SUPPORT_ENABLED=1`; from
  an ssh shell on a node of the allocation, set `FLUX_URI` to the job's
  instance (`flux uri JOBID`), else `flux run` submits a new job.

## Programs and tests

- backend/examples/bestream: STREAM copy, scale, add and triad of gpuFor
  kernels, the bandwidth ceiling the kernels can reach.
- backend/examples/bestagcg: the GPU CG against the CPU solver, GF/s and
  GB/s by a byte count, `-nb:k` k systems at once.
- backend/examples/bestagm: complete scalar and batched M solves against
  the CPU operator, `-nb`, `-masses`, `-full`, `-mixed`, `-r2req` and
  `-r2in`.  `-d:cpuOnly` builds the production reference;
  `-d:stagWorkCount` records local stencil work, CG starts and reductions.
- backend/examples/bestagres: the stopping rules of solveM and solveEE
  against stag.D, stagD2ee and the production solves: one system and several, double and mixed
  precision, fixed and atomic dot products, sources on all or on the even
  sites, masses of both signs, zero sources, too few iterations, reused
  statistics and inner tolerances below and above single precision epsilon;
  `-ipc` and `-split` as bestagcg; exits with 1 if a check fails.
  Its raw inner CG checks print CPU/GPU/batch iteration counts, recursive
  residuals, and residuals recomputed with the CPU operator in both
  precisions.  The double check promotes the rounded single precision
  links and source to double.  `-inner:0` skips those checks and `-cpu:0`
  skips the comparisons with the production double/mixed solves.
- backend/examples/berng: the generators of RNGFieldV and rng/rngGpu against
  the host RNG fields; exits with 1 if a check fails.
- tests/base/taccel: the order of queued kernels, sumFixed and gpuSum at
  the sizes where their levels change, and the neighbors GpuHaloEx brings,
  on 1, 2 and 4 ranks.
- tests/base/tgaugegpu: gauge/gaugeGpu against actionA and forceA, and
  ValueError for rect and pgm coefficients.
- tests/base/tstaglinks: adjoint consistency of 12/14/18-real caches after
  construction and gauge updates, in both cache modes; promotes FP32
  cache entries before the check to expose representation errors.
- tests/base/tstagzero: CPU zero-source solutions and finite statistics,
  including warm starts and multiple masses.
- tests/base/tstagrecon: FP64 reconstruction of FP32 links against an
  independent complex-arithmetic reference with 18-real storage; exercises
  both cache modes, mode changes and the unchanged FP64/18-real paths.
- backend/examples/bernt: the times of the rngGpu draw kernels, with 32-,
  64- and 128-bit generator loads.
- examples/staghmcgpu_sh: staghmc_sh on the GPU; `tests/extra/tstaghmc_sh/run`
  with `RUNJOB` set to the launcher.
- prod/lsd/eightFlavorSMGgpu: eightFlavorSMG on the GPU;
  `tests/extra/teightFlavorSMG/run` (needs python 3.8 or later).

## Measured

One GPU (one PVC tile), devel of September 30, the builds above.

| | H100 | B200 | MI300X | MI300A | PVC tile, OpenMP / SYCL |
|---|---|---|---|---|---|
| HBM peak (TB/s) | 3.35 | 8.0 | 5.3 | 5.3 | 1.64 |
| bestream triad (TB/s) | 3.10 | 5.84 | 4.09 | 3.35 | 0.98 / 0.98 |
| bestagcg 32^4, one system (TB/s) | 2.65 | 5.28 | 3.21 | 2.61-2.62 | 0.82 / 0.83 |
| staghmcgpu_sh 24^4, double / mixed (s per trajectory) | 0.887-0.889 / 0.679-0.680 | 0.625-0.626 / 0.528-0.529 | 0.954-0.961 / 0.770-0.772 | 1.064-1.069 / 0.866-0.885 | 2.03-2.04 / 1.68-1.69, 1.98-2.00 / 1.61-1.62 |
| eightFlavorSMGgpu 16^3x32 (s per trajectory) | 2.72-2.73 | 2.13-2.14 | 2.79-2.81 | 3.44-3.46 | 5.04-5.07 / 4.75-4.77 |

Against the code before the stopping rules, in the same jobs (runs in
the order old, new, new, old), bestagcg 32^4 takes 0.1-0.6% longer for one
system, whose CG checks its true residual at the end, and 0.9-1.6% for
three systems per hop, which finish one by one (on MI300X and MI300A those
runs spread 4%).  staghmcgpu_sh takes 2.5-3.8% longer per trajectory at
24^4 (3.3-6.3% at 16^4): its action solves, whose sources have odd sites,
take 13-15% more iterations to meet the request on M x = b, and every solve
checks b - M x.  eightFlavorSMGgpu, whose solves take the same iterations as
before, takes 0.4-2.5% longer.

## Initial CPU validation, October 1, 2026

gpu merged devel ac1a9f3b in 64caa668.  The implementation above was
tested on an Apple M4 Max with 128 GB of memory, Nim 2.3.1
(e754af2083464d791a2584a5c0239bb4b8887ca9), clang 23.1.0, QMP SINGLE,
`vlen:4`, no explicit SIMD intrinsics, and C flags
`-O3 -march=native -fno-strict-aliasing`.  Assertions were enabled in
the release builds.  The production baseline was a clean archive of
devel ac1a9f3b, built with the same settings.

With `OMP_NUM_THREADS=1`, each of these bestagres configurations reported
470 successful checks, including production solver comparisons and
raw inner CGs.  The audit below qualifies what those checks establish:

```
bin/bestagres
bin/bestagres -lat:4,4,4,8 -reals:12 -split:1
bin/bestagres -lat:4,4,4,8 -reals:18 -split:1
bin/bestagres -lat:8,8,8,8 -reals:18 -split:0
```

taccel and tgaugegpu passed.  All three tstaghmc_sh cases passed with the
production CPU program, and with staghmcgpu_sh in double and mixed
precision (the GPU version's default links of 14 reals).  Both production
eightFlavorSMG and eightFlavorSMGgpu passed the complete teightFlavorSMG
suite, including reversibility, checkpoint/resume, rejected updates, and
invalid checkpoint/configuration checks.  A separate mixed
eightFlavorSMGgpu trajectory, with both action and force tolerances set to
1e-24, matched clean devel's Hamiltonians and gauge observables at rtol
5e-11, with atol 1e-9 for the kinetic action.  Both recorded reverse dH0
= -2.91e-11, below the suite's 1e-8 bound.

The added loose-inner-tolerance cases caught a progress check against the
zero guess in the batched mixed solver.  A first reconstructed solution
can increase |b - M x| even while reducing the even-system residual.
The progress check now starts after that first reconstructed solution,
as correctM does for one system; later corrections must reduce the full
residual.

### Inner single precision CG

With the same gauge and source, r2req = epsilon(float32) = 1.1920929e-7,
links of 18 reals and fixed sums, the CPU and GPU raw inner solves gave:

| Lattice | Mass | CPU / GPU / batch iterations | True r2/b2, CPU / GPU |
|---|---|---|---|
| 4^4 | 0.1 | 108 / 108 / 108 | 1.190e-7 / 1.176e-7 |
| 4^4 | 0.01 | 628 / 625 / 625 | 4.720e-7 / 4.261e-7 |
| 4^4 | 0.001 | 643 / 642 / 642 | 1.421e-6 / 1.790e-6 |
| 4^3 x 8 | 0.001 | 1319 / 1322 / 1322 | 3.852e-4 / 5.536e-4 |
| 8^4 | 0.001 | 9969 / 9958 / 9958 | 2.168e-2 / 2.237e-2 |

All recursive residuals met the requested inner target.  The true residuals
in this table were recomputed in double with the rounded single precision
links and source.  Their drift for light masses occurs in both codes;
the final mixed solves check and correct the double precision equation.
The GPU single and batched solutions matched bit for bit in these CPU
backend runs, with both fixed and atomic sums.

The compressed case deserves separate attention: at 4^3 x 8, m = 0.001
and the same epsilon target, 12 reals gave 2453 GPU/batch iterations
versus 1319 CPU iterations, with true r2/b2 of 1.360e-3 versus 3.852e-4.
With 18 reals the GPU took 1322 iterations.  The GPU reconstructs the
third row in single precision; the CPU stores that row rounded from
double.  All converged cases in both configurations met their double
precision requests, including 1e-24.

A controlled repeat identified the adjoint inconsistency of separately
compressed forward and backward links.  Reconstructing U and U^+ from
their own rounded first two rows gives different matrices:
R(U^+) != R(U)^+ in single precision.  In this case their largest element
discrepancy was 1.23e-7.  The linear even operator assembled in double
from those reconstructed single precision links had
|A - A^+|_F / |A|_F = 3.80e-8.  Using the same forward cache for both
directions removed the link discrepancy and reduced that operator
defect to about 1.5e-17, while keeping the 12-real representation.

| GPU links | Backward hop | Inner iterations |
|---|---|---|
| 18 reals | separate adjoint cache | 1322 |
| 18 reals | adjoint of forward cache | 1322 |
| 12 reals | separate adjoint cache | 2453 |
| 12 reals | adjoint of forward cache | 1321 |

The 12-real separate-cache recurrence first reached r2/b2 = 1.205e-7 at
iteration 1314, just above the 1.192e-7 stopping threshold.  It then rose
to 2.317e-5 at iteration 1449 and crossed the threshold at iteration
2453.  The condition number was about 1.68e6 in each case.  The extra
iterations therefore accompany the loss of adjoint consistency and the
subsequent recurrence drift in this sensitive system.  Its final true
residual was 1.360e-3 with separate caches versus 4.557e-4 with the
forward cache, evaluated by the CPU operator with rounded full links.

The existing `newStagGpu(g, float32, 12, fwd = 1)` path forms local
backward hops from the same reconstructed forward links.  These
diagnostics used one rank; remote links in lh still hold independently
compressed adjoints, so the MPI boundary case needs the same consistency
analysis.

### CPU timings against production

Times below are milliseconds per complete solveEE call, including its
setup and transfers.  The gauge, layout and random seed matched; links
had 18 reals, lanes were `[2,2,1,1]`, and r2req was 1e-16.  Four samples
followed an initial warmup for each GPU path.  Production ran before and
after it, with four samples following each warmup; the table gives their
medians.  All large compilations and other tests had finished before the
timed runs.  The last two columns use six OpenMP threads for production.

| Lattice | Mass | CPU double | GPU CPU backend double | CPU mixed | GPU CPU backend mixed | CPU double, 6 threads | CPU mixed, 6 threads |
|---|---|---|---|---|---|---|---|
| 8^4 | 0.1 | 20.47 | 40.93 | 10.77 | 65.08 | 5.18 | 4.23 |
| 8^4 | 0.01 | 202.01 | 407.52 | 102.89 | 642.11 | 45.64 | 27.20 |
| 16^4 | 0.1 | 360.41 | 703.70 | 187.32 | 1102.47 | 67.75 | 42.90 |
| 16^4 | 0.01 | 3631.44 | 7077.63 | 1804.70 | 11000.17 | 667.64 | 372.48 |

For these cases production is about twice as fast in double and six
times as fast in mixed precision on one thread.  The GPU CPU backend
runs its kernels as serial host loops; production also scales across
OpenMP threads. These initial measurements describe CPU execution. The
subsequent PVC study below covers device compilation, MPI and GPU timing.

The GPU benchmark command for each lattice and mass was:

```
OMP_NUM_THREADS=1 bin/bestagcg -lat:16,16,16,16 -rg:1,1,1,1 \
  -ig:2,2,1,1 -mass:0.01 -r2req:1e-16 -maxits:50000 \
  -ncpu:1 -ngpu:5 -mixed:1 -sloppySolve:0 -reals:18 -fixed:1
```

The clean devel benchmark used the same initialization and production
solveEE call from bestagcg, with `-ncpu:5` and `-sloppySolve:0` or `1`.
bestagcg now prints the elapsed time of each complete CPU, GPU double,
and GPU mixed solve so that this comparison includes their setup.

## Audit before implementation

The final mixed solutions and HMC comparisons above provided useful CPU
backend evidence.  The pass count alone did not establish inner FP32
accuracy, device correctness, MPI correctness, or a performance improvement
from the solver changes.  The initial study left the compressed-link
defect in the code. The following findings describe that pre-repair state;
the implementation and further validation are recorded below.

### Findings

1. **The adjoint defect is confirmed, and its repair must cover every
   backward path.** A stronger control expanded the independently
   reconstructed FP32 forward/backward matrices into 18-real storage.
   The ordinary 18-real kernel reproduced the 2453 iterations and the
   same solution bits.  Holding those forward matrices fixed and replacing
   only the backward matrices with their exact adjoints gave 1321
   iterations and the same solution bits as the 12-real forward-cache
   run.  This removes the kernel/representation confound in the earlier
   comparison.  It establishes the cause for this fixture, rather than a
   universal iteration penalty.  The defect predates this session.
   `fwd = 1` repairs local pairing only; lh, constructor packing, setLinks,
   dslash and dslashB all participate in backward-link handling.  The
   14-real representation also independently reconstructs adjoints and
   needs coverage with nontrivial unitary determinants.

2. **The inner accuracy assertion is too weak.** checkInner requires the
   recursive residual to meet the target, then bounds solution differences
   by the measured true residuals.  For a Hermitian positive reference A,
   `A(xg-xc) = rc-rg` already implies
   `16 m^4 |xg-xc|^2 <= 2 (|rc|^2 + |rg|^2)`.  This is a consistency
   identity and can hold for a bad solution with a large true residual.
   For example, a zeroed downloaded solution can satisfy this bound if
   its stale recursive residual still reports convergence.  The tests
   therefore cannot call a raw solve accurate solely because this check
   passes.  They also use identical sources and masses in all raw batch
   slots and download only slot zero; that cannot detect slot permutation
   or establish accuracy of every returned slot.  The outer batch tests
   do use different systems and check each final solution.

3. **The residual references need explicit names.** The main inner table
   evaluates the CPU operator built from all 18 original link reals
   rounded to FP32, then promoted to FP64.  A compressed GPU operator
   represents different matrices.  Its residual against that common
   reference combines representation error with iteration error.  The
   matrix diagnostic subsequently measures a second residual against
   decoded GPU link matrices, but this linear model still excludes
   rounding inside FP32 stencil arithmetic.  Future checks should report
   recursive, effective-operator, and original FP64 problem residuals
   separately.  The observed large gap remains real: the 2453-iteration
   result has r2/b2 about 1.36e-3 against both matrix references.

4. **The new default tolerance has not been justified by a performance
   study.** The old GPU r2in default was 1e-6.  Changing it to FP32 epsilon
   tightens the inner squared-residual target by 8.39, whereas the upstream
   CPU change loosened its earlier floor.  In the sensitive compressed
   fixture the raw solve takes 1310 iterations at 1e-6 and 2453 at epsilon.
   That does not establish the total cost of a final mixed solve.  The
   timing table measures solveEE against production CPU, and does not
   compare the old/new GPU solveM or isolate the default change.  New
   batched M refinement also changes FP64 work and collective counts;
   norm2EO now performs a reduction for each completed slot.  CG iteration
   counts alone cannot quantify its cost.  The chosen -O3 flags are stated
   and comparable between the programs; production tuning with -Ofast
   remains a separate validation configuration.

5. **The CPU reference has a zero-source statistics bug.** A direct run
   returned `calls=1`, `iterations=0`, `r2=nan` for both solveEE and solve,
   in double and mixed precision.  Both public wrappers record r2/b2 with
   b2=0.  The solveEE behavior is introduced by the merged wrapper; the
   full solve already had this behavior.  bestagres skips the production
   comparison for zero sources, so its GPU zero-source checks did not
   expose the reference problem.

6. **Coverage is narrower than a device regression study.** QMP SINGLE
   gives no remote neighbors, even with split=1, and CPU gpuForAsync is
   synchronous.  These runs cannot test MPI halos, IPC, overlapping
   kernels, device capture rules or register pressure.  The mixed
   eightFlavorSMG comparison was one short cold-start trajectory with
   f_tol tightened from 1e-16 to 1e-24.  It is useful accuracy evidence;
   normal production tolerances, more seeds/configurations, and direct
   deltaH/force comparisons still need testing.  Several diagnostics and
   their metadata currently live only in temporary build directories.

### Implementation order and acceptance criteria

1. **Make the failures reproducible in the repository.** Preserve an
   immutable source/configuration baseline.  Add a focused link test that
   fails for the current independently compressed adjoints.  Cover
   12/14/18 reals, FP32/FP64, both cache modes, boundary phases, constructor
   and repeated setLinks updates.  Compare decoded cached data, bilinear
   adjoint identities, and actual stencil output against an independent
   reference.  Store geometry, lanes, ranks, seeds, mass, tolerance,
   compiler flags and cache orientation with results.  Keep iteration
   counts as diagnostics rather than portable equality assertions.

2. **Use one orientation for cached links.** Keep canonical forward link
   data in lf.  Store the neighbor's same forward-oriented payload in lb
   and lh; reconstruct it and use hopA/hopAm for every backward hop.
   Update newStagGpu, setLinks, scalar and batched kernels, and storage
   documentation together.  Preserve the existing cache sizes and
   locality.  Quantization and determinant handling must yield the same
   payload as the owning site's forward cache, including across rank
   boundaries.  Acceptance: matching link/adjoint pairs, a Hermitian
   decoded even operator, and correct results for fresh and updated
   SU(3) and U(3) fields.  The old defect fixture must fail before and
   pass after this change without weakening the tolerance.

3. **Repair the numerical checks and statistics.** Give raw inner solves
   separate termination and accuracy diagnostics.  Require independently
   measured residual reduction and a justified accuracy allowance on
   fixtures with known conditioning; retain the true FP64 target for
   final mixed solves.  Use distinct sources/masses in raw batch slots,
   verify every returned vector, and exercise partial batches and slot
   reuse.  A deliberately zeroed or permuted output must fail the checks.
   Check full, even-only and odd-only sources, zero sources, iteration
   exhaustion and correction stagnation.  Fix public CPU zero-source
   behavior/statistics in a separate patch, including warm-start handling,
   and compare it with GPU zero-source behavior instead of skipping it.

4. **Choose the stopping policy from measured final solves.** Keep the
   precision floor as a lower bound and restore 1e-6 as the provisional GPU
   default until the tighter default is supported by measurements.
   Compare the original GPU implementation, the cache fix alone, the new
   M refinement at 1e-6, and the epsilon setting at the same final accuracy.
   Benchmark solveEE, full solveM, heterogeneous batches and HMC.  Record
   wall time, inner iterations, refinements, FP32/FP64 stencil counts,
   reductions, halo traffic and setup costs.  Repeat in alternating order
   and report the sample spread.  Keep the validated algorithm/statistics
   changes distinct from any optional default-tolerance change.

5. **Validate the communication and device paths before finalizing.**
   Run CPU checks first, then 2/4 MPI ranks with real remote neighbors,
   forced cache/split modes, and several layouts.  Check device builds and
   the existing supported GPU backend tests, with fixed and atomic sums
   and IPC enabled/disabled where supported.  Run HMC at normal and tight
   tolerances on more than one gauge state, checking force accuracy,
   deltaH, reversibility and checkpoint/resume.  Publish exact tested
   revisions, commands and logs with the resulting scope of validation.

## Implementation and current validation

lf, lb and lh now all store links in forward orientation.  Every backward
hop takes the adjoint after reconstruction, in the scalar and batched
kernels.  newStagGpu and setLinks use this representation for local and
remote links.  The cache sizes and the automatic choice of whether to
cache backward-neighbor links are preserved.

tstaglinks fails before the repair on the 12- and 14-real FP32 caches,
with normalized adjoint defects of 1.6e-10 to 6.8e-10 after promotion.
The repaired caches pass the same 1e-13 criterion at about 1e-17, including
repeated setLinks updates and nontrivial U(3) determinants.  The CPU
zero-source tests pass in both precisions with one and six threads.
Re-running the original 2453-iteration fixture with the repaired code
gives 1321 iterations for both 12-real cache modes, versus 1322 for both
18-real cache modes.  The 12-real modes return the same residual and
solution, with a recursive r2/b2 of 1.149e-7 and a residual of 4.557e-4
against the FP64 CPU operator made from rounded FP32 links.

The mixed GPU default is again 1e-6, bounded below by epsilon.  bestagcg,
bestagm, bestagres, staghmcgpu_sh and eightFlavorSMGgpu accept `-r2in`.
Their existing final-solve tolerances remain available.  bestagres also
accepts `-innerReq`, `-innerMax`, `-innerMatch`, `-r2req`, `-fwd` and `-seed`.
Its raw batch checks use different sources and masses, download every
slot, compare scalar and batched solutions, and verify that zeroed and
permuted outputs are rejected.  On the selected hard fixtures the default
innerMax=0.25 requires at least a factor-two reduction in the true residual
norm; it does not claim that the true FP32 residual reaches the recursive
target.  Well-conditioned fixtures use a tighter bound tied to the request.
The default innerMatch=1e-4 bounds the squared batch/scalar vector difference.

Recorded residuals are checked by applying the same operator afresh to
the returned solution.  Independent CPU residuals separately check final
accuracy.  Comparing the two residual estimates to 1% was too strict for
some 14-real solveEE cases at r2req=1e-20: both met the request, while the
different operator rounding changed the tiny residual by 1-5%.

The first PVC two-rank run exposed the corresponding cross-reference
stopping check: an atomic 14-real solveEE at m=0.001 gave 9.899e-21 with
its own operator and 1.037e-20 with the original CPU operator, for a
request of 1e-20. The test now explicitly requires its freshly recomputed
own-operator residual and recorded residual to meet the request. The
independent CPU comparison allows FP64 representation/evaluation roundoff:
`min(refFactor*req, (sqrt(req)+delta)^2)`, where
`delta = refRound*epsilon64*opNorm*||x||/||b||`. For these unitary-link
fixtures, `opNorm=|m|+4` for M and `4m^2+64` for A. Defaults refRound=128
and refFactor=4 give a conservative roundoff allowance capped at twice
the requested residual norm. The allowance is independent of measured
residuals and is printed for each check. `-refRound:0` restores a strict
CPU-reference comparison. The solver continues to enforce its original
stopping tolerance. bestagm retains its separate 1.01*r2req CPU accuracy
gate for every timing sample.

The final laptop regressions pass for 12 and 14 reals with six threads,
and for 18 reals with a second seed.  Source initialization now has an
explicit barrier before changing subsets or measuring their norms;
without it, the multithreaded test could normalize a residual by a partly
initialized source.  The three mixed staghmc regressions and the full
eightFlavorSMG regression suite also pass with six threads.

Cache consistency does not guarantee equal FP32 iteration counts.  With
the new distinct sources on 4^3 x 8, the first two compressed systems at
m=0.001 and 0.0012 take 1934 and 2209 iterations at epsilon, versus CPU
1318 and 1319 and uncompressed GPU 1313 and 1312.  The corresponding
effective-operator squared residuals are 7.07e-4 and 4.95e-4 despite
recursive residuals near 1e-7.  At the retained 1e-6 default those
compressed solves take 1312 and 1313 iterations.  All final mixed solves
meet their independently checked targets.  The cache repair removes the
identified adjoint inconsistency; FP32 recurrence and compressed-operator
sensitivity still require measurements of complete refined solves.

The final full-M laptop benchmark uses 16^4, three independent sources
with masses 0.05, 0.06 and 0.07, r2req=1e-16, r2in=1e-6, 14 reals,
fixed sums and one host thread. Each implementation runs four times,
then the order is reversed; each process's first sample is discarded.
All returned vectors meet the CPU residual check. Production is clean
devel plus the same bestagm driver built with cpuOnly. Timings cover the
complete solve call, with gauge/source setup and residual verification
outside the timed region.

| Precision | Production CPU median (range), s | GPU CPU backend median (range), s | Total iterations, CPU / GPU |
|---|---|---|---|
| Double | 1.832 (1.809–1.856) | 8.525 (8.389–8.649) | 1308 / 1308 |
| Mixed | 0.964 (0.961–0.968) | 7.098 (7.075–7.127) | 1340 / 1336 |

Hot-start HMC uses seed 314159265 for the parallel RNG, 271828182 for
the serial RNG, and two trajectories of test_input.xml with its usual
force tolerance 1e-16. On the laptop the production CPU reversal error
reaches 2.81e-8; GPU mixed reaches 1.01e-7 at r2in=1e-6 and 1.05e-7 at
epsilon. The largest CPU/GPU trajectory deltaH difference is 1.97e-6.
These trajectories reject, so identical final gauge observables would
not establish agreement of their integration. Tightening only the force
tolerance to 1e-24 gives CPU reversal errors below 1.46e-10, GPU mixed
below 2.92e-11, and CPU/GPU deltaH differences below 2.04e-10.

The device study therefore uses an explicit 1e-6 reversal budget and
1e-5 CPU/GPU deltaH budget for this normal-tolerance hot fixture, then
requires 1e-8 for both at force tolerance 1e-24. Initial Hamiltonians
must agree within 1e-8. These budgets are configurable and their values
are saved beside each result. The cold suite retains its existing limits.
Force RMS and variance agree at their printed precision when all forces
are sampled. Printed force Inf is unsuitable for a CPU/GPU comparison
across thread counts: the current CPU implementation sums thread maxima.

The study runner `solver_study.py` preserves commands, statuses, logs and
benchmark samples. Its scripts and local archives are stored outside the
repository at `/private/tmp/qex-gpu-20261001/gpu-study`.
The Sunspot study uses isolated source/build directories under
`~/W/qex_solver_fix_20261001`, with the existing OpenMP/SYCL presets,
allocation Catalyst and queue workq.  Its manifest records source hashes
for the original implementation, the cache repair alone, the complete
repair and clean devel CPU reference.

## PVC results, October 1–2, 2026

Sunspot jobs 12480027, 12480031 and 12480034 completed with exit status 0
in workq under Catalyst. The user's existing checkout and build directories
were preserved; the study is in `/home/xyjin/W/qex_solver_fix_20261001`.
Builds used Nim 2.3.1 at e754af20, Intel oneAPI 2026.1.0.20260617,
MPICH 5.0.0.aurora_test.87e2045, VLEN=8, `-O3 -march=native`, assertions,
and the existing PVC OpenMP/SYCL presets. Runtime uses eight host threads
per rank, the site's GPU tile mapping wrapper, and mandatory OpenMP
offload. Python 3.12.12 runs the test scripts.

Both backends passed the following, with fixed and atomic reductions
where the solver test selects them:

- One, two and four ranks on one node: link consistency, repeated link
  updates, halo primitives, gauge forces, zero-source CPU solves, raw
  single precision CG, full-M/EE solves, heterogeneous batches and slot
  reuse, iteration limits, double/mixed HMC, and the complete
  eightFlavorSMG checkpoint/configuration suite.
- Two ranks on separate nodes: the same link, halo, gauge and solver
  checks, including forced cache/split modes and IPC disabled.
- Hot-start mixed HMC at normal and tight force tolerances, compared with
  clean devel CPU on the same rank count and initial state.

There are 50 completed check commands plus eight off-node commands per
backend, containing 7,998 printed solver comparisons per backend. The
largest logged CPU-reference M residual is at the requested target to
the printed precision. The largest A-reference ratio is 1.037, the
14-real OpenMP case discussed above; the solver's own residual satisfies
the strict request. The maximum hot-start GPU reversal error is 4.51e-7
at force tolerance 1e-16, with CPU/GPU deltaH differences at most 2.01e-6.
At force tolerance 1e-24 those maxima are 2.92e-11 and 3.21e-10.

### Complete-solve costs

The main matrix uses 24^4, masses 0.05/0.06/0.07 for three systems,
14-real links, automatic cache selection, fixed sums, full or even-only
sources, one or three systems, one/two/four ranks, and inner tolerances
1e-6 and epsilon. Every run requests r2req=1e-16 and independently checks
every returned vector with the original CPU operator at 1.01*r2req.
There are 144 benchmark commands per backend, each with four solves;
the first is discarded. Original/cache/fixed order is reversed for the
second cycle, giving six retained samples per comparison.

The table gives full-M, three-system medians at the default inner
tolerance, in milliseconds. Parentheses give the range of six samples
for the original and repaired versions.

| Backend | Ranks / tiles | Original | Cache repair only | Complete repair |
|---|---:|---:|---:|---:|
| OpenMP | 1 | 199.94 (192.74–204.77) | 198.04 | 194.60 (187.16–199.43) |
| OpenMP | 2 | 129.31 (128.78–129.51) | 130.17 | 128.19 (127.90–128.96) |
| OpenMP | 4 | 139.66 (139.52–140.26) | 138.69 | 136.22 (135.98–136.57) |
| SYCL | 1 | 202.13 (187.51–204.63) | 193.87 | 190.18 (174.99–192.61) |
| SYCL | 2 | 121.66 (121.38–122.45) | 121.51 | 120.43 (119.74–120.74) |
| SYCL | 4 | 131.67 (131.56–131.90) | 131.74 | 129.81 (129.67–130.02) |

Across all main cases the repaired/original median-time ratio ranges
from 0.973 to 0.997 for OpenMP and 0.941 to 1.004 for SYCL. These are
modest changes with overlapping ranges in several cases, especially at
one tile. Four tiles do not improve this lattice over two tiles.

The representative full-M batch changes from 1,371 to 1,340 total inner
iterations. At four ranks it changes stencil calls from 1,110 to 1,074
and reduction calls from 557 to 548, while FP64 stencil sites per rank
increase from 1,492,992 to 1,990,656. The final CPU squared residual is
about 4.89e-17 originally and 9.88e-17 after repair. Both meet the same
request; some savings come from avoiding unnecessary extra accuracy.
These counters describe stencil calls/sites and global-reduction calls,
not all device kernels or measured halo bytes.

The clean production CPU reference, one rank with eight host threads,
takes 7.15–7.20 seconds in double and 2.79–2.80 seconds in mixed precision
for the same full-M batch. These are complete API calls; GPU gauge/source
setup and output validation are outside the timed region, while setup
performed inside the production CPU solve is included.

### Light-mass sensitivity remains measurable

A separate 4^3 x 8 study at mass 0.001, 12 reals and forced separate
caches uses the same accuracy gate and six retained samples. Correcting
the link representation changes FP32 convergence and can increase cost.
For the three-system even-source batch:

| Backend | Inner r2 | Original, ms | Cache repair only, ms | Complete repair, ms |
|---|---:|---:|---:|---:|
| OpenMP | 1e-6 | 243.62 | 264.60 | 258.55 |
| OpenMP | epsilon | 219.23 | 303.39 | 297.04 |
| SYCL | 1e-6 | 179.91 | 210.73 | 216.92 |
| SYCL | epsilon | 191.59 | 221.14 | 212.83 |

The cache-only control reproduces the slowdown. In the OpenMP epsilon
case, total iterations rise from 11,801 to 13,446 after the complete
repair, and stencil calls from 7,909 to 11,205. The larger increase in
calls reflects unequal progress among batch slots. Final CPU accuracy
passes throughout. The adjoint repair is a correctness improvement; it
does not guarantee fewer iterations for every source and tolerance.

An additional study of the repaired code on that same batch compares
inner tolerances 1e-4, 1e-5, 1e-6 and epsilon with 12 and 18 reals.
At 1e-6, changing to 18 reals reduces OpenMP time from 254.10 to 209.58 ms
and SYCL from 210.74 to 175.37 ms, about 17% on each backend. These
numbers come from the same tuning job; their small difference from the
preceding table illustrates timing variation. Looser inner targets do
not improve this case. The measured policy is to retain the 1e-6 default
and use the existing storage/tolerance controls to study sensitive
workloads. The default compression mode has not been changed from this
single case.

### Evidence and limitations

The source/configuration manifests, job statuses, scripts, timing samples,
work counters and compressed logs are in
`/private/tmp/qex-gpu-20261001/gpu-study/results/sunspot-20261001`.
Regenerate the summary with:

```
python3 /private/tmp/qex-gpu-20261001/gpu-study/summarize_study.py \
  /private/tmp/qex-gpu-20261001/gpu-study/results/sunspot-20261001
```

The laptop archive includes the original
compression control and deliberately unattainable-tolerance runs:
scalar double/mixed and a mixed four-system batch stop on stagnation
well before 100,000 iterations, report residuals near 6e-30, and reject
the requested 1e-40 benchmark accuracy.

The record preserves superseded attempts: an inherited 208-thread
runtime setting, Python 3.6 incompatibilities in the scripts, and the
overly strict 14-real reference comparison. They are excluded from
accepted timing samples. Performance evidence covers the stated M
workloads and PVC backends; HMC results establish the tested numerical
behavior without isolating old/new HMC performance. CUDA and HIP were
not exercised in this study.

## FP64 reconstruction experiment, October 2, 2026

The optional reconstruction precision was checked with tstagrecon on the
laptop using one and six host threads. The point-source outputs exactly
match the independent host reconstruction rounded to FP32, for 12 and
14 reals with either cache mode. The test includes a nontrivial U(3)
phase for 14 reals. Changing the option affects the compressed FP32
operator; the 18-real and FP64 outputs are unchanged. The bestagres
checks also pass with FP64 reconstruction in both cache modes, including
raw scalar/batched FP32 CG and complete mixed solves.

On the original 4^3 x 8, mass 0.001 fixture, the raw inner CG on the CPU
backend gives the following results at the FP32 epsilon stopping target.
The last column recomputes the residual with the CPU operator containing
all 18 reals per link rounded to FP32, applied in FP64. It is a common
reference for the three representations.

| Storage / reconstruction | Iterations | Recursive r2/b2 | CPU rounded-link r2/b2 |
|---|---:|---:|---:|
| 12 FP32 / FP32 | 1321 | 1.149e-7 | 4.557e-4 |
| 12 FP32 / FP64 | 1316 | 1.053e-7 | 4.924e-4 |
| 18 FP32 | 1322 | 1.114e-7 | 5.536e-4 |

The distinct sources in bestagres show that the iteration effect can
have either sign: at the same epsilon target, the mass 0.001 source
changes from 1934 to 2169 iterations, while mass 0.0012 changes from
2209 to 1313. FP64 reconstruction reduces arithmetic error in the
reconstructed row, but cannot recover information lost when storing the
first two rows in FP32. Fewer iterations or a smaller row error do not
by themselves establish a more accurate or faster complete solve.

On Sunspot PVC, job 12480059 in workq with allocation Catalyst built
OpenMP and SYCL on two compute nodes. The source snapshot starts at
a8969b7f with the reconstruction option added. Builds use the previous
study's VLEN 8 configurations, oneAPI 2026.1.0, assertions and work
counters. Runs use eight host threads and one MPI rank per GPU tile;
OpenMP target offload is mandatory.

The same raw probe at the epsilon target is more expensive with FP64
reconstruction on both PVC backends. All recursive residuals meet that
target. The CPU residual columns use the common rounded-link reference
described above:

| Backend | 12 / FP32 iterations | 12 / FP64 iterations | 18 iterations | CPU r2/b2, 12 / FP32 | CPU r2/b2, 12 / FP64 |
|---|---:|---:|---:|---:|---:|
| OpenMP | 1314 | 2114 | 1315 | 4.979e-4 | 1.147e-3 |
| SYCL | 1314 | 2457 | 1319 | 5.419e-4 | 1.129e-3 |

At the default inner target of 1e-6, the corresponding counts are 1310,
1311 and 1310 on each backend. The cross product's arithmetic accuracy
alone therefore does not predict the FP32 CG's convergence near its
stopping threshold.

Both PVC backends pass tstagrecon on one, two and four ranks on one node,
and on two ranks across nodes. Both also pass the bestagres checks with
12, 14 and 18 reals, both cache modes, split and unsplit hops, fixed and
atomic reductions, peer IPC and MPI fallback. The inter-node CG case and
bestagcg with `-mixed:1 -recon64:1` pass as well. These checks include the
raw FP32 scalar/batched results and independent CPU checks of the final
mixed solutions.

The complete-solve comparison uses bestagm with seed 987654321, fixed
sums, r2req = 1e-16 and one or three systems on one or two tiles. Small
cases use 4^3 x 8, mass 0.001, 12 reals, forced separate caches and even
or full sources; larger cases use 24^4, mass 0.05, 12 or 14 reals, automatic
cache selection and full sources. Light cases test inner targets 1e-6 and
epsilon; larger cases use 1e-6. Batch masses are m, 1.2m and 1.4m. Each
configuration runs four solves, drops the first, then repeats in reverse
variant order, retaining six samples. Gauge and
source upload and independent CPU residual evaluation are outside the
timed solve. Every returned solution must pass the same CPU squared
residual limit of 1.01e-16.

For one tile and three systems, with inner target 1e-6, median milliseconds
and total inner iterations in parentheses are:

| Backend | Case | FP32 reconstruction | FP64 reconstruction | 18-real control |
|---|---|---:|---:|---:|
| OpenMP | Light, even source, 12 reals | 255.56 (12702) | 218.89 (11787) | 210.39 (11786) |
| OpenMP | Light, full source, 12 reals | 267.63 (14086) | 282.89 (14241) | 304.18 (14963) |
| OpenMP | 24^4, 12 reals | 207.83 (1340) | 218.32 (1340) | 209.85 (1340) |
| OpenMP | 24^4, 14 reals | 197.18 (1340) | 212.74 (1340) | 210.30 (1340) |
| SYCL | Light, even source, 12 reals | 210.83 (12789) | 188.31 (11789) | 175.30 (11776) |
| SYCL | Light, full source, 12 reals | 234.89 (14501) | 255.72 (14725) | 247.83 (14971) |
| SYCL | 24^4, 12 reals | 198.81 (1340) | 209.18 (1340) | 198.13 (1340) |
| SYCL | 24^4, 14 reals | 190.51 (1340) | 204.45 (1339) | 197.74 (1340) |

The light even-source batch improves by 14.3% on OpenMP and 10.7% on SYCL.
Its full-source counterpart slows by 5.7% and 8.9%. Across the eight
larger-lattice configurations per backend, the median time increases by
3.7-12.6% on OpenMP and 3.2-13.0% on SYCL. Their iteration counts are
unchanged on OpenMP and differ by at most one on SYCL. The epsilon target
can increase complete-solve time too: the one-system even-source OpenMP
case changes from 176.68 to 211.45 ms. The option therefore stays off by
default.

The short control against the previous executable initially gives a 4.2%
SYCL timing difference with the option disabled, while successive samples
are still settling. Job 12480060 repeats that control with 12 solves per
invocation, discards the first four, and retains 16 samples per executable
in old/new/new/old order. Old/new medians are 199.96/200.51 ms for OpenMP
and 191.90/192.28 ms for SYCL, differences of 0.27% and 0.20%, with
overlapping ranges. Both executables take 1340 iterations and return
identical CPU residuals on each backend. This control uses the one-tile,
three-system 24^4 case with 14 reals and full sources.

All 336 commands in the main PVC study complete successfully. Each timing
configuration repeats its iteration count across the six retained samples.
The largest CPU squared residual in these benchmarks is 9.9831e-17.
Scripts, source hashes, configurations, logs and all timings are stored outside the tree
in `/private/tmp/qex-recon64-20261002`, with the remote study in
`/home/xyjin/W/qex_recon64_20261002`. Regenerate the main timing summary with:

```
python3 /private/tmp/qex-recon64-20261002/summarize.py \
  /private/tmp/qex-recon64-20261002/pvc/logs
```
