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
all but the B200, which ran the code before them: both CI suites,
tests/base/taccel and tgaugegpu, backend/examples/bestagres (also with
`-split:0` and `-split:1`, links of 12 reals and an 8^4 lattice, and with
`-ipc:0` where MPI reads device memory, on Sunspot and Tuolumne), berng and
staghmcgpu `-check:1`, on 1, 2 and 4 GPUs and on 12 PVC tiles.
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
does; it records |b - M x|^2/|b|^2 in the SolverParams statistics, 0 for
b = 0.  The true |q - A x_e| can't drop below about eps |A| |x_e|, which for
small m and r2req lies above the CG target (m = 0.001 and r2req = 1e-24 on
a 4^3x8 lattice): rounding x_e by d changes b - M x by A d/(4m) through x_o,
while the correction y is small and so is its rounding.  The restarts in
double of the mixed precision solver also stop once the residual no longer
drops.  In solveM of several systems, a system whose slot stops takes these
steps on its own before the slot takes the next system.  A solve out of
iterations stops there and records the residual it reached.

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
- backend/examples/bestagres: the stopping rules of solveM and solveEE
  against stag.D and stagD2ee: one system and several, double and mixed
  precision, fixed and atomic dot products, sources on all or on the even
  sites, masses of both signs, zero sources and too few iterations;
  `-ipc` and `-split` as bestagcg; exits with 1 if a check fails.
- backend/examples/berng: the generators of RNGFieldV and rng/rngGpu against
  the host RNG fields; exits with 1 if a check fails.
- tests/base/taccel: the order of queued kernels, sumFixed and gpuSum at
  the sizes where their levels change, and the neighbors GpuHaloEx brings,
  on 1, 2 and 4 ranks.
- tests/base/tgaugegpu: gauge/gaugeGpu against actionA and forceA, and
  ValueError for rect and pgm coefficients.
- backend/examples/bernt: the times of the rngGpu draw kernels, with 32-,
  64- and 128-bit generator loads.
- examples/staghmcgpu_sh: staghmc_sh on the GPU; `tests/extra/tstaghmc_sh/run`
  with `RUNJOB` set to the launcher.
- prod/lsd/eightFlavorSMGgpu: eightFlavorSMG on the GPU;
  `tests/extra/teightFlavorSMG/run` (needs python 3.8 or later).

## Measured

One GPU (one PVC tile), devel of September 30, the builds above; the B200
ran the code before the stopping rules.

| | H100 | B200 | MI300X | MI300A | PVC tile, OpenMP / SYCL |
|---|---|---|---|---|---|
| HBM peak (TB/s) | 3.35 | 8.0 | 5.3 | 5.3 | 1.64 |
| bestream triad (TB/s) | 3.10 | 5.83 | 4.09 | 3.35 | 0.98 / 0.98 |
| bestagcg 32^4, one system (TB/s) | 2.65 | 5.30 | 3.21 | 2.61-2.62 | 0.82 / 0.83 |
| staghmcgpu_sh 24^4, double / mixed (s per trajectory) | 0.887-0.889 / 0.679-0.680 | 0.605 / 0.509 | 0.954-0.961 / 0.770-0.772 | 1.064-1.069 / 0.866-0.885 | 2.03-2.04 / 1.68-1.69, 1.98-2.00 / 1.61-1.62 |
| eightFlavorSMGgpu 16^3x32 (s per trajectory) | 2.72-2.73 | 2.11 | 2.79-2.81 | 3.44-3.46 | 5.04-5.07 / 4.75-4.77 |

Against the code before the stopping rules, in the same jobs (runs in
the order old, new, new, old), bestagcg 32^4 takes 0.1-0.6% longer for one
system, whose CG checks its true residual at the end, and 0.9-1.6% for
three systems per hop, which finish one by one (on MI300X and MI300A those
runs spread 4%).  staghmcgpu_sh takes 2.5-3.8% longer per trajectory at
24^4 (3.3-6.3% at 16^4): its action solves, whose sources have odd sites,
take 13-15% more iterations to meet the request on M x = b, and every solve
checks b - M x.  eightFlavorSMGgpu, whose solves take the same iterations as
before, takes 0.4-2.5% longer.
