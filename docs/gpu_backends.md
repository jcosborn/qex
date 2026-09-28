# GPU backends: configuring the compilers

The GPU code (backend/accel gpuFor kernels, physics/stagGpu, gauge/hypGpu,
gauge/gaugeGpu, comms/halogpu, hmc/hmcActionGpu) builds with four backends,
selected by `-d:Backend=...` in `nimargs`:

| Backend | Compiler | Kernels | Halo stores into peers on the node |
|---|---|---|---|
| OpenMP | icx (oneAPI), OpenMP offload | `target teams distribute parallel for` | Level Zero IPC (comms/zeipc) |
| SYCL | icpx (oneAPI) | `parallel_for` lambdas | Level Zero IPC (comms/zeipc) |
| CUDA | clang -x cuda | device lambdas of the qexFor<<<>>> template | CUDA IPC (comms/cudaipc) |
| HIP | amdclang -x hip | device lambdas of the qexFor<<<>>> template | HIP IPC (comms/hipipc) |

Kernels (gpuFor bodies) call templates and the procs marked gpuInline
(forced inline, as alwaysInline) or gpuCall (inlined as Nim's inline procs)
in base/basicOps.  CUDA and HIP make these host device through QEX_HD
(backend/cuda/nimbase.h, backend/hip/nimbase.h) and reject a kernel call of
any other proc.  OpenMP and SYCL compile the inline procs a kernel calls
without annotations, as Nim writes them into each C++ file that uses them.

The settings below go into `qexconfig.nims` of the build directory (the
`configure` options of the same names set them); only the lines that
differ from the defaults are shown.  Every build directory needs its own
`nimcache`: builds sharing one corrupt each other.

## OpenMP offload, Intel PVC (Aurora, Sunspot)

Nim generates C, compiled by the MPI wrapper around icx with the offload
flags in OMPFLAG; the device code is compiled ahead of time at link time.

    ccType = "clang"
    ccDef = "cc"
    cc = "mpicc"
    cflagsSpeed = "-O3 -march=native"
    ldflags = "-g -O3 -Xs '-device pvc' -ftarget-register-alloc-mode=pvc:large -ldl"
    cpp = "mpicxx"
    ldppflags = "-g -O3 -Xs '-device pvc' -ftarget-register-alloc-mode=pvc:large -ldl"
    simd = "auto"
    vlen = 8
    envs = @["OMPFLAG=-fiopenmp -fopenmp-targets=spir64_gen"]
    nimargs = @["-d:Backend=OpenMP"]

The link needs `-O3`: icx treats `-g` alone as `-O0` for the device code.

## SYCL, Intel PVC

Nim generates C++, compiled by the MPI wrapper around icpx.

    ccType = "clang"
    ccDef = "cpp"
    cc = "mpicc"
    cflagsSpeed = "-O3 -march=native"
    cpp = "mpic++"
    cppflagsAlways = "-g -fsycl -fsycl-targets=spir64_gen"
    cppflagsSpeed = "-O3 -march=native"
    ldppflags = "-g -O3 -fsycl -fsycl-targets=spir64_gen -Xsycl-target-backend '-device pvc' -ftarget-register-alloc-mode=pvc:large -ldl"
    simd = "auto"
    vlen = 8
    envs = @["OMPFLAG=-fiopenmp"]
    nimargs = @["-d:Backend=SYCL"]

For both Intel backends `-ftarget-register-alloc-mode=pvc:auto` lets the
compiler choose 128 or 256 registers per kernel: mixed precision
staghmcgpu_sh runs 5-7% faster, double precision within 2%.

Runs: one rank per tile (gpu_tile_compact.sh), 8 cores per rank,
`MPIR_CVAR_CH4_IPC_GPU_P2P_THRESHOLD=0`, `MPIR_CVAR_CH4_OFI_ENABLE_HMEM=1`
off the node, `OMP_STACKSIZE=256M` for large local volumes.

## CUDA, NVIDIA H100

clang compiles the Nim generated C++ as CUDA.  backend/cuda/nimbase.h,
found first through `-iquote`, defines QEX_HD and the qexFor kernel
template.

    ccType = "clang"
    ccDef = "cpp"
    cc = "mpicc"
    cflagsSpeed = "-O3 -march=native"
    cpp = "mpicxx"
    cppflagsAlways = "-g -x cuda --cuda-gpu-arch=sm_90 -Xarch_device -mllvm=-disable-machine-sink"
    cppflagsSpeed = "-O3 -march=native"
    ldppflags = "-g -no-pie -L$CUDA/lib64 -Wl,-rpath,$CUDA/lib64 -lcudart -ldl"
    simd = "SSE,AVX,AVX512"
    vlen = 8
    envs = @["OMPFLAG=-fopenmp"]
    nimargs = @["-d:Backend=CUDA"]

with the MPI wrappers pointed at clang (`OMPI_CC=clang OMPI_CXX=clang++`
for OpenMPI) and a CUDA toolkit clang supports (JLSE: llvm 22.1.8 with
CUDA 12.9.1).  The device flag `-disable-machine-sink` and the
`[[clang::always_inline]]` call of qexFor belong together: LLVM's machine
sinking moves the arithmetic of a hop past the loads of the later
directions, and LLVM leaves large lambdas uninlined, so the batched hop of
3 systems spilled (255 registers, up to 4 KB of stack; either change alone
does not help).  With both, 164 registers and no stack, as with nvcc, and
the default nBatch 3 fits.  To link on a node without the NVIDIA driver,
add `-L$CUDA/lib64/stubs` (the stub libcuda; the driver's at run time).

One H100 (JLSE, 2026-09-26; 3.35 TB/s HBM, stream triad 3.09 TB/s):
bestagcg 32^4 2.66 TB/s for one system (80% of the peak), 2.45 TB/s for 3
systems sharing a hop, 0.49 ns per site, system and iteration against 0.84
alone.  Per trajectory (2026-09-27, the second and third): staghmcgpu_sh
24^4, 3 terms, 0.86 s in double, 0.66 mixed; eightFlavorSMGgpu 16^3x32
2.73 s, 24^3x48 on 4 H100s 4.1 s.

## CUDA with nvcc (branch gpu-nvcc)

build/nvcc.sh turns the compile commands Nim writes for clang into nvcc
commands (`-x cu -arch=sm_XY --extended-lambda`, the gcc options through
`-Xcompiler`) for the files with kernels and gives the other files to the
host compiler.  nvcc compiles every host device function for the device and
rejects host globals there, which is why only the procs that kernels call
are host device.

    ccType = "clang"
    ccDef = "cpp"
    cc = "mpicc"
    cflagsSpeed = "-O3 -march=native"
    cpp = "<qex>/build/nvcc.sh"
    cppflagsAlways = "-g -x cuda --cuda-gpu-arch=sm_90"
    cppflagsSpeed = "-O3 -march=native"
    ldpp = "mpicxx"
    ldppflags = "-g -no-pie -L$CUDA/lib64 -Wl,-rpath,$CUDA/lib64 -lcudart -ldl"
    simd = "SSE,AVX,AVX512"
    vlen = 8
    envs = @["OMPFLAG=-fopenmp"]
    nimargs = @["-d:Backend=CUDA"]

with `CUDA_HOME` the toolkit (JLSE: CUDA 13.3.1), `QEX_HOSTCXX` the host
compiler (g++ 14), `CPATH` the MPI headers (for the host compiler and
nvcc's preprocessing) and the MPI wrappers on the host compiler
(`OMPI_CC=gcc-14 OMPI_CXX=g++-14`).  On one H100 the nvcc build runs within
1-6% of the clang build (JLSE, 2026-09-26): the same registers in the
batched hop, 2.72 TB/s for one system of bestagcg 32^4, 2.41 TB/s for 3.

## HIP, AMD MI300A

amdclang compiles the Nim generated C++ as HIP.  backend/hip/nimbase.h
includes hip_runtime.h (for `__launch_bounds__`) and defines QEX_HD and the
qexFor kernel template.  All
modules must be C++ (`ccDef = "cpp"`): in C, `__host__ __device__` does not
compile.

    ccType = "gcc"
    ccDef = "cpp"
    cc = "mpicc"
    cflagsSpeed = "-O3 -march=native"
    cpp = "mpicxx"
    cppflagsSpeed = "-O3 -march=native -x hip --offload-arch=gfx942"
    ldppflags = cppflagsAlways & " -ldl --hip-link --offload-arch=gfx942"
    simd = "SSE,AVX,AVX512"
    vlen = 8
    envs = @["STATIC_UNROLL=1"]
    nimargs = @["-d:Backend=HIP"]

with the rocmcc and cray-mpich modules (Tuolumne: rocm/10.0).  Runs:
`flux run -N1 -n4 -c24` with `MPICH_GPU_SUPPORT_ENABLED=1`; from an ssh
shell on a node of the allocation, set `FLUX_URI` to the job's instance
(`flux uri JOBID`), else `flux run` submits a new job.  The Nim phase of the
large programs takes about 13 minutes and 30 GB there.

## Programs and tests

- backend/examples/bestream: STREAM copy, scale, add and triad of gpuFor
  kernels, the bandwidth ceiling the kernels can reach.
- backend/examples/bestagcg: the GPU CG against the CPU solver, GF/s and
  GB/s by a byte count, `-nb:k` k systems at once.
- examples/staghmcgpu_sh: staghmc_sh on the GPU; `tests/extra/tstaghmc_sh/run`
  with `RUNJOB` set to the launcher.
- prod/lsd/eightFlavorSMGgpu: eightFlavorSMG on the GPU;
  `tests/extra/teightFlavorSMG/run` (needs python 3.8 or later).
