/* Nim's nimbase.h for the CUDA backend (clang -x cuda or nvcc through
   build/nvcc.sh, found first through -iquote).  QEX_HD makes the procs that
   kernels call, gpuInline and gpuCall of base/basicOps, host device.  Other
   procs stay host only: nvcc compiles every host device function for the
   device and rejects host globals there.  qexFor runs the gpuFor lambdas,
   one index per thread. */
#include_next "nimbase.h"
#ifdef __CUDACC__
#define QEX_HD __host__ __device__
/* CUDA defines __noinline__ as a macro */
#undef N_NOINLINE
#define N_NOINLINE(rettype, name) rettype __attribute__((noinline)) name
template <int T, typename F> __global__ void __launch_bounds__(T) qexFor(long n, F f) {
  long i = (long)blockIdx.x * T + threadIdx.x;
#if defined(__clang__) && !defined(__NVCC__)  /* clang as the CUDA compiler: nvcc defines __clang__ for a clang host compiler */
  if (i < n) [[clang::always_inline]] f(i);  // LLVM leaves large lambdas uninlined
#else
  if (i < n) f(i);
#endif
}
#else
#define QEX_HD  /* host compiles of build/nvcc.sh */
#endif
