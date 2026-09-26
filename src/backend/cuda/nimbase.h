/* Nim's nimbase.h for the CUDA backend (clang -x cuda or nvcc through
   build/nvcc.sh, found first through -iquote).  With clang, Nim inline procs
   compile for the device too, so that gpuFor kernels may call them; clang
   emits them for the device only when a kernel does.  nvcc and nvc++
   compile every host device function for the device and reject host
   globals there, so for them only the gpuInline procs of base/basicOps
   (QEX_HD) are host device.  qexFor runs the gpuFor lambdas, one index per
   thread. */
#include_next "nimbase.h"
#ifdef __CUDACC__
#define QEX_HD __host__ __device__
#if defined(__clang__) && !defined(__NVCC__)  /* nvcc defines __clang__ for a clang host compiler */
#undef N_INLINE
#define N_INLINE(rettype, name) __host__ __device__ inline rettype name
#endif
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
