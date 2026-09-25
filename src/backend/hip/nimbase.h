/* Nim's nimbase.h for the HIP backend (clang -x hip, found first through
   -iquote): Nim inline procs compile for the device too, so that gpuFor
   kernels may call them; clang emits them for the device only when a kernel
   does.  qexFor runs the gpuFor lambdas, one index per thread. */
#include_next "nimbase.h"
#ifdef __HIP__
#include <hip/hip_runtime.h>  /* __launch_bounds__ */
#undef N_INLINE
#define N_INLINE(rettype, name) __host__ __device__ inline rettype name
/* HIP defines __noinline__ as a macro */
#undef N_NOINLINE
#define N_NOINLINE(rettype, name) rettype __attribute__((noinline)) name
template <int T, typename F> __global__ void __launch_bounds__(T) qexFor(long n, F f) {
  long i = (long)blockIdx.x * T + threadIdx.x;
  if (i < n) f(i);
}
#endif
