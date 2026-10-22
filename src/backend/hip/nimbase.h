/* Nim's nimbase.h for the HIP backend (clang -x hip, found first through
   -iquote).  QEX_HD makes the procs that kernels call, gpuInline and gpuCall
   of base/basicOps, host device.  qexFor runs the gpuFor lambdas, one index
   per thread. */
#include_next "nimbase.h"
#ifdef __HIP__
#include <hip/hip_runtime.h>  /* __launch_bounds__ */
#define QEX_HD __host__ __device__
/* HIP defines __noinline__ as a macro */
#undef N_NOINLINE
#define N_NOINLINE(rettype, name) rettype __attribute__((noinline)) name
template <int T, typename F> __global__ void __launch_bounds__(T) qexFor(long n, F f) {
  long i = (long)blockIdx.x * T + threadIdx.x;
  if (i < n) f(i);
}
#else
#define QEX_HD
#endif
