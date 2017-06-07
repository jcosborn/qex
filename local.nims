import strutils
const xlbgq = "/soft/compilers/ibmcmp-may2016/"
const llbgq = [
  "/soft/libraries/alcf/current/xl/LAPACK/lib/liblapack.a",
  "/soft/libraries/alcf/current/xl/BLAS/lib/libblas.a",
  "-L"&xlbgq&"xlf/bg/14.1/bglib64",
  "-L"&xlbgq&"xlsmp/bg/3.1/bglib64",
  "-lxlf90_r -lxlfmath -lxlsmp -lxl" ]

######################################################################
# Configurations
# Comment starts with `#'
# Extra quotes may be required for nimble doesn't quote properly.
const
  qexDir = "."   # Path to qex source directory
  extraSrcDir = @["."]   # Extra paths to search for build targets.
  #primmeDir = "$HOME/pkg/src/primme"
  #qmpDir = "$HOME/pkg/qmp"
  #qioDir = "$HOME/pkg/qio"
  primmeDir = "$HOME/pkgs/src/primme"
  qmpDir = "$HOME/lqcd/install/qmp"
  qioDir = "$HOME/lqcd/install/qio"
  #primmeDir = "$HOME/pkgs/build/primme"
  #qmpDir = "/home/osborn/lqcd/install/qmp"
  #qioDir = "/home/osborn/lqcd/install/qio"
  qudaDir = "$HOME/lqcdM/build/quda"
  cudaLibDir = "/usr/local/cuda/lib64"
  #lapackLib = "'-Wl,-framework -Wl,Accelerate'"
  lapackLib = "'$HOME/pkg/lib/libopenblas.a -fopenmp -lm -lgfortran'"
  #lapackLib = "'"&llbgq.join" "&"'"
  ccType = "gcc"
  cc = "mpicc"
  #cc = "/home/xyjin/pkgs/src/qex/build/mpixlc2"
  cflagsAlways = "'-Wall -std=gnu99 -march=native -ldl -Wa,-q'"
  #cflagsAlways = "'-Wall -std=gnu99 -march=native -ldl'"
  cflagsDebug = "'-g3 -O0'"
  cflagsSpeed = "'-g -Ofast'"
  ompFlags = "-fopenmp"
  #cflagsAlways = "'-qinfo=pro'"
  #cflagsAlways = "'-qinfo=pro -qstrict=operationprecision'"
  #cflagsDebug = "'-g3 -O0'"
  #cflagsSpeed = "'-g -O3'"
  #ompFlags = "-qsmp=omp"
  ld = cc
  ldflags = "'-Wall -std=gnu99 -march=native -ldl'"
  #ldflags = cflagsAlways
  nimcache = "nimcache"
  verbosity = 1
  simd = "SSE AVX"
  vlen = "8"
  #simd = "QPX"
  #vlen = "4"
  extraDef = ["STATIC_UNROLL=1"]
# End of configurations
######################################################################
