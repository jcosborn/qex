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
const
  qexDir = "."   # Path to qex source directory
  extraSrcDir = @["."]   # Extra paths to search for build targets.
  #primmeDir = "$HOME/pkg/src/primme"
  #qmpDir = "$HOME/pkg/qmp"
  #qioDir = "$HOME/pkg/qio"
  qmpDir = "$HOME/lqcd/install/qmp"
  qioDir = "$HOME/lqcd/install/qio"
  ccType = "gcc"
  cc = "mpicc"
  #cc = qexDir & "/src/backend/util/ccwrapper"
  cflagsAlways = "-Wall -std=gnu11 -march=native -ldl -fno-strict-aliasing"
  #cflagsAlways = "-Wall -std=gnu11 -march=native -ldl -Wa,-q"
  #cflagsAlways = "-x cu -std=c++11 -Xcompiler -mcpu=native"
  cflagsDebug = "-g3 -O0"
  cflagsSpeed = "-g -Ofast"
  #cflagsSpeed = "-O3 --use_fast_math -Xcompiler -Ofast" #"-g -Ofast"
  ompFlags = "-fopenmp"
  #cflagsAlways = "-qinfo=pro"
  #cflagsAlways = "-qinfo=pro -qstrict=operationprecision"
  #cflagsDebug = "-g3 -O0"
  #cflagsSpeed = "-g -O3"
  #ompFlags = "-qsmp=omp"
  ld = cc
  ldflags = "-Wall -std=gnu11 -march=native -ldl"
  #ldflags = cflagsAlways
  nimcache = "nimcache"
  verbosity = 1
  simd = "SSE AVX"
  vlen = "8"
  #simd = "QPX"
  #vlen = "4"
  extraDef = ["STATIC_UNROLL=1"]
  ########################################
  # Optional dependencies
  #primmeDir = "$HOME/pkgs/src/primme"
  #primmeDir = "$HOME/pkgs/build/primme"
  #lapackLib = "'-Wl,-framework -Wl,Accelerate'"
  #lapackLib = "'$HOME/pkg/lib/libopenblas.a -fopenmp -lm -lgfortran'"
  #lapackLib = "'"&llbgq.join" "&"'"
  #qudaDir = "$HOME/lqcdM/build/quda"
  #cudaLibDir = "/usr/local/cuda/lib64"
  #cudaArch = "sm_37"
  #cudaCCBIN = "gcc"
  #cudaNVCC = "nvcc"
# End of configurations
######################################################################
