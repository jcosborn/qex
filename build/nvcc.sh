#!/bin/bash
# nvcc for the C++ compile commands Nim writes for clang CUDA (cpp in
# qexconfig.nims, see docs/gpu_backends.md): -x cuda --cuda-gpu-arch=sm_XY
# become -x cu -arch=sm_XY, -std=gnu++NN -std=c++NN, -iquote -I, and the
# gcc options go through -Xcompiler.  Files without kernels go to the host
# compiler, unless QEX_NVCC_ALL is set.
#   CUDA_HOME       the CUDA toolkit (default: from nvcc in PATH)
#   QEX_HOSTCXX     the host compiler (default g++; an MPI wrapper works)
#   QEX_NVCC_FLAGS  more nvcc options, e.g. -Xptxas -v
cuda=${CUDA_HOME:-$(dirname $(dirname $(command -v nvcc)))}
hcxx=${QEX_HOSTCXX:-g++}
for a in "$@"; do [ "$a" = -E ] && exec $hcxx "$@"; done  # the compiler probe of the build
src=${@: -1}
if [ -z "$QEX_NVCC_ALL" ] && [[ $src == *.cpp ]] && ! grep -q -E "qexFor<|__global__" "$src"; then
  args=()
  while (($#)); do
    case "$1" in
      -x) shift ;;
      --cuda-gpu-arch=*) ;;
      *) args+=("$1") ;;
    esac
    shift
  done
  exec $hcxx -I$cuda/include "${args[@]}"
fi
args=(); xc=()
while (($#)); do
  case "$1" in
    -x) shift; args+=(-x cu) ;;
    --cuda-gpu-arch=*) args+=(-arch=${1#--cuda-gpu-arch=}) ;;
    -std=gnu++*) args+=(-std=c++${1#-std=gnu++}) ;;
    -iquote) shift; args+=(-I "$1") ;;
    -f*|-march=*|-mtune=*|-pthread|-W*) xc+=("$1") ;;
    *) args+=("$1") ;;
  esac
  shift
done
x=$(IFS=,; echo "${xc[*]}")
exec $cuda/bin/nvcc -ccbin $hcxx --extended-lambda $QEX_NVCC_FLAGS ${x:+-Xcompiler "$x"} "${args[@]}"
