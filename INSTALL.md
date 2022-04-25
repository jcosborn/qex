## Quick guide

Clone github repo (devel branch recommended for now).

Create a separate build directory (optional but recommended).

From the build directory run the "configure" script found with the source.

```
QIODIR='/path/to/qio' QMPDIR='/path/to/qmp' \
QUDADIR='/path/to/quda' CUDADIR='/path/to/cuda/lib' \
path/to/qex/configure
```

This will create `Makefile` and `qexconfig.nims` in the build directory.
It will also create symlinks `qex` and `qex.nimble`.
Check the resulting `qexconfig.nims` and edit if necessary.

Try compiling a simple example:
```
make testStagProp
```
The resulting binary will be in the directory `./bin`.

Then run it
```
bin/testStagProp
```


## Nim installation

First you need [Nim](https://nim-lang.org).

You can install it either by using the script "installNim"
in this repo, or from the instructions given here:
http://nim-lang.org/download.html

You can also skip this now, and the configure script will install Nim
using the "installNim" script if it can't find the "nim" executable.


## Required dependencies


## Optional dependencies

Chroma
Grid
Primme
QUDA

## Compiler and configuration options


## Configuration files

The variables you may need to change are:

```
qexdir: root directory containing QEX code (where this README.md is)
qmpdir, qiodir, qudadir: installation directories for respective codes
cudadir: directory containing cuda runtime libraries
cc: C compiler to use
ccType: What compiler dialect the Nim code generator should use.
         Common options are: gcc, clang, icl (Intel),
         ucc (generic unix cc).
         The full list of known compilers is at the bottom of this page:
         https://github.com/nim-lang/Nim/wiki/Consts-defined-by-the-compiler
cflagsAlways: CFLAGS that are always used
cflagsDebug: extra CFLAGS used for a debug build (make debug ...)
cflagsSpeed: extra CFLAGS used for a release build (default)
verbosity: Nim compiler verbosity
simd: SIMD extensions supported (SSE,AVX,AVX512)
vlen: Default SIMD vector length to use
```

## Building using ``make``

## Building using ``nimble``

You can build under the source directory.  Modify the file
`local.nims` to suit your system.

If you wish to build in a separate directory, copy the file
`local.nims` to your build directory, and modify it there.  Copy
`qex.nimble` to your build directory.

Run `nimble tasks` for available tasks, and `nimble help` for
help building your executables.

```
To build nim files:
  nimble make [debug] [FlagsToNim] [Name=Definition] Target [MoreTargets]

`debug' will make debug build.
`Target' can be any file name, optionally with partial directory names.
The produced executables will be under `bin/'.

Examples:
  nimble make debug test0
  nimble make example/testStagProp
```

## Configuration examples

The configure command can be run from a build directory that is separate
from the source directory, or from the source directory.

``<configure>`` is the configure script, including path, in the QEX source directory.

### AXV2 using mpicc/mpicxx set to use gcc

```
export QMPDIR="$HOME/lqcd/install/qmp"
export QIODIR="$HOME/lqcd/install/qio"
<configure> \
  cc:"mpicc" \
  cflagsspeed:"-Ofast -march=native -ffast-math" \
  cpp:"mpicxx" \
  cppflagsspeed:"-Ofast -march=native -ffast-math" \
  simd:"SSE,AVX" vlen:"8"
```

### AXV2 using OpenMPI mpicc/mpicxx and specifying clang/clang++

```
export QMPDIR="$HOME/lqcd/install/qmp"
export QIODIR="$HOME/lqcd/install/qio"
<configure> \
  cctype:"clang" \
  cc:"mpicc" \
  env:"OMPI_CC=clang" \
  cflagsspeed:"-Ofast -march=native -ffast-math" \
  cpp:"mpicxx" \
  env:"OMPI_CXX=clang++" \
  cppflagsspeed:"-Ofast -march=native -ffast-math" \
  simd:"SSE,AVX" vlen:"8"
```
