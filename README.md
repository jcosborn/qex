## QEX: Quantum EXpressions lattice field theory framework

### Notice:

This code is still undergoing major development, but should be ready for
production use soon.

### Installation:

First you need [Nim](https://nim-lang.org).

You can install it either by using the script "installNim"
in this repo, or from the instructions given here:
http://nim-lang.org/download.html

You can also skip this now, and the configure script will install Nim
using the "installNim" script if it can't find the "nim" executable.

Create a separate build directory (optional but recommended).

From the build directory run the "configure" script found with the source.

```
QIODIR='/path/to/qio' QMPDIR='/path/to/qmp' \
QUDADIR='/path/to/quda' CUDADIR='/path/to/cuda/lib' \
./configure
```

This will create `Makefile`, `Makefile.nims`, `config.nims` in the build directory.
Check the resulting `config.nims` and edit if necessary.

The variables you may need to change are:

```
qexdir: root directory containing QEX code (where this README.md is)
qmpdir, qiodir, qudadir: installation directories for respective codes
cudadir: directory containing cuda runtime libraries
cc: C compiler to use
ccType: What compiler dialect the Nim code generator should use.
         Common options are: gcc, clang, llvm_gcc, icl (Intel),
         ucc (generic unix cc).
         The full list of known compilers is at the bottom of this page:
         https://github.com/nim-lang/Nim/wiki/Consts-defined-by-the-compiler
cflagsAlways: CFLAGS that are always used
cflagsDebug: extra CFLAGS used for a debug build (make debug=1 ...)
cflagsSpeed: extra CFLAGS used for a release build (default)
verbosity: Nim compiler verbosity
simd: SIMD extensions supported (QPX,SSE,AVX,AVX512)
vlen: Simd vector length to use
```

Try compiling a simple example:
```
make testStagProp
```
The resulting binary will be in the directory `./bin`.

Then run it
```
./bin/testStagProp
```
