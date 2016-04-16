# QEX: Quantum EXpressions lattice field theory framework

This is an early version of the code.  The internals will likely go through
several revisions, though the high level code written on top of it will
likely not need much (if any) revisions.

### Installation:

First you need Nim.  I recommend installing it from the instructions
given in "Installation from github" at the bottom of this page:
http://nim-lang.org/download.html

(optional) Copy "Makefile.template" to a separate build directory.

Copy "Makefile.template" to "Makefile" and edit the top section
to point to the QEX source, set compiler, etc.
(the first section is for BG/Q, the second for x86)

The variables you'll likely need to change are:

```
QEXDIR: root directory containing QEX code (where this README.md is)
QMPDIR, QIODIR: installation directories for respective codes
CC: C compiler to use
CC_TYPE: What compiler dialect the Nim code generator should use
         common options are: gcc, clang, llvm_gcc, icl (Intel),
         ucc (generic unix cc)
         The full list of known compilers is at the bottom of this page
         https://github.com/nim-lang/Nim/wiki/Consts-defined-by-the-compiler
CFLAGS_ALWAYS: CFLAGS that are always used
CFLAGS_DEBUG: extra CFLAGS used for a debug build (make debug=1 ...)
CFLAGS_SPEED: extra CFLAGS used for a release build (default)
VERBOSITY: Nim compiler verbosity
SIMD: SIMD extensions supported (QPX,SSE,AVX,AVX512)
VLEN: Simd vector length to use
```

Try compiling a simple example:
```
make testStagProp
```

Then run it
```
./testStagProp
```
