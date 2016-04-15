# QEX
##Quantum EXpressions lattice field theory framework

This is an early version of the code.  The internals will likely go through
several revisions, though the high level code written on top of it will
likely not need much (if any) revisions.

### Installation:

First you need Nim.  I recommend installing it from the instructions
given in "Installation from github" at the bottom of this page:
http://nim-lang.org/download.html

(optional) Copy "Makefile.template" to a separate build directory.

Rename "Makefile.template" to "Makefile" and edit the top section
to point to the QEX source, and set compiler, etc.
(the first section is for BG/Q, the second for x86)

The variables you'll likely need to change are:

| QEXDIR | root directory containing QEX code (where this README is) |
| QMPDIR, QIODIR | installation directories for respective codes |
| CC | C compiler to use |

CC_TYPE: What compiler dialect the Nim code generator should use
         common options are: gcc, clang, llvm_gcc, icl (Intel),
         ucc (generic unix cc)
CFLAGS_ALWAYS: CFLAGS that are always used
CFLAGS_DEBUG: CFLAGS used for a debug build (make debug=1 ...)
CFLAGS_SPEED: CFLAGS used for a release build (default)
VERBOSITY: Nim compiler verbosity
SIMD: SIMD extensions supported (QPX,SSE,AVX,AVX512)
VLEN: Simd vector length to use

Try compiling simple example:
make test2

then run it
