## QEX: Quantum EXpressions lattice field theory framework

### Notice:

This is an early version of the code, and is not fully functional.  It
is not intended for general use yet.  The internals will likely go
through several revisions, though high level code written on top
of it will likely not need much revision.  The high level framework is
incomplete, so some operations still require using lower level
constructs which may change.

### Installation:

First you need `Nim<https://nim-lang.org>`_.

You can install it either by using the script "installNim"
in this repo, or from the instructions given here:
http://nim-lang.org/download.html

Create a separate build directory (optional but recommended).

From the build directory run the "configure" script found with the source.
This will create a "Makefile" in the build directory.
Check the resulting Makefile and edit if necessary.

The variables you may need to change are:

```
QEXDIR: root directory containing QEX code (where this README.md is)
QMPDIR, QIODIR: installation directories for respective codes
CC: C compiler to use
CC_TYPE: What compiler dialect the Nim code generator should use.
         Common options are: gcc, clang, llvm_gcc, icl (Intel),
         ucc (generic unix cc).
         The full list of known compilers is at the bottom of this page:
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
