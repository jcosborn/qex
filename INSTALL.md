## Quick guide

Clone the github repo (devel branch recommended for now).

Create a separate build directory (optional but recommended).

From the build directory run the `configure` script found with the source.

```
QMPDIR='/path/to/qmp' \
QIODIR='/path/to/qio' \
QUDADIR='/path/to/quda' \
CUDADIR='/path/to/cuda/lib' \
path/to/qex/configure
```

More options to `configure` are
[here](INSTALL.md#compiler-and-configuration-options).

If the `nim` executable isn't found in your path, it will be installed
as described [here](INSTALL.md#nim-installation).

The `configure` command will create the files `Makefile` and `qexconfig.nims`
in the current directory.
It will also create the symlinks `qex` and `qex.nimble`.

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

QEX requires the Nim compiler version 1.4 or later
(1.6.4 or later recommended).

The `configure` script will search for the Nim executable in the following
places in order (using the [findNim](build/findNim) script)
- `$PATH` (using `which nim`)
- `$HOME/bin/nim`
- `$HOME/bin/nim-[0-9]*`
- `$HOME/bin/nim-*`
- `$HOME/nim/Nim/bin/nim`
- `$HOME/nim/Nim-[0-9]*/bin/nim`
- `$HOME/nim/Nim-*/bin/nim`

If Nim is not found, it will be installed using the
[installNim](build/installNim) script.

You can install Nim yourself either by using the `installNim` script directly,
or from the official instructions given
[here](http://nim-lang.org/download.html).

`installNim` defaults to downloading and building Nim into `$HOME/nim`
with a symlink added to `$HOME/bin`.
These can be overridden by setting the environment variables
- `NIMDIR` directory to place Nim source and build (default `$HOME/nim`)
- `BINDIR` directory to place symlinks (default `$HOME/bin`)

This will be necessary to set if on a system with a shared home directory
across different host (build) architectures.

By default `installNim` installs the latest stable version of Nim.
You can also specify a version, e.g. `installNim 1.6.4` for version 1.6.4, or
`installNim devel` to install the current devel branch.

You can also switch the symlinks between different (previously installed)
versions of Nim with e.g. `installNim default 1.6.4` or `installNim default devel`.
You should be able to check which version of Nim is default with `nim -v`.


## Required dependencies

QMP and QIO are currently required and need to be installed separately.
The [bootstrap-travis](bootstrap-travis) script can also be used to install them.

The QMP and QIO directories can be set in `configure` with the environment variables
`QMPDIR` and `QIODIR`.
One can also specify them as arguments to `configure`
(see [here](INSTALL.md#compiler-and-configuration-options) ).

The other required dependencies can be installed by running
`nimble install -dy`.
[Nimble](https://github.com/nim-lang/nimble) is the Nim package manager
which is used to install Nim dependencies.
Nimble is also installed (with a symlink in `$BINDIR`) by the `installNim` script.
The list of required packages that Nimble installs
is in the [qex.nimble](qex.nimble) file.


## Optional dependencies

For Chroma set
- CHROMADIR

For Grid set
- GRIDDIR

For QUDA set
- QUDADIR
- CUDALIBDIR

For Primme set
- ???

## Compiler and configuration options

The available options and their default settings can be found in
[build/configDefault.nims](build/configDefault.nims).


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

## Build guide

Details on building applications can be found in [BUILD.md](BUILD.md).

The test suite can be built with `make tests`.

This creates a `testscript.sh` in the current directory which will
run the tests.
The comments in that file explain what environment variables can be set
to run it with a parallel job launcher (i.e. mpirun).


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
