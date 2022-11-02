## Quick installation guide

Clone the github repo (devel branch recommended for now).

Create a separate build directory (optional but recommended).

From the build directory run the `configure` script found with the source.

```
NIM=/path/to/nim \
path/to/qex/configure \
  qmpdir:/path/to/qmp \
  qiodir:/path/to/qio
```

More options for `configure` are
[below](INSTALL.md#compiler-and-configuration-options).

If the `nim` executable isn't found in your path,
or specified in the `NIM` environment variable,
it will be installed as described
[below](INSTALL.md#nim-installation).

The `configure` command will create the files `Makefile` and `qexconfig.nims`
in the current directory.
It will also create the symlinks `qex` and `qex.nimble`.

Check the Nim executable path in `Makefile` and the
settings in `qexconfig.nims`, and edit if necessary.

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
[on the Nim website](http://nim-lang.org/download.html).

`installNim` defaults to downloading and building Nim into `$HOME/nim`
with a symlink added to `$HOME/bin`.
These can be overridden by setting the environment variables
- `NIMDIR` directory to place Nim source and build (default `$HOME/nim`)
- `BINDIR` directory to place symlinks (default `$HOME/bin`)

These will be necessary to set when on a system with a shared home directory
across different host (build) architectures to keep separate copies of Nim.

By default `installNim` installs the latest stable version of Nim.
You can also specify a version, e.g. `installNim 1.6.4` for version 1.6.4, or
`installNim devel` to install the current devel branch.

You can also switch the symlinks between different (previously installed)
versions of Nim with e.g. `installNim default 1.6.4` or `installNim default devel`.
You should be able to check which version of Nim is default with `nim -v`.

`installNim` with no options will also give a help message.


## Required dependencies

QMP and QIO are currently required and need to be installed separately.
The [bootstrap-travis](bootstrap-travis) script can also be used to install them.

The QMP and QIO directories can be set in `configure` with the environment variables
`QMPDIR` and `QIODIR`.
One can also specify them as arguments to `configure`
(see [below](INSTALL.md#compiler-and-configuration-options) ).

The other required dependencies can be installed by running
`nimble install -dy`.
[Nimble](https://github.com/nim-lang/nimble) is the Nim package manager
which is used to install Nim dependencies.
Nimble is also installed (with a symlink in `$BINDIR`) by the `installNim` script.
The list of required packages that Nimble installs
is in the [qex.nimble](qex.nimble#L25) file.

## Optional dependencies

For Chroma set
- chromaDir

For Grid set
- gridDir

For QUDA set
- qudaDir
- cudaLibDir

For Primme set
- primmeDir


## Compiler and configuration options

Compiler, linker and other configuration options can be passed
as arguments to `configure` like:
```
path/to/qex/configure \
  qmpdir:/path/to/qmp \
  qiodir:/path/to/qio \
  cc:"mpicc" \
  cflagsspeed:"-Ofast -march=skylake-avx512 -ffast-math" \
  cpp:"mpicxx" \
  cppflagsspeed:"-Ofast -march=skylake-avx512 -ffast-math"
```

Note that the option names are case insensitive
(qmpdir and qmpDir both work fine).

The options in the
[default QEX configuration file](build/configDefault.nims)
can be passed in this way, which will modify their values
in the generated `qexconfig.nims` file in the build directory.

Environment variables to be set during compile time
(set in the `envs` variable in the
[config file](build/configDefault.nims#L70) ),
can be passed one at a time using `env:FOO=BAR`.
See the [examples](#configuration-examples) below.
The full set can also be passed as a Nim seqeunce of strings using the
syntax `envs:'@["FOO=BAR","FOO2=BAR2","FOO3_WITH_SPACE_QUOTES=BAR3 \"IN QUOTES\""]'`.
The third one defines an environment variable with name `FOO3_WITH_SPACE_QUOTES` and
its content `BAR3 "IN QUOTES"`.
The single quotes above prevent shell from interpreting the Nim expression,
which uses double quotes to denote Nim strings.
And the backslashes preserve the double quotes in Nim strings.
For the single arguments (`env:FOO=BAR`) the double quotes will be added by Nim
so are unnecessary.

Details on the Nim compiler options can be found
[here](https://nim-lang.org/docs/nimc.html).


## Configuration files

The available options and their default settings can be found in
[build/configDefault.nims](build/configDefault.nims).

This default configuration file will be copied to the build directory
and renamed `qexconfig.nims`,
with the appropriate substitutions specified on the command line.


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

`<configure>` is the configure script, including path, in the QEX source directory.

### AXV2 using mpicc/mpicxx set to use gcc

```
<configure> \
  qmpdir:"$HOME/lqcd/install/qmp" \
  qiodir:"$HOME/lqcd/install/qio" \
  cc:"mpicc" \
  cflagsspeed:"-Ofast -march=native -ffast-math" \
  cpp:"mpicxx" \
  cppflagsspeed:"-Ofast -march=native -ffast-math" \
  simd:"SSE,AVX" vlen:8
```

### AXV512 using OpenMPI mpicc/mpicxx and specifying clang/clang++

```
<configure> \
  qmpdir:"$HOME/lqcd/install/qmp" \
  qiodir:"$HOME/lqcd/install/qio" \
  cctype:"clang" \
  cc:"mpicc" \
  env:"OMPI_CC=clang" \
  cflagsspeed:"-Ofast -march=native -ffast-math" \
  cpp:"mpicxx" \
  env:"OMPI_CXX=clang++" \
  cppflagsspeed:"-Ofast -march=native -ffast-math" \
  simd:"SSE,AVX,AVX512" vlen:16
```
