# QEX configuration file

# location of temporary build files (generated C/C++ source, objects, etc.)
nimcache = getCurrentDir() / "nimcache"

# Nim compiler build verbosity (0-3, 1 is Nim default)
buildVerbosity = 1

# type of compiler, used for Nim generated compile flags
# typically "gcc" or "clang" will work fine for most modern compilers
# https://github.com/nim-lang/Nim/wiki/Consts-defined-by-the-compiler
ccType = "gcc"

# default language backend, "cc" (C) or "cpp" (C++)
# "cpp" is required for Chroma and Grid
ccDef = "cc"

# options for C compiler
# cc            C compiler executable
# cflagsAlways  C flags for all builds (debug and release)
# cflagsDebug   C flags for debug build
# cflagsSpeed   C flags for release build
# ld            C linker executable (typically same as C compiler)
# ldflags       C linker flags
cc = "mpicc"
cflagsAlways = ""
cflagsDebug = "-g"
cflagsSpeed = "-Ofast -march=native"
ld = cc
ldflags = cflagsAlways & " -ldl"

# options for C++ compiler
# cpp             C++ compiler executable
# cppflagsAlways  C++ flags for all builds (debug and release)
# cppflagsDebug   C++ flags for debug build
# cppflagsSpeed   C++ flags for release build
# ldpp            C++ linker executable (typically same as C++ compiler)
# ldppflags       C++ linker flags
cpp = "mpicxx"
cppflagsAlways = ""
cppflagsDebug = "-g"
cppflagsSpeed = "-Ofast -march=native"
ldpp = cpp
ldppflags = cppflagsAlways & " -ldl"

# SIMD intrinsics to use, comma separated list of SSE, AVX, AXV512
#   "SSE,AVX" for up to AVX2
#   "SSE,AVX,AVX512" for up to AVX512
#   "" for no explicit intrinsics
simd = ""

# default inner (SIMD) vector length (for both single and double)
#   this is independent of hardware SIMD length
#   (i.e. 16 is valid on any platform)
vlen = 8

# required libraries
qmpDir = getHomeDir()/"lqcd/install/qmp"
qioDir = getHomeDir()/"lqcd/install/qio"

# optional libraries
qudaDir = ""
cudaLibDir = ""
nvhpcDir = ""
primmeDir = ""
chromaDir = ""
gridDir = ""

# seq of extra environment variables to define during build
#   e.g. @["OMPI_CXX=foo","BAR=1","ENV_WITH_SPACE_QUOTES=SPACE \"WITH QUOTES\""]
envs = @[]

# seq of extra arguments for Nim during build
#   e.g. @["--listCmd","-d:defPrec:S","-d:nc=4"]
nimargs = @[]
