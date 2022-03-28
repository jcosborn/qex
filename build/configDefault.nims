nimcache = getCurrentDir() / "nimcache"
buildVerbosity = 0

ccType = "gcc"

cc = "mpicc"
cflagsAlways = "-ldl"
cflagsDebug = "-g"
cflagsSpeed = "-Ofast -march=native"
ld = cc
ldflags = cflagsAlways

cpp = "mpicxx"
cppflagsAlways = "-ldl"
cppflagsDebug = "-g"
cppflagsSpeed = "-Ofast -march=native"
ldpp = cpp
ldppflags = cppflagsAlways

simd = ""
vlen = "8"

qmpDir = getHomeDir()/"lqcd/install/qmp"
qioDir = getHomeDir()/"lqcd/install/qio"
qudaDir = ""
cudaLibDir = ""
primmeDir = ""
chromaDir = ""
gridDir = ""
