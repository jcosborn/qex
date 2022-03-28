nimcache = getCurrentDir() / "nimcache"
buildVerbosity = 0

ccType = "gcc"

cc = "mpicc"
cflagsAlways = ""
cflagsDebug = "-g"
cflagsSpeed = "-Ofast -march=native"
ld = cc
ldflags = cflagsAlways

cpp = "mpicxx"
cppflagsAlways = ""
cppflagsDebug = "-g"
cppflagsSpeed = "-Ofast -march=native"
ldpp = cpp
ldppflags = cppflagsAlways

simd = ""
vlen = "8"

qmpDir = getHomeDir()/"lqcd/install/qmp"
qioDir = getHomeDir()/"lqcd/install/qio"
qudaDir = ""
cudaDir = "/usr/local/cuda"
chromaDir = ""
