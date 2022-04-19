nimcache = getCurrentDir() / "nimcache"
buildVerbosity = 0

ccType = "gcc"

cc = "mpicc"
cflagsAlways = ""
cflagsDebug = "-g"
cflagsSpeed = "-Ofast -march=native"
ld = cc
ldflags = cflagsAlways & " -ldl"

cpp = "mpicxx"
cppflagsAlways = ""
cppflagsDebug = "-g"
cppflagsSpeed = "-Ofast -march=native"
ldpp = cpp
ldppflags = cppflagsAlways & " -ldl"

simd = ""
vlen = "8"

qmpDir = getHomeDir()/"lqcd/install/qmp"
qioDir = getHomeDir()/"lqcd/install/qio"
qudaDir = ""
cudaLibDir = ""
primmeDir = ""
chromaDir = ""
gridDir = ""

envs = @[]
