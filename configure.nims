import osPaths
import strUtils

var params = newSeq[string](0)
template set(key,val: string) =
  params.insert val
  params.insert key
template `~`(key,val: untyped) =
  set(astToStr(key), val)
var envs = newSeq[string](0)

let nim = paramStr(0)
NIM ~ nim

let script = paramStr(1)
#let qexdir = thisDir()
let (qexdir, _, _) = splitFile(script)
echo "Using QEXDIR ", qexdir
set "QEXDIR", qexdir

let home = getHomeDir()

var qmpdir = home / "lqcd/install/qmp"
echo "Using QMPDIR ", qmpdir
set "QMPDIR", qmpdir

var qiodir = home / "lqcd/install/qio"
echo "Using QIODIR ", qiodir
set "QIODIR", qiodir

var machine = ""
CC ~ "mpicc"
CC_TYPE ~ "gcc"
CFLAGS_ALWAYS ~ "-Wall -std=gnu99 -march=native -ldl"
CFLAGS_DEBUG ~ "-g3 -O0"
CFLAGS_SPEED ~ "-g -O3"
OMPFLAGS ~ "-fopenmp"
LD ~ ( "$CC" % params )
LDFLAGS ~ ( "$CFLAGS_ALWAYS" % params )
VERBOSITY ~ "1"
SIMD ~ ""
VLEN ~ "1"

if dirExists "/bgsys":
  machine = "Blue Gene/Q"
  CC ~ qexdir / "build/mpixlc2"
  #CC ~ "mpixlc2"
  CFLAGS_ALWAYS ~ "-qinfo=pro"
  CFLAGS_DEBUG ~ "-g3 -O0"
  CFLAGS_SPEED ~ "-g -O3"
  OMPFLAGS ~ "-qsmp=omp"
  LD ~ ( "$CC" % params )
  LDFLAGS ~ ( "$CFLAGS_ALWAYS" % params )
  #VERBOSITY ~ "3"
  SIMD ~ "QPX"
  VLEN ~ "4"
  envs.add "STATIC_UNROLL=1"

let uname = staticExec "uname"

if machine=="" and uname=="Darwin":
  machine = "macOS"
  let features = staticExec "sysctl -n machdep.cpu.features"
  var simd = newSeq[string](0)
  var vlen = 1
  if features.contains("SSE4"):
    simd.add "SSE"
    vlen = 4
  if features.contains("AVX"):
    simd.add "AVX"
    vlen = 8
  if vlen>1:
    SIMD ~ join(simd,",")
    VLEN ~ $vlen

if machine=="" and fileExists "/proc/cpuinfo":
  # assume compute nodes are same as build nodes
  machine = "linux"
  SIMD ~ "SSE,AVX"
  VLEN ~ "8"

# cray/modules
# KNL
# linux (/proc/cpuinfo)
# check on linux/mac/vesta/cooley/theta
# gcc/icc opt level

var es = ""
for e in envs:
  es.add("envs.add \"" & e & "\"\n")
ENVS ~ es

proc confFile(fn: string) =
  FILE ~ thisDir() / fn
  var f = readFile(qexdir / fn & ".in")
  f = replace(f, "$", "!DOLLAR!")
  f = replace(f, "#", "!HASH!")
  f = replace(f, "@@", "$")
  f = f % params
  f = replace(f, "!HASH!", "#")
  f = replace(f, "!DOLLAR!", "$")
  #echo f
  writeFile(fn, f)

confFile("Makefile")
confFile("Makefile.nims")
confFile("config.nims")
