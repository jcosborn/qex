# Configurations
const
  # Extra quotes may be required for nimble doesn't quote properly.
  primmeDir {.strdefine.} = "$HOME/pkg/src/primme"
  #primmeDir {.strdefine.} = "$HOME/pkgs/src/primme"
  qexDir {.strdefine.} = "."
  qmpDir {.strdefine.} = "$HOME/pkg/qmp"
  qioDir {.strdefine.} = "$HOME/pkg/qio"
  #qmpDir {.strdefine.} = "$HOME/lqcd/install/qmp"
  #qioDir {.strdefine.} = "$HOME/lqcd/install/qio"
  qudaDir {.strdefine.} = "$HOME/lqcdM/build/quda"
  cudaLibDir {.strdefine.} = "/usr/local/cuda/lib64"
  lapackLib {.strdefine.} = "'-Wl,-framework -Wl,Accelerate'"
  #lapackLib {.strdefine.} = "'$HOME/pkg/lib/libopenblas.a -fopenmp -lm -lgfortran'"
  ccType = "gcc"
  cc = "mpicc"
  cflagsAlways = "'-Wall -std=gnu99 -march=native -ldl -Wa,-q'"
  #cflagsAlways = "'-Wall -std=gnu99 -march=native -ldl'"
  cflagsDebug = "'-g3 -O0'"
  cflagsSpeed = "'-g -Ofast'"
  ompFlags = "-fopenmp"
  ld = "mpicc"
  ldflags = "'-Wall -std=gnu99 -march=native -ldl'"
  nimcache = "nimcache"
  verbosity = 1
  simd = "SSE AVX"
  vlen = "8"

# Package

version       = "0.0.0"
author        = "James C. Osborn"
description   = "Quantum EXpressions lattice field theory framework"
license       = "MIT"
srcDir        = "src"

# Dependencies

requires "nim >= 0.16.0"
when declared(primmeDir):
  requires "primme >= 0.0.0"

const paths = @["src", "tests"]

import ospaths, sequtils, strutils

type NamePath = tuple[n,p:string]
proc targets(p:string):seq[NamePath] =
  p.listFiles.filterIt(it.splitFile.ext==".nim").mapIt((it.splitFile.name,it))
proc recTargets(ps:seq[string]):seq[NamePath]
proc recTargets(p:string):seq[NamePath] = p.targets & p.listDirs.recTargets
proc recTargets(ps:seq[string]):seq[NamePath] =
  result = @[]
  for p in ps: result &= p.recTargets
let cmd = getCommand()
echo cmd

template set(k:string, v:untyped) =
  when compiles(type(v)):
    when type(v) is string:
      let s = v
    elif type(v) is int:
      let s = $v
    else:
      let s = astToStr(v)
  else:
    let s = astToStr(v)
  switch(k,s)
template `~`(k,v:untyped) = astToStr(k).set v
template `!`(k,v:untyped) = (ccType&"."&astToStr(k)).set v
template def(v:untyped) = define ~ (astToStr(v)&"="&v)
template define(v:untyped) = define ~ v

path ~ src
cc ~ ccType
exe ! cc
linkerexe ! ld
options.always ! cflagsAlways
options.debug ! cflagsDebug
options.speed ! cflagsSpeed
options.linker ! ldflags
putenv ~ ("OMPFLAG=" & ompFlags)
putenv ~ ("QMPDIR=" & qmpDir)
putenv ~ ("QIODIR=" & qioDir)
putenv ~ ("QUDADIR=" & qudaDir)
putenv ~ ("CUDADIR=" & cudaLibDir)
putenv ~ ("VLEN=" & $vlen)
threads ~ on
tlsEmulation ~ off
verbosity ~ verbosity
nimcache ~ nimcache
warning[SmallLshouldNotBeUsed] ~ off

when declared(primmeDir):
  def primmeDir
  def lapackLib

task make, "compile and link":
  const c = paramCount()
  when c > 2:
    for i in 2..c:
      exec("nimble make "&paramStr(i))
  else:
    for s in simd.split(" "): define s
    when declared(debug):
      echo "debug build"
    else:
      define ~ release
      obj_checks ~ off
      field_checks ~ off
      range_checks ~ off
      bound_checks ~ off
      overflow_checks ~ off
      assertions ~ off
      stacktrace ~ off
      linetrace ~ off
      debugger ~ off
      line_dir ~ off
      dead_code_elim ~ on
      opt ~ speed
    setCommand("c",paramStr(2))
