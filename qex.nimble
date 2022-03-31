import os, macros

var nim = selfExe()
var qexDir = thisDir()
if dirExists(qexDir/"qex"): qexDir = qexDir/"qex"
var nimArgs = newSeq[string]()
echo "Nim: ", nim
echo "QEX dir: ", qexDir

# workaround for limitation of include
macro incl(s: static string): untyped = quote do: include `s`

include "build/configBase.nims"
const qc = getCurrentDir() / "qexconfig.nims"
when fileExists(qc): incl qc
include "build/buildTasks.nims"

# Package

version       = "0.0.0"
author        = "James C. Osborn"
description   = "Quantum EXpressions lattice field theory framework"
license       = "MIT"
srcDir        = "qex/src"

# Dependencies

requires "nim >= 1.0.2"
requires "chebyshev >= 0.2.1"
requires "mdevolve >= 1.0.0"
if primmeDir != "":
  requires "primme >= 3.0.0"

# Tasks

task make, "compile, link, and put executables in `bin' (TODO)":
  echo "make"

task clean, "remove temporary build files (TODO)":
  runClean()

task tests, "Build tests":
  buildTests()

task show, "Show Nim compile flags":
  setNimFlags()
  echo "Nim flags:"
  echo join(nimFlags," ")

task targets, "List available targets":
  const c = paramCount()
  var i = 1
  while paramStr(i) != "targets": inc i
  var f = ""
  if i < c:
    f = paramStr(i+1)
    echo "Searching for targets matching: ", f
  runTargets(f)

task help, "print out usage information":
  echo "----------------------------------------------------------------------"
  echo "To build nim files:"
  echo "  nimble make [debug] [FlagsToNim] [Name=Definition] Target [MoreTargets]"
  echo ""
  echo "`debug' will make debug build."
  echo "`Target' can be any file name, optionally with partial directory names."
  echo "The produced executables will be under `bin/'."
  echo ""
  echo "Examples:"
  echo "  nimble make debug test0"
  echo "  nimble make example/testStagProp"
