import os, macros

var nim = selfExe()
var qexDir = thisDir()
if dirExists(qexDir/"qex"): qexDir = qexDir/"qex"
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

requires "nim >= 1.4.0"
requires "chebyshev >= 0.2.1"
requires "mdevolve >= 1.0.0"
if primmeDir != "":
  requires "primme >= 3.0.0"

requires "https://github.com/usqcd-software/qiolite"

# Helpers

proc getExtraArgs(task: string): seq[string] =
  result.newSeq(0)
  const c = paramCount()
  var i = 1
  while paramStr(i) != task: inc i
  inc i
  while i <= c:
    result.add paramStr(i)
    inc i

var nimuserargs = newSeq[string](0)
proc getNimUserArgs(args: seq[string]): seq[string] =
  result.newSeq(0)
  for a in args:
    if a[0]=='-':
      nimuserargs.add a
    elif a[0]==':':
      let t = a[1..^1]
      if t[0]=='-':
        nimuserargs.add t
      else:
        nimuserargs.add "-d:" & t
    else:
      if '=' in a:
        nimuserargs.add "-d:" & a
      else:
        result.add a

# Tasks

task usage, "\n" & """
----------------------------------------------------------------------
QEX Nimble script usage:
  nimble <command> [build option | Nim option]... [path]...
----------------------------------------------------------------------
Commands:""":
  discard

let tClean = getTask "clean"
task clean, tClean.desc:
  runTask tClean

let tShow = getTask "show"
task show, tShow.desc:
  let ex = getExtraArgs("show")
  let optargs = getNimUserArgs(ex)
  setUserNimFlags(nimuserargs)
  let cmdargs = parseOpts(optargs)
  runTask tShow

let tTargets = getTask "targets"
task targets, tTargets.desc:
  let ex = getExtraArgs("targets")
  runTask tTargets, ex

let tTests = getTask "tests"
task tests, tTests.desc:
  let ex = getExtraArgs("tests")
  let optargs = getNimUserArgs(ex)
  setUserNimFlags(nimuserargs)
  runTask tTests

let tMake = getTask "make"
task make, tMake.desc:
  let ex = getExtraArgs("make")
  let optargs = getNimUserArgs(ex)
  setUserNimFlags(nimuserargs)
  let cmdargs = parseOpts(optargs)
  runMake(cmdargs)

let helpText =
  @[sepHelp,buildOptionsHelp,sepHelp,nimOptionsHelp,sepHelp,pathHelp,sepHelp].join("\n")

let exampleHelp = """
Examples:
  nimble make debug puregaugehmc
  nimble make bench/benchStagProp
""" & sepHelp

let tHelp = getTask "help"
task help, tHelp.desc & "\n" & helpText & "\n" & exampleHelp:
  discard
