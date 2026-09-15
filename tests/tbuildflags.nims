import os, macros

let nim = paramStr(0)
const qexDir = thisDir().parentDir
macro incl(s: static string): untyped = quote do: include `s`
include "../build/configBase.nims"
include "../build/buildTasks.nims"

const
  gccFlag = "-fno-allow-store-data-races"
  clangFlag = "-fno-slp-vectorize"
let (tmp, code) = gorgeEx("mktemp -d")
doAssert code == 0, tmp
let dir = tmp.strip
let cmd = "exec " & findExe("cc").quoteShell & " -undef "
var num = 0

proc script(body: string): string =
  inc num
  result = dir / ("compiler " & $num)
  writeFile(result, "#!/bin/sh\n" & body & "\n")
  exec "chmod +x " & result.quoteShell

proc check(expected = "", args: seq[string] = @[]) =
  setUserNimFlags(args)
  nimFlags.setLen(0)
  setNimFlags()
  for flag in [gccFlag, clangFlag]:
    doAssert ("--passC:" & flag in nimFlags) == (flag == expected), nimCmdArgs
  if expected != "": doAssert nimFlags[^1] == "--passC:" & expected, nimCmdArgs

try:
  # Use the real preprocessor with controlled compiler identity macros.
  var gcc: array[3, string]
  cpp = dir / "unused compiler"
  for i, v in [14, 15, 16]:
    gcc[i] = script(cmd & "-D__GNUC__=" & $v & " \"$@\"")
    cc = gcc[i]
    check(if v >= 15: gccFlag else: "")
    doAssert compilerInfo(nimFlags) == ("gcc", v)
  for mac in ["__INTEL_COMPILER", "__INTEL_LLVM_COMPILER", "__NVCOMPILER", "__PGI"]:
    cc = script(cmd & "-D__GNUC__=99 -D__clang__=1 -D__clang_major__=19 -D" & mac & "=1 \"$@\"")
    check()
  var clang: array[4, string]
  for i, v in [18, 19, 20, 23]:
    clang[i] = script(cmd & "-D__GNUC__=99 -D__clang__=1 -D__clang_major__=" & $v & " \"$@\"")
    cc = clang[i]
    check(if v == 19: clangFlag else: "")
    doAssert compilerInfo(nimFlags) == ("clang", v)

  let mpi = script("exec \"$OMPI_CC\" \"$@\"")
  cflagsSpeed = "-Ofast -march=native"
  cppflagsSpeed = "-O3 -march=native"
  cflagsDebug = "-Og"
  for (affected, unaffected, flag) in [(gcc[1], gcc[0], gccFlag), (clang[1], clang[2], clangFlag)]:
    cc = unaffected
    cpp = affected
    check()
    ccDef = "cpp"
    check(flag)
    ccDef = "cc"

    let enabled = flag.replace("-fno-", "-f")
    nimargs = @["--passC:" & enabled]
    for typ in ["gcc", "clang"]:
      check(flag, @["--cc:" & typ, "--" & typ & ".exe:" & affected.quoteShell])
    cc = affected
    check("", @["--gcc.exe:" & unaffected.quoteShell])
    nimargs.setLen(0)

    let extra = " " & flag
    check(flag)
    doAssert "--gcc.options.speed:" & ("-Ofast -march=native" & extra).quoteShell in nimFlags
    doAssert "--gcc.options.debug:" & ("-Og" & extra).quoteShell in nimFlags
    check(flag, @["--gcc.options.speed:" & ("-O2 " & enabled).quoteShell, "--gcc.options.size:'-Oz'"])
    doAssert "--gcc.options.speed:" & ("-O2 " & enabled & extra).quoteShell in nimFlags
    doAssert "--gcc.options.size:" & ("-Oz" & extra).quoteShell in nimFlags
    ccDef = "cpp"
    check(flag)
    doAssert "--gcc.cpp.options.speed:" & ("-O3 -march=native" & extra).quoteShell in nimFlags
    cpp = unaffected
    check()
    ccDef = "cc"

    cc = mpi
    envs = @["OMPI_CC=" & affected]
    check(flag)
    check("", @["--putenv:" & ("OMPI_CC=" & unaffected).quoteShell])
    envs.setLen(0)

    cc = affected
    fo.debug = true
    check(flag)
    fo.debug = false

  cc = script("exit 1")
  var failed = false
  try:
    check()
  except IOError:
    failed = true
  doAssert failed
finally:
  rmDir(dir)

echo "Build flag checks passed."
