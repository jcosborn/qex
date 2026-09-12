import os, macros

let nim = paramStr(0)
const qexDir = thisDir().parentDir
macro incl(s: static string): untyped = quote do: include `s`
include "../build/configBase.nims"
include "../build/buildTasks.nims"

const flag = "--passC:-fno-allow-store-data-races"
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

proc check(want: bool) =
  nimFlags.setLen(0)
  setNimFlags()
  doAssert (flag in nimFlags) == want, nimCmdArgs
  if want: doAssert nimFlags[^1] == flag

try:
  # Use the real preprocessor with controlled compiler identity macros.
  var gcc: array[3, string]
  cpp = dir / "unused compiler"
  for i, v in [14, 15, 16]:
    gcc[i] = script(cmd & "-D__GNUC__=" & $v & " \"$@\"")
    cc = gcc[i]
    check(v >= 15)
    doAssert compilerInfo(nimFlags) == ("gcc", v)
  for mac in ["__INTEL_COMPILER", "__INTEL_LLVM_COMPILER", "__NVCOMPILER", "__PGI"]:
    cc = script(cmd & "-D__GNUC__=99 -D" & mac & "=1 \"$@\"")
    check(false)
  cc = script(cmd & "-D__GNUC__=99 -D__clang__=1 -D__clang_major__=23 \"$@\"")
  check(false)
  doAssert compilerInfo(nimFlags) == ("clang", 23)

  cc = gcc[0]
  cpp = gcc[1]
  check(false)
  ccDef = "cpp"
  check(true)
  ccDef = "cc"

  nimargs = @["--passC:-fallow-store-data-races"]
  setUserNimFlags(@["--gcc.exe:" & gcc[1].quoteShell])
  check(true)
  setUserNimFlags(@["--cc:clang", "--clang.exe:" & gcc[1].quoteShell])
  check(true)
  nimargs.setLen(0)
  setUserNimFlags(@[])

  cc = script("exec \"$OMPI_CC\" \"$@\"")
  envs = @["OMPI_CC=" & gcc[1]]
  check(true)
  setUserNimFlags(@["--putenv:" & ("OMPI_CC=" & gcc[0]).quoteShell])
  check(false)
  envs.setLen(0)
  setUserNimFlags(@[])

  cc = gcc[1]
  fo.debug = true
  check(true)
  cc = script("exit 1")
  var failed = false
  try:
    check(true)
  except IOError:
    failed = true
  doAssert failed
finally:
  rmDir(dir)

echo "Build flag checks passed."
