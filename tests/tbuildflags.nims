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

proc mpich(def: string): string =
  # MPICH wrappers run the compiler of -cc= (mpicc) or -cxx= (mpicxx), else their default.
  script("c=" & def.quoteShell & "; for a; do shift; case $a in -cc=*|-cxx=*) c=${a#*=};; " &
         "*) set -- \"$@\" \"$a\";; esac; done; exec \"$c\" \"$@\"")

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
    doAssert (compiler.name, compiler.major) == ("gcc", v)
  for mac in ["__INTEL_COMPILER", "__INTEL_LLVM_COMPILER", "__NVCOMPILER", "__PGI"]:
    cc = script(cmd & "-D__GNUC__=99 -D__clang__=1 -D__clang_major__=19 -D" & mac & "=1 \"$@\"")
    check()
  var clang: array[4, string]
  for i, v in [18, 19, 20, 23]:
    clang[i] = script(cmd & "-D__GNUC__=99 -D__clang__=1 -D__clang_major__=" & $v & " \"$@\"")
    cc = clang[i]
    check(if v == 19: clangFlag else: "")
    doAssert (compiler.name, compiler.major) == ("clang", v)
  # Offload drivers read the input in a device pass before the host pass.
  cc = script(findExe("cc").quoteShell & " -undef \"$@\" > /dev/null; " & cmd & "-D__GNUC__=15 \"$@\"")
  check(gccFlag)

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

    cc = mpich(unaffected)
    cpp = cc
    check()
    cflagsAlways = "-cc=" & affected.quoteShell & " -g"
    check(flag)
    check("", @["--gcc.options.always:" & ("-cc=" & unaffected.quoteShell).quoteShell])
    ccDef = "cpp"
    check()
    # GPU builds compile C++ as CUDA or HIP; the probe still reads host C++.
    cppflagsAlways = "-cxx=" & affected.quoteShell & " -x cuda"
    check(flag)
    ccDef = "cc"
    cflagsAlways = ""
    cppflagsAlways = ""

    cc = affected
    fo.debug = true
    check(flag)
    fo.debug = false

  # simd = "auto" follows the target macros of the active compile options.
  let native = script("case \"$*\" in *-march=native*) f='-D__AVX__=1 -D__AVX512F__=1 -D__AVX512DQ__=1';; *) f='';; esac; " &
                      cmd & "-D__GNUC__=14 $f \"$@\"")
  cc = native
  cflagsSpeed = "-O2 -march=native"
  cflagsDebug = "-Og"
  simd = "auto"
  check()
  doAssert compiler.simd == "SSE,AVX,AVX512", nimCmdArgs
  for d in ["--d:SSE", "--d:AVX", "--d:AVX512"]:
    doAssert d in nimFlags, nimCmdArgs
  # User flags follow the auto defines, so -u:AVX512 still removes a level.
  check("", @["-u:AVX512"])
  doAssert nimFlags.find("-u:AVX512") > nimFlags.find("--d:AVX512"), nimCmdArgs
  fo.debug = true
  check()
  doAssert "--d:AVX" notin nimFlags, nimCmdArgs
  fo.debug = false
  # The always options choose the compiler and can set the target.
  cc = mpich(gcc[1])
  check(gccFlag)
  doAssert compiler.simd == "", nimCmdArgs
  cflagsAlways = "-cc=" & native.quoteShell
  check()
  doAssert compiler.simd == "SSE,AVX,AVX512", nimCmdArgs
  cflagsSpeed = "-O2"
  cflagsAlways &= " -march=native"
  check()
  doAssert compiler.simd == "SSE,AVX,AVX512", nimCmdArgs
  cflagsAlways = ""
  cc = script(cmd & "-D__GNUC__=14 -D__AVX__=1 \"$@\"")
  check()
  doAssert "--d:AVX" in nimFlags and "--d:AVX512" notin nimFlags, nimCmdArgs
  cc = script(cmd & "-D__GNUC__=14 \"$@\"")
  check()
  doAssert "--d:SSE" notin nimFlags, nimCmdArgs
  simd = "SSE"
  check()
  doAssert "--d:SSE" in nimFlags and "--d:AVX" notin nimFlags, nimCmdArgs
  simd = ""

  cc = script("exit 1")
  var failed = false
  try:
    check()
  except IOError:
    failed = true
  doAssert failed
except Exception:
  # In a finally, rmDir would hide the failure behind the VM assertion
  # described in runBuilds.
  rmDir(dir)
  raise
rmDir(dir)

echo "Build flag checks passed."
