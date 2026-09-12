# sets up build tasks, included from either build.nims or qex.nimble
# requires variables 'nim', 'qexDir' and 'nimArgs' to be declared before including
import strFormat, tables

type
  Task* = tuple[cmd:string,desc:string,f:proc(){.nimcall.}]
template newTask(c,d: string, fn: typed): untyped = (cmd:c, desc:d, f:fn)
proc emptyTask(): Task = result.cmd = ""

var configTasks = newSeq[Task](0)
template configTask(name: untyped; description: string; body: untyped) =
  proc `name CTask`*() = body
  configTasks.add newTask(astToStr(name), description, `name CTask`)

var buildTasks = newSeq[Task](0)
template buildTask(name: untyped; description: string; body: untyped) =
  proc `name BTask`*() = body
  buildTasks.add newTask(astToStr(name), description, `name BTask`)

var currentArg = ""  # set from parent files while iterating over arguments
var remainingArgs: seq[string]
proc getInt(): int =
  let t = split(currentArg,":")
  if t.len>=2: result = parseInt(t[1])
proc getString(): string =
  let t = split(currentArg,":")
  if t.len>=2: result = t[1]
  else: result = ""

var fo = newFlagsOpts()
var userNimFlags: seq[string] = @[]
proc setUserNimFlags(x: seq[string]) =
  userNimFlags = x
var nimFlags: seq[string] = @[]
var nimCmdArgs = ""
#var extraFlags = ""
type Compiler = tuple[name: string, major: int]

proc compilerInfo(flags: seq[string]): Compiler =
  # The merged flags contain configuration settings followed by user overrides.
  var cfg = initTable[string, string]()
  for arg in flags:
    let s = arg.split(':', 1)
    if s.len < 2: continue
    let key = s[0].nimIdentNormalize
    let val = parseCmdLine(s[1]).join("")
    cfg[key] = val
    if key == "--putenv":
      let env = val.split('=', 1)
      putEnv(env[0], env[1])

  let typ = cfg.getOrDefault("--cc", ccType).nimIdentNormalize
  if typ notin ["gcc", "clang"]: return
  let pre = "--" & typ & (if ccDef == "cpp": ".cpp" else: "")
  let def = if ccDef == "cpp": (if typ == "gcc": "g++" else: "clang++") else: typ
  let exe = cfg.getOrDefault(pre & ".exe", def)
  let dir = cfg.getOrDefault(pre & ".path", cfg.getOrDefault("--" & typ & ".path"))
  let cmd = (if dir.len > 0: dir / exe else: exe).quoteShell
  let lang = if ccDef == "cpp": "c++" else: "c"
  # Probe the compiler behind MPI wrappers after applying environment/flag overrides.
  # Other compilers also define __GNUC__; their own macros distinguish them from GCC.
  let src = """
#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER) || defined(__NVCOMPILER) || defined(__PGI)
#elif defined(__clang__)
qex_cc clang __clang_major__
#elif defined(__GNUC__)
qex_cc gcc __GNUC__
#endif
"""
  let (outp, code) = gorgeEx(cmd & " -E -P -x " & lang & " -", src)
  if code != 0:
    raise newException(IOError, "Compiler probe failed for " & exe & ":\n" & outp)
  for line in outp.splitLines:
    let s = line.splitWhitespace
    if s.len == 3 and s[0] == "qex_cc":
      return (name: s[1], major: parseInt(s[2]))

proc compilerFlags(c: Compiler): seq[string] =
  if c.name == "gcc" and c.major >= 15:
    # Under -Ofast, GCC can turn omp master assignments into stores by every thread.
    # Disable -fallow-store-data-races so workers cannot overwrite updates with stale values.
    result.add "--passC:-fno-allow-store-data-races"

proc setNimFlags() =
  if nimFlags.len == 0:
    nimFlags = getNimFlags(fo)
    nimFlags.add userNimFlags
    nimFlags.add compilerFlags(compilerInfo(nimFlags))
  nimCmdArgs = join(nimArgs," ") & " " & join(nimFlags," ")
  #if extraFlags != "":
  #  nimCmdArgs &= " " & extraFlags

var run = false
var runArgs = ""
#var verbosity = -1
var bindir = "bin"
var srcPaths = @[".", "qex/src", "qex/tests"]  # use relative paths for convenience
if getCurrentDir() == qexDir: srcPaths = @["."]

proc findSrc(g: string): tuple[files:seq[string],dirs:seq[string]] =
  var fs = newSeq[string]()
  var ds = newSeq[string]()
  let d = getCurrentDir()
  for p in srcPaths:
    let c = &"cd {d}; ( find {p} -type f -ipath '*{g}'; find {p} -type f -ipath '*{g}.nim' ) |sort -u"
    let f = staticExec c
    if f != "":
      for t in f.splitLines:
        if t.endswith(".nim"):
          fs.add t
    let d = staticExec &"cd {d}; find {p} -type d -ipath '*{g}' |sort"
    if d != "":
      ds.add d.splitLines
  result = (files:fs, dirs:ds)

# return true if failed
proc buildFile(f: string, outfile=""): bool =
  setNimFlags()
  var tool = ""
  #tool = "valgrind "
  var nimcmd = tool & nim & " " & nimCmdArgs
  if run: nimcmd &= " -r "
  var (dir, name, ext) = splitFile(f)
  if outfile!="": name = outfile
  else:
    if not dirExists(bindir):
      mkDir(bindir)
    name = bindir / name
  #let cc = if usecpp: "cpp" else: "c"
  let cc = ccDef
  let s = nimcmd & " " & cc & " -o:" & name & " " & f & runArgs
  echo "running: ", s
  try:
    exec s
  except:
    echo "failed: ", s
    quit(-1)
  return false

# return true if failed
proc tryBuildSource(g: string): bool =
  result = true
  var s = findSrc(g)
  if ("."/g) in s.files:
    s.dirs = @[]
    s.files = @["."/g]
  let n = s.dirs.len + s.files.len
  if n > 1:
    echo "  Error: multiple targets match:"
    if s.dirs.len > 0:
      echo "    Directories:"
      for d in s.dirs:
        echo "      ", d
    if s.files.len > 0:
      echo "    Files:"
      for f in s.files:
        echo "      ", f
    return true
  if s.dirs.len == 1:
    echo "Processing directory: ", s.dirs[0]
    for f in listFiles(s.dirs[0]):
      if f.endsWith(".nim"):
        echo "Building source: ", f
        discard buildFile(f)
    return false
  if s.files.len == 1:
    echo "Building source: ", s.files[0]
    return buildFile(s.files[0])


# === Config Tasks ===

configTask cc, "compile in C mode":
  ccDef = "cc"

configTask cpp, "compile in C++ mode":
  ccDef = "cpp"

configTask debug, "set debug build":
  fo.debug = true

configTask run, "run executable after building":
  run = true

configTask verb, "set build verbosity to N (verb:N), N in 0,1,2,3":
  buildVerbosity = getInt()


# === Build Tasks ===

proc formatCmds(tasklist: seq): string =
  var s = newSeq[string](0)
  var clen = 0
  for t in tasklist:
    clen = max(clen, t.cmd.len)
  for t in tasklist:
    var first = true
    for l in t.desc.splitLines:
      if first:
        first = false
        let c = t.cmd & " ".repeat(clen-t.cmd.len)
        s.add &"  {c}  {l.strip}"
      else:
        s.add " ".repeat(clen+4) & l.strip
  s.join("\n")

let sepHelp = '-'.repeat(72)

let buildOptionsHelp = """
build options:
""" & formatCmds(configTasks)

let nimOptionsHelp = """
Nim options:
  -<option>   Passes '-<option>' to Nim compiler
              (may need to precede with '--' so make doesn't parse it).
  :-<option>  Passes '-<option>' to Nim compiler
              (avoids issues with make trying to parse it).
  :foo        Sets Nim define 'foo'
              (equivalent to '-d:foo').
  :foo=bar    Sets Nim define 'foo' to value 'bar'
              (equivalent to '-d:foo=bar')."""

var pathHelp = """
path:
  foo.nim  Search for file matching `*foo.nim' in source paths
           (including subdirectories, but not following links)
  foo      Search for both `*foo.nim' and `*foo',
           if a directory matches compile all `*.nim' in it
    Note:  only one match is allowed,
           specify part of path to resolve ambiguity
source paths:
"""
pathHelp &=  "  " & srcPaths.join("\n  ")

buildTask help, "   Show this help message":
  echo sepHelp
  echo "QEX build script usage:"
  echo "  make [command] [build option | Nim option]... [path]..."
  echo sepHelp
  echo "commands:"
  echo formatCmds(buildTasks)
  echo "           (command make is default and can be skipped)"
  echo sepHelp
  echo buildOptionsHelp
  echo sepHelp
  echo nimOptionsHelp
  echo sepHelp
  echo pathHelp
  echo sepHelp

buildTask depends, "Install Nimble dependencies":
  exec "nimble install -dy"

buildTask show, "   Show Nim compile flags":
  setNimFlags()
  echo "Nim compile command: ", ccDef
  echo "Nim flags:"
  echo join(nimFlags," ")

proc runTargets(f: string) =
  echo "Searching for targets matching: ", f
  if f == "":
    let d = getCurrentDir()
    for p in srcPaths:
      echo "targets in path: ", p
      let r = staticExec &"cd {d}; find {p} -name \\*.nim |sort"
      for l in r.splitLines:
        echo "  ", l
  else:
    let s = findSrc(f)
    if s.dirs.len == 0:
      echo "  No matching directories found"
    else:
      echo "  Directories:"
      for d in s.dirs:
        echo "    ", d
    if s.files.len == 0:
      echo "  No matching files found"
    else:
      echo "  Nim files:"
      for f in s.files:
        echo "    ", f

let targetsDesc = """Show available build targets
               targets <name> will search for targets matching <name>
               (can include standard shell wildcards)"""
buildTask targets, targetsDesc:
  var f = getString()
  if f == "" and remainingArgs.len>0: f = remainingArgs[0]
  runTargets(f)


let cleanDesc = """  Remove contents of nimcache directory
               ("""&nimcache&")"
buildTask clean, cleanDesc:
  echo "Cleaning nimcache directory: ", nimcache
  for f in nimcache.listFiles:
    #echo f
    #if f.endsWith(".o") or f.endsWith(".c") or f.endsWith(".cpp"):
    rmFile f

#let extraTests = [
#  "gauge/wflow.nim",
#  "examples/staghmc_sh.nim",
#]
#let extraTests = readFile(qexDir/"tests"/"extra"/"extra.txt").splitLines
#  .filterIt((it.len>0) and (not it.startsWith("#")))
#echo "Extra tests: ", extraTests2
var extraTests = newSeq[string]()
var extraTestArgs = newSeq[string]()
template extraTest(f:string, a="") =
  extraTests.add f
  extraTestArgs.add a
incl qexDir/"tests"/"extra"/"extra.nims"
#echo "Extra tests: ", extraTests
#echo "Extra args: ", extraTestArgs

var extraArgs = ""
proc addTest(runscript:var seq[string], f, outdir:string) =
  let name = f.splitFile.name
  var rj = gorge("awk '$1==\"#RUNCMD\"{$1=\"\";print}' "&f).strip
  if rj == "": rj = "$RUNJOB"
  let exe = outdir/name
  discard buildFile(f, exe)
  let runner = qexDir/"tests/extra"/"t"&name/"run"
  if fileExists(runner): rj = runner
  var args = exe
  if extraArgs != "": args &= " " & extraArgs
  runscript.add("echo Running: "&args)
  runscript.add(rj&" ./"&args&" || failed=\"$failed "&name&"\"")  # use " ./"&exe for /usr/bin/env

proc buildTests(scope = "") =
  var dirs: seq[tuple[src, dst: string]]
  let optional = scope != ""
  let script = if optional: "testscript-experimental.sh" else: "testscript.sh"
  if optional:
    for d in listDirs(qexDir/"src/experimental"):
      let name = splitPath(d)[1]
      if scope == "experimental" or scope == "experimental/" & name:
        let src = d/"tests"
        if dirExists(src):
          dirs.add (src, "tests/experimental"/name)
    if dirs.len == 0:
      echo "Error: no experimental test group matches: ", scope
      quit(1)
  else:
    for d in listDirs(qexDir/"tests"):
      dirs.add (d, "tests"/splitPath(d)[1])
  var runscript = @["#!/bin/sh",
                    "# Runs QEX tests and reports on failed tests",
                    "# Environment variables that can affect this script:",
                    "#   SETUPJOBS    commands to be run once at beginning of script",
                    "#   CLEANUPJOBS  commands to be run once at end of script",
                    "#   RUNJOB       command to launch test (can be multiple ranks)",
                    "#   RUN1         command to launch test on 1 rank",
                    "$SETUPJOBS","failed=''"]
  var dorun = run
  run = false
  if not dirExists("tests"):
    mkDir("tests")
  if optional and not dirExists("tests/experimental"):
    mkDir("tests/experimental")
  var count = 0
  for (d, outdir) in dirs:
    if not dirExists(outdir):
      mkDir(outdir)
    for f in listFiles(d):
      #echo f
      let (dir, name, ext) = splitFile(f)
      #echo dir, " ", name, " ", ext
      if name[0]=='t' and ext==".nim":
        inc count
        runscript.addTest(f, outdir)
  if optional and count == 0:
    echo "Error: no experimental tests found for: ", scope
    quit(1)
  if not optional:
    for i in 0..<extraTests.len:
      let f = extraTests[i]
      extraArgs = extraTestArgs[i]
      let outdir = bindir
      if not dirExists(outdir):
        mkDir outdir
      runscript.addTest(qexDir/"src"/f, outdir)
      extraArgs = ""
  #echo runscript.join("\n")
  runscript.add("$CLEANUPJOBS")
  runscript.add("if [ X != \"X$failed\" ];then echo Failed tests: $failed;exit 1;fi")
  runscript.add("echo $0: All tests passed")
  writeFile(script, runscript.join("\n"))
  exec("chmod 755 " & script)
  if dorun:
    exec "./" & script

let testsDesc = """  Build default tests and create `testscript.sh'
               tests experimental builds src/experimental/*/tests/t*.nim
               tests experimental/<group> selects one experimental group
               Experimental selections create `testscript-experimental.sh'"""
buildTask tests, testsDesc:
  if remainingArgs.len > 1:
    echo "Usage: tests [experimental[/<group>]]"
    quit(1)
  buildTests(if remainingArgs.len == 0: "" else: remainingArgs[0])

proc runMake(args: seq[string]) =
  for a in args:
    let failed = tryBuildSource(a)
    if failed:
      echo "Error: invalid source arg: ", a
      quit(1)

let makeDesc = """   Search for each [path]... as described below,
               compile, link, and put executables in `bin'"""
buildTask make, makeDesc:
  runMake(remainingArgs)

buildTask doc, "build inline docs":
  run = false
  let ccDefSave = ccDef
  ccDef = "doc --project --index:on --outdir:htmldocs"
  let file = "qex/src/qex.nim"
  #nim doc --project --index:on --git.url:<url> --git.commit:<tag> --outdir:htmldocs <main_filename>.nim
  #runArgs = " --project --index:on --outdir:htmldocs "
  discard buildFile(file, "htmldocs")

buildTask nbook, "build nimibook docs":
  setNimFlags()
  cd("qex/nbook")
  #let d = getCurrentDir()
  #nbook = d & "/qex/nbook/nbook.nim"
  #let nbook = d & "/nbook.nim"
  let nbook = "nbook.nim"
  run = true
  runArgs = " init"
  discard buildFile(nbook)
  runArgs = " build -d:nimibParallelBuild=false " & nimCmdArgs
  discard buildFile(nbook)

########

# parses options and returns command args
proc parseOpts(args: seq[string]): seq[string] =
  result.newSeq(0)
  var iarg = 0
  while iarg<args.len:
    currentArg = args[iarg]
    var found = false
    for t in configTasks:
      #echo t.cmd
      if currentArg.len>=t.cmd.len and currentArg[0..(t.cmd.len-1)] == t.cmd:
        echo "Processing config arg: ", currentArg
        found = true
        t.f()
    #if not found: break  # assume it is a build arg
    if not found:  # assume it is a build arg
      result.add args[iarg]
    inc iarg
  #while iarg<args.len:
  #  result.add args[iarg]
  #  inc iarg

proc getTask(name: string): Task =
  result = emptyTask()
  for t in buildTasks:
    if t.cmd == name:
      result = t

proc runTask(t: Task) =
  echo "Processing build task: ", t.cmd
  t.f()

proc runTask(name: string) =
  let t = getTask(name)
  if t.cmd == name:
    runTask(t)

proc runTask(t: Task, args: seq[string]) =
  echo "Processing build task: ", t.cmd
  remainingArgs = args
  t.f()

proc runTask(name: string, args: seq[string]) =
  let t = getTask(name)
  if t.cmd == name:
    runTask(t, args)
