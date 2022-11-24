# sets up build tasks, included from either build.nims or qex.nimble
# requires variables 'nim', 'qexDir' and 'nimArgs' to be declared before including
import strFormat

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
proc setNimFlags() =
  if nimFlags.len == 0:
    nimFlags = getNimFlags(fo)
    nimFlags.add userNimFlags

var run = false
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
  var nimcmd = tool & nim & " " & join(nimArgs," ") & " " & join(nimFlags," ")
  if run: nimcmd &= " -r "
  var (dir, name, ext) = splitFile(f)
  if outfile!="": name = outfile
  else:
    if not dirExists(bindir):
      mkDir(bindir)
    name = bindir / name
  #let cc = if usecpp: "cpp" else: "c"
  let cc = ccDef
  let s = nimcmd & " " & cc & " -o:" & name & " " & f
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
  let s = findSrc(g)
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
              (may need to proceed with '--' so make doesn't parse it).
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

let extraTests = [
  "gauge/wflow.nim",
  "examples/staghmc_sh.nim",
]

proc addTest(runscript:var seq[string], f, outdir:string) =
  let name = f.splitFile.name
  var rj = gorge("awk '$1==\"#RUNCMD\"{$1=\"\";print}' "&f).strip
  if rj == "": rj = "$RUNJOB"
  let exe = outdir/name
  discard buildFile(f, exe)
  let runner = qexDir/"tests/extra"/name/"run"
  if fileExists(runner): rj = runner
  runscript.add("echo Running: "&exe)
  runscript.add(rj&" "&exe&" || failed=\"$failed "&name&"\"")

proc buildTests() =
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
  for d in listDirs(qexDir/"tests"):
    let outdir = "tests"/splitPath(d)[1]
    if not dirExists(outdir):
      mkDir(outdir)
    for f in listFiles(d):
      #echo f
      let (dir, name, ext) = splitFile(f)
      #echo dir, " ", name, " ", ext
      if name[0]=='t' and ext==".nim":
        runscript.addTest(f, outdir)
  for f in extraTests:
    let outdir = bindir
    if not dirExists(outdir):
      mkDir outdir
    runscript.addTest(qexDir/"src"/f, outdir)
  #echo runscript.join("\n")
  runscript.add("$CLEANUPJOBS")
  runscript.add("if [ X != \"X$failed\" ];then echo Failed tests: $failed;exit 1;fi")
  runscript.add("echo $0: All tests passed")
  writeFile("testscript.sh", runscript.join("\n"))
  exec("chmod 755 testscript.sh")
  if dorun:
    exec "./testscript.sh"

let testsDesc = "  Build tests and create `testscript.sh' test runner"
buildTask tests, testsDesc:
  buildTests()

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
