# sets up build tasks, included from either build.nims or qex.nimble
# requires variables 'nim', 'qexDir' and 'nimArgs' to be declared before including
import strFormat

var configTasks = newSeq[tuple[cmd:string,desc:string,f:proc(){.nimcall.}]](0)
template configTask(name: untyped; description: string; body: untyped) =
  proc `name CTask`*() = body
  configTasks.add((cmd: astToStr(name), desc: description, f: `name CTask`))

var buildTasks = newSeq[tuple[cmd:string,desc:string,f:proc(){.nimcall.}]](0)
template buildTask(name: untyped; description: string; body: untyped) =
  proc `name BTask`*() = body
  buildTasks.add((cmd: astToStr(name), desc: description, f: `name BTask`))

var currentArg = ""  # set from parent files while iterating over arguments
proc getInt(): int =
  let t = split(currentArg,":")
  if t.len>=2: result = parseInt(t[1])
proc getString(): string =
  let t = split(currentArg,":")
  if t.len>=2: result = t[1]
  else: result = ""

var nimFlags: seq[string] = @[]
proc setNimFlags() =
  if nimFlags.len == 0:
    nimFlags = getNimFlags()

var debug = false
var run = false
var usecpp = false
var verbosity = -1
var bindir = "bin"
var srcPaths = @[".", "qex/src", "qex/tests"]  # use relative paths for convenience
if getCurrentDir() == qexDir: srcPaths = @["."]

proc findSrc(g: string): tuple[files:seq[string],dirs:seq[string]] =
  var fs = newSeq[string]()
  var ds = newSeq[string]()
  let d = getCurrentDir()
  for p in srcPaths:
    let f = staticExec &"cd {d}; ( find {p} -type f -ipath {p}*{g}; find {p} -type f -ipath {p}*{g}.nim ) |sort -u"
    if f != "":
      for t in f.splitLines:
        if t.endswith(".nim"):
          fs.add t
    let d = staticExec &"cd {d}; find {p} -type d -ipath *{g} |sort"
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
  let cc = if usecpp: "cpp" else: "c"
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
  let s = findSrc(g)
  if s.dirs.len != 0:
    if s.dirs.len > 1:
      echo "  Error: multiple directories match:"
      for d in s.dirs:
        echo "    ", d
      return true
    #echo "  process dir: ", s.dirs[0]
    for f in listFiles(s.dirs[0]):
      #echo "file: ", f
      if f.endsWith(".nim"):
        discard buildFile(f)
    return false
  if s.files.len != 0:
    if s.files.len > 1:
      echo "  Error: multiple files match:"
      for f in s.files:
        echo "    ", f
      return true
    return buildFile(s.files[0])


# === Config Tasks ===

configTask debug, "set debug build":
  debug = true

configTask run, "run executable after building":
  run = true

configTask cpp, "compile in cpp mode":
  usecpp = true

configTask verb, "set build verbosity to N (verb:N)":
  verbosity = getInt()


# === Build Tasks ===

proc echoCmds(tasklist: seq) =
  var clen = 0
  for t in tasklist:
    clen = max(clen, t.cmd.len)
  for t in tasklist:
    var first = true
    for l in t.desc.splitLines:
      if first:
        first = false
        let c = t.cmd & " ".repeat(clen-t.cmd.len)
        echo &"  {c}  {l}"
      else:
        echo " ".repeat(clen+4) & l

buildTask help, "show this help message":
  let s = '-'.repeat(70)
  echo s
  echo "QEX build script usage:"
  echo "  make [config commands] [build command | Nim source]"
  echo s
  echo "config commands:"
  echoCmds(configTasks)
  echo s
  echo "build commands:"
  echoCmds(buildTasks)
  echo s
  echo "Nim source:"
  echo "  foo.nim  Search for 'foo.nim' in source paths"
  echo "           (including subdirectories, but not following links),"
  echo "           compile and put binary in 'bin' directory."
  echo "  foo      First search for 'foo.nim' and compile if found."
  echo "           If foo.nim is not found, search for directory foo"
  echo "           and compile all *.nim files in it."
  #echo "           Specify part of path to resolve ambiguity."
  echo "  source paths:"
  for p in srcPaths:
    echo "    ", p
  quit 0

buildTask show, "show Nim compile flags":
  setNimFlags()
  echo "Nim flags:"
  echo join(nimFlags," ")

proc runTargets(f: string) =
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

buildTask targets, "show available build targets\ntargets:<name> will search for targets matching <name>\n(can include standard shell wildcards)":
  # targets
  # targets:foo
  # targets:foo*
  let f = getString()
  runTargets(f)


proc runClean() =
  echo "Cleaning nimcache directory: ", nimcache
  for f in nimcache.listFiles:
    #echo f
    #if f.endsWith(".o") or f.endsWith(".c") or f.endsWith(".cpp"):
    rmFile f

buildTask clean, "remove contents of nimcache directory\n("&nimcache&")":
  runClean()

let extraTests = [
  "gauge/wflow.nim",
  "examples/staghmc_sh.nim",
]

proc addTest(runscript:var seq[string], f, outdir:string) =
  let name = f.splitFile.name
  var rj = gorge("awk '$1==\"#RUNCMD\"{$1=\"\";print}' "&f)
  if rj == "": rj = "$RUNJOB"
  let exe = outdir/name
  discard buildFile(f, exe)
  let runner = qexDir/"tests/extra"/name/"run"
  if fileExists(runner): rj = runner
  runscript.add("echo Running: "&exe)
  runscript.add(rj&" "&exe&" || failed=\"$failed "&name&"\"")

proc buildTests() =
  var runscript = @["#!/bin/sh","$SETUPJOBS","failed=''"]
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

buildTask tests, "build tests":
  buildTests()
