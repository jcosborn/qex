include "local.nims"

import ospaths, sequtils, strutils

# Package

version       = "0.0.0"
author        = "James C. Osborn"
description   = "Quantum EXpressions lattice field theory framework"
license       = "MIT"
srcDir        = qexDir/"src"

# Dependencies

requires "nim >= 0.16.0"
when declared(primmeDir):
  requires "primme >= 0.1.1"

type NamePath = tuple[n,p:string]
proc targets(p:string):seq[NamePath] =
  p.listFiles.filterIt(it.splitFile.ext==".nim").mapIt((it.splitFile.name,it))
proc recTargets(ps:seq[string]):seq[NamePath]
proc recTargets(p:string):seq[NamePath] = p.targets & p.listDirs.recTargets
proc recTargets(ps:seq[string]):seq[NamePath] =
  result = @[]
  for p in ps: result &= p.recTargets
proc findTarget(ts:seq[NamePath], t:string):NamePath =
  var i = 0
  let tp = t.splitPath
  if tp.head.len == 0:
    let n = t.splitFile.name
    while i < ts.len:
      if ts[i].n == n: break
      inc i
  else:
    while i < ts.len:
      if t.endsWith(".nim") and ts[i].p.endsWith(t): break
      if ts[i].p.endsWith(t&".nim"): break
      inc i
  if i < ts.len: return ts[i]
  else:
    echo "Error: cannot find target: `", t, "'"
    quit QuitFailure

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
template `~`(k,v:untyped) = set(astToStr(k), v)
template `!`(k,v:untyped) = set((ccType&"."&astToStr(k)), v)
template def(v:untyped) = define ~ (astToStr(v)&"="&v)
template optDef(v:untyped) =
  when declared(v): def v

proc setup =
  path ~ srcDir
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
  for s in simd.split(" "): define ~ s
  for d in extraDef:
    putenv ~ d  # We'll convert env to def in src soon.
    define ~ d
  # Here are optional external dependencies.
  optDef primmeDir
  optDef lapackLib
  optDef qudaDir     # This and cuda below are not in effect.
  optDef cudaLibDir  # Will convert them later.

proc setupRelease =
  define ~ "release"
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

task make, "compile, link, and put executables in `bin'":
  const c = paramCount()
  var debug = false
  var
    ts = newseq[string]()
    args = newseq[string]()
    defs = newseq[string]()
  for i in 1..c:
    let pi = paramStr i
    if pi[0] == '-': args.add pi
    elif ts.len == 0 and pi == "make": continue
    elif pi == "debug": debug = true
    elif '=' in pi: defs.add pi
    else: ts.add pi
  if ts.len == 0:
      exec(paramStr(0)&" help")
  elif ts.len > 1:
    for i in 0..<ts.len:
      exec(paramStr(0)&" make "&args.join(" ")&" "&ts[i]&" "&defs.join(" "))
  else:
    let (name,target) = (extraSrcDir&qexDir).recTargets.findTarget ts[0]
    setup()
    for d in defs: define ~ d
    for a in args:
      var i = 0
      while a[i] == '-': inc i
      var j = i
      while j < a.len and a[j] != ':': inc j
      if j < a.len: a[i..<j].switch a[j+1..^1]
      else: a[i..<j].switch
    if not dirExists("bin"): mkDir"bin"
    "out".set("bin/"&name)
    if debug: echo "debug build"
    else: setupRelease()
    setCommand "c", target

task targets, "List available targets":
  let ts = (extraSrcDir&qexDir).recTargets
  for t in ts: echo t.n,spaces(32-t.n.len),"\t",t.p

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

task clean, "remove temporary build files":
  rmDir "nimcache"
