import os, macros

var nim = paramStr(0)
const qexDir = thisDir().parentDir
macro incl(s: static string): untyped = quote do: include `s`
include "../build/configBase.nims"
include "../build/buildTasks.nims"

let (tmp, code) = gorgeEx("mktemp -d")
doAssert code == 0, tmp
let dir = tmp.strip
let log = dir / "log"

try:
  nim = dir / "fake nim"
  writeFile(nim, """#!/bin/sh
set -eu
cache=
out=
run=0
for arg do
  case "$arg" in
    --nimcache:*) cache=${arg#--nimcache:} ;;
    -o:*) out=${arg#-o:} ;;
    -r) run=1 ;;
  esac
done
name=${out##*/}
mkdir -p "$cache"
mkdir "$cache/active"
trap 'rmdir "$cache/active"' 0
printf 'start|%s|%s\n' "$name" "$cache" >> "$QEX_BUILD_LOG"
touch "$QEX_BUILD_SYNC/$name"
if [ "$QEX_BUILD_PAIR" = 1 ] && [ "$name" != c ]; then
  # Both jobs must start before either finishes.
  i=0
  until [ -f "$QEX_BUILD_SYNC/a" ] && [ -f "$QEX_BUILD_SYNC/b" ]; do
    i=$((i + 1))
    [ "$i" -lt 100 ] || exit 5
    sleep 0.05
  done
fi
if [ -n "$QEX_BUILD_FAIL" ] && [ "$name" != "$QEX_BUILD_FAIL" ]; then
  sleep 0.1
fi
printf 'end|%s|%s\n' "$name" "$cache" >> "$QEX_BUILD_LOG"
[ "$name" != "$QEX_BUILD_FAIL" ] || exit 3
touch "$out"
if [ "$run" = 1 ]; then
  printf 'run|%s|%s\n' "$name" "$cache" >> "$QEX_BUILD_LOG"
fi
""")
  exec "chmod +x " & nim.quoteShell
  nimcache = dir / "cache with spaces"
  cc = findExe("cc")
  bindir = dir / "bin"
  srcPaths = @[dir / "src"]
  mkDir(srcPaths[0])
  for name in ["a", "b", "c"]:
    writeFile(srcPaths[0] / name & ".nim", "discard\n")
  putEnv("QEX_BUILD_LOG", log)
  putEnv("QEX_BUILD_SYNC", dir / "sync")

  for (n, dorun, fail, opt) in [(1, false, "", ""), (2, false, "", ""),
                               (2, false, "a", ""), (2, false, "b", ""),
                               (2, false, "c", ""), (2, true, "", ""),
                               (2, false, "", "--nimcache:"), (2, false, "", "--nimcache=")]:
    discard parseOpts(@["jobs:" & $n])
    run = dorun
    if dirExists(bindir): rmDir(bindir)
    if dirExists(dir / "sync"): rmDir(dir / "sync")
    mkDir(dir / "sync")
    writeFile(log, "")
    putEnv("QEX_BUILD_PAIR", if n > 1 and not run: "1" else: "0")
    putEnv("QEX_BUILD_FAIL", fail)
    let cache = if opt == "": nimcache else: dir / "override cache"
    setUserNimFlags(if opt == "": @[] else: @[opt & cache.quoteShell])
    nimFlags.setLen(0)
    var failed = false
    try:
      runMake(@["a.nim", "b.nim", "c.nim"])
    except OSError:
      failed = true
    doAssert failed == (fail != ""), readFile(log)
    var active, peak, count, ran: int
    var caches: seq[string]
    for line in readFile(log).splitLines:
      if line == "": continue
      let s = line.split('|')
      case s[0]
      of "start":
        inc active
        inc count
        peak = max(peak, active)
        if n > 1 and not run:
          doAssert s[2] == cache / (if s[1] == "b": "job-1" else: "job-0"), line
        else:
          doAssert s[2] == cache, line
        doAssert s[2] notin caches, line
        caches.add s[2]
      of "end":
        dec active
        caches.delete(caches.find(s[2]))
      of "run": inc ran
      else: doAssert false, line
    doAssert active == 0
    doAssert peak == (if n > 1 and not run: 2 else: 1), readFile(log)
    doAssert count == (if fail in ["a", "b"]: 2 else: 3), readFile(log)
    doAssert ran == (if run: 3 else: 0), readFile(log)
    for name in ["a", "b", "c"]:
      let built = name != fail and (name != "c" or fail notin ["a", "b"])
      doAssert fileExists(bindir / name) == built, name

  for opt in ["jobs", "jobs:0", "jobs:-1"]:
    var failed = false
    try:
      discard parseOpts(@[opt])
    except ValueError:
      failed = true
    doAssert failed, opt
  runTask("clean")
  doAssert listDirs(nimcache).len == 0
  doAssert listFiles(nimcache).len == 0
finally:
  rmDir(dir)

echo "Build job checks passed."
