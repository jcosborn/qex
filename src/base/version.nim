import strutils, times

const sep = '='.repeat(78)
const head = "QEX Compilation information:"
const comptime = "  Time: " & staticExec("date")
const compmach = "  Host: " & staticExec("uname -a")
const gitlog = "  Log: " & staticExec("git log -1").indent(7).strip
const gitstat = "  Status: " & staticExec("git status -uno").indent(10).strip
const buildInfo = [sep,head,comptime,compmach,gitlog,gitstat,sep].join("\n")
const gitrev = staticExec("git rev-parse HEAD")

static: echo buildInfo

proc getBuildInfo*(): string =
  buildInfo

proc getVersion*(): string =
  gitrev

when isMainModule:
  echo getBuildInfo()
  echo getVersion()
