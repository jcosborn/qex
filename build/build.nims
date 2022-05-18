import os, strUtils

var nim = paramStr(0)
var nimuserargs = newSeq[string](0)
var script = ""
var args = newSeq[string](0)
var qexDir = thisDir().parentDir

for i in 1..paramCount():
  let p = paramStr(i)
  if p[0]=='-':
    nimuserargs.add p
  elif p[0]==':':
    let t = p[1..^1]
    if t[0]=='-':
      nimuserargs.add t
    else:
      nimuserargs.add "-d:" & t
  else:
    if script=="": script = p
    else: args.add p

echo "Using Nim: ", nim
echo "Nim user args: ", nimuserargs
echo "Makefile script: ", script
echo "Script args: ", args
echo "QEX dir: ", qexDir

include "buildTasks.nims"

setUserNimFlags(nimuserargs)
let cmdargs = parseOpts(args)

if cmdargs.len == 0:
  runTask("help")
else:
  var t = getTask(cmdargs[0])
  if t.cmd == "":
    runMake(cmdargs)
  else:
    let cargs = cmdargs[1..^1]
    runTask(t, cargs)
