import os, strUtils

var nim = paramStr(0)
var nimargs = newSeq[string](0)
var script = ""
var args = newSeq[string](0)
var qexDir = thisDir().parentDir

for i in 1..paramCount():
  let p = paramStr(i)
  if p[0]=='-':
    nimargs.add p
  else:
    if script=="": script = p
    else: args.add p

echo "Using Nim: ", nim
echo "Nim user args: ", nimargs
echo "Makefile script: ", script
echo "Script args: ", args
echo "QEX dir: ", qexDir

include "buildTasks.nims"

var iarg = 0
while iarg<args.len:
  currentArg = args[iarg]
  var found = false
  for t in configTasks:
    #echo t.cmd
    if args[iarg].len>=t.cmd.len and args[iarg][0..(t.cmd.len-1)] == t.cmd:
      echo "Processing config arg: ", currentArg
      found = true
      t.f()
  if not found: break  # assume it is a build arg
  inc iarg

while iarg<args.len:
  currentArg = args[iarg]
  echo "Processing build arg: ", currentArg
  var found = false
  for t in buildTasks:
    if args[iarg].len>=t.cmd.len and args[iarg][0..(t.cmd.len-1)] == t.cmd:
      found = true
      t.f()
  if not found: # check if source
    let failed = tryBuildSource(currentArg)
    if failed:
      echo "Error: invalid build arg: ", currentArg
      quit(1)
  inc iarg
