import os, strUtils, strFormat

var nim = paramStr(0)
#echo "Running: ", paramStr(0), " ", paramStr(1)

var qexDir = paramStr(1).splitfile[0].parentDir
echo "QEX directory: ", qexDir

var cd = getCurrentDir()
echo "Build directory: ", cd

proc makeSymLink(name,target: string) =
  var islink = staticExec(&"cd {cd}; [ -L {name} ] && echo true")
  if islink == "true":
    echo "Removing existing symlink: ", name
    var err = staticExec(&"cd {cd}; rm {name}")
    if err != "":
      echo err
      quit(1)
  if fileExists(name):
    #echo &"ERROR: file '{name}' exists"
    #echo "  please remove so a symlink can be created"
    #quit(1)
    echo &"WARNING: file '{name}' exists, skipping creation of symlink"
    return
  if dirExists(name):
    echo &"ERROR: directory '{name}' exists"
    echo "  please remove so a symlink can be created"
    quit(1)
  echo &"Creating symlink: {name} -> {target}"
  var err = staticExec(&"cd {cd}; ln -s {target} {name}")
  if err != "":
    echo err
    quit(1)

# create symlink 'qex' to QEX directory
if qexDir != "qex":
  makeSymLink("qex", qexDir)

# create symlink to qex.nimble
makeSymLink("qex.nimble", "qex/qex.nimble")

# create Makfile and Makefile.nims
var params = @[ "NIM", nim ]

proc confFile(fn: string) =
  var f = readFile(qexDir / "build" / fn & ".in")
  f = replace(f, "$", "!DOLLAR!")
  f = replace(f, "#", "!HASH!")
  f = replace(f, "@@", "$")
  f = f % params
  f = replace(f, "!HASH!", "#")
  f = replace(f, "!DOLLAR!", "$")
  echo "Creating file: ", fn
  writeFile(fn, f)

confFile("Makefile")
#confFile("Makefile.nims")
#copyFile("qex.nimble",qexDir/"qex.nimble")

# create qexconfig.nims
include "genconfig.nims"
