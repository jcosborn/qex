#import os
import strFormat, macros

# workaround for limitation of include
macro incl(s: static string): untyped = quote do: include `s`

echo &"Including: {thisDir()}/configBase.nims"
include "configBase.nims"

const qc = getCurrentDir() / "qexconfig.nims"
when fileExists(qc):
  echo "Including: ", qc
  incl qc
else:
  echo "Not found: ", qc

echo &"Including: {thisDir()}/build.nims"
include "build.nims"
