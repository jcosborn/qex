import os
import strUtils
import stdUtils
#const profileEqns* = true
const profileEqns* = false

when existsEnv("VLEN"):
  const VLEN* = getEnv("VLEN").parseInt
else:
  const VLEN* = 8

echoImm "VLEN: ", VLEN
