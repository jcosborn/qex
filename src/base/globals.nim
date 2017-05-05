import os
import strUtils
import stdUtils
import macros
#const profileEqns* = true
const profileEqns* = false

#var forceInline* {.compiletime.} = false
var forceInline* {.compiletime.} = true
macro setForceInline*(x:static[bool]):auto =
  forceInline = x
  result = newEmptyNode()

var staticUnroll* {.compiletime.} = false
#var staticUnroll* {.compiletime.} = true
macro setStaticUnroll*(x:static[bool]):auto =
  staticUnroll = x
  result = newEmptyNode()
when existsEnv("STATIC_UNROLL"):
  when getEnv("STATIC_UNROLL")=="1":
    setStaticUnroll(true)

var noAlias* {.compiletime.} = false
#var noAlias* {.compiletime.} = true
macro setNoAlias*(x:static[bool]):auto =
  noAlias = x
  result = newEmptyNode()

when existsEnv("VLEN"):
  const VLEN* = getEnv("VLEN").parseInt
else:
  const VLEN* = 8

static:
  echo "VLEN: ", VLEN
