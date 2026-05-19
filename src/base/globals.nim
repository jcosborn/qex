import os
import strUtils
#import stdUtils
import macros

const profileEqnsInt {.intdefine.} = 1
when profileEqnsInt == 0:
  const profileEqns* = false
else:
  const profileEqns* = true

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

var VLENmax {.compileTime.} = 256
var VLEN* {.compileTime.} = min(VLENmax, 8)
macro setVLEN*(n: static[int]): auto =
  VLEN = min(VLENmax, n)
  echo "VLEN: ", VLEN
  result = newEmptyNode()
macro setVLENmax*(n: static[int]): auto =
  VLENmax = n
  VLEN = min(VLENmax, VLEN)
  echo "VLENmax: ", VLENmax
  echo "VLEN: ", VLEN
  result = newEmptyNode()
when existsEnv("VLEN"):
  setVLEN(getEnv("VLEN").parseInt)
else:
  static:
    echo "VLEN: ", VLEN

var defaultNc {.compiletime.} = 3
macro setDefaultNc*(n: static[int]): untyped =
  defaultNc = n
  result = newEmptyNode()
macro getDefaultNc*(): untyped =
  return newLit(defaultNc)

template getDefPrecStr:string =
  const defPrec {.strdefine.} = "D"
  defPrec

var defPrec* {.compiletime.} = getDefPrecStr()
macro setDefaultSingle* =
  defPrec = "S"
  echo "Default precision: ", defPrec
static: echo "Default precision: ", defPrec
