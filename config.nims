import strUtils
import osPaths

# "key ~ val" sets "key" to "val"
template `~`(key,val:untyped):untyped =
  var v = astToStr(val)
  when compiles(val.type):
    when val.type is string:
      v = val
  switch(astToStr(key), v)
# "key ! val" sets "CC_TYPE.key" to "val"
template `!`(key,val:untyped):untyped =
  var v = astToStr(val)
  when compiles(val.type):
    when val.type is string:
      v = val
  switch(CC_TYPE & "." & astToStr(key), v)
proc def(s:string):bool =
  let t = "-d:" & s
  for i in 0..<paramCount():
    if paramStr(i) == t:
      result = true
template setVar(v,d:untyped):untyped =
  var v = d
  if existsEnv(v.astToStr):
    v = getEnv(v.astToStr)
  echo "using " & v.astToStr & ": ", v

setVar(QEXDIR, ".")
switch("path", QEXDIR / "src")

#let home = getEnv("HOME")
#setVar(QMPDIR, home & "/lqcd/install/qmp")
#setVar(QIODIR, home & "/lqcd/install/qio")

setVar(CC_TYPE, "gcc")
setVar(CC, "mpicc")
setVar(LD, CC)
setVar(CFLAGS_ALWAYS, "-Wall -std=gnu99 -D_C99 -D_HAS_C9X")
setVar(CFLAGS_DEBUG, "-g3 -O0")
setVar(CFLAGS_SPEED, "-O3 -march=native")
setVar(VERBOSITY, "1")
setVar(SIMD, "")
setVar(LDFLAGS, CFLAGS_ALWAYS)

var vlenS = 1
var vlenD = 1
if SIMD.contains("QPX"):
  vlenS = 4
  vlenD = 4
if SIMD.contains("SSE"):
  vlenS = 4
  vlenD = 2
if SIMD.contains("AVX"):
  vlenS = 8
  vlenD = 4
if SIMD.contains("AVX512"):
  vlenS = 16
  vlenD = 8
setVar(VLENS, $vlenS)
setVar(VLEND, $vlenD)
switch("putenv", "VLENS=" & VLENS)
switch("putenv", "VLEND=" & VLEND)
setVar(VLEN, $vlenS)
switch("putenv", "VLEN=" & VLEN)
#switch("warning[SmallLshouldNotBeUsed]", "off")

cc ~ CC_TYPE
exe ! CC
linkerexe ! LD
options.always ! CFLAGS_ALWAYS
options.debug ! CFLAGS_DEBUG
options.speed ! CFLAGS_SPEED
options.linker ! LDFLAGS

threads ~ on
tlsEmulation ~ off
verbosity ~ VERBOSITY

when not def "debug":
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
else:
  echo "debug build"

let ss = SIMD.split(',')
if ss.len>0:
  for s in items(ss):
    case s
    of "QPX":
      d ~ QPX
    of "SSE":
      d ~ SSE
    of "AVX":
      d ~ AVX
    of "AVX512":
      d ~ AVX512
    else: discard

discard """
d ~ simdS1
d ~ simdS2
d ~ simdS4
d ~ simdS8
d ~ simdS16
d ~ simdD1
d ~ simdD2
d ~ simdD4
d ~ simdD8
d ~ simdD16
"""
