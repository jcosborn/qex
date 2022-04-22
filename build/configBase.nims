import os, strUtils

var
  nimcache = getCurrentDir() / "nimcache"
  buildVerbosity = 0

  ccType = "gcc"

  cc = "gcc"
  cflagsAlways = ""
  cflagsDebug = ""
  cflagsSpeed = ""
  ld = cc
  ldflags = cflagsAlways

  cpp = "g++"
  cppflagsAlways = ""
  cppflagsDebug = ""
  cppflagsSpeed = ""
  ldpp = cpp
  ldppflags = cppflagsAlways

#  ompflags = ""

  simd = ""
  vlen = "8"

  qmpDir = ""
  qioDir = ""
  qudaDir = ""
  cudaLibDir = ""
  primmeDir = ""
  chromaDir = ""
  gridDir = ""

  envs = newSeq[string]()

  FUELCompat = false


type
  flagsOpts* = object
    debug*: bool

proc getNimFlags*(fo: flagsOpts): seq[string] =
  var defargs = newSeq[string](0)
  # "set(key, val)" sets "key" to "val"
  template set(key: string, val: untyped) =
    var v = astToStr(val)
    when compiles(val.type):
      when val.type is string:
        v = val
      when val.type is int:
        v = $val
    echo "setting: ", key, " <- \"", v, "\""
    when true: # fixes Nim 1.5.1 regression  #getCommand()=="e":
      if v=="":
        defargs.add "--" & key
      else:
        defargs.add "--" & key & ":\"" & v & "\""
    else:
      if key.len>=7 and key[0..6]!="warning":  # warnings don't seem to work here
        switch(key, v)
  # "key ~ val" sets "key" to "val"
  template `~`(key,val: untyped) =
    set(astToStr(key), val)
  # "key ! val" sets "ccType.key" to "val"
  template `!`(key,val: untyped) =
    set(ccType & "." & astToStr(key), val)

  path ~ "qex/src"

  cc ~ ccType
  exe ! cc
  linkerexe ! ld
  options.always ! cflagsAlways
  options.debug ! cflagsDebug
  options.speed ! cflagsSpeed
  options.linker ! ldflags
  cpp.exe ! cpp
  cpp.linkerexe ! ldpp
  cpp.options.always ! cppflagsAlways
  cpp.options.debug ! cppflagsDebug
  cpp.options.speed ! cppflagsSpeed
  cpp.options.linker ! ldppflags

  #putenv ~ ("OMPFLAG=" & ompflags)
  if qmpDir != "":
    d ~ ("qmpDir:" & qmpDir)
  if qioDir != "":
    d ~ ("qioDir:" & qioDir)
  if qudaDir != "":
    d ~ ("qudaDir:" & qudaDir)
    d ~ ("cudaLibDir:" & cudaLibDir)
  if primmeDir != "":
    d ~ ("primmeDir:" & primmeDir)
  if chromaDir != "":
    d ~ ("chromaDir:" & chromaDir)
  if gridDir != "":
    d ~ ("gridDir:" & gridDir)

  if existsenv("FUELCompat") and getenv("FUELCompat")!="0":
    d ~ "FUELCompat"

  threads ~ on
  tlsEmulation ~ off
  verbosity ~ buildVerbosity
  nimcache ~ nimcache
  warning[SmallLshouldNotBeUsed] ~ off
  embedsrc ~ ""

  if not fo.debug:
    d ~ "release"
    d ~ "danger"
    #obj_checks ~ off
    #field_checks ~ off
    #range_checks ~ off
    #bound_checks ~ off
    #overflow_checks ~ off
    #nilchecks ~ off
    #assertions ~ off
    #stacktrace ~ off
    #linetrace ~ off
    #debugger ~ off
    #line_dir ~ off
    #dead_code_elim ~ on
    #panics ~ on
    opt ~ speed
  else:
    echo "debug build"

  let ss = simd.split(',')
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

  putenv ~ ("VLEN=" & vlen)
  for e in envs:
    putenv ~ e
    let t = e.split("=")
    putenv(t[0],join(t[1..^1]))
  return defargs

#  echo "Finished config file: ", thisDir(), "/config.nims"
