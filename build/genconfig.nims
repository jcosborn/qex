import osPaths, strUtils, strformat, tables, macros

var args = initTable[string,string]()
var envargs = newSeq[string]()
proc setKey(k: string, v: string) =
  if k == "env":
    envargs.add v
  else:
    args[k] = v

proc fromEnv(key: string, env: string) =
  if existsEnv(env):
    args[key] = getEnv(env)

fromEnv("qmpdir", "QMPDIR")
fromEnv("qiodir", "QIODIR")
fromEnv("qudadir", "QUDADIR")
fromEnv("cudalibdir", "CUDALIBDIR")
fromEnv("chromadir", "CHROMADIR")
fromEnv("griddir", "GRIDDIR")

for i in 2..paramCount():
  let p = paramStr(i)
  let s = p.split(':')
  let k = nimIdentNormalize s[0]
  case s.len
  of 1:
    setKey(k, "true")
  of 2:
    setKey(k, s[1])
  else:
    setKey(k, join(s[1..s.high],":"))
  #echo p, " ", k, " ", args[k]

#var envarg = envargs.join(" ")
#if envarg != "":
#  args["envs"]
if envargs.len > 0:
  args["envs"] = $envargs

include "configBase.nims"
include "configDefault.nims"
const cf = "configDefault.nims"
const ls = cf.staticRead
var c = newSeq[string](0)

template kv(k: string, v: string): untyped =
  k & " = " & v
template kv(k: string, v: string, p: int): untyped =
  k & " = " & v
template kv(k: string, v: string, p: string): untyped =
  k & " = \"" & v & "\""
template kv(k: string, v: string, p: seq[string]): untyped =
  k & " = " & v

# check if symbol was initialized from another symbol
macro implSym(x: typed): untyped =
  let i = x.getImpl
  #echo i.kind
  #echo i[2].kind
  return newLit i[2].kind == nnkSym

template process1(k: string, def: string, p: typed, l: string): untyped =
  let n = k.nimIdentNormalize
  if n in args:
    #echo def
    kv(k, args[n], p)
  else:
    if p.implSym:
      kv(k, def)
    else:
      kv(k, $p, p)

macro process(ls: static string): untyped =
  result = newStmtList()
  for l in ls.splitLines:
    #echo l
    let s = l.split('=')
    if s.len == 1:
      result.add quote do:
        c.add `l`
    else:
      let k = s[0].strip
      let v = join(s[1..^1],"=").strip
      let p = parseExpr(k)
      result.add quote do:
        c.add process1(`k`, `v`, `p`, `l`)

process(ls)

writeFile("qexconfig.nims", c.join "\n")
