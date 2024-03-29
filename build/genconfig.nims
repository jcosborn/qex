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

template kvc(k: string, v: string, c: string): untyped =
  k & " = " & v & c
template kvc(k: string, v: string, p: int, c: string): untyped =
  k & " = " & v & c
template kvc(k: string, v: string, p: string, c: string): untyped =
  k & " = \"" & v & "\"" & c
template kvc(k: string, v: string, p: seq[string], c: string): untyped =
  k & " = " & v & c

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

template process1(k: string, def: string, p: typed, l: string, c: string): untyped =
  let n = k.nimIdentNormalize
  if n in args:
    #echo def
    kvc(k, args[n], p, c)
  else:
    if p.implSym:
      kvc(k, def, c)
    else:
      kvc(k, $p, p, c)

macro process(ls: static string): untyped =
  result = newStmtList()
  for l in ls.splitLines:
    #echo l
    let s = l.split('=')
    if s.len == 1 or s[0][0] == '#':
      result.add quote do:
        c.add `l`
    else:
      #echo l
      let k = s[0].strip  # key
      let v0 = join(s[1..^1],"=").strip  # value (with possible comment)
      let vs = v0.split('#')
      let v = vs[0].strip  # value (without comment)
      let p = parseExpr(k)
      if vs.len == 1:
        result.add quote do:
          c.add process1(`k`, `v`, `p`, `l`)
      else:
        let cmt = "  # " & join(vs[1..^1],"#").strip  # comment
        result.add quote do:
          c.add process1(`k`, `v`, `p`, `l`, `cmt`)

process(ls)

writeFile("qexconfig.nims", c.join "\n")
