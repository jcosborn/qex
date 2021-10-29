import macros
import strUtils
#import metaUtils
import os

var paramNames = newSeq[string](0)
var paramValues = newSeq[string](0)

proc addParam(s,r: string) =
  var i = paramNames.find(s)
  if i>=0:
    paramValues[i] = r
  else:
    paramNames.add s
    paramValues.add r

template echoParams*() =
  mixin echo
  for i in 0..<paramNames.len:
    echo paramNames[i], ": ", paramValues[i]

template cnvnone(x:typed):untyped = x
template makeTypeParam(name,typ,deflt,cnvrt: untyped): untyped {.dirty.} =
  proc name*(s: string, d=deflt): typ =
    result = d
    let n = paramCount()
    for i in 1..n:
      let p = paramstr(i)
      if p.startsWith('-'&s&':'):
        let ll = s.len + 2
        result = cnvrt(p[ll..^1])
    addParam(s, $result)

makeTypeParam(intParam, int, 0, parseInt)
makeTypeParam(floatParam, float, 0.0, parseFloat)
makeTypeParam(strParam, string, "", cnvnone)
template stringParam*(x,y: untyped): untyped = strParam(x,y)

proc intSeqParam*(s: string, d: seq[int] = @[]): seq[int] =
  result = d
  let n = paramCount()
  for i in 1..n:
    let p = paramstr(i)
    if p.startsWith('-'&s&':'):
      result.setLen(0)
      let ll = s.len + 2
      for c in split(p[ll..^1], ','):
        if c.len > 0:
          result.add parseInt(c)
  addParam(s, join(result," "))

proc floatSeqParam*(s: string, d: seq[float] = @[]): seq[float] =
  result = d
  let n = paramCount()
  for i in 1..n:
    let p = paramstr(i)
    if p.startsWith('-'&s&':'):
      result.setLen(0)
      let ll = s.len + 2
      for c in split(p[ll..^1], ','):
        if c.len > 0:
          result.add parseFloat(c)
  addParam(s, join(result," "))

template setParam*(s:string, d:string):string = strParam(s,d)
template setParam*(s:string, d:int):int = intParam(s,d)
template setParam*(s:string, d:float):float = floatParam(s,d)
template setParam*(s:string, d:seq[int]):seq[int] = intSeqParam(s,d)
template setParam*(s:string, d:seq[float]):seq[float] = floatSeqParam(s,d)

macro letParam*(decls:untyped):auto =
  #echo decls.treerepr
  result = newNimNode(nnkLetSection, decls)
  for decl in decls:
    if decl.kind == nnkAsgn:
      result.add newIdentDefs(decl[0], newEmptyNode(), newCall("setParam", newLit($decl[0]), decl[1]))
    elif decl.kind in CallNodes and decl.len == 2 and
        decl[1].kind == nnkStmtList and decl[1].len == 1 and
        decl[1][0].kind == nnkAsgn:
      result.add newIdentDefs(decl[0], newEmptyNode(),
        newCall(decl[1][0][0], newCall("setParam", newLit($decl[0]), decl[1][0][1])))
    elif decl.kind == nnkCommentStmt:
      result.add decl
    else:
      let li = decl.lineInfoObj
      error("letParam: syntax error: " &
        li.filename & ":" & $li.line & ":" & $li.column & "\n" & decl.lisprepr)
  #echo result.repr

template CLIset*(p:typed, n:untyped, prefix:string, runifset:untyped) =
  mixin echo
  let
    o = p.n
    s = prefix & astToStr(n)
  when compiles(strParam(s, p.n)):
    p.n = type(p.n)strParam(s, p.n)
  elif compiles(intParam(s, p.n)):
    p.n = type(p.n)intParam(s, p.n)
  elif compiles(floatParam(s, p.n)):
    p.n = type(p.n)floatParam(s, p.n)
  elif compiles(intSeqParam(s, p.n)):
    p.n = type(p.n)intSeqParam(s, p.n)
  elif compiles(floatSeqParam(s, p.n)):
    p.n = type(p.n)floatSeqParam(s, p.n)
  else:
    {.fatal:"Cannot set argument "&s&" of "&astToStr(p)&" for command line.".}
  if o != p.n:
    runifset
    echo "Customize $# : $# -> $#"%[s, $o, $p.n]
template CLIset*(p:typed, n:untyped, prefix = "") =
  p.CLIset n, prefix:
    discard
