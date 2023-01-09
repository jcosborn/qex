import macros
import strUtils
#import metaUtils
import os
#import qexInternal

var paramNames = newSeq[string](0)
var paramValues = newSeq[string](0)
var paramComments = newSeq[string](0)

proc addParam(s,r: string, c: string = "") =
  let i = paramNames.find(s)
  if i>=0:
    paramValues[i] = r
  else:
    paramNames.add s
    paramValues.add r
    paramComments.add c

proc addComment(s,c:string):string =
  result = s
  if c.len>0:
    let
      spc = "                                 ## "
      m = min(s.len, spc.len-8)
    result &= spc[m..^1] & c

template echoParams*(warnUnknown=false) =
  mixin echo
  bind addComment
  for i in 0..<paramNames.len:
    echo addComment(paramNames[i] & ": " & paramValues[i], paramComments[i])
  for i in 1..paramCount():
    var p = paramstr(i)
    if p[0] == '-':
      p = p[1..^1]
    let c = p.find(':')
    if c>=0:
      p = p[0..<c]
    if paramNames.find(p)<0:
      echo "Unknown argument: ",paramstr(i)

proc paramHelp*(p:string = ""):string =
  result = "Usage:\n  " & getAppFileName()
  let i = paramNames.find(p)
  if i>=0:
    result &= addComment(" -" & p & ":" & paramValues[i] & " (current value)", paramComments[i])
  else:
    result &= " -OPTION:VALUE ...\nAvailable OPTIONs and current VALUEs:"
    let spc = "                "
    for i in 0..<paramNames.len:
      let nm = paramNames[i]
      result &= "\n    " & (nm & spc[min(spc.len-1,nm.len)..^1] & " : " & paramValues[i]).addComment(paramComments[i])

template cnvnone(x:typed):untyped = x
template makeTypeParam(name,typ,deflt,cnvrt: untyped): untyped {.dirty.} =
  proc name*(s: string, d=deflt, c=""): typ =
    result = d
    let n = paramCount()
    for i in 1..n:
      let p = paramstr(i)
      if p.startsWith('-'&s&':'):
        let ll = s.len + 2
        result = cnvrt(p[ll..^1])
    addParam(s, $result, c)

makeTypeParam(intParam, int, 0, parseInt)
makeTypeParam(floatParam, float, 0.0, parseFloat)
makeTypeParam(strParam, string, "", cnvnone)
template stringParam*(x,y: untyped, c=""): untyped = strParam(x,y,c)

proc boolParam*(s: string, d = false, c=""): bool =
  result = d
  let n = paramCount()
  for i in 1..n:
    let p = paramstr(i)
    if p == '-'&s:
      result = true
    elif p.startsWith('-'&s&':'):
      let ll = s.len + 2
      let val = tolowerAscii(p[ll..^1])
      result = case val
        of "t","true","yes","y","on": true
        else: false
  addParam(s, $result, c)

proc intSeqParam*(s: string, d: seq[int] = @[], c=""): seq[int] =
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
  addParam(s, join(result," "), c)

proc floatSeqParam*(s: string, d: seq[float] = @[], c=""): seq[float] =
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
  addParam(s, join(result," "), c)

template setParam*(s:string, d:string, c:string=""):string = strParam(s,d,c)
template setParam*(s:string, d:int, c:string=""):int = intParam(s,d,c)
template setParam*(s:string, d:float, c:string=""):float = floatParam(s,d,c)
template setParam*(s:string, d:bool, c:string=""):bool = boolParam(s,d,c)
template setParam*(s:string, d:seq[int], c:string=""):seq[int] = intSeqParam(s,d,c)
template setParam*(s:string, d:seq[float], c:string=""):seq[float] = floatSeqParam(s,d,c)

macro letParam*(decls:untyped):auto =
  var
    empty = newStrLitNode("")
    comm = empty
  # echo decls.treerepr
  result = newNimNode(nnkLetSection, decls)
  for decl in decls:
    if decl.kind == nnkAsgn:
      result.add newIdentDefs(decl[0], newEmptyNode(), newCall("setParam", newLit($decl[0]), decl[1], comm))
      comm = empty
    elif decl.kind in CallNodes and decl.len == 2 and
        decl[1].kind == nnkStmtList and decl[1].len == 1 and
        decl[1][0].kind == nnkAsgn:
      result.add newIdentDefs(decl[0], newEmptyNode(),
        newCall(decl[1][0][0], newCall("setParam", newLit($decl[0]), decl[1][0][1], comm)))
      comm = empty
    elif decl.kind == nnkCommentStmt:
      comm = newStrLitNode($decl)
      result.add decl
    else:
      let li = decl.lineInfoObj
      error("letParam: syntax error: " &
        li.filename & ":" & $li.line & ":" & $li.column & "\n" & decl.lisprepr)
  # echo result.repr

template installHelpParam*(p="h") =
  if setParam(p, false, "Print the help message"):
    echo paramHelp()
    qexExit()

template assertParam*(p:auto, f:auto) =
  if not f p:
    qexError("assertion failure: " & astToStr(f(p)) & "\n" & paramHelp(astToStr p))

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
