import macros, os, strUtils, tables, algorithm, seqUtils

# TODO
# include CT in saveParams RT
# check used
# only write from rank 0
# only read from rank 0
# add way to represent empty seq/string (versus unset/unknown)
# sort paramHelp output by file and option name

type
  ParamListT = object  # used to store set parameters from cmd line or file
    value: string
    used: bool
    source: string  # command line or file
  ParamObj = object  # used to store CT/RT parameter info
    name: string
    value: string
    comment: string
    file: string
    line: int
  IiType = typeof instantiationInfo()

var inited = false
var loadParamsCommand {.compileTime.} = "loadParams"  # normally not changed
var params: Table[string,ParamListT]  # store params during setup
# param info generated at compile time, as setParam commands are compiled
var paramInfoCT {.compileTime.} = newSeq[ParamObj](0)
#var paramInfoCTRT = newSeq[ParamObj](0)  # conversion of paramInfoCT to RT variable
var paramInfoCTRT: seq[ParamObj]  # conversion of paramInfoCT to RT variable
# param info generated at run time, as setParam commands executed (using set values)
var paramInfo = newSeq[ParamObj](0)

proc getParamValue(s: string): string =
  params[s].used = true
  params[s].value

proc loadParams(fn: string) =
  mixin echo
  echo "Loading params from file: ", fn
  for l in fn.lines:
    #echo l
    let s = l.split('#',maxsplit=1)[0].strip
    if s != "":
      #echo s
      let kv = s.splitWhitespace(maxsplit=1)
      let k = kv[0]
      let v = if kv.len>1: kv[1] else: ""
      params[k] = ParamListT(value:v, used:false, source:fn)

proc setupParams =
  let n = paramCount()
  for i in 1..n:
    let p = paramstr(i)
    let s = p.split(':', maxsplit=1)
    if s[0][0] == '-':  # only parse options starting with -
      let k = s[0][1..^1]
      let v = if s.len>1: s[1] else: ""
      if k == loadParamsCommand:
        loadParams(v)
      else:
        params[k] = ParamListT(value:v, used:false, source:"cmdline")

template setup =
  if not inited:
    inited = true
    setupParams()

proc set(x: var seq[ParamObj], s,d,c,f: string, l: int) =
  let t = ParamObj(name:s,value:d,comment:c,file:f,line:l)
  var i = 0
  while(i<x.len and x[i].name!=s): inc i
  if i<x.len:
    x[i] = t
  else:
    x.add t

proc find(x: seq[ParamObj], s: string): int =
  result = -1
  for i in 0..<x.len:
    if x[i].name == s:
      result = i
      break

macro addParamCTX(s0: string, d0: auto, c0: string, ii: static IiType): auto =
  #echo s0.treerepr
  var s = ""
  if s0.kind == nnkStrLit:
    s = s0.strVal
  elif s0.kind == nnkSym and s0.getImpl.kind == nnkStrLit:
    s = s0.getImpl.strVal
  var d = ""
  if d0.kind == nnkStrLit:
    d = d0.strVal
  var c = ""
  if c0.kind == nnkStrLit:
    c = c0.strVal
  #echo s, " ", d, " ", c
  if s != "":
    paramInfoCT.set(s, d, c, ii.filename, ii.line)
  result = newEmptyNode()
macro defString(x0: bool): auto =
  #echo "bool: ", x0.treerepr
  result = newLit "#false"
  var x = x0
  while x.kind == nnkSym:
    #echo "bool: ", x.treerepr
    case x.repr
    of "true": result = newLit "true"; break
    of "false": result = newLit "false"; break
    else: x = x.getImpl[2]  # getImpl returns IdentDefs
macro defString(x: static bool): auto =
  result = newLit $x
macro defString(x: string): auto =
  #echo "string: ", x.treerepr
  result = if x.kind == nnkStrLit: x else: newLit "#string"
macro defString(x: static string): auto =
  result = newLit x
macro defString(x: int): auto =
  result = if x.kind == nnkIntLit: newLit x.repr else: newLit"#0"
macro defString(x: static int): auto =
  result = newLit $x
macro defString(x: float): auto =
  #echo "float param: ", x.treerepr
  result = if x.kind == nnkFloatLit: newLit x.repr else: newLit"#0.0"
macro defString(x: static float): auto =
  result = newLit $x
macro defString(x: seq[int]): auto =
  #echo "seq[int]: ", x.treerepr
  result = newLit"#0,0,..."
  if x.kind == nnkPrefix and x[0].repr == "@" and x[1].kind == nnkBracket:
    if x[1].len > 1:
      var s = x[1][0].repr
      for i in 1..<x[1].len:
        s &= "," & x[1][i].repr
      result = newLit s
macro defString(x: static seq[int]): auto =
  let s = x.mapIt($it).join(",")
  return newLit s
macro defString(x: seq[float]): auto =
  #echo "non static"
  #echo "seq[float]: ", x.treerepr
  result = newLit"#0.0,0.0,..."
  if x.kind == nnkPrefix and x[0].repr == "@" and x[1].kind == nnkBracket:
    if x[1].len > 1:
      var s = x[1][0].repr
      for i in 1..<x[1].len:
        s &= "," & x[1][i].repr
      result = newLit s
macro defString(x: static seq[float]): auto =
  let s = x.mapIt($it).join(",")
  return newLit s
template addParamCT(s: string, d: auto, c: string, ii: static IiType) =
  addParamCTX(s, defString(d), c, ii)

proc format(x: seq[ParamObj], linePrefix = "#line "): string =
  var nameLen, valueLen, commLen, fileLen, lineLen = 0
  template maxeq(r: var int, x: int) =  r = max(r,x)
  for t in x:
    nameLen.maxeq t.name.len
    valueLen.maxeq t.value.len
    commLen.maxeq t.comment.len
    fileLen.maxeq t.file.len
    lineLen.maxeq (linePrefix & $t.line).len
  result = ""
  template pad(s: string, n: int) =
    if n != 0:
      result &= s & spaces(n+2-s.len)
  for t in x:
    pad t.name, nameLen
    pad t.value, valueLen
    pad t.comment, commLen
    pad t.file, fileLen
    #pad linePrefix & $t.line, lineLen
    #result &= "\n"
    result &= linePrefix & $t.line & "\n"

proc fmtParams(x: seq[ParamObj], thisfile: string): string =
  result = ""
  # find all files with params
  var fileList = newSeq[string](0)
  for i in 0..<x.len:
    let f = x[i].file
    if f.notIn fileList:
      fileList.add f
  # move thisfile last
  var thisi = 0
  while thisi < fileList.len and fileList[thisi] != thisfile: inc thisi
  while thisi+1 < fileList.len:
    swap fileList[thisi], fileList[thisi+1]
    inc thisi
  if thisi < fileList.len: fileList[thisi] = thisfile
  # output one file at a time
  for f in fileList:
    var paramList = newSeq[ParamObj](0)
    for i in 0..<x.len:
      var t = x[i]
      if t.file == f:
        t.comment = "#" & t.comment
        t.file = ""
        paramList.add t
    paramList.sort proc(x, y: ParamObj):int = cmp(x.line, y.line)
    result &= "\n# Params from file " & f & "\n"
    result &= format paramList

proc saveParams(x: seq[ParamObj], fn: string, loc: string, thisfile: string) =
  mixin echo
  var o = "# QEX compile-time generated parameter file\n"
  o &= "# generated from: " & loc & "\n"
  o &= "# use -" & loadParamsCommand & ":<filename> to read in\n"
  o &= "# it can also be converted to command line arguments with:\n"
  o &= r"""#  awk '/^[a-zA-Z]/{split($2,a,"#");if(length(a[1])>0)print("-"$1":"a[1])}' <file>""" & "\n"
  o &= fmtParams(x, thisfile)
  echo "Writing parameter file to ", fn
  writeFile(fn, o)

proc addParam(s,r: string, c: string = "", ii: IiType) =
  paramInfo.set(s, r, c, ii.filename, ii.line)

let spcName = "                "
let spcValue = "                                   #"
let spcLoc = "                                            ## "
proc addComment(s,c,spc:string):string =
  result = s
  if c.len>0:
    let
      m = min(s.len, spc.len-8)
    result &= spc[m..^1] & c
proc addComment(s,c:string):string = addComment(s,c,spcValue)

proc echoParamsX*(warnUnknown=false):string =
  result = "Params:\n"
  for i in 0..<paramInfo.len:
    let t = paramInfo[i]
    result &= "  " & addComment(t.name & ": " & t.value, t.comment) & "\n"
  for p in params.keys:
    let i = paramInfoCTRT.find p
    if i < 0:
      let j = paramInfo.find p
      if j < 0:
        result &= " Unknown argument: '" & $p & "'\n"
  result.removeSuffix '\n'
template echoParams*(warnUnknown=false) =
  mixin echo
  echo echoParamsX(warnUnknown)

proc paramHelp*(p:string = ""):string =
  #echo paramInfo
  #echo paramInfoCTRT
  result = "Usage:\n  " & getAppFileName()
  var t: ParamObj
  var i = paramInfo.find p
  if i >= 0:
    t = paramInfo[i]
  else:
    i = paramInfoCTRT.find p
    if i >= 0:
      t = paramInfoCTRT[i]
  if i >= 0:
    result &= addComment(" -" & p & ":" & t.value & " (current value)", t.comment)
  else:
    result &= " -OPTION:VALUE ...\nAvailable OPTIONs and current VALUEs:"
    template p(t) =
      let nm = t.name
      var s = nm & spcName[min(spcName.len-1,nm.len)..^1] & " : " & t.value
      let f = t.file.splitfile[1]
      let c = f & "," & $t.line
      s = s.addComment(c,spcValue)
      s = s.addComment(t.comment,spcLoc)
      result &= "\n  " & s
    var paramnames = newSeq[string](0)
    for t in paramInfo:
      paramnames.add t.name
      p(t)
    for t in paramInfoCTRT:
      if t.name notin paramnames:  # FIXME: should check location too
        p(t)

template cnvnone(x:typed):untyped = x
template makeTypeParam(name,typ,deflt,cnvrt: untyped): untyped {.dirty.} =
  proc `name X`*(s: string, d: typ, c: string, ii: tuple): typ =
    setup()
    result = d
    if params.hasKey(s):
      result = cnvrt(getParamValue(s))
    addParam(s, $result, c, ii)
    #cho ii
  template name*(s: string, d=deflt, c="", index= -1): typ =
    addParamCT(s, d, c, instantiationInfo(index, fullPaths=true))
    `name X`(s, d, c, instantiationInfo(index, fullPaths=true))

makeTypeParam(intParam, int, 0, parseInt)
makeTypeParam(floatParam, float, 0.0, parseFloat)
makeTypeParam(strParam, string, "", cnvnone)
template stringParam*(x,y: untyped, c="", index= -1): untyped = strParam(x,y,c,index)

proc boolParamX*(s: string, d: bool, c: string, ii: tuple): bool =
  setup()
  result = d
  if params.hasKey(s):
    let val = tolowerAscii(getParamValue(s))
    result = case val
             of "","t","true","yes","y","on": true
             else: false
  addParam(s, $result, c, ii)
template boolParam*(s: string, d = false, c = "", index = -1): bool =
  addParamCT(s, d, c, instantiationInfo(index, fullPaths=true))
  boolParamX(s, d, c, instantiationInfo(index, fullPaths=true))

proc intSeqParamX*(s: string, d: seq[int], c: string, ii: tuple): seq[int] =
  setup()
  result = d
  if params.hasKey(s):
    result.setLen(0)
    for c in split(getParamValue(s), ','):
      if c.len > 0:
        result.add parseInt(c)
  addParam(s, join(result,","), c, ii)
#template intSeqParam*(s: string, d: seq[int] = @[], c = ""): seq[int] =
template intSeqParam*(s: string, d = newSeq[int](), c = ""): seq[int] =
  addParamCT(s, d, c, instantiationInfo(fullPaths=true))
  intSeqParamX(s, d, c, instantiationInfo(fullPaths=true))

proc floatSeqParamX*(s: string, d: seq[float], c: string, ii: tuple): seq[float] =
  setup()
  result = d
  if params.hasKey(s):
    result.setLen(0)
    for c in split(getParamValue(s), ','):
      if c.len > 0:
        result.add parseFloat(c)
  addParam(s, join(result,","), c, ii)
#template floatSeqParam*(s: string, d: seq[float] = @[], c = ""): seq[float] =
template floatSeqParam*(s: string, d = newSeq[float](), c = ""): seq[float] =
  addParamCT(s, d, c, instantiationInfo(fullPaths=true))
  floatSeqParamX(s, d, c, instantiationInfo(fullPaths=true))

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

var helpValue = false
template installHelpParam*(p: static string = "h", index= -2) =
  helpValue = boolParam(p, false, "Print the help message", index)
template processHelpParam*() =
  mixin echo
  if helpValue:
    echo paramHelp()
    qexExit()

# gets processed by init, no need for separate process call
template installLoadParams*(p: static string = loadParamsCommand, index= -2) =
  discard strParam(p, "", "Load params from file", index)
  static:
    if p != "": loadParamsCommand = p

var saveValue = ""
template installSaveParams*(p: static string = "saveParams", index= -2) =
  saveValue = strParam(p, "", "Save params to file", index)
template processSaveParams*(index = -1) =
  if saveValue != "":
    let ii = instantiationInfo(index, fullPaths=true)
    saveParams(paramInfo, saveValue, ii.filename & ":" & $ii.line, ii.filename)

var makeParamInfoCTRT: proc()
template processMakeParamInfoCTRT* =
  if not isNil makeParamInfoCTRT: makeParamInfoCTRT()

template installStandardParams*(idx= -3) =
  installSaveParams(index=idx)
  installLoadParams(index=idx)
  installHelpParam(index=idx)
  #processMakeParamInfoCTRT()

# write param file at compile time
macro writeParamFileX*(filename: static string, ii: static IiType) =
  let thisfile = ii.filename
  var fn = filename
  if fn == "":
    var base = thisfile.lastPathPart
    base.removeSuffix(".nim")
    fn = base & ".qexin.ctsample"
  saveParams(paramInfoCT, fn, thisfile & ":" & $ii.line, thisfile)
# convert paramInfoCT to RT
macro makeParamInfoCTRTX*(): auto =
  result = newStmtList()
  for t in paramInfoCT:
    let n = newLit t.name
    let v = newLit t.value
    let c = newLit t.comment
    let f = newLit t.file
    let l = newLit t.line
    result.add quote do:
      paramInfoCTRT.add ParamObj(name:`n`,value:`v`,comment:`c`,file:`f`,line:`l`)
template finalizeParams* =
  proc makeParamInfoCTRTImpl(): bool =
    makeParamInfoCTRTX()
    true
  #makeParamInfoCTRT = makeParamInfoCTRTImpl
  let dum {.global.} = makeParamInfoCTRTImpl()
template writeParamFile*(filename: static string = "") =
  writeParamFileX(filename, instantiationInfo(fullPaths=true))

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

when isMainModule:
  import qex
  qexInit()
  installHelpParam("h")
  processHelpParam()

  letParam:
    bf = false
    bt = true
    bx = if true: true else: false
    i0 = 0
    i1 = 1
    ix = if true: 2 else: 3
    f0 = 0.0
    f1 = 1.0
    fx = if true: 2.0 else: 3.0
    s0 = "foo0"
    s1 = "foo1"
    sx = if true: "foo2" else: "foo3"
    ia0 = @[0,0,0,0]
    ia1 = @[1,1,1,1]
    iax = if true: @[2,2,2,2] else: @[3,3,3,3]
    fa0 = @[0.0,0,0,0]
    fa1 = @[1.0,1,1,1]
    fax = if true: @[2.0,2,2,2] else: @[3.0,3,3,3]

  echoParams()

  defaultSetup()
  proc paramTest* =
    var i = intParam("i", 0, "Test intParam")
    var f = floatParam("f", 0.0, "Test floatParam")
    echo i, f
  paramTest()

  writeParamFile("")
  qexFinalize()
  finalizeParams()
