import threading
export threading
import comms/commsEcho, stdUtils, base/[basicOps,params]
import os, strutils, sequtils, std/monotimes, std/tables, std/algorithm, strformat
export monotimes
getOptimPragmas()

const noTicToc {.boolDefine.} = false

type PerfInfo* = object
  count*: int
  flops*: float
  secs*: float
template clear*(pi: var PerfInfo) =
  pi.count = 0
  pi.flops = 0.0
  pi.secs = 0.0
template `+=`*(r: var PerfInfo, x: PerfInfo) =
  r.count += x.count
  r.flops += x.flops
  r.secs += x.secs
proc `$`*(pi: PerfInfo): string =
  result = system.`$`(pi)
  result &= &" {1e-9*pi.flops/pi.secs:.2f} Gflops"

type TicType* = distinct int64
template getTics*: TicType = TicType(getMonoTime().ticks)
template nsec*(t:TicType):int64 = int64(t)
template seconds*(t:TicType):float = 1e-9 * float(t)
template ticDiffSecs*(x,y: TicType): float = 1e-9 * float(x.int64 - y.int64)
template `-`*(x,y: TicType): TicType = TicType(x.int64 - y.int64)

var
  ## Drop children timers if the proportion of their overhead is larger than this.
  DropWasteTimerRatio* = floatParam("dropRatio", 0.05)
  VerboseTimer* = false  ## If true print out all the timers during execution.

##[

Each tic() starts a local timer.  Each toc() records the time difference
of the current time with the one in the local timer visible in the
scope, and then update the timer with the current time.

A tree structure saves the time durations, and represents the
hierarchical runtimes of each code segments between tic/toc's.

tic() injects symbols:
  localTimer
  prevRTI
  localTimerStart
  localTic

]##

type
  II = ptr typeof(instantiationInfo())
  RTInfo = distinct int
  RTInfoObj = object
    nsec: int64
    flops: float
    overhead: int64
    childrenOverhead: int64
    count: uint32
    tic, prev, curr: CodePoint
    children: RTInfoObjList
  CodePoint = distinct int32
  CodePointObj = object
    toDropTimer: bool
    #nsec: int64
    #overhead: int64
    count: uint32
    dropcount: uint32
    name: CStr
    loc: II
  SString = static[string] | string
  List[T] = object  # No GCed type allowed
    len,cap:int32
    data:ptr UncheckedArray[T]
  RTInfoObjList = distinct List[RTInfoObj]
  CStr = object
    p, s: int32

var
  rtiListLength:int32 = 0
  rtiListLengthMax:int32 = 0
template listChangeLen[T](n:int32) =
  when T is RTInfoObj:
    rtiListLength += n
    if n>0 and rtiListLength>rtiListLengthMax:
      rtiListLengthMax = rtiListLength

proc newList[T](len:int32 = 0):List[T] {.noinit.} =
  var cap = if len == 0: 0.int32 else: 1.int32
  while cap < len: cap *= 2
  result.len = len
  result.cap = cap
  if cap > 0:
    result.data = cast[ptr UncheckedArray[T]](allocShared(sizeof(T)*cap))
    listChangeLen[T](cap)
  else:
    result.data = nil
proc newListOfCap[T](cap:int32):List[T] {.noinit.} =
  result.len = 0
  result.cap = cap
  if cap > 0:
    result.data = cast[ptr UncheckedArray[T]](allocShared(sizeof(T)*cap))
    listChangeLen[T](cap)
  else:
    result.data = nil
proc len[T](ls:List[T]):int32 = ls.len
proc setLen[T](ls:var List[T], len:int32) =
  if len > ls.cap:
    let cap0 = ls.cap
    var cap = cap0
    if cap0 == 0:
      cap = 1
      while cap < len: cap *= 2
      ls.data = cast[ptr UncheckedArray[T]](allocShared(sizeof(T)*cap.int))
    else:
      while cap < len: cap *= 2
      ls.data = cast[ptr UncheckedArray[T]](reallocShared(ls.data, sizeof(T)*cap.int))
    ls.cap = cap
    listChangeLen[T](int32(cap-cap0))
  ls.len = len
proc free[T](ls:var List[T]) =
  if ls.cap > 0:
    deallocShared(ls.data)
  listChangeLen[T](-ls.cap)
  ls.len = 0
  ls.cap = 0
  ls.data = nil
template `[]`[T](ls:List[T], n:int):untyped = ls.data[n]
proc add[T](ls:var List[T], x:T) =
  let n = ls.len
  ls.setLen(n+1)
  ls.data[n] = x
iterator items[T](ls:List[T]):T =
  for i in 0..<ls.len:
    yield ls[i]

proc setLen(ls:var RTInfoObjList, len:int32) {.borrow.}
proc free(ls:var RTInfoObjList) {.borrow.}
#proc add(ls:var RTInfoObjList, x:RTInfoObj) {.borrow.}
proc add(ls:var RTInfoObjList, x:RTInfoObj) = add(List[RTInfoObj](ls), x)
template len(ls:RTInfoObjList):int32 = List[RTInfoObj](ls).len
template `[]`(ls:RTInfoObjList, n:int32):untyped = List[RTInfoObj](ls)[n]
iterator mitems(ls:RTInfoObjList):var RTInfoObj =
  for i in 0..<ls.len:
    yield ls[i]

func isNil(x:RTInfo):bool = x.int<0
func isNil(x:CodePoint):bool = x.int<0

const defaultCStrPoolCap {.intDefine.} = 512
type CStrAtom = array[16,char]
var cstrpool = newListOfCap[CStrAtom](defaultCStrPoolCap)

proc len(s:CStr):int = int(s.s)
proc newCStr(t:string):CStr =
  const a = int32(sizeof(CStrAtom))
  let p = cstrpool.len
  let s = int32(t.len)
  let n = (s+a-1) div a
  cstrpool.setlen(p+n)
  var k = 0
  var j = 0
  for i in 0..<s:
    cstrpool[p+k][j] = t[i]
    inc j
    if j == a:
      j = 0
      inc k
  CStr(p:p, s:s)
proc equal(s:CStr, t:string):bool =
  let len = s.s
  if len != t.len:
    return false
  var k = 0
  var j = 0
  for i in 0..<len:
    if t[i] != cstrpool[s.p+k][j]:
      return false
    inc j
    if j == sizeof(CStrAtom):
      j = 0
      inc k
  return true
proc `$`(s:CStr):string =
  result = newString(s.s)
  var k = 0
  var j = 0
  for i in 0..<s.s:
    result[i] = cstrpool[s.p+k][j]
    inc j
    if j == sizeof(CStrAtom):
      j = 0
      inc k

func `$`(x:II):string =
  x.filename & ":" & $x.line & ":" & $x.column

func `$`(x:List[RTInfo]):string =
  result = "RTInfo[ "
  for i in x: result &= $i.int & " "
  result &= "]"
func `$`(x:List[CodePoint]):string =
  result = "CodePoint[ "
  for i in x: result &= $i.int & " "
  result &= "]"
func `$`(x:RTInfo):string = "RTInfo(" & $x.int & ")"
func `$`(x:CodePoint):string = "CodePonit(" & $x.int & ")"

func `$`(x:List[RTInfoObj]):string  # declaration
func `$`(x:RTInfoObjList):string {.borrow.}
func `$`(x:RTInfoObj):string =
  if x.prev.isNil:
    result = "RTInfoObj:TIC CodePoint(" & $x.curr.int & ")"
  else:
    result = "RTInfoObj(" & $x.nsec & " + " & $x.overhead & " + " & $x.childrenOverhead & " ns "
    result &= $x.flops & " f " & $x.count
    result &= " CodePoint(" & $x.tic.int & ":" & $x.prev.int & ".." & $x.curr.int & ") "
    result &= indent($x.children, 3)[3..^1] & ")"
func `$`(x:List[RTInfoObj]):string =
  result = "RTInfoObj[\n"
  for i in 0..<x.len:
    result &= "   " & $i & ":" & $x[i] & "\n"
  result[^1] = ']'

proc `$`(x:seq[CodePointObj]):string =
  result = "CodePointObj[\n"
  for i in 0..<x.len:
    result &= "   " & $i & ":" & $x[i] & "\n"
  result[^1] = ']'

proc `==`(x,y:RTInfo):bool = x.int==y.int
proc `==`(x,y:CodePoint):bool = x.int==y.int

const
  defaultRTICap {.intDefine.} = 512
  defaultLocalRTICap {.intDefine.} = 0

var
  rtiStack = newListOfCap[RTInfoObj](defaultRTICap)
  #cpHeap = newSeqOfCap[CodePointObj](defaultRTICap)
  cpHeap = newListOfCap[CodePointObj](defaultRTICap)
  frozenTimers = false

proc timersFrozen*:bool = frozenTimers
proc thawTimers* = frozenTimers = false
proc freezeTimers* = frozenTimers = true

proc newCodePoint(ii:II, s:SString = ""):CodePoint =
  let n = cpHeap.len
  cpHeap.setlen(n+1)
  cpHeap[n].toDropTimer = false
  cpHeap[n].count = 0
  cpHeap[n].dropcount = 0
  cpHeap[n].name = newCStr(s)
  cpHeap[n].loc = ii
  CodePoint(n)

template name(x:CodePoint):CStr = cpHeap[x.int].name
template loc(x:CodePoint):II = cpHeap[x.int].loc
template count(x:CodePoint):uint32 = cpHeap[x.int].count
template dropcount(x:CodePoint):uint32 = cpHeap[x.int].dropcount
template toDropTimer(x:CodePoint):bool = cpHeap[x.int].toDropTimer
template dropTimer(x:CodePoint) =
  toDropTimer(x) = true
  #echo "dropTimer: ", cpHeap[x.int].loc, " ", cpHeap[x.int].name

template overhead(x:RTInfoObj):untyped = x.overhead
template childrenOverhead(x:RTInfoObj):untyped = x.childrenOverhead
template istic(x:RTInfoObj):bool = x.prev.isNil
template isnottic(x:RTInfoObj):bool = not x.prev.isNil
template toDropTimer(x:RTInfoObj):bool = toDropTimer(x.curr)
template dropTimer(x:RTInfoObj) = dropTimer(x.curr)
proc dropTimerRecursive(x:RTInfoObj) =
  #echo "dropTimerRecursive: ", x.prev.loc, " ", x.curr.loc, " ", x.curr.name
  dropTimer(x)
  for i in 0..<x.children.len:
    dropTimer(x.children[i].tic)
    dropTimer(x.children[i].prev)
    dropTimer(x.children[i].curr)
    dropTimerRecursive(x.children[i])
proc dropTimerChildren(x:RTInfoObj) =
  #echo "dropTimerRecursive: ", x.prev.loc, " ", x.curr.loc, " ", x.curr.name
  #dropTimer(x)
  for i in 0..<x.children.len:
    #dropTimer(x.children[i].tic)
    #echo "dropTimer: ", x.children[i].tic.name, " ", x.children[i].curr.name, " ", x.children[i].curr.loc
    dropTimer(x.children[i].prev)
    #dropTimer(x.children[i].curr)
    dropTimerChildren(x.children[i])

template toDropTimer(x:RTInfo):bool = toDropTimer(rtiStack[x.int])
template dropTimer(x:RTInfo) = dropTimer(rtiStack[x.int])
#template dropTimerRecursive(x:RTInfo) = dropTimerRecursive(rtiStack[x.int])
template dropTimerChildren(x:RTInfo) = dropTimerChildren(rtiStack[x.int])

template identical(x,y:RTInfoObj):bool =
  x.tic == y.tic and x.prev == y.prev and x.curr == y.curr

proc combineList(acc:var RTInfoObjList, xs:RTInfoObjList)  # forward declaration
template combine(acc:var RTInfoObjList, x:var RTInfoObj) =
  if isnottic(x):
    let
      nsec = x.nsec
      flops = x.flops
      count = x.count
      overhead = x.overhead
      childrenOverhead = x.childrenOverhead
      tic = x.tic
      prev = x.prev
      curr = x.curr
    var children = x.children
    var ci:int32 = -1
    for i in 0..<acc.len:
      if identical(x,acc[i]):
        ci = i
        break
    if ci < 0:
      # new
      let n = acc.len
      acc.setlen(n+1)
      acc[n].nsec = nsec
      acc[n].flops = flops
      acc[n].count = count
      acc[n].overhead = overhead
      acc[n].childrenOverhead = childrenOverhead
      acc[n].tic = tic
      acc[n].prev = prev
      acc[n].curr = curr
      acc[n].children = children
    else:
      # add to ci
      acc[ci].nsec += nsec
      acc[ci].flops += flops
      acc[ci].count += count
      acc[ci].overhead += overhead
      acc[ci].childrenOverhead += childrenOverhead
      combineList(acc[ci].children, children)
      children.free

proc combineList(acc:var RTInfoObjList, xs:RTInfoObjList) =
  for i in 0..<xs.len:
    combine(acc, xs[i])

proc record(tic:RTInfo, prev:RTInfo, curr:CodePoint, t:TicType, f:float):RTInfo =
  var
    children = RTInfoObjList(newListOfCap[RTInfoObj](defaultLocalRTICap))
    oh:int64 = 0
  for i in prev.int+1..<rtiStack.len:
    # need to consolidate the stack
    combine(children, rtiStack[i])
    oh += rtiStack[i].overhead
    oh += rtiStack[i].childrenOverhead

  # Now for current RTInfo
  let n = prev.int32+1
  rtiStack.setlen(n+1)
  rtiStack[n].nsec = nsec(t)
  rtiStack[n].flops = f
  rtiStack[n].count = 1
  # Assign overhead later.
  rtiStack[n].childrenOverhead = oh
  rtiStack[n].tic = rtiStack[tic.int].curr
  rtiStack[n].prev = rtiStack[prev.int].curr
  rtiStack[n].curr = curr
  rtiStack[n].children = children
  RTInfo(n)

proc recordTic(this:CodePoint):RTInfo =
  # only needs prev and curr for tic
  # childrenOverhead for others
  # Assign overhead later.
  let n = rtiStack.len
  rtiStack.setlen(n+1)
  rtiStack[n].prev = CodePoint(-1)
  rtiStack[n].curr = this
  rtiStack[n].childrenOverhead = 0
  RTInfo(n)

proc echoTic*(s: string, ii: II) =
  echo "tic ",s,ii

proc ticSet(localTimer:var TicType, prevRTI:var RTInfo, restartTimer:var bool,
            thisCode:var CodePoint, s:SString,ii:II, localCodePtr:auto) {.alwaysInline.} =
  #echo "#### begin tic ",ii
  if unlikely VerboseTimer: echoTic(s,ii)
  if not timersFrozen():
    let theTime = getTics()
    when localCodePtr isnot bool:
      for c in items(localCodePtr[]):
        if cpHeap[c.int].name.equal(s):
          thisCode = c
          break
    if thisCode.isNil:
      thisCode = newCodePoint(ii, s)
      when localCodePtr isnot bool:
        localCodePtr[].add thisCode
    prevRTI = recordTic(thisCode)
    if toDropTimer(thisCode):
      freezeTimers()
      restartTimer = true
    localTimer = getTics()
    rtiStack[prevRTI.int].overhead = nsec(localTimer-theTime)
  else:
    localTimer = getTics()
  #echo "#### end tic ",ii

template ticI(n = -1; s:SString = "") =
  bind items
  const
    cname = compiles(static[string](s))
  var ii {.global.} = instantiationInfo(n)
  var
    localTimer {.inject.}: TicType
    prevRTI {.inject.} = RTInfo(-1)
    restartTimer {.inject.} = false
  when cname:
    var thisCode {.global.} = CodePoint(-1)
  else:
    var
      localCode {.global.} = newList[CodePoint]()
      thisCode = CodePoint(-1)
  if threadNum==0:
    when false:
      #echo "#### begin tic ",ii
      if unlikely VerboseTimer: echoTic(s,ii)
      if not timersFrozen():
        let theTime = getTics()
        when not cname:
          for c in items(localCode):
            if cpHeap[c.int].name.equal(s):
              thisCode = c
              break
        if thisCode.isNil:
          thisCode = newCodePoint(ii.addr, s)
          when not cname:
            localCode.add thisCode
        prevRTI = recordTic(thisCode)
        if toDropTimer(thisCode):
          freezeTimers()
          restartTimer = true
        localTimer = getTics()
        rtiStack[prevRTI.int].overhead = nsec(localTimer-theTime)
      #echo "#### end tic ",ii
    else:
      when cname:
        ticSet(localTimer,prevRTI,restartTimer,thisCode,s,ii.addr,false)
      else:
        ticSet(localTimer,prevRTI,restartTimer,thisCode,s,ii.addr,addr localCode)
  let
    localTimerStart {.inject, used.} = localTimer
    localTic {.inject, used.} = prevRTI

when noTicToc:
  template tic0 =
    var localTimerStart {.inject,used.} = getTics()
  template tic*() = tic0
  template tic*(n: int) = tic0
  template tic*(s: SString) = tic0
  template tic*(n: int; s: SString) = tic0
else:
  template tic*(n = -1; s:SString = "") = ticI(n-1,s)
  template tic*(s:SString = "") = ticI(-2,s)

proc echoToc*(s: string, ii: II) =
  echo "toc ",s,ii

proc tocSet(localTimer:var TicType, prevRTI:var RTInfo, restartTimer:var bool,
            thisCode:var CodePoint, f:SomeNUmber, s:SString, ii:II,
            localTic:RTInfo, localCodePtr:auto) {.alwaysInline.} =
  #echo "==== begin toc ",s," ",ii
  #echo "     rtiStack: ",indent($rtiStack,5)
  #echo "     cpHeap: ",indent($cpHeap,5)
  if unlikely VerboseTimer: echoToc(s,ii)
  if prevRTI.int32 >= 0:
    if restartTimer:
      thawTimers()
      restartTimer = false
    if not timersFrozen():
      let theTime = getTics()
      when localCodePtr isnot bool:
        for c in items(localCodePtr[]):
          if cpHeap[c.int].name.equal(s):
            thisCode = c
            break
      if thisCode.isNil:
        thisCode = newCodePoint(ii, s)
        when localCodePtr isnot bool:
          localCodePtr[].add thisCode
      let
        ns = theTime-localTimer
        thisRTI = record(localTic, prevRTI, thisCode, ns, float(f))
      var oh = rtiStack[thisRTI.int].childrenOverhead
      let c = rtiStack[thisRTI.int].children
      for i in 0..<c.len:
        if toDropTimer(c[i].prev):
          oh -= c[i].childrenOverhead
      inc thisCode.count
      if oh.float > ns.float*DropWasteTimerRatio:
      #if not toDropTimer(prevRTI) and oh.float > ns.float*DropWasteTimerRatio:
        inc thisCode.dropcount
        #if ii.filename != "scg.nim":
        #  echo "drop timer: ", oh.float, "/", ns.float, "=", oh.float / ns.float
        #  echo "  ", prevRTI.int, " ", thisRTI.int, " ", ii, " ", s
        # Signal stop if the overhead is too large.
        if thisCode.dropcount > 10 and thisCode.dropcount*10 > thisCode.count:
          #echo "dropTimer: ", rtiStack[thisRTI.int].tic.name, " ", thisCode.name, " ", thisCode.loc
          dropTimer(prevRTI)
          dropTimerChildren(thisRTI)
          #dropTimer(thisCode)
      if toDropTimer(thisCode):
        freezeTimers()
        restartTimer = true
      localTimer = getTics()
      rtiStack[thisRTI.int].overhead = nsec(localTimer-theTime)
      prevRTI = thisRTI
    else:
      localTimer = getTics()
  #echo "==== end toc ",s," ",ii

template tocI(f: SomeNumber; s:SString = ""; n = -1) =
  bind items
  const
    cname = compiles(static[string](s))
  var ii {.global.} = instantiationInfo(n)
  when cname:
    var thisCode {.global.} = CodePoint(-1)
  else:
    var
      localCode {.global.} = newList[CodePoint]()
      thisCode = CodePoint(-1)
  if threadNum==0:
    when cname:
      tocSet(localTimer,prevRTI,restartTimer,thisCode,f,s,ii.addr,localTic,false)
    else:
      tocSet(localTimer,prevRTI,restartTimer,thisCode,f,s,ii.addr,localTic,addr localCode)

when noTicToc:
  template toc*() = discard
  template toc*(n:int) = discard
  template toc*(s:SString) = discard
  template toc*(s:SString, n:int) = discard
  template toc*(n:int, flops:SomeNumber) = discard
  template toc*(s:SString = "", n = -1, flops:SomeNumber) = discard
else:
  template toc*(s:SString = ""; n = -1; flops:SomeNumber) = tocI(flops, s, n-1)
  template toc*(n = -1; flops:SomeNumber) = tocI(flops, "", n-1)
  template toc*(s:SString; n:int) = tocI(0, s, n-1)
  template toc*(s:SString) = tocI(0, s, -2)
  template toc*(n:int) = tocI(0, "", n-1)
  template toc*() = tocI(0, "", -2)

#when noTicToc:
#  template getElapsedTime*: float = 0.0
#else:
template getElapsedTime*: float = ticDiffSecs(getTics(), localTimerStart)

proc reset(x:var RTInfoObj) =
  x.nsec = 0
  x.flops = 0.0
  x.count = 0
  x.overhead = 0
  x.childrenOverhead = 0
  toDropTimer(x.curr) = false
  if x.isnottic():  # The children list is not initialized for tics.
    for c in mitems(x.children):
      reset c

template resetTimers* =
  ## Reset timers in the local scope, starting from the local tic, and below.
  ## Do nothing if localTic is uninitialized, as is the case with frozenTimers.
  when declared(localTic):
    let go = localTic.int32 >= 0
  else:
    const go = true
  if threadNum==0 and go:
    when declared(localTic):
      let p = localTic.int32
    else:
      let p:int32 = 0
    for j in p..<len(rtiStack):
      reset rtiStack[j]

template aggregateTimers* =
  ## Aggregate timers in the local scope, starting from the local tic, and below.
  ## It scrambles the local timer branches, and the direct children of the current timer.
  ## Do nothing if localTic is uninitialized, as is the case with frozenTimers.
  when declared(localTic):
    let go = localTic.int32 >= 0
  else:
    const go = true
  if threadNum==0 and go:
    when declared(localTic):
      let
        p = localTic.int32+1
        theTime = getTics()
    else:
      let p:int32 = 0
    if p<rtiStack.len-2:
      var rs = RTInfoObjList(newListOfCap[RTInfoObj](defaultLocalRTICap))
      when declared(localTic):
        var oh:int64 = 0
      for i in p..<rtiStack.len:
        combine(rs, rtiStack[i])
        when declared(localTic):
          if istic(rtiStack[i]):
            # Combine ignores the tics, so we do it here.
            # The overhead counts are lost if we don't have a localTic.
            oh += overhead(rtiStack[i])
            oh += childrenOverhead(rtiStack[i])
      when declared(localTic):
        childrenOverhead(rtiStack[localTic.int]) += oh
      let nl = rs.len
      if nl!=rtiStack.len-p:
        when declared(prevRTI):
          if prevRTI.int > p:  # It's fine if prevRTI.int==p.
            let pr = prevRTI.int32-p
            if not(pr<nl and identical(rs[pr], rtiStack[prevRTI.int32])):
              # Need to make sure prevRTI is still correct.
              # Unlike `record`, here we aggregate from local tic instead of prevRTI.
              # `combine` considers identical timers if they have the save tic/curr/prev,
              # and add local and children tocs together.
              # In general, for a rtiStack:
              # ... localTic [... A B C ... A ...] [X] [... B D E ...] X(=prevRTI) [... C D F ...]
              # we get
              # ... localTic [... A B C ...] X(=prevRTI) [... D E ... ] [... F ...]
              for i in 1..<nl:
                if identical(rs[i], rtiStack[prevRTI.int]):
                  prevRTI = RTInfo(p+i)
                  break
        copyMem(rtiStack[p].addr, rs[0].addr, nl*sizeof(RTInfoObj))
        rtiStack.setlen(p+nl)
      free(rs)
    when declared(localTic):
      childrenOverhead(rtiStack[localTic.int]) += nsec(getTics()-theTime)

type
  Tstr = tuple
    label: string
    stats: string
func markMissing(p:bool,str:string):string =
  if p: "[" & str & "]"
  else: str
template ppT(ts: RTInfoObjList, prefix = "-", total = 0'i64, overhead = 0'i64,
            count = 0'u32, initIx = 0, showAbove = 0.0, showDropped = true): seq[Tstr] =
  ppT(List[RTInfoObj](ts), prefix, total, overhead, count, initIx, showAbove, showDropped)
proc ppT(ts: List[RTInfoObj], prefix = "-", total = 0'i64, overhead = 0'i64,
        count = 0'u32, initIx = 0, showAbove = 0.0, showDropped = true): seq[Tstr] =
  var
    sub:int64 = 0
    subo:int64 = 0
    pre = prefix
  for j in initIx..<ts.len:
    let nc = ts[j].count
    if ts[j].istic or nc==0: continue
    if j==ts.len-1 and prefix.len>1:
      pre[^3] = '`'
    let
      f0 = splitFile(ts[j].prev.loc.filename)[1]
      l0 = ts[j].prev.loc.line
      f = splitFile(ts[j].curr.loc.filename)[1]
      l = ts[j].curr.loc.line
      drop = ts[j].prev.toDropTimer
      coh = ts[j].childrenOverhead
      soh = ts[j].overhead
      nsec = ts[j].nsec
      ns = nsec - coh
      oh = soh + coh
      nf = ts[j].flops
      tn = ts[j].tic.name
      pn = ts[j].prev.name
      st = ns div 1000
      ot = oh div 1000
      sc = ns div nc.int64
      oc = oh div nc.int64
      mf = nf*1e3 / ns.float
      small = total!=0 and ns.float/total.float<showAbove
      noexpand = drop or (small and ts[j].children.len>0)
      loc = pre & markMissing(noexpand, f0 & "(" & $l0 & "-" & (if f==f0:"" else:f) & $l & ")")
      nm = pre & markMissing(noexpand, (if tn.len==0:"" else: $tn & ":") & (if pn.len==0:"" else: $pn & "-") & $ts[j].curr.name)
    if total!=0:
      let
        cent = 1e2 * ns.float / total.float
        ohcent = 1e2 * oh.float / total.float
      result.add (loc, cent|(6,-1) & ohcent|(6,-1) & st|12 & ot|8 & " /" & nc|7 & " =" & sc|12 & oc|8 & int(mf)|8 & " " & nm)
    else:
      result.add (loc, ""|6 & ""|6 & st|12 & ot|8 & " /" & nc|7 & " =" & sc|12 & oc|8 & int(mf)|8 & " " & nm)
    if ts[j].children.len>0 and
        (not small) and
        ((not drop) or showDropped):
      let newprefix =
        if prefix.len==1: "|--"
        elif j<ts.len-1: prefix[0..^4] & "| |--"
        else: prefix[0..^4] & "  |--"
      result.add ppT(ts[j].children, newprefix, nsec+soh, soh, nc, 0, showAbove, showDropped)
    sub += ns
    subo += oh
  if total!=0 and count!=0:
    let
      ns = total-overhead-sub-subo
      st = ns div 1000
      ot = overhead div 1000
      sc = ns div count.int64
      oc = overhead div count.int64
      cent = 1e2 * ns.float / total.float
      ohcent = 1e2 * overhead.float / total.float
    result = @[(prefix & "#me",
      cent|(6,-1) & ohcent|(6,-1) & st|12 & ot|8 & " /" & count|7 & " =" & sc|12 & oc|8 & ""|8 & " " & prefix & "#me")] & result

proc totalTime(ts:List[RTInfoObj], initIx = 0):int64 =
  for j in initIx..<ts.len:
    if ts[j].istic:
      result += ts[j].overhead
    else:
      result += ts[j].nsec + ts[j].overhead

proc totalOverhead(ts:List[RTInfoObj], initIx = 0):int64 =
  for j in initIx..<ts.len:
    result += ts[j].overhead + ts[j].childrenOverhead

template echoTimers*(expandAbove = 0.0, expandDropped = true, aggregate = true) =
  ## Echo timers in the local scope, starting from the local tic, and below.
  ## Expand children timers if the proportion of the current work time is more than expandAbove.
  if threadNum==0:
    const width = 104
    when declared(localTic):
      let p = if localTic.int<0: 0 else: localTic.int
    else:
      let p = 0
    if aggregate: aggregateTimers()
    let
      pp = ppT(rtiStack, initIx = p, showAbove = expandAbove, showDropped = expandDropped)
      tt = totalTime(rtiStack, initIx = p)
      oh = totalOverhead(rtiStack, initIx = p)
    var n = 24
    for (s,_) in pp:
      if n<s.len: n = s.len
    inc n
    let notshowing = if expandAbove>0.0: ", not expanding contributions less than " & $(1e2*expandAbove) & " %" else:""
    echo "Timer total ",(tt.float*1e-6)|(0,-3)," ms, overhead ",(oh.float*1e-6)|(0,-3)," ms ~ ",(1e2*oh.float/tt.float)|(0,-1)," %, runtime info ",rtiListLength*sizeof(RTInfoObj)," B, max ",rtiListLengthMax*sizeof(RTInfoObj)," B, string ",cstrpool.len*sizeof(cstrpool[0])," B",notshowing
    echo '='.repeat(width)
    echo "file(lines)"|(-n), "%"|6, "OH%"|6, "microsecs"|12, "OH"|8, "count"|9, "ns/count"|14, "OH/c"|8, "mf"|8, " label"
    echo '='.repeat(width)
    for (s,t) in pp:
      echo s|(-n,'.'), t
    echo '='.repeat(width)

proc echoTimersRaw* =
  if threadNum==0:
    echo cpHeap
    echo rtiStack

proc getName(t: ptr RTInfoObj): string =
  let tn = t.tic.name
  let pn = t.prev.name
  let name = (if tn.len==0:"" else: $tn & ":") & (if pn.len==0:"" else: $pn & "-") & $t.curr.name
  if t.prev.toDropTimer:
    result = "{" & name & "}"
  else: result = name

var hs = initTable[string,RTInfoObj]()
proc makeHotspotTable(lrti: List[RTInfoObj]): tuple[ns:int64,oh:int64] =
  var nstot = int64 0
  var ohtot = int64 0
  for ri in lrti:
    let nc = ri.count
    if ri.istic or nc==0: continue
    let
      f0 = splitFile(ri.prev.loc.filename)[1]
      l0 = ri.prev.loc.line
      f = splitFile(ri.curr.loc.filename)[1]
      l = ri.curr.loc.line
      loc = f0 & "(" & $l0 & "-" & (if f==f0:"" else:f) & $l & ") " & getName(unsafeAddr ri)
      coh = ri.childrenOverhead
      soh = ri.overhead
      nsec = ri.nsec
      ns = nsec - coh
      oh = soh + coh
    nstot += ns
    ohtot += oh
    hs.withValue(loc, t): # t is found value for loc
      t.nsec += ri.nsec
      t.flops += ri.flops
      t.overhead += ri.overhead
      t.childrenOverhead += ri.childrenOverhead
      t.count += ri.count
      for i in 0..<ri.children.len:
        var j = 0
        while j < t.children.len:
          if ri.children[i].curr.loc != t.children[int32 j].curr.loc: inc j; continue
          if ri.children[i].tic.name != t.children[int32 j].tic.name: inc j; continue
          break
        if j < t.children.len:
          t.children[int32 j].nsec += ri.children[i].nsec
          t.children[int32 j].overhead += ri.children[i].overhead
        else:
          t.children.add ri.children[i]
    do: # loc not found
      hs[loc] = ri
    let tot = makeHotSpotTable(List[RTInfoObj](ri.children))
  return (nstot, ohtot)

proc echoHotspots* =
  hs.clear
  let tot = makeHotspotTable(rtiStack)
  #let nstot = tot.ns
  let ohtot = tot.oh
  let nloc = 24
  var skeys = newSeq[tuple[ns:int64,loc:string]]()
  var maxcount = 0
  var nstot = 0.0
  for k,v in hs:
    if v.children.len>0:
      let ins = v.nsec - v.childrenOverhead
      skeys.add (ns: ins, loc: "incl" & k)
    #if v.children.len>0:
    var sns = v.nsec
    if v.children.len==0: sns -= v.childrenOverhead
    for i in 0..<v.children.len:
      sns -= v.children[i].nsec + v.children[i].overhead
    skeys.add (ns: sns, loc: "self" & k)
    maxcount = max(maxcount, int v.count)
    nstot += sns
  skeys.sort(system.cmp, Sortorder.Descending)
  let countdigits = max(5, 1 + int log10(float maxcount))
  echo "Profile time: ", (1e-9*nstot)|(0,6), " overhead: ", (1e-9*ohtot)|(0,6), " total: ", (1e-9*(nstot+ohtot))|(0,6)
  echo "% time  % self count   Mf/s  #child  location                name  (#child: I=inclusive, S=self)"
  var tsns = 0.0
  for nk in skeys:
    let incl = nk.loc.startsWith("incl")
    let key = nk.loc[4..^1]
    hs.withValue(key, t):
      let
        pct = 100.0 * nk.ns / nstot
        count = t.count|countdigits
        #coh = t.childrenOverhead
        #soh = t.overhead
        #oh = soh + coh
        #ohpct = 100.0 * oh / nstot
        mfs = 1000.0 * t.flops / nk.ns.float
        mf = (if t.flops>0.0: mfs.int|7 else: "       ")
        nc = (if t.children.len>0: t.children.len|4 else: "    ")
        f0 = splitFile(t.prev.loc.filename)[1]
        l0 = t.prev.loc.line
        f = splitFile(t.curr.loc.filename)[1]
        l = t.curr.loc.line
        loc = f0 & "(" & $l0 & "-" & (if f==f0:"" else:f) & $l & ")"
        lc = loc|(-nloc,'.')
        nm = getName(t)
      #echo &"{pct:6.3f} {ohpct:6.3f} {count} {mf} {nc} {lc} {nm}"
      if incl:
        echo &"{pct:6.3f}         {count} {mf} {nc} I {lc} {nm}"
      else:
        tsns += nk.ns
        let tsnspct = 100.0 * tsns / nstot
        echo &"{pct:6.3f} {tsnspct:7.3f} {count} {mf} {nc} S {lc} {nm}"

proc echoProf*(def = 0) =
  case intParam("prof",def)
  of 1: echoHotspots()
  of 2: echoTimers()
  else: discard

when isMainModule:
  import os
  proc test =
    tic("test")
    sleep(100)
    toc("sleep 100")
    block:
      tic("block")
      sleep(1000)
      toc("sleep 1000")
    toc("sleep block")
  proc f(n = 1) =
    #echo "*** f"
    tic("f" & $n)
    sleep(10*n)
    toc("end")
  proc g(n = 5) =
    #echo "*** g"
    tic("g" & $n)
    sleep(1*n)
    toc("end")
  proc r(n = 2)
  proc s(n = 2) =
    #echo "*** s ",n
    tic("s " & $n)
    if n > 0:
      toc("s work")
      r(n-1)
      toc("s r")
    toc("s end")
  proc r(n = 2) =
    #echo "*** r ",n
    tic("r " & $n)
    if n > 0:
      sleep 20
      toc("r work")
      s(n)
      toc("r s")
    toc("r end")
  proc loop =
    tic("loop")
    for i in 0..3:
      sleep 40
      toc("no local tic")
    for i in 0..4:
      f()
      g(i)
    toc("f g")
    for i in 0..3:
      tic()
      sleep 15
      toc("sleep")
    toc("end")
  proc longloop =
    tic("longloop")
    for i in 1..10000:
      tic()
      f(0)
      toc("f")
      g(0)
      toc("g")
    toc("end")
  proc test2 =
    tic("test2")
    sleep(100)
    toc()
    f()
    toc("f 1")
    block:
      tic()
      sleep(100)
      toc("sleep 2")
      var runf {.global.} = 2
      if runf > 0:
        f()
        toc("f 2")
        dec runf
    toc("block")
  DropWasteTimerRatio = 0.10
  test()
  echoTimers(aggregate=false)
  test()
  echoTimers(aggregate=false)
  resetTimers()
  test()
  echoTimers(aggregate=false)
  test()
  test2()
  echoTimers()
  tic()
  test2()
  toc("test2 1")
  for i in 0..1:
    test2()
  toc("test2 loop")
  test2()
  toc("test2 2")
  echoTimers(0.1)
  resetTimers()
  r()
  echoTimers()
  loop()
  echoTimers()
  for i in 0..4:
    tic("llloop")
    longloop()
    toc("one")
  toc("end")
  echoTimers()
  echoTimersRaw()
  echoHotspots()
