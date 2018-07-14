import threading
export threading
import comms/comms
import stdUtils
import times
import strUtils
import macros
import os
import algorithm
#include system/timers

when true:
  type
    TicType* = distinct float
  template getTics*(): untyped = TicType(epochTime())
  #template toSeconds*(x: TicType): untyped = 1e-9 * x.float
  template `-`*(x,y: TicType): untyped = TicType(float(x) - float(y))
  template ticDiffSecs*(x,y: untyped): untyped = float(x - y)
else:
  type
    TicType* = distinct Ticks
  template getTics*(): untyped = TicType(getTicks())
  #template toSeconds*(x: TicType): untyped = 1e-9 * x.float
  template `-`*(x,y: TicType): untyped = TicType(Ticks(x) - Ticks(y))
  template ticDiffSecs*(x,y: untyped): untyped = 1e-9 * (x - y).float

type
  TimerInfo* = object
    seconds*: float
    flops*: float
    count*: int
    name*: string
    timerP*: ptr TicType
    start*: ptr TimerInfo
    prev*: ptr TimerInfo
    location*: type(instantiationInfo())

var ticSeq* = newSeq[ptr TimerInfo](0)
var tocSeq* = newSeq[ptr TimerInfo](0)

#var cFunc {.importC: "__func__", noDecl.}: cstring
#template getFunctionName(): untyped = $cFunc
template getFunctionName(): untyped = ""

#macro getFunctionName:auto =
  #let s = callsite()
  #echo s.repr
  #return newLit("")

#{.emit:"#include <stdio.h>".}
#template ctrace(s:cstring):untyped =
#  {.emit:"""#ifdef __func__
#printf("%s(%i): %s\n", __func__, __LINE__, `s`);
#endif""".}
#macro ctrace(s:string):auto =

proc initTic*(ti: var TimerInfo, s: string, timer: var TicType, inst: tuple) =
  ti.name = s
  ti.timerP = timer.addr
  ti.location = inst
  #echo "tic", ti.location, cast[ByteAddress](ti.timerP)
  ticSeq.add ti.addr
template tic*(n= -1): untyped =
  var ti {.global.}: TimerInfo
  var timer {.global,inject.}: TicType
  if threadNum==0:
    if ti.timerP==nil:
      ti.initTic(getFunctionName(), timer, instantiationInfo(n))
    inc ti.count
    timer = getTics()

proc initToc*(ti: var TimerInfo, s: string, timer: var TicType, inst: tuple) =
  ti.name = s
  ti.timerP = timer.addr
  ti.location = inst
  #echo "toc", ti.location, cast[ByteAddress](ti.timerP)
  tocSeq.add ti.addr
template tocI(f: untyped; s=""; n= -1): untyped =
  var ti {.global.}: TimerInfo
  if threadNum==0:
    if ti.timerP==nil:
      ti.initToc(s, timer, instantiationInfo(n))
    ti.flops += f.float
    inc ti.count
    ti.seconds += ticDiffSecs(getTics(), timer)
template toc*(s=""; n= -1; flops:untyped):untyped = tocI(flops, s, n-1)
template toc*(n= -1; flops:untyped):untyped = tocI(flops, "", n-1)
template toc*(s:string; n:int):untyped = tocI(0, s, n-1)
template toc*(s:string):untyped = tocI(0, s, -2)
template toc*(n:int):untyped = tocI(0, "", n-1)
template toc*():untyped = tocI(0, "", -2)

template getElapsedTime*(): untyped {.dirty.} =
  mixin timer
  ticDiffSecs(getTics(), timer)

proc resetTimers* =
  for ti in tocSeq:
    ti.seconds = 0.0
    ti.count = 0
    ti.flops = 0.0

proc echoTimers* =
  for ti in tocSeq:
    for t in ticSeq:
      if t.timerP == ti.timerP: ti.start = t; break
  tocSeq.sort do (x, y: ptr TimerInfo) -> int:
    #result = cmp(x.location.filename, y.location.filename)
    result = cmp(x.timerP.ptrInt, y.timerP.ptrInt)
    #if result==0:
    #  result = cmp(x.start.location.line, y.start.location.line)
    if result==0:
      result = cmp(x.location.line, y.location.line)
  var prev:ptr TimerInfo = nil
  for ti in tocSeq:
    if prev!=nil and ti.timerP==prev.timerP:
      ti.prev = prev
    prev = ti
    #echo ti.location.line, " ", ti.timerP.ptrInt
  tocSeq.sort do (x, y: ptr TimerInfo) -> int:
    result = cmp(x.location.filename, y.location.filename)
    if result==0:
      result = cmp(x.location.line, y.location.line)
    if result==0:
      result = cmp(x.timerP.ptrInt, y.timerP.ptrInt)

  echo '='.repeat(76)
  echo "file(lines)"|(-24), "microsecs"|10, "count"|8, "ns/count"|12, "mf"|8, " label"
  echo '='.repeat(76)
  #var tot = 0.0
  for ti in tocSeq:
    let tc = ti.start
    let f = splitFile($ti.location.filename)[1]
    let l = $ti.location.line
    var l0 = tc.location.line
    var secs = ti.seconds
    if ti.prev!=nil:
      l0 = ti.prev.location.line
      secs -= ti.prev.seconds
    let loc = f & "(" & $l0 & "-" & $l & ")"
    let st = (secs*1e6+0.5).int | 10
    let c = ti.count | 6
    let sc = (secs*1e9/(ti.count.float+1e-9)).int | 10
    var mf = (ti.flops*1e-6/(secs+1e-9)).int | 8
    if ti.flops<0 or ti.count==0: mf = "-"|8
    echo loc|(-24,'.'), st, " /", c, " =", sc, mf, " ", ti.name
    #tot += secs
  #echo "total"|(-24,'.'), (tot*1e6).int | 10
  echo '='.repeat(76)

when isMainModule:
  import os
  proc test =
    tic()
    sleep(100)
    toc()
    block:
      tic()
      sleep(1000)
      toc()
    toc("sleep2")
  test()
  echoTimers()
  test()
  echoTimers()
  resetTimers()
  echoTimers()
  test()
  echoTimers()
