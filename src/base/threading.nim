import times
import strUtils
import stdUtils
import macros
import omp
import metaUtils

type
  ThreadShare* = object
    p*:pointer
    counter*:int
  ThreadObj* = object
    threadNum*:int
    numThreads*:int
    share*:ptr cArray[ThreadShare]

var threadNum*{.threadvar.}:int
var numThreads*{.threadvar.}:int
var threadLocals*{.threadvar.}:ThreadObj
var inited = false

template initThreadLocals*(ts:seq[ThreadShare]):untyped =
  threadLocals.threadNum = threadNum
  threadLocals.numThreads = numThreads
  threadLocals.share = cast[ptr cArray[ThreadShare]](ts[0].addr)
  threadLocals.share[threadNum].p = nil
  threadLocals.share[threadNum].counter = 0
proc init =
  inited = true
  threadNum = 0
  numThreads = 1
  var ts = newSeq[ThreadShare](numThreads)
  initThreadLocals(ts)
template threadsInit* =
  if not inited:
    init()
template checkInit* =
  threadsInit()
  #if not inited:
    #let ii = instantiationInfo()
    #let ln = ii.line
    #let fn = ii.filename[0 .. ^5]
    #echo format("error: $#($#): threads not initialized",fn,ln)
    #quit(-1)

template threads*(body:untyped):untyped =
  checkInit()
  let tidOld = threadNum
  let nidOld = numThreads
  let tlOld = threadLocals
  #proc tproc2{.genSym,inline.} =
  #  body
  proc tproc{.genSym.} =
    var ts:seq[ThreadShare]
    ompParallel:
      threadNum = ompGetThreadNum()
      numThreads = ompGetNumThreads()
      if threadNum==0: ts.newSeq(numThreads)
      threadBarrierO()
      initThreadLocals(ts)
      #echoAll threadNum, " s: ", ptrInt(threadLocals.share)
      body
      #tproc2()
      threadBarrierO()
  tproc()
  threadNum = tidOld
  numThreads = nidOld
  threadLocals = tlOld
template threads*(x0:untyped;body:untyped):untyped =
  checkInit()
  let tidOld = threadNum
  let nidOld = numThreads
  let tlOld = threadLocals
  proc tproc(xx:var type(x0)) {.genSym.} =
    var ts:seq[ThreadShare]
    ompParallel:
      threadNum = ompGetThreadNum()
      numThreads = ompGetNumThreads()
      if threadNum==0: ts.newSeq(numThreads)
      threadBarrierO()
      initThreadLocals(ts)
      #echoAll threadNum, " s: ", ptrInt(threadLocals.share)
      subst(x0,xx):
        body
  tproc(x0)
  threadNum = tidOld
  numThreads = nidOld
  threadLocals = tlOld

template getMaxThreads*() = ompGetMaxThreads()
template threadBarrierO* = ompBarrier
template threadMaster*(x:untyped) = ompMaster(x)
template threadSingle*(x:untyped) = ompSingle(x)
template threadCritical*(x:untyped) = ompCritical(x)

template threadDivideLow*(x,y: untyped): untyped =
  x + (threadNum*(y-x)) div numThreads
template threadDivideHigh*(x,y: untyped): untyped =
  x + ((threadNum+1)*(y-x)) div numThreads


proc tForX*(index,i0,i1,body:NimNode):NimNode =
  return quote do:
    let d = 1+`i1` - `i0`
    let ti0 = `i0` + (threadNum*d) div numThreads
    let ti1 = `i0` + ((threadNum+1)*d) div numThreads
    for `index` in ti0 ..< ti1:
      `body`
macro tFor*(index,i0,i1: untyped; body: untyped): untyped =
  result = tForX(index, i0, i1, body)
macro tFor*(index: untyped; slice: Slice; body: untyped): untyped =
  #echo index.treeRepr
  #echo treeRepr(slice)
  var i0,i1: NimNode
  #echo slice.kind
  if slice.kind == nnkStmtListExpr:
    i0 = slice[1][1]
    i1 = slice[1][2]
  else:
    i0 = slice[1]
    i1 = slice[2]
  result = tForX(index, i0, i1, body)

discard """
iterator `.|`*[S, T](a: S, b: T): T {.inline.} =
  mixin threadNum
  var d = b - T(a)
  var res = T(a) + (threadNum*d) div numThreads
  var bb = T(a) + ((threadNum+1)*d) div numThreads
  while res <= bb:
    yield res
    inc(res)
"""

template t0wait* = threadBarrier()
template t0waitX* =
  if threadNum==0:
    inc threadLocals.share[0].counter
    let tbar0 = threadLocals.share[0].counter
    for b in 1..<numThreads:
      let p{.volatile.} = threadLocals.share[b].counter.addr
      while true:
        if p[] >= tbar0: break
  else:
    inc threadLocals.share[threadNum].counter
    #fence()

template twait0* = threadBarrier()
template twait0X* =
  if threadNum==0:
    inc threadLocals.share[0].counter
    #fence()
  else:
    inc threadLocals.share[threadNum].counter
    let tbar0 = threadLocals.share[threadNum].counter
    let p{.volatile.} = threadLocals.share[0].counter.addr
    while true:
      if p[] >= tbar0: break

template threadBarrier* =
  #t0wait
  #twait0
  ompBarrier

macro threadSum*(a:varargs[untyped]):auto =
  #echo a.treeRepr
  result = newNimNode(nnkStmtList)
  var sum = newNimNode(nnkStmtList)
  let tid = ident("threadNum")
  let nid = ident("numThreads")
  let p = newLit(1)
  for i in 0..<a.len:
    let gi = !("g" & $i)
    let ai = a[i]
    result.add(quote do:
      var `gi`{.global.}:array[`p`*512,type(`ai`)]
      #`gi`[`p`*`tid`] = `ai`
      deepCopy(`gi`[`p`*`tid`], `ai`)
      )
    let s = quote do:
      `ai` = `gi`[0]
      for i in 1..<`nid`:
        `ai` += `gi`[`p`*i]
    sum.add(s)
  let m = quote do:
    threadBarrier()
    `sum`
    threadBarrier()
  result.add(m)
  result = newBlockStmt(result)
  #echo result.treeRepr
macro threadSum2*(a:varargs[untyped]):auto =
  #echo a.treeRepr
  result = newNimNode(nnkStmtList)
  var g0 = newNimNode(nnkStmtList)
  var gp = newNimNode(nnkStmtList)
  var a0 = newNimNode(nnkStmtList)
  for i in 0..<a.len:
    let gi = !("g" & $i)
    let ai = a[i]
    let t = quote do:
      var `gi`{.global.}:type(`ai`)
    result.add(t[0])
    let x0 = quote do:
      `gi` = `ai`
    g0.add(x0[0])
    #echo g0.treeRepr
    let xp = quote do:
      `gi` += `ai`
    gp.add(xp[0])
    #echo gp.treeRepr
    let ax = quote do:
      `ai` = `gi`
    a0.add(ax[0])
    #echo a0.treeRepr
  #echo result.treeRepr
  let m = quote do:
    if threadNum==0:
      `g0`
      threadBarrier()
      threadBarrier()
    else:
      threadBarrier()
      {.emit:"#pragma omp critical"}
      block:
        `gp`
      threadBarrier()
    `a0`
  result.add(m)
  #echo result.treeRepr

when isMainModule:
  threadsInit()
  echo threadNum, "/", numThreads
  threads:
    echo threadNum, "/", numThreads
    let n = numThreads
    let s = (n*(n-1)) div 2
    var x = threadNum
    threadSum(x)
    echo threadNum, ": ", x, "  ", s
    threadSum(x)
    echo threadNum, ": ", x, "  ", n*s

    let nrep = 1000

    threadBarrier()
    var t0 = epochTime()
    for i in 1..nrep:
      threadBarrier()
    var t1 = epochTime()
    echo "threadBarrier time: ", int(1e9*(t1-t0)/nrep.float), " ns"

    var f = 0.1
    threadBarrier()
    t0 = epochTime()
    for i in 1..nrep:
      threadSum(f)
    t1 = epochTime()
    echo "threadSum(float) time: ", int(1e9*(t1-t0)/nrep.float), " ns"
