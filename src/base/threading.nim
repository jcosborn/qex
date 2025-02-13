import times
import strUtils
import stdUtils
import macros
import omp
import metaUtils
import base/basicOps
import bitops
getOptimPragmas

type
  ThreadShare* = object
    p*:pointer
    counter*:int
    extra*:int
  ThreadObj* = object
    threadNum*:int
    numThreads*:int
    #counter*: int
    share*:ptr cArray[ThreadShare]

var threadNum*{.threadvar.}:int
var numThreads*{.threadvar.}:int
var threadLocals*{.threadvar.}:ThreadObj
var inited = false
var ts: pointer = nil
var nts = 0

proc allocTs* {.alwaysInline.} =
  if numThreads > nts and threadNum == 0:
    if ts == nil:
      ts = allocShared(numThreads*sizeof(ThreadShare))
    else:
      ts = reallocShared(ts, numThreads*sizeof(ThreadShare))
    nts = numThreads

#template initThreadLocals*(ts:seq[ThreadShare]) =
template initThreadLocals* =
  bind ts
  threadLocals.threadNum = threadNum
  threadLocals.numThreads = numThreads
  #threadLocals.share = cast[ptr cArray[ThreadShare]](ts[0].addr)
  threadLocals.share = cast[ptr cArray[ThreadShare]](ts)
  threadLocals.share[threadNum].p = nil
  threadLocals.share[threadNum].counter = 0
proc init =
  inited = true
  threadNum = 0
  numThreads = 1
  #var ts = newSeq[ThreadShare](numThreads)
  #initThreadLocals(ts)
  allocTs()
  initThreadLocals()
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

macro emitStackTraceX(x: typed): untyped =
  template est(x) =
    {.emit: "// instantiationInfo: " & x.}
  let ii = x.repr.replace("\n","")
  result = getAst(est(ii))

template emitStackTrace: untyped =
  emitStackTraceX(instantiationInfo(-1))
  emitStackTraceX(instantiationInfo(-2))
  emitStackTraceX(instantiationInfo(-3))

template threads*(body:untyped):untyped =
  checkInit()
  doAssert(numThreads==1)
  #let tidOld = threadNum
  #let nidOld = numThreads
  #let tlOld = threadLocals
  #proc tproc2{.genSym,inline.} =
  #  body
  proc tproc{.genSym.} =
    emitStackTrace()
    #var ts:seq[ThreadShare]
    ompParallel:
      threadNum = ompGetThreadNum()
      numThreads = ompGetNumThreads()
      #if threadNum==0: ts.newSeq(numThreads)
      allocTs()
      threadBarrierO()
      #initThreadLocals(ts)
      initThreadLocals()
      threadBarrierO()
      #echoAll threadNum, " s: ", ptrInt(threadLocals.share)
      body
      #tproc2()
      threadBarrierO()
  tproc()
  #threadNum = tidOld
  #numThreads = nidOld
  #threadLocals = tlOld
  threadNum = 0
  numThreads = 1
  initThreadLocals()
template threads*(x0:untyped;body:untyped):untyped =
  checkInit()
  #let tidOld = threadNum
  #let nidOld = numThreads
  #let tlOld = threadLocals
  proc tproc(xx:var type(x0)) {.genSym.} =
    #var ts:seq[ThreadShare]
    ompParallel:
      threadNum = ompGetThreadNum()
      numThreads = ompGetNumThreads()
      #if threadNum==0: ts.newSeq(numThreads)
      allocTs()
      threadBarrierO()
      #initThreadLocals(ts)
      initThreadLocals()
      threadBarrierO()
      #echoAll threadNum, " s: ", ptrInt(threadLocals.share)
      subst(x0,xx):
        body
      threadBarrierO()
  tproc(x0)
  #threadNum = tidOld
  #numThreads = nidOld
  #threadLocals = tlOld
  threadNum = 0
  numThreads = 1
  initThreadLocals()

template nothreads*(body: untyped): untyped =
  ## convenient way to turn off threading
  block:
    body


template getMaxThreads*() = ompGetMaxThreads()
template threadBarrierO* = ompBarrier
template threadMaster*(x:untyped) = ompMaster(x)
template threadSingle*(x:untyped) = ompSingle(x)
template threadCritical*(x:untyped) = ompCritical(x)
template threadFlush* = ompFlush
template threadFlushRelease* = ompFlushRelease
template threadFlushAcquire* = ompFlushAcquire
template threadFlushSeqCst* = ompFlushSeqCst
template threadAtomicRead*(body:typed) = ompAtomicRead(body)
template threadAtomicWrite*(body:typed) = ompAtomicWrite(body)

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

template t0waitB* = threadBarrier()
#template t0waitO* = t0waitB()
template t0waitA* =
  if numThreads > 1:
    if threadNum==0:
      inc threadLocals.share[0].counter
      let tbar0 = threadLocals.share[0].counter
      for b in 1..<numThreads:
        let p = threadLocals.share[b].counter.addr
        while true:
          #fence()
          #ompAcquire
          #if p[] >= tbar0: break
          var t {.noInit.}: type(p[])
          ompAtomicRead: t = p[]
          if t >= tbar0: break
    else:
      #inc threadLocals.share[threadNum].counter
      #fence()
      #ompRelease
      let t = threadLocals.share[threadNum].counter + 1
      ompAtomicWrite:
        threadLocals.share[threadNum].counter = t
template t0waitC* =
  if numThreads > 1:
    if threadNum==0:
      inc threadLocals.share[0].counter
      let tbar0 = threadLocals.share[0].counter
      var left = numThreads - 1
      for i in 1..<numThreads: threadLocals.share[i].extra = 1
      while left > 0:
        for i in 1..<numThreads:
          if threadLocals.share[i].extra > 0:
            let p = threadLocals.share[i].counter.addr
            var t {.noInit.}: type(p[])
            ompAtomicRead: t = p[]
            if t >= tbar0:
              threadLocals.share[i].extra = 0
              dec left
    else:
      let t = threadLocals.share[threadNum].counter + 1
      ompAtomicWrite:
        threadLocals.share[threadNum].counter = t
#template t0wait* = t0waitA()
#template t0wait* = t0waitB()
template t0wait* = t0waitC()

template twait0B* = threadBarrier()
template twait0A* =
  if numThreads > 1:
    if threadNum==0:
      #inc threadLocals.share[0].counter
      #fence()
      #ompRelease
      let t = threadLocals.share[0].counter + 1
      ompAtomicWrite:
        threadLocals.share[0].counter = t
    else:
      inc threadLocals.share[threadNum].counter
      let tbar0 = threadLocals.share[threadNum].counter
      let p = threadLocals.share[0].counter.addr
      while true:
        #fence()
        #ompAcquire
        #if p[] >= tbar0: break
        var t {.noInit.}: type(p[])
        ompAtomicRead: t = p[]
        if t >= tbar0: break
template twait0* = twait0A()
#template twait0* = twait0B()

template threadBarrierA* =
  threadFlushRelease
  t0waitA
  twait0A
  threadFlushAcquire
template threadBarrier* = threadBarrierA
#template threadBarrier* = ompBarrier

template threadSum01A0*[T](a: T) =
  ## sum value with result on thread 0, atomic version
  if threadNum==0:
    tic("threadSum01A0")
    t0wait()
    toc("t0wait")
    for i in 1..<numThreads:
      var p{.noInit.}: pointer
      threadAtomicRead:
        p = threadLocals.share[i].p
      a += cast[ptr T](p)[]
    toc("sum")
    twait0()
    toc("twait0")
  else:
    threadAtomicWrite:
      threadLocals.share[threadNum].p = a.addr
    t0wait()
    twait0()

template threadSum01A*[T](a: T) =
  ## sum value with result on thread 0, atomic version
  if threadNum==0:
    tic("threadSum01A")
    inc threadLocals.share[0].counter
    let tbar0 = threadLocals.share[0].counter
    for b in 1..<numThreads:
      let pc = threadLocals.share[b].counter.addr
      while true:
        var t {.noInit.}: type(pc[])
        ompAtomicRead: t = pc[]
        if t >= tbar0: break
      var p{.noInit.}: pointer
      threadAtomicRead:
        p = threadLocals.share[b].p
      a += cast[ptr T](p)[]
    toc("sum")
    twait0()
    toc("twait0")
  else:
    threadAtomicWrite:
      threadLocals.share[threadNum].p = a.addr
    let t = threadLocals.share[threadNum].counter + 1
    threadAtomicWrite:
      threadLocals.share[threadNum].counter = t
    twait0()

template threadSum01B*[T](a: T) =
  ## sum value with result on thread 0, barrier version
  block:
    tic("threadSum01B")
    if threadNum!=0:
      threadLocals.share[threadNum].p = a.addr
    threadBarrier()
    toc("threadBarrier first")
    if threadNum==0:
      for i in 1..<numThreads:
        a += cast[ptr T](threadLocals.share[i].p)[]
    toc("sum")
    threadBarrier()
    toc("threadBarrier last")

template threadSum01T*[T](a: T) =
  ## sum value with result on thread 0, tree version
  threadAtomicWrite:
    threadLocals.share[threadNum].p = a.addr
  var c = threadLocals.share[threadNum].counter
  var done = false
  var b = 1
  while b < numThreads:
    inc c
    if not done:
      let o = bitxor(threadNum, b)
      if bitand(threadNum, b) == 0:
        if o < numThreads:
          while true:
            var d{.noInit.}: int
            threadAtomicRead:
              d = threadLocals.share[o].counter
            if d >= c: break
          var p{.noInit.}: pointer
          threadAtomicRead:
            p = threadLocals.share[o].p
          a += cast[ptr T](p)[]
        threadAtomicWrite:
          threadLocals.share[threadNum].counter = c
      else:
        threadAtomicWrite:
          threadLocals.share[threadNum].counter = c
        while true:
          var d{.noInit.}: int
          threadAtomicRead:
            d = threadLocals.share[o].counter
          if d >= c: break
        done = true
    else:
      threadLocals.share[threadNum].counter = c
    b *= 2

template threadSum01*(a: auto) = threadSum01A(a)
#template threadSum01*(a: auto) = threadSum01B(a)
#template threadSum01*(a: auto) = threadSum01T(a)
template threadSum0*(a: auto) = threadSum01(a)

# threadMax0 FIXME

template threadBroadcast1A*[T](a: T) =
  if threadNum==0:
    tic("threadBroadcast1A")
    threadAtomicWrite:
      threadLocals.share[0].p = a.addr
    twait0()
    toc("twait0")
    t0wait()
    toc("t0wait")
  else:
    twait0()
    var p{.noInit.}: pointer
    threadAtomicRead:
      p = threadLocals.share[0].p
    a = cast[ptr T](p)[]
    t0wait()
template threadBroadcast1*(a: auto) = threadBroadcast1A(a)
template threadBroadcast*(a: auto) = threadBroadcast1(a)

macro threadSum*(a:varargs[untyped]):auto =
  #echo a.treeRepr
  result = newNimNode(nnkStmtList)
  var sum = newNimNode(nnkStmtList)
  let tid = ident("threadNum")
  let nid = ident("numThreads")
  let p = newLit(1)
  for i in 0..<a.len:
    let gi = ident("g" & $i)
    let ai = a[i]
    result.add(quote do:
      var `gi`{.global.}:array[`p`*512,type(`ai`)]
      `gi`[`p`*`tid`] = `ai`
      #deepCopy(`gi`[`p`*`tid`], `ai`)
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
    let gi = ident("g" & $i)
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

type ThreadSingle*[T] = distinct T
template `[]`*[T](x: ThreadSingle[T]): auto = T(x)
template `[]=`*[T](x: ThreadSingle[T], y: auto): auto =
  T(x) = y
template `=`*(x: var ThreadSingle, y: auto) =
  threadSingle:
    x[] = y
template `+=`*(x: var ThreadSingle, y: auto) =
  threadSingle:
    x[] += y
converter fromThreadSingle*[T](x: ThreadSingle[T]): T = T(x)

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
