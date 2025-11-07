import base/[threading,metaUtils,stdUtils,basicOps,profile]
import times
import os
import macros
#import strUtils
import comms/comms
getOptimPragmas()

template sum*(c:Comm, v:var SomeNumber) = c.allReduce(v)
template sum*(c:Comm, v:ptr float32, n:int) = c.allReduce(v,n)
template sum*(c:Comm, v:ptr float64, n:int) = c.allReduce(v,n)
#template sum*(c:Comm, v:ptr array, n:int) =
#  c.allReduce(v[][0].addr, n*v[].len)
template sum*[R,T](c:Comm, v:ptr array[R,T], n:int) =
  #treerep:
  #  v[][rangeLow(R)]
  c.allReduce(v[][rangeLow(R)].addr, n*rangeLen(R))
#template qmpSum*(v:ptr tuple, n:int):untyped =
#  qmpSum(v[][0].addr, n*(sizeOf(v) div sizeOf(v[0])))
#template qmpSum*(v:ptr object, n:int):untyped =
#  qmpSum(v[][].addr, n)
#template qmpSum*(v: object) =
#  qmpSum(asNumberPtr(v), numNumbers(v))
#template qmpSum*(v:ptr typed, n:int):untyped =
#  qmpSum(v[][].addr, n)
#template QmpSum(v:array[int,int]):untyped =
#  var tQmpSumDoubleArray(v)
#template qmpSum*[I,T](v:array[I,T]):untyped =
#  qmpSum(v[0].addr, v.len)
#template qmpSum*(v:openArray[float64]):untyped =
#  QmpSumDoubleArray(v[0].addr,v.len.cint)
template sum*[T](c:Comm, v:seq[T]):untyped =
  sum(c, v[0].addr, v.len)
#template qmpSum*[I,T](v:seq[array[I,T]]):untyped =
#  qmpSum(v[0][0].addr, v.len.cint*sizeOf(v[0]))
#template qmpSum*(v:openArray[array]):untyped =
#  qmpSum(v[0][0].addr, v.len.cint*sizeOf(v[0]))
#template qmpSum*(v:tuple):untyped =
#  qmpSum(v[0].addr, sizeOf(v) div sizeOf(v[0]))
#template qmpSum*[T](v:T):untyped =
#template qmpSum*(v:typed):untyped =
#  qmpSum(v[])
#template qmpSum*[T](v:T):untyped =
#  qmpSum(v[])
template sum*(c:Comm, v:typed) =
  when numberType(v) is float64:
    sum(c, cast[ptr float64](addr v), sizeof(v) div sizeof(float64))
  elif numberType(v) is float32:
    sum(c, cast[ptr float32](addr v), sizeof(v) div sizeof(float32))
  else:
    sum(v[])

#template qmpMax*(v:float32):untyped = QmpMaxFloat(v.addr)
#template qmpMax*(v:float64):untyped = QmpMaxDouble(v.addr)
#template qmpMin*(v:float32):untyped = QmpMinFloat(v.addr)
#template qmpMin*(v:float64):untyped = QmpMinDouble(v.addr)

proc unwrap(x:NimNode):seq[NimNode] =
  result = @[]
  let t = x.getType
  #echo x.treeRepr
  #echo t.treeRepr
  #echo t.typekind
  if t.typekind==ntyTuple:
    let n = t.len - 1
    for i in 0..<n:
      let id = newLit(i)
      result.add(quote do:
        `x`[`id`]
        )
  else:
    result.add(quote do:
      `x`[]
      )
  #echo result.repr

macro rankSumN*(comm:Comm, a:varargs[typed]):auto =
  #echo "rankSum: ", a.repr
  #echo a.treeRepr
  var i0 = 0
  let t0 = a[0].getType
  if a.len==1:
    let a0 = a[0]
    result = quote do:
      if threadNum==0:
        `comm`.sum(`a0`)
    return result
  #echo t0.repr
  #echo t0.typekind
  if t0.typekind==ntyFloat32 or t0.typekind==ntyFloat:
    #echo "got float"
    i0 = -1
    for i in 1..<a.len:
      #echo a[i].getType.repr
      if a[i].getType.repr != t0.repr:
        if a[i].getType is float32|float64:
          quit("can't mix float types in rankSum")
        i0 = i
        break
  if i0<0:
    let b = newNimNode(nnkBracket)
    var s = newNimNode(nnkStmtList)
    let t = ident("t")
    for i in 0..<a.len:
      b.add a[i]
      let ai = a[i]
      let x = quote do:
        `ai` = `t`[`i`]
      s.add x
    result = quote do:
      if threadNum==0:
        var `t` = `b`
        `comm`.sum(`t`)
        `s`
  else:
    result = newCall("rankSum")
    for i in 0..<a.len:
      if i==i0:
        let ai = unwrap(a[i])
        for j in 0..<ai.len:
          result.add(ai[j])
      else:
        result.add(a[i])
  #echo result.repr
macro rankSum*(comm: Comm, a:varargs[untyped]):auto =
  if a.len==1:
    let a0 = a[0]
    result = quote do:
      if threadNum==0:
        `comm`.sum(`a0`)
  else:
    result = newCall(ident("rankSumN"),comm)
    for v in a: result.add v

template rankSum*(a:varargs[untyped]) =
  let comm = getDefaultComm()
  comm.rankSum(a)

#[
#var count = 0
template threadRankSum1*[T](comm: Comm, a: T) =
  mixin rankSum
  #[
  #if threadNum==0: inc count
  #threadBarrier()
  threadLocals.share[threadNum].p = a.addr
  #echoAll count, " ", myrank, " ", threadNum, " v: ", cast[ByteAddress](a.addr)
  #echoAll count, " ", myrank, " ", threadNum, " s: ", ptrInt(threadLocals.share)
  if threadNum==0:
    #threadBarrier()
    t0wait()
    for i in 1..<numThreads:
      #echo "test1"
      #echo count, " ", i, " ", cast[ByteAddress](threadLocals.share[i].p)
      a += cast[ptr type(a)](threadLocals.share[i].p)[]
      #echo "test2"
    rankSum(a)
    #threadBarrier()
    twait0()
    #threadBarrier()
    t0wait()
  else:
    #threadBarrier()
    t0wait()
    #threadBarrier()
    twait0()
    a = cast[ptr type(a)](threadLocals.share[0].p)[]
    #threadBarrier()
    t0wait()
  ]#
  var ta{.global.}: type(a)
  #var ta2{.global.}:array[512,type(a)]
  if threadNum==0:
    tic("threadRankSum1")
    t0wait()
    toc("t0wait")
    for i in 1..<numThreads:
      a += cast[ptr type(a)](threadLocals.share[i].p)[]
    toc("sum")
    rankSum(comm,a)
    toc("rankSum")
    ta = a
    twait0()
    toc("twait0")
  else:
    threadAtomicWrite:
      threadLocals.share[threadNum].p = a.addr
    #ta2[threadNum] = a
    t0wait()
    twait0()
    a = ta
]#

#template threadRankSum1x*[T](comm: Comm, a: T) =
proc threadRankSum1x*[T](comm: Comm, a: var T) =
  mixin rankSum
  tic("threadRankSum1x")
  threadSum0(a)
  toc("threadSum0")
  if threadNum==0:
    rankSum(comm,a)
  toc("rankSum")
  threadBroadcast(a)
  toc("threadBroadcast")

#[
proc threadRankSumN*(comm: Comm, a: NimNode): auto =
  echo a.treeRepr
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
    comm.rankSum(`gi`)
  result.add(m)
  #echo result.treeRepr
]#

#[
macro threadRankSum*(comm: Comm, a0:varargs[untyped]):auto =
  var a = a0
  #let ta = a.getType()
  #echo a.treerepr
  #echo ta.treerepr
  #echo a[0].treerepr
  #echo a[0].getType.treerepr
  #if ta.kind == nnkBracketExpr and ta[0].repr == "varargs":
  #  a = a[0][1]
  #echo a.treerepr
  #echo a.getType.treerepr
  if a.len==1:
    #template trs1(comm:Comm,x:untyped):untyped = threadRankSum1(comm, x)
    #result = getAst(trs1(comm,a[0]))
    #result = getAst(threadRankSum1(comm,a[0]))
    result = newCall(bindSym("threadRankSum1"), comm, a[0])
    echo result.repr
  else:
    echo "threadRankSumN not implemented"
    doAssert(false)
    #result = threadRankSumN(comm,a)

template threadRankSum*(a:varargs[untyped]) =
  static: echo a.getType.treerepr
  when a.len==1:
    let comm = getDefaultComm()
    threadRankSum1(comm, a[0])
  else:
    echo "threadRankSumN not implemented"
    doAssert(false)
    #result = threadRankSumN(comm,a)
]#
#template threadRankSum*(a:typed):auto =

template threadRankSum*(c: Comm, a: typed) =
  threadRankSum1x(c, a)
template threadRankSum*(a: typed) =
  let comm = getDefaultComm()
  threadRankSum1x(comm, a)

macro rankMax*(a:varargs[untyped]):auto =
  if a.len==1:
    let a0 = a[0]
    result = quote do:
      if threadNum==0:
        qmpMax(`a0`)
  else:
    error("rankMax not imlemented for multiple arguments.")
    #result = newCall(ident("rankMaxN"))
    #for v in a: result.add v

template threadRankMax1*(a:untyped):untyped =
  var ta{.global.}:type(a)
  if threadNum==0:
    t0wait()
    #threadBarrier()
    for i in 1..<numThreads:
      let c = cast[ptr type(a)](threadLocals.share[i].p)[]
      if a < c: a = c
    rankMax(a)
    ta = a
    twait0()
    #threadBarrier()
  else:
    threadLocals.share[threadNum].p = a.addr
    t0wait()
    #threadBarrier()
    twait0()
    #threadBarrier()
    a = ta

macro threadRankMax*(a:varargs[untyped]):auto =
  if a.len==1:
    template trm1(x:untyped):untyped = threadRankMax1(x)
    result = getAst(trm1(a[0]))
  else:
    error("threadRankMax not imlemented for multipel arguments.")
    #result = threadRankMaxN(a)

when isMainModule:
  commsInit()
  echo "rank ", myRank, "/", nRanks
  printf("rank %i/%i\n", myRank, nRanks)
  threads:
    echo threadNum, "/", numThreads
    let n = nRanks * numThreads
    let s = (n*(n-1)) div 2
    var x = myRank*numThreads + threadNum
    threadRankSum(x)
    echo threadNum, ": ", x, "  ", s
    threadRankSum(x)
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

    f = 0.1
    threadBarrier()
    t0 = epochTime()
    for i in 1..nrep:
      threadRankSum(f)
    t1 = epochTime()
    echo "threadRankSum(float) time: ", int(1e9*(t1-t0)/nrep.float), " ns"

    f = 0.1
    threadBarrier()
    if threadNum==0:
      t0 = epochTime()
      for i in 1..nrep:
        rankSum(f)
      t1 = epochTime()
      echo "rankSum(float) time: ", int(1e9*(t1-t0)/nrep.float), " ns"
    threadBarrier()

  commsFinalize()
