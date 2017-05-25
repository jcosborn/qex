import base/threading
import times
import os
import macros
import strUtils
import qmp

var myRank* = 0
var nRanks* = 1

proc commsInit* =
  var argc {.importc:"cmdCount", global.}:cint
  var argv {.importc:"cmdLine", global.}:ptr cstring
  var prv = QMP_THREAD_FUNNELED
  #var prv = QMP_THREAD_SERIALIZED
  let err = QMP_init_msg_passing(argc.addr, argv.addr, prv, prv.addr)
  myRank = int(QMP_get_node_number())
  nRanks = int(QMP_get_number_of_nodes())
proc commsFinalize* =
  QMP_finalize_msg_passing()
proc commsAbort*(status = -1) =
  QMP_abort(status.cint)

proc evalArgs*(call:var NimNode; args:NimNode):NimNode =
  result = newStmtList()
  for i in 0..<args.len:
    let t = genSym()
    let a = args[i]
    result.add(quote do:
      let `t` = `a`
      )
    call.add(t)
proc cprintf*(fmt:cstring){.importc:"printf",varargs,header:"<stdio.h>".}
#proc printfOrdered(
macro printf*(fmt:string; args:varargs[untyped]):auto =
  var call = newCall(ident("cprintf"), fmt)
  result = evalArgs(call, args)
  result.add(quote do:
    if myRank==0 and threadNum==0:
      `call`
    )
proc echoRaw*(x: varargs[typed, `$`]) {.magic: "Echo".}
macro echoAll*(args:varargs[untyped]):auto =
  var call = newCall(bindSym"echoRaw")
  result = evalArgs(call, args)
  result.add(quote do:
    `call`
    )
macro echoRank*(args:varargs[untyped]):auto =
  var call = newCall(bindSym"echoRaw")
  call.add ident"myRank"
  call.add newLit"/"
  call.add ident"nRanks"
  call.add newLit": "
  result = evalArgs(call, args)
  template f(x:untyped):untyped =
    if threadNum==0: x
  result.add getAst(f(call))
macro echo0*(args: varargs[untyped]): auto =
  var call = newCall(bindSym"echoRaw")
  result = evalArgs(call, args)
  result.add(quote do:
    if myRank==0 and threadNum==0:
      `call`
    )
macro makeEchos(n:static[int]):auto =
  template ech(x,y: untyped): untyped =
    template echo*(): untyped =
      when nimvm:
        x
      else:
        y
  result = newStmtList()
  var er = newCall(bindSym"echoRaw")
  var e0 = newCall(bindSym"echo0")
  var ea = newSeq[NimNode](0)
  for i in 1..n:
    let ai = ident("a" & $i)
    er.add ai
    e0.add ai
    ea.add newNimNode(nnkIdentDefs).add(ai).add(ident"untyped").add(newEmptyNode())
    var t = getAst(ech(er,e0))
    #echo t.treerepr
    for j in 0..<i:
      t[0][3].add ea[j]
    result.add t
  #echoAll result.repr
makeEchos(10)

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

macro rankSumN*(a:varargs[typed]):auto =
  #echo "rankSum: ", a.repr
  #echo a.treeRepr
  var i0 = 0
  let t0 = a[0].getType
  if a.len==1:
    let a0 = a[0]
    result = quote do:
      if threadNum==0:
        qmpSum(`a0`)
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
    var s = newNimNode(nnkStmtList)
    let t = !"t"
    for i in 0..<a.len:
      let ai = a[i]
      let x = quote do:
        `ai` = `t`[`i`]
      s.add x[0]
    result = quote do:
      if threadNum==0:
        var `t` = `a`
        qmpSum(`t`)
        `s`
  else:
    result = newCall(!"rankSum")
    for i in 0..<a.len:
      if i==i0:
        let ai = unwrap(a[i])
        for j in 0..<ai.len:
          result.add(ai[j])
      else:
        result.add(a[i])
  #echo result.repr
macro rankSum*(a:varargs[untyped]):auto =
  if a.len==1:
    let a0 = a[0]
    result = quote do:
      if threadNum==0:
        qmpSum(`a0`)
  else:
    result = newCall(ident("rankSumN"))
    for v in a: result.add v

#var count = 0
template threadRankSum1*(a:untyped):untyped =
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
  var ta{.global.}:type(a)
  #var ta2{.global.}:array[512,type(a)]
  if threadNum==0:
    t0wait()
    for i in 1..<numThreads:
      a += cast[ptr type(a)](threadLocals.share[i].p)[]
      #a += ta2[threadNum]
    rankSum(a)
    ta = a
    twait0()
  else:
    threadLocals.share[threadNum].p = a.addr
    #ta2[threadNum] = a
    t0wait()
    twait0()
    a = ta

proc threadRankSumN*(a:NimNode):auto =
  echo a.treeRepr
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
  result.add(m)
  #echo result.treeRepr
macro threadRankSum*(a:varargs[untyped]):auto =
  if a.len==1:
    template trs1(x:untyped):untyped = threadRankSum1(x)
    result = getAst(trs1(a[0]))
  else:
    result = threadRankSumN(a)


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
