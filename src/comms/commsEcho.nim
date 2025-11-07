import macros
import base/[threading,metaUtils]
import commsTypes

proc evalArgs*(call:var NimNode; args:NimNode):NimNode =
  result = newStmtList()
  for i in 0..<args.len:
    let t = genSym()
    let a = args[i]
    result.add(quote do:
      when `a` is openarray:
        let `t` = $`a`
      else:
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
macro echo0*(args: varargs[untyped]): untyped =
  var call = newCall(bindSym"echoRaw")
  result = evalArgs(call, args)
  result.add(quote do:
    bind myRank
    if myRank==0 and threadNum==0:
      `call`
    )
  #echo result.repr
macro makeEchos(n:static[int]): untyped =
  template ech(x,y: untyped) =
    template echo* =
      when nimvm:
        x
      else:
        y
  result = newStmtList()
  for i in 1..n:
    var er = newCall(bindSym"echoRaw")
    var e0 = newCall(bindSym"echo0")
    var ea = newSeq[NimNode](0)
    for j in 1..i:
      let ai = ident("a" & $j)
      er.add ai
      e0.add ai
      ea.add newNimNode(nnkIdentDefs).add(ai).add(ident"untyped").add(newEmptyNode())
    var t = getAst(ech(er,e0)).peelStmt
    #echo t.treerepr
    for j in 0..<i: t[3].add ea[j]
    result.add t
  #echo result.treerepr
makeEchos(64)
#[
template echo*(a1: untyped) =
  when nimvm:
    echoRaw(a1)
  else:
    echo0(a1)
template echo*(a1,a2: untyped) =
  when nimvm:
    echoRaw(a1,a2)
  else:
    echo0(a1,a2)
]#
