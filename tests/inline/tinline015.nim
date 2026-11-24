import base/[metaUtils, basicOps]
import macros
getOptimPragmas()

# inlineProcs copies the body of a proc into the code, except for procs
# marked for kernels (gpuInline, gpuCall), which stay calls; it stops on a
# proc whose copy would go wrong.

proc twice(x: float): float {.gpuInline.} = 2.0*x
proc thrice(x: float): float {.gpuCall.} = 3.0*x
proc plus1(x: float): float = x + 1.0
proc at(x: var array[3, float]; i: int): var float = x[i]
proc cadd(x, y: cint): cint =
  {.emit: "result = x + y;".}

macro calls(name: static string; body: typed): untyped =
  ## the calls of name in body, with body
  var n = 0
  proc count(x: NimNode) =
    if x.kind == nnkTypeOfExpr: return
    if x.kind in CallNodes and x[0].kind == nnkSym and $x[0] == name: inc n
    for c in x: count(c)
  when defined(dumpCalls): echo body.treeRepr
  count(body)
  result = newStmtList(body, newLit(n))

proc run =
  var r = 0.0
  let kept = calls("twice"):
    inlineProcs:
      r = thrice(twice(plus1(1.0)))
  echo "* marked procs stay calls: ", r, " ", kept
  doAssert r == 12.0 and kept == 1
  let copied = calls("plus1"):
    inlineProcs:
      r = plus1(r)
  echo "* other procs are copied: ", r, " ", copied
  doAssert r == 13.0 and copied == 0
  var a = [1.0, 2.0, 3.0]
  doAssert not compiles((inlineProcs: (at(a, 1) = 5.0)))
  doAssert not compiles((inlineProcs: (discard cadd(1, 2))))
  echo "* a var result and C text naming a parameter stop the compile"

run()
