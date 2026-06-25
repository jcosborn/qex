import std/tables
import std/sets
import base, traverse

proc updated*(x: Gvalue) =
  x.staticZeroLeaf = false
  let grt = x.runtime
  inc grt.graphEpochCounter
  x.epoch = grt.graphEpochCounter

proc debugEval(node: Gvalue) =
  if not node.runtime.graphDebug:
    return
  var s = "[graph/core] eval: " & node.nodeRepr
  node.walkInputView(iwmEval, proc(input: Gvalue) =
    s &= "\n  " & input.nodeRepr
  )
  echo s

proc eval*[T: Gvalue](v: T): T {.discardable.} =
  var seen = initHashSet[NodeKey]()
  var active = initHashSet[NodeKey]()

  proc walkNode(node: Gvalue) =
    let key = node.nodeKey
    if key in active:
      raiseError("cycle detected while evaluating graph:\n" & node.nodeRepr)
    if not seen.markSeenNode(node):
      return
    active.incl key
    defer:
      active.excl key

    var maxep = 0
    node.walkInputView(iwmEval, proc(input: Gvalue) =
      walkNode(input)
      if maxep < input.epoch:
        maxep = input.epoch)
    if node.epoch < maxep:
      let f = node.gfunc
      node.debugEval
      let forward =
        if f == nil:
          nil
        else:
          f.forward
      if forward != nil:
        forward node
        if node.epoch < maxep:
          node.epoch = maxep
        inc node.runtime.runCountsByNode.mgetOrPut(node.stableNodeId, 0)
      else:
        raiseError(
          "inputs.len: " & $node.inputs.len &
          ", but no forward function defined for:\n" &
          node.nodeRepr)

  walkNode(v)
  v
