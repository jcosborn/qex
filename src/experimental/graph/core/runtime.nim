import std/tables
import std/sets
import base, traverse

proc updated*(x: Gvalue) =
  let grt = x.runtime
  inc grt.graphEpochCounter
  x.epoch = grt.graphEpochCounter

proc debugEval(node: Gvalue) =
  if not node.runtime.graphDebug:
    return
  var s = "[graph/core] eval: " & node.nodeRepr
  node.walkDeps(iwmEval, proc(input: Gvalue) =
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
    node.walkDeps(iwmEval, proc(input: Gvalue) =
      walkNode(input)
      if maxep < input.epoch:
        maxep = input.epoch)
    if node.epoch < maxep:
      let f = node.gfunc
      node.debugEval
      let forward = f.forward
      if forward != nil:
        forward node
        if node.epoch < maxep:
          node.epoch = maxep
        let nodeId = node.stableNodeId
        node.runtime.runCountsByNode[nodeId] =
          node.runtime.runCountsByNode.getOrDefault(nodeId) + 1
      else:
        raiseError(
          "inputs.len: " & $node.inputs.len &
          ", but no forward function defined for:\n" &
          node.nodeRepr)

  walkNode(v)
  v
