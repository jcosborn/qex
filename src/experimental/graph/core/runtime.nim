import base, traverse

proc updated*(x: Gvalue) =
  let grt = x.runtime
  inc grt.graphEpochCounter
  x.epoch = grt.graphEpochCounter

proc maxInputEpoch(inputs: openArray[Gvalue]): int =
  for input in inputs:
    if result < input.epoch:
      result = input.epoch

proc debugEval(node: Gvalue) =
  if not node.runtime.graphDebug:
    return
  var s = "[graph/core] eval: " & node.nodeRepr
  node.walkEvalDeps(proc(input: Gvalue) =
    s &= "\n  " & input.nodeRepr
  )
  echo s

proc evaluated*(x: Gvalue) =
  ## Marks a node current with respect to raw `inputs` only.
  ## This ignores custom eval deps and hidden deps; use it only when the
  ## caller owns the node's full forward freshness contract.
  x.epoch = x.inputs.maxInputEpoch

proc eval*(v: Gvalue): Gvalue {.discardable.} =
  var seen = initNodeSet()

  proc walkNode(node: Gvalue) =
    if node == nil:
      raiseError("eval traversal encountered nil node")
    if not seen.markSeenNode(node):
      return

    var maxep = 0
    node.walkEvalDeps(proc(input: Gvalue) =
      walkNode(input)
      if maxep < input.epoch:
        maxep = input.epoch)
    if node.epoch < maxep:
      let f = node.gfunc
      node.debugEval
      if f.forward != nil:
        f.forward node
        node.epoch = maxep
        # `f` is checked above; `node.runtime` is non-nil by construction.
        node.runtime.recordRun(f)
      else:
        raiseError(
          "inputs.len: " & $node.inputs.len &
          ", but no forward function defined for:\n" &
          node.nodeRepr)

  walkNode(v)
  v
