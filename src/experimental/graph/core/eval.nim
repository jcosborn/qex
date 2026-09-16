import std/tables
import std/sets
import base, traverse
import base/profile  # qex timing facility: getTics/ticDiffSecs

proc updated*(x: Gvalue) =
  let grt = x.runtime
  if grt.evalFrame == nil and (x.staticZeroLeaf or x.restoreValue != nil or
      (x.gfunc != nil and not x.valueOverride)):
    inc grt.boundaryRevision
  x.staticZeroLeaf = false
  x.restoreValue = nil
  x.valueReady = true
  x.valueOverride = true
  x.stale = false
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
  let grt = v.runtime
  if grt.evalFrame != nil:
    grt.evalFrame(v)
    return v
  var seen = initHashSet[NodeKey]()
  var active = initHashSet[NodeKey]()

  proc walkNode(node: Gvalue) =
    let key = node.nodeKey
    if key in active:
      raiseError("cycle detected while evaluating graph:\n" & node.nodeRepr)
    if node.stale:
      raiseError("published plan result is invalid; evaluate its plan successfully and reacquire the result:\n" & node.nodeRepr)
    if key in seen and node.valueReady and node.hasStorage:
      return
    seen.incl key
    active.incl key
    defer:
      active.excl key

    var maxep = 0
    var hasInputs = false
    node.walkInputView(iwmEval, proc(input: Gvalue) =
      hasInputs = true
      walkNode(input)
      if maxep < input.epoch:
        maxep = input.epoch)
    if not node.valueReady or not node.hasStorage or node.epoch < maxep:
      let wasOverride = node.valueOverride
      node.valueReady = false
      var ready = false
      defer:
        node.valueReady = ready
      node.ensureStorage
      let f = node.gfunc
      node.debugEval
      let forward = if f == nil: nil else: f.forward
      if forward != nil:
        if f.bufferMode == bmZero and node.bufferBytes > 0:
          node.clearBuffer
        let t0 = getTics()
        forward node
        let secs = ticDiffSecs(getTics(), t0)
        if wasOverride:
          # Replacing an override changes the logical value, even when every
          # input is older. Ordinary storage restoration keeps its prior epoch.
          inc grt.boundaryRevision
          inc grt.graphEpochCounter
          node.epoch = grt.graphEpochCounter
        node.valueOverride = false
        grt.runStatsByNode.mgetOrPut(node.stableNodeId, RunStat())
          .record(secs, f.name)
      elif hasInputs:
        raiseError("inputs.len: " & $node.inputs.len &
          ", but no forward function defined for:\n" & node.nodeRepr)
      if node.epoch < maxep:
        node.epoch = maxep
      ready = true

  walkNode(v)
  v
