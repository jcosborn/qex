import base

proc walkRawInputs*(v: Gvalue, visit: GnodeVisit) =
  for j in 0..<v.inputs.len:
    let input = v.inputs[j]
    if input == nil:
      raiseError("node has nil input at index " & $j & ":\n" & v.nodeRepr)
    visit input

proc modeLabel(mode: InputWalkMode): string =
  case mode
  of iwmEval:
    "eval"
  of iwmGradSignature:
    "signature"
  of iwmDepend:
    "dependency"

proc walkDepsForMode(v: Gvalue,
                     mode: InputWalkMode,
                     visit: GnodeVisit) =
  if v == nil:
    raiseError(mode.modeLabel & " input walk received nil node")

  let spec = v.gfunc.depWalkForMode(mode)

  case spec.kind
  of dwkRawInputs:
    v.walkRawInputs(visit)
  of dwkNoInputs:
    discard
  of dwkCustomInputs:
    discard v.runDepWalkInputs(spec, visit)

  v.walkHiddenDeps(mode, visit)

proc walkEvalDeps*(v: Gvalue, visit: GnodeVisit) =
  v.walkDepsForMode(iwmEval, visit)

proc walkGradSignatureDeps*(v: Gvalue, visit: GnodeVisit) =
  v.walkDepsForMode(iwmGradSignature, visit)

proc walkDependDeps*(v: Gvalue, visit: GnodeVisit) =
  v.walkDepsForMode(iwmDepend, visit)

proc walkDeps*(v: Gvalue,
               mode: InputWalkMode,
               visit: GnodeVisit) =
  case mode
  of iwmEval:
    v.walkEvalDeps(visit)
  of iwmGradSignature:
    v.walkGradSignatureDeps(visit)
  of iwmDepend:
    v.walkDependDeps(visit)

proc collectNodeInputs*(v: Gvalue,
                        mode: InputWalkMode): seq[Gvalue] =
  if v == nil:
    return @[]
  var children: seq[Gvalue] = @[]
  v.walkDeps(mode, proc(child: Gvalue) =
    children.add child)
  result = children

proc treeRepr*(v: Gvalue, mode: InputWalkMode): string =
  ## Renders the graph through the same dependency surface used by `mode`.
  treeReprImpl(v, proc(node: Gvalue): seq[Gvalue] = node.collectNodeInputs(mode))

proc sameNode*(a: Gvalue, b: Gvalue): bool =
  a.nodeKey == b.nodeKey

proc walkDeferredEvalGraph*(root: Gvalue,
                            isDeferredNode: proc(v: Gvalue): bool {.closure.},
                            visit: GnodeVisit = nil) =
  ## Eval traversal that skips through deferred nodes by following their
  ## dependency view until a concrete value is needed. Non-deferred nodes also
  ## expose depend-mode hidden deps so closure boundaries remain reachable.
  var seen = initNodeSet()

  proc walkNode(node: Gvalue) =
    if node == nil:
      raiseError("eval traversal encountered nil node")
    if not seen.markSeenNode(node):
      return
    if isDeferredNode != nil and isDeferredNode(node):
      node.walkDependDeps(walkNode)
      return
    node.walkEvalDeps(walkNode)
    node.walkHiddenDeps(iwmDepend, walkNode)
    if visit != nil:
      visit node

  walkNode(root)
