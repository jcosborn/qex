from strutils import join
import std/sets
import base

proc walkInputView*(v: Gvalue,
               mode: InputWalkMode,
               visit: GnodeVisit) =
  proc checkedVisit(child: Gvalue) =
    if child == nil:
      raiseValueError("input view produced nil dependency:\n" & v.nodeRepr)
    if child.runtime != v.runtime:
      raiseValueError("input view dependency mixes graph runtimes:\n" & v.nodeRepr)
    discard child.stableNodeId
    visit child

  let walk =
    if v.gfunc == nil:
      nil
    else:
      v.gfunc.inputView
  if walk == nil:
    for input in v.inputs:
      visit input
  else:
    walk(v, mode, checkedVisit)

proc collectInputView*(v: Gvalue,
                        mode: InputWalkMode): seq[Gvalue] =
  var children: seq[Gvalue] = @[]
  v.walkInputView(mode, proc(child: Gvalue) =
    children.add child)
  result = children

proc reaches*(v, target: Gvalue,
              mode: InputWalkMode,
              stop: Gvalue = nil): bool =
  ## Whether `target` is on `v`'s `mode` dependency surface (`v` counts), cut at `stop`.
  var seen = initHashSet[NodeKey]()
  var found = false
  proc visit(node: Gvalue) =
    if found or (stop != nil and node.nodeKey == stop.nodeKey):
      return
    if node.nodeKey == target.nodeKey:
      found = true
      return
    if not seen.markSeenNode(node):
      return
    node.walkInputView(mode, visit)
  visit(v)
  found

proc treeRepr*(v: Gvalue, mode: InputWalkMode): string =
  ## Renders the graph through the same input view used by `mode`.
  var shared = newseq[Gvalue]()
  var rendered = initHashSet[NodeKey]()

  proc render(x: Gvalue): seq[string] =
    let sharedIndex = shared.find x
    result = @[x.nodeRepr]
    if x.nodeKey in rendered:
      if sharedIndex >= 0:
        result[0] &= " #" & $sharedIndex
      return
    if sharedIndex >= 0:
      result[0] &= " #" & $sharedIndex & "#"
    rendered.incl x.nodeKey
    for child in x.collectInputView(mode):
      for renderedChild in render(child):
        result.add("  " & renderedChild)

  var seen = initHashSet[NodeKey]()
  proc collectShared(node: Gvalue) =
    if not seen.markSeenNode(node):
      if shared.find(node) < 0:
        shared.add node
      return
    for input in node.collectInputView(mode):
      collectShared(input)

  collectShared(v)

  result = render(v).join "\n"

proc treeRepr*(v: Gvalue): string =
  v.treeRepr(iwmReachable)
