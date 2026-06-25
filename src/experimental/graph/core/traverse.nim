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
