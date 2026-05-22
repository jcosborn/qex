import base

proc walkRawInputs*(v: Gvalue, visit: GnodeVisit) =
  for j in 0..<v.inputs.len:
    visit v.inputs[j]

proc walkDeps*(v: Gvalue,
               mode: InputWalkMode,
               visit: GnodeVisit) =
  proc checkedVisit(child: Gvalue) =
    if child == nil:
      raiseValueError("dependency walk produced nil dependency" & v.nodeContext)
    visit child

  let walk = v.gfunc.depWalkForMode(mode)
  if walk == nil:
    v.walkRawInputs(checkedVisit)
  else:
    walk(v, checkedVisit)

  v.walkHiddenDeps(mode, checkedVisit)

proc collectNodeInputs*(v: Gvalue,
                        mode: InputWalkMode): seq[Gvalue] =
  var children: seq[Gvalue] = @[]
  v.walkDeps(mode, proc(child: Gvalue) =
    children.add child)
  result = children

proc treeRepr*(v: Gvalue, mode: InputWalkMode): string =
  ## Renders the graph through the same dependency surface used by `mode`.
  treeReprImpl(v, proc(node: Gvalue): seq[Gvalue] = node.collectNodeInputs(mode))

proc sameNode*(a: Gvalue, b: Gvalue): bool =
  a.nodeKey == b.nodeKey
