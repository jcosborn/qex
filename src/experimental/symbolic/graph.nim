## lazy generic computational graph

#[

requires: `--multimethods:on`

We want the symbolic graph nodes to be type generic.  Therefore we
need dynamic dispatch based on object types at runtime.  This
implementation uses Nim's builtin multimethods for this purpose.

Nim's multimethods may slow down single method dispatch time, which
would affect the performance of the comms module.  We need to measure
it to understand the impact.

The functions, ident and add, are the same as a variable argument
function, in terms of symbolic formulae, but we treat functions
with different number of arguments differently, because the
implementations of the functions would be different.  It also avoids
increasing dynamic dispatch overhead.

Typed values are enclosed in derived types of `SymNodeValueConcrete`
and referenced from nodes, which have a single type.  Since Nim
doesn't allow mixing generic functions with methods, we need to
define a method for each and every combination of concrete types
we use, and wrap our existing generic functions in each method.

Both SymNodeValueConcrete and SymNode are references.  In the graph
we build, a shared node means the same variable.  We create new
nodes that equal to the existing nodes with the ident function to
create new referenced node objects, in order to avoid false sharing.
We use copySymNodeValue to create new referenced value objects,
such that different nodes refer to distinct value objects.

We use the tag sntVisited to avoid repeatedly traverse shared nodes.
The recursive graph traversal function all ends with Rec, just to
remind us to call `tagClearRec(z, sntVisited)` in the top level.

Further optimizations only possible after building all the graphs:
- Remove ident nodes
- Analyze and reuse allocations when possible

]#

#
# basic type support
#

type
  SymNodeTag = enum
    sntVisited, sntAssigned, sntNeedUpdate, sntNeedGradient, sntFixedGradient
    # sntReusable, ...
  SymNodeTags = set[SymNodeTag]
  SymNodeValue* = ref object of RootObj  ## Represent unallocated symbolic value
  SymNodeValueConcrete* = ref object of SymNodeValue  ## For any concrete values
  SymNodeGradient = object
    ## for a particular variable v
    dependent: SymNode  ## a variable v that depends on this node, x
    gradient: SymNode  ## dv/dx
  SymNode* = ref object
    # This object can be cyclic, because gradients refer to ancestor nodes
    value*: SymNodeValue
    inputs*: seq[SymNode]
    forward: proc(z: SymNode)  ## runs the actual compute
    arg: SymNodeValue  ## extra argument forward/backward uses
    runCount: int
    allocateValue: proc(z: SymNode)
    backward: proc(z: SymNode, i: int, dep: SymNode): SymNode  ## create graphs
    gradients: seq[SymNodeGradient]  ## saved gradient graphs
    name: string
    tag: SymNodeTags
    id: int
  SymNodeError = object of Defect
  SymNodeValueError = object of Defect

template raiseError*(msg: string) =
  raise newException(SymNodeError, msg)

template raiseValueError*(msg: string) =
  raise newException(SymNodeValueError, msg)

template raiseErrorBaseMethod*(msg: string) =
  raise newException(
    SymNodeError,
    "Base method invoked: " & msg &
      "\nMake sure to pass `--multimethods:on` and check there is a custom method for each derived type.")

method `$`*(v: SymNodeValue): string {.base.} = "SymNodeValue"

func `$`*(z: SymNode): string =
  z.name & "#" & $z.id

func nodeRepr*(z: SymNode): string
func treeRepr*(z: SymNode): string

method copySymNodeValue*(v: SymNodeValue): SymNodeValue {.base.} =
  ## nothing to copy
  v

method copySymNodeValue*(v: SymNodeValueConcrete): SymNodeValueConcrete =
  raiseValueError("Custom method required for concrete value: " & $v)

proc newSymNode*(
    value = SymNodeValue(),
    inputs: seq[SymNode] = @[],
    forward: proc(z: SymNode) = nil,
    arg: SymNodeValue = nil,
    runCount: int = 0,
    allocateValue: proc(z: SymNode) = nil,
    backward: proc(z: SymNode, i: int, dep: SymNode): SymNode = nil,
    gradients: seq[SymNodeGradient] = @[],
    name: string = "",
    tag: SymNodeTags = {}): SymNode =
  ## Create new SymNode with a unique id.
  var id {.global.} = 0
  result = SymNode(value: value, inputs: inputs, forward: forward, arg: arg, runCount: runCount,
    allocateValue: allocateValue, backward: backward, gradients: gradients, name: name, tag: tag, id: id)
  id.inc

proc copySymNode*(z: SymNode): SymNode =
  newSymNode(value = z.value.copySymNodeValue, inputs = z.inputs, forward = z.forward,
    arg = if z.arg != nil: z.arg.copySymNodeValue else: nil, runCount = z.runCount,
    allocateValue = z.allocateValue, backward = z.backward, gradients = z.gradients,
    name = z.name, tag = z.tag)

proc assignSymNode*(z: SymNode, x: SymNode) =
  z.value = x.value.copySymNodeValue
  z.inputs = x.inputs
  z.forward = x.forward
  if x.arg != nil:
    z.arg = x.arg.copySymNodeValue
  z.runCount = x.runCount
  z.allocateValue = x.allocateValue
  z.backward = x.backward
  z.gradients = x.gradients
  z.name = x.name
  z.tag = x.tag

proc gradientDependentOrNil*(z: SymNode, dep: SymNode): SymNode =
  ## May return nil if not exists.
  for g in z.gradients:
    if dep == g.dependent:
      return g.gradient
  # We don't have a matching dependent variable.
  return nil

proc gradientDependentAssign(z: SymNode, dep: SymNode, grad: SymNode) =
  ## Replace if exists, otherwise add to the list.
  for g in z.gradients.mitems:
    if dep == g.dependent:
      g.gradient = grad
      return
  z.gradients.add SymNodeGradient(dependent: dep, gradient: grad)

proc assign*(z: SymNode, v: SymNodeValueConcrete) =
  z.value = v
  z.tag.incl sntAssigned

#
# generic symbol support
#

proc newSym*(s: string): SymNode =
  newSymNode(name = s)

method identSymNodeValue*(z: SymNodeValue, x: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.repr & "\n  " & x.repr)

method identAllocateSymNodeValue*(z: SymNode, x: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.nodeRepr & "\n  " & x.repr)

method iaddSymNodeValue*(z: SymNodeValue, x: SymNodeValue, y: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.repr & "\n  " & x.repr & "\n  " & y.repr)

method iaddAllocateSymNodeValue*(z: SymNode, x: SymNodeValue, y: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.nodeRepr & "\n  " & x.repr & "\n  " & y.repr)

#
# float support
#

type SymNodeValueFloat* = ref object of SymNodeValueConcrete
  floatValue*: float

method `$`*(v: SymNodeValueFloat): string = $v.floatValue

proc assign*(z: SymNode, v: float) =
  z.assign SymNodeValueFloat(floatValue: v)

method copySymNodeValue*(v: SymNodeValueFloat): SymNodeValueFloat =
  SymNodeValueFloat(floatValue: v.floatValue)

method identSymNodeValue*(z: SymNodeValueFloat, x: SymNodeValueFloat) =
  z.floatValue = x.floatValue

method identAllocateSymNodeValue*(z: SymNode, x: SymNodeValueFloat) =
  z.value = SymNodeValueFloat()

method iaddSymNodeValue*(z: SymNodeValueFloat, x: SymNodeValueFloat, y: SymNodeValueFloat) =
  z.floatValue = x.floatValue + y.floatValue

method iaddAllocateSymNodeValue*(z: SymNode, x: SymNodeValueFloat, y: SymNodeValueFloat) =
  z.value = SymNodeValueFloat()

#
# minimum algebra for the nodes
#

proc ident*(x:SymNode): SymNode
proc add*(x: SymNode, y: SymNode): SymNode

proc identForward(z: SymNode) =
  identSymNodeValue(z.value, z.inputs[0].value)

proc identAllocate(z: SymNode) =
  identAllocateSymNodeValue(z, z.inputs[0].value)

proc identBackward(z: SymNode, i: int, dep: SymNode): SymNode =
  if i != 0:
    raiseError("ident has 1 operand, got i = " & $i)
  let g = z.gradientDependentOrNil dep
  if g == nil:
    return newSymNode(value = SymNodeValueFloat(floatValue: 1.0), name = "One[ident]")
  else:
    return g.ident

proc ident*(x:SymNode): SymNode =
  newSymNode(
    inputs = @[x],
    forward = identForward,
    allocateValue = identAllocate,
    backward = identBackward,
    name = "ident")

proc addForward(z: SymNode) =
  iaddSymNodeValue(z.value, z.inputs[0].value, z.inputs[1].value)

proc addAllocate(z: SymNode) =
  iaddAllocateSymNodeValue(z, z.inputs[0].value, z.inputs[1].value)

proc addBackward(z: SymNode, i: int, dep: SymNode): SymNode =
  if i != 0 and i != 1:
    raiseError("add has 2 operands, got i = " & $i)
  let g = z.gradientDependentOrNil dep
  if g == nil:
    return newSymNode(value = SymNodeValueFloat(floatValue: 1.0), name = "One[add]")
  else:
    return g.ident

proc add*(x: SymNode, y: SymNode): SymNode =
  newSymNode(
    inputs = @[x, y],
    forward = addForward,
    allocateValue = addAllocate,
    backward = addBackward,
    name = "add")

#
# graph traversal and evaluation
#

proc tagClearRec(z: SymNode, tag: SymNodeTag) =
  ## This does not use sntVisited, so it will repeat on shared nodes.
  if tag in z.tag:
    z.tag.excl tag
  for i in z.inputs:
    i.tagClearRec tag

proc allocateRec(z: SymNode) =
  if sntVisited notin z.tag:
    z.tag.incl sntVisited
    for i in z.inputs:
      i.allocateRec
    if not (z.value of SymNodeValueConcrete):
      if z.allocateValue == nil:
        raiseError("undefined allocateValue for node: " & z.nodeRepr)
      z.allocateValue z

proc allocate*(z: SymNode) =
  z.allocateRec
  z.tagClearRec sntVisited

proc tagUpdate*(z: SymNode) =
  ## call this for newly created inner nodes in *Backward
  z.tag.incl sntNeedUpdate

proc tagUpdateRec(z: SymNode) =
  if sntVisited in z.tag:
    return
  z.tag.incl sntVisited
  if sntAssigned in z.tag:
    z.tag.excl sntAssigned
  else:
    var needupdate = false
    for i in z.inputs:
      needupdate = needupdate or sntAssigned in i.tag
      i.tagUpdateRec
      needupdate = needupdate or sntNeedUpdate in i.tag
    if needupdate:
      z.tag.incl sntNeedUpdate

proc evalRec(z: SymNode) =
  if sntVisited in z.tag:
    if sntNeedUpdate in z.tag:
      raiseError "cycle detected"
  elif sntNeedUpdate in z.tag:
    z.tag.incl sntVisited
    for i in z.inputs:
      i.evalRec
    if z.forward != nil:
      z.forward z
      z.runCount.inc
    elif z.inputs.len > 0:
      raiseError("inputs.len: " & $z.inputs.len & ", but no forward function defined for:\n" & z.nodeRepr)
    z.tag.excl sntNeedUpdate

proc eval*(z: SymNode) =
  z.tagUpdateRec
  z.tagClearRec sntVisited
  z.evalRec
  z.tagClearRec sntVisited

proc tagUpdateNeedGradientRec(z: SymNode) =
  if sntVisited in z.tag:
    return
  z.tag.incl sntVisited
  var needgradient = false
  for i in z.inputs:
    i.tagUpdateNeedGradientRec
    needgradient = needgradient or sntNeedGradient in i.tag
  if needgradient and sntNeedGradient notin z.tag:
    z.tag.incl sntNeedGradient

proc gradientRec(z: SymNode, dep: SymNode) =
  ## gradient of dep with respect to z
  # We tag newly created nodes from z.backward(z, i, dep), with needUpdate.
  for i in 0..<z.inputs.len:
    let input = z.inputs[i]
    if sntNeedGradient in input.tag:
      if z.backward == nil:
        raiseError(z.nodeRepr & ":" & $i & ":" & input.nodeRepr & ": backward undefined")

      # Start making gradient.
      let grad = input.gradientDependentOrNil dep

      if sntVisited in input.tag:
        if grad != nil:
          if sntFixedGradient notin input.tag:
            let childGrad = z.backward(z, i, dep)
            # We always request update for the newly created gradient node.
            childGrad.tag.incl sntNeedUpdate
            # We need to combine the gradient without breaking the existing graph.
            # Because the previous built graph may have a reference of this node, `grad`,
            # our new node has to reuse `grad`. We use a copy of `grad` and assign back.
            grad.assignSymNode(grad.copySymNode.add childGrad)
        else:
          raiseError(z.nodeRepr & ":" & $i & ":" & input.nodeRepr & ": visited but no gradient")
      else:
        # Not visited this time.  Construct gradient if needed.
        input.tag.incl sntVisited
        if grad == nil:
          let childGrad = z.backward(z, i, dep)
          # We always request update for the newly created gradient node.
          childGrad.tag.incl sntNeedUpdate
          gradientDependentAssign(input, dep, childGrad)
        else:
          # Existent gradient means it was setup previously outside of this gradientRec(z, dep),
          # by gradient calls of the same dep with respective to other variables,
          # and we trust the construction.
          input.tag.incl sntFixedGradient
        input.gradientRec dep

proc gradient*(dep: SymNode, x: SymNode): SymNode =
  if sntNeedGradient notin x.tag:
    x.tag.incl sntNeedGradient
  dep.tagUpdateNeedGradientRec
  dep.tagClearRec sntVisited
  dep.gradientRec dep
  dep.tagClearRec sntVisited
  dep.tagClearRec sntNeedGradient
  dep.tagClearRec sntFixedGradient
  x.gradientDependentOrNil dep

proc optimize*(output: seq[SymNode], variables: seq[SymNode]) =
  # TODO
  discard

proc sharedNodesRec(shared: var seq[SymNode], z: SymNode) =
  if sntVisited in z.tag:
    var found = false
    for n in shared:
      if n == z:
        found = true
        break
    if not found:
      shared.add z
  else:
    z.tag.incl sntVisited
    for i in z.inputs:
      shared.sharedNodesRec i

#
# to string
#

func nodeRepr*(z: SymNode): string =
  result = $z & $z.tag & ": " & $z.value & ", run: " & $z.runCount
  if z.arg != nil:
    result &= ", arg: " & $z.arg
  if z.inputs.len > 0:
    result &= ", inputs: {"
    for i in 0..<z.inputs.len:
      if i > 0:
        result &= ", "
      result &= "[" & $i & "]: " & $z.inputs[i]
    result &= "}"
  if z.gradients.len > 0:
    result &= ", gradients: {"
    for i in 0..<z.gradients.len:
      if i > 0:
        result &= ", "
      result &= "[" & $i & "]: " & $z.gradients[i].dependent & " -> " & $z.gradients[i].gradient
    result &= "}"
  if z.value != nil:
    result &= ": " & $z.value

func toStringRec(z: SymNode, pre: string, shared: seq[SymNode]): string =
  result = pre & z.nodeRepr
  if sntVisited in z.tag:
    result &= " [shared]"
  else:
    z.tag.incl sntVisited
    for zid in 0..<shared.len:
      if z == shared[zid]:
        result &= " (shared)"
        break
    for i in z.inputs:
      result &= "\n" & i.toStringRec(pre & "  ", shared)

func treeRepr*(z: SymNode): string =
  var shared = newseq[SymNode]()
  shared.sharedNodesRec z
  z.tagClearRec sntVisited
  result = z.toStringRec("", shared)
  z.tagClearRec sntVisited

#
# example
#

when isMainModule:
  let x = newSym("x")
  echo x.nodeRepr
  let y = newSym("y")
  echo y.nodeRepr
  let z = x.add y
  z.name = "z"
  echo z.nodeRepr
  x.assign 0.1
  y.assign 0.2
  z.allocate
  echo "\nafter assign & allocate\n", z.nodeRepr
  z.eval
  echo "\nafter z.eval\n", z.nodeRepr
  echo "z.value = ", z.value
  let dzdx = z.gradient x
  dzdx.name = "dzdx"
  echo "\n", dzdx.nodeRepr
  dzdx.eval
  echo "\nafter dzdx.eval:\n", dzdx.nodeRepr
  echo "dzdx.value = ", dzdx.value
  let w = z.add z.add x
  w.name = "w"
  echo "\nw = z.add z.add x:\n", w.treeRepr
  let dwdx = w.gradient x
  dwdx.name = "dwdx"
  echo "\nafter w.gradient x:\n", w.treeRepr
  echo "\ndwdx.treeRepr:\n", dwdx.treeRepr
  dwdx.allocate
  echo "\nafter dwdx.allocate:\n", dwdx.treeRepr
  dwdx.eval
  echo "\nafter dwdx.eval:\n", dwdx.treeRepr
  echo "dwdx.value = ", dwdx.value
  let dwdy = w.gradient y
  dwdy.name = "dwdy"
  echo "\nafter w.gradient y:\n", w.treeRepr
  echo "\n", dwdy.treeRepr
  dwdy.allocate
  dwdy.eval
  echo "dwdy.value = ", dwdy.value
  echo "\n## Current graphs"
  echo "\n", z.treeRepr
  echo "\n", w.treeRepr
  echo "\n", dzdx.treeRepr
  echo "\n", dwdx.treeRepr
  echo "\n", dwdy.treeRepr
  x.assign 0.3
  echo "\n## x.assign 0.3\n"
  dzdx.eval
  echo dzdx.treeRepr
  dwdx.eval
  echo dwdx.treeRepr
  w.allocate
  w.eval
  echo w.treeRepr
