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
create new referenced node objects, in order to avoid false sharing.  TODO: when do we really need ident?
We use copySymNodeValue to create new referenced value objects,
such that different nodes refer to distinct value objects.

We use the tag sntVisited to avoid repeatedly traverse shared nodes.
The recursive graph traversal function all ends with Rec, just to
remind us to call `tagClearRec(z, sntVisited)` in the top level.

Further optimizations only possible after building all the graphs:
- Remove ident nodes
- Analyze and reuse allocations when possible

TEMP NOTES

Make backward dispatch based on types
- all nodes have to already have the value type
- only need to call the method using the values, which allow dispatch at runtime
- need to create typed values when creating the node
- so the node creating functions must dispatch based on the value types
- everything needs to be dispatched with methods
- should I just make graph nodes with type instead?

]#

from math import exp

#
# basic type support
#

type
  SymNodeTag = enum
    sntVisited, sntNeedGradient, sntFixedGradient
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
    arg: SymNodeValue  ## extra argument forward/backward uses, must be immutable and can be shared, use getArg/setArg
    runCount: int
    epoch: int  ## for resolving dependency, tracks update
    allocateValue: proc(z: SymNode)
    backward: proc(z: SymNode, i: int, dep: SymNode): SymNode  ## create graphs
    gradients: seq[SymNodeGradient]  ## saved gradient graphs
    name: string
    tag: SymNodeTags
    id: int
    refCount: int
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

method isSymNodeValueConcrete*(v: SymNodeValue): bool {.base.} = false
method isSymNodeValueConcrete*(v: SymNodeValueConcrete): bool = true

func getArg*(z: SymNode): SymNodeValue = z.arg
proc setArg*(z: SymNode, v: SymNodeValue) =
  if z.arg != nil and z.arg.isSymNodeValueConcrete:
    raiseValueError("Cannot set z.arg, which is a concrete value: " & $z.arg)
  else:
    z.arg = v

func nodeRepr*(z: SymNode): string
func treeRepr*(z: SymNode): string

method copySymNodeValue*(v: SymNodeValue): SymNodeValue {.base.} =
  ## nothing to copy
  v

method copySymNodeValue*(v: SymNodeValueConcrete): SymNodeValue =
  raiseValueError("Custom method required for concrete value: " & $v)

proc newSymNode*(
    value = SymNodeValue(),
    inputs: seq[SymNode] = @[],
    forward: proc(z: SymNode) = nil,
    arg: SymNodeValue = nil,
    runCount: int = 0,
    epoch: int = 0,
    allocateValue: proc(z: SymNode) = nil,
    backward: proc(z: SymNode, i: int, dep: SymNode): SymNode = nil,
    gradients: seq[SymNodeGradient] = @[],
    name: string = "",
    refCount: int = 0,
    tag: SymNodeTags = {}): SymNode =
  ## Create new SymNode with a unique id.
  var id {.global.} = 0
  result = SymNode(value: value, inputs: inputs, forward: forward, arg: arg, runCount: runCount,
    epoch: epoch, allocateValue: allocateValue, backward: backward, gradients: gradients,
    name: name, tag: tag, id: id, refCount: refCount)
  id.inc
  for i in result.inputs:
    i.refCount.inc

proc copySymNode*(z: SymNode): SymNode =
  newSymNode(value = z.value.copySymNodeValue, inputs = z.inputs, forward = z.forward,
    arg = z.arg, runCount = z.runCount, epoch = z.epoch,
    allocateValue = z.allocateValue, backward = z.backward, gradients = z.gradients,
    name = z.name, tag = z.tag, refCount = z.refCount)

proc assignSymNode*(z: SymNode, x: SymNode) =
  z.value = x.value.copySymNodeValue
  z.inputs = x.inputs
  z.forward = x.forward
  z.arg = x.arg
  z.runCount = x.runCount
  z.epoch = x.epoch
  z.allocateValue = x.allocateValue
  z.backward = x.backward
  z.gradients = x.gradients
  z.name = x.name
  z.tag = x.tag
  z.refCount = x.refCount

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

proc updated*(z: SymNode) =
  ## Tag this SymNode after updating its value.
  ## Graph dependency resolution depends on this.
  var epoch {.global.} = 0
  epoch.inc
  z.epoch = epoch

proc assign*(z: SymNode, v: SymNodeValueConcrete) =
  z.value = v
  z.updated

#
# generic symbol support
#

proc newSym*(s: string): SymNode =
  newSymNode(name = s)

method identAllocateSymNodeValue*(x: SymNodeValue): SymNodeValue {.base.} =
  raiseErrorBaseMethod("args:\n  " & x.repr)

method identSymNodeValue*(z: SymNodeValue, x: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.repr & "\n  " & x.repr)

method zeroAllocateSymNodeValue*(x: SymNodeValue): SymNodeValue {.base.} =
  raiseErrorBaseMethod("args:\n  " & x.repr)

method zeroSymNodeValue*(z: SymNodeValue, x: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.repr & "\n  " & x.repr)

method indexAllocateSymNodeValue*(x: SymNodeValue, i: SymNodeValue): SymNodeValue {.base.} =
  raiseErrorBaseMethod("args:\n  " & x.repr & "\n  " & i.repr)

method indexSymNodeValue*(z: SymNodeValue, x: SymNodeValue, i: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.repr & "\n  " & x.repr & "\n  " & i.repr)

method indexUpdateAllocateSymNodeValue*(x: SymNodeValue, i: SymNodeValue, y: SymNodeValue): SymNodeValue {.base.} =
  raiseErrorBaseMethod("args:\n  " & x.repr & "\n  " & i.repr & "\n  " & y.repr)

method indexUpdateSymNodeValue*(z: SymNodeValue, x: SymNodeValue, i: SymNodeValue, y: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.repr & "\n  " & x.repr & "\n  " & i.repr & "\n  " & y.repr)

method negAllocateSymNodeValue*(x: SymNodeValue): SymNodeValue {.base.} =
  raiseErrorBaseMethod("args:\n  " & x.repr)

method negSymNodeValue*(z: SymNodeValue, x: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.repr & "\n  " & x.repr)

method invAllocateSymNodeValue*(x: SymNodeValue): SymNodeValue {.base.} =
  raiseErrorBaseMethod("args:\n  " & x.repr)

method invSymNodeValue*(z: SymNodeValue, x: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.repr & "\n  " & x.repr)

method addAllocateSymNodeValue*(x: SymNodeValue, y: SymNodeValue): SymNodeValue {.base.} =
  raiseErrorBaseMethod("args:\n  " & x.repr & "\n  " & y.repr)

method addSymNodeValue*(z: SymNodeValue, x: SymNodeValue, y: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.repr & "\n  " & x.repr & "\n  " & y.repr)

method mulAllocateSymNodeValue*(x: SymNodeValue, y: SymNodeValue): SymNodeValue {.base.} =
  raiseErrorBaseMethod("args:\n  " & x.repr & "\n  " & y.repr)

method mulSymNodeValue*(z: SymNodeValue, x: SymNodeValue, y: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.repr & "\n  " & x.repr & "\n  " & y.repr)

method subAllocateSymNodeValue*(x: SymNodeValue, y: SymNodeValue): SymNodeValue {.base.} =
  raiseErrorBaseMethod("args:\n  " & x.repr & "\n  " & y.repr)

method subSymNodeValue*(z: SymNodeValue, x: SymNodeValue, y: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.repr & "\n  " & x.repr & "\n  " & y.repr)

method divideAllocateSymNodeValue*(x: SymNodeValue, y: SymNodeValue): SymNodeValue {.base.} =
  raiseErrorBaseMethod("args:\n  " & x.repr & "\n  " & y.repr)

method divideSymNodeValue*(z: SymNodeValue, x: SymNodeValue, y: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.repr & "\n  " & x.repr & "\n  " & y.repr)

method expAllocateSymNodeValue*(x: SymNodeValue): SymNodeValue {.base.} =
  raiseErrorBaseMethod("args:\n  " & x.repr)

method expSymNodeValue*(z: SymNodeValue, x: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.repr & "\n  " & x.repr)

#
# support of generic multiple values per node
#

type SymNodeValueTuple* = ref object of SymNodeValue
  tupleValue*: seq[SymNodeValue]

method getTuple*(v: SymNodeValue): seq[SymNodeValue] {.base.} =
  raiseValueError("Custom method required for value: " & $v)

method getTuple*(v:SymNodeValueTuple): seq[SymNodeValue] = v.tupleValue

method `$`*(v: SymNodeValueTuple): string = $v.tupleValue

type SymNodeValueInt* = ref object of SymNodeValue
  intValue*: int

method getInt*(v: SymNodeValue): int {.base.} =
  raiseValueError("Custom method required for value: " & $v)

method getInt*(v:SymNodeValueInt): int = v.intValue

method `$`*(v: SymNodeValueInt): string = $v.intValue

method identAllocateSymNodeValue*(x: SymNodeValueTuple): SymNodeValue =
  var vs = newseq[SymNodeValue](x.tupleValue.len)
  for i in 0..<vs.len:
    vs[i] = identAllocateSymNodeValue(x.tupleValue[i])
  SymNodeValueTuple(tupleValue: vs)

method identSymNodeValue*(z: SymNodeValueTuple, x: SymNodeValueTuple) =
  for i in 0..<x.tupleValue.len:
     identSymNodeValue(z.tupleValue[i], x.tupleValue[i])

method zeroAllocateSymNodeValue*(x: SymNodeValueTuple): SymNodeValue =
  identAllocateSymNodeValue(x)

method zeroSymNodeValue*(z: SymNodeValueTuple, x: SymNodeValueTuple) =
  for i in 0..<x.tupleValue.len:
     zeroSymNodeValue(z.tupleValue[i], x.tupleValue[i])

method indexAllocateSymNodeValue*(x: SymNodeValueTuple, i: SymNodeValueInt): SymNodeValue =
  # reuse ident with the indexed value
  identAllocateSymNodeValue(x.tupleValue[i.intValue])

method indexSymNodeValue*(z: SymNodeValue, x: SymNodeValueTuple, i: SymNodeValueInt) =
  # reuse ident with the indexed value
  identSymNodeValue(z, x.tupleValue[i.intValue])

method indexUpdateAllocateSymNodeValue*(x: SymNodeValueTuple, i: SymNodeValueInt, y: SymNodeValue): SymNodeValue =
  # reuse ident with the original value
  identAllocateSymNodeValue(x)

method indexUpdateSymNodeValue*(z: SymNodeValueTuple, x: SymNodeValueTuple, i: SymNodeValueInt, y: SymNodeValue) =
  for k in 0..<x.tupleValue.len:
    if k == i.intValue:
      identSymNodeValue(z.tupleValue[k], y)
    else:
      identSymNodeValue(z.tupleValue[k], x.tupleValue[k])

#
# float support
#

type SymNodeValueFloat* = ref object of SymNodeValueConcrete
  floatValue*: float

proc newSymNodeFloat*(floatValue: float, name = ""): SymNode =
  result = newSymNode(value = SymNodeValueFloat(floatValue: floatValue), name = name)
  result.updated

converter toSymNode*(x: float): SymNode =
  newSymNodeFloat(x, "autoConv")

method getFloat*(v: SymNodeValue): float {.base.} =
  raiseValueError("Custom method required for value: " & $v)

method getFloat*(v:SymNodeValueFloat): float = v.floatValue

method `$`*(v: SymNodeValueFloat): string = $v.floatValue

converter toSymNodeValueFloat*(x: float): SymNodeValueFloat = SymNodeValueFloat(floatValue: x)
proc assign*(z: SymNode, v: float) =
  z.assign SymNodeValueFloat(floatValue: v)

method copySymNodeValue*(v: SymNodeValueFloat): SymNodeValue =
  SymNodeValueFloat(floatValue: v.floatValue)

method identAllocateSymNodeValue*(x: SymNodeValueFloat): SymNodeValue =
  SymNodeValueFloat()

method identSymNodeValue*(z: SymNodeValueFloat, x: SymNodeValueFloat) =
  z.floatValue = x.floatValue

method zeroAllocateSymNodeValue*(x: SymNodeValueFloat): SymNodeValue =
  SymNodeValueFloat()

method zeroSymNodeValue*(z: SymNodeValueFloat, x: SymNodeValueFloat) =
  z.floatValue = 0.0

method negAllocateSymNodeValue*(x: SymNodeValueFloat): SymNodeValue =
  SymNodeValueFloat()

method negSymNodeValue*(z: SymNodeValueFloat, x: SymNodeValueFloat) =
  z.floatValue = -x.floatValue

method invAllocateSymNodeValue*(x: SymNodeValueFloat): SymNodeValue =
  SymNodeValueFloat()

method invSymNodeValue*(z: SymNodeValueFloat, x: SymNodeValueFloat) =
  z.floatValue = 1.0/x.floatValue

method addAllocateSymNodeValue*(x: SymNodeValueFloat, y: SymNodeValueFloat): SymNodeValue =
  SymNodeValueFloat()

method addSymNodeValue*(z: SymNodeValueFloat, x: SymNodeValueFloat, y: SymNodeValueFloat) =
  z.floatValue = x.floatValue + y.floatValue

method mulAllocateSymNodeValue*(x: SymNodeValueFloat, y: SymNodeValueFloat): SymNodeValue =
  SymNodeValueFloat()

method mulSymNodeValue*(z: SymNodeValueFloat, x: SymNodeValueFloat, y: SymNodeValueFloat) =
  z.floatValue = x.floatValue * y.floatValue

method subAllocateSymNodeValue*(x: SymNodeValueFloat, y: SymNodeValueFloat): SymNodeValue =
  SymNodeValueFloat()

method subSymNodeValue*(z: SymNodeValueFloat, x: SymNodeValueFloat, y: SymNodeValueFloat) =
  z.floatValue = x.floatValue - y.floatValue

method divideAllocateSymNodeValue*(x: SymNodeValueFloat, y: SymNodeValueFloat): SymNodeValue =
  SymNodeValueFloat()

method divideSymNodeValue*(z: SymNodeValueFloat, x: SymNodeValueFloat, y: SymNodeValueFloat) =
  z.floatValue = x.floatValue / y.floatValue

method expAllocateSymNodeValue*(x: SymNodeValueFloat): SymNodeValue =
  SymNodeValueFloat()

method expSymNodeValue*(z: SymNodeValueFloat, x: SymNodeValueFloat) =
  z.floatValue = exp(x.floatValue)

#
# generic algebra
#

proc ident*(x: SymNode): SymNode
proc zero*(x: SymNode): SymNode
proc index*(x: SymNode, i: SymNodeValue): SymNode
proc indexUpdate*(x: SymNode, i: SymNodeValue, y: SymNode): SymNode
proc neg*(x: SymNode): SymNode
proc inv*(x: SymNode): SymNode
proc add*(x: SymNode, y: SymNode): SymNode
proc mul*(x: SymNode, y: SymNode): SymNode
proc sub*(x: SymNode, y: SymNode): SymNode
proc divide*(x: SymNode, y: SymNode): SymNode
proc exp*(x: SymNode): SymNode

proc `-`*(x: SymNode): SymNode = x.neg
proc `+`*(x: SymNode, y: SymNode): SymNode = x.add y
proc `*`*(x: SymNode, y: SymNode): SymNode = x.mul y
proc `-`*(x: SymNode, y: SymNode): SymNode = x.sub y
proc `/`*(x: SymNode, y: SymNode): SymNode = x.divide y

proc identForward(z: SymNode) =
  identSymNodeValue(z.value, z.inputs[0].value)

proc identAllocate(z: SymNode) =
  z.value = identAllocateSymNodeValue(z.inputs[0].value)

proc identBackward(z: SymNode, i: int, dep: SymNode): SymNode =
  if i != 0:
    raiseError("ident has 1 operand, got i = " & $i)
  let g = z.gradientDependentOrNil dep
  if g == nil:
    return newSymNodeFloat(floatValue = 1.0, name = "One[ident]")
  else:
    return g.ident

proc ident*(x:SymNode): SymNode =
  newSymNode(
    inputs = @[x],
    forward = identForward,
    allocateValue = identAllocate,
    backward = identBackward,
    name = "ident")

proc zeroAllocate(z: SymNode) =
  z.value = zeroAllocateSymNodeValue(z.inputs[0].value)

proc zeroForward(z: SymNode) =
  zeroSymNodeValue(z.value, z.inputs[0].value)

proc zeroBackward(z: SymNode, i: int, dep: SymNode): SymNode =
  if i != 0:
    raiseError("ident has 1 operand, got i = " & $i)
  return newSymNodeFloat(floatValue = 0.0, name = "Zero[zero]")

proc zero*(x:SymNode): SymNode =
  newSymNode(
    inputs = @[x],
    forward = zeroForward,
    allocateValue = zeroAllocate,
    backward = zeroBackward,
    name = "zero")

proc indexAllocate(z: SymNode) =
  z.value = indexAllocateSymNodeValue(z.inputs[0].value, z.arg)

proc indexForward(z: SymNode) =
  indexSymNodeValue(z.value, z.inputs[0].value, z.arg)

proc indexBackward(z: SymNode, i: int, dep: SymNode): SymNode =
  if i != 0:
    raiseError("index has 1 operand, got i = " & $i)
  let g = z.gradientDependentOrNil dep
  let v =
    if g == nil:
      newSymNodeFloat(floatValue = 1.0, name = "One[ident]")
    else:
      g.ident
  result = z.inputs[0].zero.indexUpdate(z.arg, v)

proc index*(x: SymNode, i: SymNodeValue): SymNode =
  newSymNode(
    inputs = @[x],
    arg = i,
    forward = indexForward,
    allocateValue = indexAllocate,
    backward = indexBackward,
    name = "index")

proc indexUpdateAllocate(z: SymNode) =
  z.value = indexUpdateAllocateSymNodeValue(z.inputs[0].value, z.arg, z.inputs[1].value)

proc indexUpdateForward(z: SymNode) =
  indexUpdateSymNodeValue(z.value, z.inputs[0].value, z.arg, z.inputs[1].value)

proc indexUpdateBackward(z: SymNode, i: int, dep: SymNode): SymNode =
  let g = z.gradientDependentOrNil dep
  if g == nil:
    raiseValueError("gradient of " & $dep & " with respect to " & $z & " does not exists")
  case i
  of 0:
    return g.indexUpdate(z.arg, newSymNodeFloat(floatValue = 0.0, name = "[" & $z.arg & "]=Zero[indexUpdate]"))
  of 1:
    return z.inputs[0].zero.indexUpdate(z.arg, g.index(z.arg))
  else:
    raiseError("indexUpdate has 2 operand, got i = " & $i)

proc indexUpdate*(x: SymNode, i: SymNodeValue, y: SymNode): SymNode =
  newSymNode(
    inputs = @[x,y],
    arg = i,
    forward = indexUpdateForward,
    allocateValue = indexUpdateAllocate,
    backward = indexUpdateBackward,
    name = "indexUpdate")

proc negAllocate(z: SymNode) =
  z.value = negAllocateSymNodeValue(z.inputs[0].value)

proc negForward(z: SymNode) =
  negSymNodeValue(z.value, z.inputs[0].value)

proc negBackward(z: SymNode, i: int, dep: SymNode): SymNode =
  if i != 0:
    raiseError("neg has 1 operand, got i = " & $i)
  let g = z.gradientDependentOrNil dep
  let v = newSymNodeFloat(floatValue = -1.0, name = "NegOne[neg]")
  if g == nil:
    return v
  else:
    return v.mul g.ident

proc neg*(x: SymNode): SymNode =
  newSymNode(
    inputs = @[x],
    forward = negForward,
    allocateValue = negAllocate,
    backward = negBackward,
    name = "neg")

proc invAllocate(z: SymNode) =
  z.value = invAllocateSymNodeValue(z.inputs[0].value)

proc invForward(z: SymNode) =
  invSymNodeValue(z.value, z.inputs[0].value)

proc invBackward(z: SymNode, i: int, dep: SymNode): SymNode =
  if i != 0:
    raiseError("inv has 1 operand, got i = " & $i)
  let g = z.gradientDependentOrNil dep
  let v = newSymNodeFloat(floatValue = -1.0, name = "NegOne[neg]")
  # TODO this is incorrect for matrices
  if g == nil:
    return v.mul(z).mul(z)
  else:
    return g.mul(z).mul(z)

proc inv*(x: SymNode): SymNode =
  newSymNode(
    inputs = @[x],
    forward = invForward,
    allocateValue = invAllocate,
    backward = invBackward,
    name = "inv")

proc addAllocate(z: SymNode) =
  z.value = addAllocateSymNodeValue(z.inputs[0].value, z.inputs[1].value)

proc addForward(z: SymNode) =
  addSymNodeValue(z.value, z.inputs[0].value, z.inputs[1].value)

proc addBackward(z: SymNode, i: int, dep: SymNode): SymNode =
  if i != 0 and i != 1:
    raiseError("add has 2 operands, got i = " & $i)
  let g = z.gradientDependentOrNil dep
  if g == nil:
    return newSymNodeFloat(floatValue = 1.0, name = "One[add]")
  else:
    return g.ident

proc add*(x: SymNode, y: SymNode): SymNode =
  newSymNode(
    inputs = @[x, y],
    forward = addForward,
    allocateValue = addAllocate,
    backward = addBackward,
    name = "add")

proc mulAllocate(z: SymNode) =
  z.value = mulAllocateSymNodeValue(z.inputs[0].value, z.inputs[1].value)

proc mulForward(z: SymNode) =
  mulSymNodeValue(z.value, z.inputs[0].value, z.inputs[1].value)

proc mulBackward(z: SymNode, i: int, dep: SymNode): SymNode =
  let g = z.gradientDependentOrNil dep
  # TODO this is incorrect for matrices
  case i
  of 0:
    if g == nil:
      return z.inputs[1]
    else:
      return g.mul z.inputs[1]
  of 1:
    if g == nil:
      return z.inputs[0]
    else:
      return g.mul z.inputs[0]
  else:
    raiseError("mul has 2 operands, got i = " & $i)

proc mul*(x: SymNode, y: SymNode): SymNode =
  newSymNode(
    inputs = @[x, y],
    forward = mulForward,
    allocateValue = mulAllocate,
    backward = mulBackward,
    name = "mul")

proc subAllocate(z: SymNode) =
  z.value = subAllocateSymNodeValue(z.inputs[0].value, z.inputs[1].value)

proc subForward(z: SymNode) =
  subSymNodeValue(z.value, z.inputs[0].value, z.inputs[1].value)

proc subBackward(z: SymNode, i: int, dep: SymNode): SymNode =
  let g = z.gradientDependentOrNil dep
  case i
  of 0:
    if g == nil:
      return newSymNodeFloat(floatValue = 1.0, name = "One[sub]")
    else:
      return g
  of 1:
    let v = newSymNodeFloat(floatValue = -1.0, name = "NegOne[sub]")
    if g == nil:
      return v
    else:
      return v.mul g
  else:
    raiseError("mul has 2 operands, got i = " & $i)

proc sub*(x: SymNode, y: SymNode): SymNode =
  newSymNode(
    inputs = @[x, y],
    forward = subForward,
    allocateValue = subAllocate,
    backward = subBackward,
    name = "sub")

proc divideAllocate(z: SymNode) =
  z.value = divideAllocateSymNodeValue(z.inputs[0].value, z.inputs[1].value)

proc divideForward(z: SymNode) =
  divideSymNodeValue(z.value, z.inputs[0].value, z.inputs[1].value)

proc divideBackward(z: SymNode, i: int, dep: SymNode): SymNode =
  let g = z.gradientDependentOrNil dep
  # TODO this is incorrect for matrices
  case i
  of 0:
    if g == nil:
      return z.inputs[1].inv
    else:
      return g.divide z.inputs[1]
  of 1:
    let v = newSymNodeFloat(floatValue = -1.0, name = "NegOne[neg]")
    if g == nil:
      return v.mul(z.divide z.inputs[1])
    else:
      return v.mul(g.mul(z.divide z.inputs[1]))
  else:
    raiseError("divide has 2 operands, got i = " & $i)

proc divide*(x: SymNode, y: SymNode): SymNode =
  newSymNode(
    inputs = @[x, y],
    forward = divideForward,
    allocateValue = divideAllocate,
    backward = divideBackward,
    name = "divide")

proc expAllocate(z: SymNode) =
  z.value = expAllocateSymNodeValue(z.inputs[0].value)

proc expForward(z: SymNode) =
  expSymNodeValue(z.value, z.inputs[0].value)

proc expBackward(z: SymNode, i: int, dep: SymNode): SymNode =
  if i != 0:
    raiseError("exp has 1 operand, got i = " & $i)
  let g = z.gradientDependentOrNil dep
  # TODO this is incorrect for matrices
  if g == nil:
    return z
  else:
    return g.mul z

proc exp*(x: SymNode): SymNode =
  newSymNode(
    inputs = @[x],
    forward = expForward,
    allocateValue = expAllocate,
    backward = expBackward,
    name = "exp")

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

proc evalRec(z: SymNode) =
  if sntVisited in z.tag:
    return
  z.tag.incl sntVisited  # if there's a cycle, it won't get evaluated
  var highEpoch = 0
  for i in z.inputs:
    i.evalRec
    if highEpoch < i.epoch:
      highEpoch = i.epoch
  if z.epoch < highEpoch:  # TODO: leaf nodes with functions not handled properly
    if z.forward != nil:
      z.forward z
      z.runCount.inc
    elif z.inputs.len > 0:
      raiseError("inputs.len: " & $z.inputs.len & ", but no forward function defined for:\n" & z.nodeRepr)
    z.epoch = highEpoch

proc eval*(z: SymNode) =
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
            # We need to combine the gradient without breaking the existing graph.
            # Because the previous built graph may have a reference of this node, `grad`,
            # our new node has to reuse `grad`. We use a copy of `grad` and assign back.
            #echo "TODO: recombine grad: when do we need to reuse the node? Reference counting?"
            # In this gradient call during the previous traversal, the grad node was used by the children of this node.
            #echo "      ",grad.nodeRepr
            grad.assignSymNode(grad.copySymNode.add childGrad)
        else:
          raiseError(z.nodeRepr & ":" & $i & ":" & input.nodeRepr & ": visited but no gradient")
      else:
        # Not visited this time.  Construct gradient if needed.
        input.tag.incl sntVisited
        if grad == nil:
          let childGrad = z.backward(z, i, dep)
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
  result = $z & $z.tag & ": " & $z.value & ", run: " & $z.runCount & ", epoch: " & $z.epoch & ", ref: " & $z.refCount
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
  var failed = 0
  proc assertRunCount(z: SymNode, c: int) =
    if z.runCount != c:
      failed.inc
      echo "Failed: expect run count: ",c," but got node: ",z.nodeRepr

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

  w.assertRunCount 1
  z.assertRunCount 2
  dwdx.assertRunCount 1

  x.assign 0.7
  y.assign 0.5
  echo "\n## x.assign 0.7 & y.assign 0.5\n"
  dzdx.eval
  echo dzdx.treeRepr
  dwdx.eval
  echo dwdx.treeRepr
  w.eval
  echo w.treeRepr

  w.assertRunCount 2
  z.assertRunCount 3
  dwdx.assertRunCount 1

  if failed > 0:
    quit 1
