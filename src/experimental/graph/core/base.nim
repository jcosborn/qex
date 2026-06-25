from strutils import toHex, strip
import std/[sets, tables]

type
  NodeId* = uint64
  NodeKey* = pointer
  NodeTable*[T] = Table[NodeKey, T]
  NodeSet* = HashSet[NodeKey]
  Gvalue* {.acyclic.} = ref object of RootObj
    ## A value participates in the graph only through traversal-visible dependencies.
    ## Construct graph nodes through `graphNode` or `newMultiOutputNode` so runtime
    ## and dependency invariants stay checked.
    inputs*: seq[Gvalue]
    gfunc*: Gfunc
    epoch*: int
    staticZeroLeaf*: bool
    stableId: NodeId
    runtime*: GraphRuntime
  GnodeVisit* = proc(n: Gvalue) {.closure.}
  InputWalkMode* = enum
    iwmEval, iwmReachable, iwmBackward
  ## Recomputes `z`'s cached value. A forward hook may copy concrete value
  ## storage into `z`, but must not change canonical topology: `z.inputs`,
  ## `z.gfunc`, or unrelated nodes.
  GforwardHook* = proc(z: Gvalue)
  ## `zb` is the upstream gradient; nil means `z` is the root of `grad(dep, x)`.
  ## The hook returns this node's contribution to backward input `i`, the actual
  ## backward dependency `input` (which need not line up with raw input order).
  GbackwardHook* = proc(
    zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue
  ## Emits this node's dependency surface for one traversal mode.
  ## A nil hook means "walk raw inputs" for every mode.
  GinputViewHook* = proc(z: Gvalue, mode: InputWalkMode, visit: GnodeVisit)
  Gfunc* {.acyclic.} = ref object of RootObj
    ## Represents a graph operation: inputs -> output.
    ## Dependency exposure contract:
    ## - plain dependencies live in `z.inputs`
    ## - non-raw eval/reachable/backward surfaces live in `inputView`
    forward*: GforwardHook
    backward*: GbackwardHook
    inputView*: GinputViewHook
    name*: string
  GradCacheEntry* = ref object
    revision*: uint64
    ## Contains only adjoints from completed gradient builds.
    adjoints*: Table[NodeId, Gvalue]
  GradCacheStats* = object
    revisionHits*: int
    revisionMisses*: int
    directHits*: int
    invalidations*: int
  ApplyCacheKey* = object
    selectedLambdaId*: NodeId
    revision*: uint64
  ApplyCacheEntry* = ref object
    key*: ApplyCacheKey
    instantiated*: Gvalue
  ApplyCacheStats* = object
    instantiationHits*: int
    instantiationMisses*: int
  FunctionalRuntimeState* = object
    applyCacheByNode*: Table[NodeId, ApplyCacheEntry]
    applyCacheStats*: ApplyCacheStats
  GraphRuntime* = ref object
    nextStableNodeId*: NodeId
    symbolicRevision*: uint64
    graphEpochCounter*: int
    graphDebug*: bool
    functional*: FunctionalRuntimeState
    gradCacheByOutput*: Table[NodeId, GradCacheEntry]
    gradCacheStats*: GradCacheStats
    runCountsByNode*: Table[NodeId, int]

type
  GraphError* = object of CatchableError
  GraphValueError* = object of GraphError

template raiseError*(msg: string) =
  raise newException(GraphError, msg)

template raiseValueError*(msg: string) =
  raise newException(GraphValueError, msg)

template raiseErrorBaseMethod*(msg: string) =
  raiseError(
    "Base method invoked: " & msg &
    "\nCheck that graph-building code keeps concrete Gvalue subtypes, and cast erased inputs back to their expected type in backward builders.")

proc initGraphRuntime*(): GraphRuntime =
  GraphRuntime(
    functional: FunctionalRuntimeState(
      applyCacheByNode: initTable[NodeId, ApplyCacheEntry]()),
    gradCacheByOutput: initTable[NodeId, GradCacheEntry](),
    runCountsByNode: initTable[NodeId, int]())

proc runCount*(x: Gvalue): int =
  if x.stableId == 0:
    return 0
  x.runtime.runCountsByNode.getOrDefault(x.stableId)

proc `$`*(x: Gfunc): string = x.name

method `$`*(x: Gvalue): string {.base.} =
  let f = x.gfunc
  result = "Gvalue(" & $x.epoch & ")"
  if f != nil:
    result &= " " & $f

proc nodeRepr*(x: Gvalue): string =
  let f = x.gfunc
  result = $x & " (" & $x.epoch & ")" & "@0X" &
    strip(toHex(cast[int](x)), trailing = false, chars = {'0'})
  if f != nil:
    result &= " " & $f & "@0X" &
      strip(toHex(cast[int](f)), trailing = false, chars = {'0'})

proc sharedGraphRuntime*(values: openArray[Gvalue],
                         label = "graph value"): GraphRuntime =
  ## The single runtime-identity check: every value must be a constructed graph
  ## value (non-nil, runtime-bearing) and all must share one runtime, which is
  ## returned. Empty input returns nil, which graph-node construction reads as
  ## "no input edges". `[a, b]`/`[a, b, c]` covers the fixed-arity callers, so
  ## there is no separate two-value variant.
  for value in values:
    if value == nil:
      raiseValueError(label & " is nil")
    if value.runtime == nil:
      raiseValueError(label & " has nil runtime")
    if value.stableId == 0:
      raiseValueError(label & " has no stable node id")
    if result == nil:
      result = value.runtime
    elif result != value.runtime:
      raiseValueError(label & " mixes multiple graph runtimes")

proc assignStableNodeId*[T: Gvalue](x: T): T {.discardable.} =
  # Callers are always a graphNode result (pre-validated) or a non-nil object
  # literal with a non-nil runtime, so a malformed value crashes here as user error.
  let grt = x.runtime
  if x.stableId == 0:
    inc grt.nextStableNodeId
    x.stableId = grt.nextStableNodeId
  x

proc stableNodeId*(x: Gvalue): NodeId =
  if x == nil:
    raiseValueError("stable node id target is nil")
  if x.stableId == 0:
    raiseValueError("graph value has no stable node id:\n" & x.nodeRepr)
  x.stableId

proc nodeKey*(x: Gvalue): NodeKey {.inline.} =
  cast[NodeKey](x)

proc markSeenNode*(seen: var NodeSet, x: Gvalue): bool {.inline.} =
  ## Returns true only for the first visit to a node key.
  if x.nodeKey in seen:
    return false
  seen.incl x.nodeKey
  true

method newOneOf*(x: Gvalue): Gvalue {.base.} =
  raiseErrorBaseMethod("newOneOf(" & $x & ")")
method valCopy*(z: Gvalue, x: Gvalue) {.base.} =
  raiseErrorBaseMethod("valCopy(" & $z & "," & $x & ")")
method copyCompatible*(prototype: Gvalue, value: Gvalue): bool {.base.} =
  false
method zeroLike*(x: Gvalue): Gvalue {.base.} = x.newOneOf
method oneLike*(x: Gvalue): Gvalue {.base.} =
  raiseErrorBaseMethod("oneLike(" & $x & ")")
method rootGradientSeed*(x: Gvalue): Gvalue {.base.} =
  x.oneLike

proc rootedUpstream*(upstream: Gvalue, node: Gvalue): Gvalue =
  ## Effective upstream adjoint inside a backward hook: the given `upstream`, or
  ## `node`'s own root seed when this node is the root of the backward build
  ## (`upstream == nil`). Structural nodes (cond, multi, apply) seed from
  ## `rootGradientSeed`; leaf ops use `scaledUpstreamOr` instead.
  if upstream == nil: node.rootGradientSeed else: upstream

method addLike*(prototype: Gvalue, x: Gvalue, y: Gvalue): Gvalue {.base.} =
  raiseErrorBaseMethod("addLike(" & $prototype & "," & $x & "," & $y & ")")

method isZero*(x: Gvalue): bool {.base.} =
  raiseErrorBaseMethod("isZero(" & $x & ")")

proc markStaticZeroLeaf*[T: Gvalue](v: T): T {.discardable.} =
  if v.gfunc != nil or v.inputs.len != 0:
    raiseValueError("static zero must be an inputless leaf")
  if not v.isZero:
    raiseValueError("static zero leaf must hold a zero value")
  v.staticZeroLeaf = true
  v

proc isStaticZeroLeaf*(v: Gvalue): bool =
  v.staticZeroLeaf and v.isZero

method scaleLike*(contribution: Gvalue, upstream: Gvalue): Gvalue {.base.} =
  raiseErrorBaseMethod("scaleLike(" & $contribution & ", " & $upstream & ")")

proc graphNode*[T: Gvalue, I: Gvalue](node: T,
                                      inputs: openArray[I],
                                      gfunc: Gfunc,
                                      label = "graph node"): T =
  var nodeInputs: seq[Gvalue] = @[]
  for input in inputs:
    nodeInputs.add input
  if nodeInputs.len > 0 and gfunc == nil:
    raiseValueError(label & " with inputs requires a graph function")
  let erasedNode: Gvalue = node
  if erasedNode == nil:
    raiseValueError(label & " result is nil")
  if erasedNode.runtime == nil:
    raiseValueError(label & " result has nil runtime")
  let inputGrt = sharedGraphRuntime(nodeInputs, label)
  # `nil` only means there were no inputs; values themselves have runtimes.
  if inputGrt != nil and erasedNode.runtime != inputGrt:
    raiseValueError(label & " mixes multiple graph runtimes")
  node.inputs = nodeInputs
  node.gfunc = gfunc
  if nodeInputs.len > 0 or gfunc != nil:
    node.staticZeroLeaf = false
  node.assignStableNodeId
  node
