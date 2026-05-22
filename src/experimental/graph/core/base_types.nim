from strutils import toHex, strip
import std/tables
import std/sets

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
    stableId: NodeId
    runtime*: GraphRuntime
  GnodeVisit* = proc(n: Gvalue) {.closure.}
  InputWalkMode* = enum
    iwmEval, iwmGradSignature, iwmDepend
  GradSigTokenKind* = enum
    gstNode, gstInput, gstCallable, gstFunc
  GradSigToken* = object
    kind*: GradSigTokenKind
    key*: uint64
  GforwardHook* = proc(z: Gvalue)
  ## `zb` is the upstream gradient; nil means `z` is the root of `grad(dep, x)`.
  ## The hook returns this node's contribution to raw input `i`.
  GbackwardHook* = proc(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue
  ## Appends operation-specific structural tokens to the gradient cache signature.
  GsignatureHook* = proc(z: Gvalue, tokens: var seq[GradSigToken])
  GinputWalkHook* = proc(z: Gvalue, visit: GnodeVisit)
  GdepWalks* = object
    ## A nil hook means "walk raw inputs" for that mode.
    eval*: GinputWalkHook
    gradSignature*: GinputWalkHook
    depend*: GinputWalkHook
  ## Overrides raw-input propagation when a node can contribute directly to
  ## `target`. `dep` is the original output/root passed to `grad(dep, target)`.
  ## Returning nil means this node has no direct contribution to `target`;
  ## unlike raw `backward`, nil is a valid zero-contribution result here.
  GtargetBackwardHook* = proc(zb: Gvalue, z: Gvalue, target: Gvalue, dep: Gvalue): Gvalue
  ## Reduces a lazy node to the concrete value used by callable resolution.
  GcallableReduceHook* = proc(z: Gvalue): Gvalue
  Gfunc* {.acyclic.} = ref object
    ## Represents a graph operation: inputs -> output.
    ## Dependency exposure contract:
    ## - plain dependencies live in `z.inputs`
    ## - hidden dependencies must be surfaced via a mode-specific dependency walk
    ##   hook or `walkHiddenDeps`
    ## - cache-signature-only structure may be surfaced via `signature`/`appendSignatureTokens`
    forward*: GforwardHook
    backward*: GbackwardHook
    signature*: GsignatureHook
    depWalks*: GdepWalks
    backwardTarget*: GtargetBackwardHook
    reduceCallable*: GcallableReduceHook
    ## Stable structural key for static operator metadata that is not in inputs.
    signatureKey*: uint64
    name*: string
  GradSignature* = seq[GradSigToken]
  GradCacheEntry* = ref object
    hasSignature*: bool
    signature*: GradSignature
    grads*: Table[NodeId, Gvalue]
  GradCacheStats* = object
    signatureHits*: int
    signatureMisses*: int
    directHits*: int
    directMisses*: int
    invalidations*: int
  ApplyCacheKey* = object
    callableKey*: NodeId
    inputKeys*: seq[NodeId]
  ApplyCacheEntry* = ref object
    key*: ApplyCacheKey
    reduced*: Gvalue
    partials*: Table[NodeId, Gvalue]
  ApplyCacheStats* = object
    reduceHits*: int
    reduceMisses*: int
    partialHits*: int
    partialMisses*: int
  GraphRuntime* = ref object
    nextStableNodeId*: NodeId
    graphEpochCounter*: int
    graphDebug*: bool
    gradCacheByOutput*: Table[NodeId, GradCacheEntry]
    gradCacheStats*: GradCacheStats
    applyCacheByNode*: Table[NodeId, ApplyCacheEntry]
    applyCacheStats*: ApplyCacheStats
    runCountsByNode*: Table[NodeId, int]
    lambdaResolveDepthLimit*: int
    applyGradPrepareDepthLimit*: int
    applyGradPrepareDepth*: int

type
  GraphError* = object of CatchableError
  GraphValueError* = object of GraphError
  GraphUnresolvedValueError* = object of GraphValueError

template raiseError*(msg: string) =
  raise newException(GraphError, msg)

template raiseValueError*(msg: string) =
  raise newException(GraphValueError, msg)

template raiseUnresolvedValueError*(msg: string) =
  raise newException(GraphUnresolvedValueError, msg)

template raiseErrorBaseMethod*(msg: string) =
  raiseError(
    "Base method invoked: " & msg &
    "\nCheck that graph-building code keeps concrete Gvalue subtypes, and cast erased inputs back to their expected type in backward builders.")

proc initGraphRuntime*(): GraphRuntime =
  GraphRuntime(
    gradCacheByOutput: initTable[NodeId, GradCacheEntry](),
    applyCacheByNode: initTable[NodeId, ApplyCacheEntry](),
    runCountsByNode: initTable[NodeId, int](),
    lambdaResolveDepthLimit: 64,
    applyGradPrepareDepthLimit: 4096)

proc clearGradCache*(grt: GraphRuntime) =
  grt.gradCacheByOutput = initTable[NodeId, GradCacheEntry]()

proc resetGradCacheStats*(grt: GraphRuntime) =
  grt.gradCacheStats = GradCacheStats()

proc resetGradCache*(grt: GraphRuntime) =
  clearGradCache(grt)
  resetGradCacheStats(grt)

proc clearApplyCache*(grt: GraphRuntime) =
  grt.applyCacheByNode = initTable[NodeId, ApplyCacheEntry]()

proc resetApplyCacheStats*(grt: GraphRuntime) =
  grt.applyCacheStats = ApplyCacheStats()

proc resetApplyCache*(grt: GraphRuntime) =
  clearApplyCache(grt)
  resetApplyCacheStats(grt)

proc attachRuntime*[T: Gvalue](x: T, grt: GraphRuntime): T {.discardable.} =
  if x == nil:
    raiseValueError("cannot attach graph runtime to nil value")
  if grt == nil:
    raiseValueError("graph value requires non-nil runtime")
  if x.runtime != nil and x.runtime != grt:
    raiseValueError("graph value runtime is already set")
  x.runtime = grt
  x

proc newGfunc*(
    forward: GforwardHook = nil,
    backward: GbackwardHook = nil,
    signature: GsignatureHook = nil,
    depWalks: GdepWalks = GdepWalks(),
    backwardTarget: GtargetBackwardHook = nil,
    reduceCallable: GcallableReduceHook = nil,
    signatureKey: uint64 = 0,
    name: string): Gfunc =
  if name.strip.len == 0:
    raiseValueError("graph function name cannot be empty")
  if backward != nil and backwardTarget != nil:
    raiseValueError(name & " cannot define both backward and backwardTarget")
  Gfunc(
    forward: forward,
    backward: backward,
    signature: signature,
    depWalks: depWalks,
    backwardTarget: backwardTarget,
    reduceCallable: reduceCallable,
    signatureKey: signatureKey,
    name: name)

proc signatureKey*(f: Gfunc): uint64 =
  if f == nil:
    return 0
  if f.signatureKey != 0:
    return f.signatureKey
  cast[uint64](f)

proc nodeKey*(x: Gvalue): NodeKey {.inline.}
proc stableNodeId*(x: Gvalue): NodeId

proc runCount*(x: Gvalue): int =
  if x.gfunc == nil or x.stableId == 0:
    return 0
  x.runtime.runCountsByNode.getOrDefault(x.stableId)

proc `$`*(x: Gfunc): string

method `$`*(x: Gvalue): string {.base.} =
  let f = x.gfunc
  result = "Gvalue(" & $x.epoch & ")"
  if f != nil:
    result &= " " & $f

proc `$`*(x: Gfunc): string = x.name

proc nodeRepr*(x: Gvalue): string =
  let f = x.gfunc
  result = $x & " (" & $x.epoch & ")" & "@0X" &
    strip(toHex(cast[int](x)), trailing = false, chars = {'0'})
  if f != nil:
    result &= " " & $f & "@0X" &
      strip(toHex(cast[int](f)), trailing = false, chars = {'0'})

proc sharedGraphRuntime*(values: openArray[Gvalue]): GraphRuntime =
  for value in values:
    let grt = value.runtime
    if result == nil:
      result = grt
    elif result != grt:
      raiseValueError("graph node mixes multiple graph runtimes")

proc stableNodeId*(x: Gvalue): NodeId =
  let grt = x.runtime
  if x.stableId == 0:
    inc grt.nextStableNodeId
    x.stableId = grt.nextStableNodeId
  x.stableId

proc depWalkForMode*(f: Gfunc,
                     mode: InputWalkMode): GinputWalkHook =
  if f == nil:
    return nil
  case mode
  of iwmEval:
    result = f.depWalks.eval
  of iwmGradSignature:
    result = f.depWalks.gradSignature
  of iwmDepend:
    result = f.depWalks.depend

proc nodeKey*(x: Gvalue): NodeKey {.inline.} =
  cast[NodeKey](x)

proc markSeenNode*(seen: var NodeSet, x: Gvalue): bool {.inline.} =
  ## Returns true only for the first visit to a node key.
  if x.nodeKey in seen:
    return false
  seen.incl x.nodeKey
  true
