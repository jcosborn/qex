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
    ## Fields other than `inputs` are ignored unless a custom walk/signature hook exposes them.
    inputs*: seq[Gvalue]
    gfunc*: Gfunc
    epoch*: int
    stableId: NodeId
    runtime*: GraphRuntime
  GnodeVisit* = proc(n: Gvalue) {.closure.}
  InputWalkMode* = enum
    iwmEval, iwmGradSignature, iwmDepend
  DepWalkKind* = enum
    dwkRawInputs, dwkNoInputs, dwkCustomInputs
  DeferredApplyKind* = enum
    dakNone, dakApply, dakApplyPartial
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
  GdepWalkSpec* = object
    kind*: DepWalkKind
    inputs*: GinputWalkHook
  GdepWalks* = object
    eval*: GdepWalkSpec
    gradSignature*: GdepWalkSpec
    depend*: GdepWalkSpec
  ## Overrides raw-input propagation when a node can contribute directly to `target`.
  ## `dep` is the original output/root passed to `grad(dep, target)`.
  GtargetBackwardHook* = proc(zb: Gvalue, z: Gvalue, target: Gvalue, dep: Gvalue): Gvalue
  ## Reduces a lazy node to the concrete value used by callable resolution.
  GreduceValueHook* = proc(z: Gvalue): Gvalue
  Gfunc* {.acyclic.} = ref object
    ## Represents a graph operation: inputs -> output.
    ## Dependency exposure contract:
    ## - plain dependencies live in `z.inputs`
    ## - hidden dependencies must be surfaced via a mode-specific dependency walk
    ##   spec or `walkHiddenDeps`
    ## - cache-signature-only structure may be surfaced via `signature`/`appendSignatureTokens`
    forward*: GforwardHook
    backward*: GbackwardHook
    signature*: GsignatureHook
    depWalks*: GdepWalks
    backwardTarget*: GtargetBackwardHook
    reduceValue*: GreduceValueHook
    deferredApplyKind*: DeferredApplyKind
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
    hasKey*: bool
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
    runCountsByFunc*: Table[uint64, int]
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
    "\nMake sure to pass `--multimethods:on` and check there is a custom method for each derived type.")

proc initGraphRuntime*(): GraphRuntime =
  GraphRuntime(
    gradCacheByOutput: initTable[NodeId, GradCacheEntry](),
    applyCacheByNode: initTable[NodeId, ApplyCacheEntry](),
    runCountsByFunc: initTable[uint64, int](),
    lambdaResolveDepthLimit: 64,
    applyGradPrepareDepthLimit: 4096)

const
  rawInputDepWalkSpec* = GdepWalkSpec(kind: dwkRawInputs)
  noInputDepWalkSpec* = GdepWalkSpec(kind: dwkNoInputs)

proc rawInputDepWalk*(): GdepWalkSpec =
  rawInputDepWalkSpec

proc noInputDepWalk*(): GdepWalkSpec =
  noInputDepWalkSpec

proc walkedInputs*(inputs: GinputWalkHook): GdepWalkSpec =
  if inputs == nil:
    raiseValueError("walkedInputs requires a callback")
  GdepWalkSpec(
    kind: dwkCustomInputs,
    inputs: inputs)

proc newDepWalks*(
    eval: GdepWalkSpec = rawInputDepWalkSpec,
    gradSignature: GdepWalkSpec = rawInputDepWalkSpec,
    depend: GdepWalkSpec = rawInputDepWalkSpec): GdepWalks =
  GdepWalks(
    eval: eval,
    gradSignature: gradSignature,
    depend: depend)

proc newGfunc*(
    forward: GforwardHook = nil,
    backward: GbackwardHook = nil,
    signature: GsignatureHook = nil,
    depWalks: GdepWalks = newDepWalks(),
    backwardTarget: GtargetBackwardHook = nil,
    reduceValue: GreduceValueHook = nil,
    deferredApplyKind: DeferredApplyKind = dakNone,
    signatureKey: uint64 = 0,
    name: string): Gfunc =
  Gfunc(
    forward: forward,
    backward: backward,
    signature: signature,
    depWalks: depWalks,
    backwardTarget: backwardTarget,
    reduceValue: reduceValue,
    deferredApplyKind: deferredApplyKind,
    signatureKey: signatureKey,
    name: name)

proc signatureKey*(f: Gfunc): uint64 =
  if f == nil:
    return 0
  if f.signatureKey != 0:
    return f.signatureKey
  cast[uint64](f)

proc gfuncKey(f: Gfunc): uint64 {.inline.} =
  cast[uint64](f)

proc recordRun*(grt: GraphRuntime, f: Gfunc) =
  ## Run counts are keyed by `Gfunc` identity, not by node instance.
  ## `grt` is non-nil by node construction; eval records only runnable funcs.
  let key = f.gfuncKey
  grt.runCountsByFunc[key] = grt.runCountsByFunc.getOrDefault(key) + 1

proc runCount*(grt: GraphRuntime, f: Gfunc): int =
  ## Returns executions recorded for this function object in the runtime.
  ## `grt` is non-nil by node construction; callers pass a checked function.
  grt.runCountsByFunc.getOrDefault(f.gfuncKey)

proc runCount*(x: Gvalue): int =
  if x == nil or x.gfunc == nil:
    return 0
  x.runtime.runCount(x.gfunc)

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

proc epochOf*(x: Gvalue): int =
  if x == nil:
    return 0
  x.epoch

proc sharedGraphRuntime*(values: openArray[Gvalue]): GraphRuntime =
  for value in values:
    if value == nil:
      continue
    let grt = value.runtime
    if result == nil:
      result = grt
    elif result != grt:
      raiseValueError("graph node mixes multiple graph runtimes")

proc stableNodeId*(x: Gvalue): NodeId =
  if x == nil:
    return 0
  let grt = x.runtime
  if x.stableId == 0:
    inc grt.nextStableNodeId
    x.stableId = grt.nextStableNodeId
  x.stableId

proc depWalkForMode*(f: Gfunc,
                     mode: InputWalkMode): GdepWalkSpec =
  if f == nil:
    return rawInputDepWalkSpec
  case mode
  of iwmEval:
    result = f.depWalks.eval
  of iwmGradSignature:
    result = f.depWalks.gradSignature
  of iwmDepend:
    result = f.depWalks.depend

proc runDepWalkInputs*(x: Gvalue,
                       spec: GdepWalkSpec,
                       visit: GnodeVisit): bool =
  if x == nil or spec.kind != dwkCustomInputs or spec.inputs == nil:
    return false
  spec.inputs(x, visit)
  true

proc appendGfuncSignature*(x: Gvalue,
                           tokens: var seq[GradSigToken]) =
  if x == nil or x.gfunc == nil or x.gfunc.signature == nil:
    return
  x.gfunc.signature(x, tokens)

proc hasTargetBackward*(x: Gvalue): bool {.inline.} =
  x != nil and x.gfunc != nil and x.gfunc.backwardTarget != nil

proc runTargetBackward*(x: Gvalue,
                        zb: Gvalue,
                        target: Gvalue,
                        dep: Gvalue): Gvalue =
  if not x.hasTargetBackward:
    return nil
  x.gfunc.backwardTarget(zb, x, target, dep)

proc hasReduceValue*(x: Gvalue): bool {.inline.} =
  x != nil and x.gfunc != nil and x.gfunc.reduceValue != nil

proc runReduceValue*(x: Gvalue): Gvalue =
  if not x.hasReduceValue:
    return nil
  x.gfunc.reduceValue(x)

proc nodeKey*(x: Gvalue): NodeKey {.inline.} =
  cast[NodeKey](x)

proc nodeFromKey*(key: NodeKey): Gvalue {.inline.} =
  cast[Gvalue](key)

proc initNodeTable*[T](): NodeTable[T] {.inline.} =
  initTable[NodeKey, T]()

proc initNodeSet*(): NodeSet {.inline.} =
  initHashSet[NodeKey]()

proc hasNode*[T](t: NodeTable[T], x: Gvalue): bool {.inline.} =
  t.hasKey(x.nodeKey)

proc getNode*[T](t: NodeTable[T], x: Gvalue): T {.inline.} =
  t[x.nodeKey]

proc putNode*[T](t: var NodeTable[T], x: Gvalue, value: T) {.inline.} =
  t[x.nodeKey] = value

proc delNode*[T](t: var NodeTable[T], x: Gvalue) {.inline.} =
  t.del(x.nodeKey)

proc nodeOrDefault*[T](t: NodeTable[T], x: Gvalue, default: T): T {.inline.} =
  t.getOrDefault(x.nodeKey, default)

proc containsNode*(s: NodeSet, x: Gvalue): bool {.inline.} =
  x.nodeKey in s

proc inclNode*(s: var NodeSet, x: Gvalue) {.inline.} =
  s.incl x.nodeKey

proc exclNode*(s: var NodeSet, x: Gvalue) {.inline.} =
  s.excl x.nodeKey

proc markSeenNode*(seen: var NodeSet, x: Gvalue): bool {.inline.} =
  ## Returns true only for the first visit to a node key.
  if seen.containsNode(x):
    return false
  seen.inclNode x
  true
