import core, scalar
import std/tables
import std/sets

# Section: Functional Node Types and Runtime State

type
  Glocal* {.final.} = ref object of Gvalue
    retProto: Gvalue
    bound: Gvalue
  Gcallable* {.final.} = ref object of Gvalue
    retProto: Gvalue
    bound: Gvalue
  Glambda* {.final.} = ref object of Gvalue
    param: Gvalue
    body: Gvalue
    envParams: seq[Gvalue]
    envValues: seq[Gvalue]
    resolved: bool
  CallableResolution = object
    value: Gvalue
    fn: Glambda
  CallableResolveMode = enum
    crmShallow, crmReduced
  ApplySnapshot = object
    fun: Gvalue
    arg: Gvalue
    deps: seq[Gvalue]
    fn: Glambda
    callableKey: NodeKey
    freshEpoch: int
    hasReductionState: bool
  ApplyPartialView = object
    base: Gvalue
    targets: seq[Gvalue]
  ApplySignature = object
    callableKey: NodeKey
    inputKeys: seq[NodeKey]
  ApplyCacheEntry = ref object
    hasSignature: bool
    signature: ApplySignature
    freshEpoch: int
    reduced: Gvalue
    partials: NodeTable[Gvalue]
  DependencyEpochContext = object
    seenCallable: NodeSet
    seenValues: NodeSet
  CallableDepWalkMode = enum
    cdwmEpoch, cdwmValueDeps
  ApplyCacheStats* = object
    reduceHits*: int
    reduceMisses*: int
    partialHits*: int
    partialMisses*: int
  NodeBindings = NodeTable[Gvalue]

var applyCache = initNodeTable[ApplyCacheEntry]()
var applyCacheStats*: ApplyCacheStats
var gapplyDeferred: Gfunc
var gapplyPartialDeferred: Gfunc

# Section: Private Forward Declarations

# Lambda cloning and closure conversion.
proc cloneWithSubst(v: Gvalue, subst: NodeBindings,
                    memo: var NodeBindings): Gvalue
proc cloneResolvedLambdaWithSubst(fn: Glambda,
                                  subst: NodeBindings,
                                  memo: var NodeBindings): Gvalue
proc recloseLambda(fn: Glambda)
proc instantiateLambdaBody(fn: Glambda, x: Gvalue): Gvalue

# Callable resolution and apply normalization.
proc resolveCallable(fun: Gvalue, mode: CallableResolveMode): CallableResolution
proc resolveDirectLambda(fun: Gvalue): Glambda
proc resolveLambda(fun: Gvalue): Glambda
proc buildApplySnapshot(v: Gvalue): ApplySnapshot
proc parseApplyPartial(v: Gvalue): ApplyPartialView
proc prepareApplyReduction(snapshot: var ApplySnapshot)
proc appendApplySignature(v: Gvalue, tokens: var seq[GradSigToken])
proc collectCallableValueDeps(roots: openArray[Gvalue],
                              seeded: openArray[Gvalue]): seq[Gvalue]
proc inputDependencyEpoch(snapshot: ApplySnapshot): int
proc callableShellRepr(label: string, bound: Gvalue): string
proc appendCallableToken(v: Gvalue, tokens: var seq[GradSigToken])
proc isCallableBoundary(v: Gvalue, mode: CallableDepWalkMode): bool
proc callableBoundaryDeps(v: Gvalue,
                          mode: CallableDepWalkMode): seq[Gvalue]

# Deferred apply-partial node construction.
proc deferredTargetNode(baseInputs: seq[Gvalue],
                        targets: seq[Gvalue],
                        gfunc: Gfunc,
                        label: string): Gvalue
proc applyPartialDeferredNode(baseInputs: seq[Gvalue], targets: seq[Gvalue]): Gvalue
proc applyPartialDeferred(dep: Gvalue, target: Gvalue): Gvalue
proc applyPartialDeferredAppendTarget(z: Gvalue, target: Gvalue): Gvalue
proc applyPartialDeferredContrib(zb: Gvalue, z: Gvalue, target: Gvalue): Gvalue
proc dependencyEpoch(ctx: var DependencyEpochContext, v: Gvalue): int
proc dependencyEpoch(v: Gvalue): int
proc inputDependencyEpoch(ctx: var DependencyEpochContext,
                          inputs: openArray[Gvalue]): int
proc inputDependencyEpoch(inputs: openArray[Gvalue]): int
proc ensureApplyReduction(v: Gvalue): ApplyCacheEntry
proc requireApplyReduction(v: Gvalue): ApplyCacheEntry
proc setApplyReduction(entry: ApplyCacheEntry,
                       snapshot: ApplySnapshot,
                       reduced: Gvalue)
proc freshCallableBound(v: Gvalue): Gvalue

# Section: Runtime Limits and Stats Reset

var lambdaResolveDepthLimit* = 64
var applyGradPrepareDepthLimit* = 4096
var applyGradPrepareDepth = 0

proc resetApplyCacheStats*() =
  applyCache = initNodeTable[ApplyCacheEntry]()
  applyCacheStats = ApplyCacheStats()

# Section: Small Shared Helpers

proc copySeq[T](src: seq[T]): seq[T] =
  result = newseq[T](src.len)
  for i in 0..<src.len:
    result[i] = src[i]

proc initNodeBindings(): NodeBindings =
  initNodeTable[Gvalue]()

proc updateMax(value: var int, candidate: int) =
  if value < candidate:
    value = candidate

proc initDependencyEpochContext(): DependencyEpochContext =
  result.seenCallable = initNodeSet()
  result.seenValues = initNodeSet()

proc markSeenNode(seen: var NodeSet, v: Gvalue): bool =
  if seen.containsNode(v):
    return false
  seen.inclNode v
  true

proc bindSubst(subst: var NodeBindings, key, value: Gvalue) =
  subst.putNode(key, value)

proc dropSubst(subst: var NodeBindings, key: Gvalue) =
  if key != nil and subst.hasNode(key):
    subst.delNode(key)

proc scanDependencyInputs(ctx: var DependencyEpochContext,
                          value: var int,
                          inputs: openArray[Gvalue]) =
  for input in inputs:
    value.updateMax ctx.dependencyEpoch(input)

# Section: Callable Resolution and Apply Signatures

proc callableKey(fn: Glambda): NodeKey =
  if fn == nil:
    return nil
  fn.nodeKey

proc applySnapshotInputCount(snapshot: ApplySnapshot): int =
  2 + snapshot.deps.len

proc visitApplySnapshotInputs(snapshot: ApplySnapshot,
                              visit: proc(n: Gvalue) {.closure.}) =
  visit snapshot.fun
  visit snapshot.arg
  for dep in snapshot.deps:
    visit dep

proc applySignatureMatches(signature: ApplySignature,
                           snapshot: ApplySnapshot): bool =
  if signature.callableKey != snapshot.callableKey:
    return false
  if signature.inputKeys.len != snapshot.applySnapshotInputCount:
    return false
  if signature.inputKeys[0] != snapshot.fun.nodeKey:
    return false
  if signature.inputKeys[1] != snapshot.arg.nodeKey:
    return false
  for i in 0..<snapshot.deps.len:
    if signature.inputKeys[i + 2] != snapshot.deps[i].nodeKey:
      return false
  true

proc setApplySignature(entry: ApplyCacheEntry, snapshot: ApplySnapshot) =
  entry.signature.callableKey = snapshot.callableKey
  entry.signature.inputKeys = newseq[NodeKey](snapshot.applySnapshotInputCount)
  entry.signature.inputKeys[0] = snapshot.fun.nodeKey
  entry.signature.inputKeys[1] = snapshot.arg.nodeKey
  for i in 0..<snapshot.deps.len:
    entry.signature.inputKeys[i + 2] = snapshot.deps[i].nodeKey

proc applyCacheEntry(v: Gvalue): ApplyCacheEntry =
  if not applyCache.hasNode(v):
    applyCache.putNode(v, ApplyCacheEntry())
  applyCache.getNode(v)

proc applyEntryMatches(entry: ApplyCacheEntry,
                       snapshot: ApplySnapshot): bool =
  if not entry.hasSignature or entry.reduced == nil:
    return false
  if not entry.signature.applySignatureMatches(snapshot):
    return false
  entry.freshEpoch >= snapshot.freshEpoch

proc maxInputEpoch(v: Gvalue): int =
  if v == nil:
    return 0
  for input in v.inputs:
    if input != nil and result < input.epochOf:
      result = input.epochOf

proc freshCallableBound(v: Gvalue): Gvalue =
  ## Reuse an evaluated callable wrapper only while none of its inputs are newer.
  if v == nil or not (v of Gcallable):
    return nil
  let cv = Gcallable(v)
  if cv.bound == nil:
    return nil
  var maxep = v.maxInputEpoch
  var ctx = initDependencyEpochContext()
  let boundEpoch = ctx.dependencyEpoch(cv.bound)
  if maxep < boundEpoch:
    maxep = boundEpoch
  if v.epochOf < maxep:
    return nil
  cv.bound

proc callableResultProto(v: Gvalue): Gvalue =
  if v == nil:
    return nil
  let fn = resolveDirectLambda(v)
  if fn != nil:
    return fn.body
  if v of Gcallable:
    return Gcallable(v).retProto
  if v of Glocal:
    return Glocal(v).retProto
  if v of Glambda:
    return Glambda(v).body
  nil

proc nextCallableBinding(v: Gvalue): Gvalue =
  if v == nil:
    return nil
  if v of Glocal:
    let lf = Glocal(v)
    if lf.bound != nil:
      return lf.bound
  let boundCallable = v.freshCallableBound
  if boundCallable != nil:
    return boundCallable
  if v of Gcallable:
    let cv = Gcallable(v)
    if cv.bound != nil and v.inputs.len == 0:
      return cv.bound
  nil

proc resolveCallableValue(v: Gvalue, mode: CallableResolveMode): Gvalue =
  var current = v
  var depth = 0
  while depth < lambdaResolveDepthLimit:
    if current == nil:
      return nil
    let next = current.nextCallableBinding
    if next != nil:
      current = next
      inc depth
      continue
    if mode == crmReduced and current.gfunc == gapplyDeferred and current.inputs.len >= 2:
      try:
        let entry = current.ensureApplyReduction
        current = entry.reduced
      except GraphUnresolvedValueError:
        return nil
      inc depth
      continue
    return current
  raiseValueError(
    "callable resolution exceeded depth limit " & $lambdaResolveDepthLimit &
    ":\n" & v.nodeRepr)

proc resolveCallable(fun: Gvalue, mode: CallableResolveMode): CallableResolution =
  result.value = fun.resolveCallableValue(mode)
  if result.value == nil or not (result.value of Glambda):
    return
  let fn = Glambda(result.value)
  if fn.resolved:
    result.fn = fn

proc resolveDirectLambda(fun: Gvalue): Glambda =
  fun.resolveCallable(crmShallow).fn

proc symbolicCallableToken(v: Gvalue): pointer =
  let resolved = v.resolveCallable(crmShallow)
  if resolved.value == nil:
    return nil
  if resolved.fn != nil:
    return callableKey(resolved.fn)
  if resolved.value != v:
    return resolved.value.nodeKey
  nil

proc appendApplySignature(v: Gvalue, tokens: var seq[GradSigToken]) =
  let snapshot = v.buildApplySnapshot
  for input in [snapshot.fun, snapshot.arg]:
    let token = symbolicCallableToken(input)
    if token != nil:
      tokens.add GradSigToken(kind: gstCallable, nodePtr: token)
  for dep in snapshot.deps:
    tokens.add GradSigToken(kind: gstInput, nodePtr: dep.nodeKey)

proc isCallableLike(v: Gvalue): bool =
  v.callableResultProto != nil

proc isCallableBoundary(v: Gvalue, mode: CallableDepWalkMode): bool =
  if v == nil:
    return false
  case mode
  of cdwmEpoch:
    (v of Gcallable) or
      (v of Glocal and Glocal(v).retProto != nil) or
      (v of Glambda and Glambda(v).resolved)
  of cdwmValueDeps:
    (v of Gcallable) or
      (v of Glocal and Glocal(v).retProto != nil) or
      (v of Glambda)

proc callableBoundaryDeps(v: Gvalue,
                          mode: CallableDepWalkMode): seq[Gvalue] =
  var deps: seq[Gvalue] = @[]

  proc addDep(dep: Gvalue) =
    if dep != nil:
      deps.add dep

  case mode
  of cdwmEpoch:
    if v of Glambda:
      for ev in Glambda(v).envValues:
        addDep(ev)
      return deps
    if v of Gcallable:
      let cv = Gcallable(v)
      addDep(cv.bound)
      for input in v.inputs:
        addDep(input)
      return deps
    if v of Glocal:
      v.walkSymbolicDeps(proc(n: Gvalue) = addDep(n))
      for input in v.inputs:
        addDep(input)
      return deps
  of cdwmValueDeps:
    if v of Glambda:
      for ev in Glambda(v).envValues:
        addDep(ev)
      return deps
    if v.gfunc != nil:
      v.prepareNode
    if v of Gcallable:
      let boundCallable = v.freshCallableBound
      if boundCallable != nil:
        addDep(boundCallable)
        return deps
    v.walkSymbolicDeps(proc(n: Gvalue) = addDep(n))
    for input in v.inputs:
      addDep(input)
    return deps

# Section: Public Local and Prototype Constructors

proc local*(): Gvalue =
  result = Glocal()
  result.updated
proc local*(retPrototype: Gvalue): Gvalue =
  result = Glocal(retProto: retPrototype)
  result.updated
proc localValue*(prototype: Gvalue): Gvalue = prototype.newOneOf
proc localScalar*(): Gscalar =
  result = Gscalar()
  result.updated
proc localInt*(): Gint =
  result = Gint()
  result.updated

# Section: Functional Node Methods

method newOneOf*(x: Glocal): Gvalue = Glocal(retProto: x.retProto)
method valCopy*(z: Glocal, x: Glocal) =
  z.retProto = x.retProto
  z.bound = x.bound
  z.updated
method valCopy*(z: Glocal, x: Gvalue) =
  z.bound = x
  z.updated
method `$`*(x: Glocal): string = callableShellRepr("local", x.bound)

method newOneOf*(x: Gcallable): Gvalue = Gcallable(retProto: x.retProto)
method valCopy*(z: Gcallable, x: Gcallable) =
  z.retProto = x.retProto
  z.bound = x.bound
method valCopy*(z: Gcallable, x: Gvalue) = z.bound = x
method `$`*(x: Gcallable): string = callableShellRepr("callable", x.bound)

method newOneOf*(x: Glambda): Gvalue =
  Glambda(
    param: x.param,
    body: x.body,
    envParams: x.envParams,
    envValues: x.envValues,
    resolved: false)
method `$`*(x: Glambda): string =
  result = "lambda(" & $x.param & " -> " & $x.body & ")"
  if x.envValues.len > 0:
    result &= "[env:" & $x.envValues.len & "]"

method walkSymbolicDeps*(v: Glocal, visit: proc(n: Gvalue) {.closure.}) =
  if v.bound != nil:
    visit v.bound

method appendSignatureTokens*(v: Glocal, tokens: var seq[GradSigToken]) =
  appendCallableToken(v, tokens)

method walkSymbolicDeps*(v: Gcallable, visit: proc(n: Gvalue) {.closure.}) =
  let boundCallable = freshCallableBound(v)
  if boundCallable != nil:
    visit boundCallable
  elif v.inputs.len == 0 and v.bound != nil:
    visit v.bound

method appendSignatureTokens*(v: Gcallable, tokens: var seq[GradSigToken]) =
  appendCallableToken(v, tokens)

method walkSymbolicDeps*(v: Glambda, visit: proc(n: Gvalue) {.closure.}) =
  for ev in v.envValues:
    if ev != nil:
      visit ev

# Section: Lambda Closure Conversion and Cloning

proc collectCaptureValues(v: Gvalue, bound: var NodeSet,
                          seenCaps: var NodeSet,
                          seenNodes: var NodeSet,
                          caps: var seq[Gvalue]) =
  if v == nil:
    return
  if seenNodes.containsNode(v):
    return
  seenNodes.inclNode v
  if bound.containsNode(v):
    return

  if v of Glambda and Glambda(v).resolved:
    let fn = Glambda(v)
    for ev in fn.envValues:
      collectCaptureValues(ev, bound, seenCaps, seenNodes, caps)
    return

  let boundCallable = v.freshCallableBound
  if boundCallable != nil:
    collectCaptureValues(boundCallable, bound, seenCaps, seenNodes, caps)
    return

  if v.gfunc == nil and v.inputs.len == 0:
    if not seenCaps.containsNode(v) and not (v of Glambda):
      seenCaps.inclNode v
      caps.add v
    return

  for i in v.inputs:
    collectCaptureValues(i, bound, seenCaps, seenNodes, caps)

proc trimSubstForLambda(subst: NodeBindings, fn: Glambda): NodeBindings =
  result = subst
  result.dropSubst(fn.param)
  for ep in fn.envParams:
    result.dropSubst(ep)

proc cloneResolvedLambdaWithSubst(fn: Glambda,
                                  subst: NodeBindings,
                                  memo: var NodeBindings): Gvalue =
  result = Glambda(param: fn.param, resolved: fn.resolved)
  memo.putNode(fn, result)

  let r = Glambda(result)
  if fn.envParams.len > 0:
    r.envParams = fn.envParams.copySeq
  if fn.envValues.len > 0:
    r.envValues = newseq[Gvalue](fn.envValues.len)
    for i in 0..<fn.envValues.len:
      r.envValues[i] = cloneWithSubst(fn.envValues[i], subst, memo)
  if fn.body != nil:
    let innerSubst = subst.trimSubstForLambda(fn)
    var innerMemo = initNodeBindings()
    r.body = cloneWithSubst(fn.body, innerSubst, innerMemo)

proc cloneWithFreshMemo(v: Gvalue, subst: NodeBindings): Gvalue =
  var memo = initNodeBindings()
  cloneWithSubst(v, subst, memo)

proc bindLambdaEnvValues(subst: var NodeBindings, fn: Glambda) =
  if fn.envParams.len != fn.envValues.len:
    raiseValueError("lambda env arity mismatch during instantiation")
  for i in 0..<fn.envParams.len:
    subst.bindSubst(fn.envParams[i], fn.envValues[i])

proc cloneWithSubst(v: Gvalue, subst: NodeBindings,
                    memo: var NodeBindings): Gvalue =
  if v == nil:
    return nil

  if subst.hasNode(v):
    return subst.getNode(v)
  if memo.hasNode(v):
    return memo.getNode(v)

  if v of Glambda and Glambda(v).resolved:
    return Glambda(v).cloneResolvedLambdaWithSubst(subst, memo)

  let boundCallable = v.freshCallableBound
  if boundCallable != nil:
    result = cloneWithSubst(boundCallable, subst, memo)
    memo.putNode(v, result)
    return result

  if v.gfunc == nil and v.inputs.len == 0:
    return v

  result = v.newOneOf
  memo.putNode(v, result)
  if v.inputs.len > 0:
    result.inputs = newseq[Gvalue](v.inputs.len)
    for i in 0..<v.inputs.len:
      result.inputs[i] = cloneWithSubst(v.inputs[i], subst, memo)
  result.gfunc = v.gfunc
  if v.locals.len > 0:
    result.locals = newseq[Gvalue](v.locals.len)

proc openLambdaBody(fn: Glambda): Gvalue =
  if fn.envParams.len == 0:
    return fn.body
  var subst = initNodeBindings()
  subst.bindLambdaEnvValues(fn)
  fn.body.cloneWithFreshMemo(subst)

proc collectLambdaCaptures(fn: Glambda): seq[Gvalue] =
  var bound = initNodeSet()
  bound.inclNode fn.param
  var seenCaps = initNodeSet()
  var seenNodes = initNodeSet()
  collectCaptureValues(fn.body, bound, seenCaps, seenNodes, result)

proc rewriteLambdaCaptures(fn: Glambda, captures: seq[Gvalue]) =
  fn.envValues = captures
  fn.envParams.setLen(0)
  if captures.len == 0:
    return

  fn.envParams = newseq[Gvalue](captures.len)
  var subst = initNodeBindings()
  for i in 0..<captures.len:
    let ep = localValue(captures[i])
    fn.envParams[i] = ep
    subst.bindSubst(captures[i], ep)
  fn.body = fn.body.cloneWithFreshMemo(subst)

proc recloseLambda(fn: Glambda) =
  if fn == nil or not fn.resolved:
    return
  if fn.body == nil:
    raiseValueError("resolved lambda body cannot be nil")

  fn.body = openLambdaBody(fn)
  let captures = fn.collectLambdaCaptures
  fn.rewriteLambdaCaptures(captures)

proc instantiateLambdaBody(fn: Glambda, x: Gvalue): Gvalue =
  var subst = initNodeBindings()
  subst.bindSubst(fn.param, x)
  subst.bindLambdaEnvValues(fn)
  fn.body.cloneWithFreshMemo(subst)

proc resolveLambda(fun: Gvalue): Glambda =
  fun.resolveCallable(crmReduced).fn

# Section: Dependency Epoch Tracking

proc dependencyEpoch(ctx: var DependencyEpochContext, v: Gvalue): int =
  if v == nil:
    return 0

  if v.isCallableBoundary(cdwmEpoch):
    if not ctx.seenCallable.markSeenNode(v):
      return 0
    var maxEpoch = 0
    if v of Glocal:
      maxEpoch.updateMax v.epochOf
    for dep in v.callableBoundaryDeps(cdwmEpoch):
      maxEpoch.updateMax ctx.dependencyEpoch(dep)
    return maxEpoch

  if not ctx.seenValues.markSeenNode(v):
    return 0
  result = v.epochOf
  ctx.scanDependencyInputs(result, v.inputs)

proc dependencyEpoch(v: Gvalue): int =
  if v == nil:
    return 0
  var ctx = initDependencyEpochContext()
  result = ctx.dependencyEpoch(v)

proc inputDependencyEpoch(ctx: var DependencyEpochContext,
                          inputs: openArray[Gvalue]): int =
  for input in inputs:
    if input != nil:
      result.updateMax ctx.dependencyEpoch(input)

proc inputDependencyEpoch(inputs: openArray[Gvalue]): int =
  var ctx = initDependencyEpochContext()
  result = ctx.inputDependencyEpoch(inputs)

# Section: Public Lambda Construction and Apply Result Prototypes

proc lambda*(param: Gvalue, body: Gvalue): Gvalue =
  if param == nil:
    raiseValueError("lambda parameter cannot be nil")
  if body == nil:
    raiseValueError("lambda body cannot be nil")
  let fn = Glambda(param: param, body: body, resolved: true)
  recloseLambda(fn)
  result = fn
  result.updated

proc materializeApplyResultProto(proto: Gvalue): Gvalue =
  if proto == nil:
    return nil
  let nestedProto = proto.callableResultProto
  if nestedProto != nil:
    return Gcallable(retProto: materializeApplyResultProto(nestedProto))
  proto.newOneOf

proc applyResultProto(fun: Gvalue): Gvalue =
  let retProto = fun.callableResultProto
  if retProto == nil:
    return nil
  materializeApplyResultProto(retProto)

# Section: Deferred Apply-Partial Node Construction

proc deferredTargetNode(baseInputs: seq[Gvalue],
                        targets: seq[Gvalue],
                        gfunc: Gfunc,
                        label: string): Gvalue =
  if baseInputs.len == 0:
    raiseValueError(label & " requires at least one base input")
  if targets.len == 0:
    raiseValueError(label & " requires at least one target")
  for j in 0..<baseInputs.len:
    if baseInputs[j] == nil:
      raiseError(label & " base input cannot be nil")
  for target in targets:
    if target == nil:
      raiseError(label & " target cannot be nil")
  result = targets[^1].newOneOf
  result.inputs = newseq[Gvalue](baseInputs.len + targets.len)
  for j in 0..<baseInputs.len:
    result.inputs[j] = baseInputs[j]
  for j in 0..<targets.len:
    result.inputs[baseInputs.len + j] = targets[j]
  result.gfunc = gfunc

proc applyPartialDeferredNode(baseInputs: seq[Gvalue], targets: seq[Gvalue]): Gvalue =
  if baseInputs.len != 1:
    raiseValueError("applyPartialDeferred base input length mismatch: expected 1, got: " &
      $baseInputs.len)
  deferredTargetNode(baseInputs, targets, gapplyPartialDeferred, "applyPartialDeferred")

proc applyPartialDeferred(dep: Gvalue, target: Gvalue): Gvalue =
  if dep == nil or target == nil:
    raiseError("applyPartialDeferred has nil input")
  applyPartialDeferredNode(@[dep], @[target])

proc parseApplyPartial(v: Gvalue): ApplyPartialView =
  if v.inputs.len < 2:
    raiseValueError("applyPartialDeferred node requires a base apply and at least one target")
  result.base = v.inputs[0]
  if result.base == nil:
    raiseError("applyPartialDeferred node has nil base input:\n" & v.nodeRepr)
  result.targets = newseq[Gvalue](v.inputs.len - 1)
  for j in 1..<v.inputs.len:
    let target = v.inputs[j]
    if target == nil:
      raiseError("applyPartialDeferred target cannot be nil:\n" & v.nodeRepr)
    result.targets[j - 1] = target

proc applyPartialDeferredAppendTarget(z: Gvalue, target: Gvalue): Gvalue =
  if target == nil:
    raiseError("applyPartialDeferred append target cannot be nil")
  if z.gfunc != gapplyPartialDeferred:
    raiseValueError("applyPartialDeferred append target expects applyPartialDeferred node")
  let view = z.parseApplyPartial
  var targets = view.targets
  targets.add target
  applyPartialDeferredNode(@[view.base], targets)

proc applyPartialDeferredContrib(zb: Gvalue, z: Gvalue, target: Gvalue): Gvalue =
  result = applyPartialDeferredAppendTarget(z, target)
  if zb != nil:
    result = zb * result

# Section: Callable Dependency Collection and Apply Views

proc collectCallableValueDeps(roots: openArray[Gvalue],
                              seeded: openArray[Gvalue]): seq[Gvalue] =
  var seenCallable = initNodeSet()
  var seenValues = initNodeSet()
  var deps: seq[Gvalue] = @[]

  for dep in seeded:
    if dep != nil:
      seenValues.inclNode dep
  for root in roots:
    if root != nil:
      seenValues.inclNode root

  proc collect(v: Gvalue) =
    if v == nil:
      return
    if v.isCallableBoundary(cdwmValueDeps):
      if seenCallable.containsNode(v):
        return
      seenCallable.inclNode v
      for dep in v.callableBoundaryDeps(cdwmValueDeps):
        collect(dep)
      return
    if seenValues.containsNode(v):
      return
    seenValues.inclNode v
    deps.add v

  for root in roots:
    collect(root)
  result = deps

proc buildApplySnapshot(v: Gvalue): ApplySnapshot =
  if v.inputs.len < 2:
    raiseValueError("apply node requires at least two inputs")
  result.fun = v.inputs[0]
  result.arg = v.inputs[1]
  if result.fun == nil:
    raiseError("apply node has nil function input:\n" & v.nodeRepr)
  if result.arg == nil:
    raiseError("apply node has nil argument input:\n" & v.nodeRepr)

  let fn = resolveDirectLambda(result.fun)
  if fn != nil:
    result.deps = fn.envValues.copySeq
  for dep in collectCallableValueDeps([result.fun, result.arg], result.deps):
    result.deps.add dep

proc prepareApplyReduction(snapshot: var ApplySnapshot) =
  if snapshot.hasReductionState:
    return
  snapshot.fn = resolveLambda(snapshot.fun)
  if snapshot.fn != nil:
    snapshot.callableKey = callableKey(snapshot.fn)
  snapshot.freshEpoch = snapshot.inputDependencyEpoch
  snapshot.hasReductionState = true

proc inputDependencyEpoch(snapshot: ApplySnapshot): int =
  var ctx = initDependencyEpochContext()
  var maxEpoch = 0
  snapshot.visitApplySnapshotInputs(proc(n: Gvalue) =
    if n != nil:
      maxEpoch.updateMax ctx.dependencyEpoch(n))
  result = maxEpoch

proc visitApplySnapshotInputs(v: Gvalue, visit: proc(n: Gvalue) {.closure.}) =
  v.buildApplySnapshot.visitApplySnapshotInputs(visit)

proc emptyLike(x: Gvalue): Gvalue =
  ## Type-preserving neutral constructor; callable-like nodes may not support update(0).
  x.newOneOf

proc gradOrEmpty(expr: Gvalue, target: Gvalue): Gvalue =
  result = expr.gradIsolated(target)
  if result == nil:
    result = emptyLike(target)

proc resolvedCallableValue(v: Gvalue): Gvalue =
  ## Collapse freshly-bound callable wrappers to their callable descriptor after eval.
  let direct = v.resolveCallable(crmShallow)
  if direct.value == nil:
    return nil
  let reduced = direct.value.resolveCallable(crmReduced)
  if reduced.value != nil:
    return reduced.value
  direct.value

proc callableShellRepr(label: string, bound: Gvalue): string =
  if bound == nil:
    return label
  label & "(" & $bound & ")"

proc appendCallableToken(v: Gvalue, tokens: var seq[GradSigToken]) =
  let token = symbolicCallableToken(v)
  if token != nil:
    tokens.add GradSigToken(kind: gstCallable, nodePtr: token)

# Section: Apply Cache and Lazy Reduction

proc setApplyReduction(entry: ApplyCacheEntry,
                       snapshot: ApplySnapshot,
                       reduced: Gvalue) =
  entry.hasSignature = true
  entry.setApplySignature(snapshot)
  entry.freshEpoch = snapshot.freshEpoch
  entry.reduced = reduced
  entry.partials = initNodeTable[Gvalue]()

proc ensureApplyReduction(v: Gvalue): ApplyCacheEntry =
  var snapshot = v.buildApplySnapshot
  snapshot.prepareApplyReduction
  if snapshot.fn == nil:
    raiseUnresolvedValueError("deferred apply unresolved at eval: " & snapshot.fun.nodeRepr)

  result = v.applyCacheEntry
  if result.applyEntryMatches(snapshot):
    applyCacheStats.reduceHits.inc
    return
  applyCacheStats.reduceMisses.inc

  let reduced = instantiateLambdaBody(snapshot.fn, snapshot.arg)
  if reduced of Glambda and Glambda(reduced).resolved:
    recloseLambda(Glambda(reduced))

  result.setApplyReduction(snapshot, reduced)

proc requireApplyReduction(v: Gvalue): ApplyCacheEntry =
  result = v.ensureApplyReduction
  if result.reduced == nil:
    raiseValueError("deferred apply preparation did not produce reduced body")

# Section: Lazy Apply-Partials and Runtime Hooks

proc ensureApplyPartial(z: Gvalue, target: Gvalue): Gvalue =
  if target == nil:
    raiseValueError("apply partial target cannot be nil")
  let entry = z.requireApplyReduction
  if entry.partials.hasNode(target):
    applyCacheStats.partialHits.inc
    return entry.partials.getNode(target)
  applyCacheStats.partialMisses.inc
  if isCallableLike(target):
    result = emptyLike(target)
    entry.partials.putNode(target, result)
    return
  inc applyGradPrepareDepth
  defer:
    dec applyGradPrepareDepth
  if applyGradPrepareDepth > applyGradPrepareDepthLimit:
    raiseValueError(
      "apply partial materialization exceeded depth limit " &
      $applyGradPrepareDepthLimit &
      "\napply: " & z.nodeRepr &
      "\ntarget: " & target.nodeRepr)
  result = entry.reduced.gradOrEmpty(target)
  entry.partials.putNode(target, result)

proc applyDeferredf(v: Gvalue) =
  let entry = v.requireApplyReduction
  let reduced = entry.reduced
  let callable = resolveLambda(reduced)
  if callable != nil:
    v.valCopy callable
    return
  discard reduced.eval
  v.valCopy reduced.resolvedCallableValue

proc applyDeferredContribution(zb: Gvalue, z: Gvalue, target: Gvalue): Gvalue =
  if target == nil or isCallableLike(target):
    return nil
  let partial = applyPartialDeferred(z, target)
  if zb != nil and zb.isCallableLike:
    return nil
  if zb == nil:
    return partial
  zb * partial

proc applyDeferredb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard dep
  if i < 0 or i >= z.inputs.len:
    raiseValueError("apply backward input index out of range: " & $i)
  if i == 0:
    return nil
  applyDeferredContribution(zb, z, z.inputs[i])

proc applyDeferredBackwardTarget(zb: Gvalue,
                                 z: Gvalue,
                                 target: Gvalue,
                                 dep: Gvalue): Gvalue =
  discard dep
  applyDeferredContribution(zb, z, target)

proc applyDeferredWalkInputs(v: Gvalue,
                             mode: InputWalkMode,
                             visit: proc(n: Gvalue) {.closure.},
                             onUnknown: proc(tbranch, fbranch: Gvalue) {.closure.}) =
  discard mode
  let _ = onUnknown
  v.visitApplySnapshotInputs(visit)

proc applyPartialDeferredf(v: Gvalue) =
  let view = v.parseApplyPartial
  var expr = view.base.ensureApplyPartial(view.targets[0])
  for j in 1..<view.targets.len:
    expr = expr.gradOrEmpty(view.targets[j])
  discard expr.eval
  v.valCopy expr

proc applyPartialDeferredb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard zb
  discard dep
  let view = z.parseApplyPartial
  if i < 0 or i >= z.inputs.len:
    raiseValueError("applyPartialDeferred backward input index out of range: " & $i)
  if i == 0:
    return emptyLike(view.base)
  emptyLike(view.targets[i - 1])

proc applyPartialDeferredWalkInputs(v: Gvalue,
                                    mode: InputWalkMode,
                                    visit: proc(n: Gvalue) {.closure.},
                                    onUnknown: proc(tbranch, fbranch: Gvalue) {.closure.}) =
  let view = v.parseApplyPartial
  let _ = onUnknown
  case mode
  of iwmEval:
    var seen = initNodeSet()
    proc walkEvalDeps(n: Gvalue) =
      if n == nil:
        raiseError("applyPartialDeferred eval dependency walk encountered nil node")
      if seen.containsNode(n):
        return
      seen.inclNode n
      n.prepareNode
      if n.gfunc == gapplyDeferred or n.gfunc == gapplyPartialDeferred:
        n.walkPreparedDependInputs(walkEvalDeps)
        return
      n.walkPreparedEvalInputs(walkEvalDeps)
      n.walkSymbolicDeps(walkEvalDeps)
      visit n
    walkEvalDeps(view.base)
  of iwmGradSignature, iwmDepend:
    visit(view.base)

proc applyPartialDeferredBackwardTarget(zb: Gvalue,
                                        z: Gvalue,
                                        target: Gvalue,
                                        dep: Gvalue): Gvalue =
  discard dep
  applyPartialDeferredContrib(zb, z, target)

gapplyPartialDeferred = newGfunc(
  forward = applyPartialDeferredf,
  backward = applyPartialDeferredb,
  walkInputs = applyPartialDeferredWalkInputs,
  backwardTarget = applyPartialDeferredBackwardTarget,
  name = "applyPartialDeferred")
gapplyDeferred = newGfunc(
  forward = applyDeferredf,
  backward = applyDeferredb,
  walkInputs = applyDeferredWalkInputs,
  backwardTarget = applyDeferredBackwardTarget,
  signature = appendApplySignature,
  name = "applyDeferred")

# Section: Public Apply Construction

proc apply*(fun: Gvalue, x: Gvalue): Gvalue =
  if fun == nil:
    raiseValueError("apply function cannot be nil")
  if x == nil:
    raiseValueError("apply argument cannot be nil")

  let proto = applyResultProto(fun)
  if proto == nil:
    raiseValueError("apply expects a lambda or local function placeholder, got: " & fun.nodeRepr)

  result = proto.newOneOf
  result.inputs = @[fun, x]
  result.gfunc = gapplyDeferred

# Section: Demonstrations

when isMainModule:
  block:
    let x = toGvalue(3.0)
    let y = toGvalue(2.0)
    let v = localScalar()
    let z = apply(lambda(v, v + v), x * y)
    echo "## z (before eval)"
    echo z.treeRepr
    z.eval
    echo "## z (after eval)"
    echo z.treeRepr
    echo "z = ",z
    let dzdx = z.grad x
    echo "## dzdx (before eval)"
    echo dzdx.treeRepr
    dzdx.eval
    echo "## dzdx (after eval)"
    echo dzdx.treeRepr
    echo "dzdx = ",dzdx
  block:
    let f = local(localScalar())
    let u = localScalar()
    let hof = lambda(f, lambda(u, apply(f, u) + 1.0))
    let v = localScalar()
    let a = toGvalue(2.0)
    let g = lambda(v, a * v)
    let z = apply(apply(hof, g), 3.0)
    echo "## z (before eval)"
    echo z.treeRepr
    z.eval
    echo "## z (after eval)"
    echo z.treeRepr
    echo "z = ",z
    let dzda = z.grad a
    echo "## dzda (before eval)"
    echo dzda.treeRepr
    dzda.eval
    echo "## dzda (after eval)"
    echo dzda.treeRepr
    echo "dzda = ",dzda
  block:
    # Short-depth recursion via Y combinator: F(n) = if n == 0 then base else F(n-1) + step
    let protoArg = localScalar()
    let protoRet = localScalar()
    let fnProto = lambda(protoArg, protoRet)
    let x = local(fnProto)
    let f = local(fnProto)
    let Y = lambda(f, apply(lambda(x, apply(f, apply(x, x))), lambda(x, apply(f, apply(x, x)))))

    let rf = local(localScalar())
    let u = localScalar()
    let v = localScalar()
    let base = toGvalue(1.0)
    let step = toGvalue(1.0)
    let F = lambda(rf, lambda(u,
      cond(equal(u, 0.0), base,
        apply(lambda(v, apply(rf, v) + step), u - 1.0))))

    let z = apply(apply(Y, F), 3.0)
    echo "## Y recursion z (before eval)"
    echo z.treeRepr
    z.eval
    echo "## Y recursion z (after eval)"
    echo z.treeRepr
    echo "Y recursion z = ",z
    let dzdstep = z.grad step
    echo "## Y recursion dzdstep (before eval)"
    echo dzdstep.treeRepr
    dzdstep.eval
    echo "## Y recursion dzdstep (after eval)"
    echo dzdstep.treeRepr
    echo "Y recursion dzdstep = ",dzdstep
    let dzdbase = z.grad base
    echo "## Y recursion dzdbase (before eval)"
    echo dzdbase.treeRepr
    dzdbase.eval
    echo "## Y recursion dzdbase (after eval)"
    echo dzdbase.treeRepr
    echo "Y recursion dzdbase = ",dzdbase
