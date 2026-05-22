import std/[sets, tables]
import ../../core
import ../../core/base
import ../../core/traverse
import ../../scalar
import ../callable, ../callable/walk, ../closure/normalize

proc reducedApplyExpr(v: Gvalue): Gvalue
proc isApplyPartialNode(v: Gvalue): bool
proc isDeferredApplyNode(v: Gvalue): bool

proc applyDeferredBackwardTarget(zb: Gvalue,
                                 z: Gvalue,
                                 target: Gvalue,
                                 dep: Gvalue): Gvalue

type
  ApplyDeps = object
    ## Raw fun/arg inputs plus callable-visible value deps. Dependency walks and
    ## the apply cache key both derive from this record so they cannot drift.
    fun: Gvalue
    arg: Gvalue
    callableVisibleDeps: seq[Gvalue]
  ApplyPartialView = object
    base: Gvalue
    target: Gvalue

proc requireApplyInputs(v: Gvalue): tuple[fun: Gvalue, arg: Gvalue] =
  v.requireInputCountExactly(2, "apply node")
  (
    fun: v.requireNodeInput(0, "apply node", "function"),
    arg: v.requireNodeInput(1, "apply node", "argument"))

proc requireApplyPartialView(v: Gvalue,
                             label = "applyPartialDeferred node"): ApplyPartialView =
  v.requireInputCountExactly(2, label)
  result.base = v.requireNodeInput(0, label, "base")
  result.target = v.requireNodeInput(1, label, "target")

proc callableResultCompatible(expected: Gvalue,
                              actual: Gvalue,
                              remainingDepth: int): bool =
  if remainingDepth <= 0:
    raiseValueError("callable argument compatibility exceeded depth limit")
  if expected.copyCompatible(actual):
    return true

  # Callable shells are independently materialized, so compare the value shape
  # they eventually return instead of accepting any callable-like pair.
  let expectedResult = expected.bindingResultProto
  let actualResult = actual.bindingResultProto
  if expectedResult == nil or actualResult == nil:
    return false
  callableResultCompatible(expectedResult, actualResult, remainingDepth - 1)

proc lambdaArgCompatible(param: Gvalue, arg: Gvalue): bool =
  let expectedCallableResult = param.bindingResultProto
  if expectedCallableResult != nil:
    let actualCallableResult = arg.bindingResultProto
    if actualCallableResult == nil:
      return false
    return callableResultCompatible(
      expectedCallableResult,
      actualCallableResult,
      param.runtime.lambdaResolveDepthLimit)
  param.copyCompatible(arg)

proc collectApplyDeps(v: Gvalue): ApplyDeps =
  let inputs = v.requireApplyInputs
  result.fun = inputs.fun
  result.arg = inputs.arg
  result.callableVisibleDeps = collectCallableValueDeps(
    [result.fun, result.arg],
    [result.fun, result.arg])

proc cacheInputKeys(deps: ApplyDeps): seq[NodeId] =
  result.add deps.fun.stableNodeId
  result.add deps.arg.stableNodeId
  var seen = initHashSet[NodeKey]()
  seen.incl deps.fun.nodeKey
  seen.incl deps.arg.nodeKey
  for dep in deps.callableVisibleDeps:
    if seen.markSeenNode(dep):
      result.add dep.stableNodeId

proc walkApplyDeps(v: Gvalue,
                   visit: proc(n: Gvalue) {.closure.}) =
  let deps = v.collectApplyDeps
  visit deps.fun
  visit deps.arg
  for dep in deps.callableVisibleDeps:
    visit dep

proc appendApplySignature(v: Gvalue, tokens: var seq[GradSigToken]) =
  # The signature records callable identity tokens only; value deps are already
  # exposed through `walkApplyDeps`.
  let applyInputs = v.requireApplyInputs
  for input in [applyInputs.fun, applyInputs.arg]:
    let token = symbolicCallableToken(input)
    if token != 0:
      tokens.add GradSigToken(kind: gstCallable, key: token)

proc ensureReduction(v: Gvalue): ApplyCacheEntry =
  let deps = v.collectApplyDeps
  let fn = resolveLambda(deps.fun)
  if fn == nil:
    raiseUnresolvedValueError("deferred apply unresolved at eval: " & deps.fun.nodeRepr)
  if not fn.param.lambdaArgCompatible(deps.arg):
    raiseValueError(
      "apply argument is incompatible with lambda parameter" &
      "\nparameter: " & fn.param.nodeRepr &
      "\nargument: " & deps.arg.nodeRepr &
      "\napply: " & v.nodeRepr)
  let key = ApplyCacheKey(
    callableKey: fn.stableNodeId,
    inputKeys: deps.cacheInputKeys)

  let grt = v.runtime
  let nodeId = v.stableNodeId
  if not grt.applyCacheByNode.hasKey(nodeId):
    grt.applyCacheByNode[nodeId] = ApplyCacheEntry()
  result = grt.applyCacheByNode[nodeId]
  if result.reduced != nil and result.key == key:
    grt.applyCacheStats.reduceHits.inc
    return
  grt.applyCacheStats.reduceMisses.inc

  let reduced = instantiateNormalizedBody(fn, deps.arg)
  if reduced == nil:
    raiseValueError("deferred apply reduction produced nil body")
  result.key = key
  result.reduced = reduced
  result.partials = initTable[NodeId, Gvalue]()

proc reducedApplyExpr(v: Gvalue): Gvalue =
  let entry = v.ensureReduction
  entry.reduced

proc materializeConcreteApplyPartial(z: Gvalue,
                                     reduced: Gvalue,
                                     target: Gvalue): Gvalue =
  let grt = z.runtime
  inc grt.applyGradPrepareDepth
  defer:
    dec grt.applyGradPrepareDepth
  if grt.applyGradPrepareDepth > grt.applyGradPrepareDepthLimit:
    raiseValueError(
      "apply partial materialization exceeded depth limit " &
      $grt.applyGradPrepareDepthLimit &
      "\napply: " & z.nodeRepr &
      "\ntarget: " & target.nodeRepr)
  reduced.grad(target)

proc ensurePartial(z: Gvalue, target: Gvalue): Gvalue =
  ## Resolve the target-specific partial for an apply expression. Nested
  ## apply-partial nodes chain through their base partial, then differentiate the
  ## materialized partial with respect to the requested target.
  let grt = z.runtime
  let targetGrt = target.runtime
  # Apply-partial cache keys are runtime-local stable ids.
  if grt != targetGrt:
    raiseValueError("apply partial mixes multiple graph runtimes")
  if z.isApplyPartialNode:
    let view = z.requireApplyPartialView
    return view.base.ensurePartial(view.target).grad(target)
  let entry = z.ensureReduction
  let targetId = target.stableNodeId
  if entry.partials.hasKey(targetId):
    grt.applyCacheStats.partialHits.inc
    return entry.partials[targetId]
  grt.applyCacheStats.partialMisses.inc
  result =
    if target.bindingResultProto != nil:
      target.zeroLike
    else:
      z.materializeConcreteApplyPartial(entry.reduced, target)
  entry.partials[targetId] = result

proc applyDeferredf(v: Gvalue) =
  let reduced = v.reducedApplyExpr
  let callable = resolveLambda(reduced)
  if callable != nil:
    v.valCopy callable
    return
  discard reduced.eval
  v.valCopy reduced.resolveCallableChain(crmReduced)

proc applyPartialDeferredf(v: Gvalue) =
  let view = v.requireApplyPartialView
  let expr = view.base.ensurePartial(view.target)
  discard expr.eval
  v.valCopy expr

proc walkApplyPartialMaterializationDeps(base: Gvalue,
                                         visit: proc(n: Gvalue) {.closure.}) =
  ## Apply-partial eval traversal is apply-owned: deferred nodes expose
  ## depend-mode deps and stay lazy, while concrete nodes expose eval deps plus
  ## depend-mode hidden deps for callable captures.
  var seen = initHashSet[NodeKey]()

  proc walkNode(node: Gvalue) =
    if not seen.markSeenNode(node):
      return
    if node.isDeferredApplyNode:
      node.walkDeps(iwmDepend, walkNode)
      return
    node.walkDeps(iwmEval, walkNode)
    node.walkHiddenDeps(iwmDepend, walkNode)
    visit node

  walkNode(base)

proc applyPartialDeferredWalkMaterializationInputs(
    v: Gvalue,
    visit: proc(n: Gvalue) {.closure.}) =
  walkApplyPartialMaterializationDeps(v.requireApplyPartialView.base, visit)

proc applyPartialDeferredWalkGraphInputs(v: Gvalue,
                                         visit: proc(n: Gvalue) {.closure.}) =
  visit(v.requireApplyPartialView.base)

# Partial application nodes evaluate through deferred base work, but signature
# and depend walks expose only the base. The raw target input is handled by the
# backward target hook when a target-specific partial is materialized.
let gapplyPartialDeferred = newGfunc(
  forward = applyPartialDeferredf,
  depWalks = GdepWalks(
    eval: applyPartialDeferredWalkMaterializationInputs,
    gradSignature: applyPartialDeferredWalkGraphInputs,
    depend: applyPartialDeferredWalkGraphInputs),
  backwardTarget = applyDeferredBackwardTarget,
  name = "applyPartialDeferred")
# Deferred apply exposes callable, argument, and capture deps in every traversal
# mode. `reduceCallable` materializes the resolved lambda body only on demand.
let gapplyDeferred = newGfunc(
  forward = applyDeferredf,
  depWalks = GdepWalks(
    eval: walkApplyDeps,
    gradSignature: walkApplyDeps,
    depend: walkApplyDeps),
  backwardTarget = applyDeferredBackwardTarget,
  reduceCallable = reducedApplyExpr,
  signature = appendApplySignature,
  name = "applyDeferred")

proc isApplyPartialNode(v: Gvalue): bool =
  v.gfunc == gapplyPartialDeferred

proc isDeferredApplyNode(v: Gvalue): bool =
  v.gfunc == gapplyDeferred or v.isApplyPartialNode

proc applyDeferredBackwardTarget(zb: Gvalue,
                                 z: Gvalue,
                                 target: Gvalue,
                                 dep: Gvalue): Gvalue =
  discard dep
  if target.bindingResultProto != nil:
    return nil
  if zb != nil and zb.bindingResultProto != nil:
    raiseValueError(
      "apply partial cannot scale callable upstream gradient" &
      "\napply: " & z.nodeRepr &
      "\ntarget: " & target.nodeRepr &
      "\nupstream: " & zb.nodeRepr)
  result = graphNode(
    target.newOneOf,
    @[z, target],
    gapplyPartialDeferred,
    "applyPartialDeferred")
  if zb != nil:
    result = result.scaleLike zb

proc apply*(fun: Gvalue, x: Gvalue): Gvalue =
  let proto = applyResultProto(fun)
  if proto == nil:
    raiseValueError("apply expects a callable value or callable placeholder, got: " & fun.nodeRepr)
  let checkedArg = x.requireGraphValue("apply argument")
  let direct = resolveDirectLambda(fun)
  if direct != nil and not direct.param.lambdaArgCompatible(checkedArg):
    raiseValueError(
      "apply argument is incompatible with lambda parameter" &
      "\nparameter: " & direct.param.nodeRepr &
      "\nargument: " & checkedArg.nodeRepr)

  graphNode(proto.newOneOf, @[fun, checkedArg], gapplyDeferred, "apply")

proc applyLiteralArg(fun: Gvalue,
                     value: float): Gvalue =
  let direct = resolveDirectLambda(fun)
  if direct != nil:
    return scalarLeafLike(direct.param, value)
  scalarLeafLike(fun, value)

proc applyLiteralArg(fun: Gvalue,
                     value: int): Gvalue =
  let direct = resolveDirectLambda(fun)
  if direct != nil:
    return numericLeafLike(direct.param, value)
  intLeafLike(fun, value)

proc apply*(fun: Gvalue,
            x: float): Gvalue =
  apply(fun, applyLiteralArg(fun, x))

proc apply*(fun: Gvalue,
            x: int): Gvalue =
  apply(fun, applyLiteralArg(fun, x))
