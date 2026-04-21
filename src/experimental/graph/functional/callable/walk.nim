import ../../core
import types, values

type
  CallableBoundaryMode = enum
    cbmValue, cbmEpoch
  FreshCallableBinding = proc(w: Gwrapper, v: Gvalue): Gvalue {.closure.}
  CallableBoundary = object
    matched: bool
    deps: seq[Gvalue]
  CallableWalkCtx = object
    seenCallable: NodeSet
    seenValues: NodeSet
    excludedValues: NodeSet

proc addDependency(deps: var seq[Gvalue], dep: Gvalue) {.inline.} =
  if dep != nil:
    deps.add dep

proc collectDeps(deps: var seq[Gvalue], values: openArray[Gvalue]) =
  for value in values:
    deps.addDependency(value)

proc collectLambdaEnvDeps(deps: var seq[Gvalue], fn: Glambda) =
  for binding in fn.env:
    deps.addDependency(binding.value)

proc collectInputDeps(deps: var seq[Gvalue], v: Gvalue) =
  deps.collectDeps(v.inputs)

proc collectHiddenDeps(deps: var seq[Gvalue],
                       v: Gvalue,
                       mode: InputWalkMode = iwmDepend) =
  var collected = deps
  v.walkHiddenDeps(mode, proc(dep: Gvalue) = collected.addDependency(dep))
  deps = collected

proc collectGraphDeps(v: Gvalue, deps: var seq[Gvalue]) =
  deps.collectHiddenDeps(v)
  deps.collectInputDeps(v)

proc hasBoundary(w: Gwrapper): bool =
  w != nil and (w.kind == wkCallable or w.retProto != nil)

proc initCallableBoundary(deps: seq[Gvalue] = @[]): CallableBoundary =
  CallableBoundary(matched: true, deps: deps)

proc noCallableBoundary(): CallableBoundary =
  CallableBoundary()

proc wrapperBoundary(w: Gwrapper,
                     v: Gvalue,
                     mode: CallableBoundaryMode,
                     freshBound: FreshCallableBinding): CallableBoundary =
  if not w.hasBoundary:
    return noCallableBoundary()

  var deps: seq[Gvalue] = @[]
  case mode
  of cbmValue:
    if w.kind == wkCallable and freshBound != nil:
      let boundCallable = freshBound(w, v)
      if boundCallable != nil:
        deps.addDependency(boundCallable)
        return initCallableBoundary(deps)
    v.collectGraphDeps(deps)
  of cbmEpoch:
    case w.kind
    of wkCallable:
      deps.addDependency(w.bound)
      deps.collectInputDeps(v)
    of wkLocal:
      v.collectGraphDeps(deps)
  initCallableBoundary(deps)

proc lambdaBoundary(v: Gvalue,
                    mode: CallableBoundaryMode): CallableBoundary =
  discard mode
  let fn = v.resolvedClosure
  if fn == nil:
    return noCallableBoundary()

  var deps: seq[Gvalue] = @[]
  deps.collectLambdaEnvDeps(fn)
  initCallableBoundary(deps)

proc describeCallableBoundary(v: Gvalue,
                              mode: CallableBoundaryMode,
                              freshBound: FreshCallableBinding = nil): CallableBoundary =
  result = v.lambdaBoundary(mode)
  if result.matched:
    return
  result = v.wrapper.wrapperBoundary(v, mode, freshBound)

proc updateMax(value: var int, candidate: int) {.inline.} =
  if value < candidate:
    value = candidate

proc initCallableWalkCtx(): CallableWalkCtx =
  result.seenCallable = initNodeSet()
  result.seenValues = initNodeSet()
  result.excludedValues = initNodeSet()

proc seedExcludedValues(ctx: var CallableWalkCtx,
                        values: openArray[Gvalue]) =
  for value in values:
    if value != nil:
      ctx.excludedValues.inclNode value

proc pushDependencies(stack: var seq[Gvalue], deps: openArray[Gvalue]) =
  if deps.len == 0:
    return
  for i in countdown(deps.len - 1, 0):
    stack.add deps[i]

proc visitIfBound(bound: Gvalue, visit: proc(n: Gvalue) {.closure.}) =
  if bound != nil:
    visit bound

proc maxInputEpoch(v: Gvalue): int =
  if v == nil:
    return 0
  for input in v.inputs:
    if input != nil and result < input.epochOf:
      result = input.epochOf

proc callableDependencyEpoch*(v: Gvalue): int

proc freshBound*(w: Gwrapper, v: Gvalue): Gvalue =
  ## Reuse a callable wrapper binding while its producer is fresh.
  if not w.isWrapperKind(wkCallable):
    return nil
  if w.bound == nil:
    return nil
  var maxep = v.maxInputEpoch
  let boundEpoch = callableDependencyEpoch(w.bound)
  if maxep < boundEpoch:
    maxep = boundEpoch
  if v.epochOf < maxep:
    return nil
  w.bound

proc symbolicBinding*(w: Gwrapper, v: Gvalue): Gvalue =
  if w == nil:
    return nil
  case w.kind
  of wkLocal:
    return w.bound
  of wkCallable:
    result = w.freshBound(v)
    if result == nil and w.inputs.len == 0:
      result = w.bound

proc freshCallableBound*(v: Gvalue): Gvalue =
  v.callableWrapper.freshBound(v)

proc symbolicWrapperBinding*(v: Gvalue): Gvalue =
  v.wrapper.symbolicBinding(v)

proc visitWrapperHiddenDeps*(v: Gvalue, visit: proc(n: Gvalue) {.closure.}) =
  visitIfBound(v.symbolicWrapperBinding, visit)

proc walkCallableGraph(v: Gvalue,
                       ctx: var CallableWalkCtx,
                       mode: CallableBoundaryMode,
                       visitBoundary: proc(v: Gvalue) {.closure.} = nil,
                       visitValue: proc(v: Gvalue) {.closure.} = nil) =
  ## Walk graph structure until a callable boundary, then continue through the
  ## dependencies exposed at that boundary.
  var stack: seq[Gvalue] = @[]
  if v != nil:
    stack.add v
  while stack.len > 0:
    let node = stack[^1]
    stack.setLen(stack.len - 1)
    let boundary = node.describeCallableBoundary(
      mode,
      proc(w: Gwrapper, owner: Gvalue): Gvalue =
        w.freshBound(owner))
    if boundary.matched:
      if not ctx.seenCallable.markSeenNode(node):
        continue
      if visitBoundary != nil:
        visitBoundary(node)
      stack.pushDependencies(boundary.deps)
      continue
    if not ctx.seenValues.markSeenNode(node):
      continue
    if visitValue != nil:
      visitValue(node)
    stack.pushDependencies(node.collectNodeInputs(iwmDepend))

proc walkCallableValueDeps(v: Gvalue,
                           ctx: var CallableWalkCtx,
                           visitCollectedValue: proc(v: Gvalue) {.closure.}) =
  ## Seeded values are traversed but not reported.
  let excludedValues = ctx.excludedValues
  let collectValue = proc(node: Gvalue) =
    if visitCollectedValue != nil and not excludedValues.containsNode(node):
      visitCollectedValue(node)
  walkCallableGraph(v, ctx, cbmValue, visitValue = collectValue)

proc dependencyEpoch(ctx: var CallableWalkCtx, v: Gvalue): int =
  var epoch = 0
  let markBoundary = proc(node: Gvalue) =
    if node.isLocalWrapper:
      epoch.updateMax node.epochOf
  let markValue = proc(node: Gvalue) =
    epoch.updateMax node.epochOf
  walkCallableGraph(
    v,
    ctx,
    cbmEpoch,
    visitBoundary = markBoundary,
    visitValue = markValue)
  result = epoch

proc collectCallableValueDeps*(roots: openArray[Gvalue],
                               seedValues: openArray[Gvalue]): seq[Gvalue] =
  var ctx = initCallableWalkCtx()
  var collected: seq[Gvalue] = @[]
  ctx.seedExcludedValues(seedValues)
  let collectValue = proc(node: Gvalue) =
    collected.add node
  for root in roots:
    walkCallableValueDeps(root, ctx, collectValue)
  result = collected

proc callableDependencyEpoch*(v: Gvalue): int =
  var ctx = initCallableWalkCtx()
  ctx.dependencyEpoch(v)

proc trackDependencyEpoch(ctx: var CallableWalkCtx,
                          value: Gvalue,
                          freshEpoch: var int) =
  if value == nil:
    return
  freshEpoch.updateMax ctx.dependencyEpoch(value)

proc collectDependencyEpoch*(values: openArray[Gvalue]): int =
  var ctx = initCallableWalkCtx()
  for value in values:
    ctx.trackDependencyEpoch(value, result)

method walkHiddenDeps*(v: Gwrapper,
                       mode: InputWalkMode,
                       visit: GnodeVisit) =
  if mode notin {iwmGradSignature, iwmDepend}:
    return
  Gvalue(v).visitWrapperHiddenDeps(visit)
