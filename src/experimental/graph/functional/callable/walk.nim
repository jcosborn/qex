import std/sets
import ../../core
import ../../core/base
import types

type
  CallableFreshnessCtx = object
    active: NodeSet
  CallableWalkCtx = object
    seenCallable: NodeSet
    seenValues: NodeSet
    excludedValues: NodeSet
  CallableWalkResult = object
    maxEpoch: int
    valueDeps: seq[Gvalue]

proc addCallableDep(deps: var seq[Gvalue],
                    dep: Gvalue,
                    label: string) =
  if dep == nil:
    raiseValueError(label & " produced nil dependency")
  deps.add dep

proc producerDeps(v: Gvalue): seq[Gvalue] =
  ## Collect the producer-side deps that decide whether a callable wrapper's
  ## cached binding is still fresh. Wrapper self-bindings are excluded here;
  ## including them would make stale callable captures look current.
  if not (v of Gwrapper):
    return v.collectNodeInputs(iwmDepend)

  let walk = v.gfunc.depWalkForMode(iwmDepend)
  if walk == nil:
    for input in v.inputs:
      result.addCallableDep(input, "callable wrapper raw input")
  else:
    var collected: seq[Gvalue] = @[]
    walk(v, proc(input: Gvalue) =
      collected.addCallableDep(input, "custom callable dependency walk"))
    result = collected

proc freshBound(w: Gwrapper,
                freshness: var CallableFreshnessCtx): Gvalue

proc maxCallableVisibleEpoch(roots: openArray[Gvalue],
                             freshness: var CallableFreshnessCtx): int

proc freshBound(w: Gwrapper,
                freshness: var CallableFreshnessCtx): Gvalue =
  ## Reuse a callable wrapper binding while its producer is fresh.
  if w == nil or w.kind != wkCallable or w.bound == nil:
    return nil
  let node = Gvalue(w)
  # A leaf callable wrapper has no producer graph to refresh, so its last
  # binding is its symbolic value.
  if node.gfunc == nil and node.inputs.len == 0:
    return w.bound
  let key = node.nodeKey
  if key in freshness.active:
    return nil
  freshness.active.incl key
  defer:
    freshness.active.excl key

  # Closure captures remain dependency-visible through the bound callable; they
  # do not make the producer binding itself stale.
  let maxep = maxCallableVisibleEpoch(node.producerDeps, freshness)
  if node.epoch < maxep:
    return nil
  w.bound

proc symbolicBinding(w: Gwrapper,
                     freshness: var CallableFreshnessCtx): Gvalue =
  if w == nil:
    return nil
  case w.kind
  of wkLocal:
    return w.bound
  of wkCallable:
    w.freshBound(freshness)

proc freshCallableBound*(v: Gvalue): Gvalue =
  if not (v of Gwrapper):
    return nil
  let w = Gwrapper(v)
  if w.kind != wkCallable:
    return nil
  var freshness = CallableFreshnessCtx(active: initHashSet[NodeKey]())
  w.freshBound(freshness)

proc symbolicWrapperBinding*(v: Gvalue): Gvalue =
  if not (v of Gwrapper):
    return nil
  var freshness = CallableFreshnessCtx(active: initHashSet[NodeKey]())
  Gwrapper(v).symbolicBinding(freshness)

proc walkCallableVisibleDeps(roots: openArray[Gvalue],
                             seedValues: openArray[Gvalue],
                             freshness: var CallableFreshnessCtx): CallableWalkResult =
  ## Walk graph structure until callable boundaries, then follow only the
  ## binding that is symbolically valid for the boundary's current freshness.
  ## The same walk answers both freshness and apply-cache value-dependency
  ## questions so those contracts cannot drift.
  var ctx = CallableWalkCtx(
    seenCallable: initHashSet[NodeKey](),
    seenValues: initHashSet[NodeKey](),
    excludedValues: initHashSet[NodeKey]())
  for value in seedValues:
    ctx.excludedValues.incl value.nodeKey
  var stack: seq[Gvalue] = @[]
  if roots.len > 0:
    for i in countdown(roots.len - 1, 0):
      stack.add roots[i]

  while stack.len > 0:
    let node = stack[^1]
    stack.setLen(stack.len - 1)

    if node of Glambda:
      let fn = Glambda(node)
      if fn.isResolvedClosure:
        if not ctx.seenCallable.markSeenNode(node):
          continue
        var deps: seq[Gvalue] = @[]
        for binding in fn.env:
          deps.addCallableDep(binding.value, "lambda closure binding")
        if deps.len > 0:
          for i in countdown(deps.len - 1, 0):
            stack.add deps[i]
        continue

    let w =
      if node of Gwrapper:
        Gwrapper(node)
      else:
        nil
    if w != nil:
      if not ctx.seenCallable.markSeenNode(node):
        continue
      var deps: seq[Gvalue] = @[]
      case w.kind
      of wkLocal:
        if result.maxEpoch < node.epoch:
          result.maxEpoch = node.epoch
        deps = node.collectNodeInputs(iwmDepend)
      of wkCallable:
        let boundCallable = w.freshBound(freshness)
        if boundCallable != nil:
          stack.add boundCallable
          continue
        deps = node.producerDeps
      if deps.len > 0:
        for i in countdown(deps.len - 1, 0):
          stack.add deps[i]
      continue

    if not ctx.seenValues.markSeenNode(node):
      continue
    if result.maxEpoch < node.epoch:
      result.maxEpoch = node.epoch
    if node.nodeKey notin ctx.excludedValues:
      result.valueDeps.add node
    let deps = node.collectNodeInputs(iwmDepend)
    if deps.len > 0:
      for i in countdown(deps.len - 1, 0):
        stack.add deps[i]

proc maxCallableVisibleEpoch(roots: openArray[Gvalue],
                             freshness: var CallableFreshnessCtx): int =
  let seedValues: seq[Gvalue] = @[]
  let walk = walkCallableVisibleDeps(roots, seedValues, freshness)
  walk.maxEpoch

proc collectCallableValueDeps*(roots: openArray[Gvalue],
                               seedValues: openArray[Gvalue]): seq[Gvalue] =
  ## Return apply-cache-visible ordinary value dependencies under callable
  ## boundaries. Seed values are excluded from the returned value deps; fresh
  ## callable wrappers expose their bound callable, stale wrappers expose
  ## producer deps, and custom walks must not yield nil dependencies.
  var freshness = CallableFreshnessCtx(active: initHashSet[NodeKey]())
  let walk = walkCallableVisibleDeps(roots, seedValues, freshness)
  walk.valueDeps

method walkHiddenDeps*(v: Gwrapper,
                       mode: InputWalkMode,
                       visit: GnodeVisit) =
  if mode notin {iwmGradSignature, iwmDepend}:
    return
  var freshness = CallableFreshnessCtx(active: initHashSet[NodeKey]())
  let bound = v.symbolicBinding(freshness)
  if bound != nil:
    visit bound
