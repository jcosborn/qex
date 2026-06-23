import std/[sets, tables]
import base, traverse
import cond

proc resetGradCache*(grt: GraphRuntime, stats = true) =
  grt.gradCacheByOutput = initTable[NodeId, GradCacheEntry]()
  if stats:
    grt.gradCacheStats = GradCacheStats()

proc findGrad*(input: Gvalue, output: Gvalue): Gvalue =
  let grt = sharedGraphRuntime([input, output], "findGrad")
  let entry = grt.gradCacheByOutput.getOrDefault(output.stableNodeId)
  if entry == nil:
    return nil
  if entry.revision != grt.symbolicRevision:
    return nil
  entry.adjoints.getOrDefault(input.stableNodeId)

type
  GradBuildPlan = object
    needsGradient: NodeSet
    order: seq[Gvalue]
  GradBuildContext = ref object
    dep: Gvalue
    x: Gvalue
    rootSeed: Gvalue
    cache: GradCacheEntry
    plan: GradBuildPlan
    pendingAdjoints: Table[NodeId, Gvalue]

proc requireGradResultShape(label: string, x: Gvalue, result: Gvalue) =
  if not x.copyCompatible(result):
    raiseValueError(
      label & " result is not input-shaped" &
      "\ninput: " & x.nodeRepr &
      "\nresult: " & result.nodeRepr)

proc currentAdjoint(ctx: GradBuildContext, node: Gvalue): Gvalue =
  let nodeId = node.stableNodeId
  result = ctx.pendingAdjoints.getOrDefault(nodeId)
  if result != nil:
    return
  result = ctx.cache.adjoints.getOrDefault(nodeId)

proc addGradContribution(ctx: GradBuildContext,
                         input: Gvalue,
                         contrib: Gvalue) =
  if contrib == nil:
    raiseValueError("gradient contribution cannot be nil")
  let inputId = input.stableNodeId
  let existing = ctx.currentAdjoint(input)
  if existing != nil:
    ctx.pendingAdjoints[inputId] = input.addLike(existing, contrib)
  else:
    ctx.pendingAdjoints[inputId] = contrib

proc buildGradPlan(dep: Gvalue, x: Gvalue): GradBuildPlan =
  var plan = GradBuildPlan(needsGradient: initHashSet[NodeKey]())
  var visited = initHashSet[NodeKey]()
  var active = initHashSet[NodeKey]()

  proc markGradientNode(node: Gvalue): bool =
    let key = node.nodeKey
    if key in visited:
      return key in plan.needsGradient
    if key in active:
      raiseError("cycle detected while building gradient plan:\n" & node.nodeRepr)

    active.incl key
    var need = node.nodeKey == x.nodeKey
    node.walkInputView(iwmBackward, proc(child: Gvalue) =
      if markGradientNode(child):
        need = true)
    active.excl key

    visited.incl key
    if need:
      plan.needsGradient.incl key
      plan.order.add node
    need

  discard markGradientNode(dep)
  result = plan

proc buildInputContribution(graphFunc: Gfunc,
                            upstream: Gvalue,
                            node: Gvalue,
                            inputIndex: int,
                            input: Gvalue): Gvalue =
  # A conditional upstream is split so each branch gets its own guarded backward
  # graph. Static-zero branches are skipped here, so an inactive branch never
  # builds a backward graph. `buildInputContribution` itself never returns nil,
  # so a nil branch result can only mean the static-zero short-circuit.
  if upstream != nil and upstream.isCondNode:
    let parts = upstream.condParts
    proc branchContribution(branch: Gvalue): Gvalue =
      if branch.isStaticZeroLeaf:
        return nil
      buildInputContribution(graphFunc, branch, node, inputIndex, input)
    let whenTrue = branchContribution(parts.whenTrue)
    let whenFalse = branchContribution(parts.whenFalse)
    if whenTrue == nil and whenFalse == nil:
      return input.zeroLike
    return newCondNode(
      parts.selector,
      if whenTrue == nil: input.zeroLike else: whenTrue,
      if whenFalse == nil: input.zeroLike else: whenFalse)

  let backward = graphFunc.backward
  if backward == nil:
    raiseError(
      node.nodeRepr & ":" & $inputIndex & ":" & input.nodeRepr &
      ": backward undefined")
  result = backward(upstream, node, inputIndex, input)
  if result == nil:
    raiseValueError(
      graphFunc.name & " backward returned nil for input " & $inputIndex &
      ":\n" & node.nodeRepr)

proc executeGradPlan(ctx: GradBuildContext): Gvalue =
  for j in countdown(ctx.plan.order.high, 0):
    let node = ctx.plan.order[j]
    let upstream =
      if node.nodeKey == ctx.dep.nodeKey:
        ctx.rootSeed
      else:
        let adjoint = ctx.currentAdjoint(node)
        if adjoint == nil:
          continue
        adjoint
    let graphFunc = node.gfunc
    # The gradient target leaf x is in plan.order but has no backward function.
    if graphFunc == nil:
      continue
    var i = 0
    node.walkInputView(iwmBackward, proc(input: Gvalue) =
      let inputIndex = i
      inc i
      if input.nodeKey notin ctx.plan.needsGradient:
        return
      if ctx.cache.adjoints.hasKey(input.stableNodeId):
        return
      let contribution = buildInputContribution(
        graphFunc, upstream, node, inputIndex, input)
      ctx.addGradContribution(input, contribution))
  let adjoint = ctx.currentAdjoint(ctx.x)
  if adjoint != nil:
    return adjoint
  ctx.x.zeroLike

proc buildGradAdjoint(dep: Gvalue,
                      x: Gvalue,
                      rootSeed: Gvalue,
                      cache: GradCacheEntry,
                      label: string): tuple[
                        adjoint: Gvalue,
                        pendingAdjoints: Table[NodeId, Gvalue]] =
  result.pendingAdjoints = initTable[NodeId, Gvalue]()
  if dep.nodeKey == x.nodeKey:
    result.adjoint = rootedUpstream(rootSeed, x)
  else:
    let plan = buildGradPlan(dep, x)
    if dep.nodeKey notin plan.needsGradient:
      result.adjoint = x.zeroLike
    else:
      let ctx = GradBuildContext(
        dep: dep,
        x: x,
        rootSeed: rootSeed,
        cache: cache,
        plan: plan,
        pendingAdjoints: initTable[NodeId, Gvalue]())
      result.adjoint = ctx.executeGradPlan
      result.pendingAdjoints = ctx.pendingAdjoints
  requireGradResultShape(label, x, result.adjoint)

proc gradImpl(dep: Gvalue, x: Gvalue): Gvalue =
  discard sharedGraphRuntime([dep, x], "grad")
  let grt = dep.runtime
  let revision = grt.symbolicRevision
  let depId = dep.stableNodeId
  let inputId = x.stableNodeId
  let cache =
    if grt.gradCacheByOutput.hasKey(depId):
      let entry = grt.gradCacheByOutput[depId]
      if entry.revision == revision:
        grt.gradCacheStats.revisionHits.inc
        if entry.adjoints.hasKey(inputId):
          grt.gradCacheStats.directHits.inc
          return entry.adjoints[inputId]
      else:
        grt.gradCacheStats.revisionMisses.inc
        grt.gradCacheStats.invalidations.inc
        entry.revision = revision
        entry.adjoints = initTable[NodeId, Gvalue]()
      entry
    else:
      grt.gradCacheStats.revisionMisses.inc
      let entry = GradCacheEntry(
        revision: revision,
        adjoints: initTable[NodeId, Gvalue]())
      grt.gradCacheByOutput[depId] = entry
      entry

  let built = buildGradAdjoint(dep, x, nil, cache, "grad")
  result = built.adjoint
  for nodeId, adjoint in built.pendingAdjoints:
    cache.adjoints[nodeId] = adjoint
  cache.adjoints[inputId] = result

proc gradSeeded*(dep: Gvalue, x: Gvalue, seed: Gvalue): Gvalue =
  discard sharedGraphRuntime([dep, x, seed], "gradSeeded")
  let seedProto = dep.zeroLike
  if not seedProto.copyCompatible(seed):
    raiseValueError(
      "gradSeeded seed is incompatible with output cotangent" &
      "\noutput: " & dep.nodeRepr &
      "\noutput cotangent prototype: " & seedProto.nodeRepr &
      "\nseed: " & seed.nodeRepr)

  result = buildGradAdjoint(
    dep,
    x,
    seed,
    GradCacheEntry(adjoints: initTable[NodeId, Gvalue]()),
    "gradSeeded").adjoint

proc grad*[T: Gvalue](dep: Gvalue, x: T): T =
  T(gradImpl(dep, x))
