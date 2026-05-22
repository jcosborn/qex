import std/[sets, tables]
import base, traverse
import grad_signature, grad_cache

type
  GradBuildPlan = object
    visited: NodeSet
    needsGradient: NodeSet
    order: seq[Gvalue]
  GradBuildContext = object
    dep: Gvalue
    x: Gvalue
    cache: GradCacheEntry
    plan: GradBuildPlan
    pendingAdjoints: Table[NodeId, Gvalue]
    pendingCompleteAdjoints: HashSet[NodeId]
    pendingInputContributions: Table[NodeId, Table[int, Gvalue]]
    pendingExpandedInputs: Table[NodeId, HashSet[int]]
    pendingExpandedTargets: Table[NodeId, HashSet[NodeId]]
    usedBackwardTarget: bool

proc initGradBuildPlan(): GradBuildPlan =
  result.visited = initHashSet[NodeKey]()
  result.needsGradient = initHashSet[NodeKey]()

proc initGradBuildContext(dep: Gvalue,
                          x: Gvalue,
                          cache: GradCacheEntry): GradBuildContext =
  result.dep = dep
  result.x = x
  result.cache = cache
  result.plan = initGradBuildPlan()
  result.pendingAdjoints = initTable[NodeId, Gvalue]()
  result.pendingCompleteAdjoints = initHashSet[NodeId]()
  result.pendingInputContributions = initTable[NodeId, Table[int, Gvalue]]()
  result.pendingExpandedInputs = initTable[NodeId, HashSet[int]]()
  result.pendingExpandedTargets = initTable[NodeId, HashSet[NodeId]]()

proc containsInput(table: Table[NodeId, HashSet[int]],
                   nodeId: NodeId,
                   inputIndex: int): bool =
  table.hasKey(nodeId) and inputIndex in table[nodeId]

proc containsTarget(table: Table[NodeId, HashSet[NodeId]],
                    nodeId: NodeId,
                    targetId: NodeId): bool =
  table.hasKey(nodeId) and targetId in table[nodeId]

proc markInput(table: var Table[NodeId, HashSet[int]],
               nodeId: NodeId,
               inputIndex: int) =
  var inputs =
    if table.hasKey(nodeId): table[nodeId]
    else: initHashSet[int]()
  inputs.incl inputIndex
  table[nodeId] = inputs

proc markTarget(table: var Table[NodeId, HashSet[NodeId]],
                nodeId: NodeId,
                targetId: NodeId) =
  var targets =
    if table.hasKey(nodeId): table[nodeId]
    else: initHashSet[NodeId]()
  targets.incl targetId
  table[nodeId] = targets

proc targetExpanded(ctx: GradBuildContext,
                    nodeId: NodeId,
                    targetId: NodeId): bool =
  ctx.pendingExpandedTargets.containsTarget(nodeId, targetId) or
    ctx.cache.expandedTargets.containsTarget(nodeId, targetId)

proc currentAdjoint(ctx: GradBuildContext, node: Gvalue): Gvalue =
  let nodeId = node.stableNodeId
  if ctx.pendingAdjoints.hasKey(nodeId):
    return ctx.pendingAdjoints[nodeId]
  nil

proc inputContribution(table: Table[NodeId, Table[int, Gvalue]],
                       nodeId: NodeId,
                       inputIndex: int): Gvalue =
  if not table.hasKey(nodeId):
    return nil
  if not table[nodeId].hasKey(inputIndex):
    return nil
  table[nodeId][inputIndex]

proc putInputContribution(table: var Table[NodeId, Table[int, Gvalue]],
                          nodeId: NodeId,
                          inputIndex: int,
                          contribution: Gvalue) =
  var inputContributions =
    if table.hasKey(nodeId): table[nodeId]
    else: initTable[int, Gvalue]()
  inputContributions[inputIndex] = contribution
  table[nodeId] = inputContributions

proc addGradContribution(ctx: var GradBuildContext,
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

proc nodeNeedsGradient(plan: GradBuildPlan, node: Gvalue): bool =
  node.nodeKey in plan.needsGradient

proc buildGradPlan(dep: Gvalue, x: Gvalue): GradBuildPlan =
  var plan = initGradBuildPlan()
  var active = initHashSet[NodeKey]()

  proc markGradientNode(node: Gvalue): bool =
    let key = node.nodeKey
    if key in plan.visited:
      return plan.nodeNeedsGradient(node)
    if key in active:
      raiseError("cycle detected while building gradient plan:\n" & node.nodeRepr)

    active.incl key
    var need = sameNode(node, x)
    node.walkDeps(iwmDepend, proc(child: Gvalue) =
      if markGradientNode(child):
        need = true)
    active.excl key

    plan.visited.incl key
    if need:
      plan.needsGradient.incl key
      plan.order.add node
    need

  discard markGradientNode(dep)
  result = plan

proc propagateInputGradients(ctx: var GradBuildContext,
                             node: Gvalue,
                             upstream: Gvalue) =
  let graphFunc = node.gfunc
  if graphFunc == nil:
    return
  let nodeId = node.stableNodeId
  # `backward` only propagates along raw `inputs`; hidden deps must come from
  # traversal or `backwardTarget`.
  for inputIndex in 0..<node.inputs.len:
    let input = node.inputs[inputIndex]
    if not ctx.plan.nodeNeedsGradient(input):
      continue
    if ctx.pendingExpandedInputs.containsInput(nodeId, inputIndex):
      continue
    if ctx.cache.expandedInputs.containsInput(nodeId, inputIndex):
      let contribution =
        ctx.cache.inputContributions.inputContribution(nodeId, inputIndex)
      if contribution == nil:
        raiseError(
          node.nodeRepr & ":" & $inputIndex & ":" & input.nodeRepr &
          ": cached gradient input contribution missing")
      ctx.addGradContribution(input, contribution)
      ctx.pendingExpandedInputs.markInput(nodeId, inputIndex)
      continue
    let backward = graphFunc.backward
    if backward == nil:
      raiseError(
        node.nodeRepr & ":" & $inputIndex & ":" & input.nodeRepr &
        ": backward undefined")
    let contribution = backward(upstream, node, inputIndex, ctx.dep)
    if contribution == nil:
      raiseValueError(
        graphFunc.name & " backward returned nil for input " & $inputIndex &
        ":\n" & node.nodeRepr)
    ctx.addGradContribution(input, contribution)
    ctx.pendingInputContributions.putInputContribution(
      nodeId, inputIndex, contribution)
    ctx.pendingExpandedInputs.markInput(nodeId, inputIndex)

proc accumulateNodeGradient(ctx: var GradBuildContext, node: Gvalue) =
  let isRoot = sameNode(node, ctx.dep)
  let upstream =
    if isRoot:
      nil
    else:
      let adjoint = ctx.currentAdjoint(node)
      if adjoint == nil:
        return
      adjoint
  # `backwardTarget` overrides ordinary backward propagation for this node.
  # A nil contribution is an explicit zero contribution for the target.
  let graphFunc = node.gfunc
  if graphFunc != nil and graphFunc.backwardTarget != nil:
    ctx.usedBackwardTarget = true
    let nodeId = node.stableNodeId
    let targetId = ctx.x.stableNodeId
    if ctx.targetExpanded(nodeId, targetId):
      return
    let contribution = graphFunc.backwardTarget(upstream, node, ctx.x, ctx.dep)
    if contribution != nil:
      ctx.addGradContribution(ctx.x, contribution)
    ctx.pendingExpandedTargets.markTarget(nodeId, targetId)
    return
  ctx.propagateInputGradients(node, upstream)

proc executeGradPlan(ctx: var GradBuildContext): Gvalue =
  for j in countdown(ctx.plan.order.high, 0):
    ctx.accumulateNodeGradient(ctx.plan.order[j])
  let adjoint = ctx.currentAdjoint(ctx.x)
  if adjoint != nil:
    return adjoint
  ctx.x.zeroLike

proc mergeInputs(dst: var Table[NodeId, HashSet[int]],
                 src: Table[NodeId, HashSet[int]]) =
  for nodeId, inputs in src:
    var merged =
      if dst.hasKey(nodeId): dst[nodeId]
      else: initHashSet[int]()
    for inputIndex in inputs:
      merged.incl inputIndex
    dst[nodeId] = merged

proc mergeTargets(dst: var Table[NodeId, HashSet[NodeId]],
                  src: Table[NodeId, HashSet[NodeId]]) =
  for nodeId, targets in src:
    var merged =
      if dst.hasKey(nodeId): dst[nodeId]
      else: initHashSet[NodeId]()
    for targetId in targets:
      merged.incl targetId
    dst[nodeId] = merged

proc mergeInputContributions(
    dst: var Table[NodeId, Table[int, Gvalue]],
    src: Table[NodeId, Table[int, Gvalue]]) =
  for nodeId, inputContributions in src:
    var merged =
      if dst.hasKey(nodeId): dst[nodeId]
      else: initTable[int, Gvalue]()
    for inputIndex, contribution in inputContributions:
      merged[inputIndex] = contribution
    dst[nodeId] = merged

proc markCompleteAdjoints(ctx: var GradBuildContext) =
  if ctx.usedBackwardTarget:
    ctx.pendingCompleteAdjoints.incl ctx.x.stableNodeId
    return
  for nodeId in ctx.pendingAdjoints.keys:
    ctx.pendingCompleteAdjoints.incl nodeId

proc commitGradBuild(ctx: GradBuildContext) =
  for nodeId, adjoint in ctx.pendingAdjoints:
    if nodeId in ctx.pendingCompleteAdjoints or
        nodeId notin ctx.cache.completeAdjoints:
      ctx.cache.adjoints[nodeId] = adjoint
  for nodeId in ctx.pendingCompleteAdjoints:
    ctx.cache.completeAdjoints.incl nodeId
  ctx.cache.inputContributions.mergeInputContributions(
    ctx.pendingInputContributions)
  ctx.cache.expandedInputs.mergeInputs(ctx.pendingExpandedInputs)
  ctx.cache.expandedTargets.mergeTargets(ctx.pendingExpandedTargets)

proc gradIsolated*(dep: Gvalue, x: Gvalue): Gvalue

proc gradImpl(dep: Gvalue, x: Gvalue): Gvalue =
  discard requireSameGraphRuntime(dep, x, "grad", "output", "input")
  let signature = dep.buildGradSignature
  let cache = dep.ensureGradCacheEntry()
  let grt = dep.runtime
  let inputId = x.stableNodeId

  if cache.hasSignature and cache.signature == signature:
    grt.gradCacheStats.signatureHits.inc
    if inputId in cache.completeAdjoints and cache.adjoints.hasKey(inputId):
      grt.gradCacheStats.directHits.inc
      return cache.adjoints[inputId]
    grt.gradCacheStats.directMisses.inc
  else:
    grt.gradCacheStats.signatureMisses.inc
    if cache.hasSignature:
      grt.gradCacheStats.invalidations.inc
    cache.resetGradCacheEntry(signature)

  if sameNode(dep, x):
    result = x.oneLike
    cache.adjoints[inputId] = result
    cache.completeAdjoints.incl inputId
    return

  var ctx = initGradBuildContext(dep, x, cache)
  ctx.plan = buildGradPlan(dep, x)
  if ctx.plan.nodeNeedsGradient(dep):
    result = ctx.executeGradPlan
  else:
    result = x.zeroLike
  ctx.pendingAdjoints[inputId] = result
  ctx.markCompleteAdjoints()
  ctx.commitGradBuild()

proc grad*(dep: Gvalue, x: Gvalue): Gvalue =
  gradImpl(dep, x)

proc grad*[T: Gvalue](dep: Gvalue, x: T): T =
  T(gradImpl(dep, x))

proc gradIsolated*(dep: Gvalue, x: Gvalue): Gvalue =
  let grt = requireSameGraphRuntime(dep, x, "gradIsolated", "output", "input")
  withIsolatedGradCache(grt, proc(): Gvalue =
    gradImpl(dep, x))

proc gradIsolated*[T: Gvalue](dep: Gvalue, x: T): T =
  T(gradIsolated(dep, Gvalue(x)))
