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
    plan: GradBuildPlan
    contribs: NodeTable[Gvalue]

proc initGradBuildPlan(): GradBuildPlan =
  result.visited = initHashSet[NodeKey]()
  result.needsGradient = initHashSet[NodeKey]()

proc initGradBuildContext(dep: Gvalue, x: Gvalue): GradBuildContext =
  result.dep = dep
  result.x = x
  result.plan = initGradBuildPlan()
  result.contribs = initTable[NodeKey, Gvalue]()

proc addGradContribution(ctx: var GradBuildContext,
                         input: Gvalue,
                         contrib: Gvalue) =
  if contrib == nil:
    raiseValueError("gradient contribution cannot be nil")
  let inputKey = input.nodeKey
  if ctx.contribs.hasKey(inputKey):
    ctx.contribs[inputKey] = input.addLike(ctx.contribs[inputKey], contrib)
  else:
    ctx.contribs[inputKey] = contrib

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
  # `backward` only propagates along raw `inputs`; hidden deps must come from
  # traversal or `backwardTarget`.
  for inputIndex in 0..<node.inputs.len:
    let input = node.inputs[inputIndex]
    if not ctx.plan.nodeNeedsGradient(input):
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

proc accumulateNodeGradient(ctx: var GradBuildContext, node: Gvalue) =
  let isRoot = sameNode(node, ctx.dep)
  let nodeKey = node.nodeKey
  if not isRoot and not ctx.contribs.hasKey(nodeKey):
    return
  let upstream = if isRoot: nil else: ctx.contribs[nodeKey]
  # `backwardTarget` overrides ordinary backward propagation for this node.
  # A nil contribution is an explicit zero contribution for the target.
  let graphFunc = node.gfunc
  if graphFunc != nil and graphFunc.backwardTarget != nil:
    let contribution = graphFunc.backwardTarget(upstream, node, ctx.x, ctx.dep)
    if contribution != nil:
      ctx.addGradContribution(ctx.x, contribution)
    return
  ctx.propagateInputGradients(node, upstream)

proc executeGradPlan(ctx: var GradBuildContext): Gvalue =
  for j in countdown(ctx.plan.order.high, 0):
    ctx.accumulateNodeGradient(ctx.plan.order[j])
  let xKey = ctx.x.nodeKey
  if ctx.contribs.hasKey(xKey):
    return ctx.contribs[xKey]
  ctx.x.zeroLike

proc gradIsolated*(dep: Gvalue, x: Gvalue): Gvalue

proc gradImpl(dep: Gvalue, x: Gvalue): Gvalue =
  discard requireSameGraphRuntime(dep, x, "grad", "output", "input")
  var ctx = initGradBuildContext(dep, x)
  let signature = ctx.dep.buildGradSignature
  let cache = ctx.dep.ensureGradCacheEntry()
  let grt = ctx.dep.runtime
  let inputId = ctx.x.stableNodeId

  if cache.hasSignature and cache.signature == signature:
    grt.gradCacheStats.signatureHits.inc
    if cache.grads.hasKey(inputId):
      grt.gradCacheStats.directHits.inc
      return cache.grads[inputId]
    grt.gradCacheStats.directMisses.inc
  else:
    grt.gradCacheStats.signatureMisses.inc
    if cache.hasSignature:
      grt.gradCacheStats.invalidations.inc
    cache.hasSignature = true
    cache.signature = signature
    cache.grads = initTable[NodeId, Gvalue]()

  if sameNode(ctx.dep, ctx.x):
    result = ctx.x.oneLike
    cache.grads[inputId] = result
    return

  ctx.plan = buildGradPlan(dep, x)
  if ctx.plan.nodeNeedsGradient(dep):
    result = ctx.executeGradPlan
  else:
    result = x.zeroLike
  cache.grads[inputId] = result

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
