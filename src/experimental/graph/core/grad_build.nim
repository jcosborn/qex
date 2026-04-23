import base, traverse

type
  GradBuildPlan* = object
    relevant*: NodeTable[bool]
    order*: seq[Gvalue]
  GradBuildContext* = object
    dep*: Gvalue
    x*: Gvalue
    plan*: GradBuildPlan
    contribs*: NodeTable[Gvalue]

proc initGradBuildPlan*(): GradBuildPlan =
  result.relevant = initNodeTable[bool]()

proc initGradBuildContext*(dep: Gvalue, x: Gvalue): GradBuildContext =
  result.dep = dep
  result.x = x
  result.plan = initGradBuildPlan()
  result.contribs = initNodeTable[Gvalue]()

proc addGradContribution(ctx: var GradBuildContext,
                         input: Gvalue,
                         contrib: Gvalue) =
  if input == nil or contrib == nil:
    return
  if ctx.contribs.hasNode(input):
    ctx.contribs.putNode(input, input.addLike(ctx.contribs.getNode(input), contrib))
  else:
    ctx.contribs.putNode(input, contrib)

proc nodeNeedsGradient*(plan: GradBuildPlan, node: Gvalue): bool =
  plan.relevant.nodeOrDefault(node, false)

proc buildGradPlan*(dep: Gvalue, x: Gvalue): GradBuildPlan =
  var plan = initGradBuildPlan()
  var active = initNodeSet()

  proc markGradientNode(node: Gvalue): bool =
    if node == nil:
      return false
    if plan.relevant.hasNode(node):
      return plan.nodeNeedsGradient(node)
    if active.containsNode(node):
      return false

    active.inclNode node
    var need = sameNode(node, x)
    node.walkDependDeps(proc(child: Gvalue) =
      if markGradientNode(child):
        need = true)
    active.exclNode node

    plan.relevant.putNode(node, need)
    if need:
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
    if input == nil:
      raiseError("node has nil input at index " & $inputIndex & ":\n" & node.nodeRepr)
    if not ctx.plan.nodeNeedsGradient(input):
      continue
    if graphFunc.backward == nil:
      raiseError(
        node.nodeRepr & ":" & $inputIndex & ":" & input.nodeRepr &
        ": backward undefined")
    ctx.addGradContribution(
      input,
      graphFunc.backward(upstream, node, inputIndex, ctx.dep))

proc accumulateNodeGradient(ctx: var GradBuildContext, node: Gvalue) =
  let isRoot = sameNode(node, ctx.dep)
  if not isRoot and not ctx.contribs.hasNode(node):
    return
  let upstream = if isRoot: nil else: ctx.contribs.getNode(node)
  # `backwardTarget` overrides ordinary backward propagation for this node.
  if node.hasTargetBackward:
    ctx.addGradContribution(
      ctx.x,
      node.runTargetBackward(upstream, ctx.x, ctx.dep))
    return
  ctx.propagateInputGradients(node, upstream)

proc executeGradPlan*(ctx: var GradBuildContext): Gvalue =
  for j in countdown(ctx.plan.order.high, 0):
    ctx.accumulateNodeGradient(ctx.plan.order[j])
  if ctx.contribs.hasNode(ctx.x):
    return ctx.contribs.getNode(ctx.x)
  ctx.x.zeroLike
