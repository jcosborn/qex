import ../../core
import ../callable
import clone

type
  CaptureCtx = object
    bound: NodeSet
    seenCaps: NodeSet
    seenNodes: NodeSet
    captures: seq[Gvalue]

proc initCaptureCtx(boundValue: Gvalue): CaptureCtx =
  result.bound = initNodeSet()
  result.seenCaps = initNodeSet()
  result.seenNodes = initNodeSet()
  if boundValue != nil:
    result.bound.inclNode boundValue

proc isClosureLeaf(v: Gvalue): bool =
  v != nil and v.gfunc == nil and v.inputs.len == 0

proc addCaptureLeaf(ctx: var CaptureCtx, v: Gvalue) =
  if not ctx.seenCaps.containsNode(v) and v.lambdaValue == nil:
    ctx.seenCaps.inclNode v
    ctx.captures.add v

proc collectCaptureValues(v: Gvalue, ctx: var CaptureCtx)

proc collectClosureCaptureChildren(v: Gvalue,
                                   ctx: var CaptureCtx): bool =
  let fn = v.resolvedClosure
  if fn != nil:
    for binding in fn.env:
      collectCaptureValues(binding.value, ctx)
    return true

  let boundCallable = v.freshCallableBound
  if boundCallable != nil:
    collectCaptureValues(boundCallable, ctx)
    return true

  false

proc collectCaptureValues(v: Gvalue, ctx: var CaptureCtx) =
  if v == nil:
    return
  if ctx.seenNodes.containsNode(v):
    return
  ctx.seenNodes.inclNode v
  if ctx.bound.containsNode(v):
    return

  if v.collectClosureCaptureChildren(ctx):
    return

  if v.isClosureLeaf:
    ctx.addCaptureLeaf(v)
    return

  for child in v.collectNodeInputs(iwmDepend):
    collectCaptureValues(child, ctx)

proc collectLambdaBodyCaptures(fn: Glambda): seq[Gvalue] =
  var ctx = initCaptureCtx(fn.param)
  collectCaptureValues(fn.body, ctx)
  ctx.captures

proc initLambdaSubst*(fn: Glambda, arg: Gvalue = nil): Bindings =
  result = initBindings()
  for binding in fn.env:
    result.bindNode(binding.param, binding.value)
  if arg != nil:
    result.bindNode(fn.param, arg)

proc initCaptureParamBindings*(captures: openArray[Gvalue]): tuple[
    env: seq[LambdaBinding], subst: Bindings] =
  result.env = newseq[LambdaBinding](captures.len)
  result.subst = initBindings()
  for i in 0..<captures.len:
    let param = localValue(captures[i])
    result.env[i] = LambdaBinding(param: param, value: captures[i])
    result.subst.bindNode(captures[i], param)

proc materializeLambdaEnvBody*(fn: Glambda): Gvalue =
  if fn.env.len == 0:
    return fn.body
  let subst = fn.initLambdaSubst
  fn.body.cloneWithFreshMemo(subst)

proc bindLambdaCaptureParams*(fn: Glambda, captures: seq[Gvalue]) =
  fn.env.setLen(0)
  if captures.len == 0:
    return

  let captureBindings = initCaptureParamBindings(captures)
  fn.env = captureBindings.env
  fn.body = fn.body.cloneWithFreshMemo(captureBindings.subst)

proc normalizeClosure*(fn: Glambda) =
  ## Normalize a resolved lambda into closure form.
  ## After this pass, captures live in `env` and the body can be instantiated by
  ## ordinary substitution without having to rediscover free variables.
  if not fn.isResolvedClosure:
    return

  fn.body = materializeLambdaEnvBody(fn)
  let captures = fn.collectLambdaBodyCaptures
  fn.bindLambdaCaptureParams(captures)

proc instantiateBody*(fn: Glambda, x: Gvalue): Gvalue =
  let subst = fn.initLambdaSubst(x)
  fn.body.cloneWithFreshMemo(subst)

proc instantiateNormalizedBody*(fn: Glambda, x: Gvalue): Gvalue =
  result = fn.instantiateBody(x)
  if result of Glambda and Glambda(result).isResolvedClosure:
    normalizeClosure(Glambda(result))
