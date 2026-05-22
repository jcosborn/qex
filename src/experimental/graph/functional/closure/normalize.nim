import std/[sets, tables]
import ../../core
import ../../core/base
import ../callable
import ../callable/walk
import clone

type
  CaptureCtx = object
    bound: NodeSet
    seenCaps: NodeSet
    seenNodes: NodeSet
    captures: seq[Gvalue]

proc initCaptureCtx(boundValue: Gvalue): CaptureCtx =
  result.bound = initHashSet[NodeKey]()
  result.seenCaps = initHashSet[NodeKey]()
  result.seenNodes = initHashSet[NodeKey]()
  result.bound.incl boundValue.nodeKey

proc isClosureLeaf(v: Gvalue): bool =
  v.gfunc == nil and v.inputs.len == 0

proc addCaptureLeaf(ctx: var CaptureCtx, v: Gvalue) =
  let key = v.nodeKey
  if key notin ctx.seenCaps and not (v of Glambda):
    ctx.seenCaps.incl key
    ctx.captures.add v

proc collectCaptureValues(v: Gvalue, ctx: var CaptureCtx)

proc collectClosureCaptureChildren(v: Gvalue,
                                   ctx: var CaptureCtx): bool =
  if v of Glambda:
    let fn = Glambda(v)
    if fn.isResolvedClosure:
      for binding in fn.env:
        collectCaptureValues(binding.value, ctx)
      return true

  let boundCallable = v.freshCallableBound
  if boundCallable != nil:
    collectCaptureValues(boundCallable, ctx)
    return true

  false

proc collectCaptureValues(v: Gvalue, ctx: var CaptureCtx) =
  let key = v.nodeKey
  if key in ctx.seenNodes:
    return
  ctx.seenNodes.incl key
  if key in ctx.bound:
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

proc initLambdaSubst(fn: Glambda, arg: Gvalue = nil): Bindings =
  result = initTable[NodeKey, Gvalue]()
  for binding in fn.env:
    result[binding.param.nodeKey] = binding.value
  if arg != nil:
    result[fn.param.nodeKey] = arg

proc initCaptureParamBindings(captures: openArray[Gvalue]): tuple[
    env: seq[LambdaBinding], subst: Bindings] =
  result.env = newseq[LambdaBinding](captures.len)
  result.subst = initTable[NodeKey, Gvalue]()
  for i in 0..<captures.len:
    let param = captures[i].newOneOf
    result.env[i] = LambdaBinding(param: param, value: captures[i])
    result.subst[captures[i].nodeKey] = param

proc materializeLambdaEnvBody(fn: Glambda): Gvalue =
  if fn.env.len == 0:
    return fn.body
  let subst = fn.initLambdaSubst
  fn.body.cloneWithFreshMemo(subst)

proc bindLambdaCaptureParams(fn: Glambda, captures: seq[Gvalue]) =
  if captures.len == 0:
    fn.env.setLen(0)
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
  fn.env.setLen(0)
  let captures = fn.collectLambdaBodyCaptures
  fn.bindLambdaCaptureParams(captures)

proc instantiateBody(fn: Glambda, x: Gvalue): Gvalue =
  let subst = fn.initLambdaSubst(x)
  fn.body.cloneWithFreshMemo(subst)

proc instantiateNormalizedBody*(fn: Glambda, x: Gvalue): Gvalue =
  result = fn.instantiateBody(x)
  if result of Glambda and Glambda(result).isResolvedClosure:
    normalizeClosure(Glambda(result))
