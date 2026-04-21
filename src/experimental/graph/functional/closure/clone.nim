import ../../core
import ../callable

type
  CloneCtx* = object
    subst: Bindings
    memo: Bindings

proc initCloneCtx*(subst: Bindings): CloneCtx =
  result.subst = subst
  result.memo = initBindings()

proc clone*(v: Gvalue, ctx: var CloneCtx): Gvalue

proc cloneInputs(dst: Gvalue,
                 src: Gvalue,
                 ctx: var CloneCtx) =
  if src.inputs.len > 0:
    dst.inputs = newseq[Gvalue](src.inputs.len)
    for i in 0..<src.inputs.len:
      dst.inputs[i] = clone(src.inputs[i], ctx)

proc withoutLambdaBindings*(subst: Bindings, fn: Glambda): Bindings =
  result = subst
  result.deleteBinding(fn.param)
  for binding in fn.env:
    result.deleteBinding(binding.param)

proc cloneEnv(env: openArray[LambdaBinding],
              ctx: var CloneCtx): seq[LambdaBinding] =
  if env.len == 0:
    return @[]
  result = newseq[LambdaBinding](env.len)
  for i in 0..<env.len:
    result[i] = LambdaBinding(
      param: env[i].param,
      value: clone(env[i].value, ctx))

proc cloneResolvedLambda(fn: Glambda,
                         ctx: var CloneCtx): Gvalue =
  result = Glambda(
    param: fn.param,
    runtime: fn.runtime)
  ctx.memo.putNode(fn, result)

  let r = Glambda(result)
  r.env = cloneEnv(fn.env, ctx)
  if fn.body != nil:
    var innerCtx = initCloneCtx(ctx.subst.withoutLambdaBindings(fn))
    r.body = clone(fn.body, innerCtx)

proc cloneWithFreshMemo*(v: Gvalue, subst: Bindings): Gvalue =
  var ctx = initCloneCtx(subst)
  clone(v, ctx)

proc clone*(v: Gvalue, ctx: var CloneCtx): Gvalue =
  if v == nil:
    return nil

  if ctx.subst.hasNode(v):
    return ctx.subst.getNode(v)
  if ctx.memo.hasNode(v):
    return ctx.memo.getNode(v)

  let fn = v.resolvedClosure
  if fn != nil:
    return fn.cloneResolvedLambda(ctx)

  let boundCallable = v.freshCallableBound
  if boundCallable != nil:
    result = clone(boundCallable, ctx)
    ctx.memo.putNode(v, result)
    return result

  if v.gfunc == nil and v.inputs.len == 0:
    return v

  result = v.newOneOf
  ctx.memo.putNode(v, result)
  result.cloneInputs(v, ctx)
  result.gfunc = v.gfunc
