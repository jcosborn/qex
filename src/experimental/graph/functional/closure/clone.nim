import std/tables
import ../../core
import ../../core/base
import ../callable
import ../callable/walk

type
  CloneCtx = object
    subst: Bindings
    memo: Bindings

proc initCloneCtx(subst: Bindings): CloneCtx =
  result.subst = subst
  result.memo = initTable[NodeKey, Gvalue]()

proc clone(v: Gvalue, ctx: var CloneCtx): Gvalue

proc cloneInputs(src: Gvalue,
                 ctx: var CloneCtx): seq[Gvalue] =
  if src.inputs.len == 0:
    return @[]
  result = newseq[Gvalue](src.inputs.len)
  for i in 0..<src.inputs.len:
    result[i] = clone(src.inputs[i], ctx)

proc withoutLambdaBindings(subst: Bindings, fn: Glambda): Bindings =
  result = subst
  if result.hasKey(fn.param.nodeKey):
    result.del(fn.param.nodeKey)
  for binding in fn.env:
    if result.hasKey(binding.param.nodeKey):
      result.del(binding.param.nodeKey)

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
  result = Glambda(param: fn.param).attachRuntime(fn.runtime)
  ctx.memo[fn.nodeKey] = result

  let r = Glambda(result)
  var innerCtx = initCloneCtx(ctx.subst.withoutLambdaBindings(fn))
  r.env = cloneEnv(fn.env, ctx)
  r.body = clone(fn.body, innerCtx)

proc cloneWithFreshMemo*(v: Gvalue, subst: Bindings): Gvalue =
  var ctx = initCloneCtx(subst)
  clone(v, ctx)

proc clone(v: Gvalue, ctx: var CloneCtx): Gvalue =
  let key = v.nodeKey
  if ctx.subst.hasKey(key):
    return ctx.subst[key]
  if ctx.memo.hasKey(key):
    return ctx.memo[key]

  if v of Glambda:
    let fn = Glambda(v)
    if fn.isResolvedClosure:
      return fn.cloneResolvedLambda(ctx)

  let boundCallable = v.freshCallableBound
  if boundCallable != nil:
    result = clone(boundCallable, ctx)
    ctx.memo[key] = result
    return result

  if v.gfunc == nil and v.inputs.len == 0:
    return v

  # Generic graph-node clone is valid only for nodes whose structural state is
  # fully represented by `inputs` plus `gfunc`. Nodes with extra structural
  # fields need explicit clone handling or prototype state preserved by
  # `newOneOf`.
  result = v.newOneOf
  ctx.memo[key] = result
  result.inputs = v.cloneInputs(ctx)
  result.gfunc = v.gfunc
