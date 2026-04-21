import std/tables
import base, traverse
import grad_signature, grad_cache, grad_build

proc gradIsolated*(dep: Gvalue, x: Gvalue): Gvalue

proc gradOrZero*(dep: Gvalue, x: Gvalue): Gvalue =
  if dep == nil or x == nil:
    raiseError("gradOrZero has nil input")
  result = dep.gradIsolated(x)
  if result == nil:
    result = x.zeroLike

proc prepareGradCache(ctx: GradBuildContext,
                      cache: var GradCacheEntry): Gvalue =
  let signature = ctx.dep.buildGradSignature
  cache = ctx.dep.ensureGradCacheEntry()
  let grt = ctx.dep.runtime

  if cache.hasSignature and cache.signature == signature:
    grt.gradCacheStats.signatureHits.inc
    if cache.hasCachedGrad(ctx.x):
      result = cache.cachedGrad(ctx.x)
      grt.gradCacheStats.directHits.inc
      return
    grt.gradCacheStats.directMisses.inc
  else:
    grt.gradCacheStats.signatureMisses.inc
    if cache.hasSignature:
      grt.gradCacheStats.invalidations.inc
    cache.hasSignature = true
    cache.signature = signature
    cache.grads = initTable[NodeId, Gvalue]()

  if sameNode(ctx.dep, ctx.x):
    result = ctx.x.constLike(1)
    cache.cacheGrad(ctx.x, result)

proc gradImpl(dep: Gvalue, x: Gvalue): Gvalue =
  let depGrt = dep.runtime
  let xGrt = x.runtime
  # Grad caches are runtime-owned; never build a cross-runtime derivative.
  if depGrt != xGrt:
    raiseValueError("grad mixes multiple graph runtimes")
  var ctx = initGradBuildContext(dep, x)
  var cache: GradCacheEntry
  result = ctx.prepareGradCache(cache)
  if result != nil:
    return

  ctx.plan = buildGradPlan(dep, x)
  if ctx.plan.nodeNeedsGradient(dep):
    result = ctx.executeGradPlan
  else:
    result = x.zeroLike
  cache.cacheGrad(ctx.x, result)

proc grad*(dep: Gvalue, x: Gvalue): Gvalue =
  gradImpl(dep, x)

proc gradIsolated*(dep: Gvalue, x: Gvalue): Gvalue =
  withIsolatedGradCache(dep.runtime, proc(): Gvalue =
    gradImpl(dep, x))
