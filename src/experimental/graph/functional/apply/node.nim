import ../../core
import ../../scalar
import ../callable, meta, partial, reduce_cache

proc isDeferredApplyNode(v: Gvalue): bool
proc applyDeferredBackwardTarget(zb: Gvalue,
                                 z: Gvalue,
                                 target: Gvalue,
                                 dep: Gvalue): Gvalue
proc applyPartialDeferredBackwardTarget(zb: Gvalue,
                                        z: Gvalue,
                                        target: Gvalue,
                                        dep: Gvalue): Gvalue

proc reduceApplyDeferredValue(v: Gvalue): Gvalue =
  v.reducedApplyExpr

proc applyDeferredf(v: Gvalue) =
  let reduced = v.reducedApplyExpr
  let callable = resolveLambda(reduced)
  if callable != nil:
    v.valCopy callable
    return
  discard reduced.eval
  v.valCopy reduced.resolvedCallableValue

proc applyDeferredWalkInputs(v: Gvalue,
                             visit: proc(n: Gvalue) {.closure.}) =
  v.walkApplyDeps(visit)

proc applyPartialDeferredf(v: Gvalue) =
  let expr = v.reduceApplyPartialExpr
  discard expr.eval
  v.valCopy expr

proc applyPartialDeferredWalkEvalInputs(v: Gvalue,
                                        visit: proc(n: Gvalue) {.closure.}) =
  ## Deferred partials stay lazy with respect to targets.
  walkDeferredEvalGraph(v.applyPartialBase, isDeferredApplyNode, visit)

proc applyPartialDeferredWalkGraphInputs(v: Gvalue,
                                         visit: proc(n: Gvalue) {.closure.}) =
  visit(v.applyPartialBase)

let gapplyPartialDeferred = newGfunc(
  forward = applyPartialDeferredf,
  depWalks = newDepWalks(
    eval = walkedInputs(applyPartialDeferredWalkEvalInputs),
    gradSignature = walkedInputs(applyPartialDeferredWalkGraphInputs),
    depend = walkedInputs(applyPartialDeferredWalkGraphInputs)),
  backwardTarget = applyPartialDeferredBackwardTarget,
  deferredApplyKind = dakApplyPartial,
  name = "applyPartialDeferred")
let gapplyDeferred = newGfunc(
  forward = applyDeferredf,
  depWalks = newDepWalks(
    eval = walkedInputs(applyDeferredWalkInputs),
    gradSignature = walkedInputs(applyDeferredWalkInputs),
    depend = walkedInputs(applyDeferredWalkInputs)),
  backwardTarget = applyDeferredBackwardTarget,
  reduceValue = reduceApplyDeferredValue,
  signature = appendApplySignature,
  deferredApplyKind = dakApply,
  name = "applyDeferred")

proc isDeferredApplyNode(v: Gvalue): bool =
  v != nil and v.gfunc != nil and v.gfunc.deferredApplyKind in {dakApply, dakApplyPartial}

proc applyDeferredContribution(zb: Gvalue, z: Gvalue, target: Gvalue): Gvalue =
  if not target.needsConcreteApplyPartial:
    return nil
  if not zb.acceptsApplyPartialScale:
    return nil
  let partial = newApplyPartialNode(
    gapplyPartialDeferred,
    z,
    target,
    "applyPartialDeferred")
  if zb == nil:
    return partial
  zb * partial

proc applyDeferredBackwardTarget(zb: Gvalue,
                                 z: Gvalue,
                                 target: Gvalue,
                                 dep: Gvalue): Gvalue =
  discard dep
  applyDeferredContribution(zb, z, target)

proc applyPartialDeferredBackwardTarget(zb: Gvalue,
                                        z: Gvalue,
                                        target: Gvalue,
                                        dep: Gvalue): Gvalue =
  discard dep
  contributeApplyPartialTarget(
    gapplyPartialDeferred,
    zb,
    z,
    target,
    "applyPartialDeferred")

proc directLambdaArgCompatible(param: Gvalue, arg: Gvalue): bool =
  # Callable wrapper parameters have their own binding semantics. Keep this
  # fast path to ordinary value parameters until callable shape gets a type
  # object that can describe higher-order arguments precisely.
  if param.wrapper != nil:
    return true
  param.copyCompatible(arg)

proc apply*(fun: Gvalue, x: Gvalue): Gvalue =
  if fun == nil:
    raiseValueError("apply function cannot be nil")
  if x == nil:
    raiseValueError("apply argument cannot be nil")

  let proto = applyResultProto(fun)
  if proto == nil:
    raiseValueError("apply expects a callable value or callable placeholder, got: " & fun.nodeRepr)
  let direct = resolveDirectLambda(fun)
  if direct != nil and not direct.param.directLambdaArgCompatible(x):
    raiseValueError(
      "apply argument is incompatible with lambda parameter" &
      "\nparameter: " & direct.param.nodeRepr &
      "\nargument: " & x.nodeRepr)

  graphNode(proto.newOneOf, @[fun, x], gapplyDeferred, "apply")

proc applyLiteralArg(fun: Gvalue,
                     value: float): Gvalue =
  if fun == nil:
    raiseValueError("apply function cannot be nil")
  let direct = resolveDirectLambda(fun)
  if direct != nil:
    return scalarLeafLike(direct.param, value)
  scalarLeafLike(fun, value)

proc applyLiteralArg(fun: Gvalue,
                     value: int): Gvalue =
  if fun == nil:
    raiseValueError("apply function cannot be nil")
  let direct = resolveDirectLambda(fun)
  if direct != nil:
    return numericLeafLike(direct.param, value)
  intLeafLike(fun, value)

proc apply*(fun: Gvalue,
            x: float): Gvalue =
  apply(fun, applyLiteralArg(fun, x))

proc apply*(fun: Gvalue,
            x: int): Gvalue =
  apply(fun, applyLiteralArg(fun, x))
