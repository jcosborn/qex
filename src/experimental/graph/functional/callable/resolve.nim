import ../../core
import types
import values
import walk

proc nextCallableBinding(v: Gvalue): Gvalue =
  v.symbolicWrapperBinding

proc resolveCallableChain*(v: Gvalue, mode: CallableResolveMode): Gvalue =
  let grt = v.runtime
  var current = v
  var depth = 0
  while depth < grt.lambdaResolveDepthLimit:
    if current == nil:
      return nil
    let next = current.nextCallableBinding
    if next != nil:
      current = next
      inc depth
      continue
    if mode == crmReduced and current.isDeferredApplyValue:
      try:
        current = current.reduceDeferredApplyValue
      except GraphUnresolvedValueError:
        return nil
      inc depth
      continue
    return current
  raiseValueError(
    "callable resolution exceeded depth limit " & $grt.lambdaResolveDepthLimit &
    ":\n" & v.nodeRepr)

proc inspectCallable*(v: Gvalue,
                      mode: CallableResolveMode = crmShallow): CallableInspect =
  result.directValue = v.resolveCallableChain(crmShallow)
  result.directFn = result.directValue.resolvedClosure
  if mode == crmReduced:
    result.reducedValue = v.resolveCallableChain(crmReduced)
    result.reducedFn = result.reducedValue.resolvedClosure
    result.hasReduced = true

proc resolveDirectLambda*(fun: Gvalue): Glambda =
  inspectCallable(fun).directFn

proc resolveLambda*(fun: Gvalue): Glambda =
  fun.inspectCallable(crmReduced).reducedFn

proc resolvedCallableValue*(v: Gvalue): Gvalue =
  ## Return the resolved callable after following wrappers and reducible applies.
  let inspect = inspectCallable(v, crmReduced)
  if inspect.directValue == nil:
    return nil
  if inspect.hasReduced and inspect.reducedValue != nil:
    return inspect.reducedValue
  inspect.directValue

proc resultProto*(v: Gvalue): Gvalue =
  v.bindingResultProto

proc symbolicCallableToken*(v: Gvalue): uint64 =
  let inspect = inspectCallable(v)
  if inspect.directValue == nil:
    return 0
  if inspect.directFn != nil:
    return callableKey(inspect.directFn)
  if inspect.directValue != v:
    return inspect.directValue.stableNodeId
  0

proc appendCallableToken*(v: Gvalue, tokens: var seq[GradSigToken]) =
  let token = symbolicCallableToken(v)
  if token != 0:
    tokens.add GradSigToken(kind: gstCallable, key: token)

proc isCallableLike*(v: Gvalue): bool =
  v.resultProto != nil

proc materializeApplyResultProto(proto: Gvalue): Gvalue =
  if proto == nil:
    return nil
  let nestedProto = proto.resultProto
  if nestedProto != nil:
    return callableWrapperNode(materializeApplyResultProto(nestedProto))
  proto.newOneOf

proc applyResultProto*(fun: Gvalue): Gvalue =
  let retProto = fun.resultProto
  if retProto == nil:
    return nil
  materializeApplyResultProto(retProto)

method appendSignatureTokens*(v: Gwrapper, tokens: var seq[GradSigToken]) =
  appendCallableToken(v, tokens)
