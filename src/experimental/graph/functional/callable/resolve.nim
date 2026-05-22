import ../../core
import ../../core/base
import types
import values
import walk

proc resolveCallableChain*(v: Gvalue, mode: CallableResolveMode): Gvalue =
  discard v.requireGraphValue("callable resolution")
  let grt = v.runtime
  var current = v
  var depth = 0
  while depth < grt.lambdaResolveDepthLimit:
    let next = current.symbolicWrapperBinding
    if next != nil:
      current = next
      inc depth
      continue
    let graphFunc = current.gfunc
    if mode == crmReduced and graphFunc != nil and graphFunc.reduceCallable != nil:
      let reduced = graphFunc.reduceCallable(current)
      if reduced == nil:
        raiseValueError("callable reduction returned nil:\n" & current.nodeRepr)
      current = reduced
      inc depth
      continue
    return current
  raiseValueError(
    "callable resolution exceeded depth limit " & $grt.lambdaResolveDepthLimit &
    ":\n" & v.nodeRepr)

proc resolveDirectLambda*(fun: Gvalue): Glambda =
  let resolved = fun.resolveCallableChain(crmShallow)
  if resolved of Glambda:
    result = Glambda(resolved)
    if not result.isResolvedClosure:
      result = nil

proc resolveLambda*(fun: Gvalue): Glambda =
  let resolved = fun.resolveCallableChain(crmReduced)
  if resolved of Glambda:
    result = Glambda(resolved)
    if not result.isResolvedClosure:
      result = nil

proc symbolicCallableToken*(v: Gvalue): uint64 =
  let directValue = v.resolveCallableChain(crmShallow)
  if directValue of Glambda:
    let directFn = Glambda(directValue)
    if directFn.isResolvedClosure:
      return directFn.stableNodeId
  if directValue != v:
    return directValue.stableNodeId
  0

proc materializeApplyResultProto(proto: Gvalue,
                                 remainingDepth: int): Gvalue =
  let checkedProto = proto.requireGraphValue("apply result prototype")
  let nestedProto = checkedProto.bindingResultProto
  if nestedProto != nil:
    if remainingDepth <= 0:
      raiseValueError(
        "apply result prototype exceeded depth limit " &
        $checkedProto.runtime.lambdaResolveDepthLimit &
        ":\n" & checkedProto.nodeRepr)
    return callableWrapperNode(
      materializeApplyResultProto(nestedProto, remainingDepth - 1))
  checkedProto.newOneOf

proc applyResultProto*(fun: Gvalue): Gvalue =
  let checkedFun = fun.requireGraphValue("apply result prototype")
  let retProto = checkedFun.bindingResultProto
  if retProto == nil:
    return nil
  materializeApplyResultProto(retProto, checkedFun.runtime.lambdaResolveDepthLimit)

method appendSignatureTokens*(v: Gwrapper, tokens: var seq[GradSigToken]) =
  let token = symbolicCallableToken(v)
  if token != 0:
    tokens.add GradSigToken(kind: gstCallable, key: token)
